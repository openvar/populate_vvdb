"""Tests for vvDBFlatBuildDbSNP.py (format, normalization, reader, builder, CLI)."""

from __future__ import annotations

import os

import pytest

import vvDBFlatBuildDbSNP as B


# ---------------------------------------------------------------------------
# normalization: trim_alleles
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("chrom,pos,ref,alt,expected", [
    ("1", 100, "C", "T", ("1", 100, "C", "T")),        # SNV unchanged
    ("1", 250, "AA", "GA", ("1", 250, "A", "G")),      # shared suffix trimmed
    ("1", 250, "TA", "TG", ("1", 251, "A", "G")),      # shared prefix advances POS
    ("1", 100, "TAT", "TGT", ("1", 101, "A", "G")),    # prefix + suffix
    ("1", 100, "c", "t", ("1", 100, "C", "T")),        # upper-cased
    ("1", 100, "A", "AT", ("1", 100, "A", "AT")),      # insertion stays (>=1 base)
    ("1", 100, "AT", "A", ("1", 100, "AT", "A")),      # deletion stays (>=1 base)
])
def test_trim_alleles(chrom, pos, ref, alt, expected):
    assert B.trim_alleles(chrom, pos, ref, alt) == expected


# ---------------------------------------------------------------------------
# contig maps
# ---------------------------------------------------------------------------
def test_get_contig_map():
    assert B.get_contig_map("grch38")["NC_000001.11"] == "1"
    assert B.get_contig_map("grch37")["NC_000001.10"] == "1"
    assert B.get_contig_map("hg38") is B.GRCH38
    assert B.get_contig_map("38") is B.GRCH38
    assert B.get_contig_map("hg19") is B.GRCH37
    assert B.get_contig_map("unknown") == {}


def test_contig_maps_differ_between_builds():
    # chr1 has a different accession in each build.
    assert B.GRCH38["NC_000001.11"] == "1"
    assert B.GRCH37["NC_000001.10"] == "1"
    assert "NC_000001.11" not in B.GRCH37


# ---------------------------------------------------------------------------
# small helpers
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("raw,expected", [
    ("rs6025", 6025), ("6025", 6025), (6025, 6025),
    ("rs6025;rs1", 6025), ("6025,7", 6025),
    ("", None), (".", None), (None, None), ("abc", None),
])
def test_strip_rs(raw, expected):
    assert B._strip_rs(raw) == expected


@pytest.mark.parametrize("path,expected", [
    ("x.bcf", "bcf"), ("x.tsv", "tsv"), ("x.txt", "tsv"),
    ("x.vcf", "vcf"), ("x.vcf.gz", "vcf"), ("x.vcf.bgz", "vcf"),
])
def test_detect_source_type(path, expected):
    assert B._detect_source_type(path) == expected


# ---------------------------------------------------------------------------
# on-disk format: pack/iter entries, MetaWriter/MetaIndex
# ---------------------------------------------------------------------------
def test_pack_iter_entries_roundtrip():
    pairs = [(1, 100), (2, 200), (65535, 2 ** 40)]
    data = b"".join(B.pack_entry(r, v) for r, v in pairs)
    assert len(B.pack_entry(1, 2)) == B.ENTRY_SIZE == 16
    assert list(B.iter_entries(data)) == pairs


def test_meta_roundtrip_and_candidate_block(tmp_path):
    mw = B.MetaWriter(block_size=4)
    mw.add_block(first_rsid=10, comp_offset=0, comp_len=100, n_entries=4)
    mw.add_block(first_rsid=100, comp_offset=100, comp_len=120, n_entries=4)
    mw.add_block(first_rsid=1000, comp_offset=220, comp_len=90, n_entries=2)
    p = tmp_path / "x.meta"
    with open(p, "wb") as fh:
        mw.write(fh)

    mi = B.MetaIndex.load(str(p))
    assert mi.block_size == 4
    assert mi.n_entries == 10
    assert [b.first_rsid for b in mi.blocks] == [10, 100, 1000]

    # binary search picks the rightmost block with first_rsid <= rsid
    assert mi.candidate_block(5) is None       # before every block
    assert mi.candidate_block(10) == 0
    assert mi.candidate_block(50) == 0
    assert mi.candidate_block(100) == 1
    assert mi.candidate_block(150) == 1
    assert mi.candidate_block(5000) == 2


def test_meta_bad_magic(tmp_path):
    p = tmp_path / "bad.meta"
    p.write_bytes(b"XXXX" + b"\x00" * 60)
    with pytest.raises(ValueError):
        B.MetaIndex.load(str(p))


# ---------------------------------------------------------------------------
# builder outputs
# ---------------------------------------------------------------------------
def test_build_outputs_exist(synth):
    assert set(synth["outputs"]) == {"master", "tabix", "rsid_idx", "rsid_meta"}
    for path in synth["outputs"].values():
        assert os.path.exists(path), path


def test_build_is_idempotent(synth):
    # Building again over the same output dir must not rebuild; returns same paths.
    builder = B.FlatDbSNPBuilder(out_dir=synth["out_dir"])
    again = builder.build("grch38", synth["source"], version="test")
    assert again == synth["outputs"]


# ---------------------------------------------------------------------------
# reader: DbSNPFlat, both directions
# ---------------------------------------------------------------------------
@pytest.fixture
def reader(synth):
    db = B.DbSNPFlat(synth["stem"] + ".variants.tsv.gz",
                     synth["stem"] + ".rsid.idx",
                     synth["stem"] + ".rsid.meta")
    yield db
    db.close()


@pytest.mark.parametrize("rsid", [99, 201, 202, 203, 6025, 334, 429358])
def test_reader_rsid_to_variant(reader, synth, rsid):
    chrom, pos, ref, alts = synth["landmarks"][rsid]
    groups = reader.rsid_to_variants(rsid)
    assert len(groups) == 1
    g = groups[0]
    assert (g.rsid, g.chrom, g.pos, g.ref, tuple(g.alts)) == (rsid, chrom, pos, ref, alts)


def test_reader_resolve_rsid_string(reader):
    assert reader.resolve_rsid(6025) == ["1-169549811-C-T"]
    assert reader.resolve_rsid("rs334") == ["11-5227002-T-A,C,G"]


@pytest.mark.parametrize("rsid", [99, 201, 202, 203, 6025, 334, 429358])
def test_reader_variant_to_rsid_each_alt(reader, synth, rsid):
    chrom, pos, ref, alts = synth["landmarks"][rsid]
    for alt in alts:
        assert rsid in reader.variant_to_rsids(chrom, pos, ref, alt)


def test_reader_variant_to_rsid_singular(reader):
    assert reader.variant_to_rsid("1", 169549811, "C", "T") == 6025


def test_reader_normalizes_query(reader):
    # non-minimal query for rs99 (source AA/GA) resolves via trim_alleles -> A/G
    assert 99 in reader.variant_to_rsids("1", 250, "AA", "GA")


def test_reader_unknown(reader):
    assert reader.rsid_to_variants(999_999_999) == []
    assert reader.variant_to_rsids("1", 12345, "A", "G") == []
    assert reader.variant_to_rsid("1", 12345, "A", "G") is None


def test_variant_group_to_string():
    g = B.VariantGroup(rsid=1, chrom="1", pos=5, ref="C", alts=("T", "G"))
    assert g.to_string() == "1-5-C-T,G"


# ---------------------------------------------------------------------------
# CLI (_main)
# ---------------------------------------------------------------------------
def test_cli_build_from_tsv(tmp_path):
    tsv = tmp_path / "mini.tsv"
    tsv.write_text("1\t100\trs5\tC\tT\n11\t200\trs6\tA\tG\n")
    out = tmp_path / "cli"
    rc = B._main(["--out", str(out), "--version", "clitest", "--grch38", str(tsv)])
    assert rc == 0

    stem = str(out / "clitest" / "dbsnp_grch38_clitest")
    for suffix in (".variants.tsv.gz", ".rsid.idx", ".rsid.meta"):
        assert os.path.exists(stem + suffix)

    db = B.DbSNPFlat(stem + ".variants.tsv.gz", stem + ".rsid.idx",
                     stem + ".rsid.meta")
    try:
        assert db.variant_to_rsid("1", 100, "C", "T") == 5
        assert db.resolve_rsid(6) == ["11-200-A-G"]
    finally:
        db.close()


def test_cli_requires_a_source():
    with pytest.raises(SystemExit):
        B._main([])                      # neither --grch37 nor --grch38


def test_cli_missing_source_file(tmp_path):
    with pytest.raises(SystemExit):
        B._main(["--out", str(tmp_path / "o"), "--grch38",
                 str(tmp_path / "does_not_exist.vcf")])
