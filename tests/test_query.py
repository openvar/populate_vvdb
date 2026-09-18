"""Tests for vvDBFlatQueryDbSNP.py (the backend query API)."""

from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor

import pytest

import vvDBFlatQueryDbSNP as Q


# ---------------------------------------------------------------------------
# input helpers
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("label,expected", [
    ("grch38", "grch38"), ("GRCh38", "grch38"), ("hg38", "grch38"), ("38", "grch38"),
    ("grch37", "grch37"), ("GRCh37", "grch37"), ("hg19", "grch37"), ("37", "grch37"),
    (" grch38 ", "grch38"),
])
def test_normalise_assembly_ok(label, expected):
    assert Q.normalise_assembly(label) == expected


@pytest.mark.parametrize("label", ["hg99", "b38", "", "chr1"])
def test_normalise_assembly_bad(label):
    with pytest.raises(ValueError):
        Q.normalise_assembly(label)


@pytest.mark.parametrize("raw,expected", [
    (6025, 6025), ("6025", 6025), ("rs6025", 6025), ("RS6025", 6025),
    (" rs6025 ", 6025),
])
def test_parse_rsid_ok(raw, expected):
    assert Q.parse_rsid(raw) == expected


@pytest.mark.parametrize("raw", ["abc", "rsX", "", 0, -5, True, "rs-1"])
def test_parse_rsid_bad(raw):
    with pytest.raises(ValueError):
        Q.parse_rsid(raw)


@pytest.mark.parametrize("val,expected", [("C", "C"), (" C ", "C"), ("atg", "atg")])
def test_require_allele_ok(val, expected):
    assert Q._require_allele(val, "ref") == expected


@pytest.mark.parametrize("val", ["", "   "])
def test_require_allele_bad(val):
    with pytest.raises(ValueError):
        Q._require_allele(val, "ref")


def test_store_stem():
    assert Q.store_stem("/base", "grch38", "b157") == "/base/b157/dbsnp_grch38_b157"
    assert Q.store_stem("/base", "hg38", "b157") == "/base/b157/dbsnp_grch38_b157"


# ---------------------------------------------------------------------------
# result objects (pure)
# ---------------------------------------------------------------------------
def test_variant_to_rsid_result_helpers():
    q = Q.VariantQuery("grch38", "1", 5, "C", "T")
    assert q.to_string() == "1-5-C-T"
    r = Q.VariantToRsidResult(q, None, (6025, 7))
    assert r.found and r.first == 6025 and r.labels == ("rs6025", "rs7")
    empty = Q.VariantToRsidResult(q, None, ())
    assert not empty.found and empty.first is None and empty.labels == ()


def test_variant_placement_helpers():
    p = Q.VariantPlacement(334, "11", 5227002, "T", ("A", "C", "G"))
    assert p.label == "rs334"
    assert p.to_string() == "11-5227002-T-A,C,G"
    assert p.to_vcf_tuples() == [
        ("11", 5227002, "T", "A"), ("11", 5227002, "T", "C"), ("11", 5227002, "T", "G")]

    res = Q.RsidToVariantResult("grch38", 334, (p,))
    assert res.found and res.label == "rs334"
    assert res.to_strings() == ["11-5227002-T-A,C,G"]
    assert len(res.to_vcf_tuples()) == 3


# ---------------------------------------------------------------------------
# construction
# ---------------------------------------------------------------------------
@pytest.fixture
def client(synth):
    c = Q.DbSNPFlatQuery({"grch38": synth["stem"]})
    yield c
    c.close()


def test_construct_and_introspect(client):
    assert client.assemblies == ("grch38",)
    assert client.has_assembly("grch38")
    assert client.has_assembly("hg38")
    assert not client.has_assembly("grch37")
    assert not client.has_assembly("nonsense")


def test_construct_empty_raises():
    with pytest.raises(ValueError):
        Q.DbSNPFlatQuery({})


def test_construct_missing_file_raises(tmp_path):
    with pytest.raises(FileNotFoundError):
        Q.DbSNPFlatQuery({"grch38": str(tmp_path / "missing_stem")})


def test_from_directory_picks_present_only(synth):
    c = Q.DbSNPFlatQuery.from_directory(synth["out_dir"], version="test",
                                        assemblies=("grch37", "grch38"))
    try:
        assert c.assemblies == ("grch38",)   # grch37 store absent -> skipped
    finally:
        c.close()


def test_from_directory_require_all_raises(synth):
    with pytest.raises(FileNotFoundError):
        Q.DbSNPFlatQuery.from_directory(synth["out_dir"], version="test",
                                        assemblies=("grch37",), require_all=True)


# ---------------------------------------------------------------------------
# direction 1: variant -> rsID
# ---------------------------------------------------------------------------
def test_variant_to_rsids(client, synth):
    hit = client.variant_to_rsids("grch38", "1", 169549811, "C", "T")
    assert hit.found
    assert hit.first == 6025
    assert hit.rsids == (6025,)
    assert hit.labels == ("rs6025",)


def test_variant_to_rsid_singular(client):
    assert client.variant_to_rsid("grch38", "1", 169549811, "C", "T") == 6025


def test_variant_to_rsids_normalized_field(client):
    # raw AA/GA is echoed back normalized to A/G and still resolves rs99
    hit = client.variant_to_rsids("grch38", "1", 250, "AA", "GA")
    assert hit.first == 99
    assert (hit.normalized.ref, hit.normalized.alt) == ("A", "G")


def test_variant_to_rsids_absent_is_empty(client):
    hit = client.variant_to_rsids("grch38", "1", 12345, "A", "G")
    assert not hit.found and hit.rsids == ()


def test_variant_to_rsids_empty_allele_raises(client):
    with pytest.raises(ValueError):
        client.variant_to_rsids("grch38", "1", 100, "", "G")


# ---------------------------------------------------------------------------
# direction 2: rsID -> variant
# ---------------------------------------------------------------------------
def test_rsid_to_variants_multiallelic(client):
    res = client.rsid_to_variants("grch38", "rs334")
    assert res.found
    assert res.to_strings() == ["11-5227002-T-A,C,G"]
    assert res.to_vcf_tuples() == [
        ("11", 5227002, "T", "A"), ("11", 5227002, "T", "C"), ("11", 5227002, "T", "G")]


def test_rsid_spelling_equivalence(client):
    a = client.rsid_to_variants("grch38", 6025).to_strings()
    b = client.rsid_to_variants("grch38", "6025").to_strings()
    c = client.rsid_to_variants("grch38", "rs6025").to_strings()
    assert a == b == c == ["1-169549811-C-T"]


def test_rsid_to_variants_absent_is_empty(client):
    res = client.rsid_to_variants("grch38", 999_999_999)
    assert not res.found and res.placements == ()


def test_rsid_to_variants_malformed_raises(client):
    with pytest.raises(ValueError):
        client.rsid_to_variants("grch38", "not-an-id")


def test_unconfigured_assembly_raises(client):
    with pytest.raises(KeyError):
        client.rsid_to_variants("grch37", 6025)


# ---------------------------------------------------------------------------
# lifecycle + thread safety
# ---------------------------------------------------------------------------
def test_context_manager(synth):
    with Q.DbSNPFlatQuery({"grch38": synth["stem"]}) as c:
        assert c.variant_to_rsid("grch38", "1", 169549811, "C", "T") == 6025


def test_thread_safe(client):
    def one(_):
        v = client.variant_to_rsid("grch38", "1", 169549811, "C", "T")
        r = client.rsid_to_variants("grch38", 334).to_strings()
        return v, tuple(r)

    with ThreadPoolExecutor(max_workers=8) as ex:
        results = list(ex.map(one, range(64)))
    assert set(results) == {(6025, ("11-5227002-T-A,C,G",))}
