#!/usr/bin/env python3
from __future__ import annotations
import bisect
import os
import struct
import subprocess
import tempfile
import zlib
from dataclasses import dataclass
from typing import BinaryIO, Iterable, Iterator, Optional
import pysam


# ============================================================================
# Contig maps  (RefSeq accession -> chromosome name, per assembly)
# ============================================================================
# name -> accession, per assembly (primary assembly units only)
_GRCH38_NAME_TO_ACC = {
    "1": "NC_000001.11", "2": "NC_000002.12", "3": "NC_000003.12",
    "4": "NC_000004.12", "5": "NC_000005.10", "6": "NC_000006.12",
    "7": "NC_000007.14", "8": "NC_000008.11", "9": "NC_000009.12",
    "10": "NC_000010.11", "11": "NC_000011.10", "12": "NC_000012.12",
    "13": "NC_000013.11", "14": "NC_000014.9", "15": "NC_000015.10",
    "16": "NC_000016.10", "17": "NC_000017.11", "18": "NC_000018.10",
    "19": "NC_000019.10", "20": "NC_000020.11", "21": "NC_000021.9",
    "22": "NC_000022.11", "X": "NC_000023.11", "Y": "NC_000024.10",
    "MT": "NC_012920.1",
}

_GRCH37_NAME_TO_ACC = {
    "1": "NC_000001.10", "2": "NC_000002.11", "3": "NC_000003.11",
    "4": "NC_000004.11", "5": "NC_000005.9", "6": "NC_000006.11",
    "7": "NC_000007.13", "8": "NC_000008.10", "9": "NC_000009.11",
    "10": "NC_000010.10", "11": "NC_000011.9", "12": "NC_000012.11",
    "13": "NC_000013.10", "14": "NC_000014.8", "15": "NC_000015.9",
    "16": "NC_000016.9", "17": "NC_000017.10", "18": "NC_000018.9",
    "19": "NC_000019.9", "20": "NC_000020.10", "21": "NC_000021.8",
    "22": "NC_000022.10", "X": "NC_000023.10", "Y": "NC_000024.9",
    "MT": "NC_012920.1",
}

# accession -> name (what the builder uses)
GRCH38 = {acc: name for name, acc in _GRCH38_NAME_TO_ACC.items()}
GRCH37 = {acc: name for name, acc in _GRCH37_NAME_TO_ACC.items()}


def get_contig_map(assembly: str) -> dict:
    """
    Return the accession->name map for an assembly label. Accepts 'GRCh37',
    'grch38', 'hg19', 'hg38', '37', '38', etc. Unknown labels (e.g. 'test')
    return an empty map, so plain-named inputs pass through unchanged.
    """
    a = assembly.strip().lower()
    if a in ("grch38", "hg38", "38"):
        return GRCH38
    if a in ("grch37", "hg19", "37"):
        return GRCH37
    return {}


# ============================================================================
# On-disk format
# ============================================================================
MAGIC = b"VVR1"
FORMAT_VERSION = 1
DEFAULT_BLOCK_SIZE = 4096

# One index entry: rsID (uint64), master BGZF virtual offset (uint64).
_ENTRY = struct.Struct("<QQ")
ENTRY_SIZE = _ENTRY.size  # 16

# Meta header: magic, format version, block_size, n_blocks, n_entries.
_META_HEADER = struct.Struct("<4sHIQQ")

# One meta record per block: first_rsid, comp_offset, comp_len, n_entries.
_META_RECORD = struct.Struct("<QQII")

# Master TSV columns.
C_CHROM, C_POS, C_END, C_RSID, C_REF, C_ALT = range(6)


# ============================================================================
# Normalization seam
# ============================================================================
def trim_alleles(chrom: str, pos, ref: str, alt: str):
    """
    Reference-free minimal-representation trimming: drop shared suffix bases,
    then shared prefix bases (advancing POS), always keeping >=1 base in each
    allele. This is the default normalizer and makes many equivalent VCF
    spellings match (e.g. ``1-250-AA-GA`` -> ``1-250-A-G``).

    It is NOT full left-alignment: shifting an indel left through a homopolymer
    needs the reference sequence. For that, pass VariantValidator's own
    reference-aware normalizer (hgvs / ``bcftools norm -f``) to both the builder
    and :class:`DbSNPFlat` instead - the ONLY requirement is that the same
    normalizer is used on the stored data and on the query.
    """
    pos = int(pos)
    ref, alt = ref.upper(), alt.upper()
    while len(ref) > 1 and len(alt) > 1 and ref[-1] == alt[-1]:
        ref, alt = ref[:-1], alt[:-1]
    while len(ref) > 1 and len(alt) > 1 and ref[0] == alt[0]:
        ref, alt, pos = ref[1:], alt[1:], pos + 1
    return chrom, pos, ref, alt


def pack_entry(rsid: int, voffset: int) -> bytes:
    """Serialise a single ``(rsid, voffset)`` index entry."""
    return _ENTRY.pack(rsid, voffset)


def iter_entries(raw: bytes):
    """Yield ``(rsid, voffset)`` pairs from a decompressed block."""
    for rsid, voffset in _ENTRY.iter_unpack(raw):
        yield rsid, voffset


@dataclass(frozen=True)
class BlockRecord:
    """A single signpost describing one compressed block on disk."""

    first_rsid: int
    comp_offset: int
    comp_len: int
    n_entries: int


class MetaWriter:
    """Accumulates :class:`BlockRecord` signposts and writes the ``.meta`` file."""

    def __init__(self, block_size: int = DEFAULT_BLOCK_SIZE):
        self.block_size = block_size
        self._blocks: list[BlockRecord] = []
        self._n_entries = 0

    def add_block(self, first_rsid: int, comp_offset: int, comp_len: int,
                  n_entries: int) -> None:
        self._blocks.append(
            BlockRecord(first_rsid, comp_offset, comp_len, n_entries))
        self._n_entries += n_entries

    def write(self, fh: BinaryIO) -> None:
        fh.write(_META_HEADER.pack(
            MAGIC, FORMAT_VERSION, self.block_size,
            len(self._blocks), self._n_entries))
        for b in self._blocks:
            fh.write(_META_RECORD.pack(
                b.first_rsid, b.comp_offset, b.comp_len, b.n_entries))


class MetaIndex:
    """Loads the ``.meta`` signpost table and binary-searches it by rsID."""

    def __init__(self, block_size: int, n_entries: int,
                 blocks: list[BlockRecord]):
        self.block_size = block_size
        self.n_entries = n_entries
        self.blocks = blocks
        self._first_rsids = [b.first_rsid for b in blocks]

    @classmethod
    def load(cls, path: str) -> "MetaIndex":
        with open(path, "rb") as fh:
            header = fh.read(_META_HEADER.size)
            if len(header) != _META_HEADER.size:
                raise ValueError(f"{path}: truncated meta header")
            magic, version, block_size, n_blocks, n_entries = \
                _META_HEADER.unpack(header)
            if magic != MAGIC:
                raise ValueError(f"{path}: bad magic {magic!r}")
            if version != FORMAT_VERSION:
                raise ValueError(f"{path}: unsupported format version {version}")
            blocks: list[BlockRecord] = []
            for _ in range(n_blocks):
                rec = fh.read(_META_RECORD.size)
                if len(rec) != _META_RECORD.size:
                    raise ValueError(f"{path}: truncated meta record")
                blocks.append(BlockRecord(*_META_RECORD.unpack(rec)))
        return cls(block_size, n_entries, blocks)

    def candidate_block(self, rsid: int) -> Optional[int]:
        """
        Index of the block whose range must contain the *first* occurrence of
        ``rsid`` (rightmost block with ``first_rsid <= rsid``); ``None`` when
        ``rsid`` sorts before every block.
        """
        if not self.blocks:
            return None
        i = bisect.bisect_right(self._first_rsids, rsid) - 1
        return i if i >= 0 else None


# ============================================================================
# Read-side API
# ============================================================================
@dataclass(frozen=True)
class VariantGroup:
    """One locus of one rsID: a shared REF with one or more ALT alleles."""

    rsid: int
    chrom: str
    pos: int
    ref: str
    alts: tuple[str, ...]

    def to_string(self) -> str:
        """``chr-pos-ref-alt,alt`` (the rsID -> sequence output format)."""
        return f"{self.chrom}-{self.pos}-{self.ref}-{','.join(self.alts)}"


class DbSNPFlat:
    """
    Bidirectional reader over one assembly's flat dbSNP store.

    ``normalizer`` must be the SAME callable used by the builder, or
    ``sequence -> rsID`` lookups will miss. It defaults to :func:`trim_alleles`
    (the builder's default too); pass ``None`` to disable, or a reference-aware
    normalizer for full left-alignment.
    """

    def __init__(self, master_gz: str, rsid_idx: str, rsid_meta: str,
                 normalizer=trim_alleles):
        self._meta = MetaIndex.load(rsid_meta)
        self._bgzf = pysam.BGZFile(master_gz, "rb")   # rsID -> master rows
        self._tabix = pysam.TabixFile(master_gz)      # sequence -> rsID
        self._idx_fh = open(rsid_idx, "rb")           # compressed index blocks
        self._normalizer = normalizer

    # ---- lifecycle -----------------------------------------------------
    def close(self) -> None:
        for h in (self._bgzf, self._tabix, self._idx_fh):
            try:
                h.close()
            except Exception:
                pass

    def __enter__(self) -> "DbSNPFlat":
        return self

    def __exit__(self, *exc) -> None:
        self.close()

    # ---- Direction 1:  rsID -> sequence --------------------------------
    def _read_block(self, block_idx: int) -> list[tuple[int, int]]:
        rec = self._meta.blocks[block_idx]
        self._idx_fh.seek(rec.comp_offset)
        raw = zlib.decompress(self._idx_fh.read(rec.comp_len))
        return list(iter_entries(raw))

    def _offsets_for_rsid(self, rsid: int) -> list[int]:
        """
        All master virtual offsets registered for ``rsid`` (usually one; more
        for a multi-mapping rsID). Duplicate entries are globally contiguous but
        can straddle block boundaries, so we scan the candidate block and walk
        outward in both directions while adjacent blocks can still hold it.
        """
        block_idx = self._meta.candidate_block(rsid)
        if block_idx is None:
            return []
        blocks = self._meta.blocks
        offsets: set[int] = set()

        offsets.update(v for r, v in self._read_block(block_idx) if r == rsid)

        nxt = block_idx + 1
        while nxt < len(blocks) and blocks[nxt].first_rsid == rsid:
            offsets.update(v for r, v in self._read_block(nxt) if r == rsid)
            nxt += 1

        k = block_idx
        while k > 0 and blocks[k].first_rsid == rsid:
            k -= 1
            offsets.update(v for r, v in self._read_block(k) if r == rsid)

        return sorted(offsets)

    def _read_group_at(self, voffset: int) -> Optional[VariantGroup]:
        """Read the contiguous (chrom, pos, rsid) group starting at ``voffset``."""
        self._bgzf.seek(voffset)
        chrom = pos = ref = rsid = None
        alts: list[str] = []
        while True:
            line = self._bgzf.readline()
            if not line:
                break
            f = [x.decode() for x in line.rstrip(b"\n").split(b"\t")]
            key = (f[C_CHROM], int(f[C_POS]), int(f[C_RSID]), f[C_REF])
            if chrom is None:
                chrom, pos, rsid, ref = key[0], key[1], key[2], key[3]
            elif key != (chrom, pos, rsid, ref):
                break
            alts.append(f[C_ALT])
        if chrom is None:
            return None
        return VariantGroup(rsid, chrom, pos, ref, tuple(alts))

    def rsid_to_variants(self, rsid) -> list[VariantGroup]:
        """Resolve an rsID to its variant group(s). Empty list if unknown."""
        if isinstance(rsid, str):
            rsid = int(rsid[2:] if rsid.lower().startswith("rs") else rsid)
        out = []
        for voffset in self._offsets_for_rsid(rsid):
            g = self._read_group_at(voffset)
            if g is not None:
                out.append(g)
        return out

    def resolve_rsid(self, rsid) -> list[str]:
        """rsID -> list of ``chr-pos-ref-alt,alt`` strings."""
        return [g.to_string() for g in self.rsid_to_variants(rsid)]

    # ---- Direction 2:  sequence -> rsID --------------------------------
    def variant_to_rsids(self, chrom: str, pos: int, ref: str,
                         alt: str) -> list[int]:
        """
        Map an exact ``chrom, pos, ref, alt`` to its rsID(s). The query is run
        through the same ``normalizer`` the store was built with, then matched
        exactly on REF and ALT.
        """
        if self._normalizer is not None:
            chrom, pos, ref, alt = self._normalizer(chrom, pos, ref, alt)
        pos = int(pos)
        try:
            rows = self._tabix.fetch(str(chrom), max(0, pos - 2), pos + 1)
        except (ValueError, KeyError):
            return []  # contig absent from this assembly
        found: list[int] = []
        for row in rows:
            f = row.split("\t")
            if int(f[C_POS]) == pos and f[C_REF] == ref and f[C_ALT] == alt:
                rsid = int(f[C_RSID])
                if rsid not in found:
                    found.append(rsid)
        return found

    def variant_to_rsid(self, chrom: str, pos: int, ref: str,
                        alt: str) -> Optional[int]:
        """Convenience: first matching rsID, or ``None``."""
        hits = self.variant_to_rsids(chrom, pos, ref, alt)
        return hits[0] if hits else None


# ============================================================================
# Builder
# ============================================================================
# Master TSV columns (0-based): CHROM, POS, END, rsID, REF, ALT
MASTER_SEQ_COL = 0
MASTER_START_COL = 1
MASTER_END_COL = 2

# A parsed variant row: one ALT allele.
Variant = tuple[str, int, int, str, str]  # (chrom, pos, rsid, ref, alt)


def _strip_rs(raw) -> Optional[int]:
    """'rs6025' / '6025' -> 6025 ; unusable ids -> None."""
    if raw is None:
        return None
    raw = str(raw).strip()
    if raw in ("", "."):
        return None
    if raw.lower().startswith("rs"):
        raw = raw[2:]
    raw = raw.split(";")[0].split(",")[0]  # dbSNP ids are occasionally joined
    try:
        return int(raw)
    except ValueError:
        return None


def _detect_source_type(path: str) -> str:
    low = path.lower()
    if low.endswith(".bcf"):
        return "bcf"
    if low.endswith((".tsv", ".txt")):
        return "tsv"
    return "vcf"  # .vcf, .vcf.gz, .vcf.bgz


class FlatDbSNPBuilder:
    def __init__(self, out_dir: str = "dbsnp_flat",
                 block_size: int = DEFAULT_BLOCK_SIZE,
                 compress_level: int = 6,
                 contig_map: Optional[dict] = None,
                 normalizer=trim_alleles):
        self.out_dir = os.path.abspath(out_dir)
        self.block_size = block_size
        # zlib level 0-9 (stdlib; no third-party compression dependency).
        self.compress_level = compress_level
        # Explicit accession->name override; when None, build() uses the map for
        # the assembly being built (get_contig_map). Empty dict disables mapping.
        self._contig_map_override = contig_map
        self.contig_map: dict = contig_map or {}
        # Applied to every stored variant; MUST match DbSNPFlat's normalizer.
        self.normalizer = normalizer
        os.makedirs(self.out_dir, exist_ok=True)

    # ------------------------------------------------------------------
    # Helpers preserved from the original builder
    # ------------------------------------------------------------------
    @staticmethod
    def _extract_version_tag(path: str) -> str:
        """vvdbSNP_GRCh37_b157_20241205.bcf -> b157_20241205"""
        name = os.path.basename(path).rsplit(".", 1)[0]
        parts = name.split("_")
        return "_".join(parts[-2:]) if len(parts) >= 2 else name

    # ------------------------------------------------------------------
    # Source parsing -> stream of (chrom, pos, rsid, ref, alt)
    # ------------------------------------------------------------------
    def _iter_bcf(self, path: str) -> Iterator[Variant]:
        """
        Decode a BCF/VCF with ``bcftools query`` (fast, streaming) and split
        multi-allelic records in Python. Matches the original pipeline, which
        already produces BCFs with bcftools.
        """
        cmd = ["bcftools", "query",
               "-f", "%CHROM\t%POS\t%REF\t%ALT\t%ID\n", path]
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True)
        assert proc.stdout is not None
        for line in proc.stdout:
            chrom, pos, ref, alt_field, id_field = \
                line.rstrip("\n").split("\t")
            rsid = _strip_rs(id_field)
            if rsid is None:
                continue
            chrom = self.contig_map.get(chrom, chrom)
            for alt in alt_field.split(","):
                if alt and not alt.startswith("<"):
                    yield chrom, int(pos), rsid, ref, alt
        proc.stdout.close()
        if proc.wait() != 0:
            raise subprocess.CalledProcessError(proc.returncode, cmd)

    def _iter_vcf(self, path: str) -> Iterator[Variant]:
        """Read a VCF/BCF directly with pysam (no bcftools binary needed)."""
        with pysam.VariantFile(path) as vf:
            for rec in vf:
                rsid = _strip_rs(rec.id)
                if rsid is None or rec.ref is None:
                    continue
                chrom = self.contig_map.get(rec.chrom, rec.chrom)
                for alt in (rec.alts or ()):
                    if alt and not alt.startswith("<"):
                        yield chrom, rec.pos, rsid, rec.ref, alt

    def _iter_tsv(self, path: str) -> Iterator[Variant]:
        """Minimal TSV: ``chrom  pos  rsid  ref  alt`` (alt may be comma list)."""
        with open(path) as fh:
            for line in fh:
                if not line.strip() or line.startswith("#"):
                    continue
                chrom, pos, rsid_raw, ref, alt_field = \
                    line.rstrip("\n").split("\t")[:5]
                rsid = _strip_rs(rsid_raw)
                if rsid is None:
                    continue
                chrom = self.contig_map.get(chrom, chrom)
                for alt in alt_field.split(","):
                    alt = alt.strip()
                    if alt and not alt.startswith("<"):
                        yield chrom, int(pos), rsid, ref, alt

    def _iter_source(self, source, source_type: str) -> Iterator[Variant]:
        if source_type == "bcf":
            stream = self._iter_bcf(source)
        elif source_type == "vcf":
            stream = self._iter_vcf(source)
        elif source_type == "tsv":
            stream = self._iter_tsv(source)
        else:
            raise ValueError(f"unknown source_type {source_type!r}")

        if self.normalizer is None:
            yield from stream
            return
        for chrom, pos, rsid, ref, alt in stream:
            chrom, pos, ref, alt = self.normalizer(chrom, pos, ref, alt)
            yield chrom, int(pos), rsid, ref, alt

    # ------------------------------------------------------------------
    # Master file
    # ------------------------------------------------------------------
    def _write_master(self, variants: Iterable[Variant], master_tsv: str) -> None:
        """
        Write the sorted master TSV (uncompressed), sorted by (chrom, pos, rsid)
        so all ALT rows of one rsID at one locus are contiguous - needed for the
        index build and for multi-allelic re-assembly at query time.
        """
        unsorted = master_tsv + ".unsorted"
        with open(unsorted, "w") as out:
            for chrom, pos, rsid, ref, alt in variants:
                end = pos + len(ref) - 1
                out.write(f"{chrom}\t{pos}\t{end}\t{rsid}\t{ref}\t{alt}\n")

        env = dict(os.environ, LC_ALL="C")
        with open(master_tsv, "w") as out:
            subprocess.run(
                ["sort", "-t", "\t", "-k1,1", "-k2,2n", "-k4,4n",
                 "-k5,5", "-k6,6", unsorted],
                check=True, stdout=out, env=env)
        os.unlink(unsorted)

    def _index_master(self, master_tsv: str) -> str:
        """bgzip + tabix the master TSV via pysam. Returns the .gz path."""
        return pysam.tabix_index(
            master_tsv,
            seq_col=MASTER_SEQ_COL, start_col=MASTER_START_COL,
            end_col=MASTER_END_COL, zerobased=False,
            force=True, keep_original=False)

    # ------------------------------------------------------------------
    # rsID index
    # ------------------------------------------------------------------
    def _collect_pairs(self, master_gz: str, pairs_path: str) -> None:
        """
        Walk the bgzipped master and emit one ``rsid<TAB>voffset`` pair per
        distinct (chrom, pos, rsid) group, where ``voffset`` is the BGZF virtual
        offset of that group's first row.
        """
        with pysam.BGZFile(master_gz, "rb") as bgz, \
                open(pairs_path, "w") as pairs:
            last_key = None
            while True:
                vfo = bgz.tell()
                line = bgz.readline()
                if not line:
                    break
                f = line.rstrip(b"\n").split(b"\t")
                key = (f[0], f[1], f[3], f[4])  # chrom, pos, rsid, ref
                if key != last_key:
                    pairs.write(f"{int(f[3])}\t{vfo}\n")
                    last_key = key

    def _sort_pairs(self, pairs_path: str) -> str:
        """Sort pairs by rsID (numeric). THIS is what makes the index searchable."""
        sorted_path = pairs_path + ".sorted"
        env = dict(os.environ, LC_ALL="C")
        with open(sorted_path, "w") as out:
            subprocess.run(
                ["sort", "-t", "\t", "-k1,1n", "-k2,2n", pairs_path],
                check=True, stdout=out, env=env)
        return sorted_path

    def _write_index(self, sorted_pairs: str, idx_path: str,
                     meta_path: str) -> None:
        meta = MetaWriter(block_size=self.block_size)

        def flush(buf, comp_offset):
            raw = b"".join(pack_entry(r, v) for r, v in buf)
            comp = zlib.compress(raw, self.compress_level)
            idx_out.write(comp)
            meta.add_block(buf[0][0], comp_offset, len(comp), len(buf))
            return comp_offset + len(comp)

        with open(sorted_pairs) as pairs, open(idx_path, "wb") as idx_out:
            buf, comp_offset = [], 0
            for line in pairs:
                rsid_s, voff_s = line.split("\t")
                buf.append((int(rsid_s), int(voff_s)))
                if len(buf) == self.block_size:
                    comp_offset = flush(buf, comp_offset)
                    buf = []
            if buf:
                flush(buf, comp_offset)

        with open(meta_path, "wb") as meta_out:
            meta.write(meta_out)

    # ------------------------------------------------------------------
    # Orchestration
    # ------------------------------------------------------------------
    def build(self, assembly: str, source: str,
              source_type: Optional[str] = None,
              version: Optional[str] = None) -> dict:
        """
        Build all artefacts for one assembly. ``source`` is a BCF/VCF/TSV path;
        ``source_type`` is auto-detected from the extension when omitted, and
        ``version`` defaults to the tag parsed from the filename.

        Backwards compatible with ``populate.py``'s ``builder.build(asm, bcf)``.
        """
        if source_type is None:
            source_type = _detect_source_type(source)
        if version is None:
            version = self._extract_version_tag(source)

        # Resolve the accession->name map for this assembly (unless an explicit
        # override was supplied at construction).
        if self._contig_map_override is None:
            self.contig_map = get_contig_map(assembly)
        else:
            self.contig_map = self._contig_map_override

        version_dir = os.path.join(self.out_dir, version)
        os.makedirs(version_dir, exist_ok=True)
        stem = os.path.join(version_dir, f"dbsnp_{assembly}_{version}")

        master_tsv = stem + ".variants.tsv"
        master_gz = master_tsv + ".gz"
        tbi = master_gz + ".tbi"
        idx_path = stem + ".rsid.idx"
        meta_path = stem + ".rsid.meta"

        outputs = {"master": master_gz, "tabix": tbi,
                   "rsid_idx": idx_path, "rsid_meta": meta_path}

        # Skip if a complete build already exists (idempotent, like the original).
        if all(os.path.exists(p) for p in outputs.values()):
            print(f"✔ Flat dbSNP already exists for {assembly} ({version}),"
                  f" skipping")
            return outputs

        print(f"Building flat dbSNP for {assembly} (version {version}) "
              f"from {source_type}: {source}")

        # 1. master
        self._write_master(
            self._iter_source(source, source_type), master_tsv)
        self._index_master(master_tsv)  # -> master_gz + tbi

        # 2. rsID index (sorted by rsID)
        with tempfile.TemporaryDirectory(dir=version_dir) as tmp:
            pairs = os.path.join(tmp, "pairs.tsv")
            self._collect_pairs(master_gz, pairs)
            sorted_pairs = self._sort_pairs(pairs)
            self._write_index(sorted_pairs, idx_path, meta_path)

        print(f"✔ Flat dbSNP build complete for {assembly}")
        return outputs


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def _main(argv=None) -> int:
    import argparse
    import sys

    parser = argparse.ArgumentParser(
        description="Build flat dbSNP files (master + rsID index) from "
                    "dbSNP BCF/VCF (or TSV).")
    parser.add_argument("--out", default="dbsnp_flat",
                        help="Output directory (default: dbsnp_flat)")
    parser.add_argument("--version", default=None,
                        help="Version tag (default: parsed from filename)")
    parser.add_argument("--source-type", choices=["bcf", "vcf", "tsv"],
                        default=None, help="Input type (default: by extension)")
    parser.add_argument("--block-size", type=int, default=DEFAULT_BLOCK_SIZE)
    parser.add_argument("--compress-level", type=int, default=6,
                        help="zlib compression level 0-9 (default: 6)")
    parser.add_argument("--no-normalise", action="store_true",
                        help="Disable allele trimming of stored variants")
    parser.add_argument("--no-contig-map", action="store_true",
                        help="Do not rename RefSeq accessions to chromosome names")
    parser.add_argument("--grch37", help="Path to GRCh37 dbSNP BCF/VCF/TSV")
    parser.add_argument("--grch38", help="Path to GRCh38 dbSNP BCF/VCF/TSV")
    args = parser.parse_args(argv)

    if not args.grch37 and not args.grch38:
        parser.error("Specify at least one of --grch37 / --grch38")

    builder = FlatDbSNPBuilder(
        out_dir=args.out, block_size=args.block_size,
        compress_level=args.compress_level,
        contig_map={} if args.no_contig_map else None,
        normalizer=None if args.no_normalise else trim_alleles)

    for assembly, src in (("grch37", args.grch37), ("grch38", args.grch38)):
        if not src:
            continue
        if not os.path.exists(src):
            sys.exit(f"ERROR: {assembly} source not found: {src}")
        paths = builder.build(assembly, src, source_type=args.source_type,
                              version=args.version)
        for key, path in paths.items():
            print(f"  {key:10} {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(_main())

