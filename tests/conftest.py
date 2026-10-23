"""
Shared pytest fixtures for the flat-dbSNP test suite.

The centrepiece is the session-scoped ``synth`` fixture: it writes a small,
NCBI-style source VCF (RefSeq-accession contigs, *no* ``##contig`` header lines,
bgzipped + tabix-indexed) that contains a set of KNOWN landmark variants plus
enough filler spread across three chromosomes for the validator's random
sampling to work, then builds a flat store from it with the real builder. All
three modules are exercised against this one store.
"""

from __future__ import annotations

import os
import random
import sys

import pytest
import pysam

# Make the three modules under test importable (they live one dir up).
_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

import vvDBFlatBuildDbSNP as build  # noqa: E402

# GRCh38 accession -> (plain name, length). Lengths match CONTIG_LENGTHS in the
# validator, so sample_source can pick random positions on these contigs.
_ACC = {
    "NC_000001.11": ("1", 248956422),
    "NC_000011.10": ("11", 135086622),
    "NC_000019.10": ("19", 58617616),
}

# Landmark records as they appear in the SOURCE VCF (rsID, raw REF/ALT).
_LANDMARKS_SOURCE = [
    ("NC_000001.11", 250, "rs99", "AA", "GA"),        # non-minimal SNV -> A/G
    ("NC_000001.11", 300000, "rs201", "A", "ATG"),    # insertion
    ("NC_000001.11", 400000, "rs202", "AGT", "A"),    # deletion
    ("NC_000001.11", 500000, "rs203", "AGT", "CGA"),  # MNV
    ("NC_000001.11", 169549811, "rs6025", "C", "T"),  # F5 Leiden (chr1)
    ("NC_000011.10", 5227002, "rs334", "T", "A,C,G"), # sickle-cell, multi-allelic
    ("NC_000019.10", 44908684, "rs429358", "T", "C"), # APOE (chr19)
]

# Expected STORED (normalized) form of each landmark, for assertions:
#   rsid -> (chrom, pos, ref, alts)
LANDMARKS = {
    99: ("1", 250, "A", ("G",)),
    201: ("1", 300000, "A", ("ATG",)),
    202: ("1", 400000, "AGT", ("A",)),
    203: ("1", 500000, "AGT", ("CGA",)),
    6025: ("1", 169549811, "C", ("T",)),
    334: ("11", 5227002, "T", ("A", "C", "G")),
    429358: ("19", 44908684, "T", ("C",)),
}

_LANDMARK_POSITIONS = {250, 300000, 400000, 500000, 169549811, 5227002, 44908684}


def _make_source(path: str) -> None:
    """Write a bgzipped + tabix-indexed synthetic dbSNP VCF at ``path``."""
    rng = random.Random(1234)
    bases = "ACGT"
    per_contig: dict[str, list] = {acc: [] for acc in _ACC}

    for acc, pos, rid, ref, alt in _LANDMARKS_SOURCE:
        per_contig[acc].append((pos, rid, ref, alt))

    rsid = 1_000_000
    step = 60_000
    for acc, (_name, length) in _ACC.items():
        pos = 15_000
        i = 0
        while pos < length - 10:
            rsid += 1
            i += 1
            k = i % 8
            if k in (0, 1, 2):                                  # SNV
                r = rng.choice(bases)
                a = rng.choice([b for b in bases if b != r])
                ref, alt = r, a
            elif k == 3:                                        # insertion
                r = rng.choice(bases)
                ref, alt = r, r + rng.choice(bases) + rng.choice(bases)
            elif k == 4:                                        # deletion
                r = rng.choice(bases)
                ref, alt = r + rng.choice(bases), r
            elif k == 5:                                        # MNV
                ref = "".join(rng.choice(bases) for _ in range(3))
                alt = "".join(rng.choice(bases) for _ in range(3))
                if alt == ref:
                    alt = ("A" if ref[0] != "A" else "C") + ref[1:]
            elif k == 6:                                        # multi-allelic SNV
                r = rng.choice(bases)
                a1, a2 = [b for b in bases if b != r][:2]
                ref, alt = r, f"{a1},{a2}"
            else:                                               # multi-allelic mixed
                r = rng.choice(bases)
                a1 = rng.choice([b for b in bases if b != r])
                ref, alt = r, f"{a1},{r + rng.choice(bases)}"
            if pos not in _LANDMARK_POSITIONS:
                per_contig[acc].append((pos, f"rs{rsid}", ref, alt))
            pos += step + rng.randint(-3000, 3000)

    plain = path[:-3] if path.endswith(".gz") else path + ".vcf"
    with open(plain, "w") as fh:
        fh.write("##fileformat=VCFv4.2\n")
        # NOTE: deliberately no ##contig lines, mimicking real NCBI dbSNP VCFs.
        fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for acc in _ACC:                       # contigs grouped, records sorted
            for pos, rid, ref, alt in sorted(per_contig[acc]):
                fh.write(f"{acc}\t{pos}\t{rid}\t{ref}\t{alt}\t.\t.\t.\n")
    pysam.tabix_compress(plain, path, force=True)
    pysam.tabix_index(path, preset="vcf", force=True)


@pytest.fixture(scope="session")
def synth(tmp_path_factory) -> dict:
    """
    Build the synthetic source VCF + flat store once for the whole session.

    Returns a dict with: ``source`` (vcf.gz), ``out_dir``, ``version``,
    ``assembly``, ``stem`` (store path stem), ``outputs`` (build() return),
    and ``landmarks`` (rsid -> expected stored (chrom, pos, ref, alts)).
    """
    d = tmp_path_factory.mktemp("dbsnp_flat_suite")
    source = str(d / "synth_grch38.vcf.gz")
    _make_source(source)

    out_dir = str(d / "store")
    builder = build.FlatDbSNPBuilder(out_dir=out_dir)
    outputs = builder.build("grch38", source, version="test")
    stem = os.path.join(out_dir, "test", "dbsnp_grch38_test")

    return {
        "source": source,
        "out_dir": out_dir,
        "version": "test",
        "assembly": "grch38",
        "stem": stem,
        "outputs": outputs,
        "landmarks": LANDMARKS,
    }
