"""Tests for vvDBFlatValidate.py (validation logic + reporting + CLI)."""

from __future__ import annotations

import json
import math
import os
import random

import pytest

import vvDBFlatValidate as V
import vvDBFlatQueryDbSNP as Q


# ---------------------------------------------------------------------------
# pure helpers
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("ref,alt,expected", [
    ("C", "T", "snv"),
    ("A", "ATG", "insertion"),
    ("AGT", "A", "deletion"),
    ("AGT", "CGA", "mnv"),
])
def test_classify(ref, alt, expected):
    assert V._classify(ref, alt) == expected


def test_pctl():
    vals = [float(x) for x in range(1, 11)]
    assert V._pctl(vals, 0) == 1.0
    assert V._pctl(vals, 100) == 10.0
    assert V._pctl([5.0], 50) == 5.0
    assert math.isnan(V._pctl([], 50))


def test_summ_ms():
    s = V._summ_ms([0.001, 0.002, 0.003])   # seconds in -> ms out
    assert s["n"] == 3
    assert s["min_ms"] == 1.0
    assert s["max_ms"] == 3.0
    assert s["mean_ms"] == 2.0
    assert V._summ_ms([]) == {"n": 0}


def test_rate():
    assert V._rate(3, 4) == 75.0
    assert V._rate(4, 4) == 100.0
    assert math.isnan(V._rate(0, 0))


# ---------------------------------------------------------------------------
# truth-set sampling
# ---------------------------------------------------------------------------
def test_sample_source_basic(synth):
    recs = V.sample_source(synth["source"], "grch38", per_type=5,
                           rng=random.Random(0))
    assert recs
    assert all(r.chrom in ("1", "11", "19") for r in recs)
    assert all(isinstance(r.rsid, int) for r in recs)
    # at least the common SNV bucket should have filled
    types = {V._classify(r.ref, a) for r in recs for a in r.alts}
    assert "snv" in types


def test_sample_source_reproducible(synth):
    a = V.sample_source(synth["source"], "grch38", 5, random.Random(0))
    b = V.sample_source(synth["source"], "grch38", 5, random.Random(0))
    assert [(r.rsid, r.pos) for r in a] == [(r.rsid, r.pos) for r in b]


def test_sample_source_requires_index(tmp_path):
    with pytest.raises(SystemExit):
        V.sample_source(str(tmp_path / "no_such.vcf.gz"), "grch38", 5,
                        random.Random(0))


def test_sample_store_rsids(synth):
    ids = V.sample_store_rsids(synth["stem"], 20, random.Random(0))
    assert ids and all(isinstance(x, int) for x in ids)
    with Q.DbSNPFlatQuery({"grch38": synth["stem"]}) as c:
        assert c.rsid_to_variants("grch38", ids[0]).found


# ---------------------------------------------------------------------------
# accuracy check (through the query API)
# ---------------------------------------------------------------------------
def _landmark_records():
    return [
        V.TruthRecord(6025, "1", 169549811, "C", ("T",), False),
        V.TruthRecord(334, "11", 5227002, "T", ("A", "C", "G"), True),
        V.TruthRecord(99, "1", 250, "AA", ("GA",), False),  # non-minimal source
    ]


def test_check_accuracy_pass(synth):
    with Q.DbSNPFlatQuery({"grch38": synth["stem"]}) as c:
        res = V.check_accuracy(c, "grch38", _landmark_records(), normalize=True)
    assert res["passed"]
    assert res["forward"]["pass_rate"] == 100.0
    assert res["reverse"]["pass_rate"] == 100.0
    assert res["round_trip"]["pass_rate"] == 100.0


def test_check_accuracy_detects_wrong_rsid(synth):
    # A truth record claiming the wrong rsID at rs6025's locus must FAIL.
    bad = [V.TruthRecord(424242, "1", 169549811, "C", ("T",), False)]
    with Q.DbSNPFlatQuery({"grch38": synth["stem"]}) as c:
        res = V.check_accuracy(c, "grch38", bad, normalize=True)
    assert not res["passed"]
    assert res["forward"]["pass_rate"] < 100.0


def test_check_accuracy_no_normalise_misses_nonminimal(synth):
    # With normalization off, the non-minimal rs99 (AA/GA) can't match the
    # trimmed store entry -> the check must notice.
    rec = [V.TruthRecord(99, "1", 250, "AA", ("GA",), False)]
    with Q.DbSNPFlatQuery({"grch38": synth["stem"]}, normalizer=None) as c:
        res = V.check_accuracy(c, "grch38", rec, normalize=False)
    assert not res["passed"]


# ---------------------------------------------------------------------------
# format + curated checks
# ---------------------------------------------------------------------------
def test_check_formats_pass(synth):
    recs = [
        V.TruthRecord(6025, "1", 169549811, "C", ("T",), False),
        V.TruthRecord(334, "11", 5227002, "T", ("A", "C", "G"), True),
    ]
    with Q.DbSNPFlatQuery({"grch38": synth["stem"]}) as c:
        res = V.check_formats(c, "grch38", recs, normalize=True)
    assert res["passed"]
    assert res["variant_formats"]["canonical"]["pass_rate"] == 100.0
    assert res["variant_formats"]["lowercase"]["pass_rate"] == 100.0
    assert res["variant_formats"]["nonminimal_suffix"]["pass_rate"] == 100.0
    assert res["rsid_spelling"]["pass_rate"] == 100.0


def test_check_curated_pass(synth):
    # the store contains rs6025 (chr1), rs334 (chr11), rs429358 (chr19)
    with Q.DbSNPFlatQuery({"grch38": synth["stem"]}) as c:
        res = V.check_curated(c, "grch38")
    assert res["passed"]
    assert {row["rsid"] for row in res["landmarks"]} == {6025, 334, 429358}
    assert all(row["ok"] for row in res["landmarks"])


# ---------------------------------------------------------------------------
# benchmark
# ---------------------------------------------------------------------------
def test_benchmark(synth):
    stem = synth["stem"]
    client = Q.DbSNPFlatQuery({"grch38": stem})
    try:
        out = V.benchmark(lambda: Q.DbSNPFlatQuery({"grch38": stem}), client,
                          "grch38", [("1", 169549811, "C", "T")], [6025, 334],
                          n_queries=10, rng=random.Random(0))
    finally:
        client.close()
    assert out["variant_to_rsid"]["n"] == 10
    assert out["rsid_to_variant"]["n"] == 10
    assert out["variant_to_rsid"]["median_ms"] >= 0
    assert out["rsid_to_variant"]["cold_first_query_ms"] is not None


# ---------------------------------------------------------------------------
# full per-assembly driver
# ---------------------------------------------------------------------------
def test_validate_assembly_pass(synth):
    res = V.validate_assembly("grch38", synth["stem"], synth["source"],
                              per_type=5, bench_queries=10, normalize=True, seed=1)
    assert res["passed"] is True
    assert res["accuracy"]["forward"]["pass_rate"] == 100.0
    assert res["accuracy"]["reverse"]["pass_rate"] == 100.0
    assert res["accuracy"]["round_trip"]["pass_rate"] == 100.0
    assert res["formats"]["passed"]
    assert res["curated"]["passed"]
    assert "variant_to_rsid" in res["timing"]


def test_validate_assembly_no_source(synth):
    res = V.validate_assembly("grch38", synth["stem"], None,
                              per_type=5, bench_queries=10, normalize=True, seed=1)
    assert res["accuracy"] is None
    assert res["formats"] is None
    assert res["store_round_trip"]["passed"]
    assert res["passed"] is True


# ---------------------------------------------------------------------------
# reporting
# ---------------------------------------------------------------------------
def test_write_reports(synth, tmp_path):
    res = V.validate_assembly("grch38", synth["stem"], None,
                              per_type=3, bench_queries=5, normalize=True, seed=1)
    json_path, md_path = V.write_reports([res], str(tmp_path / "rep"))
    assert os.path.exists(json_path) and os.path.exists(md_path)
    with open(json_path) as fh:
        payload = json.load(fh)
    assert payload["results"][0]["assembly"] == "grch38"
    assert "# Flat dbSNP validation report" in open(md_path).read()


# ---------------------------------------------------------------------------
# CLI (_main)
# ---------------------------------------------------------------------------
def test_main_single_assembly_pass(synth):
    rc = V._main(["--assembly", "grch38", "--store", synth["stem"],
                  "--source", synth["source"], "--sample-per-type", "5",
                  "--bench-queries", "10", "--seed", "1", "--no-report-files"])
    assert rc == 0


def test_main_writes_reports(synth, tmp_path):
    rc = V._main(["--assembly", "grch38", "--store", synth["stem"],
                  "--sample-per-type", "3", "--bench-queries", "5",
                  "--report-dir", str(tmp_path / "out")])
    assert rc == 0
    files = os.listdir(tmp_path / "out")
    assert any(f.endswith(".json") for f in files)
    assert any(f.endswith(".md") for f in files)


@pytest.mark.parametrize("argv", [
    [],                                                       # nothing
    ["--assembly", "grch38"],                                 # --assembly, no --store
    ["--store", "x"],                                         # --store, no --assembly
    ["--assembly", "grch38", "--store", "x", "--grch37-store", "y"],  # mixed interfaces
])
def test_main_argument_guards(argv):
    with pytest.raises(SystemExit):
        V._main(argv)
