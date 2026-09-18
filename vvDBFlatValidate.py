#!/usr/bin/env python3
"""
vvDBFlatValidate.py

Validate a built flat dbSNP store (see ``vvDBFlatBuildDbSNP.py``) for
**accuracy** and **query speed**, and write a human-readable report so the
results are reproducible and documentable by anyone.

What it checks
--------------
1. ACCURACY against an independent truth set - the *source* dbSNP VCF the store
   was built from (e.g. NCBI ``GCF_000001405.25.gz`` for GRCh37,
   ``GCF_000001405.40.gz`` for GRCh38). Real records are randomly sampled
   straight from that VCF, **stratified by variant type** (SNV, insertion,
   deletion, MNV) and tracking multi-allelic sites, then each is checked in
   BOTH directions:

     * sequence -> rsID   ``DbSNPFlatQuery.variant_to_rsids``  (must contain rsID)
     * rsID -> sequence   ``DbSNPFlatQuery.rsid_to_variants``  (must contain locus)

   The store is derived from dbSNP, so testing it against itself would prove
   nothing; re-parsing the *source* VCF independently is what makes this a real
   accuracy test.

   Every check runs through the **backend query API** (``DbSNPFlatQuery`` from
   ``vvDBFlatQueryDbSNP.py``) - the same functions a VariantValidator
   integration will call - so a pass certifies that querying itself works, not a
   parallel reimplementation of it.

2. FORMAT robustness - the same lookups issued in several input spellings that
   must all resolve to the same answer:

     * rsID:    ``rs6025`` / ``6025`` / ``6025`` (int)
     * variant: canonical, lower-cased alleles, and non-minimal (prefix/suffix
                padded) spellings that the normalizer must collapse.

3. ROUND-TRIP consistency - rsID -> variant -> rsID returns the original rsID.
   (Self-consistent; catches index/offset corruption. Note the normalizer is a
   *shared* seam used at build and query time, so a symmetric bug inside
   ``trim_alleles`` is invisible to any self-consistent test - the curated
   landmarks below, printed with their coordinates, are the manual cross-check
   against dbSNP for that.)

4. CURATED landmarks - a tiny, file-free smoke set of well-known variants,
   checked for the right chromosome + round-trip and printed for eyeballing.

5. TIMING - per-query latency (min/mean/median/p95/p99/max) and throughput for
   both directions, plus a cold first-query measurement.

Usage
-----
Validate ONE assembly at a time (run these as two separate commands; each
writes its own report)::

    python vvDBFlatValidate.py --assembly grch37 \
        --store  dbsnp_flat/<ver>/dbsnp_grch37_<ver> \
        --source /data/GCF_000001405.25.gz

    python vvDBFlatValidate.py --assembly grch38 \
        --store  dbsnp_flat/<ver>/dbsnp_grch38_<ver> \
        --source /data/GCF_000001405.40.gz \
        --sample-per-type 100 --bench-queries 1000 --seed 7

Or validate both in a single run with the paired flags::

    python vvDBFlatValidate.py \
        --grch37-store  dbsnp_flat/<ver>/dbsnp_grch37_<ver> \
        --grch37-source /data/GCF_000001405.25.gz \
        --grch38-store  dbsnp_flat/<ver>/dbsnp_grch38_<ver> \
        --grch38-source /data/GCF_000001405.40.gz \
        --report-dir validation_reports

Without a source VCF (no independent accuracy test; still does curated +
round-trip + timing using rsIDs harvested from the store itself)::

    python vvDBFlatValidate.py --assembly grch38 \
        --store dbsnp_flat/<ver>/dbsnp_grch38_<ver>

The store *stem* is the path without the ``.variants.tsv.gz`` / ``.rsid.idx`` /
``.rsid.meta`` suffixes (the same value ``vvDBFlatQuery.py`` takes).

Exit status is non-zero if any accuracy / format / round-trip check fails.
"""

from __future__ import annotations

import argparse
import json
import os
import random
import statistics
import sys
import time
import zlib
from dataclasses import dataclass, field, asdict
from datetime import datetime, timezone
from typing import Callable, Iterable, Optional

import pysam

from vvDBFlatBuildDbSNP import (
    MetaIndex,
    get_contig_map,
    iter_entries,
    trim_alleles,
)
# Validation drives the SAME backend query API a VariantValidator integration
# would call, so a pass genuinely certifies that querying works - not a parallel
# reimplementation of it.
from vvDBFlatQueryDbSNP import DbSNPFlatQuery

# Variant-type buckets we stratify the truth sample into.
VTYPES = ("snv", "insertion", "deletion", "mnv")

# Primary-assembly-unit lengths (bp), by plain chromosome name. NCBI dbSNP VCFs
# do NOT define contigs (or their lengths) in the header, so the validator
# carries these itself to drive weighted random-position sampling. These are the
# standard GRCh37 / GRCh38 reference lengths and are stable.
CONTIG_LENGTHS = {
    "grch38": {
        "1": 248956422, "2": 242193529, "3": 198295559, "4": 190214555,
        "5": 181538259, "6": 170805979, "7": 159345973, "8": 145138636,
        "9": 138394717, "10": 133797422, "11": 135086622, "12": 133275309,
        "13": 114364328, "14": 107043718, "15": 101991189, "16": 90338345,
        "17": 83257441, "18": 80373285, "19": 58617616, "20": 64444167,
        "21": 46709983, "22": 50818468, "X": 156040895, "Y": 57227415,
        "MT": 16569,
    },
    "grch37": {
        "1": 249250621, "2": 243199373, "3": 198022430, "4": 191154276,
        "5": 180915260, "6": 171115067, "7": 159138663, "8": 146364022,
        "9": 141213431, "10": 135534747, "11": 135006516, "12": 133851895,
        "13": 115169878, "14": 107349540, "15": 102531392, "16": 90354753,
        "17": 81195210, "18": 78077248, "19": 59128983, "20": 63025520,
        "21": 48129895, "22": 51304566, "X": 155270560, "Y": 59373566,
        "MT": 16569,
    },
}


# ---------------------------------------------------------------------------
# small helpers
# ---------------------------------------------------------------------------
def _classify(ref: str, alt: str) -> str:
    """Classify one REF/ALT allele pair by length (pre-normalization)."""
    lr, la = len(ref), len(alt)
    if lr == 1 and la == 1:
        return "snv"
    if la > lr:
        return "insertion"
    if lr > la:
        return "deletion"
    return "mnv"  # equal length, > 1 (MNV / block substitution)


def _strip_rs(raw) -> Optional[int]:
    """'rs6025' / '6025' -> 6025 ; unusable ids -> None.  (Mirrors the builder.)"""
    if raw is None:
        return None
    raw = str(raw).strip()
    if raw in ("", "."):
        return None
    if raw.lower().startswith("rs"):
        raw = raw[2:]
    raw = raw.split(";")[0].split(",")[0]
    try:
        return int(raw)
    except ValueError:
        return None


def _pctl(sorted_vals: list[float], p: float) -> float:
    """Nearest-rank percentile of an already-sorted list (p in 0..100)."""
    if not sorted_vals:
        return float("nan")
    if len(sorted_vals) == 1:
        return sorted_vals[0]
    k = int(round((p / 100.0) * (len(sorted_vals) - 1)))
    return sorted_vals[min(max(k, 0), len(sorted_vals) - 1)]


def _summ_ms(latencies_s: list[float]) -> dict:
    """Summarise a list of per-query latencies (seconds) into a ms stats dict."""
    ms = sorted(t * 1000.0 for t in latencies_s)
    n = len(ms)
    if n == 0:
        return {"n": 0}
    total_s = sum(latencies_s)
    return {
        "n": n,
        "min_ms": round(ms[0], 4),
        "mean_ms": round(statistics.fmean(ms), 4),
        "median_ms": round(_pctl(ms, 50), 4),
        "p95_ms": round(_pctl(ms, 95), 4),
        "p99_ms": round(_pctl(ms, 99), 4),
        "max_ms": round(ms[-1], 4),
        "throughput_qps": round(n / total_s, 1) if total_s > 0 else float("inf"),
    }


# ---------------------------------------------------------------------------
# sampled truth record
# ---------------------------------------------------------------------------
@dataclass
class TruthRecord:
    """One dbSNP site sampled from the source VCF (already contig-mapped)."""

    rsid: int
    chrom: str          # plain name (1..22, X, Y, MT)
    pos: int
    ref: str
    alts: tuple[str, ...]
    multiallelic: bool


def sample_source(source: str, assembly: str, per_type: int,
                  rng: random.Random, window: int = 200_000,
                  max_attempts: int = 200_000) -> list[TruthRecord]:
    """
    Randomly sample real records from the source dbSNP VCF, stratified so each
    variant type in :data:`VTYPES` reaches ``per_type`` allele-cases (MNVs are
    rare and may fall short - that is reported, not fatal).

    Sampling uses tabix random-region jumps (fast; no full scan of the ~30 GB
    file). Only primary assembly units (the contigs in the assembly's contig
    map) are sampled, so every locus has a plain chromosome name.

    NCBI dbSNP VCFs do NOT define contigs (or their lengths) in the header, so
    the contig list is taken from the tabix index and the chromosome lengths
    come from :data:`CONTIG_LENGTHS`. Records are read as raw VCF text (via
    ``TabixFile``), which also avoids the noisy INFO/RS parse warnings htslib
    emits on these files. A ``.tbi``/``.csi`` index alongside the .gz is
    required for the random-access fetch.
    """
    accmap = get_contig_map(assembly)          # accession -> plain name
    if not accmap:
        raise SystemExit(f"unknown assembly {assembly!r} (no contig map)")
    lengths = CONTIG_LENGTHS.get(assembly, {})  # plain name -> length

    try:
        tbx = pysam.TabixFile(source)
    except (OSError, ValueError) as e:
        raise SystemExit(
            f"cannot open {source} for random access ({e}). A tabix index is "
            f"required: expected {source}.tbi (download it from the same NCBI "
            f"directory as the .gz, or build it with:  python -c \"import "
            f"pysam; pysam.tabix_index('{source}', preset='vcf', force=True)\").")

    available = set(tbx.contigs)               # accession names present in index
    pool, weights = [], []
    for acc, name in accmap.items():
        if acc in available and lengths.get(name):
            pool.append((acc, name))
            weights.append(lengths[name])
    if not pool:
        tbx.close()
        raise SystemExit(
            f"none of the {assembly} primary contigs found in {source} "
            f"(index has e.g. {sorted(available)[:3]}); is --assembly correct?")

    counts = {t: 0 for t in VTYPES}
    seen_keys: set[tuple] = set()
    records: list[TruthRecord] = []
    target_total = per_type * len(VTYPES)
    attempts = 0

    def _need(vtype: str) -> bool:
        return counts[vtype] < per_type

    while attempts < max_attempts and sum(counts.values()) < target_total:
        attempts += 1
        (acc, name), = rng.choices(pool, weights=weights, k=1)
        length = lengths[name]
        start = rng.randint(1, max(1, length - 1))
        try:
            rows = tbx.fetch(acc, start - 1, min(length, start - 1 + window))
        except (ValueError, OSError):
            continue
        for line in rows:
            f = line.split("\t")
            if len(f) < 5:
                continue
            _acc, pos_s, id_field, ref, alt_field = f[0], f[1], f[2], f[3], f[4]
            rsid = _strip_rs(id_field)
            if rsid is None or not ref or ref == ".":
                continue
            alts = tuple(a for a in alt_field.split(",")
                         if a and a != "." and not a.startswith("<"))
            if not alts:
                continue
            try:
                pos = int(pos_s)
            except ValueError:
                continue
            ref = ref.upper()
            key = (name, pos, ref)
            if key in seen_keys:
                continue
            # Does this record help fill any still-open bucket?
            rec_types = [_classify(ref, a) for a in alts]
            if not any(_need(t) for t in rec_types):
                continue
            seen_keys.add(key)
            for a, t in zip(alts, rec_types):
                if counts[t] < per_type:
                    counts[t] += 1
            records.append(TruthRecord(
                rsid=rsid, chrom=name, pos=pos, ref=ref,
                alts=tuple(a.upper() for a in alts),
                multiallelic=len(alts) > 1))
            if sum(counts.values()) >= target_total:
                break
    tbx.close()
    rng.shuffle(records)
    return records


def sample_store_rsids(store_stem: str, k: int,
                       rng: random.Random) -> list[int]:
    """
    Harvest ``k`` real rsIDs straight from the store's rsID index (used for
    timing / round-trip when no source VCF is supplied). Reads random
    compressed blocks and picks random entries from them.
    """
    meta = MetaIndex.load(store_stem + ".rsid.meta")
    if not meta.blocks:
        return []
    out: list[int] = []
    with open(store_stem + ".rsid.idx", "rb") as fh:
        for _ in range(k):
            b = rng.choice(meta.blocks)
            fh.seek(b.comp_offset)
            entries = list(iter_entries(zlib.decompress(fh.read(b.comp_len))))
            if entries:
                out.append(rng.choice(entries)[0])
    return out


# ---------------------------------------------------------------------------
# curated landmarks (file-free smoke test).  Coordinates printed for a manual
# cross-check against https://www.ncbi.nlm.nih.gov/snp/ ; the automated check is
# chromosome + round-trip so it never false-fails on a coordinate revision.
# ---------------------------------------------------------------------------
CURATED = {
    "grch37": [
        {"rsid": 6025, "chrom": "1", "note": "F5 Leiden (F5), chr1"},
        {"rsid": 334, "chrom": "11", "note": "HbS sickle-cell (HBB), chr11"},
        {"rsid": 429358, "chrom": "19", "note": "APOE e4 site, chr19"},
    ],
    "grch38": [
        {"rsid": 6025, "chrom": "1", "note": "F5 Leiden (F5), chr1"},
        {"rsid": 334, "chrom": "11", "note": "HbS sickle-cell (HBB), chr11"},
        {"rsid": 429358, "chrom": "19", "note": "APOE e4 site, chr19"},
    ],
}


# ---------------------------------------------------------------------------
# the checks
# ---------------------------------------------------------------------------
def check_accuracy(client: DbSNPFlatQuery, assembly: str,
                   records: list[TruthRecord],
                   normalize: bool, max_report: int = 20) -> dict:
    """
    Both-direction accuracy + round-trip over the sampled truth records, run
    entirely through the backend query API (``DbSNPFlatQuery``): direction 1 via
    :meth:`variant_to_rsids`, direction 2 via :meth:`rsid_to_variants`.
    """
    norm = trim_alleles if normalize else (lambda c, p, r, a: (c, int(p), r, a))

    fwd_ok = fwd_tot = 0            # variant -> rsID
    rev_ok = rev_tot = 0           # rsID -> variant
    rt_ok = rt_tot = 0             # round-trip variant->rsID->...->rsID
    per_type = {t: {"fwd_ok": 0, "fwd_tot": 0, "rev_ok": 0, "rev_tot": 0}
                for t in VTYPES}
    failures: list[dict] = []

    for r in records:
        # ---- direction: sequence -> rsID (one case per ALT) --------------
        for a in r.alts:
            t = _classify(r.ref, a)
            fwd_tot += 1
            per_type[t]["fwd_tot"] += 1
            hits = client.variant_to_rsids(assembly, r.chrom, r.pos, r.ref, a).rsids
            if r.rsid in hits:
                fwd_ok += 1
                per_type[t]["fwd_ok"] += 1
            elif len(failures) < max_report:
                failures.append({
                    "direction": "variant->rsid", "type": t,
                    "query": f"{r.chrom}-{r.pos}-{r.ref}-{a}",
                    "expected_rsid": r.rsid, "got": list(hits)})

        # ---- direction: rsID -> variant (whole record at once) -----------
        rev_tot += 1
        vtype_rec = _classify(r.ref, r.alts[0])
        per_type[vtype_rec]["rev_tot"] += 1
        placements = client.rsid_to_variants(assembly, r.rsid).placements
        got = {(g.chrom, g.pos, g.ref, a) for g in placements for a in g.alts}
        # What the store *should* hold for this record: each ALT put through the
        # same normalizer the store was built with, giving (chrom, pos, ref, alt).
        expected = {norm(r.chrom, r.pos, r.ref, a) for a in r.alts}
        if expected <= got:
            rev_ok += 1
            per_type[vtype_rec]["rev_ok"] += 1
        elif len(failures) < max_report:
            failures.append({
                "direction": "rsid->variant", "type": vtype_rec,
                "query": f"rs{r.rsid}", "expected_loci": sorted(map(str, expected)),
                "got": sorted(f"{c}-{p}-{rf}-{al}" for c, p, rf, al in got)[:8]})

        # ---- round-trip: every returned variant must map back to rsID ----
        for (c, p, rf, al) in got:
            rt_tot += 1
            back = client.variant_to_rsids(assembly, c, p, rf, al).rsids
            if r.rsid in back:
                rt_ok += 1
            elif len(failures) < max_report:
                failures.append({
                    "direction": "round-trip", "query": f"{c}-{p}-{rf}-{al}",
                    "expected_rsid": r.rsid, "got": list(back)})

    return {
        "records": len(records),
        "forward": {"ok": fwd_ok, "total": fwd_tot,
                    "pass_rate": _rate(fwd_ok, fwd_tot)},
        "reverse": {"ok": rev_ok, "total": rev_tot,
                    "pass_rate": _rate(rev_ok, rev_tot)},
        "round_trip": {"ok": rt_ok, "total": rt_tot,
                       "pass_rate": _rate(rt_ok, rt_tot)},
        "by_type": {t: {**per_type[t],
                        "fwd_rate": _rate(per_type[t]["fwd_ok"],
                                          per_type[t]["fwd_tot"]),
                        "rev_rate": _rate(per_type[t]["rev_ok"],
                                          per_type[t]["rev_tot"])}
                    for t in VTYPES},
        "multiallelic_records": sum(1 for r in records if r.multiallelic),
        "failures": failures,
        "passed": (fwd_ok == fwd_tot and rev_ok == rev_tot and rt_ok == rt_tot),
    }


def check_formats(client: DbSNPFlatQuery, assembly: str,
                  records: list[TruthRecord], normalize: bool,
                  subset: int = 300, max_report: int = 20) -> dict:
    """
    Same answers must come back regardless of input spelling.

    variant -> rsID formats: canonical (the raw source spelling, which for real
    dbSNP indels is frequently non-minimal, so this alone exercises the
    normalizer), lower-case alleles, and (when the normalizer is on) an
    explicitly non-minimal shared-suffix spelling.
    rsID -> variant formats: 'rs6025' vs '6025' vs 6025.

    Note: a non-minimal *prefix* (left-context) spelling is deliberately NOT
    tested. ``trim_alleles`` is reference-free and strips shared suffixes first,
    so prepending an arbitrary base is not an equivalence-preserving transform
    for indels (e.g. T>TCT padded to AT>ATCT collapses to A>ATC). True
    left-alignment needs the reference genome; see the ``trim_alleles`` docs.
    """
    sample = records[:subset]
    variant_formats = ["canonical", "lowercase"]
    if normalize:
        variant_formats += ["nonminimal_suffix"]

    vres = {f: {"ok": 0, "tot": 0} for f in variant_formats}
    failures: list[dict] = []

    for r in sample:
        a = r.alts[0]
        forms = {
            "canonical": (r.chrom, r.pos, r.ref, a),
            "lowercase": (r.chrom, r.pos, r.ref.lower(), a.lower()),
        }
        if normalize:
            # Appending the SAME trailing base to REF and ALT is a legitimate
            # alternate spelling: trim_alleles strips the shared suffix and
            # recovers exactly trim(raw) (the AA>GA -> A>G example in its docs).
            forms["nonminimal_suffix"] = (r.chrom, r.pos, r.ref + "A", a + "A")
        for f, (c, p, rf, al) in forms.items():
            vres[f]["tot"] += 1
            if r.rsid in client.variant_to_rsids(assembly, c, p, rf, al).rsids:
                vres[f]["ok"] += 1
            elif len(failures) < max_report:
                failures.append({"format": f, "query": f"{c}-{p}-{rf}-{al}",
                                 "expected_rsid": r.rsid})

    # rsID input-format equivalence (exercises DbSNPFlatQuery's own rsID parse).
    rid_ok = rid_tot = 0
    for r in sample:
        rid_tot += 1
        a = set(client.rsid_to_variants(assembly, r.rsid).to_strings())
        b = set(client.rsid_to_variants(assembly, f"rs{r.rsid}").to_strings())
        c = set(client.rsid_to_variants(assembly, str(r.rsid)).to_strings())
        if a == b == c and a:
            rid_ok += 1
        elif len(failures) < max_report:
            failures.append({"format": "rsid-spelling", "query": f"rs{r.rsid}",
                             "note": "rs/int/str gave different results"})

    passed = (all(v["ok"] == v["tot"] for v in vres.values())
              and rid_ok == rid_tot)
    return {
        "variant_formats": {f: {**v, "pass_rate": _rate(v["ok"], v["tot"])}
                            for f, v in vres.items()},
        "rsid_spelling": {"ok": rid_ok, "total": rid_tot,
                          "pass_rate": _rate(rid_ok, rid_tot)},
        "failures": failures, "passed": passed,
    }


def check_curated(client: DbSNPFlatQuery, assembly: str) -> dict:
    """Well-known variants: right chromosome + round-trips. Soft (warn) check."""
    rows = []
    all_ok = True
    for item in CURATED.get(assembly, []):
        res = client.rsid_to_variants(assembly, item["rsid"])
        groups = res.placements
        strings = res.to_strings()
        chroms = {g.chrom for g in groups}
        chrom_ok = item["chrom"] in chroms
        rt_ok = all(item["rsid"]
                    in client.variant_to_rsids(assembly, g.chrom, g.pos,
                                               g.ref, a).rsids
                    for g in groups for a in g.alts) if groups else False
        ok = bool(groups) and chrom_ok and rt_ok
        all_ok = all_ok and ok
        rows.append({"rsid": item["rsid"], "note": item["note"],
                     "expected_chrom": item["chrom"], "resolved": strings,
                     "chrom_ok": chrom_ok, "round_trip_ok": rt_ok, "ok": ok})
    return {"landmarks": rows, "passed": all_ok}


def benchmark(make_client: Callable[[], DbSNPFlatQuery],
              client: DbSNPFlatQuery, assembly: str,
              variant_pool: list[tuple], rsid_pool: list[int],
              n_queries: int, rng: random.Random) -> dict:
    """
    Latency + throughput for both directions through the backend query API,
    plus a cold first-query time measured on a freshly opened client (its
    per-thread handle is opened lazily on that first call, so this captures the
    open+query cost a backend worker pays once).
    """
    def _cold(fn) -> Optional[float]:
        try:
            fresh = make_client()
        except Exception:
            return None
        try:
            t0 = time.perf_counter()
            fn(fresh)
            return (time.perf_counter() - t0) * 1000.0
        finally:
            fresh.close()

    out: dict = {"n_requested": n_queries}

    # sequence -> rsID
    if variant_pool:
        cold = _cold(lambda cl: cl.variant_to_rsids(assembly, *variant_pool[0]))
        picks = [rng.choice(variant_pool) for _ in range(n_queries)]
        lat = []
        for (c, p, rf, al) in picks:
            t0 = time.perf_counter()
            client.variant_to_rsids(assembly, c, p, rf, al)
            lat.append(time.perf_counter() - t0)
        out["variant_to_rsid"] = {**_summ_ms(lat), "cold_first_query_ms":
                                  round(cold, 4) if cold is not None else None}

    # rsID -> sequence
    if rsid_pool:
        cold = _cold(lambda cl: cl.rsid_to_variants(assembly, rsid_pool[0]))
        picks = [rng.choice(rsid_pool) for _ in range(n_queries)]
        lat = []
        for rid in picks:
            t0 = time.perf_counter()
            client.rsid_to_variants(assembly, rid)
            lat.append(time.perf_counter() - t0)
        out["rsid_to_variant"] = {**_summ_ms(lat), "cold_first_query_ms":
                                  round(cold, 4) if cold is not None else None}
    return out


def _rate(ok: int, tot: int) -> float:
    return round(100.0 * ok / tot, 4) if tot else float("nan")


# ---------------------------------------------------------------------------
# per-assembly driver
# ---------------------------------------------------------------------------
def validate_assembly(assembly: str, store_stem: str, source: Optional[str],
                      per_type: int, bench_queries: int, normalize: bool,
                      seed: int) -> dict:
    rng = random.Random(seed)

    for suffix in (".variants.tsv.gz", ".rsid.idx", ".rsid.meta"):
        if not os.path.exists(store_stem + suffix):
            raise SystemExit(f"missing store file: {store_stem + suffix}")

    normalizer = trim_alleles if normalize else None

    # The validation runs through the backend query API, so a pass certifies the
    # very code a VariantValidator integration will call.
    def _open() -> DbSNPFlatQuery:
        return DbSNPFlatQuery({assembly: store_stem}, normalizer=normalizer)

    client = _open()
    result: dict = {
        "assembly": assembly, "store": store_stem, "source": source,
        "normalize": normalize, "seed": seed,
        "sample_per_type": per_type, "bench_queries": bench_queries,
    }
    try:
        # --- truth sample -------------------------------------------------
        records: list[TruthRecord] = []
        if source:
            print(f"[{assembly}] sampling truth records from {source} ...")
            records = sample_source(source, assembly, per_type, rng)
            got = {t: 0 for t in VTYPES}
            for r in records:
                for a in r.alts:
                    got[_classify(r.ref, a)] += 1
            print(f"[{assembly}] sampled {len(records)} records "
                  f"({sum(1 for r in records if r.multiallelic)} multi-allelic); "
                  f"allele cases by type: {got}")
            result["accuracy"] = check_accuracy(client, assembly, records,
                                                normalize)
            result["formats"] = check_formats(client, assembly, records,
                                              normalize)
        else:
            print(f"[{assembly}] no source VCF -> skipping independent accuracy "
                  f"test (curated + round-trip + timing only)")
            result["accuracy"] = None
            result["formats"] = None

        # --- curated ------------------------------------------------------
        result["curated"] = check_curated(client, assembly)

        # --- build the query pools for the benchmark ----------------------
        if records:
            variant_pool = [(r.chrom, r.pos, r.ref, r.alts[0]) for r in records]
            rsid_pool = [r.rsid for r in records]
        else:
            rsid_pool = sample_store_rsids(store_stem, max(bench_queries, 500), rng)
            variant_pool = []
            for rid in rsid_pool[:1000]:
                ps = client.rsid_to_variants(assembly, rid).placements
                if ps:
                    g = ps[0]
                    variant_pool.append((g.chrom, g.pos, g.ref, g.alts[0]))
            # round-trip over the harvested rsIDs so there is still a
            # correctness signal without a source VCF.
            rt_ok = rt_tot = 0
            for rid in rsid_pool:
                for g in client.rsid_to_variants(assembly, rid).placements:
                    for a in g.alts:
                        rt_tot += 1
                        if rid in client.variant_to_rsids(
                                assembly, g.chrom, g.pos, g.ref, a).rsids:
                            rt_ok += 1
            result["store_round_trip"] = {
                "ok": rt_ok, "total": rt_tot, "pass_rate": _rate(rt_ok, rt_tot),
                "passed": rt_ok == rt_tot}

        # --- timing -------------------------------------------------------
        print(f"[{assembly}] benchmarking {bench_queries} queries per direction ...")
        result["timing"] = benchmark(_open, client, assembly, variant_pool,
                                     rsid_pool, bench_queries, rng)
    finally:
        client.close()

    # --- overall verdict --------------------------------------------------
    checks = []
    if result.get("accuracy") is not None:
        checks.append(result["accuracy"]["passed"])
    if result.get("formats") is not None:
        checks.append(result["formats"]["passed"])
    if "store_round_trip" in result:
        checks.append(result["store_round_trip"]["passed"])
    # curated is a soft/warn check, excluded from the hard verdict.
    result["passed"] = all(checks) if checks else None
    result["curated_passed"] = result["curated"]["passed"]
    return result


# ---------------------------------------------------------------------------
# reporting
# ---------------------------------------------------------------------------
def _fmt_timing(t: dict) -> str:
    if not t or t.get("n", 0) == 0:
        return "  (no queries)\n"
    return ("    n={n}  min={min_ms}  mean={mean_ms}  median={median_ms}  "
            "p95={p95_ms}  p99={p99_ms}  max={max_ms}  "
            "({throughput_qps} q/s)  cold_first={cold}\n").format(
        cold=t.get("cold_first_query_ms"), **t)


def print_report(results: list[dict]) -> None:
    for r in results:
        a = r["assembly"]
        print("\n" + "=" * 72)
        print(f"ASSEMBLY {a.upper()}   store={r['store']}")
        print(f"  source={r['source']}  normalize={r['normalize']}  "
              f"seed={r['seed']}")
        print("=" * 72)

        acc = r.get("accuracy")
        if acc:
            print(f"ACCURACY (independent truth from source VCF, "
                  f"{acc['records']} records, "
                  f"{acc['multiallelic_records']} multi-allelic)")
            print(f"  variant->rsID : {acc['forward']['ok']}/"
                  f"{acc['forward']['total']}  ({acc['forward']['pass_rate']}%)")
            print(f"  rsID->variant : {acc['reverse']['ok']}/"
                  f"{acc['reverse']['total']}  ({acc['reverse']['pass_rate']}%)")
            print(f"  round-trip    : {acc['round_trip']['ok']}/"
                  f"{acc['round_trip']['total']}  "
                  f"({acc['round_trip']['pass_rate']}%)")
            print("  by variant type (fwd / rev pass-rate):")
            for t in VTYPES:
                bt = acc["by_type"][t]
                print(f"    {t:10} fwd {bt['fwd_ok']}/{bt['fwd_tot']} "
                      f"({bt['fwd_rate']}%)   rev {bt['rev_ok']}/{bt['rev_tot']} "
                      f"({bt['rev_rate']}%)")
            if acc["failures"]:
                print(f"  !! {len(acc['failures'])} failure(s) shown:")
                for f in acc["failures"][:10]:
                    print(f"     {f}")
        else:
            print("ACCURACY: skipped (no source VCF)")
            if "store_round_trip" in r:
                srt = r["store_round_trip"]
                print(f"  store round-trip: {srt['ok']}/{srt['total']} "
                      f"({srt['pass_rate']}%)")

        fmt = r.get("formats")
        if fmt:
            print("FORMAT ROBUSTNESS")
            for f, v in fmt["variant_formats"].items():
                print(f"  variant/{f:14} {v['ok']}/{v['tot']} "
                      f"({v['pass_rate']}%)")
            rs = fmt["rsid_spelling"]
            print(f"  rsID spelling      {rs['ok']}/{rs['total']} "
                  f"({rs['pass_rate']}%)")

        cur = r["curated"]
        print(f"CURATED LANDMARKS (soft check) "
              f"{'OK' if cur['passed'] else 'WARN'}")
        for row in cur["landmarks"]:
            flag = "ok " if row["ok"] else "WARN"
            print(f"  [{flag}] rs{row['rsid']:<8} chr{row['expected_chrom']:<3} "
                  f"{row['note']}")
            for s in row["resolved"][:4]:
                print(f"          -> {s}")

        tim = r.get("timing", {})
        print("TIMING")
        print("  variant->rsID:")
        print(_fmt_timing(tim.get("variant_to_rsid", {})), end="")
        print("  rsID->variant:")
        print(_fmt_timing(tim.get("rsid_to_variant", {})), end="")

        verdict = ("PASS" if r["passed"] else "FAIL") if r["passed"] is not None \
            else "N/A (no hard checks run)"
        print(f"\nVERDICT [{a}]: {verdict}")


def write_reports(results: list[dict], report_dir: str) -> tuple[str, str]:
    os.makedirs(report_dir, exist_ok=True)
    ts = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    asm = "_".join(r["assembly"] for r in results)
    json_path = os.path.join(report_dir, f"validation_{asm}_{ts}.json")
    md_path = os.path.join(report_dir, f"validation_{asm}_{ts}.md")

    payload = {"generated_utc": ts, "tool": "vvDBFlatValidate.py",
               "results": results}
    with open(json_path, "w") as fh:
        json.dump(payload, fh, indent=2, default=str)

    with open(md_path, "w") as fh:
        fh.write(f"# Flat dbSNP validation report\n\n")
        fh.write(f"_Generated {ts} by `vvDBFlatValidate.py`._\n\n")
        for r in results:
            _write_md_assembly(fh, r)
    return json_path, md_path


def _md_table(fh, header: list[str], rows: list[list]) -> None:
    fh.write("| " + " | ".join(header) + " |\n")
    fh.write("| " + " | ".join("---" for _ in header) + " |\n")
    for row in rows:
        fh.write("| " + " | ".join(str(c) for c in row) + " |\n")
    fh.write("\n")


def _write_md_assembly(fh, r: dict) -> None:
    a = r["assembly"]
    verdict = ("PASS" if r["passed"] else "FAIL") if r["passed"] is not None \
        else "N/A"
    fh.write(f"## {a.upper()} — **{verdict}**\n\n")
    fh.write(f"- store: `{r['store']}`\n")
    fh.write(f"- source: `{r['source']}`\n")
    fh.write(f"- normalize: `{r['normalize']}`  seed: `{r['seed']}`  "
             f"sample/type: `{r['sample_per_type']}`  "
             f"bench queries: `{r['bench_queries']}`\n\n")

    acc = r.get("accuracy")
    if acc:
        fh.write(f"### Accuracy (independent truth, {acc['records']} records, "
                 f"{acc['multiallelic_records']} multi-allelic)\n\n")
        _md_table(fh, ["direction", "passed/total", "pass rate"],
                  [["variant → rsID",
                    f"{acc['forward']['ok']}/{acc['forward']['total']}",
                    f"{acc['forward']['pass_rate']}%"],
                   ["rsID → variant",
                    f"{acc['reverse']['ok']}/{acc['reverse']['total']}",
                    f"{acc['reverse']['pass_rate']}%"],
                   ["round-trip",
                    f"{acc['round_trip']['ok']}/{acc['round_trip']['total']}",
                    f"{acc['round_trip']['pass_rate']}%"]])
        fh.write("By variant type:\n\n")
        _md_table(fh, ["type", "variant→rsID", "rsID→variant"],
                  [[t, f"{acc['by_type'][t]['fwd_ok']}/"
                       f"{acc['by_type'][t]['fwd_tot']} "
                       f"({acc['by_type'][t]['fwd_rate']}%)",
                    f"{acc['by_type'][t]['rev_ok']}/"
                    f"{acc['by_type'][t]['rev_tot']} "
                    f"({acc['by_type'][t]['rev_rate']}%)"] for t in VTYPES])
        if acc["failures"]:
            fh.write(f"**{len(acc['failures'])} failure(s):**\n\n```\n")
            for f in acc["failures"][:20]:
                fh.write(f"{f}\n")
            fh.write("```\n\n")
    else:
        fh.write("### Accuracy\n\nSkipped (no source VCF supplied).\n\n")
        if "store_round_trip" in r:
            srt = r["store_round_trip"]
            fh.write(f"Store round-trip: {srt['ok']}/{srt['total']} "
                     f"({srt['pass_rate']}%)\n\n")

    fmt = r.get("formats")
    if fmt:
        fh.write("### Format robustness\n\n")
        rows = [[f"variant / {f}", f"{v['ok']}/{v['tot']}",
                 f"{v['pass_rate']}%"]
                for f, v in fmt["variant_formats"].items()]
        rows.append(["rsID spelling (rs/int/str)",
                     f"{fmt['rsid_spelling']['ok']}/"
                     f"{fmt['rsid_spelling']['total']}",
                     f"{fmt['rsid_spelling']['pass_rate']}%"])
        _md_table(fh, ["input format", "passed/total", "pass rate"], rows)

    cur = r["curated"]
    fh.write(f"### Curated landmarks (soft check — "
             f"{'OK' if cur['passed'] else 'WARN'})\n\n")
    _md_table(fh, ["rsID", "expected chr", "chr ok", "round-trip ok",
                   "resolved", "note"],
              [[f"rs{row['rsid']}", row["expected_chrom"], row["chrom_ok"],
                row["round_trip_ok"], "; ".join(row["resolved"][:3]) or "(none)",
                row["note"]] for row in cur["landmarks"]])

    tim = r.get("timing", {})
    fh.write("### Timing (latency per query)\n\n")
    rows = []
    for label, key in (("variant → rsID", "variant_to_rsid"),
                       ("rsID → variant", "rsid_to_variant")):
        t = tim.get(key)
        if t and t.get("n"):
            rows.append([label, t["n"], f"{t['min_ms']}", f"{t['mean_ms']}",
                         f"{t['median_ms']}", f"{t['p95_ms']}", f"{t['p99_ms']}",
                         f"{t['max_ms']}", f"{t['throughput_qps']}",
                         t.get("cold_first_query_ms")])
    if rows:
        _md_table(fh, ["direction", "n", "min ms", "mean ms", "median ms",
                       "p95 ms", "p99 ms", "max ms", "q/s", "cold 1st ms"], rows)
    fh.write("\n")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def _main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    # Simple single-assembly interface (run grch37 and grch38 as separate
    # commands, each with its own report):
    #     python vvDBFlatValidate.py --assembly grch38 --store <stem> --source <vcf>
    ap.add_argument("--assembly", choices=["grch37", "grch38"],
                    help="Validate just this assembly (use with --store/--source)")
    ap.add_argument("--store", help="Store path stem for --assembly")
    ap.add_argument("--source", help="Source dbSNP VCF (bgzip+tbi) for --assembly")
    # Paired interface (validate one or both assemblies in a single run):
    ap.add_argument("--grch37-store", help="GRCh37 store path stem")
    ap.add_argument("--grch37-source", help="GRCh37 source dbSNP VCF (bgzip+tbi)")
    ap.add_argument("--grch38-store", help="GRCh38 store path stem")
    ap.add_argument("--grch38-source", help="GRCh38 source dbSNP VCF (bgzip+tbi)")
    ap.add_argument("--sample-per-type", type=int, default=200,
                    help="Truth records to sample per variant type (default 200)")
    ap.add_argument("--bench-queries", type=int, default=2000,
                    help="Timed queries per direction (default 2000)")
    ap.add_argument("--seed", type=int, default=42,
                    help="RNG seed for reproducible sampling (default 42)")
    ap.add_argument("--no-normalise", action="store_true",
                    help="Disable normalization (must match how the store was "
                         "built)")
    ap.add_argument("--report-dir", default="validation_reports",
                    help="Where to write the JSON + Markdown report "
                         "(default: validation_reports)")
    ap.add_argument("--no-report-files", action="store_true",
                    help="Print to console only; do not write report files")
    args = ap.parse_args(argv)

    jobs = []
    if args.assembly:
        if not args.store:
            ap.error("--assembly requires --store")
        if args.grch37_store or args.grch38_store:
            ap.error("use either --assembly/--store/--source OR the "
                     "--grch37-*/--grch38-* flags, not both")
        jobs.append((args.assembly, args.store, args.source))
    else:
        if args.store or args.source:
            ap.error("--store/--source require --assembly")
        if args.grch37_store:
            jobs.append(("grch37", args.grch37_store, args.grch37_source))
        if args.grch38_store:
            jobs.append(("grch38", args.grch38_store, args.grch38_source))
    if not jobs:
        ap.error("Specify --assembly/--store (single) or "
                 "--grch37-store / --grch38-store (paired)")

    normalize = not args.no_normalise
    results = []
    for assembly, store, source in jobs:
        if source and not os.path.exists(source):
            sys.exit(f"ERROR: {assembly} source not found: {source}")
        results.append(validate_assembly(
            assembly, store, source,
            per_type=args.sample_per_type,
            bench_queries=args.bench_queries,
            normalize=normalize, seed=args.seed))

    print_report(results)

    if not args.no_report_files:
        json_path, md_path = write_reports(results, args.report_dir)
        print(f"\nReports written:\n  {json_path}\n  {md_path}")

    # Non-zero exit if any hard check failed.
    hard = [r["passed"] for r in results if r["passed"] is not None]
    return 0 if all(hard) else 1


if __name__ == "__main__":
    raise SystemExit(_main())

