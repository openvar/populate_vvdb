# `vvDBFlatValidate.py`

## What this file is for (the big picture)

After the store is built (`vvDBFlatBuildDbSNP.py`) and wrapped in a query API
(`vvDBFlatQueryDbSNP.py`), how do we *prove* the querying actually works? This
file is the **validator**: it checks the store for **accuracy** (are the answers
right?) and **speed** (how fast is each query?), and writes a human-readable
report anyone can reproduce.

Two design choices make the check trustworthy:

1. **The "truth" is the original dbSNP VCF, re-read independently.** The store is
   built *from* that VCF, so comparing the store to itself would prove nothing.
   The validator samples real records straight from the source VCF and checks the
   store round-trips them.
2. **Every check runs through the real backend API (`DbSNPFlatQuery`).** So a
   passing run certifies the *actual functions* a VariantValidator integration
   will call; not a look-alike.

> Scope note: a pass means the store is a **faithful re-indexing of dbSNP** for
> that build. It does **not** prove dbSNP is biologically/clinically correct, and
> the reference-free normalizer means some coordinate→rsID indel lookups can miss
> — see the caveats in the code comments.

---

## Map of the file (what lives where)

| Part | Functions / classes | Job |
|---|---|---|
| Build a truth set | `sample_source`, `TruthRecord`, `CONTIG_LENGTHS`, `_classify` | pull real records from the source VCF |
| Accuracy | `check_accuracy` | both directions + round-trip agree with truth |
| Formats | `check_formats` | different input spellings all work |
| Landmarks | `check_curated`, `CURATED` | famous variants land on the right chromosome |
| Speed | `benchmark`, `_summ_ms`, `_pctl` | latency + throughput per direction |
| No-source mode | `sample_store_rsids` | still test without the source VCF |
| Orchestration | `validate_assembly` | run all checks for one assembly |
| Reporting | `print_report`, `write_reports` | console + JSON + Markdown |
| Command line | `_main` | run it from a terminal |

---

## Step 1 : Build an independent "truth set" (`sample_source`)

**`sample_source`** pulls real records straight from the source dbSNP VCF (e.g.
`GCF_000001405.40.gz`), **without scanning all 30 GB**. It uses the VCF's tabix
index to jump to random windows and grabs records there.

Two real-world wrinkles it handles:

- **NCBI dbSNP VCFs don't list chromosomes (or their sizes) in the header.** So
  the validator reads the chromosome list from the **tabix index** and carries
  the chromosome lengths itself in the constant **`CONTIG_LENGTHS`** (standard
  GRCh37/GRCh38 sizes). The lengths let it pick random positions weighted by
  chromosome size.
- It reads records as **raw text** via `pysam.TabixFile` (which also avoids noisy
  parser warnings).

It keeps sampling until it has a **balanced mix of variant types**, classified by
**`_classify`** into the four buckets in `VTYPES`:

| Type | Meaning | Example |
|---|---|---|
| `snv` | one base → one base | `C → T` |
| `insertion` | ALT longer than REF | `A → ATG` |
| `deletion` | REF longer than ALT | `AGT → A` |
| `mnv` | equal length > 1 | `AGT → CGA` |

Each kept record becomes a **`TruthRecord`** : the ground truth, already
chromosome-mapped to plain names:

```python
TruthRecord(rsid=6025, chrom="1", pos=169549811, ref="C",
            alts=("T",), multiallelic=False)
```

A fixed random seed makes the sample reproducible.

---

## Step 2 : Accuracy: both directions + round-trip (`check_accuracy`)

**`check_accuracy`** takes the truth records and, for each one, asks the backend
API three questions.

Take the truth record `rs6025 = 1-169549811-C-T`:

**(a) variant → rsID** : feed the coordinates, expect the rsID back:
```python
client.variant_to_rsids("grch38", "1", 169549811, "C", "T").rsids   # must contain 6025
```

**(b) rsID → variant** : feed the rsID, expect the coordinates back. For a
multi-allelic truth record like `rs334 = 11-5227002-T-{A,C,G}`, it checks the
store returns **all** three alternates:
```python
client.rsid_to_variants("grch38", 334).placements   # must cover 11-5227002-T-A/C/G
```
(It compares against the truth run through the *same normalizer* the store uses,
so non-minimal spellings line up.)

**(c) round-trip** : take what direction (b) returned and feed it back into
direction (a); the original rsID must come out. This catches index/offset
corruption.

The result records pass-rates **overall and per variant type**, plus a few
example failures. A perfect store scores 100% on all three. Helper **`_rate`**
turns "132 of 132" into "100.0%".

---

## Step 3 : Format robustness (`check_formats`)

A clinician (or another pipeline) may spell the same query differently. This step
checks the answer doesn't change. For each sampled record it re-issues the
**variant → rsID** query in several spellings:

| Format | Example for `1-169549811-C-T` | Should still find rs6025? |
|---|---|---|
| `canonical` | `C / T` (raw source spelling) | yes |
| `lowercase` | `c / t` | yes |
| `nonminimal_suffix` | `CA / TA` (shared trailing base) | yes (normalizer trims it) |

For the **rsID → variant** direction it checks the three rsID spellings —
`6025`, `"rs6025"`, `"6025"` — all return the **same** answer (this exercises the
client's `parse_rsid`).

> Note in the code: a non-minimal *prefix* spelling is deliberately **not**
> tested, because reference-free `trim_alleles` can't treat it as equivalent for
> indels ; that needs the genome sequence. 

---

## Step 4 : Curated landmarks (`check_curated`)

**`check_curated`** looks up a tiny hand-picked list (**`CURATED`**) of famous
variants : `rs6025` (Factor V Leiden, chr1), `rs334` (sickle-cell, chr11),
`rs429358` (APOE, chr19). For each it checks the store returns something **on the
expected chromosome** and that it **round-trips**, and it prints the actual
coordinates so a human can eyeball them against the NCBI website.

This is a **soft check** (labelled `WARN`, not counted in the pass/fail verdict):
it's a sanity smoke test and a documentation aid, not a strict gate ; a synthetic
test store won't contain `rs6025`, and that shouldn't fail the run.

---

## Step 5 : Speed (`benchmark`)

**`benchmark`** times many random queries in **each direction separately** and
summarises the latencies with **`_summ_ms`** (which uses **`_pctl`** for
percentiles). It reports:

```
variant→rsID:  n=2000  min  mean  median  p95  p99  max  (… queries/sec)  cold_first=…
rsID→variant:  n=2000  …
```

- `median` / `p95` / `p99` tell you the typical and worst-case latency.
- `cold_first` is measured on a **freshly opened client**, capturing the one-time
  "open the files" cost a backend worker pays on its first query.

---

## Step 6 : Running without the source VCF (`sample_store_rsids`)

If you don't have the 30 GB source handy, `validate_assembly` still does
something useful. **`sample_store_rsids`** harvests real rsIDs directly from the
store's own index (reading random compressed blocks), and the validator then runs
**round-trip consistency** (rsID → variant → rsID) plus the timing benchmark on
them. You lose the *independent* accuracy check (there's no external truth), but
you still catch index corruption and get performance numbers.

---

## Step 7 : Putting one assembly together (`validate_assembly`)

**`validate_assembly`** is the conductor for a single genome build. It:

1. Opens the backend client: `DbSNPFlatQuery({assembly: store_stem}, normalizer=…)`.
2. If a source VCF was given: `sample_source` → `check_accuracy` → `check_formats`.
   Otherwise: the no-source path (Step 6).
3. Always: `check_curated`.
4. Always: `benchmark`.
5. Computes an overall **verdict**: pass only if accuracy, formats, and
   round-trip all pass (curated is excluded, being soft).

Everything flows through the **same `DbSNPFlatQuery`** object, so the validation
truly exercises the backend query path.

---

## Step 8 : Reporting (`print_report`, `write_reports`)

- **`print_report`** prints the sectioned summary to the console.
- **`write_reports`** writes a machine-readable **`.json`** and a
  documentation-friendly **`.md`** (tables built by `_md_table`) into the report
  directory, timestamped, e.g. `validation_grch38_20260901T…Z.md`.

---

## Step 9 : The command line (`_main`)

**`_main`** is the terminal entry point. Two ways to call it:

**One assembly at a time** (each writes its own report):
```bash
python vvDBFlatValidate.py --assembly grch38 \
    --store  /data/dbsnp_flat/b157/dbsnp_grch38_b157 \
    --source /data/GCF_000001405.40.gz
```

**Both at once** with the paired flags:
```bash
python vvDBFlatValidate.py \
    --grch37-store … --grch37-source … \
    --grch38-store … --grch38-source …
```

Useful options: `--sample-per-type` (how many of each variant type to test),
`--bench-queries` (how many timed queries), `--seed` (reproducibility),
`--no-normalise` (only if the store was built that way), `--report-dir`,
`--no-report-files`.

**Exit code:** `0` if every hard check passed, non-zero if any accuracy / format
/ round-trip check failed — so it can gate a CI pipeline.

---

## What a passing run looks like

```
ACCURACY (independent truth from source VCF, 800 records, 96 multi-allelic)
  variant->rsID : 946/946  (100.0%)
  rsID->variant : 800/800  (100.0%)
  round-trip    : 946/946  (100.0%)
FORMAT ROBUSTNESS
  variant/canonical         800/800 (100.0%)
  variant/lowercase         800/800 (100.0%)
  variant/nonminimal_suffix 800/800 (100.0%)
  rsID spelling             800/800 (100.0%)
TIMING
  variant->rsID:  median≈0.02 ms  …
  rsID->variant:  median≈0.09 ms  …
VERDICT [grch38]: PASS
```

---

