# VariantValidator populate_vdb

Repo that populates the validator database and also provides a flat dbSNP store for fast rsID ↔ variant lookups.

# Flat dbSNP store: build, query, validate

Fast **two-way conversion between dbSNP rsIDs and genomic variants**
(`rs6025 ⇄ 1-169549811-C-T`), without a database server and without scanning the
~30 GB dbSNP VCF at query time.

From the official dbSNP VCF for a genome build, `vvDBFlatBuildDbSNP.py` produces a
small, indexed **flat store** that answers both directions in well under a
millisecond:

- **rsID → variant** ("where is `rs6025`?")
- **variant → rsID** ("does `1-169549811-C-T` have an rsID?")

## The three modules

| File | Role |
|---|---|
| [`vvDBFlatBuildDbSNP.py`](vvDBFlatBuildDbSNP.py) | On-disk format + reader + **builder** (turns a dbSNP VCF into the store) |
| [`vvDBFlatQueryDbSNP.py`](vvDBFlatQueryDbSNP.py) | **Backend query API** (importable library, both directions, thread-safe) |
| [`vvDBFlatValidate.py`](vvDBFlatValidate.py) | **Validator** (proves the store round-trips dbSNP, both directions, + timing) |

Full, step-by-step explanations of each — with worked examples and references to
every function/class — are in the **[`docs/`](docs/)** directory:

- [`docs/vvDBFlatBuildDbSNP.md`](docs/vvDBFlatBuildDbSNP.md)
- [`docs/vvDBFlatQueryDbSNP.md`](docs/vvDBFlatQueryDbSNP.md)
- [`docs/vvDBFlatValidate.md`](docs/vvDBFlatValidate.md)

---

## Requirements : runs in the VariantValidator venv, no new dependencies

The only third-party dependency is **`pysam`**, which VariantValidator already
uses. Everything else is the Python standard library plus the system `sort`
command (present on macOS/Linux). So this project runs **inside the existing
VariantValidator virtual environment (`vvenv`) as-is** — activating that env is
enough, nothing new to install.

```bash
conda activate vvenv     # or however you enter the VariantValidator env
```

(`bcftools` is optional and only used if you build directly from a `.bcf`; the
default VCF path uses `pysam` and needs no external binary. `pytest` is needed
only to run the test suite.)

---

## Running the test suite

The suite is **self-contained**: it builds a tiny synthetic store in a fixture,
so it needs no network and no downloads.

```bash
# from this directory:
python3 -m pytest
```

Useful variants:

```bash
python3 -m pytest tests/test_query.py          # one module
python3 -m pytest tests/test_build.py -v       # verbose
python3 -m pytest tests/test_build.py::test_trim_alleles -v
```

What it contains (**127 tests**):

| File | Tests | Covers |
|---|---|---|
| `tests/test_build.py` | 51 | normalization, contig maps, on-disk format, the builder, the `DbSNPFlat` reader (both directions), the CLI |
| `tests/test_query.py` | 50 | input parsing, result objects, `DbSNPFlatQuery` construction/errors, both directions, thread-safety, lifecycle |
| `tests/test_validate.py` | 26 | helpers, source/store sampling, accuracy/format/curated checks (incl. negative tests), benchmark, `validate_assembly`, reporting, CLI |

`tests/conftest.py` writes an NCBI-style source VCF (RefSeq-accession contigs, no
`##contig` header, bgzipped + tabix-indexed) with known landmark variants, then
runs the real builder on it once per session. See
[`tests/README.md`](tests/README.md).

---

## Step 1 : Download the dbSNP VCF and build the store

### Download

dbSNP publishes one whole-genome VCF per genome build. The current build (**b157**,
which is also `latest_release`) lives at NCBI:

| Assembly | File | Approx size |
|---|---|---|
| GRCh37 | `GCF_000001405.25.gz` | ~28 GB |
| GRCh38 | `GCF_000001405.40.gz` | ~30 GB |

Download the `.gz`, its **`.tbi`** index (the validator needs it for random
access), and the `.md5` to verify. Example for GRCh38:

```bash
BASE=https://ftp.ncbi.nih.gov/snp/latest_release/VCF
curl -O $BASE/GCF_000001405.40.gz
curl -O $BASE/GCF_000001405.40.gz.tbi
curl -O $BASE/GCF_000001405.40.gz.md5
md5sum -c GCF_000001405.40.gz.md5     # (macOS: md5 GCF_000001405.40.gz, compare by eye)
```

(Swap `.40` → `.25` for GRCh37. To pin a specific build, use
`https://ftp.ncbi.nih.gov/snp/archive/b157/VCF/…` instead of `latest_release`.)

### Build

```bash
python3 vvDBFlatBuildDbSNP.py --out dbsnp_flat --version b157 \
    --grch38 /path/GCF_000001405.40.gz \
    --grch37 /path/GCF_000001405.25.gz
```

You can pass just one assembly. The build streams the VCF once, renames RefSeq
accessions to plain chromosome names (`NC_000001.11 → 1`), normalizes alleles,
and writes the store. (Runtime is roughly tens of minutes per assembly.)

### Outputs

For each assembly, four files under `dbsnp_flat/<version>/`:

```
dbsnp_flat/b157/
├── dbsnp_grch38_b157.variants.tsv.gz       # MASTER: sorted, bgzipped rows
├── dbsnp_grch38_b157.variants.tsv.gz.tbi   #   its tabix index  → variant → rsID
├── dbsnp_grch38_b157.rsid.idx              # rsID index (compressed blocks)
└── dbsnp_grch38_b157.rsid.meta             #   its signpost table → rsID → variant
```

The shared prefix `dbsnp_flat/b157/dbsnp_grch38_b157` is the **store "stem"** you
pass to the query API and the validator. Building is idempotent — re-running skips
an assembly whose four files already exist.

---

## Step 2 : Validate the store

```bash
python3 vvDBFlatValidate.py --assembly grch38 \
    --store  dbsnp_flat/b157/dbsnp_grch38_b157 \
    --source /path/GCF_000001405.40.gz
```

Run GRCh37 as a separate command with its own `--store`/`--source`, or validate
both at once with the paired flags (`--grch37-store/--grch37-source/…`).

### The data and truth set it uses

- **Data under test:** the built flat store (the four files above), queried
  **through the real backend API** (`DbSNPFlatQuery`) ; so a pass certifies the
  exact code a VariantValidator integration will call.
- **Truth set:** the **original dbSNP source VCF, re-read independently**
  (`--source`). The validator randomly samples real records straight from that
  VCF — stratified across variant types (SNV, insertion, deletion, MNV) and
  including multi-allelic sites — and checks the store round-trips each one in
  **both directions**. Because the store is *built* from that VCF, comparing it
  to an independent re-parse of the same VCF is what makes this a real accuracy
  test rather than a tautology. (The source needs its `.tbi` for the random
  sampling. Without `--source` it still runs curated + round-trip + timing using
  rsIDs harvested from the store itself.)

### Outputs

- A **console report**: accuracy per direction and per variant type, format
  robustness, curated landmarks, timing (latency percentiles + throughput), and a
  PASS/FAIL verdict.
- **Report files** in `validation_reports/` (override with `--report-dir`;
  suppress with `--no-report-files`): a machine-readable
  `validation_<assembly>_<timestamp>.json` and a documentation-friendly
  `.md` table.
- **Exit code** `0` on pass, non-zero if any accuracy / format / round-trip check
  fails — so it can gate CI.

> **Scope note.** A pass means the store is a **faithful, lossless re-indexing of
> dbSNP** for that build. It does *not* certify dbSNP itself as biologically or
> clinically correct, and the default reference-free normalizer is not full
> left-alignment (some coordinate→rsID indel lookups can miss). See the caveats in
> [`docs/vvDBFlatValidate.md`](docs/vvDBFlatValidate.md).

---

## Step 3 : Integrating the query API into the VariantValidator backend

`vvDBFlatQueryDbSNP.py` is a **library, not a CLI** ; import it and call the two
directions per request.

**Once, at app/process startup** (paths from your config):

```python
from vvDBFlatQueryDbSNP import DbSNPFlatQuery

# discover stores under the builder's layout …
DBSNP = DbSNPFlatQuery.from_directory("/data/dbsnp_flat", version="b157")

# … or pass explicit stems:
DBSNP = DbSNPFlatQuery({
    "grch38": "/data/dbsnp_flat/b157/dbsnp_grch38_b157",
    "grch37": "/data/dbsnp_flat/b157/dbsnp_grch37_b157",
})
```

Construct it **once** and reuse it; it checks the files exist up front, and it is
**thread-safe** (each worker thread lazily gets its own file handles), so it fits
a request-serving backend with a worker pool.

**Then, per request** : the two directions are separate methods returning
structured objects (not printed strings):

```python
# variant → rsID
hit = DBSNP.variant_to_rsids("grch38", "1", 169549811, "C", "T")
if hit.found:
    rsid = hit.first            # 6025   (hit.rsids / hit.labels for all matches)

# rsID → variant  (accepts 6025 / "6025" / "rs6025")
res = DBSNP.rsid_to_variants("grch38", "rs6025")
for chrom, pos, ref, alt in res.to_vcf_tuples():
    ...                          # ('1', 169549811, 'C', 'T'), plus any other ALTs
```

Integration notes for the backend developer:

- A variant/rsID that simply isn't present returns an **empty result**
  (`found == False`) , not an exception. Malformed input (bad rsID, empty allele,
  unknown assembly) raises, so callers can distinguish "not found" from "bad
  request".
- The **`normalizer` must match how the store was built** (default
  `trim_alleles`). This is the seam for plugging in VariantValidator's own
  reference-aware normalizer later ; pass it to *both* the builder and
  `DbSNPFlatQuery` so build-time and query-time representations agree.
- Results are frozen dataclasses (`VariantToRsidResult`, `RsidToVariantResult`,
  `VariantPlacement`) ; easy to log or serialize to JSON in an endpoint.

Full API reference and behaviour: [`docs/vvDBFlatQueryDbSNP.md`](docs/vvDBFlatQueryDbSNP.md).
