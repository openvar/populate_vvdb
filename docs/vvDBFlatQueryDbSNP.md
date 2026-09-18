# `vvDBFlatQueryDbSNP.py` : explained step by step

## What this file is for (the big picture)

`vvDBFlatBuildDbSNP.py` already contains a reader (`DbSNPFlat`) that can query the
flat store. So why another file?

Because a raw reader is not quite what a **backend** (like VariantValidator)
needs. This file is a **thin, safe wrapper**. This file is a *library, not a command-line
tool*, meant to be imported and called once per request. On top of the raw
reader it adds the four things a production backend needs:

1. **Manages both genome builds** (GRCh37 and GRCh38) behind one object.
2. **Cleans up the inputs** : accepts an rsID as `6025`, `"rs6025"`, or `"6025"`,
   and an assembly as `"grch38"`, `"hg38"`, or `"38"`.
3. **Returns tidy result objects** (not printed text) that backend code can read
   fields off, or turn into JSON.
4. **Is safe to use from many request threads at once** (the raw reader is not).

It exposes the two directions as **two separate methods**, exactly as asked:

```
variant_to_rsids(assembly, chrom, pos, ref, alt)   →  sequence → rsID
rsid_to_variants(assembly, rsid)                    →  rsID → sequence
```

---

## Map of the file (what lives where)

| Part | Functions / classes | Job |
|---|---|---|
| Input cleaning | `normalise_assembly`, `parse_rsid`, `_require_allele`, `store_stem` | validate and standardise inputs/paths |
| Result objects | `VariantQuery`, `VariantToRsidResult`, `VariantPlacement`, `RsidToVariantResult` | structured answers |
| The client | `DbSNPFlatQuery` | holds the stores, runs the two lookups, thread-safe |

---

## Step 1 : Create the client once (at startup)

You build **one `DbSNPFlatQuery`** when the app starts and reuse it for every
request. Two ways to construct it:

```python
from vvDBFlatQueryDbSNP import DbSNPFlatQuery

# A) point it at the builder's folder layout:
DBSNP = DbSNPFlatQuery.from_directory("/data/dbsnp_flat", version="b157")

# B) or give explicit path stems per assembly:
DBSNP = DbSNPFlatQuery({
    "grch38": "/data/dbsnp_flat/b157/dbsnp_grch38_b157",
    "grch37": "/data/dbsnp_flat/b157/dbsnp_grch37_b157",
})
```

- **`DbSNPFlatQuery.from_directory`** works out each stem with the helper
  **`store_stem`** (`<base>/<version>/dbsnp_<assembly>_<version>`) and skips any
  assembly whose files are missing (unless `require_all=True`).
- **`DbSNPFlatQuery.__init__`** canonicalises each label with
  **`normalise_assembly`** and **checks the three files exist right now**, so a
  bad path fails loudly at startup instead of mid-request.

A "**stem**" is the shared path prefix *without* the `.variants.tsv.gz` /
`.rsid.idx` / `.rsid.meta` endings (those three suffixes live in the constant
`_STORE_SUFFIXES`).

### The normalizer must match how the store was built

`DbSNPFlatQuery` takes a `normalizer` (default **`trim_alleles`**, imported from
the build module). This **must** be the same normalizer the store was built
with, or `variant → rsID` lookups will silently miss. Pass `normalizer=None`
only if the store was built with `--no-normalise`. (This is the "normalization
seam": to get full reference-aware left-alignment later, plug VariantValidator's
own normalizer in here *and* in the builder.)

---

## Step 2 : Direction 1: variant → rsID (`variant_to_rsids`)

**Question:** *"Does `1-169549811-C-T` have an rsID in GRCh38?"*

```python
hit = DBSNP.variant_to_rsids("grch38", "1", 169549811, "C", "T")

hit.found      # True
hit.first      # 6025          (handy for the usual 1-answer case)
hit.rsids      # (6025,)       (all matches, as a tuple)
hit.labels     # ('rs6025',)   (same, with the 'rs' prefix)
hit.normalized # VariantQuery(... "1", 169549811, "C", "T")
```

What happens inside **`DbSNPFlatQuery.variant_to_rsids`**:

1. `normalise_assembly("grch38")` picks the right store; the chromosome is
   passed through as-is; `pos` is forced to `int`; **`_require_allele`** checks
   REF and ALT are non-empty.
2. **`_handle(assembly)`** hands back this thread's `DbSNPFlat` (see Step 4).
3. It records what the query looks like **after normalization** (the
   `.normalized` field) — useful for logging/debugging.
4. It calls the raw reader's `variant_to_rsids`, which normalizes the query the
   same way and matches exactly on REF/ALT.
5. It wraps the answer in a **`VariantToRsidResult`**.

**Empty vs. error:** a variant that simply isn't in dbSNP returns a result with
`found == False` (not an exception). A genuinely bad input — an empty allele —
raises `ValueError`. This distinction lets backend code treat "not found" as a
normal outcome.

There is also a shortcut, **`variant_to_rsid`** (singular), returning just the
first rsID or `None`.

---

## Step 3 — Direction 2: rsID → variant (`rsid_to_variants`)

**Question:** *"Where is `rs334` in GRCh38?"* (a multi-allelic site)

```python
res = DBSNP.rsid_to_variants("grch38", "rs334")

res.found                 # True
res.to_strings()          # ['11-5227002-T-A,C,G']
res.to_vcf_tuples()       # [('11',5227002,'T','A'), ('11',5227002,'T','C'), ('11',5227002,'T','G')]
for p in res.placements:  # each p is a VariantPlacement
    p.chrom, p.pos, p.ref, p.alts   # '11', 5227002, 'T', ('A','C','G')
```

What happens inside **`DbSNPFlatQuery.rsid_to_variants`**:

1. **`parse_rsid`** turns `6025` / `"rs6025"` / `"6025"` into the number `6025`.
   A malformed id (`"abc"`, `0`, `-5`, `True`) raises `ValueError`; a valid id
   that isn't present just yields an empty result.
2. **`_handle(assembly)`** gets this thread's reader.
3. It calls the raw reader's `rsid_to_variants` and wraps each returned group in
   a **`VariantPlacement`**, collected into a **`RsidToVariantResult`**.

A `VariantPlacement` is "one place one rsID sits": a shared REF with one or more
ALTs. Its helpers:

- `VariantPlacement.to_string()` → `"11-5227002-T-A,C,G"`
- `VariantPlacement.to_vcf_tuples()` → one `(chrom, pos, ref, alt)` per ALT.

If an rsID maps to **more than one location** (it happens), `placements` has more
than one entry, and `RsidToVariantResult.to_vcf_tuples()` flattens them all.

---

## Step 4 : Why it's thread-safe (`_handle` + `threading.local`)

The raw `DbSNPFlat` keeps live file handles (BGZF/tabix) and is **not**
thread-safe ; two threads sharing one would corrupt each other's reads.

**`DbSNPFlatQuery._handle`** solves this with **`threading.local`**: each thread
gets **its own `DbSNPFlat` per assembly**, opened lazily the first time that
thread asks for it, then reused. This matches the reader's own advice ("open one
per thread/process") and fits a web backend with a fixed pool of worker threads:
each worker opens its handles once and reuses them for the life of the process.

```
Thread A ─┐                          ┌─ its own DbSNPFlat(grch38)
          ├─ DBSNP.variant_to_rsids ─┤
Thread B ─┘                          └─ a separate DbSNPFlat(grch38)
```

The **first** query on a fresh thread pays a small one-time "open the files"
cost; every query after that is fast.

---

## Step 5 : Cleanup (`close`, context manager)

**`DbSNPFlatQuery.close`** closes the **calling thread's** handles. Handles held
by other live worker threads are released when the process exits. You can also
use it as a context manager for short-lived scripts:

```python
with DbSNPFlatQuery.from_directory("/data/dbsnp_flat", "b157") as q:
    print(q.rsid_to_variants("grch38", "rs6025").to_strings())
# handles for this thread are closed on exit
```

Small extras: `DBSNP.assemblies` lists what it can serve (e.g. `('grch38',)`),
and `DBSNP.has_assembly("grch37")` checks availability. Asking for an assembly
that wasn't configured raises `KeyError` with a helpful message.

---

## The result objects at a glance

| Class | Key fields / methods | Represents |
|---|---|---|
| `VariantToRsidResult` | `.found`, `.first`, `.rsids`, `.labels`, `.normalized` | answer to variant → rsID |
| `RsidToVariantResult` | `.found`, `.placements`, `.to_strings()`, `.to_vcf_tuples()` | answer to rsID → variant |
| `VariantPlacement` | `.chrom`, `.pos`, `.ref`, `.alts`, `.to_vcf_tuples()` | one locus of one rsID |
| `VariantQuery` | `.chrom`, `.pos`, `.ref`, `.alt`, `.to_string()` | the (normalized) input, echoed back |

All are frozen dataclasses 

---

## End-to-end example

```python
from vvDBFlatQueryDbSNP import DbSNPFlatQuery

DBSNP = DbSNPFlatQuery.from_directory("/data/dbsnp_flat", version="b157")

# variant → rsID
h = DBSNP.variant_to_rsids("grch38", "1", 169549811, "C", "T")
if h.found:
    print("rsID:", h.first)                 # rsID: 6025

# rsID → variant
r = DBSNP.rsid_to_variants("grch38", "rs6025")
for chrom, pos, ref, alt in r.to_vcf_tuples():
    print(chrom, pos, ref, alt)             # 1 169549811 C T  (and any other ALTs)
```

---
