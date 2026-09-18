# `vvDBFlatBuildDbSNP.py` 

## What this file is for 

dbSNP is a giant catalogue of human genetic variants. Each variant has an
**rsID** (e.g. `rs6025`) and a **position + alleles** (e.g. chromosome `1`,
position `169549811`, reference `C`, alternate `T`). People need to convert
between the two directions:

- **rsID → variant**  ("what/where is `rs6025`?")
- **variant → rsID**  ("does `1-169549811-C-T` have an rsID?")

The official dbSNP file for one genome build is a ~30 GB VCF. Searching that
directly is slow, and running a database server is heavy. This file builds a
small, fast, **flat-file store** that answers both questions quickly, and it
also contains the **reader** used to query that store.

It produces **two artefacts per genome build (assembly)**:

| Artefact | Files | Answers |
|---|---|---|
| **MASTER** | `*.variants.tsv.gz` + `.tbi` | variant → rsID (and holds the actual rows) |
| **rsID INDEX** | `*.rsid.idx` + `*.rsid.meta` | rsID → variant |

So this one file plays **three roles**: it defines the on-disk **format**, it
provides the **reader** (`DbSNPFlat`), and it provides the **builder**
(`FlatDbSNPBuilder`) that creates the files.

---

## Map of the file (what lives where)

| Part | Functions / classes | Job |
|---|---|---|
| Contig maps | `get_contig_map`, `GRCH37`, `GRCH38` | RefSeq accession (`NC_000001.11`) → plain name (`1`) |
| On-disk format | `MAGIC`, `pack_entry`, `iter_entries`, `BlockRecord`, `MetaWriter`, `MetaIndex` | the byte layout of the rsID index |
| Normalization | `trim_alleles` | make equivalent allele spellings identical |
| Reader | `VariantGroup`, `DbSNPFlat` | query the built store, both directions |
| Builder | `FlatDbSNPBuilder`, `_strip_rs`, `_detect_source_type` | turn a source VCF into the store |
| Command line | `_main` | run the builder from a terminal |

---

## Part A: Building the store, step by step

Everything below is orchestrated by **`FlatDbSNPBuilder.build()`**. We'll follow
a tiny example all the way through.

### Our tiny example input

Imagine a mini dbSNP VCF for GRCh38 with three records (chromosomes are written
as RefSeq **accessions**, exactly like the real dbSNP file):

```
NC_000001.11   169549811   rs6025   C   T
NC_000011.10   5227002      rs334    T   A,C,G
NC_000001.11   10000        rs99     AA  GA
```

- `rs6025` : a simple SNV (one base → one base).
- `rs334` : **multi-allelic**: one position, three alternates.
- `rs99` : a **non-minimal** spelling (`AA→GA`) that should really be `A→G`.

### Step 1 : Detect the input type

**`_detect_source_type(path)`** looks at the file extension and returns
`"bcf"`, `"vcf"`, or `"tsv"`. `build()` uses this to pick the right reader.

### Step 2 : Read the source into a stream of variants

Depending on the type, one of these runs:

- **`_iter_bcf`** : decodes a `.bcf`/`.vcf` using the `bcftools` program.
- **`_iter_vcf`** : reads a VCF/BCF directly with `pysam` (no external program).
- **`_iter_tsv`** : reads a simple tab-separated test file.

**`_iter_source`** picks one of them and wraps the result. Each reader does two
important things per record:

1. **Extract the rsID as a number** with **`_strip_rs`**: `"rs6025" → 6025`.
   Records without a usable rsID are dropped.
2. **Rename the chromosome** using the contig map (Step 2a below).
3. **Split multi-allelic records** : `rs334`'s `A,C,G` becomes three separate
   `(chrom, pos, rsid, ref, alt)` rows.

After this step our stream is:

```
("NC_000001.11" → "1", 169549811, 6025, "C", "T")
("NC_000011.10" → "11", 5227002, 334, "T", "A")
("NC_000011.10" → "11", 5227002, 334, "T", "C")
("NC_000011.10" → "11", 5227002, 334, "T", "G")
("NC_000001.11" → "1", 10000, 99, "AA", "GA")
```

#### Step 2a : Chromosome renaming (`get_contig_map`)

dbSNP names chromosomes as accessions like `NC_000001.11`. VariantValidator
wants plain names like `1`. **`get_contig_map("grch38")`** returns the lookup
`{"NC_000001.11": "1", ...}`. Note the accessions **differ between builds**
(`NC_000001.10` is chr1 in GRCh37, `NC_000001.11` in GRCh38) — that's why there
are two maps, `GRCH37` and `GRCH38`.

### Step 3 : Normalize the alleles (`trim_alleles`)

Still inside `_iter_source`, every variant passes through the **normalizer**
(default **`trim_alleles`**). This trims shared bases so that equivalent
spellings become identical:

```
trim_alleles("1", 10000, "AA", "GA")  →  ("1", 10000, "A", "G")
```

It removed the shared trailing `A`. (SNVs like `C→T` are already minimal and
pass through unchanged.) This matters because the **same normalizer must be used
when querying**, so both sides agree. `trim_alleles` is reference-*free*: it does
not do full left-alignment (that needs the genome sequence).

### Step 4 : Write the sorted MASTER file (`_write_master`)

**`_write_master`** writes every variant as a row of **six columns** and sorts
them. The columns are:

```
CHROM   POS   END   rsID   REF   ALT
```

`END = POS + len(REF) - 1`. For our example the sorted master (uncompressed) is:

```
1    10000       10000       99     A   G
1    169549811   169549811   6025   C   T
11   5227002     5227002     334    T   A
11   5227002     5227002     334    T   C
11   5227002     5227002     334    T   G
```

The sort (by chromosome, then position, then rsID, then alleles) guarantees that
**all alternates of one variant sit next to each other** — needed later to
re-assemble multi-allelic sites.

### Step 5 : Compress + index the master (`_index_master`)

**`_index_master`** bgzips the file to `*.variants.tsv.gz` and builds a
**tabix** index (`.tbi`). Tabix lets you jump straight to "chromosome 11 around
position 5227002" without scanning the whole file. This is what makes the
**variant → rsID** direction fast.

### Step 6 : Collect (rsID → location) pairs (`_collect_pairs`)

Now the builder prepares the **rsID → variant** direction. **`_collect_pairs`**
walks the compressed master and, for each distinct variant group, records:

```
rsID   →   voffset
```

`voffset` is the **BGZF virtual offset** : a bookmark pointing at that group's
first row inside the compressed master. For our example:

```
99     →  <offset of the rs99 row>
6025   →  <offset of the rs6025 row>
334    →  <offset of the first rs334 row>
```

(One pair per variant, even though `rs334` has three alternate rows ; they share
one starting offset.)

### Step 7 : Sort the pairs **by rsID** (`_sort_pairs`) ← the critical step

**`_sort_pairs`** sorts those pairs by the **rsID number**. Our three rsIDs
`6025`, `334`, `99` become, in numeric order:

```
99     →  <offset>
334    →  <offset>
6025   →  <offset>
```

**Why this matters:** the index is later searched with **binary search** on the
rsID. Binary search only works on sorted data. An earlier prototype left the
pairs in genome-position order — which is effectively random rsIDs — so lookups
silently fell back to scanning everything. This sort is the fix. (The module
docstring/`build()` call this out explicitly.)

### Step 8 : Write the compressed index + signposts (`_write_index`)

**`_write_index`** turns the sorted pairs into the two index files:

- **`*.rsid.idx`** : the pairs packed into fixed 16-byte records
  (`pack_entry` writes each `(rsid, voffset)`), grouped into **blocks** of
  `block_size` entries (default 4096), and each block **zlib-compressed**.
- **`*.rsid.meta`** : a tiny **signpost table** built by **`MetaWriter`**. For
  every block it stores `(first_rsid, where the block is, how big, how many
  entries)` via `MetaWriter.add_block`, then `MetaWriter.write`.

Think of `.meta` as the **table of contents** and `.idx` as the **chapters**.
To find an rsID you read the small table of contents, jump to the one block that
could contain it, and decompress only that block.

`build()` also **skips work if all four output files already exist** (it is
idempotent), and it lays them out as
`<out>/<version>/dbsnp_<assembly>_<version>.*`.

---

## Part B : Reading the store (how a lookup actually works)

The reader is **`DbSNPFlat`**. You open it with the three store files:

```python
db = DbSNPFlat("….variants.tsv.gz", "….rsid.idx", "….rsid.meta")
```

Its constructor opens: the `.meta` signpost table (`MetaIndex.load`), the
compressed master for random reads (`pysam.BGZFile`), the tabix handle
(`pysam.TabixFile`), and the raw `.idx` file.

### Direction 1 : variant → rsID (`variant_to_rsids`)

Example: *"what rsID is `11-5227002-T-C`?"*

1. `variant_to_rsids` runs the query through the **same normalizer** (`trim_alleles`).
2. It uses the **tabix** index to fetch only the master rows near
   `11:5227002`.
3. It keeps rows whose POS, REF and ALT match exactly → returns `[334]`.

Because the store is normalized, a non-minimal query also works: asking for
`1-10000-AA-GA` gets trimmed to `1-10000-A-G` and finds `rs99`. There is also
`variant_to_rsid` (singular) that returns just the first hit.

### Direction 2 : rsID → variant (`rsid_to_variants` / `resolve_rsid`)

Example: *"where is `rs334`?"*

1. **`MetaIndex.candidate_block`** binary-searches the signpost table for the
   block whose range covers `334`.
2. **`_offsets_for_rsid`** decompresses that block (`_read_block`), finds the
   `voffset(s)` stored for `334`, and (for safety) checks neighbouring blocks in
   case duplicates straddle a boundary.
3. **`_read_group_at`** seeks the master to that offset and reads the contiguous
   rows for the variant, re-assembling the multi-allelic site into one
   **`VariantGroup`**: `rsid=334, chrom="11", pos=5227002, ref="T",
   alts=("A","C","G")`.
4. `resolve_rsid` is a convenience that returns the group as a string via
   `VariantGroup.to_string()` → `"11-5227002-T-A,C,G"`.

An rsID that maps to more than one place returns more than one group.

---

## Part C : The command line (`_main`)

**`_main`** is the terminal entry point. It parses flags (`--grch37`,
`--grch38`, `--out`, `--version`, `--source-type`, `--block-size`,
`--compress-level`, `--no-normalise`, `--no-contig-map`) and calls
`FlatDbSNPBuilder.build()` for each assembly. Example:

```bash
python vvDBFlatBuildDbSNP.py --out dbsnp_flat --version b157 \
    --grch38 /data/GCF_000001405.40.gz
```

That single command runs Steps 1–8 above and writes the four store files under
`dbsnp_flat/b157/`.

---