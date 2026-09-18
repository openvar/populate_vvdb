# Test suite for the flat-dbSNP tools

Covers the three modules:

| Test file | Module under test | What it checks |
|---|---|---|
| `test_build.py` | `vvDBFlatBuildDbSNP.py` | `trim_alleles`, contig maps, `_strip_rs`/`_detect_source_type`, the on-disk format (`pack_entry`/`iter_entries`, `MetaWriter`/`MetaIndex`), the builder, the `DbSNPFlat` reader (both directions), and the CLI |
| `test_query.py` | `vvDBFlatQueryDbSNP.py` | input parsing, result objects, `DbSNPFlatQuery` construction/errors, both direction methods, thread-safety, lifecycle |
| `test_validate.py` | `vvDBFlatValidate.py` | helpers, source/store sampling, accuracy/format/curated checks, benchmark, `validate_assembly`, reporting, and the CLI |

## How it works

`conftest.py` builds a **synthetic store once per session**: it writes a small
NCBI-style source VCF (RefSeq-accession contigs, no `##contig` header lines,
bgzipped + tabix-indexed) with known landmark variants plus filler, then runs the
real builder on it. Every test queries that store, so the suite needs no network
and no 30 GB download.

## Requirements

- Python 3.9+
- `pytest`
- `pysam` (already required by the tools themselves)

```bash
pip install pytest pysam
```

## Running

From the project root (`Ultimate_ME_populate_vvdb/`):

```bash
# everything
python3 -m pytest

# one module's tests
python3 -m pytest tests/test_query.py

# a single test, verbose
python3 -m pytest tests/test_build.py::test_trim_alleles -v

# quieter, with a short summary of any failures
python3 -m pytest -q
```

`pytest.ini` sets `pythonpath = .` so the modules import without installation,
and `testpaths = tests` so a bare `pytest` finds them.
