#!/usr/bin/env python3
"""
vvDBFlatQueryDbSNP.py

Backend query API for the flat dbSNP store (built by ``vvDBFlatBuildDbSNP.py``).

This is a **library, not a CLI**. It exposes the two lookup directions as two
separate calls, meant to be imported into the VariantValidator backend and
called per request, per variant:

    variant_to_rsids(assembly, chrom, pos, ref, alt)   sequence -> rsID
    rsid_to_variants(assembly, rsid)                   rsID     -> sequence

The query logic mirrors ``vvDBFlatValidate.py`` (both wrap :class:`DbSNPFlat`
from ``vvDBFlatBuildDbSNP.py``). On top of that raw reader this module adds the
things a backend integration needs and a validation script does not:

  * per-assembly store management (GRCh37 and/or GRCh38), opened lazily;
  * input parsing/validation (rsID as ``6025`` / ``"rs6025"`` / ``"6025"``);
  * structured, JSON-friendly return types instead of printed strings;
  * thread-safety: :class:`DbSNPFlat` is stateful and NOT thread-safe, so each
    thread gets its own handle (its own guidance: "open one per thread/process").

Normalization seam
------------------
``DbSNPFlat`` runs every ``sequence -> rsID`` query through the SAME normalizer
the store was built with, then matches exactly. The default here is
:func:`trim_alleles` (the builder's default too). It is reference-*free*
minimal-representation trimming, NOT reference-aware left-alignment, so a
coordinate query that needs left-alignment (an indel in a repeat/homopolymer)
can miss. To close that gap, build the store AND construct this client with
VariantValidator's own reference-aware normalizer - the only rule is that the
two match. Pass ``normalizer=None`` only if the store was built with
``--no-normalise``.

Integration sketch
------------------
    from vvDBFlatQueryDbSNP import DbSNPFlatQuery

    # once, at process/app startup (paths from config):
    DBSNP = DbSNPFlatQuery.from_directory("/data/dbsnp_flat", version="b157")
    # or:  DbSNPFlatQuery({"grch38": "/data/.../dbsnp_grch38_b157", ...})

    # then, per request:
    hit = DBSNP.variant_to_rsids("grch38", "1", 169549811, "C", "T")
    if hit.found:
        rs = hit.first                      # -> 6025

    res = DBSNP.rsid_to_variants("grch38", "rs6025")
    for p in res.placements:
        chrom, pos, ref, alt = p.to_vcf_tuples()[0]
"""

from __future__ import annotations

import os
import threading
from dataclasses import dataclass, field
from typing import Callable, Optional

from vvDBFlatBuildDbSNP import DbSNPFlat, trim_alleles

# The three files that make up one assembly's store, appended to a path *stem*.
_STORE_SUFFIXES = (".variants.tsv.gz", ".rsid.idx", ".rsid.meta")

# A normalizer is (chrom, pos, ref, alt) -> (chrom, pos, ref, alt); or None.
Normalizer = Optional[Callable[[str, int, str, str], tuple]]


# ---------------------------------------------------------------------------
# input helpers
# ---------------------------------------------------------------------------
def normalise_assembly(label) -> str:
    """
    Canonicalise an assembly label to ``"grch37"`` or ``"grch38"``.

    Accepts the same spellings as ``vvDBFlatBuildDbSNP.get_contig_map``
    ('GRCh38', 'hg38', '38', 'GRCh37', 'hg19', '37', ...). Raises ``ValueError``
    for anything else (a backend should fail loudly on an unknown assembly).
    """
    a = str(label).strip().lower()
    if a in ("grch38", "hg38", "38"):
        return "grch38"
    if a in ("grch37", "hg19", "37"):
        return "grch37"
    raise ValueError(f"unknown assembly {label!r}; expected grch37 or grch38")


def parse_rsid(rsid) -> int:
    """
    Parse an rsID given as ``6025`` (int), ``"rs6025"`` or ``"6025"`` -> ``6025``.

    Raises ``ValueError`` on a malformed id (that is caller error, distinct from
    a well-formed id that simply isn't in the store, which yields an empty
    result rather than an exception).
    """
    if isinstance(rsid, bool):  # bool is an int subclass; reject it explicitly
        raise ValueError(f"invalid rsID {rsid!r}")
    if isinstance(rsid, int):
        if rsid <= 0:
            raise ValueError(f"invalid rsID {rsid!r}")
        return rsid
    s = str(rsid).strip()
    if s.lower().startswith("rs"):
        s = s[2:]
    try:
        val = int(s)
    except (TypeError, ValueError):
        raise ValueError(f"invalid rsID {rsid!r}; expected e.g. 6025 or 'rs6025'")
    if val <= 0:
        raise ValueError(f"invalid rsID {rsid!r}")
    return val


def _require_allele(value, which: str) -> str:
    """Validate a REF/ALT allele string is present and non-empty."""
    s = str(value).strip()
    if not s:
        raise ValueError(f"empty {which} allele")
    return s


def store_stem(base_dir: str, assembly: str, version: str) -> str:
    """
    Path stem for one assembly's store under the builder's layout
    (``<base_dir>/<version>/dbsnp_<assembly>_<version>``). The three store files
    are this stem + :data:`_STORE_SUFFIXES`.
    """
    asm = normalise_assembly(assembly)
    return os.path.join(base_dir, version, f"dbsnp_{asm}_{version}")


# ---------------------------------------------------------------------------
# result types  (structured; safe to serialise / log / return from an endpoint)
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class VariantQuery:
    """The (already assembly-scoped) variant that was looked up."""

    assembly: str
    chrom: str
    pos: int
    ref: str
    alt: str

    def to_string(self) -> str:
        return f"{self.chrom}-{self.pos}-{self.ref}-{self.alt}"


@dataclass(frozen=True)
class VariantToRsidResult:
    """
    Result of a ``sequence -> rsID`` lookup.

    ``rsids`` are the dbSNP ids whose stored REF/ALT match the (normalized)
    query at that locus - usually one, occasionally more.
    """

    query: VariantQuery
    normalized: Optional[VariantQuery]   # query after the store's normalizer
    rsids: tuple[int, ...]

    @property
    def found(self) -> bool:
        return bool(self.rsids)

    @property
    def labels(self) -> tuple[str, ...]:
        return tuple(f"rs{r}" for r in self.rsids)

    @property
    def first(self) -> Optional[int]:
        """The first matching rsID, or ``None`` - convenience for the 1:1 case."""
        return self.rsids[0] if self.rsids else None


@dataclass(frozen=True)
class VariantPlacement:
    """One locus of one rsID: a shared REF with one or more ALT alleles."""

    rsid: int
    chrom: str
    pos: int
    ref: str
    alts: tuple[str, ...]

    @property
    def label(self) -> str:
        return f"rs{self.rsid}"

    def to_string(self) -> str:
        """``chrom-pos-ref-alt,alt`` (multi-allelic ALTs comma-joined)."""
        return f"{self.chrom}-{self.pos}-{self.ref}-{','.join(self.alts)}"

    def to_vcf_tuples(self) -> list[tuple[str, int, str, str]]:
        """Expand to one ``(chrom, pos, ref, alt)`` per ALT allele."""
        return [(self.chrom, self.pos, self.ref, a) for a in self.alts]


@dataclass(frozen=True)
class RsidToVariantResult:
    """
    Result of an ``rsID -> sequence`` lookup.

    ``placements`` is usually a single locus; a multi-mapping rsID yields more.
    Empty when the rsID is not in this assembly's store.
    """

    assembly: str
    rsid: int
    placements: tuple[VariantPlacement, ...]

    @property
    def label(self) -> str:
        return f"rs{self.rsid}"

    @property
    def found(self) -> bool:
        return bool(self.placements)

    def to_strings(self) -> list[str]:
        return [p.to_string() for p in self.placements]

    def to_vcf_tuples(self) -> list[tuple[str, int, str, str]]:
        """All ``(chrom, pos, ref, alt)`` tuples across every placement/ALT."""
        out: list[tuple[str, int, str, str]] = []
        for p in self.placements:
            out.extend(p.to_vcf_tuples())
        return out


# ---------------------------------------------------------------------------
# the client
# ---------------------------------------------------------------------------
class DbSNPFlatQuery:
    """
    Bidirectional query client over one or both assemblies' flat dbSNP stores.

    Construct once (at app/process startup) and reuse for every request. The two
    directions are separate methods:

        variant_to_rsids(assembly, chrom, pos, ref, alt) -> VariantToRsidResult
        rsid_to_variants(assembly, rsid)                 -> RsidToVariantResult

    Thread-safety
    -------------
    ``DbSNPFlat`` holds stateful BGZF/tabix handles and is NOT thread-safe, so
    this client keeps **one handle per thread per assembly** (opened lazily on
    first use in each thread). That fits a request-serving backend with a
    bounded worker pool: each worker opens its handles once and reuses them.
    :meth:`close` closes the *calling* thread's handles; handles held by other
    live threads are released when the process exits. (If you prefer a single
    shared handle instead, guard every call with your own lock and open one
    ``DbSNPFlat`` yourself - but per-thread handles avoid that contention.)
    """

    def __init__(self, stores: dict, normalizer: Normalizer = trim_alleles):
        """
        ``stores`` maps assembly label -> store path *stem* (the path without the
        ``.variants.tsv.gz`` / ``.rsid.idx`` / ``.rsid.meta`` suffixes). At least
        one assembly is required; the files must exist (checked now, fail fast).

        ``normalizer`` MUST match how the store was built (default
        :func:`trim_alleles`; ``None`` if built with ``--no-normalise``).
        """
        if not stores:
            raise ValueError("no stores given; expected at least one assembly")
        self._stems: dict[str, str] = {}
        for label, stem in stores.items():
            asm = normalise_assembly(label)
            for suffix in _STORE_SUFFIXES:
                path = stem + suffix
                if not os.path.exists(path):
                    raise FileNotFoundError(
                        f"{asm}: missing store file {path}")
            self._stems[asm] = stem
        self._normalizer = normalizer
        self._local = threading.local()

    # -- construction helpers -------------------------------------------
    @classmethod
    def from_directory(cls, base_dir: str, version: str,
                       assemblies=("grch37", "grch38"),
                       normalizer: Normalizer = trim_alleles,
                       require_all: bool = False) -> "DbSNPFlatQuery":
        """
        Build a client from the builder's on-disk layout
        (``<base_dir>/<version>/dbsnp_<assembly>_<version>``). Assemblies whose
        files are absent are skipped unless ``require_all`` is set.
        """
        stores: dict[str, str] = {}
        for asm in assemblies:
            stem = store_stem(base_dir, asm, version)
            present = all(os.path.exists(stem + s) for s in _STORE_SUFFIXES)
            if present:
                stores[normalise_assembly(asm)] = stem
            elif require_all:
                raise FileNotFoundError(
                    f"{asm}: store not found at stem {stem}")
        if not stores:
            raise FileNotFoundError(
                f"no stores found under {base_dir} for version {version}")
        return cls(stores, normalizer=normalizer)

    # -- introspection --------------------------------------------------
    @property
    def assemblies(self) -> tuple[str, ...]:
        """Assemblies this client can serve (canonical labels)."""
        return tuple(sorted(self._stems))

    def has_assembly(self, assembly) -> bool:
        try:
            return normalise_assembly(assembly) in self._stems
        except ValueError:
            return False

    # -- per-thread handle management -----------------------------------
    def _handle(self, assembly: str) -> DbSNPFlat:
        asm = normalise_assembly(assembly)
        stem = self._stems.get(asm)
        if stem is None:
            raise KeyError(
                f"no store configured for {asm}; have {self.assemblies}")
        cache = getattr(self._local, "handles", None)
        if cache is None:
            cache = self._local.handles = {}
        db = cache.get(asm)
        if db is None:
            db = DbSNPFlat(stem + ".variants.tsv.gz",
                           stem + ".rsid.idx",
                           stem + ".rsid.meta",
                           normalizer=self._normalizer)
            cache[asm] = db
        return db

    # -- Direction 1:  sequence -> rsID ---------------------------------
    def variant_to_rsids(self, assembly, chrom, pos, ref, alt
                         ) -> VariantToRsidResult:
        """
        Map an exact ``chrom, pos, ref, alt`` to its dbSNP rsID(s) for one
        assembly. The query is run through the store's normalizer, then matched
        exactly on REF and ALT (this is the same call the validator checks).

        ``chrom`` is a plain chromosome name as stored in the master file
        (``1``..``22``, ``X``, ``Y``, ``MT``); ``pos`` is 1-based. Returns a
        :class:`VariantToRsidResult` (``.found`` / ``.first`` / ``.rsids`` /
        ``.labels``); the result is empty (not an error) when nothing matches.
        """
        asm = normalise_assembly(assembly)
        chrom = str(chrom).strip()
        if not chrom:
            raise ValueError("empty chrom")
        pos = int(pos)
        ref = _require_allele(ref, "ref")
        alt = _require_allele(alt, "alt")

        db = self._handle(asm)
        query = VariantQuery(asm, chrom, pos, ref, alt)

        normalized: Optional[VariantQuery] = None
        if self._normalizer is not None:
            nc, npos, nref, nalt = self._normalizer(chrom, pos, ref, alt)
            normalized = VariantQuery(asm, nc, int(npos), nref, nalt)

        rsids = tuple(db.variant_to_rsids(chrom, pos, ref, alt))
        return VariantToRsidResult(query, normalized, rsids)

    def variant_to_rsid(self, assembly, chrom, pos, ref, alt) -> Optional[int]:
        """Convenience: the first matching rsID for the 1:1 case, or ``None``."""
        return self.variant_to_rsids(assembly, chrom, pos, ref, alt).first

    # -- Direction 2:  rsID -> sequence ---------------------------------
    def rsid_to_variants(self, assembly, rsid) -> RsidToVariantResult:
        """
        Resolve an rsID to its variant placement(s) for one assembly. ``rsid``
        may be an int, ``"rs6025"`` or ``"6025"``. Returns a
        :class:`RsidToVariantResult` (``.found`` / ``.placements`` /
        ``.to_vcf_tuples()`` / ``.to_strings()``); empty when the rsID is not in
        this assembly's store. Raises ``ValueError`` on a malformed rsID.
        """
        asm = normalise_assembly(assembly)
        rid = parse_rsid(rsid)
        db = self._handle(asm)
        placements = tuple(
            VariantPlacement(g.rsid, g.chrom, g.pos, g.ref, tuple(g.alts))
            for g in db.rsid_to_variants(rid))
        return RsidToVariantResult(asm, rid, placements)

    # -- lifecycle ------------------------------------------------------
    def close(self) -> None:
        """Close the calling thread's open handles (see thread-safety note)."""
        cache = getattr(self._local, "handles", None)
        if cache:
            for db in cache.values():
                try:
                    db.close()
                except Exception:
                    pass
            cache.clear()

    def __enter__(self) -> "DbSNPFlatQuery":
        return self

    def __exit__(self, *exc) -> None:
        self.close()


