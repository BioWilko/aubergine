"""Per-sample event index and local haplotype reconstruction.

`apply` walks a reference window and emits the sample's sequence over it, with the origin
of every base. It is used for disagreement analysis, VCF alleles, CDS projection and the
round-trip test.
"""

from __future__ import annotations

from bisect import bisect_left, bisect_right
from collections.abc import Callable, Iterable
from dataclasses import dataclass, field

import numpy as np

from .call import Ambiguous, Call, CallSet, Variant
from .fasta import ACGT

REF, SUB, INS, MISS = "ref", "sub", "ins", "missing"


@dataclass(frozen=True)
class Base:
    base: str
    kind: str  # ref | sub | ins | missing
    pos: int  # reference position (for inserted bases: the anchor p)

    @property
    def from_sample(self) -> bool:
        return self.kind in (SUB, INS)


@dataclass
class SampleEvents:
    """One sample's events on one contig, indexed for state and haplotype queries."""

    sample: str
    contig: str
    length: int
    keys: set[str] = field(default_factory=set)
    snps: dict[int, tuple[str, str]] = field(default_factory=dict)  # pos -> (alt, key)
    ambiguous: dict[int, str] = field(default_factory=dict)  # pos -> code
    deletions: list[tuple[int, int, str]] = field(default_factory=list)
    insertions: dict[int, tuple[str, str]] = field(default_factory=dict)  # p -> (raw, key)
    missing: list[tuple[int, int]] = field(default_factory=list)
    missing_insertions: dict[int, int] = field(default_factory=dict)  # p -> length
    covered: tuple[int, int] | None = None

    def __post_init__(self):
        self.reindex()

    def reindex(self) -> None:
        L = self.length
        miss = np.zeros(L + 2, dtype=np.int64)
        for s, e in self.missing:
            miss[s : e + 1] = 1
        dele = np.zeros(L + 2, dtype=np.int64)
        for s, e, _ in self.deletions:
            dele[s : e + 1] = 1
        self._missing = miss.astype(bool)
        self._deleted = dele.astype(bool)
        self._miss_cum = np.cumsum(miss)
        self._del_cum = np.cumsum(dele)
        self._snp_pos = sorted(self.snps)
        self._amb_pos = sorted(self.ambiguous)

    @staticmethod
    def _range_any(cum: np.ndarray, lo: int, hi: int) -> bool:
        if hi < lo:
            return False
        return bool(cum[hi] - cum[lo - 1] > 0)

    def missing_in(self, lo: int, hi: int) -> bool:
        return self._range_any(self._miss_cum, lo, hi)

    def deleted_in(self, lo: int, hi: int) -> bool:
        return self._range_any(self._del_cum, lo, hi)

    def is_missing(self, q: int) -> bool:
        return bool(self._missing[q])

    def is_deleted(self, q: int) -> bool:
        return bool(self._deleted[q])

    def ambiguous_in(self, lo: int, hi: int) -> bool:
        return bisect_right(self._amb_pos, hi) > bisect_left(self._amb_pos, lo)

    def snps_in(self, lo: int, hi: int) -> list[int]:
        return self._snp_pos[bisect_left(self._snp_pos, lo) : bisect_right(self._snp_pos, hi)]

    def add_variant(self, v: Variant, raw: str = "") -> None:
        """Add one called variant (no reindex; call reindex() after a batch)."""
        self.keys.add(v.key)
        if v.type == "snp":
            self.snps[v.start] = (v.alt, v.key)
        elif v.type == "del":
            self.deletions.append((v.start, v.end, v.key))
        else:
            self.insertions[v.start] = (raw or v.ins_seq, v.key)


def build_sample_events(cs: CallSet) -> dict[tuple[str, str], SampleEvents]:
    """Index every (sample, contig) of a call set."""
    out: dict[tuple[str, str], SampleEvents] = {}
    for s in cs.samples:
        for contig in cs.reference.ids:
            se = SampleEvents(s.id, contig, cs.reference.length(contig))
            se.missing = [(a, b) for c, a, b in s.missing if c == contig]
            se.missing_insertions = {p: n for c, p, n in s.missing_insertions if c == contig}
            cov = [(a, b) for c, a, b in s.covered if c == contig]
            se.covered = cov[0] if cov else None
            out[(s.id, contig)] = se
    for c in cs.calls:
        v = cs.variants[c.key]
        out[(c.sample, v.contig)].add_variant(v, c.ins_seq_raw)
    for a in cs.ambiguous:
        out[(a.sample, a.contig)].ambiguous[a.pos] = a.code
    for se in out.values():
        se.deletions.sort()
        se.reindex()
    return out


def subset_events(
    template: SampleEvents,
    variants: Iterable[tuple[Variant, str]],
    ambiguous: Iterable[Ambiguous] = (),
) -> SampleEvents:
    """A SampleEvents containing only the given (variant, raw insertion) pairs."""
    se = SampleEvents(template.sample, template.contig, template.length)
    for v, raw in variants:
        se.add_variant(v, raw)
    for a in ambiguous:
        se.ambiguous[a.pos] = a.code
    se.deletions.sort()
    se.reindex()
    return se


def calls_as_events(calls: Iterable[Call], cs: CallSet) -> list[tuple[Variant, str]]:
    return [(cs.variants[c.key], c.ins_seq_raw) for c in calls]


def apply(
    ref_seq: str,
    se: SampleEvents,
    start: int,
    end: int,
    missing_as_n: bool = False,
    insertion_ok: Callable[[int], bool] | None = None,
) -> list[Base]:
    """Emit the sample's bases over reference window [start, end] (1-based, inclusive).

    An insertion at p is emitted directly after position p, if start <= p <= end, or
    before position 1 if p == 0 and start == 1. With `insertion_ok`, only insertion points
    it accepts are emitted, and missing insertions are never emitted.
    """
    out: list[Base] = []

    def insertion(p: int):
        if insertion_ok is not None and not insertion_ok(p):
            return
        if p in se.insertions:
            raw = se.insertions[p][0]
            out.extend(Base(b, INS, p) for b in raw)
        elif missing_as_n and insertion_ok is None and p in se.missing_insertions:
            out.extend(Base("N", MISS, p) for _ in range(se.missing_insertions[p]))

    if start == 1:
        insertion(0)
    for q in range(start, end + 1):
        if se.is_deleted(q):
            pass
        elif missing_as_n and se.is_missing(q):
            out.append(Base("N", MISS, q))
        elif q in se.snps:
            out.append(Base(se.snps[q][0], SUB, q))
        elif q in se.ambiguous:
            out.append(Base(se.ambiguous[q], SUB, q))
        else:
            out.append(Base(ref_seq[q - 1], REF, q))
        insertion(q)
    return out


def as_string(bases: list[Base]) -> str:
    return "".join(b.base for b in bases)


def window_haplotype(ref_seq: str, se: SampleEvents, start: int, end: int) -> tuple[str, str]:
    """Classify a sample's haplotype over a locus window.

    Returns (status, sequence) with status `missing`, `ambiguous` or `callable`.
    """
    if se.missing_in(start, end):
        return "missing", ""
    for p in se.missing_insertions:
        if start - 1 <= p <= end and (p + 1 <= end or p == 0):
            return "missing", ""
    bases = apply(ref_seq, se, start, end)
    seq = as_string(bases)
    if any(b.from_sample and b.base not in ACGT for b in bases):
        return "ambiguous", seq
    return "callable", seq
