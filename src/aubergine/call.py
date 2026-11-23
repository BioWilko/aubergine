"""Variant calling from a reference-anchored MSA.

Coordinates are 1-based on the ungapped reference. Internally an event is described by
`(type, start, end)`: an SNP at `start == end`, a deletion of bases `start..end`, or an
insertion between `start` (= p, 0..L) and `end` (= p + 1). VCF-style anchored POS/REF/ALT
and the variant key are derived from that.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from . import __version__
from .fasta import ACGT, GAP, ContigAlignment, N, Reference, SetAlignment

# Column classes
SKIP, MATCH, SNP, AMBIG, MISSING, DEL, INS, OVERHANG, UNCALLABLE = range(9)

# Per-reference-position status of a sample
PLAIN, S_SNP, S_AMBIG, S_MISSING, S_DELETED, S_UNCALLABLE = range(6)
_CLASS_TO_STATUS = np.zeros(9, dtype=np.int8)
_CLASS_TO_STATUS[SNP] = S_SNP
_CLASS_TO_STATUS[AMBIG] = S_AMBIG
_CLASS_TO_STATUS[MISSING] = S_MISSING
_CLASS_TO_STATUS[DEL] = S_DELETED
_CLASS_TO_STATUS[UNCALLABLE] = S_UNCALLABLE

_ACGT_CODES = np.frombuffer(ACGT.encode(), dtype=np.uint8)

TYPE_NAMES = {"snp": "snp", "ins": "insertion", "del": "deletion"}


class CallError(ValueError):
    pass


# --------------------------------------------------------------------------- keys


def make_key(vtype: str, contig: str, pos: int, ref: str, alt: str) -> str:
    return f"{vtype}:{contig}:{pos}:{ref}:{alt}"


def parse_key(key: str) -> tuple[str, str, int, str, str]:
    """Inverse of make_key; a contig ID containing ':' round-trips."""
    vtype, rest = key.split(":", 1)
    contig, pos, ref, alt = rest.rsplit(":", 3)
    return vtype, contig, int(pos), ref, alt


# --------------------------------------------------------------------------- model


@dataclass
class Variant:
    contig: str
    type: str  # snp | ins | del
    start: int
    end: int
    pos: int
    ref: str
    alt: str
    ins_seq: str = ""  # inserted sequence as written in ALT (IUPAC -> N)
    hp_len: int = 0
    hp_run: tuple[str, int, int] | None = None
    slide: tuple[int, int] = (0, 0)
    cds: list[dict] = field(default_factory=list)  # [{gene, delta_cds}]

    @property
    def key(self) -> str:
        return make_key(self.type, self.contig, self.pos, self.ref, self.alt)

    @property
    def length(self) -> int:
        return len(self.alt) - len(self.ref)

    @property
    def len_mod3(self) -> int:
        return abs(self.length) % 3

    @property
    def is_indel(self) -> bool:
        return self.type != "snp"

    def callable_span(self, contig_length: int) -> tuple[int, int]:
        """Reference bases whose state decides whether this variant can be assessed."""
        if self.type == "ins":
            return max(self.start, 1), min(self.end, contig_length)
        return self.start, self.end


@dataclass
class Call:
    sample: str
    key: str
    seq_start: int
    ins_seq_raw: str = ""
    flags: list[str] = field(default_factory=list)
    frameshift: str = "."


@dataclass
class Ambiguous:
    sample: str
    contig: str
    pos: int
    ref: str
    code: str


@dataclass
class SampleResult:
    id: str
    records: dict[str, str] = field(default_factory=dict)  # contig -> record id
    covered: list[tuple[str, int, int]] = field(default_factory=list)
    missing: list[tuple[str, int, int]] = field(default_factory=list)
    missing_insertions: list[tuple[str, int, int]] = field(default_factory=list)  # contig, p, len
    overhang_bases: int = 0
    summary: dict = field(default_factory=dict)


@dataclass
class CallSet:
    reference: Reference
    set_name: str
    params: dict
    samples: list[SampleResult] = field(default_factory=list)
    variants: dict[str, Variant] = field(default_factory=dict)
    calls: list[Call] = field(default_factory=list)
    ambiguous: list[Ambiguous] = field(default_factory=list)
    alignments: list[tuple[str, str]] = field(default_factory=list)
    annotation: dict | None = None
    cds_status: list[dict] = field(default_factory=list)
    tool_version: str = __version__

    @property
    def sample_ids(self) -> list[str]:
        return [s.id for s in self.samples]

    def sort_key(self, key: str):
        v = self.variants[key]
        return (self.reference.index(v.contig), v.pos, v.type, v.ref, v.alt)

    def sorted_keys(self) -> list[str]:
        return sorted(self.variants, key=self.sort_key)

    def calls_by_sample(self) -> dict[str, list[Call]]:
        out: dict[str, list[Call]] = {s: [] for s in self.sample_ids}
        for c in self.calls:
            out[c.sample].append(c)
        return out


# --------------------------------------------------------------------------- reference context


class RefContext:
    """Per-contig reference helpers: 1-based base access and homopolymer runs."""

    def __init__(self, contig: str, seq: str):
        self.contig = contig
        self.seq = seq
        self.length = len(seq)
        L = self.length
        self.run_start = np.zeros(L + 2, dtype=np.int64)
        self.run_end = np.zeros(L + 2, dtype=np.int64)
        q = 1
        while q <= L:
            r = q
            while r < L and seq[r] == seq[q - 1]:
                r += 1
            self.run_start[q : r + 1] = q
            self.run_end[q : r + 1] = r
            q = r + 1

    def base(self, q: int) -> str:
        return self.seq[q - 1]

    def span(self, s: int, e: int) -> str:
        return self.seq[s - 1 : e]

    def run_at(self, q: int) -> tuple[int, int] | None:
        if not 1 <= q <= self.length or self.base(q) not in ACGT:
            return None
        return int(self.run_start[q]), int(self.run_end[q])

    def homopolymer(
        self, positions: range, bases: set[str] | None = None
    ) -> tuple[int, tuple[str, int, int] | None]:
        """Longest ACGT run containing any of `positions` (optionally restricted by base)."""
        best = None
        for q in positions:
            run = self.run_at(q)
            if run is None or (bases is not None and self.base(q) not in bases):
                continue
            if best is None or run[1] - run[0] > best[1] - best[0]:
                best = run
        if best is None:
            return 0, None
        return best[1] - best[0] + 1, (self.contig, best[0], best[1])


def build_variant(ctx: RefContext, vtype: str, start: int, end: int, ins_seq: str = "") -> Variant:
    """Anchor an (already normalised) event VCF-style and annotate hp context and slide span."""
    L = ctx.length
    if vtype == "snp":
        pos, ref, alt = start, ctx.base(start), ins_seq
        v = Variant(ctx.contig, "snp", start, end, pos, ref, alt)
        v.hp_len, v.hp_run = ctx.homopolymer(range(start - 1, start + 2))
        v.slide = (start, start)
        return v
    if vtype == "del":
        if start > 1:
            pos, ref, alt = start - 1, ctx.span(start - 1, end), ctx.base(start - 1)
        else:
            if end >= L:
                raise CallError(f"{ctx.contig}: deletion of the entire contig cannot be anchored")
            pos, ref, alt = 1, ctx.span(1, end + 1), ctx.base(end + 1)
        v = Variant(ctx.contig, "del", start, end, pos, ref, alt)
        v.hp_len, v.hp_run = ctx.homopolymer(range(start - 1, end + 2))
        v.slide = deletion_slide(ctx, start, end)
        return v
    # insertion between p=start and p+1=end
    p = start
    if p >= 1:
        pos, ref, alt = p, ctx.base(p), ctx.base(p) + ins_seq
    else:
        pos, ref, alt = 1, ctx.base(1), ins_seq + ctx.base(1)
    v = Variant(ctx.contig, "ins", p, p + 1, pos, ref, alt, ins_seq=ins_seq)
    v.hp_len, v.hp_run = ctx.homopolymer(range(p, p + 2), set(ins_seq))
    v.slide = insertion_slide(ctx, p, ins_seq)
    return v


def deletion_slide(ctx: RefContext, s: int, e: int) -> tuple[int, int]:
    ls, le = s, e
    while ls > 1 and ctx.base(ls - 1) == ctx.base(le):
        ls, le = ls - 1, le - 1
    rs, re_ = s, e
    while re_ < ctx.length and ctx.base(rs) == ctx.base(re_ + 1):
        rs, re_ = rs + 1, re_ + 1
    return ls, re_


def insertion_slide(ctx: RefContext, p: int, seq: str) -> tuple[int, int]:
    lp, ls = p, seq
    while lp >= 1 and ctx.base(lp) == ls[-1]:
        ls = ls[-1] + ls[:-1]
        lp -= 1
    rp, rs = p, seq
    while rp + 1 <= ctx.length and ctx.base(rp + 1) == rs[0]:
        rs = rs[1:] + rs[0]
        rp += 1
    return max(lp, 1), min(rp + 1, ctx.length)


# --------------------------------------------------------------------------- calling


@dataclass
class _Event:
    type: str
    start: int
    end: int
    seq_start: int
    raw: str = ""  # SNP alt base or raw inserted sequence


def classify_columns(ca: ContigAlignment, i: int) -> np.ndarray:
    R = ca.ref_row
    S = ca.rows[i]
    ncol = R.size
    rgap = R == GAP
    sgap = S == GAP
    cls = np.full(ncol, MATCH, dtype=np.int8)
    cov = ca.covered_columns(i)
    inside = np.zeros(ncol, dtype=bool)
    if cov is not None:
        inside[cov[0] : cov[1] + 1] = True
    rb = ~rgap
    cls[rgap & sgap] = SKIP
    cls[rgap & ~sgap & inside] = INS
    cls[rgap & ~sgap & ~inside] = OVERHANG
    cls[rb & sgap & ~inside] = MISSING
    cls[rb & sgap & inside] = DEL
    s_n = S == N
    cls[rb & s_n] = MISSING
    diff = rb & ~sgap & ~s_n & (S != R)
    racgt = np.isin(R, _ACGT_CODES)
    sacgt = np.isin(S, _ACGT_CODES)
    cls[diff & racgt & sacgt] = SNP
    cls[diff & racgt & ~sacgt] = AMBIG
    cls[diff & ~racgt] = UNCALLABLE
    return cls


def _intervals(positions: np.ndarray) -> list[tuple[int, int]]:
    if positions.size == 0:
        return []
    breaks = np.flatnonzero(np.diff(positions) != 1)
    starts = np.concatenate(([positions[0]], positions[breaks + 1]))
    ends = np.concatenate((positions[breaks], [positions[-1]]))
    return [(int(s), int(e)) for s, e in zip(starts, ends, strict=True)]


def _normalise(events: list[_Event], status: np.ndarray, ins_points: set[int], ctx: RefContext):
    """Left-normalise in place, never sliding over the sample's own non-reference state."""
    for ev in sorted(events, key=lambda x: (x.start, x.type)):
        if ev.type == "del":
            s, e = ev.start, ev.end
            while (
                s > 1
                and ctx.base(s - 1) == ctx.base(e)
                and status[s - 1] == PLAIN
                and (s - 1) not in ins_points
            ):
                status[s - 1] = S_DELETED
                status[e] = PLAIN
                s, e = s - 1, e - 1
                ev.seq_start -= 1
            ev.start, ev.end = s, e
        elif ev.type == "ins":
            p, seq = ev.start, ev.raw
            while (
                p >= 1
                and ctx.base(p) == seq[-1]
                and status[p] == PLAIN
                and (p - 1) not in ins_points
            ):
                ins_points.discard(p)
                ins_points.add(p - 1)
                seq = seq[-1] + seq[:-1]
                p -= 1
                ev.seq_start -= 1
            ev.start, ev.end, ev.raw = p, p + 1, seq


def _alt_seq(raw: str) -> str:
    return "".join(c if c in ACGT else "N" for c in raw)


def _column_events(
    S: np.ndarray, cls: np.ndarray, ref_pos: np.ndarray
) -> tuple[list[_Event], list[tuple[int, int]]]:
    """Group classified columns into SNP / deletion / insertion events.

    Each column is classified independently and runs are grouped afterwards, so an event can
    start in the column right after another closes. Skip columns are transparent.
    Returns the events and the all-N insertions as (p, length).
    """
    rank = np.cumsum(cls != SKIP) - 1
    chars_before = np.cumsum(S != GAP) - (S != GAP)
    events: list[_Event] = []
    missing_ins: list[tuple[int, int]] = []

    def flush(kind, cols):
        cols = np.array(cols)
        seq_start = int(chars_before[cols[0]]) + 1
        if kind == DEL:
            positions = ref_pos[cols]
            events.append(_Event("del", int(positions.min()), int(positions.max()), seq_start))
            return
        raw = bytes(S[cols]).decode()
        p = int(ref_pos[cols[0]])
        if set(raw) == {"N"}:
            missing_ins.append((p, len(raw)))
        else:
            events.append(_Event("ins", p, p + 1, seq_start, raw))

    run_kind, run_cols, run_last = None, [], -2
    for c in np.flatnonzero((cls == SNP) | (cls == DEL) | (cls == INS)):
        k = int(cls[c])
        r = int(rank[c])
        if run_kind == k and r == run_last + 1:
            run_cols.append(c)
            run_last = r
            continue
        if run_kind is not None:
            flush(run_kind, run_cols)
            run_kind = None
        if k == SNP:
            q = int(ref_pos[c])
            events.append(_Event("snp", q, q, int(chars_before[c]) + 1, chr(S[c])))
        else:
            run_kind, run_cols, run_last = k, [c], r
    if run_kind is not None:
        flush(run_kind, run_cols)
    return events, missing_ins


def call_contig(ca: ContigAlignment, hp_min: int):
    """Call one contig. Returns per-sample results, variants, calls and ambiguous calls."""
    ctx = RefContext(ca.contig, ca.ref_seq)
    L = ctx.length
    rb = ca.ref_row != GAP
    ref_cols = np.flatnonzero(rb)
    ref_pos = ca.ref_pos
    variants: dict[str, Variant] = {}
    calls: list[Call] = []
    ambiguous: list[Ambiguous] = []
    per_sample: list[SampleResult] = []

    for i, sid in enumerate(ca.sample_ids):
        S = ca.rows[i]
        cls = classify_columns(ca, i)
        status = np.zeros(L + 2, dtype=np.int8)
        status[ref_pos[rb]] = _CLASS_TO_STATUS[cls[rb]]
        status[0] = status[L + 1] = S_MISSING

        events, missing_ins = _column_events(S, cls, ref_pos)

        ins_points = {ev.start for ev in events if ev.type == "ins"}
        _normalise(events, status, ins_points, ctx)

        missing_positions = np.flatnonzero(status[1 : L + 1] == S_MISSING) + 1
        missing_ivs = _intervals(missing_positions)
        missing_mask = status == S_MISSING

        for ev in events:
            if ev.type == "snp":
                v = build_variant(ctx, "snp", ev.start, ev.end, ev.raw)
            elif ev.type == "del":
                v = build_variant(ctx, "del", ev.start, ev.end)
            else:
                v = build_variant(ctx, "ins", ev.start, ev.end, _alt_seq(ev.raw))
            key = v.key
            variants.setdefault(key, v)
            flags = []
            if v.hp_len >= hp_min:
                flags.append("homopolymer")
            if ev.type == "ins" and set(ev.raw) - set(ACGT):
                flags.append("contains_ambiguity")
            flank = (v.start, v.end) if v.type == "ins" else (v.start - 1, v.end + 1)
            if any(1 <= q <= L and missing_mask[q] for q in flank):
                flags.append("near_missing")
            calls.append(
                Call(
                    sample=sid,
                    key=key,
                    seq_start=ev.seq_start,
                    ins_seq_raw=ev.raw if ev.type == "ins" else "",
                    flags=flags,
                )
            )

        for q in np.flatnonzero(status[1 : L + 1] == S_AMBIG) + 1:
            col = int(ref_cols[q - 1])
            ambiguous.append(Ambiguous(sid, ca.contig, int(q), ctx.base(int(q)), chr(S[col])))

        cov = ca.covered_columns(i)
        sr = SampleResult(id=sid, records={ca.contig: sid})
        if cov is not None:
            sr.covered.append((ca.contig, int(ref_pos[cov[0]]), int(ref_pos[cov[1]])))
        sr.missing = [(ca.contig, s, e) for s, e in missing_ivs]
        sr.missing_insertions = [(ca.contig, p, n) for p, n in missing_ins]
        sr.overhang_bases = int(np.count_nonzero(cls == OVERHANG))
        per_sample.append(sr)

    return per_sample, variants, calls, ambiguous


def call_set(alignment: SetAlignment, reference: Reference, set_name: str, params: dict) -> CallSet:
    cs = CallSet(reference=reference, set_name=set_name, params=dict(params))
    by_id: dict[str, SampleResult] = {}
    for contig in reference.ids:
        samples, variants, calls, ambiguous = call_contig(alignment[contig], params["hp_min"])
        for sr in samples:
            agg = by_id.setdefault(sr.id, SampleResult(id=sr.id))
            agg.records.update(sr.records)
            agg.covered += sr.covered
            agg.missing += sr.missing
            agg.missing_insertions += sr.missing_insertions
            agg.overhang_bases += sr.overhang_bases
        cs.variants.update(variants)
        cs.calls += calls
        cs.ambiguous += ambiguous
    cs.samples = list(by_id.values())
    order = {s: i for i, s in enumerate(cs.sample_ids)}
    cs.calls.sort(key=lambda c: (order[c.sample], cs.sort_key(c.key)))
    cs.ambiguous.sort(key=lambda a: (order[a.sample], reference.index(a.contig), a.pos))
    summarise(cs)
    return cs


def summarise(cs: CallSet) -> None:
    """(Re)compute per-sample summary counts."""
    total_len = sum(cs.reference.length(c) for c in cs.reference.ids)
    per = {s.id: {"n_snp": 0, "n_ins": 0, "n_del": 0, "n_ambiguous": 0} for s in cs.samples}
    for c in cs.calls:
        per[c.sample]["n_" + cs.variants[c.key].type] += 1
    for a in cs.ambiguous:
        per[a.sample]["n_ambiguous"] += 1
    cds_counts = {s.id: dict.fromkeys(_CDS_COUNTS, 0) for s in cs.samples}
    for st in cs.cds_status:
        for col, (field_name, value) in _CDS_COUNTS.items():
            if (st["stop"]["stop_status"] if field_name == "stop" else st[field_name]) == value:
                cds_counts[st["sample"]][col] += 1
    for s in cs.samples:
        missing = sum(e - st + 1 for _, st, e in s.missing)
        s.summary = {
            **per[s.id],
            "missing_bases": missing,
            "covered_fraction": round((total_len - missing) / total_len, 4) if total_len else 0,
            "overhang_bases": s.overhang_bases,
            **(cds_counts[s.id] if cs.annotation else dict.fromkeys(_CDS_COUNTS, ".")),
        }


_CDS_COUNTS = {
    "frameshifted_cds": ("status", "frameshift"),
    "restored_cds": ("status", "restored"),
    "premature_stop_cds": ("stop", "premature_stop"),
    "stop_lost_cds": ("stop", "stop_lost"),
}
