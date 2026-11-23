"""Disagreement analysis: loci, local haplotypes, pair reasons and discordance outputs."""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path

from .call import CallSet
from .haplotype import SampleEvents, window_haplotype
from .io import write_tsv
from .params import HAPLOTYPE_DISPLAY_FLANK, HAPLOTYPE_DISPLAY_MAX

# Pair reasons in precedence order (first match wins)
REASONS = [
    "missing_data",
    "ambiguity",
    "equivalent_representation",
    "paired_indel",
    "indel_length",
    "indel_sequence",
    "indel_presence",
    "different_allele",
    "snp_presence",
    "mixed",
]

# Locus ranking: substantive callable differences before artefacts and missing data
LOCUS_RANKING = [
    "paired_indel",
    "indel_length",
    "indel_sequence",
    "indel_presence",
    "different_allele",
    "snp_presence",
    "mixed",
    "equivalent_representation",
    "ambiguity",
    "missing_data",
]

FLAGS = ["near_end", "homopolymer", "frameshift"]


@dataclass
class Locus:
    id: str
    contig: str
    start: int
    end: int
    keys: list[str] = field(default_factory=list)


@dataclass
class PairResult:
    a: str
    b: str
    hap_a: str
    hap_b: str
    reason: str
    only_a: list[str]
    only_b: list[str]
    flags: list[str]


@dataclass
class LocusResult:
    locus: Locus
    labels: dict[str, str]  # sample -> H0/H1/.../missing/ambiguous
    haplotypes: dict[str, str]  # label -> sequence
    pairs: list[PairResult]
    haplotype_concordant: bool
    callable_discordant: bool
    genes: list[str]

    @property
    def discordant(self) -> bool:
        return bool(self.pairs)

    @property
    def reasons(self) -> list[str]:
        present = {p.reason for p in self.pairs}
        return [r for r in REASONS if r in present]

    @property
    def primary_reason(self) -> str | None:
        present = {p.reason for p in self.pairs}
        return next((r for r in LOCUS_RANKING if r in present), None)

    @property
    def flags(self) -> list[str]:
        present = {f for p in self.pairs for f in p.flags}
        return [f for f in FLAGS if f in present]

    @property
    def odd_samples(self) -> list[str]:
        """Callable samples differing from the unique majority haplotype ([] on a tie)."""
        counts = Counter(lab for lab in self.labels.values() if lab.startswith("H"))
        top = counts.most_common(2)
        if len(top) < 2 or top[0][1] == top[1][1]:
            return []
        majority = top[0][0]
        return [s for s, lab in self.labels.items() if lab.startswith("H") and lab != majority]

    def pair(self, a: str, b: str) -> PairResult | None:
        for p in self.pairs:
            if {p.a, p.b} == {a, b}:
                return p
        return None


def build_loci(cs: CallSet, locus_join_gap: int) -> list[Locus]:
    spans = sorted(
        (cs.reference.index(v.contig), v.slide[0], v.slide[1], key)
        for key, v in cs.variants.items()
    )
    loci: list[Locus] = []
    cur = None
    for ci, s, e, key in spans:
        contig = cs.reference.ids[ci]
        if cur is not None and cur.contig == contig and s - cur.end - 1 <= locus_join_gap:
            cur.end = max(cur.end, e)
            cur.keys.append(key)
        else:
            cur = Locus("", contig, s, e, [key])
            loci.append(cur)
    for loc in loci:
        loc.id = f"L{cs.reference.index(loc.contig) + 1}:{loc.start}-{loc.end}"
        loc.keys.sort(key=cs.sort_key)
    return loci


def pair_reason(status_a, seq_a, status_b, seq_b, da, db, variants) -> str:
    if "missing" in (status_a, status_b):
        return "missing_data"
    if "ambiguous" in (status_a, status_b):
        return "ambiguity"
    if seq_a == seq_b:
        return "equivalent_representation"
    d = da | db

    def ins_and_del(keys):
        types = {variants[k].type for k in keys}
        return "ins" in types and "del" in types

    if len(seq_a) == len(seq_b) and (ins_and_del(da) or ins_and_del(db)):
        return "paired_indel"
    if all(variants[k].is_indel for k in d):
        if da and db:
            return "indel_length" if len(seq_a) != len(seq_b) else "indel_sequence"
        return "indel_presence"
    if all(variants[k].type == "snp" for k in d):
        pos_a = {variants[k].start for k in da}
        pos_b = {variants[k].start for k in db}
        return "different_allele" if pos_a & pos_b else "snp_presence"
    return "mixed"


def _near_end(se: SampleEvents, start: int, end: int, end_dist: int) -> bool:
    if se.covered is None:
        return False
    cs_, ce = se.covered
    return start - cs_ <= end_dist or ce - end <= end_dist


def analyse(
    cs: CallSet, sample_events: dict[tuple[str, str], SampleEvents], thresholds: dict
) -> list[LocusResult]:
    """Analyse every locus (discordant or not), in reference order."""
    samples = cs.sample_ids
    sample_keys = {s: set() for s in samples}
    frameshifted = set()
    for c in cs.calls:
        sample_keys[c.sample].add(c.key)
        if "fs" in c.frameshift:
            frameshifted.add((c.sample, c.key))
    cds = (cs.annotation or {}).get("cds", [])

    results = []
    for loc in build_loci(cs, thresholds["locus_join_gap"]):
        ref_seq = cs.reference.seq(loc.contig)
        loc_keys = set(loc.keys)
        keysets = {s: frozenset(sample_keys[s] & loc_keys) for s in samples}
        haps = {
            s: window_haplotype(ref_seq, sample_events[(s, loc.contig)], loc.start, loc.end)
            for s in samples
        }
        ref_window = ref_seq[loc.start - 1 : loc.end]
        counts = Counter(seq for st, seq in haps.values() if st == "callable" and seq != ref_window)
        ordered = sorted(counts, key=lambda x: (-counts[x], x))
        label_of_seq = {ref_window: "H0", **{seq: f"H{i + 1}" for i, seq in enumerate(ordered)}}
        labels = {
            s: (label_of_seq[seq] if st == "callable" else st) for s, (st, seq) in haps.items()
        }
        haplotypes = {
            lab: seq
            for seq, lab in sorted(label_of_seq.items(), key=lambda x: int(x[1][1:]))
            if lab in labels.values()
        }

        pairs = []
        for i, a in enumerate(samples):
            for b in samples[i + 1 :]:
                if keysets[a] == keysets[b]:
                    continue
                da, db = set(keysets[a] - keysets[b]), set(keysets[b] - keysets[a])
                reason = pair_reason(*haps[a], *haps[b], da, db, cs.variants)
                flags = []
                sea, seb = sample_events[(a, loc.contig)], sample_events[(b, loc.contig)]
                ed = thresholds["end_dist"]
                if _near_end(sea, loc.start, loc.end, ed) or _near_end(seb, loc.start, loc.end, ed):
                    flags.append("near_end")
                if any(cs.variants[k].hp_len >= thresholds["hp_min"] for k in da | db):
                    flags.append("homopolymer")
                if any((a, k) in frameshifted for k in da) or any(
                    (b, k) in frameshifted for k in db
                ):
                    flags.append("frameshift")
                pairs.append(
                    PairResult(
                        a,
                        b,
                        labels[a],
                        labels[b],
                        reason,
                        sorted(da, key=cs.sort_key),
                        sorted(db, key=cs.sort_key),
                        flags,
                    )
                )

        callable_labels = {lab for lab in labels.values() if lab.startswith("H")}
        genes = [
            c["gene"]
            for c in cds
            if c["contig"] == loc.contig
            and any(s <= loc.end and e >= loc.start for s, e in c["segments"])
        ]
        results.append(
            LocusResult(
                locus=loc,
                labels=labels,
                haplotypes=haplotypes,
                pairs=pairs,
                haplotype_concordant=len(callable_labels) <= 1,
                callable_discordant=len(callable_labels) >= 2,
                genes=genes,
            )
        )
    return results


CDS_EFFECTS = ["fs", "fs_compensated", "premature_stop"]


def cds_effects(cs: CallSet, results: list[LocusResult]) -> dict[str, list[str]]:
    """locus id -> ['gene:effect(sample,...)', ...] for frameshifts and premature stops.

    A frameshift is attributed to the locus holding the frameshifting call. A premature stop
    is attributed to a locus when one of the sample's calls there overlaps any base of the
    stop codon; a stop caused downstream by a frameshift is therefore not attributed (the
    locus already shows `fs`).
    """
    order = {s: i for i, s in enumerate(cs.sample_ids)}
    calls_at: dict[str, list] = {}
    for c in cs.calls:
        calls_at.setdefault(c.key, []).append(c)
    stops = [
        st
        for st in cs.cds_status
        if st["stop"]["stop_status"] == "premature_stop" and st["stop"].get("stop_codon")
    ]
    out = {}
    for r in results:
        found: dict[tuple[str, str], set[str]] = {}
        for key in r.locus.keys:
            for c in calls_at.get(key, []):
                for entry in c.frameshift.split("|"):
                    gene, _, effect = entry.rpartition(":")
                    if effect in ("fs", "fs_compensated"):
                        found.setdefault((gene, effect), set()).add(c.sample)
        for st in stops:
            lo, hi = st["stop"]["stop_codon"]
            if st["contig"] != r.locus.contig or hi < r.locus.start or lo > r.locus.end:
                continue
            for key in r.locus.keys:
                v = cs.variants[key]
                if (
                    v.start <= hi
                    and v.end >= lo
                    and any(c.sample == st["sample"] for c in calls_at.get(key, []))
                ):
                    found.setdefault((st["gene"], "premature_stop"), set()).add(st["sample"])
                    break
        out[r.locus.id] = [
            f"{gene}:{effect}({','.join(sorted(samples, key=order.__getitem__))})"
            for (gene, effect), samples in sorted(
                found.items(), key=lambda x: (x[0][0], CDS_EFFECTS.index(x[0][1]))
            )
        ]
    return out


def _display(seq: str) -> str:
    if len(seq) <= HAPLOTYPE_DISPLAY_MAX:
        return seq or "-"
    f = HAPLOTYPE_DISPLAY_FLANK
    return f"{seq[:f]}…{seq[-f:]}(len={len(seq)})"


DISCORDANCE_HEADER = [
    "set",
    "locus",
    "contig",
    "start",
    "end",
    "keys",
    "haplotype_concordant",
    "callable_discordant",
    "primary_reason",
    "reasons",
    "flags",
    "genes",
    "haplotypes",
    "odd_samples",
    "cds_effect",
]

PAIRS_HEADER = [
    "set",
    "locus",
    "sample_a",
    "sample_b",
    "hap_a",
    "hap_b",
    "reason",
    "only_a",
    "only_b",
    "flags",
]

SUMMARY_HEADER = [
    "set",
    "reason",
    "n_loci_primary",
    "n_loci_callable_primary",
    "n_loci_any",
    "n_pairs",
]


def write_outputs(cs: CallSet, results: list[LocusResult], paths: dict[str, Path]) -> None:
    samples = cs.sample_ids
    disc = [r for r in results if r.discordant]
    effects = cds_effects(cs, disc)
    rows = []
    for r in disc:
        loc = r.locus
        haps = ";".join(f"{lab}={_display(seq)}" for lab, seq in r.haplotypes.items())
        rows.append(
            [
                cs.set_name,
                loc.id,
                loc.contig,
                loc.start,
                loc.end,
                loc.keys,
                r.haplotype_concordant,
                r.callable_discordant,
                r.primary_reason,
                r.reasons,
                r.flags,
                r.genes,
                haps,
                r.odd_samples,
                effects[loc.id],
            ]
            + [r.labels[s] for s in samples]
        )
    write_tsv(paths["discordance"], DISCORDANCE_HEADER + samples, rows)

    pair_rows = []
    for r in disc:
        loc = r.locus
        for p in r.pairs:
            pair_rows.append(
                [
                    cs.set_name,
                    loc.id,
                    p.a,
                    p.b,
                    p.hap_a,
                    p.hap_b,
                    p.reason,
                    p.only_a,
                    p.only_b,
                    p.flags,
                ]
            )
    write_tsv(paths["discordance_pairs"], PAIRS_HEADER, pair_rows)

    summary = []
    for reason in REASONS:
        primary = [r for r in disc if r.primary_reason == reason]
        summary.append(
            [
                cs.set_name,
                reason,
                len(primary),
                sum(1 for r in primary if r.callable_discordant),
                sum(1 for r in disc if reason in r.reasons),
                sum(1 for r in disc for p in r.pairs if p.reason == reason),
            ]
        )
    write_tsv(paths["discordance_summary"], SUMMARY_HEADER, summary)


NOT_REAL = {"missing_data", "ambiguity", "equivalent_representation"}

SET_SUMMARY_HEADER = [
    "set",
    "n_samples",
    "n_variant_keys",
    "discordant_loci",
    "discordant_positions",
    "callable_discordant_loci",
    "uncallable_only_loci",
    *REASONS,
    *FLAGS,
    "fs_loci",
    "premature_stop_loci",
    "odd_one_out",
    "most_discordant_pair",
]


def pair_differences(disc: list[LocusResult], a: str, b: str) -> tuple[int, int]:
    """(real differences, equivalent-representation loci) between two samples."""
    real = equiv = 0
    for r in disc:
        p = r.pair(a, b)
        if p is None:
            continue
        if p.reason == "equivalent_representation":
            equiv += 1
        elif p.reason not in NOT_REAL:
            real += 1
    return real, equiv


def set_summary_row(cs: CallSet, results: list[LocusResult]) -> list:
    """One row describing the disagreement between the samples of one set."""
    samples = cs.sample_ids
    disc = [r for r in results if r.discordant]
    n_callable = sum(r.callable_discordant for r in disc)
    primary = Counter(r.primary_reason for r in disc)
    flags = Counter(f for r in disc for f in r.flags)
    effects = cds_effects(cs, disc).values()
    odd = Counter(s for r in disc for s in r.odd_samples)
    pairs = [
        (pair_differences(disc, a, b), a, b)
        for i, a in enumerate(samples)
        for b in samples[i + 1 :]
    ]
    top = max(pairs, key=lambda x: x[0][0], default=None)
    return [
        cs.set_name,
        len(samples),
        len(cs.variants),
        len(disc),
        sum(r.locus.end - r.locus.start + 1 for r in disc),
        n_callable,
        len(disc) - n_callable,
        *(primary[reason] for reason in REASONS),
        *(flags[f] for f in FLAGS),
        *(
            (
                sum(any(":fs" in e for e in eff) for eff in effects),
                sum(any(":premature_stop" in e for e in eff) for eff in effects),
            )
            if cs.annotation
            else (None, None)  # not evaluated without --annotation
        ),
        [f"{s}:{odd[s]}" for s in samples if odd[s]],
        f"{top[1]},{top[2]}:{top[0][0]}" if top and top[0][0] else None,
    ]
