"""Genotype states, wide variants table, pairwise concordance and truth metrics."""

from __future__ import annotations

from pathlib import Path

from .call import TYPE_NAMES, CallSet, Variant
from .discord import LocusResult
from .haplotype import SampleEvents
from .io import write_tsv

STATES = ["alt", "missing", "ambiguous", "deleted", "other_allele", "ref"]
UNASSESSABLE = {"missing", "ambiguous"}
TYPES = ["snp", "insertion", "deletion"]


class CompareError(ValueError):
    pass


def state(v: Variant, se: SampleEvents) -> str:
    """State of variant v in one sample (first matching rule, §7.1)."""
    if v.key in se.keys:
        return "alt"
    lo, hi = v.callable_span(se.length)
    if se.missing_in(lo, hi) or (v.type == "ins" and v.start in se.missing_insertions):
        return "missing"
    if se.ambiguous_in(lo, hi):
        return "ambiguous"
    if se.deleted_in(lo, hi):
        return "deleted"
    if se.snps_in(lo, hi):
        return "other_allele"
    if v.type == "ins":
        if v.start in se.insertions:
            return "other_allele"
    elif any(lo <= p and p + 1 <= hi for p in se.insertions):
        return "other_allele"
    return "ref"


def state_table(
    cs: CallSet, sample_events: dict[tuple[str, str], SampleEvents]
) -> dict[str, dict[str, str]]:
    """key -> sample -> state, for every variant key in the set."""
    return {
        key: {s: state(v, sample_events[(s, v.contig)]) for s in cs.sample_ids}
        for key, v in cs.variants.items()
    }


VARIANTS_HEADER = [
    "key",
    "contig",
    "type",
    "pos",
    "ref",
    "alt",
    "length",
    "hp_len",
    "cds",
    "locus",
]


def write_variants_wide(cs: CallSet, states, loci: list[LocusResult], path: Path) -> None:
    """One row per variant key with its attributes, locus and state in every sample."""
    locus_of = {k: r for r in loci for k in r.locus.keys}
    rows = []
    for key in cs.sorted_keys():
        v = cs.variants[key]
        rows.append(
            [
                key,
                v.contig,
                v.type,
                v.pos,
                v.ref,
                v.alt,
                v.length,
                v.hp_len,
                "|".join(f"{x['gene']}:{x['delta_cds']}" for x in v.cds) or ".",
                locus_of[key].locus.id,
            ]
            + [states[key][s] for s in cs.sample_ids]
        )
    write_tsv(path, VARIANTS_HEADER + cs.sample_ids, rows)


PAIRWISE_HEADER = [
    "sample_a",
    "sample_b",
    "type",
    "shared",
    "only_a",
    "only_b",
    "comparable",
    "excluded",
    "jaccard",
]


def _ratio(num: int, den: int) -> str:
    return f"{num / den:.4f}" if den else "NA"


def pairwise_rows(cs: CallSet, states) -> list[list]:
    samples = cs.sample_ids
    keys_by_type = {t: [] for t in TYPES}
    for key, v in cs.variants.items():
        keys_by_type[TYPE_NAMES[v.type]].append(key)
    rows = []
    for i, a in enumerate(samples):
        for b in samples[i + 1 :]:
            totals = [0] * 5
            for t in TYPES + ["all"]:
                if t == "all":
                    shared, only_a, only_b, comparable, excluded = totals
                else:
                    shared = only_a = only_b = comparable = excluded = 0
                    for key in keys_by_type[t]:
                        sa, sb = states[key][a], states[key][b]
                        if sa in UNASSESSABLE or sb in UNASSESSABLE:
                            excluded += 1
                            continue
                        comparable += 1
                        if sa == "alt" and sb == "alt":
                            shared += 1
                        elif sa == "alt":
                            only_a += 1
                        elif sb == "alt":
                            only_b += 1
                    for j, x in enumerate((shared, only_a, only_b, comparable, excluded)):
                        totals[j] += x
                rows.append(
                    [
                        a,
                        b,
                        t,
                        shared,
                        only_a,
                        only_b,
                        comparable,
                        excluded,
                        _ratio(shared, shared + only_a + only_b),
                    ]
                )
    return rows


def write_pairwise(cs: CallSet, states, path: Path) -> None:
    write_tsv(path, PAIRWISE_HEADER, pairwise_rows(cs, states))


# --------------------------------------------------------------------------- truth

TRUTH_METRICS_HEADER = [
    "set",
    "truth",
    "sample",
    "type",
    "tp",
    "fp",
    "fn",
    "uncalled",
    "truth_uncallable",
    "equivalent",
    "precision",
    "recall",
    "f1",
]

TRUTH_DETAIL_HEADER = [
    "set",
    "truth",
    "sample",
    "key",
    "type",
    "class",
    "query_state",
    "truth_state",
    "locus",
    "reason",
]


def truth_tables(cs: CallSet, states, loci: list[LocusResult], truth: str):
    if truth not in cs.sample_ids:
        raise CompareError(
            f"--truth {truth!r} is not a sample in this set (samples: {', '.join(cs.sample_ids)})"
        )
    locus_of = {k: r for r in loci for k in r.locus.keys}
    metrics, detail = [], []
    for q in cs.sample_ids:
        if q == truth:
            continue
        counts = {
            t: dict.fromkeys(("TP", "FP", "FN", "UNCALLED", "TRUTH_UNCALLABLE", "eq"), 0)
            for t in TYPES
        }
        for key in cs.sorted_keys():
            v = cs.variants[key]
            t = TYPE_NAMES[v.type]
            ts, qs = states[key][truth], states[key][q]
            if ts == "alt":
                cls = "TP" if qs == "alt" else "UNCALLED" if qs in UNASSESSABLE else "FN"
            elif qs == "alt":
                cls = "TRUTH_UNCALLABLE" if ts in UNASSESSABLE else "FP"
            else:
                continue
            r = locus_of[key]
            reason = None
            if cls != "TP":
                pair = r.pair(truth, q)
                reason = pair.reason if pair else None
            counts[t][cls] += 1
            if cls in ("FP", "FN") and reason == "equivalent_representation":
                counts[t]["eq"] += 1
            detail.append([cs.set_name, truth, q, key, t, cls, qs, ts, r.locus.id, reason])
        total = dict.fromkeys(counts["snp"], 0)
        for t in TYPES + ["all"]:
            c = total if t == "all" else counts[t]
            if t != "all":
                for k in total:
                    total[k] += c[k]
            tp, fp, fn = c["TP"], c["FP"], c["FN"]
            metrics.append(
                [
                    cs.set_name,
                    truth,
                    q,
                    t,
                    tp,
                    fp,
                    fn,
                    c["UNCALLED"],
                    c["TRUTH_UNCALLABLE"],
                    c["eq"],
                    _ratio(tp, tp + fp),
                    _ratio(tp, tp + fn),
                    _ratio(2 * tp, 2 * tp + fp + fn),
                ]
            )
    return metrics, detail


def write_truth(cs: CallSet, states, loci, truth: str, paths: dict[str, Path]) -> None:
    metrics, detail = truth_tables(cs, states, loci, truth)
    write_tsv(paths["truth_metrics"], TRUTH_METRICS_HEADER, metrics)
    write_tsv(paths["truth_detail"], TRUTH_DETAIL_HEADER, detail)
