"""Human-readable discordance report: the first place to look."""

from __future__ import annotations

from collections import Counter
from pathlib import Path

from .call import CallSet
from .discord import LOCUS_RANKING, LocusResult, _display, cds_effects

REPORT_MAX_LOCI = 50
REPORT_MAX_MATRIX_SAMPLES = 12
REPORT_TOP_PAIRS = 20
NOT_REAL = {"missing_data", "ambiguity", "equivalent_representation"}


def _pair_counts(disc: list[LocusResult], a: str, b: str) -> tuple[int, int]:
    """(real differences, equivalent-representation loci) for one sample pair."""
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


def _pair_cell(real: int, equiv: int) -> str:
    return f"{real} (+{equiv} equiv)" if equiv else str(real)


def _who_carries_what(r: LocusResult, samples: list[str]) -> str:
    counts = Counter(r.labels.values())
    haps = sorted(
        (lab for lab in counts if lab.startswith("H")), key=lambda x: (-counts[x], int(x[1:]))
    )
    big = None
    if haps and counts[haps[0]] > 1 and (len(haps) == 1 or counts[haps[0]] > counts[haps[1]]):
        big = haps[0]
    parts = []
    for lab in haps + [x for x in ("missing", "ambiguous") if x in counts]:
        if lab.startswith("H"):
            name = _display(r.haplotypes[lab]) + (" (ref)" if lab == "H0" else "")
        else:
            name = lab
        if lab == big:
            who = f"{counts[lab]} samples"
        else:
            who = ", ".join(s for s in samples if r.labels[s] == lab)
        parts.append(f"{name}: {who}")
    return " / ".join(parts)


def _cell(items: list[str]) -> str:
    return ";".join(items) if items else "."


def write_report(cs: CallSet, results: list[LocusResult], path: Path) -> None:
    samples = cs.sample_ids
    disc = [r for r in results if r.discordant]
    n_callable = sum(r.callable_discordant for r in disc)
    effects = cds_effects(cs, disc)
    lines: list[str] = []
    w = lines.append

    w(f"# Discordance report: {cs.set_name}\n")
    w(
        f"{len(samples)} samples, {len(cs.variants)} variant keys, "
        f"**{len(disc)} discordant loci** ({n_callable} where callable samples carry different "
        f"haplotypes; {len(disc) - n_callable} explained only by missing data or ambiguity "
        "codes).\n"
    )
    w(
        "Files: `discordance_summary.tsv` (counts per reason), `discordance.tsv` (one row per "
        "locus), `discordance_pairs.tsv` (one row per locus and sample pair), `variants.tsv` "
        "(per-variant states), `calls.tsv` (per-sample calls).\n"
    )

    w("## Samples\n")
    w(
        "`odd one out`: discordant loci where the sample's haplotype differs from the majority "
        "haplotype. `missing`/`ambiguous`: discordant loci where the sample cannot be called. "
        "`covered` and `ambiguity codes` are genome-wide.\n"
    )
    w("| sample | covered | ambiguity codes (genome) | odd one out | missing | ambiguous |")
    w("|---|--:|--:|--:|--:|--:|")
    for s in cs.samples:
        odd = sum(s.id in r.odd_samples for r in disc)
        miss = sum(r.labels[s.id] == "missing" for r in disc)
        amb = sum(r.labels[s.id] == "ambiguous" for r in disc)
        w(
            f"| {s.id} | {s.summary['covered_fraction']:.1%} | {s.summary['n_ambiguous']} | "
            f"{odd} | {miss} | {amb} |"
        )

    w("\n## Sample pairs\n")
    w(
        "Loci where the two samples carry different sequences (missing data, ambiguity and "
        "equivalent representations excluded); `+N equiv` = loci where they differ only in how "
        "the same sequence is written as variant keys.\n"
    )
    if len(samples) <= REPORT_MAX_MATRIX_SAMPLES:
        short = {s: f"S{i + 1}" for i, s in enumerate(samples)}
        w("| | " + " | ".join(short[s] for s in samples) + " |")
        w("|---|" + "--:|" * len(samples))
        for i, a in enumerate(samples):
            cells = [
                _pair_cell(*_pair_counts(disc, a, b)) if j > i else ""
                for j, b in enumerate(samples)
            ]
            w(f"| {short[a]} | " + " | ".join(cells) + " |")
        w("\n" + "; ".join(f"{short[s]} = {s}" for s in samples) + "\n")
    else:
        pairs = [
            (a, b, *_pair_counts(disc, a, b))
            for i, a in enumerate(samples)
            for b in samples[i + 1 :]
        ]
        pairs = [p for p in pairs if p[2] or p[3]]
        pairs.sort(key=lambda p: (-p[2], -p[3]))
        n = len(pairs)
        shown = f"top {REPORT_TOP_PAIRS} of {n}" if n > REPORT_TOP_PAIRS else f"all {n}"
        w(f"Pairs with the most differences ({shown} pairs with any difference).\n")
        w("| sample a | sample b | differences |")
        w("|---|---|--:|")
        for a, b, real, equiv in pairs[:REPORT_TOP_PAIRS]:
            w(f"| {a} | {b} | {_pair_cell(real, equiv)} |")
        w("")

    w("## Discordant loci\n")
    ranked = sorted(
        disc,
        key=lambda r: (
            LOCUS_RANKING.index(r.primary_reason),
            cs.reference.index(r.locus.contig),
            r.locus.start,
        ),
    )
    capped = (
        f"; first {REPORT_MAX_LOCI} of {len(ranked)} shown, the rest are in `discordance.tsv`"
        if len(ranked) > REPORT_MAX_LOCI
        else ""
    )
    w(f"Most substantive first{capped}.\n")
    w(
        "Who carries what: each haplotype's sequence over the locus (`ref` = reference) and "
        "its samples; the largest group is given as a count.\n"
    )
    w("| locus | genes | primary reason | flags | who carries what | CDS effect |")
    w("|---|---|---|---|---|---|")
    for r in ranked[:REPORT_MAX_LOCI]:
        w(
            f"| {r.locus.id} | {_cell(r.genes)} | {r.primary_reason} | {_cell(r.flags)} | "
            f"{_who_carries_what(r, samples)} | {_cell(effects[r.locus.id])} |"
        )
    if not ranked:
        w("| . | . | . | . | no discordant loci | . |")

    with open(path, "w") as fh:
        fh.write("\n".join(lines) + "\n")
