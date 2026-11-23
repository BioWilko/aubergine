"""Command-line interface: align / call / compare / run."""

from __future__ import annotations

import argparse
import hashlib
import json
import logging
import sys
from pathlib import Path

from . import __version__
from .call import CallError, CallSet, call_set
from .cds import AnnotationError, annotate, parse_gff3, write_cds_tsv
from .compare import CompareError, state_table, write_pairwise, write_truth, write_variants_wide
from .discord import analyse, write_outputs
from .fasta import FastaError, Reference, has_content, load_alignment, read_fasta
from .haplotype import build_sample_events
from .io import (
    SUMMARY_HEADER,
    InputError,
    effective_thresholds,
    read_callset,
    summary_rows,
    write_calls_tsv,
    write_callset,
    write_summary_tsv,
    write_tsv,
)
from .mafft import MafftError, align, mafft_version
from .params import FIXED_BEHAVIOURS, THRESHOLDS
from .report import write_report
from .vcf import write_vcf

log = logging.getLogger("aubergine")

FASTA_SUFFIXES = (".fa", ".fasta", ".fna")
ERRORS = (FastaError, CallError, AnnotationError, CompareError, InputError, MafftError, OSError)


class UsageError(ValueError):
    pass


# --------------------------------------------------------------------------- helpers


def _threshold_type(name: str):
    minimum = THRESHOLDS[name][1]

    def parse(text: str) -> int:
        try:
            value = int(text)
        except ValueError:
            raise argparse.ArgumentTypeError(f"must be an integer, got {text!r}") from None
        if value < minimum:
            raise argparse.ArgumentTypeError(f"must be >= {minimum}, got {value}")
        return value

    return parse


def _add_thresholds(p: argparse.ArgumentParser, with_defaults: bool) -> None:
    g = p.add_argument_group("thresholds")
    help_text = {
        "hp_min": "an event is 'homopolymer' if its hp_len >= N",
        "locus_join_gap": "slide spans separated by <= N reference bases join one locus",
        "end_dist": "distance (bp) from a sequence end that sets the near_end flag",
    }
    for name, (default, minimum) in THRESHOLDS.items():
        suffix = f"default {default}" if with_defaults else "default: value in variants.json"
        g.add_argument(
            "--" + name.replace("_", "-"),
            dest=name,
            type=_threshold_type(name),
            default=default if with_defaults else None,
            metavar="N",
            help=f"{help_text[name]} (>= {minimum}; {suffix})",
        )


def _thresholds(args) -> dict:
    return {name: getattr(args, name) for name in THRESHOLDS}


def load_reference(path: Path) -> Reference:
    ref = Reference.from_fasta(path)
    if len(ref.ids) > 1:
        raise UsageError(
            f"{path}: multi-record references are not supported yet "
            f"({len(ref.ids)} records: {', '.join(ref.ids[:5])})"
        )
    return ref


def set_name_from(path: Path) -> str:
    name = Path(path).name
    if name.endswith(".gz"):
        name = name[:-3]
    for suffix in FASTA_SUFFIXES:
        if name.endswith(suffix):
            name = name[: -len(suffix)]
            break
    if name.endswith(".aligned"):
        name = name[: -len(".aligned")]
    return name


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


# --------------------------------------------------------------------------- stages

CALL_OUTPUTS = ["variants.json", "calls.tsv", "vcf", "summary.tsv", "cds.tsv"]
COMPARE_OUTPUTS = {
    "variants": "variants.tsv",
    "pairwise": "pairwise.tsv",
    "discordance_report": "discordance_report.md",
    "discordance": "discordance.tsv",
    "discordance_pairs": "discordance_pairs.tsv",
    "discordance_summary": "discordance_summary.tsv",
    "truth_metrics": "truth_metrics.tsv",
    "truth_detail": "truth_detail.tsv",
}


def run_call(
    reference: Reference,
    msa: Path,
    set_name: str,
    thresholds: dict,
    annotation: Path | None,
    out: dict[str, Path],
) -> CallSet:
    alignment = load_alignment(msa, reference)
    cs = call_set(alignment, reference, set_name, thresholds)
    cs.alignments = [(c, str(msa)) for c in reference.ids]
    if annotation is not None:
        cds_list = parse_gff3(annotation, reference)
        annotate(cs, cds_list, build_sample_events(cs), str(annotation))
        write_cds_tsv(cs, out["cds.tsv"])
    sample_events = build_sample_events(cs)
    write_callset(cs, out["variants.json"])
    write_calls_tsv(cs, out["calls.tsv"])
    if out.get("summary.tsv") is not None:
        write_summary_tsv(cs, out["summary.tsv"])
    write_vcf(cs, state_table(cs, sample_events), sample_events, out["vcf"])
    log.info("%s: %d samples, %d variants", set_name, len(cs.samples), len(cs.variants))
    return cs


def run_compare(cs: CallSet, thresholds: dict, truth: str | None, out: dict[str, Path]) -> None:
    sample_events = build_sample_events(cs)
    states = state_table(cs, sample_events)
    loci = analyse(cs, sample_events, thresholds)
    write_variants_wide(cs, states, loci, out["variants"])
    write_pairwise(cs, states, out["pairwise"])
    write_outputs(cs, loci, out)
    write_report(cs, loci, out["discordance_report"])
    if truth is not None:
        write_truth(cs, states, loci, truth, out)
    n_disc = sum(1 for r in loci if r.discordant)
    log.info("%s: %d loci, %d discordant", cs.set_name, len(loci), n_disc)


# --------------------------------------------------------------------------- commands


def cmd_align(args) -> None:
    reference = load_reference(args.reference)
    align(reference, args.set, args.output, args.threads)
    log.info("wrote %s", args.output)


def cmd_call(args) -> None:
    reference = load_reference(args.reference)
    prefix = str(args.output)
    out = {name: Path(f"{prefix}.{name}") for name in CALL_OUTPUTS}
    run_call(
        reference,
        args.msa,
        args.set_name or set_name_from(args.msa),
        _thresholds(args),
        args.annotation,
        out,
    )


def cmd_compare(args) -> None:
    cs = read_callset(args.variants)
    thresholds = effective_thresholds(cs, _thresholds(args), log)
    prefix = str(args.output)
    out = {name: Path(f"{prefix}.{fn}") for name, fn in COMPARE_OUTPUTS.items()}
    run_compare(cs, thresholds, args.truth, out)


def discover_inputs(inputs: list[Path]) -> list[Path]:
    files: list[Path] = []
    for p in inputs:
        if p.is_dir():
            found = sorted(
                f
                for f in p.iterdir()
                if f.is_file()
                and any(f.name.endswith(s) or f.name.endswith(s + ".gz") for s in FASTA_SUFFIXES)
            )
            if not found:
                log.warning("%s: no FASTA files found", p)
            files += found
        elif p.is_file():
            files.append(p)
        else:
            raise UsageError(f"input not found: {p}")
    if not files:
        raise UsageError("no input FASTA files")
    names = [set_name_from(f) for f in files]
    dups = sorted({n for n in names if names.count(n) > 1})
    if dups:
        raise UsageError(f"several inputs give the same set name: {', '.join(dups)}")
    return files


def cmd_run(args, argv: list[str]) -> None:
    reference = load_reference(args.reference)
    thresholds = _thresholds(args)
    files = discover_inputs(args.inputs)
    outdir = Path(args.output)
    outdir.mkdir(parents=True, exist_ok=True)
    run_rows = []
    skipped = []
    mafft_cmd: list[str] = []
    for f in files:
        name = set_name_from(f)
        records = read_fasta(f)
        n_content = sum(has_content(r.seq) for r in records)
        if n_content < 2:
            log.warning(
                "set %s: skipped, only %d of %d records have bases other than N",
                name,
                n_content,
                len(records),
            )
            skipped.append(
                {"set": name, "input": str(f), "records": len(records), "with_content": n_content}
            )
            continue
        setdir = outdir / name
        setdir.mkdir(exist_ok=True)
        log.info("set %s: aligning %s", name, f)
        aligned = setdir / "aligned.fasta"
        mafft_cmd = align(reference, f, aligned, args.threads)
        call_out = {n: setdir / n for n in CALL_OUTPUTS}
        call_out["vcf"] = setdir / "variants.vcf"
        del call_out["summary.tsv"]  # written once, at the top level
        cs = run_call(reference, aligned, name, thresholds, args.annotation, call_out)
        run_compare(cs, thresholds, None, {n: setdir / fn for n, fn in COMPARE_OUTPUTS.items()})
        run_rows += summary_rows(cs)
    write_tsv(outdir / "summary.tsv", SUMMARY_HEADER, run_rows)
    info = {
        "tool_version": __version__,
        "command_line": ["aubergine", *argv],
        "mafft_version": mafft_version(),
        "mafft_command": mafft_cmd,
        "reference": {
            "path": str(args.reference),
            "sha256": _sha256(args.reference),
            "contigs": [{"id": c.id, "sha256": c.sha256} for c in reference.contigs.values()],
        },
        "annotation": str(args.annotation) if args.annotation else None,
        "inputs": [str(f) for f in files],
        "skipped": skipped,
        "thresholds": thresholds,
        "fixed_behaviours": FIXED_BEHAVIOURS,
    }
    with open(outdir / "run_info.json", "w") as fh:
        json.dump(info, fh, indent=1)
        fh.write("\n")


# --------------------------------------------------------------------------- parser


def build_parser() -> argparse.ArgumentParser:
    common = argparse.ArgumentParser(add_help=False)
    common.add_argument("-v", "--verbose", action="store_true", help="log debug messages")
    common.add_argument("-q", "--quiet", action="store_true", help="log only errors")
    common.add_argument("--version", action="version", version=f"aubergine {__version__}")

    parser = argparse.ArgumentParser(
        prog="aubergine",
        description="MSA-based variant calling and disagreement analysis for consensus genomes.",
        parents=[common],
    )
    sub = parser.add_subparsers(dest="command", required=True)

    p = sub.add_parser("align", parents=[common], help="align a set of sequences to the reference")
    p.add_argument("-r", "--reference", type=Path, required=True)
    p.add_argument("-o", "--output", type=Path, required=True, help="aligned FASTA to write")
    p.add_argument("--threads", type=int, default=1)
    p.add_argument("set", type=Path, help="FASTA of sequences to align")

    p = sub.add_parser("call", parents=[common], help="call variants from a reference-anchored MSA")
    p.add_argument("-r", "--reference", type=Path, required=True)
    p.add_argument("-o", "--output", required=True, help="output prefix")
    p.add_argument("--annotation", type=Path, help="GFF3 with CDS features")
    p.add_argument("--set-name", help="set name (default: from the MSA file name)")
    _add_thresholds(p, with_defaults=True)
    p.add_argument("msa", type=Path)

    p = sub.add_parser(
        "compare", parents=[common], help="compare samples and explain disagreements"
    )
    p.add_argument("-o", "--output", required=True, help="output prefix")
    p.add_argument("--truth", metavar="SAMPLE_ID", help="sample in the set to treat as truth")
    _add_thresholds(p, with_defaults=False)
    p.add_argument("variants", type=Path, help="variants.json written by `aubergine call`")

    p = sub.add_parser("run", parents=[common], help="align, call and compare one or more sets")
    p.add_argument("-r", "--reference", type=Path, required=True)
    p.add_argument("-o", "--output", type=Path, required=True, help="output directory")
    p.add_argument("--threads", type=int, default=1)
    p.add_argument("--annotation", type=Path, help="GFF3 with CDS features")
    _add_thresholds(p, with_defaults=True)
    p.add_argument("inputs", type=Path, nargs="+", help="FASTA files and/or directories")
    return parser


def _configure_logging(level: int) -> None:
    for h in list(log.handlers):
        log.removeHandler(h)
    handler = logging.StreamHandler(sys.stderr)
    handler.setFormatter(logging.Formatter("%(levelname)s: %(message)s"))
    log.addHandler(handler)
    log.setLevel(level)


def main(argv: list[str] | None = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    args = build_parser().parse_args(argv)
    _configure_logging(
        logging.DEBUG if args.verbose else logging.ERROR if args.quiet else logging.INFO
    )
    try:
        if args.command == "align":
            cmd_align(args)
        elif args.command == "call":
            cmd_call(args)
        elif args.command == "compare":
            cmd_compare(args)
        else:
            cmd_run(args, argv)
    except (UsageError, *ERRORS) as e:
        log.error("%s", e)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
