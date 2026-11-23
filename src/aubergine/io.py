"""Canonical JSON (variants.json) read/write and TSV writers."""

from __future__ import annotations

import json
from collections.abc import Iterable, Sequence
from pathlib import Path

from .call import Ambiguous, Call, CallSet, SampleResult, Variant
from .fasta import Contig, Reference
from .params import THRESHOLDS

FORMAT_VERSION = 1


class InputError(ValueError):
    pass


def _variant_to_json(v: Variant) -> dict:
    return {
        "key": v.key,
        "contig": v.contig,
        "type": v.type,
        "pos": v.pos,
        "ref": v.ref,
        "alt": v.alt,
        "start": v.start,
        "end": v.end,
        "length": v.length,
        "len_mod3": v.len_mod3,
        "hp_len": v.hp_len,
        "hp_run": list(v.hp_run) if v.hp_run else None,
        "slide": list(v.slide),
        "cds": v.cds,
    }


def _variant_from_json(d: dict) -> Variant:
    ins_seq = ""
    if d["type"] == "ins":
        ins_seq = d["alt"][1:] if d["start"] >= 1 else d["alt"][:-1]
    return Variant(
        contig=d["contig"],
        type=d["type"],
        start=d["start"],
        end=d["end"],
        pos=d["pos"],
        ref=d["ref"],
        alt=d["alt"],
        ins_seq=ins_seq,
        hp_len=d["hp_len"],
        hp_run=tuple(d["hp_run"]) if d["hp_run"] else None,
        slide=tuple(d["slide"]),
        cds=d.get("cds", []),
    )


def callset_to_json(cs: CallSet) -> dict:
    return {
        "format_version": FORMAT_VERSION,
        "tool_version": cs.tool_version,
        "params": cs.params,
        "reference": {
            "path": cs.reference.path,
            "contigs": [
                {"id": c.id, "length": c.length, "sha256": c.sha256, "sequence": c.seq}
                for c in cs.reference.contigs.values()
            ],
        },
        "set": cs.set_name,
        "alignments": [{"contig": c, "path": p} for c, p in cs.alignments],
        "annotation": cs.annotation,
        "samples": [
            {
                "id": s.id,
                "records": s.records,
                "summary": s.summary,
                "covered": [list(x) for x in s.covered],
                "missing": [list(x) for x in s.missing],
                "missing_insertions": [list(x) for x in s.missing_insertions],
                "overhang_bases": s.overhang_bases,
            }
            for s in cs.samples
        ],
        "variants": [_variant_to_json(cs.variants[k]) for k in cs.sorted_keys()],
        "calls": [
            {
                "sample": c.sample,
                "key": c.key,
                "seq_start": c.seq_start,
                "ins_seq_raw": c.ins_seq_raw,
                "flags": c.flags,
                "frameshift": c.frameshift,
            }
            for c in cs.calls
        ],
        "ambiguous": [
            {"sample": a.sample, "contig": a.contig, "pos": a.pos, "ref": a.ref, "code": a.code}
            for a in cs.ambiguous
        ],
        "cds_status": cs.cds_status,
    }


def write_callset(cs: CallSet, path: Path) -> None:
    with open(path, "w") as fh:
        json.dump(callset_to_json(cs), fh, indent=1)
        fh.write("\n")


def read_callset(path: Path) -> CallSet:
    try:
        with open(path) as fh:
            d = json.load(fh)
    except (OSError, json.JSONDecodeError) as e:
        raise InputError(f"{path}: cannot read variants JSON: {e}") from e
    if d.get("format_version") != FORMAT_VERSION:
        raise InputError(f"{path}: unsupported format_version {d.get('format_version')!r}")
    ref = Reference(
        {c["id"]: Contig(c["id"], c["sequence"]) for c in d["reference"]["contigs"]},
        d["reference"].get("path"),
    )
    cs = CallSet(reference=ref, set_name=d["set"], params=d.get("params", {}))
    cs.tool_version = d["tool_version"]
    cs.alignments = [(a["contig"], a["path"]) for a in d.get("alignments", [])]
    cs.annotation = d.get("annotation")
    for s in d["samples"]:
        cs.samples.append(
            SampleResult(
                id=s["id"],
                records=s["records"],
                covered=[tuple(x) for x in s["covered"]],
                missing=[tuple(x) for x in s["missing"]],
                missing_insertions=[tuple(x) for x in s["missing_insertions"]],
                overhang_bases=s["overhang_bases"],
                summary=s["summary"],
            )
        )
    for v in d["variants"]:
        var = _variant_from_json(v)
        cs.variants[var.key] = var
    cs.calls = [Call(**c) for c in d["calls"]]
    cs.ambiguous = [Ambiguous(**a) for a in d["ambiguous"]]
    cs.cds_status = d.get("cds_status", [])
    return cs


def effective_thresholds(cs: CallSet, overrides: dict, log) -> dict:
    """CLI value if given, else the value recorded in variants.json (missing -> error)."""
    out = {}
    for name in THRESHOLDS:
        recorded = cs.params.get(name)
        given = overrides.get(name)
        if recorded is None:
            raise InputError(f"variants.json has no recorded value for threshold {name!r}")
        if given is not None and given != recorded:
            log.warning(
                "%s: using command-line value %s instead of %s recorded in variants.json",
                name,
                given,
                recorded,
            )
        out[name] = given if given is not None else recorded
    log.info("effective thresholds: %s", ", ".join(f"{k}={v}" for k, v in out.items()))
    return out


def fmt(value) -> str:
    if value is None or value == "":
        return "."
    if isinstance(value, bool):
        return "yes" if value else "no"
    if isinstance(value, float):
        return f"{value:.4f}"
    if isinstance(value, list | tuple):
        return ";".join(fmt(v) for v in value) if value else "."
    return str(value)


def write_tsv(path: Path, header: Sequence[str], rows: Iterable[Sequence]) -> None:
    with open(path, "w") as fh:
        fh.write("\t".join(header) + "\n")
        for row in rows:
            fh.write("\t".join(fmt(v) for v in row) + "\n")


CALLS_TSV_HEADER = ["sample", "key", "seq_start", "frameshift", "flags"]


def write_calls_tsv(cs: CallSet, path: Path) -> None:
    """One row per call plus one per ambiguity code; variant attributes live in variants.tsv."""
    order = {s: i for i, s in enumerate(cs.sample_ids)}
    rows = []
    for c in cs.calls:
        sk = (order[c.sample], *cs.sort_key(c.key))
        rows.append((sk, [c.sample, c.key, c.seq_start, c.frameshift, c.flags]))
    for a in cs.ambiguous:
        sk = (order[a.sample], cs.reference.index(a.contig), a.pos, "ambiguous", a.ref, a.code)
        key = f"ambiguous:{a.contig}:{a.pos}:{a.ref}:{a.code}"
        rows.append((sk, [a.sample, key, None, None, None]))
    rows.sort(key=lambda r: r[0])
    write_tsv(path, CALLS_TSV_HEADER, (r[1] for r in rows))


SUMMARY_HEADER = [
    "set",
    "sample",
    "n_snp",
    "n_ins",
    "n_del",
    "n_ambiguous",
    "missing_bases",
    "covered_fraction",
    "overhang_bases",
    "frameshifted_cds",
    "restored_cds",
    "premature_stop_cds",
    "stop_lost_cds",
]


def summary_rows(cs: CallSet) -> list[list]:
    return [[cs.set_name, s.id] + [s.summary[h] for h in SUMMARY_HEADER[2:]] for s in cs.samples]


def write_summary_tsv(cs: CallSet, path: Path) -> None:
    write_tsv(path, SUMMARY_HEADER, summary_rows(cs))
