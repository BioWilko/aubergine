"""GFF3 CDS parsing, frameshift detection and simple premature stop detection."""

from __future__ import annotations

import re
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from urllib.parse import unquote

from .call import CallSet, summarise
from .fasta import ACGT, Reference
from .haplotype import INS, Base, SampleEvents, apply
from .io import write_tsv

STOP_CODONS = {"TAA", "TAG", "TGA"}  # NCBI table 1
_COMPLEMENT = str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN")


class AnnotationError(ValueError):
    pass


@dataclass
class CDS:
    id: str
    gene: str
    contig: str
    strand: str
    segments: list[tuple[int, int]]  # transcription order
    transl_table: str = "1"
    offsets: list[int] = field(init=False)
    coords: dict[int, list[int]] = field(init=False)  # ref pos -> CDS coords

    def __post_init__(self):
        self.offsets = []
        self.coords = defaultdict(list)
        o = 0
        for s, e in self.segments:
            self.offsets.append(o)
            for q in range(s, e + 1):
                self.coords[q].append(self._coord(q, s, e, o))
            o += e - s + 1

    def _coord(self, q, s, e, o):
        return o + (q - s + 1 if self.strand == "+" else e - q + 1)

    @property
    def length(self) -> int:
        return sum(e - s + 1 for s, e in self.segments)

    @property
    def end_pos(self) -> int:
        """Reference position of the last coding base in transcription order."""
        s, e = self.segments[-1]
        return e if self.strand == "+" else s

    def overlaps(self, lo: int, hi: int) -> bool:
        return any(s <= hi and e >= lo for s, e in self.segments)

    def deletion_delta(self, s: int, e: int) -> tuple[int, bool]:
        """(-number of CDS coordinates removed, overlaps first/last codon)."""
        removed = [c for q in range(s, e + 1) for c in self.coords.get(q, [])]
        L = self.length
        boundary = any(c <= 3 or c >= L - 2 for c in removed)
        return -len(removed), boundary

    def insertion_coding(self, p: int) -> bool:
        """True if p and p+1 are consecutive CDS coordinates (in transcription order)."""
        before, after = (p, p + 1) if self.strand == "+" else (p + 1, p)
        nxt = set(self.coords.get(after, []))
        return any(c + 1 in nxt for c in self.coords.get(before, []))

    def to_json(self) -> dict:
        return {
            "id": self.id,
            "gene": self.gene,
            "contig": self.contig,
            "strand": self.strand,
            "segments": [list(x) for x in self.segments],
            "transl_table": self.transl_table,
        }


def _attrs(col: str) -> dict[str, str]:
    out = {}
    for part in col.strip().split(";"):
        if "=" in part:
            k, v = part.split("=", 1)
            out[k.strip()] = unquote(v.strip())
    return out


def _sanitise(label: str) -> str:
    return re.sub(r"[\s,;:|=]", "_", label)


def parse_gff3(path: Path, reference: Reference) -> list[CDS]:
    rows = defaultdict(list)
    meta = {}
    bad_seqids = set()
    with open(path) as fh:
        for n, line in enumerate(fh, 1):
            if line.startswith("##FASTA"):
                break
            if line.startswith("#") or not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) != 9:
                raise AnnotationError(f"{path}:{n}: expected 9 tab-separated columns")
            if cols[2] != "CDS":
                continue
            seqid, start, end, strand = cols[0], int(cols[3]), int(cols[4]), cols[6]
            if seqid not in reference.contigs:
                bad_seqids.add(seqid)
                continue
            if strand not in "+-" or not strand:
                raise AnnotationError(f"{path}:{n}: CDS strand must be + or -")
            a = _attrs(cols[8])
            cid = a.get("ID") or a.get("Parent") or f"cds{n}"
            rows[cid].append((start, end))
            prev = meta.setdefault(
                cid,
                {
                    "contig": seqid,
                    "strand": strand,
                    "label": a.get("gene") or a.get("Name") or cid,
                    "transl_table": a.get("transl_table", "1"),
                },
            )
            if prev["contig"] != seqid or prev["strand"] != strand:
                raise AnnotationError(f"{path}: CDS {cid!r} spans contigs or strands")
    if bad_seqids:
        raise AnnotationError(
            f"{path}: GFF seqid(s) not in reference: {', '.join(sorted(bad_seqids))} "
            f"(reference contigs: {', '.join(reference.ids)})"
        )
    labels = [m["label"] for m in meta.values()]
    out = []
    for cid, m in meta.items():
        label = m["label"] if labels.count(m["label"]) == 1 else f"{m['label']}[{cid}]"
        segs = sorted(rows[cid])
        if m["strand"] == "-":
            segs = sorted(segs, key=lambda x: (-x[1], -x[0]))
        for s, e in segs:
            if not 1 <= s <= e <= reference.length(m["contig"]):
                raise AnnotationError(f"{path}: CDS {cid!r} segment {s}-{e} is outside the contig")
        out.append(CDS(cid, _sanitise(label), m["contig"], m["strand"], segs, m["transl_table"]))
    out.sort(key=lambda c: (reference.index(c.contig), min(s for s, _ in c.segments), c.id))
    return out


# --------------------------------------------------------------------------- projection


def project(ref_seq: str, cds: CDS, se: SampleEvents | None) -> list[tuple[Base, int]]:
    """Sample (or reference, if se is None) sequence of a CDS, with CDS coordinates.

    Returns (base, cds_coord) pairs in transcription order; an inserted base carries the
    coordinate of the base preceding it in transcription order.
    """
    out = []
    for (s, e), o in zip(cds.segments, cds.offsets, strict=True):
        if se is None:
            bases = [Base(ref_seq[q - 1], "ref", q) for q in range(s, e + 1)]
        else:
            bases = apply(ref_seq, se, s, e, missing_as_n=True, insertion_ok=cds.insertion_coding)
        if cds.strand == "-":
            bases = [Base(b.base.translate(_COMPLEMENT), b.kind, b.pos) for b in reversed(bases)]
        for b in bases:
            if b.kind == INS:
                anchor = b.pos if cds.strand == "+" else b.pos + 1
            else:
                anchor = b.pos
            out.append((b, cds._coord(anchor, s, e, o)))
    return out


def _codons(seq: list[tuple[Base, int]]):
    for i in range(0, len(seq) - 2, 3):
        yield i // 3, seq[i : i + 3]


def validate_reference(ref_seq: str, cds: CDS) -> str | None:
    if cds.transl_table != "1":
        return "unsupported_transl_table"
    if cds.length % 3:
        return "ref_length_not_mod3"
    codons = ["".join(b.base for b, _ in c) for _, c in _codons(project(ref_seq, cds, None))]
    if codons[-1] not in STOP_CODONS:
        return "ref_no_terminal_stop"
    if any(c in STOP_CODONS for c in codons[:-1]):
        return "ref_internal_stop"
    return None


def stop_status(
    ref_seq: str, cds: CDS, se: SampleEvents, note: str | None, out_of_frame, missing_ins: bool
) -> dict:
    Lc = cds.length
    res = {
        "stop_status": "not_evaluated",
        "incomplete": False,
        "note": note,
        "stop_aa_pos": None,
        "protein_len": None,
        "ref_protein_len": Lc // 3 - 1 if not Lc % 3 else None,
        "stop_ref_pos": None,
        "stop_codon": None,
        "stop_frameshift": None,
        "trailing_bases": None,
    }
    if note:
        return res
    seq = project(ref_seq, cds, se)
    incomplete = missing_ins
    first_stop = None
    n_codons = 0
    for i, codon in _codons(seq):
        n_codons = i + 1
        bases = "".join(b.base for b, _ in codon)
        if any(c not in ACGT for c in bases):
            incomplete = True
            continue
        if bases in STOP_CODONS:
            first_stop = (i, codon)
            break
    res["trailing_bases"] = len(seq) % 3
    res["incomplete"] = incomplete
    if first_stop is None:
        res["stop_status"] = "unknown" if incomplete else "stop_lost"
        res["protein_len"] = n_codons
        if res["stop_status"] == "stop_lost":
            res["stop_frameshift"] = bool(out_of_frame) and out_of_frame[-1][2]
        return res
    i, codon = first_stop
    base, coord = codon[0]
    premature = coord <= Lc - 3 if base.kind == INS else coord < Lc - 2
    res["stop_status"] = "premature_stop" if premature else "normal"
    res["stop_aa_pos"] = i + 1
    res["protein_len"] = i
    res["stop_ref_pos"] = f"{cds.contig}:{base.pos}"
    # reference span of all three stop-codon bases (an inserted base spans its anchor p..p+1)
    res["stop_codon"] = [
        min(b.pos for b, _ in codon),
        max(b.pos + (b.kind == INS) for b, _ in codon),
    ]
    if premature:
        res["stop_frameshift"] = any(lo <= base.pos <= hi for lo, hi, _ in out_of_frame)
    return res


# --------------------------------------------------------------------------- annotation


def annotate(
    cs: CallSet,
    cds_list: list[CDS],
    sample_events: dict[tuple[str, str], SampleEvents],
    path: str,
) -> None:
    """Add CDS effects to variants and calls, and per sample x CDS status, in place."""
    cs.annotation = {"path": path, "cds": [c.to_json() for c in cds_list]}

    boundary_keys = set()
    for v in cs.variants.values():
        v.cds = []
        if not v.is_indel:
            continue
        for c in cds_list:
            if c.contig != v.contig:
                continue
            if v.type == "del":
                delta, boundary = c.deletion_delta(v.start, v.end)
                if boundary and delta:
                    boundary_keys.add(v.key)
            else:
                delta = len(v.ins_seq) if c.insertion_coding(v.start) else 0
            if delta:
                v.cds.append({"gene": c.gene, "delta_cds": delta})

    calls_by_sample = cs.calls_by_sample()
    ref_notes = {c.id: validate_reference(cs.reference.seq(c.contig), c) for c in cds_list}
    cs.cds_status = []
    for s in cs.samples:
        calls = calls_by_sample[s.id]
        for c in cds_list:
            ref_seq = cs.reference.seq(c.contig)
            se = sample_events[(s.id, c.contig)]
            coding = []
            for call in calls:
                v = cs.variants[call.key]
                for x in v.cds:
                    if x["gene"] == c.gene:
                        coding.append((v, x["delta_cds"]))
            net = sum(d for _, d in coding)
            for call in calls:
                v = cs.variants[call.key]
                vals = []
                for x in v.cds:
                    if x["gene"] != c.gene:
                        continue
                    if x["delta_cds"] % 3 == 0:
                        vals.append("inframe")
                    else:
                        vals.append("fs_compensated" if net % 3 == 0 else "fs")
                if vals:
                    entry = f"{c.gene}:{vals[0]}"
                    call.frameshift = (
                        entry if call.frameshift == "." else f"{call.frameshift}|{entry}"
                    )
            if not coding:
                status = "intact"
            elif net % 3:
                status = "frameshift"
            elif any(d % 3 for _, d in coding):
                status = "restored"
            else:
                status = "inframe"
            oof = _out_of_frame(c, coding)
            missing_ins = any(c.insertion_coding(p) for p in se.missing_insertions)
            incomplete = any(c.overlaps(lo, hi) for lo, hi in se.missing) or missing_ins
            stop = stop_status(ref_seq, c, se, ref_notes[c.id], oof, missing_ins)
            cs.cds_status.append(
                {
                    "sample": s.id,
                    "contig": c.contig,
                    "gene": c.gene,
                    "id": c.id,
                    "strand": c.strand,
                    "net_delta": net,
                    "status": status,
                    "out_of_frame": [[lo, hi] for lo, hi, _ in oof],
                    "incomplete": incomplete,
                    "stop": stop,
                }
            )

    for call in cs.calls:
        if call.key in boundary_keys and "cds_boundary" not in call.flags:
            call.flags.append("cds_boundary")
    summarise(cs)


def _out_of_frame(cds: CDS, coding) -> list[tuple[int, int, bool]]:
    """Reference spans (lo, hi, runs_to_cds_end) where the sample's frame is shifted."""

    def order(item):
        v, _ = item
        coords = [x for q in range(v.start, v.end + 1) for x in cds.coords.get(q, [])]
        return min(coords) if coords else 0

    spans = []
    cum = 0
    open_pos = None
    for v, d in sorted(coding, key=order):
        prev = cum
        cum += d
        here = v.start if cds.strand == "+" else v.end
        if prev % 3 == 0 and cum % 3:
            open_pos = here
        elif prev % 3 and cum % 3 == 0 and open_pos is not None:
            spans.append((min(open_pos, here), max(open_pos, here), False))
            open_pos = None
    if open_pos is not None:
        spans.append((min(open_pos, cds.end_pos), max(open_pos, cds.end_pos), True))
    return spans


CDS_TSV_HEADER = [
    "sample",
    "contig",
    "gene",
    "strand",
    "net_delta",
    "status",
    "out_of_frame",
    "stop_status",
    "stop_frameshift",
    "stop_aa_pos",
    "protein_len",
    "ref_protein_len",
    "stop_ref_pos",
    "incomplete",
    "stop_incomplete",
    "note",
]


def write_cds_tsv(cs: CallSet, path: Path) -> None:
    rows = []
    for st in cs.cds_status:
        stop = st["stop"]
        rows.append(
            [
                st["sample"],
                st["contig"],
                st["gene"],
                st["strand"],
                st["net_delta"],
                st["status"],
                [f"{lo}-{hi}" for lo, hi in st["out_of_frame"]],
                stop["stop_status"],
                stop["stop_frameshift"],
                stop["stop_aa_pos"],
                stop["protein_len"],
                stop["ref_protein_len"],
                stop["stop_ref_pos"],
                st["incomplete"],
                stop["incomplete"],
                stop["note"],
            ]
        )
    write_tsv(path, CDS_TSV_HEADER, rows)
