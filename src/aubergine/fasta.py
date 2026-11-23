"""FASTA reading/writing, validation, and the contig-aware reference and MSA models."""

from __future__ import annotations

import gzip
import hashlib
from collections.abc import Iterable
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

ACGT = "ACGT"
IUPAC = "RYSWKMBDHVN"
ALPHABET = set(ACGT + IUPAC + "-")
GAP = ord("-")
N = ord("N")

_ACGT_CODES = np.frombuffer(ACGT.encode(), dtype=np.uint8)


class FastaError(ValueError):
    pass


@dataclass
class Record:
    id: str
    description: str
    seq: str


def _open(path: Path):
    path = Path(path)
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return open(path)


def normalise_seq(seq: str) -> str:
    """Upper-case, map U->T and '.'->'-', strip whitespace."""
    return "".join(seq.split()).upper().replace("U", "T").replace(".", "-")


def has_content(seq: str) -> bool:
    """True if the sequence has at least one A, C, G or T (not only N, ambiguity codes or gaps)."""
    return any(b in seq for b in ACGT)


def read_fasta(path: Path, allow_gaps: bool = True) -> list[Record]:
    """Read a (optionally gzipped) FASTA file, normalising and validating sequences."""
    records: list[Record] = []
    header = None
    chunks: list[str] = []

    def flush():
        if header is None:
            return
        rid = header.split()[0] if header.split() else ""
        if not rid:
            raise FastaError(f"{path}: record with empty ID")
        seq = normalise_seq("".join(chunks))
        bad = set(seq) - ALPHABET
        if not allow_gaps:
            bad |= set(seq) & {"-"}
        if bad:
            raise FastaError(
                f"{path}: record {rid!r} contains invalid characters: {''.join(sorted(bad))}"
            )
        records.append(Record(rid, header, seq))

    with _open(path) as fh:
        for line in fh:
            line = line.rstrip("\r\n")
            if line.startswith(">"):
                flush()
                header = line[1:].strip()
                chunks = []
            elif line.strip():
                if header is None:
                    raise FastaError(f"{path}: sequence data before first header")
                chunks.append(line)
        flush()

    ids = [r.id for r in records]
    dups = sorted({i for i in ids if ids.count(i) > 1})
    if dups:
        raise FastaError(f"{path}: duplicate record IDs: {', '.join(dups)}")
    return records


def write_fasta(path: Path, records: Iterable[Record], width: int = 0) -> None:
    with open(path, "w") as fh:
        for r in records:
            fh.write(f">{r.description or r.id}\n")
            if width:
                for i in range(0, len(r.seq), width):
                    fh.write(r.seq[i : i + width] + "\n")
            else:
                fh.write(r.seq + "\n")


@dataclass
class Contig:
    id: str
    seq: str

    @property
    def length(self) -> int:
        return len(self.seq)

    @property
    def sha256(self) -> str:
        return hashlib.sha256(self.seq.encode()).hexdigest()


@dataclass
class Reference:
    """Ordered mapping contig_id -> sequence (FASTA order)."""

    contigs: dict[str, Contig] = field(default_factory=dict)
    path: str | None = None

    @classmethod
    def from_fasta(cls, path: Path) -> Reference:
        records = read_fasta(path, allow_gaps=False)
        if not records:
            raise FastaError(f"{path}: no records in reference")
        return cls({r.id: Contig(r.id, r.seq) for r in records}, str(path))

    @property
    def ids(self) -> list[str]:
        return list(self.contigs)

    def index(self, contig: str) -> int:
        return self.ids.index(contig)

    def seq(self, contig: str) -> str:
        return self.contigs[contig].seq

    def length(self, contig: str) -> int:
        return self.contigs[contig].length


@dataclass
class ContigAlignment:
    """One contig's MSA: the reference row plus sample rows, with coordinate maps."""

    contig: str
    ref_row: np.ndarray  # uint8, alignment columns
    sample_ids: list[str]
    rows: np.ndarray  # uint8, shape (n_samples, n_columns)
    ref_seq: str  # ungapped reference
    ref_pos: np.ndarray = field(init=False)  # 1-based last ref base at or before column

    def __post_init__(self):
        self.ref_pos = np.cumsum(self.ref_row != GAP).astype(np.int64)

    @property
    def n_columns(self) -> int:
        return int(self.ref_row.shape[0])

    @property
    def ref_is_acgt(self) -> np.ndarray:
        return np.isin(self.ref_row, _ACGT_CODES)

    def covered_columns(self, i: int) -> tuple[int, int] | None:
        """First/last column where sample i has a character at a reference base."""
        row = self.rows[i]
        mask = (self.ref_row != GAP) & (row != GAP)
        idx = np.flatnonzero(mask)
        if idx.size == 0:
            return None
        return int(idx[0]), int(idx[-1])


SetAlignment = dict  # contig_id -> ContigAlignment


def _to_array(seq: str) -> np.ndarray:
    return np.frombuffer(seq.encode(), dtype=np.uint8).copy()


def load_alignment(path: Path, reference: Reference) -> SetAlignment:
    """Load an MSA containing the reference record(s) and samples.

    The reference row is found by the reference record ID; its ungapped sequence must
    equal the reference. All other records are samples of that contig.
    """
    records = read_fasta(path)
    if not records:
        raise FastaError(f"{path}: alignment is empty")
    lengths = {len(r.seq) for r in records}
    if len(lengths) > 1:
        raise FastaError(f"{path}: alignment rows have different lengths: {sorted(lengths)}")
    by_id = {r.id: r for r in records}

    out: SetAlignment = {}
    ref_ids = set(reference.ids)
    for contig in reference.ids:
        if contig not in by_id:
            raise FastaError(
                f"{path}: reference record {contig!r} not found in alignment "
                f"(IDs present: {', '.join(list(by_id)[:10])}{'...' if len(by_id) > 10 else ''})"
            )
        ref_row = by_id[contig].seq
        if ref_row.replace("-", "") != reference.seq(contig):
            raise FastaError(
                f"{path}: ungapped reference row {contig!r} does not match the reference sequence"
            )
        samples = [r for r in records if r.id not in ref_ids]
        out[contig] = ContigAlignment(
            contig=contig,
            ref_row=_to_array(ref_row),
            sample_ids=[r.id for r in samples],
            rows=(
                np.vstack([_to_array(r.seq) for r in samples])
                if samples
                else np.zeros((0, len(ref_row)), dtype=np.uint8)
            ),
            ref_seq=reference.seq(contig),
        )
    return out
