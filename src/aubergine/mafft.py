"""mafft wrapper: align a set of sequences to a reference with --addfragments."""

from __future__ import annotations

import os
import shutil
import subprocess
import tempfile
from pathlib import Path

from .fasta import FastaError, Record, Reference, read_fasta, write_fasta
from .params import MAFFT_ARGS


class MafftError(RuntimeError):
    pass


def find_mafft() -> str:
    exe = shutil.which("mafft")
    if exe is None:
        raise MafftError(
            "mafft not found on PATH; install it with `conda install -c bioconda mafft`"
        )
    return exe


def mafft_version() -> str:
    proc = subprocess.run([find_mafft(), "--version"], capture_output=True, text=True, check=False)
    return (
        (proc.stderr or proc.stdout).strip().splitlines()[0] if (proc.stderr or proc.stdout) else ""
    )


def _tail(text: str, n: int = 20) -> str:
    return "\n".join(text.strip().splitlines()[-n:])


def align_contig(ref: Record, samples: list[Record], out_path: Path, threads: int = 1) -> list[str]:
    """Align samples to one reference record; write the MSA atomically. Returns the command."""
    exe = find_mafft()
    out_path = Path(out_path)
    out_dir = out_path.parent
    out_dir.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=out_dir, prefix=".aubergine-mafft-") as tmp:
        ref_fa = Path(tmp, "reference.fasta")
        set_fa = Path(tmp, "set.fasta")
        tmp_out = Path(tmp, "aligned.fasta")
        write_fasta(ref_fa, [Record(ref.id, ref.id, ref.seq)])
        write_fasta(set_fa, [Record(r.id, r.id, r.seq.replace("-", "")) for r in samples])
        cmd = [exe, "--thread", str(threads), *MAFFT_ARGS, str(set_fa), str(ref_fa)]
        with open(tmp_out, "w") as fh:
            proc = subprocess.run(cmd, stdout=fh, stderr=subprocess.PIPE, text=True, check=False)
        if proc.returncode != 0:
            raise MafftError(f"mafft failed (exit {proc.returncode}):\n{_tail(proc.stderr)}")
        if tmp_out.stat().st_size == 0:
            raise MafftError(f"mafft produced no output:\n{_tail(proc.stderr)}")
        try:
            records = read_fasta(tmp_out)
        except FastaError as e:
            raise MafftError(f"mafft output is not a valid alignment: {e}") from e
        ids = {r.id for r in records}
        if ref.id not in ids:
            raise MafftError(f"mafft output lacks the reference record {ref.id!r}")
        aligned_ref = next(r for r in records if r.id == ref.id)
        if aligned_ref.seq.replace("-", "") != ref.seq:
            raise MafftError("mafft output: ungapped reference row differs from the reference")
        missing = [r.id for r in samples if r.id not in ids]
        if missing:
            raise MafftError(f"mafft output lacks input records: {', '.join(missing)}")
        os.replace(tmp_out, out_path)
    return cmd


def align(reference: Reference, set_path: Path, out_path: Path, threads: int = 1) -> list[str]:
    """Align one set. v1: the reference has one contig, and every record is a sample of it."""
    samples = read_fasta(set_path)
    if not samples:
        raise MafftError(f"{set_path}: no sequences to align")
    clash = [r.id for r in samples if r.id in reference.contigs]
    if clash:
        raise MafftError(f"{set_path}: sample IDs clash with reference IDs: {', '.join(clash)}")
    cmd: list[str] = []
    for contig in reference.ids:
        cmd = align_contig(
            Record(contig, contig, reference.seq(contig)), samples, out_path, threads
        )
    return cmd
