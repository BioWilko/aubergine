import random

import pytest
from conftest import keys_of, require_tool

from aubergine import mafft
from aubergine.call import call_set
from aubergine.fasta import Reference, load_alignment
from aubergine.params import default_thresholds


def _reference(n=300, seed=7):
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(n))


def _clean_site(ref, q):
    """Move q until ref[q-1], ref[q], ref[q+1] are all different (so no normalisation shift)."""
    while len({ref[q - 2], ref[q - 1], ref[q]}) < 3:
        q += 1
    return q


def test_align_then_call_recovers_known_edits(tmp_path):
    require_tool("mafft")
    ref = _reference()
    snp = _clean_site(ref, 100)
    dele = _clean_site(ref, 150)
    ins = _clean_site(ref, 200)
    snp2 = _clean_site(ref, 250)
    other = {"A": "C", "C": "G", "G": "T", "T": "A"}
    ins_seq = next(b for b in "ACGT" if b not in (ref[ins - 1], ref[ins]))

    def mutate(edits):
        out = list(ref)
        for kind, q in sorted(edits, key=lambda e: -e[1]):
            if kind == "snp":
                out[q - 1] = other[ref[q - 1]]
            elif kind == "del":
                del out[q - 1]
            else:
                out.insert(q, ins_seq * 2)
        return "".join(out)

    (tmp_path / "ref.fasta").write_text(f">ref\n{ref}\n")
    (tmp_path / "set.fasta").write_text(
        f">q1\n{mutate([('snp', snp)])}\n"
        f">q2\n{mutate([('del', dele)])}\n"
        f">q3\n{mutate([('ins', ins), ('snp', snp2)])}\n"
    )
    reference = Reference.from_fasta(tmp_path / "ref.fasta")
    mafft.align(reference, tmp_path / "set.fasta", tmp_path / "set.aligned.fasta")
    cs = call_set(
        load_alignment(tmp_path / "set.aligned.fasta", reference),
        reference,
        "set",
        default_thresholds(),
    )
    assert keys_of(cs, "q1") == [f"snp:{snp}:{ref[snp - 1]}:{other[ref[snp - 1]]}"]
    assert keys_of(cs, "q2") == [f"del:{dele - 1}:{ref[dele - 2 : dele]}:{ref[dele - 2]}"]
    assert keys_of(cs, "q3") == [
        f"ins:{ins}:{ref[ins - 1]}:{ref[ins - 1]}{ins_seq * 2}",
        f"snp:{snp2}:{ref[snp2 - 1]}:{other[ref[snp2 - 1]]}",
    ]


def test_failing_mafft_gives_clear_error(tmp_path, monkeypatch):
    require_tool("mafft")
    monkeypatch.setattr(mafft, "MAFFT_ARGS", ["--no-such-option", "--addfragments"])
    (tmp_path / "ref.fasta").write_text(">ref\nACGTACGTAC\n")
    (tmp_path / "set.fasta").write_text(">a\nACGTACGTAC\n")
    reference = Reference.from_fasta(tmp_path / "ref.fasta")
    with pytest.raises(mafft.MafftError, match="mafft failed"):
        mafft.align(reference, tmp_path / "set.fasta", tmp_path / "out.fasta")
    assert not (tmp_path / "out.fasta").exists()


def test_missing_mafft_suggests_bioconda(monkeypatch):
    monkeypatch.setattr(mafft.shutil, "which", lambda _: None)
    with pytest.raises(mafft.MafftError, match="bioconda"):
        mafft.find_mafft()
