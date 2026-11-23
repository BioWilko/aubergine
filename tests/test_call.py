"""Calling fixtures F1-F9 (one per confirmed legacy bug) and the round-trip invariant.

Reference used by most fixtures (20 bp):
    A1 T2 G3 C4 A5 T6 C7 G8 T9 A10 G11 C12 T13 A14 G15 G16 C17 T18 A19 C20
"""

import pytest
from conftest import keys_of

from aubergine.haplotype import REF, apply, build_sample_events
from aubergine.io import write_calls_tsv

F1 = {  # deletion immediately followed by an insertion (B1)
    "ref": "ATGCATCG--TAGCTAGGCTAC",
    "s1": "ATGCAT--AATAGCTAGGCTAC",
    "s2": "ATGCAT--A-TAGCTAGGCTAC",
}
F2 = {  # insertion immediately followed by a deletion (B1/B2)
    "ref": "ATGCAT--CGTAGCTAGGCTAC",
    "s1": "ATGCATGG--TAGCTAGGCTAC",
}
F3 = {  # insertion whose closing column is an SNP (B2)
    "ref": "ATGCAT--CGTAGCTAGGCTAC",
    "s1": "ATGCATGGAGTAGCTAGGCTAC",
}
F4 = {  # SNP directly after a deletion and directly after an insertion (control)
    "ref": "ATGCATCGTAGC--TAGGCTAC",
    "s1": "ATGC--AGTAGCAACAGGCTAC",
}
F5 = {  # deletion spanning a both-gap column opened by another sample's insertion (B4)
    "ref": "ATGCATC-GTAGCTAGGCTAC",
    "s1": "ATGCATCAGTAGCTAGGCTAC",
    "s2": "ATGCA----TAGCTAGGCTAC",
}
F6 = {  # terminal gaps (B5), last run flushed at the end (B3), reference overhang
    "ref": "--ATGCATCGTAGCTAGGCTAC--",
    "s1": "----GCATCGTAGCTAGGCT----",
    "s2": "GGATGCATCGTAGCTAGGCTACTT",
    "s3": "--ATGCATCGTAGCTAGG--AC--",
}
F7 = {  # insertion starting with N; all-N insertion; IUPAC in insertion (B6)
    "ref": "ATGCAT---CGTAGCTAGGCTAC",
    "s1": "ATGCATNAGCGTAGCTAGGCTAC",
    "s2": "ATGCATNNNCGTAGCTAGGCTAC",
    "s3": "ATGCATRA-CGTAGCTAGGCTAC",
}
F8 = {  # IUPAC at a reference base; N run; reference-N column (B6/B7)
    "ref": "ATGCATCGTANCTAGGCTAC",
    "s1": "ATGCRTCGTAGCTANNNTAC",
}
F9 = {  # identical to the reference (B10)
    "ref": "ATGCATCGTAGCTAGGCTAC",
    "s1": "ATGCATCGTAGCTAGGCTAC",
    "s2": "ATGCATCGTAGCTAGGCTCC",
}
ALL = {"F1": F1, "F2": F2, "F3": F3, "F4": F4, "F5": F5, "F6": F6, "F7": F7, "F8": F8, "F9": F9}


def test_f1_deletion_then_insertion(mk):
    cs = mk(F1)
    assert keys_of(cs, "s1") == ["del:6:TCG:T", "ins:8:G:GAA"]
    assert keys_of(cs, "s2") == ["del:6:TCG:T", "ins:8:G:GA"]
    assert all(v.type != "snp" for v in cs.variants.values())


def test_f2_insertion_then_deletion(mk):
    cs = mk(F2)
    assert keys_of(cs, "s1") == ["del:6:TCG:T", "ins:6:T:TGG"]


def test_f3_insertion_closed_by_snp(mk):
    cs = mk(F3)
    assert keys_of(cs, "s1") == ["ins:6:T:TGG", "snp:7:C:A"]


def test_f4_snp_after_indels(mk):
    cs = mk(F4)
    assert keys_of(cs, "s1") == ["del:4:CAT:C", "snp:7:C:A", "ins:12:C:CAA", "snp:13:T:C"]


def test_f5_deletion_across_both_gap_column(mk):
    cs = mk(F5)
    assert keys_of(cs, "s1") == ["ins:7:C:CA"]
    assert keys_of(cs, "s2") == ["del:5:ATCG:A"]


def test_f6_terminal_gaps_overhang_and_final_run(mk):
    cs = mk(F6)
    s1, s2, s3 = cs.samples
    assert keys_of(cs, "s1") == []
    assert s1.missing == [("ref", 1, 2), ("ref", 19, 20)]
    assert s1.covered == [("ref", 3, 18)]
    assert keys_of(cs, "s2") == []
    assert s2.overhang_bases == 4
    assert s2.summary["overhang_bases"] == 4
    assert keys_of(cs, "s3") == ["del:16:GCT:G"]
    assert s3.missing == []


def test_f7_insertions_with_n_and_iupac(mk):
    cs = mk(F7)
    calls = {c.sample: c for c in cs.calls}
    assert keys_of(cs, "s1") == ["ins:6:T:TNAG"]
    assert calls["s1"].ins_seq_raw == "NAG"
    assert "contains_ambiguity" in calls["s1"].flags
    assert keys_of(cs, "s2") == []
    assert cs.samples[1].missing_insertions == [("ref", 6, 3)]
    assert keys_of(cs, "s3") == ["ins:6:T:TNA"]
    assert calls["s3"].ins_seq_raw == "RA"
    assert "contains_ambiguity" in calls["s3"].flags


def test_f8_ambiguity_missing_and_reference_n(mk):
    cs = mk(F8)
    assert keys_of(cs, "s1") == []
    assert [(a.pos, a.ref, a.code) for a in cs.ambiguous] == [(5, "A", "R")]
    assert cs.samples[0].missing == [("ref", 15, 17)]
    assert cs.samples[0].summary["n_ambiguous"] == 1
    assert cs.samples[0].summary["missing_bases"] == 3


def test_calls_tsv_keeps_ambiguity_rows(mk, tmp_path):
    cs = mk(F8)
    write_calls_tsv(cs, tmp_path / "calls.tsv")
    lines = (tmp_path / "calls.tsv").read_text().splitlines()
    assert lines[1:] == ["s1\tambiguous:ref:5:A:R\t.\t.\t."]


def test_f9_reference_identical_sample_is_reported(mk):
    cs = mk(F9)
    assert cs.sample_ids == ["s1", "s2"]
    assert keys_of(cs, "s1") == []
    assert cs.samples[0].summary["n_snp"] == 0
    assert cs.samples[0].summary["covered_fraction"] == 1.0


def test_seq_start(mk):
    cs = mk(F4)
    starts = {c.key.split(":", 2)[0] + ":" + c.key.split(":")[2]: c.seq_start for c in cs.calls}
    # del of A5,T6 -> next sample base is the 5th; SNP at ref 7 is sample base 5
    assert starts == {"del:4": 5, "snp:7": 5, "ins:12": 11, "snp:13": 13}


def test_near_missing_flag(mk):
    cs = mk({"ref": "ATGCATCGTAGCTAGGCTAC", "s1": "ATGCATCGTNACTAGGCTAC"})
    (call,) = cs.calls
    assert call.key == "snp:ref:11:G:A"
    assert "near_missing" in call.flags


@pytest.mark.parametrize("name", list(ALL))
def test_round_trip(mk, name):
    """apply() over each sample's covered span reproduces its ungapped sequence."""
    rows = ALL[name]
    cs = mk(rows)
    events = build_sample_events(cs)
    ref_seq = rows["ref"].replace("-", "")
    ref_row = rows["ref"]
    for s in cs.samples:
        row = rows[s.id]
        cols = [
            i for i, (r, c) in enumerate(zip(ref_row, row, strict=True)) if r != "-" and c != "-"
        ]
        expected = row[cols[0] : cols[-1] + 1].replace("-", "")
        _, lo, hi = s.covered[0]
        got = apply(ref_seq, events[(s.id, "ref")], lo, hi, missing_as_n=True)
        assert len(got) == len(expected), (s.id, "".join(b.base for b in got), expected)
        for b, e in zip(got, expected, strict=True):
            uncallable_ref = b.kind == REF and b.base not in "ACGT"
            assert b.base == e or uncallable_ref, (s.id, "".join(x.base for x in got), expected)
