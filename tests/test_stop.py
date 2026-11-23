import pytest
from cds_helpers import R1, R2, edit, gff, revcomp


def _stop(cs):
    st = cs.cds_status[0]
    return st, st["stop"]


def test_snp_creates_stop(mk):
    st, stop = _stop(mk(edit(R2, subs={6: "T"}), gff=gff(3, 29)))
    assert stop["stop_status"] == "premature_stop"
    assert stop["stop_frameshift"] is False
    assert (stop["stop_aa_pos"], stop["protein_len"], stop["ref_protein_len"]) == (2, 1, 8)
    assert stop["stop_ref_pos"] == "ref:6"


def test_one_base_deletion_creates_stop(mk):
    st, stop = _stop(mk(edit(R1, dels={8}), gff=gff(3, 20)))
    assert st["status"] == "frameshift"
    assert stop["stop_status"] == "premature_stop"
    assert stop["stop_frameshift"] is True
    assert (stop["stop_aa_pos"], stop["stop_ref_pos"]) == (3, "ref:10")


def test_frameshift_without_downstream_stop_is_stop_lost(mk):
    st, stop = _stop(mk(edit(R2, dels={10}), gff=gff(3, 29)))
    assert st["status"] == "frameshift"
    assert stop["stop_status"] == "stop_lost"
    assert stop["stop_frameshift"] is True


def test_restored_frameshift_with_stop_inside(mk):
    st, stop = _stop(mk(edit(R1, dels={8}, ins={16: "T"}), gff=gff(3, 20)))
    assert st["status"] == "restored"
    assert stop["stop_status"] == "premature_stop"
    assert stop["stop_frameshift"] is True


def test_deleting_reference_stop_codon(mk):
    st, stop = _stop(mk(edit(R1, dels={18, 19, 20}), gff=gff(3, 20)))
    assert stop["stop_status"] == "stop_lost"
    assert stop["stop_frameshift"] is False


@pytest.mark.parametrize(
    "rows",
    [edit(R1, subs={20: "G"}), edit(R1, dels={12, 13, 14})],
    ids=["TAA_to_TAG", "inframe_deletion"],
)
def test_normal(mk, rows):
    st, stop = _stop(mk(rows, gff=gff(3, 20)))
    assert stop["stop_status"] == "normal"
    assert stop["incomplete"] is False


@pytest.mark.parametrize(
    "subs", [{9: "T", 10: "N"}, {9: "T", 10: "A", 11: "R"}], ids=["TNA", "TAR"]
)
def test_unknown_codon_is_not_a_stop(mk, subs):
    st, stop = _stop(mk(edit(R1, subs=subs), gff=gff(3, 20)))
    assert stop["stop_status"] == "normal"
    assert stop["incomplete"] is True


def test_missing_stop_codon_is_unknown(mk):
    st, stop = _stop(mk(edit(R1, subs={18: "N", 19: "N", 20: "N"}), gff=gff(3, 20)))
    assert stop["stop_status"] == "unknown"
    assert stop["incomplete"] is True


def test_minus_strand_snp_creates_stop(mk):
    # CDS base 4 (A of codon 2) is reference position 26 on the minus strand
    st, stop = _stop(mk(edit(revcomp(R2), subs={26: "A"}), gff=gff(3, 29, strand="-")))
    assert stop["stop_status"] == "premature_stop"
    assert (stop["stop_aa_pos"], stop["stop_ref_pos"]) == (2, "ref:26")


@pytest.mark.parametrize(
    "ref, annotation, note",
    [
        ("CCATGTAAAAATAACC", gff(3, 14), "ref_internal_stop"),
        (R1, gff(3, 19), "ref_length_not_mod3"),
        (R1, gff(3, 20, attrs="ID=c1;gene=g1;transl_table=11"), "unsupported_transl_table"),
    ],
)
def test_not_evaluated(mk, ref, annotation, note):
    st, stop = _stop(mk(edit(ref, subs={1: "G"}), gff=annotation))
    assert stop["stop_status"] == "not_evaluated"
    assert stop["note"] == note


def test_summary_counts(mk):
    cs = mk(edit(R2, subs={6: "T"}), gff=gff(3, 29))
    assert cs.samples[0].summary["premature_stop_cds"] == 1
    assert cs.samples[0].summary["stop_lost_cds"] == 0
