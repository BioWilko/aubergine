import pytest
from cds_helpers import R2, edit, gff, revcomp

from aubergine.cds import AnnotationError, parse_gff3
from aubergine.fasta import Reference


def _status(cs, sample="s1"):
    return next(s for s in cs.cds_status if s["sample"] == sample)


def _fs(cs):
    return [c.frameshift for c in cs.calls]


def test_one_base_deletion_is_frameshift(mk):
    cs = mk(edit(R2, dels={10}), gff=gff(3, 29))
    assert _fs(cs) == ["g1:fs"]
    st = _status(cs)
    assert (st["net_delta"], st["status"]) == (-1, "frameshift")
    assert st["out_of_frame"] == [[9, 29]]
    assert cs.variants[cs.calls[0].key].cds == [{"gene": "g1", "delta_cds": -1}]


def test_three_base_deletion_is_inframe(mk):
    cs = mk(edit(R2, dels={12, 13, 14}), gff=gff(3, 29))
    assert _fs(cs) == ["g1:inframe"]
    assert _status(cs)["status"] == "inframe"


def test_compensating_indels_are_restored(mk):
    cs = mk(edit(R2, dels={6}, ins={20: "T"}), gff=gff(3, 29))
    assert _fs(cs) == ["g1:fs_compensated", "g1:fs_compensated"]
    st = _status(cs)
    assert (st["net_delta"], st["status"], st["out_of_frame"]) == (0, "restored", [[6, 20]])
    assert cs.samples[0].summary["restored_cds"] == 1


def test_minus_strand(mk):
    cs = mk(edit(revcomp(R2), dels={10}), gff=gff(3, 29, strand="-"))
    assert _fs(cs) == ["g1:fs"]
    st = _status(cs)
    assert st["status"] == "frameshift"
    assert st["out_of_frame"][0][0] == 3  # runs to the CDS end, which is at the left


def test_two_segments_and_slippage_base(mk):
    segs = [(3, 15), (15, 29)]
    junction = mk(edit(R2, dels={14, 15, 16}), gff=gff(0, 0, segments=segs))
    assert junction.variants[junction.calls[0].key].cds == [{"gene": "g1", "delta_cds": -4}]
    slip = mk(edit(R2, dels={15}), gff=gff(0, 0, segments=segs))
    assert slip.variants[slip.calls[0].key].cds == [{"gene": "g1", "delta_cds": -2}]


def test_insertion_before_start_codon_is_not_coding(mk):
    cs = mk(edit(R2, ins={2: "G"}), gff=gff(3, 29))
    assert _fs(cs) == ["."]
    assert _status(cs)["status"] == "intact"


def test_deletion_over_stop_codon_flags_boundary(mk):
    cs = mk(edit(R2, dels={28}), gff=gff(3, 29))
    assert "cds_boundary" in cs.calls[0].flags


def test_missing_data_inside_cds_is_incomplete(mk):
    cs = mk(edit(R2, subs={10: "N"}), gff=gff(3, 29))
    assert _status(cs)["incomplete"] is True


def test_no_annotation(mk):
    cs = mk(edit(R2, dels={10}))
    assert _fs(cs) == ["."]
    assert cs.cds_status == []
    assert cs.samples[0].summary["frameshifted_cds"] == "."


def test_unknown_seqid_is_an_error(tmp_path):
    (tmp_path / "r.fa").write_text(f">ref\n{R2}\n")
    (tmp_path / "a.gff3").write_text(gff(3, 29, seqid="chrX"))
    with pytest.raises(AnnotationError, match="chrX"):
        parse_gff3(tmp_path / "a.gff3", Reference.from_fasta(tmp_path / "r.fa"))


def test_duplicate_gene_labels_are_disambiguated(tmp_path):
    (tmp_path / "r.fa").write_text(f">ref\n{R2}\n")
    (tmp_path / "a.gff3").write_text(
        gff(3, 29, attrs="ID=a;gene=orf1") + gff(3, 14, attrs="ID=b;gene=orf1")
    )
    cds = parse_gff3(tmp_path / "a.gff3", Reference.from_fasta(tmp_path / "r.fa"))
    assert sorted(c.gene for c in cds) == ["orf1[a]", "orf1[b]"]
