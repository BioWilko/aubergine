import gzip

import pytest

from aubergine.fasta import FastaError, Reference, load_alignment, read_fasta


def test_read_normalises_case_u_and_dots(tmp_path):
    p = tmp_path / "x.fa"
    p.write_text(">a desc\nacgu\n..n\n>b\nACGT-N\n")
    recs = read_fasta(p)
    assert [(r.id, r.seq) for r in recs] == [("a", "ACGT--N"), ("b", "ACGT-N")]
    assert recs[0].description == "a desc"


def test_read_gz(tmp_path):
    p = tmp_path / "x.fa.gz"
    with gzip.open(p, "wt") as fh:
        fh.write(">a\nACGT\n")
    assert read_fasta(p)[0].seq == "ACGT"


def test_invalid_character(tmp_path):
    p = tmp_path / "x.fa"
    p.write_text(">a\nACGTX\n")
    with pytest.raises(FastaError, match="invalid characters: X"):
        read_fasta(p)


def test_duplicate_ids(tmp_path):
    p = tmp_path / "x.fa"
    p.write_text(">a\nACGT\n>a other\nACGT\n")
    with pytest.raises(FastaError, match="duplicate record IDs: a"):
        read_fasta(p)


def test_reference_rejects_gaps(tmp_path):
    p = tmp_path / "r.fa"
    p.write_text(">r\nAC-GT\n")
    with pytest.raises(FastaError):
        Reference.from_fasta(p)


def _ref(tmp_path, seq="ACGTACGT"):
    p = tmp_path / "r.fa"
    p.write_text(f">ref\n{seq}\n")
    return Reference.from_fasta(p)


def test_ragged_alignment(tmp_path):
    ref = _ref(tmp_path)
    p = tmp_path / "a.fa"
    p.write_text(">ref\nACGTACGT\n>s\nACGTACG\n")
    with pytest.raises(FastaError, match="different lengths"):
        load_alignment(p, ref)


def test_missing_reference_id(tmp_path):
    ref = _ref(tmp_path)
    p = tmp_path / "a.fa"
    p.write_text(">other\nACGTACGT\n>s\nACGTACGT\n")
    with pytest.raises(FastaError, match="reference record 'ref' not found"):
        load_alignment(p, ref)


def test_reference_row_must_match(tmp_path):
    ref = _ref(tmp_path)
    p = tmp_path / "a.fa"
    p.write_text(">ref\nACGTACGA\n>s\nACGTACGT\n")
    with pytest.raises(FastaError, match="does not match"):
        load_alignment(p, ref)


def test_reference_not_first_and_lowercase(tmp_path):
    ref = _ref(tmp_path)
    p = tmp_path / "a.fa"
    p.write_text(">s\nacgtacgt\n>ref\nacgt-acgt\n".replace("acgtacgt\n>ref", "acgt-acgt\n>ref"))
    aln = load_alignment(p, ref)["ref"]
    assert aln.sample_ids == ["s"]
    assert list(aln.ref_pos) == [1, 2, 3, 4, 4, 5, 6, 7, 8]
