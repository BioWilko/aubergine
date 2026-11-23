import pytest

from aubergine.call import make_key, parse_key
from aubergine.cli import main


@pytest.mark.parametrize("contig", ["MN908947.3", "chr:1", "a:b:c"])
def test_key_round_trip(contig):
    key = make_key("del", contig, 10, "GT", "G")
    assert key == f"del:{contig}:10:GT:G"
    assert parse_key(key) == ("del", contig, 10, "GT", "G")


def test_two_record_reference_is_rejected(tmp_path, caplog):
    (tmp_path / "ref.fasta").write_text(">seg1\nACGT\n>seg2\nTTGA\n")
    (tmp_path / "aln.fasta").write_text(">seg1\nACGT\n>s\nACGT\n")
    rc = main(
        [
            "call",
            "-r",
            str(tmp_path / "ref.fasta"),
            "-o",
            str(tmp_path / "x"),
            str(tmp_path / "aln.fasta"),
        ]
    )
    assert rc == 1
    assert "multi-record references are not supported yet" in caplog.text
