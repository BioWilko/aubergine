import pytest
from conftest import FIXTURES, analysis

from aubergine.cli import main
from aubergine.compare import CompareError, truth_tables

EXAMPLE = FIXTURES / "truth_example"


def _tsv(path):
    lines = path.read_text().splitlines()
    header = lines[0].split("\t")
    return [dict(zip(header, ln.split("\t"), strict=True)) for ln in lines[1:]]


def test_worked_example(tmp_path):
    out = tmp_path / "ex"
    assert (
        main(
            [
                "-q",
                "call",
                "-r",
                str(EXAMPLE / "ref.fasta"),
                "-o",
                str(out),
                str(EXAMPLE / "example.aligned.fasta"),
            ]
        )
        == 0
    )
    assert main(["-q", "compare", "--truth", "T", "-o", str(out), f"{out}.variants.json"]) == 0
    metrics = {r["type"]: r for r in _tsv(tmp_path / "ex.truth_metrics.tsv")}
    assert {r["sample"] for r in metrics.values()} == {"Q"}

    def nums(t):
        r = metrics[t]
        return [r[k] for k in ("tp", "fp", "fn", "uncalled", "precision", "recall", "f1")]

    assert nums("snp") == ["1", "1", "1", "0", "0.5000", "0.5000", "0.5000"]
    assert nums("insertion") == ["0", "0", "0", "1", "NA", "NA", "NA"]
    assert nums("deletion") == ["1", "0", "0", "0", "1.0000", "1.0000", "1.0000"]
    assert nums("all") == ["2", "1", "1", "1", "0.6667", "0.6667", "0.6667"]
    assert metrics["all"]["equivalent"] == "0"

    detail = {r["key"]: r for r in _tsv(tmp_path / "ex.truth_detail.tsv")}
    assert detail["del:ref:10:GT:G"]["class"] == "TP"
    assert detail["del:ref:10:GT:G"]["reason"] == "."
    assert detail["snp:ref:3:G:A"]["class"] == "FN"
    assert detail["snp:ref:3:G:A"]["reason"] == "snp_presence"
    assert detail["ins:ref:15:C:CTT"]["class"] == "UNCALLED"
    assert detail["ins:ref:15:C:CTT"]["reason"] == "missing_data"
    # The SNP at 18 shares locus L1:11-18 with the insertion Q cannot see (N at 15-16),
    # so the pair reason at that locus is missing_data (locus_join_gap = 2).
    assert detail["snp:ref:18:T:C"]["class"] == "FP"
    assert detail["snp:ref:18:T:C"]["locus"] == "L1:11-18"
    assert detail["snp:ref:18:T:C"]["reason"] == "missing_data"


def test_truth_uncallable(mk):
    cs = mk({"ref": "ATGCATCGTAGC", "T": "ATGCANCGTAGC", "Q": "ATGCAACGTAGC"})
    states, loci, _ = analysis(cs)
    metrics, detail = truth_tables(cs, states, loci, "T")
    assert [d[5] for d in detail] == ["TRUTH_UNCALLABLE"]
    snp = next(m for m in metrics if m[3] == "snp")
    assert (snp[5], snp[8]) == (0, 1)  # fp, truth_uncallable


def test_equivalent_counts_decomposition_difference(mk):
    cs = mk({"ref": "CCATGCC", "T": "CCA-CCC", "Q": "CCAC-CC"})
    states, loci, _ = analysis(cs)
    metrics, detail = truth_tables(cs, states, loci, "T")
    assert {d[5] for d in detail} == {"FP", "FN"}
    assert {d[9] for d in detail} == {"equivalent_representation"}
    total = next(m for m in metrics if m[3] == "all")
    assert total[9] == 4  # equivalent: 2 FN + 2 FP


def test_unknown_truth_and_exclusion(mk):
    cs = mk({"ref": "ATGCATCGTAGC", "T": "ATGCAACGTAGC", "Q": "ATGCAACGTAGC"})
    states, loci, _ = analysis(cs)
    with pytest.raises(CompareError, match="not a sample"):
        truth_tables(cs, states, loci, "nope")
    metrics, _ = truth_tables(cs, states, loci, "T")
    assert {m[2] for m in metrics} == {"Q"}
