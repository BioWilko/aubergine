from conftest import FIXTURES, analysis

from aubergine.cli import main
from aubergine.compare import pairwise_rows

# A1 T2 G3 C4 A5 T6 C7 G8 T9 A10 [ins ins] G11 C12 T13 A14 G15 G16 C17 T18 A19 C20
STATES = {
    "ref": "ATGCATCGTA--GCTAGGCTAC",
    "s1": "ATGCATAGTATTGC-AGGCTAC",  # snp 7 C>A, ins TT after 10, del T13
    "s2": "ATGCATNGTAG-G---GGCTAC",  # N at 7, ins G after 10, del 12-14
    "s3": "ATGCATRGTA--GCNAGGCTAC",  # R at 7, N at 13
    "s4": "ATGCA---TA--GCTAGGCTAC",  # del 6-8
    "s5": "ATGCATGGTA--GCGAGGCTAC",  # snp 7 C>G, snp 13 T>G
    "s6": "ATGCATCGTA--GCKAGGCTAC",  # K at 13
    "s7": "ATGCATCGTANNGCTAGGCTAC",  # all-N insertion after 10
}


def test_every_state(mk):
    cs = mk(STATES)
    states, _, _ = analysis(cs)
    samples = ["s1", "s2", "s3", "s4", "s5", "s6", "s7"]
    row = {k.replace(":ref:", ":"): [states[k][s] for s in samples] for k in states}
    assert row["snp:7:C:A"] == [
        "alt",
        "missing",
        "ambiguous",
        "deleted",
        "other_allele",
        "ref",
        "ref",
    ]
    assert row["del:12:CT:C"] == [
        "alt",
        "deleted",
        "missing",
        "ref",
        "other_allele",
        "ambiguous",
        "ref",
    ]
    assert row["ins:10:A:ATT"] == [
        "alt",
        "other_allele",
        "ref",
        "ref",
        "ref",
        "ref",
        "missing",
    ]


def test_pairwise_counts_and_na_jaccard(mk):
    cs = mk(STATES)
    states, _, _ = analysis(cs)
    rows = {(r[0], r[1], r[2]): r[3:] for r in pairwise_rows(cs, states)}
    # s1 vs s5 SNPs: s1 has 7A; s5 has 7G and 13G -> 0 shared, 1 only_a, 2 only_b
    assert rows[("s1", "s5", "snp")] == [0, 1, 2, 3, 0, "0.0000"]
    # s6 and s7 carry no SNPs: nothing shared or private -> NA
    snp_67 = rows[("s6", "s7", "snp")]
    assert snp_67[:3] == [0, 0, 0] and snp_67[-1] == "NA"
    # s2 is missing at 7 -> that key is excluded for any pair involving s2
    assert rows[("s1", "s2", "snp")][4] >= 1
    all_row = rows[("s1", "s5", "all")]
    parts = [rows[("s1", "s5", t)] for t in ("snp", "insertion", "deletion")]
    assert all_row[:5] == [sum(p[i] for p in parts) for i in range(5)]


def test_wide_variants_columns(tmp_path):
    out = tmp_path / "ex"
    src = FIXTURES / "discord_example"
    assert (
        main(
            [
                "-q",
                "call",
                "-r",
                str(src / "ref.fasta"),
                "-o",
                str(out),
                str(src / "example.aligned.fasta"),
            ]
        )
        == 0
    )
    assert main(["-q", "compare", "-o", str(out), f"{out}.variants.json"]) == 0
    lines = [ln.split("\t") for ln in (tmp_path / "ex.variants.tsv").read_text().splitlines()]
    header = lines[0]
    cols = {r[0]: dict(zip(header, r, strict=True)) for r in lines[1:]}
    assert cols["del:ref:5:CA:C"]["locus"] == "L1:6-11"
    assert "locus_reason" not in header
    assert cols["snp:ref:22:T:C"]["locus"] == "L1:22-22"
    assert cols["snp:ref:22:T:C"]["s3"] == "missing"
