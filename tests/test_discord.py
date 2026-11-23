import random

import pytest
from cds_helpers import gff
from conftest import FIXTURES, analysis

from aubergine.cli import main
from aubergine.discord import (
    REASONS,
    SET_SUMMARY_HEADER,
    Locus,
    LocusResult,
    build_loci,
    cds_effects,
    set_summary_row,
)
from aubergine.io import read_callset

EXAMPLE = FIXTURES / "discord_example"


def _run_example(tmp_path, extra=()):
    out = tmp_path / "example"
    call = ["-q", "call", "-r", str(EXAMPLE / "ref.fasta"), "-o", str(out)]
    assert main([*call, str(EXAMPLE / "example.aligned.fasta")]) == 0
    assert main(["-q", "compare", "-o", str(out), *extra, f"{out}.variants.json"]) == 0
    return out


def _tsv(path):
    lines = path.read_text().splitlines()
    header = lines[0].split("\t")
    return [dict(zip(header, ln.split("\t"), strict=True)) for ln in lines[1:]]


def test_worked_example_end_to_end(tmp_path):
    out = _run_example(tmp_path)
    for name in ("discordance", "discordance_pairs", "discordance_summary"):
        got = (tmp_path / f"example.{name}.tsv").read_text()
        assert got == (EXAMPLE / "expected" / f"{name}.tsv").read_text(), name
    report = (tmp_path / "example.discordance_report.md").read_text()
    assert report == (EXAMPLE / "expected" / "discordance_report.md").read_text()
    pw = {
        (r["sample_a"], r["sample_b"], r["type"]): r
        for r in _tsv(tmp_path / "example.pairwise.tsv")
    }
    row = pw[("s1", "s2", "all")]
    assert (row["shared"], row["only_a"], row["only_b"], row["jaccard"]) == (
        "1",
        "3",
        "3",
        "0.1429",
    )
    assert out.with_suffix(".variants.json").exists()


def _reasons(cs, **kw):
    _, loci, _ = analysis(cs, **kw)
    return {(p.a, p.b): p.reason for r in loci for p in r.pairs}, loci


# Minimal fixtures, one per reason code. Reference: A1 T2 G3 C4 A5 T6 C7 G8 T9 A10 G11 C12
REF12 = "ATGCATCGTAGC"
REASON_FIXTURES = {
    "missing_data": {"ref": REF12, "s1": "ATGCAACGTAGC", "s2": "ATGCANCGTAGC"},
    "ambiguity": {"ref": REF12, "s1": "ATGCAACGTAGC", "s2": "ATGCAWCGTAGC"},
    # del of C7 + ins C after G8 vs two SNPs C7G, G8C
    "equivalent_representation": {
        "ref": "ATGCATCG-TAGC",
        "s1": "ATGCAT-GCTAGC",
        "s2": "ATGCATGC-TAGC",
    },
    "paired_indel": {"ref": "ATGCATCG-TAGC", "s1": "ATGCAT-GCTAGC", "s2": "ATGCATCG-TAGC"},
    "indel_length": {"ref": "ATGCATATATCGC", "s1": "ATGC--ATATCGC", "s2": "ATGC----ATCGC"},
    "indel_sequence": {"ref": "ATGCAT--CGTAGC", "s1": "ATGCATGACGTAGC", "s2": "ATGCATCACGTAGC"},
    "indel_presence": {"ref": REF12, "s1": "ATGCA-CGTAGC", "s2": REF12},
    "different_allele": {"ref": REF12, "s1": "ATGCAACGTAGC", "s2": "ATGCAGCGTAGC"},
    "snp_presence": {"ref": REF12, "s1": "ATGCAACGTAGC", "s2": REF12},
    "mixed": {"ref": REF12, "s1": "ATGCAACGTAGC", "s2": "ATGCAT-GTAGC"},
}


@pytest.mark.parametrize("reason", REASONS)
def test_each_reason_code(mk, reason):
    reasons, _ = _reasons(mk(REASON_FIXTURES[reason]))
    assert reasons == {("s1", "s2"): reason}


def test_indel_length_homopolymer_flag(mk):
    _, outside = _reasons(mk(REASON_FIXTURES["indel_length"]))
    assert "homopolymer" not in outside[0].flags
    _, inside = _reasons(mk({"ref": "ATGCAAAAAACGC", "s1": "ATGC-AAAAACGC", "s2": "ATGC--AAAACGC"}))
    assert inside[0].pairs[0].reason == "indel_length"
    assert "homopolymer" in inside[0].flags


def test_precedence(mk):
    rows = {
        "ref": "ATGCATCG-TAGC",
        "s1": "ATGCAT-GCTAGC",  # del + ins
        "s2": "ATGCATGC-TAGC",  # same haplotype as s1, as two SNPs
        "s3": "ATGCATGC-NAGC",  # s2's SNPs plus N inside the locus window (7-9)
        "s4": "ATGCATGC-WAGC",  # s2's SNPs plus an ambiguity inside the window
    }
    reasons, loci = _reasons(mk(rows))
    assert len(loci) == 1
    assert reasons[("s1", "s2")] == "equivalent_representation"  # not paired_indel
    assert reasons[("s1", "s3")] == "missing_data"  # missing beats equivalent
    assert reasons[("s1", "s4")] == "ambiguity"  # ambiguity beats equivalent


def test_equivalent_del_snp_vs_snp_del(mk):
    # C1 C2 A3 T4 G5 C6 C7: del T4 + G5C  vs  T4C + del G5
    reasons, loci = _reasons(mk({"ref": "CCATGCC", "s1": "CCA-CCC", "s2": "CCAC-CC"}))
    assert reasons == {("s1", "s2"): "equivalent_representation"}
    assert loci[0].haplotype_concordant and not loci[0].callable_discordant


def test_equivalent_when_normalisation_blocked_in_one_sample(mk):
    reasons, _ = _reasons(mk({"ref": "ACGTTTACG", "s1": "ACGTC-ACG", "s2": "ACGT-CACG"}))
    assert reasons == {("s1", "s2"): "equivalent_representation"}


def _loci(cs, gap):
    return [(loc.start, loc.end) for loc in build_loci(cs, gap)]


def test_locus_join_gap(mk):
    # SNPs at 5 and 8 (2 bases between) join; 5 and 9 (3 bases between) do not
    two = mk({"ref": REF12, "s1": "ATGCCTCTTAGC"})
    assert _loci(two, 2) == [(5, 8)]
    cs = mk({"ref": REF12, "s1": "ATGCCTCGAAGC"})
    assert _loci(cs, 2) == [(5, 5), (9, 9)]
    assert _loci(cs, 3) == [(5, 9)]


def test_large_deletion_absorbs_nearby_snps(mk):
    cs = mk({"ref": "ATGCATCGTAGCTAGG", "s1": "AT----------TCGG", "s2": "ATGCATCGTAGCTCGG"})
    assert _loci(cs, 2) == [(3, 14)]


def test_flags_follow_thresholds(mk):
    cs = mk(
        {
            "ref": "ACGTCAAAAAAGTCCGATGCATGC",
            "s1": "ACGTC-AAAAAGTCCGATGCATGC",
            "s2": "ACGTC--AAAAGTCCGATGCATGC",
        }
    )
    _, loci = _reasons(cs)
    assert loci[0].flags == ["near_end", "homopolymer"]
    _, loci = _reasons(cs, end_dist=0, hp_min=7)
    assert loci[0].flags == []


def test_frameshift_flag(mk):
    gff = "ref\t.\tCDS\t3\t11\t.\t+\t0\tID=c1;gene=g\n"
    _, loci = _reasons(
        mk({"ref": "CCATGAAATAACC", "s1": "CCATG-AATAACC", "s2": "CCATGAAATAACC"}, gff=gff)
    )
    assert "frameshift" in loci[0].flags


def test_equal_keysets_not_reported_and_single_sample(mk, tmp_path):
    reasons, _ = _reasons(
        mk({"ref": REF12, "s1": "ATGCAACGTAGC", "s2": "ATGCAACGTAGC", "s3": REF12})
    )
    assert ("s1", "s2") not in reasons
    one = tmp_path / "one"
    one.mkdir()
    (one / "r.fa").write_text(f">ref\n{REF12}\n")
    (one / "a.fa").write_text(f">ref\n{REF12}\n>s1\nATGCAACGTAGC\n")
    assert (
        main(["-q", "call", "-r", str(one / "r.fa"), "-o", str(one / "x"), str(one / "a.fa")]) == 0
    )
    assert main(["-q", "compare", "-o", str(one / "x"), str(one / "x.variants.json")]) == 0
    assert len((one / "x.discordance.tsv").read_text().splitlines()) == 1
    assert len((one / "x.discordance_pairs.tsv").read_text().splitlines()) == 1
    summary = _tsv(one / "x.discordance_summary.tsv")
    assert [r["reason"] for r in summary] == REASONS
    assert all(r["n_pairs"] == "0" for r in summary)
    report = (one / "x.discordance_report.md").read_text()
    assert "**0 discordant loci**" in report
    assert "no discordant loci" in report


def _locus_result(labels):
    return LocusResult(Locus("L1:1-1", "ref", 1, 1), labels, {}, [], False, False, [])


@pytest.mark.parametrize(
    "labels, odd",
    [
        ({"a": "H0", "b": "H0", "c": "H1", "d": "missing"}, ["c"]),
        ({"a": "H0", "b": "H1"}, []),  # tie: no majority
        ({"a": "H1", "b": "H1", "c": "missing"}, []),  # one callable haplotype
        ({"a": "H0", "b": "H2", "c": "H2", "d": "H1", "e": "ambiguous"}, ["a", "d"]),
    ],
)
def test_odd_samples(labels, odd):
    assert _locus_result(labels).odd_samples == odd


# CDS 3..26: ATG TGG AAA CCC GGG TTT AAA TAA
REF_CDS = "CC" + "ATGTGGAAACCCGGGTTTAAATAA" + "CC"


def test_cds_effects(mk):
    s1 = REF_CDS[:6] + "A" + REF_CDS[7:]  # TGG -> TAG at codon position 2 (ref 7)
    s2 = REF_CDS[:12] + "-" + REF_CDS[13:]  # one C of CCC deleted: frameshift
    cs = mk({"ref": REF_CDS, "s1": s1, "s2": s2, "s3": REF_CDS}, gff=gff(3, 26))
    stop = next(st["stop"] for st in cs.cds_status if st["sample"] == "s1")
    assert stop["stop_status"] == "premature_stop"
    assert stop["stop_ref_pos"] == "ref:6" and stop["stop_codon"] == [6, 8]
    _, results, _ = analysis(cs)
    effects = cds_effects(cs, results)
    by_start = {r.locus.start: effects[r.locus.id] for r in results}
    assert by_start[7] == ["g1:premature_stop(s1)"]
    assert by_start[12] == ["g1:fs(s2)"]


def test_shuffled_sample_order_changes_only_columns(tmp_path):
    rows = (EXAMPLE / "example.aligned.fasta").read_text().split(">")[1:]
    rows = [">" + r for r in rows]
    shuffled = [rows[0]] + random.Random(3).sample(rows[1:], 3)
    assert shuffled[1:] != rows[1:]
    (tmp_path / "shuf.aligned.fasta").write_text("".join(shuffled))

    def labels(path):
        return {
            r["locus"]: {
                **{k: r[k] for k in ("s1", "s2", "s3", "primary_reason", "reasons")},
                "odd_samples": set(r["odd_samples"].split(";")),
            }
            for r in _tsv(path)
        }

    a = _run_example(tmp_path)
    out = tmp_path / "shuf"
    assert (
        main(
            [
                "-q",
                "call",
                "-r",
                str(EXAMPLE / "ref.fasta"),
                "-o",
                str(out),
                str(tmp_path / "shuf.aligned.fasta"),
            ]
        )
        == 0
    )
    assert main(["-q", "compare", "-o", str(out), f"{out}.variants.json"]) == 0
    assert labels(a.with_suffix(".discordance.tsv")) == labels(out.with_suffix(".discordance.tsv"))


def test_set_summary_row_worked_example(tmp_path):
    out = _run_example(tmp_path)
    cs = read_callset(out.with_suffix(".variants.json"))
    _, results, _ = analysis(cs)
    row = dict(zip(SET_SUMMARY_HEADER, set_summary_row(cs, results), strict=True))
    reasons = {r: row[r] for r in REASONS if row[r]}
    assert reasons == {"paired_indel": 1, "indel_length": 1, "missing_data": 1}
    assert (row["n_samples"], row["n_variant_keys"]) == (3, 7)
    # L1:6-11 (6 bp) + L1:16-18 (3 bp) + L1:22-22 (1 bp)
    assert (row["discordant_loci"], row["discordant_positions"]) == (3, 10)
    assert (row["callable_discordant_loci"], row["uncallable_only_loci"]) == (2, 1)
    assert (row["near_end"], row["homopolymer"], row["frameshift"]) == (3, 1, 0)
    assert (row["fs_loci"], row["premature_stop_loci"]) == (None, None)  # no annotation
    assert row["odd_one_out"] == ["s3:1"]
    assert row["most_discordant_pair"] == "s1,s3:2"
