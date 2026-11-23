from aubergine.haplotype import INS, MISS, REF, SUB, apply, build_sample_events, window_haplotype


def _events(cs, sample="s1"):
    return build_sample_events(cs)[(sample, "ref")]


def test_insertion_at_position_zero(mk):
    cs = mk({"ref": "A-CGTACG", "s1": "AACGTACG"})
    bases = apply("ACGTACG", _events(cs), 1, 3)
    assert [(b.base, b.kind, b.pos) for b in bases] == [
        ("A", INS, 0),
        ("A", REF, 1),
        ("C", REF, 2),
        ("G", REF, 3),
    ]


def test_insertion_directly_before_deletion_and_origins(mk):
    # ins GG after T6, then C7 G8 deleted, then T9 -> C (SNP)
    cs = mk({"ref": "ATGCAT--CGTAGC", "s1": "ATGCATGG--CAGC"})
    bases = apply("ATGCATCGTAGC", _events(cs), 5, 10)
    assert "".join(b.base for b in bases) == "ATGGCA"
    assert [b.kind for b in bases] == [REF, REF, INS, INS, SUB, REF]
    assert [b.pos for b in bases] == [5, 6, 6, 6, 9, 10]


def test_ambiguous_code_is_a_sample_base(mk):
    cs = mk({"ref": "ATGCATCG", "s1": "ATGCRTCG"})
    bases = apply("ATGCATCG", _events(cs), 4, 6)
    assert [(b.base, b.kind) for b in bases] == [("C", REF), ("R", SUB), ("T", REF)]
    assert window_haplotype("ATGCATCG", _events(cs), 4, 6) == ("ambiguous", "CRT")


def test_missing_classification(mk):
    cs = mk({"ref": "ATGCATCG", "s1": "ATGNNTCG"})
    se = _events(cs)
    assert window_haplotype("ATGCATCG", se, 2, 4)[0] == "missing"
    assert window_haplotype("ATGCATCG", se, 6, 8) == ("callable", "TCG")
    assert [b.kind for b in apply("ATGCATCG", se, 3, 5, missing_as_n=True)] == [REF, MISS, MISS]


def test_missing_insertion_inside_window(mk):
    cs = mk({"ref": "ATGC--ATCG", "s1": "ATGCNNATCG"})
    se = _events(cs)
    assert se.missing_insertions == {4: 2}
    assert window_haplotype("ATGCATCG", se, 3, 6)[0] == "missing"
    assert window_haplotype("ATGCATCG", se, 6, 8)[0] == "callable"
    assert "".join(b.base for b in apply("ATGCATCG", se, 1, 8, missing_as_n=True)) == "ATGCNNATCG"


def test_reference_n_stays_callable(mk):
    cs = mk({"ref": "ATGNATCG", "s1": "ATGAATCG"})
    assert window_haplotype("ATGNATCG", _events(cs), 3, 5) == ("callable", "GNA")
