from conftest import keys_of


def test_same_deletion_in_homopolymer_gives_same_key(mk):
    cs = mk({"ref": "ACGTTTACG", "s1": "ACG-TTACG", "s2": "ACGTT-ACG"})
    assert keys_of(cs, "s1") == keys_of(cs, "s2") == ["del:3:GT:G"]


def test_same_insertion_in_repeat_gives_same_key(mk):
    cs = mk(
        {
            "ref": "GTCA--CACA--GT",
            "s1": "GTCACACACA--GT",
            "s2": "GTCA--CACACAGT",
        }
    )
    assert keys_of(cs, "s1") == keys_of(cs, "s2") == ["ins:2:T:TCA"]


def test_shift_blocked_by_own_snp(mk):
    cs = mk({"ref": "ACGTTTACG", "s1": "ACGTC-ACG", "s2": "ACGTT-ACG"})
    assert keys_of(cs, "s1") == ["del:5:TT:T", "snp:5:T:C"]
    assert keys_of(cs, "s2") == ["del:3:GT:G"]


def test_shift_blocked_by_n(mk):
    cs = mk({"ref": "ACGTTTACG", "s1": "ACGTN-ACG"})
    assert keys_of(cs, "s1") == ["del:5:TT:T"]


def test_insertion_shift_blocked_by_own_snp(mk):
    cs = mk({"ref": "ACGTTT-ACG", "s1": "ACGTCTTACG"})
    # inserted T slides from after T6 to after T5, then stops at the sample's own SNP at 5
    assert keys_of(cs, "s1") == ["ins:5:T:TT", "snp:5:T:C"]


def test_deletion_at_position_one_is_right_anchored(mk):
    cs = mk({"ref": "AACGTACG", "s1": "A-CGTACG"})
    assert keys_of(cs, "s1") == ["del:1:AA:A"]
    v = next(iter(cs.variants.values()))
    assert (v.start, v.end) == (1, 1)


def test_insertion_at_position_zero(mk):
    cs = mk({"ref": "A-CGTACG", "s1": "AACGTACG"})
    assert keys_of(cs, "s1") == ["ins:1:A:AA"]
    v = next(iter(cs.variants.values()))
    assert (v.start, v.end, v.ins_seq) == (0, 1, "A")


def test_slide_spans(mk):
    cs = mk(
        {
            "ref": "ACGTTTT-ACG",
            "s1": "ACGTTT--ACG",
            "s2": "ACGTTTTTACG",
            "s3": "ACCTTTT-ACG",
        }
    )
    by_type = {v.type: v for v in cs.variants.values()}
    assert by_type["del"].slide == (4, 7)
    assert by_type["ins"].slide == (3, 8)
    assert by_type["snp"].slide == (3, 3)


def test_slide_clipped_at_contig_start(mk):
    cs = mk({"ref": "A-AACG", "s1": "AAAACG"})
    (v,) = cs.variants.values()
    assert v.slide == (1, 4)


def test_hp_len_at_reference_start_and_end(mk):
    cs = mk(
        {
            "ref": "AAAAAACGTCGTTTTTT",
            "s1": "AAAAA-CGTCGTTTTTT",
            "s2": "AAAAAACGTCGT-TTTT",
        }
    )
    hp = {k.split(":")[0] + k.split(":")[2]: v for k, v in cs.variants.items()}
    assert cs.variants["del:ref:1:AA:A"].hp_len == 6
    assert cs.variants["del:ref:1:AA:A"].hp_run == ("ref", 1, 6)
    assert cs.variants["del:ref:11:GT:G"].hp_len == 6
    assert cs.variants["del:ref:11:GT:G"].hp_run == ("ref", 12, 17)
    assert len(hp) == 2


def test_hp_len_next_to_insertion(mk):
    cs = mk(
        {
            "ref": "CGAAAAA-CGT",
            "s1": "CGAAAAAACGT",
            "s2": "CGAAAAA-CGT",
            "s3": "CGAAAAAGCGT",
        }
    )
    assert cs.variants["ins:ref:2:G:GA"].hp_len == 5
    assert cs.variants["ins:ref:7:A:AG"].hp_len == 0
    assert cs.variants["ins:ref:7:A:AG"].hp_run is None


def test_homopolymer_flag_follows_hp_min(tmp_path):
    from conftest import make_callset

    rows = {"ref": "CGAAAAAACGT", "s1": "CGAAAAA-CGT"}
    (tmp_path / "a").mkdir()
    (tmp_path / "b").mkdir()
    low = make_callset(tmp_path / "a", rows, hp_min=5)
    high = make_callset(tmp_path / "b", rows, hp_min=7)
    assert "homopolymer" in low.calls[0].flags
    assert "homopolymer" not in high.calls[0].flags
    assert low.calls[0].key == high.calls[0].key
