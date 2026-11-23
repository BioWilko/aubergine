import json
import logging

import pytest
from cds_helpers import R2, edit, gff
from conftest import FIXTURES, require_tool

from aubergine import __version__
from aubergine.cli import main

EXAMPLE = FIXTURES / "discord_example"


def _call(out, *extra):
    return main(
        [
            "-q",
            "call",
            "-r",
            str(EXAMPLE / "ref.fasta"),
            "-o",
            str(out),
            *extra,
            str(EXAMPLE / "example.aligned.fasta"),
        ]
    )


def test_version(capsys):
    with pytest.raises(SystemExit) as e:
        main(["--version"])
    assert e.value.code == 0
    assert __version__ in capsys.readouterr().out


def test_call_and_compare_outputs(tmp_path):
    out = tmp_path / "ex"
    assert _call(out) == 0
    for name in ("variants.json", "calls.tsv", "vcf", "summary.tsv"):
        assert (tmp_path / f"ex.{name}").stat().st_size > 0
    assert (tmp_path / "ex.summary.tsv").read_text().startswith("set\tsample\t")
    assert not (tmp_path / "ex.cds.tsv").exists()
    assert main(["-q", "compare", "--truth", "s1", "-o", str(out), f"{out}.variants.json"]) == 0
    for name in (
        "variants.tsv",
        "pairwise.tsv",
        "discordance_report.md",
        "discordance.tsv",
        "discordance_pairs.tsv",
        "discordance_summary.tsv",
        "truth_metrics.tsv",
        "truth_detail.tsv",
    ):
        assert (tmp_path / f"ex.{name}").exists(), name


def test_call_with_annotation(tmp_path):
    rows = edit(R2, dels={10})
    (tmp_path / "ref.fasta").write_text(f">ref\n{R2}\n")
    (tmp_path / "aln.fasta").write_text("".join(f">{k}\n{v}\n" for k, v in rows.items()))
    (tmp_path / "ann.gff3").write_text(gff(3, 29))
    assert (
        main(
            [
                "-q",
                "call",
                "-r",
                str(tmp_path / "ref.fasta"),
                "--annotation",
                str(tmp_path / "ann.gff3"),
                "-o",
                str(tmp_path / "x"),
                str(tmp_path / "aln.fasta"),
            ]
        )
        == 0
    )
    lines = (tmp_path / "x.cds.tsv").read_text().splitlines()
    row = dict(zip(lines[0].split("\t"), lines[1].split("\t"), strict=True))
    assert (row["status"], row["stop_status"], row["stop_frameshift"]) == (
        "frameshift",
        "stop_lost",
        "yes",
    )
    assert "FS=g1" in (tmp_path / "x.vcf").read_text()


def test_thresholds_recorded_and_used(tmp_path, caplog):
    out = tmp_path / "ex"
    assert _call(out, "--hp-min", "7", "--locus-join-gap", "4", "--end-dist", "0") == 0
    params = json.loads((tmp_path / "ex.variants.json").read_text())["params"]
    assert params == {"hp_min": 7, "locus_join_gap": 4, "end_dist": 0}

    caplog.set_level(logging.INFO, logger="aubergine")
    assert main(["compare", "-o", str(out), f"{out}.variants.json"]) == 0
    assert "hp_min=7, locus_join_gap=4, end_dist=0" in caplog.text
    assert "WARNING" not in caplog.text
    flags = (tmp_path / "ex.discordance.tsv").read_text()
    assert "homopolymer" not in flags and "near_end" not in flags

    caplog.clear()
    assert main(["compare", "--hp-min", "5", "-o", str(out), f"{out}.variants.json"]) == 0
    assert "using command-line value 5 instead of 7" in caplog.text
    assert "homopolymer" in (tmp_path / "ex.discordance.tsv").read_text()


def test_json_without_threshold_is_an_error(tmp_path, caplog):
    out = tmp_path / "ex"
    assert _call(out) == 0
    path = tmp_path / "ex.variants.json"
    d = json.loads(path.read_text())
    del d["params"]["end_dist"]
    path.write_text(json.dumps(d))
    assert main(["compare", "-o", str(out), str(path)]) == 1
    assert "no recorded value for threshold 'end_dist'" in caplog.text


@pytest.mark.parametrize(
    "bad", [["--hp-min", "1"], ["--locus-join-gap", "-1"], ["--end-dist", "x"]]
)
def test_invalid_threshold_values(tmp_path, bad):
    with pytest.raises(SystemExit) as e:
        _call(tmp_path / "ex", *bad)
    assert e.value.code == 2


def test_run(tmp_path):
    require_tool("mafft")
    sets = tmp_path / "sets"
    sets.mkdir()
    ref = R2
    (tmp_path / "ref.fasta").write_text(f">ref\n{ref}\n")
    (sets / "a.fasta").write_text(f">a1\n{ref}\n>a2\n{ref[:9]}{ref[10:]}\n")
    (sets / "b.fa").write_text(f">b1\n{ref[:5]}T{ref[6:]}\n>b2\n{ref}\n")
    (tmp_path / "ann.gff3").write_text(gff(3, 29))
    out = tmp_path / "out"
    argv = [
        "-q",
        "run",
        "-r",
        str(tmp_path / "ref.fasta"),
        "--annotation",
        str(tmp_path / "ann.gff3"),
        "--hp-min",
        "6",
        "-o",
        str(out),
        str(sets),
    ]
    assert main(argv) == 0
    for name in ("a", "b"):
        for f in (
            "aligned.fasta",
            "variants.json",
            "calls.tsv",
            "variants.vcf",
            "cds.tsv",
            "variants.tsv",
            "pairwise.tsv",
            "discordance_report.md",
            "discordance.tsv",
            "discordance_pairs.tsv",
            "discordance_summary.tsv",
        ):
            assert (out / name / f).exists(), (name, f)
        assert not (out / name / "summary.tsv").exists()
    summary = (out / "summary.tsv").read_text().splitlines()
    assert [ln.split("\t")[:2] for ln in summary[1:]] == [
        ["a", "a1"],
        ["a", "a2"],
        ["b", "b1"],
        ["b", "b2"],
    ]
    info = json.loads((out / "run_info.json").read_text())
    assert info["thresholds"] == {"hp_min": 6, "locus_join_gap": 2, "end_dist": 10}
    assert info["mafft_version"]
    assert info["fixed_behaviours"]["genetic_code"] == "NCBI table 1"
    assert json.loads((out / "a" / "variants.json").read_text())["params"]["hp_min"] == 6
    cds = (out / "b" / "cds.tsv").read_text()
    assert "premature_stop" in cds
    disc = (out / "b" / "discordance.tsv").read_text()
    assert "g1:premature_stop(b1)" in disc


def test_run_skips_sets_without_two_informative_records(tmp_path, caplog):
    require_tool("mafft")
    sets = tmp_path / "sets"
    sets.mkdir()
    ref = R2
    (tmp_path / "ref.fasta").write_text(f">ref\n{ref}\n")
    (sets / "ok.fasta").write_text(f">o1\n{ref}\n>o2\n{ref[:5]}T{ref[6:]}\n>o3\n{'N' * len(ref)}\n")
    (sets / "one.fasta").write_text(f">x1\n{ref}\n>x2\n{'N' * len(ref)}\n>x3\nNNRN--\n")
    (sets / "none.fasta").write_text(f">y1\n{'N' * len(ref)}\n")
    out = tmp_path / "out"
    assert main(["run", "-r", str(tmp_path / "ref.fasta"), "-o", str(out), str(sets)]) == 0
    assert (out / "ok" / "discordance_report.md").exists()
    assert not (out / "one").exists() and not (out / "none").exists()
    assert "set one: skipped, only 1 of 3 records have bases other than N" in caplog.text
    summary = (out / "summary.tsv").read_text().splitlines()[1:]
    assert {ln.split("\t")[0] for ln in summary} == {"ok"}
    info = json.loads((out / "run_info.json").read_text())
    assert [(s["set"], s["records"], s["with_content"]) for s in info["skipped"]] == [
        ("none", 1, 0),
        ("one", 3, 1),
    ]
