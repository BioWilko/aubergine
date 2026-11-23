import re
import subprocess

from conftest import make_callset, require_tool

from aubergine.compare import state_table
from aubergine.haplotype import build_sample_events
from aubergine.vcf import write_vcf

# A1 T2 G3 C4 A5 T6 C7 G8 T9 A10 [ins ins] G11 C12 T13 A14 G15 G16 C17 T18 A19 C20
ROWS = {
    "ref": "ATGCATCGTA--GCTAGGCTAC",
    "s1": "ATGCATAGTGTTGCTAGGCTAC",  # C7A; A10G plus insertion TT after 10
    "s2": "ATGCATGGTA--GCTAGGCTAC",  # C7G
    "s3": "ATGC----TA--GCTAGGCTAC",  # deletion of 5-8, spanning the record at 7
    "s4": "ATGCATCGTA--GCTAGGCTAC",  # reference
}


def _vcf(tmp_path, rows):
    cs = make_callset(tmp_path, rows)
    cs.reference.path = str(tmp_path / "ref.fasta")
    se = build_sample_events(cs)
    path = tmp_path / "out.vcf"
    write_vcf(cs, state_table(cs, se), se, path)
    records = {}
    for line in path.read_text().splitlines():
        if line.startswith("#"):
            continue
        f = line.split("\t")
        info = dict(x.split("=") for x in f[7].split(";"))
        records[int(f[1])] = {"ref": f[3], "alt": f[4], "info": info, "gt": f[9:]}
    return path, records


def test_records(tmp_path):
    _, rec = _vcf(tmp_path, ROWS)
    assert sorted(rec) == [4, 7, 10]
    assert rec[7]["ref"] == "C" and rec[7]["alt"] == "A,G"
    assert rec[7]["gt"] == ["1:alt", "2:alt", ".:deleted", "0:ref"]
    assert rec[7]["info"]["AC"] == "1,1" and rec[7]["info"]["AN"] == "3"
    assert rec[10]["alt"] == "GTT"
    assert rec[10]["info"]["TYPE"] == "complex"
    assert rec[10]["info"]["LEN"] == "2"
    assert rec[4]["ref"] == "CATCG" and rec[4]["alt"] == "C"
    assert rec[4]["info"]["TYPE"] == "del"
    assert rec[4]["gt"][2] == "1:alt"


def test_deletion_at_position_one(tmp_path):
    _, rec = _vcf(tmp_path, {"ref": "AACGTACG", "s1": "A-CGTACG", "s2": "AACGTACG"})
    assert rec[1]["ref"] == "AA" and rec[1]["alt"] == "A"
    assert rec[1]["gt"] == ["1:alt", "0:ref"]


def test_header(tmp_path):
    path, _ = _vcf(tmp_path, ROWS)
    text = path.read_text()
    assert text.startswith("##fileformat=VCFv4.3\n")
    assert "##contig=<ID=ref,length=20>" in text
    assert "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\ts3\ts4" in text


def test_bcftools_accepts_and_finds_nothing_to_realign(tmp_path):
    bcftools = require_tool("bcftools")
    path, _ = _vcf(tmp_path, ROWS)
    view = subprocess.run([bcftools, "view", str(path)], capture_output=True, text=True)
    assert view.returncode == 0, view.stderr
    norm = subprocess.run(
        [bcftools, "norm", "-f", str(tmp_path / "ref.fasta"), "-c", "e", str(path)],
        capture_output=True,
        text=True,
    )
    assert norm.returncode == 0, norm.stderr
    # e.g. "Lines   total/split/joined/realigned/removed/skipped:\t3/0/0/0/0/0"
    m = re.search(r"Lines\s+(\S+):\s+(\S+)", norm.stderr)
    assert m, norm.stderr
    counts = dict(zip(m.group(1).split("/"), m.group(2).split("/"), strict=True))
    assert counts["total"] == "3"
    assert counts["realigned"] == "0"
