import os
import shutil
from pathlib import Path

import pytest

from aubergine.call import call_set
from aubergine.cds import annotate, parse_gff3
from aubergine.compare import state_table
from aubergine.discord import analyse
from aubergine.fasta import Reference, load_alignment
from aubergine.haplotype import build_sample_events
from aubergine.params import default_thresholds

FIXTURES = Path(__file__).parent / "fixtures"


def require_tool(name: str) -> str:
    exe = shutil.which(name)
    if exe is None:
        if os.environ.get("AUBERGINE_REQUIRE_TOOLS") == "1":
            pytest.fail(f"{name} is required (AUBERGINE_REQUIRE_TOOLS=1) but not on PATH")
        pytest.skip(f"{name} not on PATH")
    return exe


def write_msa(path: Path, ref_id: str, rows: dict[str, str]) -> Path:
    path.write_text("".join(f">{k}\n{v}\n" for k, v in rows.items()))
    return path


def make_callset(tmp_path: Path, rows: dict[str, str], ref_id="ref", gff=None, **thresholds):
    """Call an MSA given as {id: aligned row}; the row named ref_id is the reference."""
    ref_seq = rows[ref_id].replace("-", "")
    ref_fa = tmp_path / "ref.fasta"
    ref_fa.write_text(f">{ref_id}\n{ref_seq}\n")
    msa = write_msa(tmp_path / "set.aligned.fasta", ref_id, rows)
    reference = Reference.from_fasta(ref_fa)
    params = {**default_thresholds(), **thresholds}
    cs = call_set(load_alignment(msa, reference), reference, "set", params)
    if gff is not None:
        gff_path = tmp_path / "ann.gff3"
        gff_path.write_text(gff)
        annotate(cs, parse_gff3(gff_path, reference), build_sample_events(cs), str(gff_path))
    return cs


def keys_of(cs, sample: str) -> list[str]:
    """Sample's variant keys with the contig stripped (type:POS:REF:ALT)."""
    out = []
    for c in cs.calls:
        if c.sample == sample:
            t, rest = c.key.split(":", 1)
            out.append(t + ":" + rest.split(":", 1)[1])
    return out


def analysis(cs, **overrides):
    thresholds = {**cs.params, **overrides}
    se = build_sample_events(cs)
    return state_table(cs, se), analyse(cs, se, thresholds), se


@pytest.fixture
def mk(tmp_path):
    def _mk(rows, **kw):
        return make_callset(tmp_path, rows, **kw)

    return _mk
