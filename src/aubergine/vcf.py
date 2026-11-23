"""Multi-sample, haploid VCF 4.3 writer."""

from __future__ import annotations

from collections import defaultdict
from pathlib import Path

from . import __version__
from .call import CallSet
from .compare import STATES
from .fasta import ACGT
from .haplotype import apply, as_string, subset_events

HEADER_LINES = [
    '##INFO=<ID=TYPE,Number=A,Type=String,Description="Allele type: snp, ins, del or complex">',
    '##INFO=<ID=LEN,Number=A,Type=Integer,Description="Signed length change of the allele">',
    '##INFO=<ID=LEN_MOD3,Number=A,Type=Integer,Description="abs(LEN) modulo 3">',
    "##INFO=<ID=HP_LEN,Number=A,Type=Integer,"
    'Description="Longest reference homopolymer run containing or flanking the allele">',
    "##INFO=<ID=FS,Number=A,Type=String,"
    'Description="|-joined CDS names in which this allele alone shifts the frame, else .">',
    '##INFO=<ID=AC,Number=A,Type=Integer,Description="Number of samples carrying the allele">',
    '##INFO=<ID=AN,Number=1,Type=Integer,Description="Number of samples with a called genotype">',
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Haploid genotype">',
    "##FORMAT=<ID=ST,Number=1,Type=String,"
    'Description="Sample state: alt, ref, missing, ambiguous, deleted or other_allele">',
]

_STATE_RANK = {s: i for i, s in enumerate(STATES)}


def _vcf_base(b: str) -> str:
    return b if b in ACGT else "N"


def write_vcf(cs: CallSet, states, sample_events, path: Path) -> None:
    samples = cs.sample_ids
    calls_by = defaultdict(list)  # (sample, contig, pos) -> [call]
    for c in cs.calls:
        v = cs.variants[c.key]
        calls_by[(c.sample, v.contig, v.pos)].append(c)

    with open(path, "w") as fh:
        fh.write("##fileformat=VCFv4.3\n")
        fh.write(f"##source=aubergine {__version__}\n")
        if cs.reference.path:
            fh.write(f"##reference={cs.reference.path}\n")
        for contig in cs.reference.contigs.values():
            fh.write(f"##contig=<ID={contig.id},length={contig.length}>\n")
        for line in HEADER_LINES:
            fh.write(line + "\n")
        fh.write(
            "\t".join(
                ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"] + samples
            )
            + "\n"
        )

        groups = defaultdict(list)
        for key, v in cs.variants.items():
            groups[(cs.reference.index(v.contig), v.pos)].append(key)

        for ci, pos in sorted(groups):
            contig = cs.reference.ids[ci]
            ref_seq = cs.reference.seq(contig)
            keys = groups[(ci, pos)]
            ref_len = max(len(cs.variants[k].ref) for k in keys)
            ref = ref_seq[pos - 1 : pos - 1 + ref_len]
            end = pos + ref_len - 1

            sample_allele: dict[str, str] = {}
            allele_info: dict[str, dict] = {}
            for s in samples:
                calls = calls_by.get((s, contig, pos), [])
                if not calls:
                    continue
                template = sample_events[(s, contig)]
                se = subset_events(template, [(cs.variants[c.key], c.ins_seq_raw) for c in calls])
                allele = "".join(_vcf_base(b) for b in as_string(apply(ref_seq, se, pos, end)))
                if allele == ref:
                    continue
                sample_allele[s] = allele
                if allele not in allele_info:
                    vs = [cs.variants[c.key] for c in calls]
                    deltas = defaultdict(int)
                    for v in vs:
                        for x in v.cds:
                            deltas[x["gene"]] += x["delta_cds"]
                    allele_info[allele] = {
                        "TYPE": vs[0].type if len(vs) == 1 else "complex",
                        "HP_LEN": max(v.hp_len for v in vs),
                        "FS": "|".join(g for g, d in deltas.items() if d % 3) or ".",
                    }

            alts = sorted(allele_info, key=lambda a: (len(a), a))
            index = {a: i + 1 for i, a in enumerate(alts)}
            gts = []
            ac = dict.fromkeys(alts, 0)
            an = 0
            for s in samples:
                if s in sample_allele:
                    a = sample_allele[s]
                    ac[a] += 1
                    an += 1
                    gts.append(f"{index[a]}:alt")
                    continue
                st = min((states[k][s] for k in keys), key=_STATE_RANK.__getitem__)
                if st == "ref":
                    an += 1
                    gts.append("0:ref")
                else:
                    gts.append(f".:{st}")

            lens = [len(a) - len(ref) for a in alts]
            per_allele = {
                "TYPE": [allele_info[a]["TYPE"] for a in alts],
                "LEN": lens,
                "LEN_MOD3": [abs(x) % 3 for x in lens],
                "HP_LEN": [allele_info[a]["HP_LEN"] for a in alts],
                "FS": [allele_info[a]["FS"] for a in alts],
                "AC": [ac[a] for a in alts],
            }
            info = ";".join(
                [f"{k}={','.join(map(str, v))}" for k, v in per_allele.items()] + [f"AN={an}"]
            )
            fh.write(
                "\t".join(
                    [contig, str(pos), ".", ref, ",".join(alts), ".", ".", info, "GT:ST"] + gts
                )
                + "\n"
            )
