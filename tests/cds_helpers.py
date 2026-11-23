"""Coding references and an edit-to-alignment helper for the CDS and stop tests."""

_COMP = str.maketrans("ACGT", "TGCA")

# CDS = ATG AAA CCC GGG TTT AAA CCC GGG TAA at reference 3..29 (reference = CDS + 2)
R2 = "CC" + "ATGAAACCCGGGTTTAAACCCGGGTAA" + "CC"
# CDS = ATG GCC CTA AGG GCA TAA at reference 3..20
R1 = "CC" + "ATGGCCCTAAGGGCATAA" + "CC"


def revcomp(seq: str) -> str:
    return seq.translate(_COMP)[::-1]


def gff(start, end, strand="+", attrs="ID=c1;gene=g1", seqid="ref", segments=None):
    segs = segments or [(start, end)]
    return "".join(f"{seqid}\t.\tCDS\t{s}\t{e}\t.\t{strand}\t0\t{attrs}\n" for s, e in segs)


def edit(ref: str, subs=None, dels=(), ins=None) -> dict[str, str]:
    """Reference and one sample ('s1') as an alignment, from 1-based edits.

    subs: {pos: base}; dels: positions; ins: {p: seq} inserted after position p.
    """
    subs, ins = subs or {}, ins or {}
    ref_row, row = [], []
    for q, b in enumerate(ref, 1):
        ref_row.append(b)
        row.append("-" if q in dels else subs.get(q, b))
        if q in ins:
            ref_row.append("-" * len(ins[q]))
            row.append(ins[q])
    return {"ref": "".join(ref_row), "s1": "".join(row)}
