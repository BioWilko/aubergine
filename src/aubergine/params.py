"""Single source of defaults: tunable thresholds and fixed behaviours."""

HP_MIN = 5
LOCUS_JOIN_GAP = 2
END_DIST = 10

# name -> (default, minimum allowed value)
THRESHOLDS = {
    "hp_min": (HP_MIN, 2),
    "locus_join_gap": (LOCUS_JOIN_GAP, 0),
    "end_dist": (END_DIST, 0),
}


def default_thresholds() -> dict[str, int]:
    return {name: default for name, (default, _) in THRESHOLDS.items()}


MAFFT_ARGS = ["--6merpair", "--addfragments"]

FIXED_BEHAVIOURS = {
    "mafft_mode": " ".join(MAFFT_ARGS),
    "iupac": "non-N IUPAC codes are ambiguous; N and terminal gaps are missing",
    "ins_touching_del": "two events; SNPs stay single-base",
    "normalisation": "always left-normalise",
    "genetic_code": "NCBI table 1",
}

# Sequences longer than this are abbreviated in discordance.tsv
HAPLOTYPE_DISPLAY_MAX = 50
HAPLOTYPE_DISPLAY_FLANK = 20
