# aubergine

Variant calling and **disagreement analysis** for sets of consensus genomes of the same sample
(for example one consensus per workflow), from a reference-anchored multiple sequence alignment.

For every place where the genomes in a set disagree, `aubergine` reports which genomes disagree
and why: missing data, an ambiguity code, the same sequence written differently, an indel of a
different length (in a homopolymer or not), a different allele, and so on. It also writes
per-sample calls (TSV, JSON, VCF), CDS-aware frameshift and premature-stop calls, and,
optionally, truth metrics against one genome in the set.

## Install

`aubergine` needs `mafft` (and `bcftools` for the test suite), which are not on PyPI, so it is
installed into a conda environment. Bioconda is the intended distribution channel; until a
recipe exists, install from a checkout:

```bash
mamba env create -f environment.yml     # or: conda env create -f environment.yml
mamba activate aubergine
pip install -e .
aubergine --version
```

Development: `pip install -e .[dev]`, then `pytest` and `ruff check . && ruff format --check .`.
Tests that need `mafft`/`bcftools` are skipped when the tools are missing, unless
`AUBERGINE_REQUIRE_TOOLS=1` is set (as in CI), in which case they fail.

## Usage

```
aubergine align   -r REF.fasta -o SET.aligned.fasta [--threads N] SET.fasta
aubergine call    -r REF.fasta [--annotation CDS.gff3] [THRESHOLDS] -o PREFIX MSA.fasta
aubergine compare [--truth SAMPLE_ID] [THRESHOLDS] -o PREFIX PREFIX.variants.json
aubergine run     -r REF.fasta -o OUTDIR [--threads N] [--annotation CDS.gff3] [THRESHOLDS] INPUT...

THRESHOLDS = [--hp-min N] [--locus-join-gap N] [--end-dist N]
```

- **`align`** runs `mafft --thread N --6merpair --addfragments SET REF`. The reference is the
  anchor, so insertions appear as reference-gap columns. The output is checked (exit status,
  reference row present and unchanged, every input ID present) and written atomically.
- **`call`** accepts any MSA containing the reference record. The reference row is found by the
  ID of the record in `-r` (a missing ID is an error); its ungapped sequence must equal the
  reference. Outputs: `PREFIX.variants.json` (canonical, embeds the reference and thresholds),
  `PREFIX.calls.tsv` (one row per sample call or ambiguity code: `sample, key, seq_start,
  frameshift, flags`), `PREFIX.vcf`, `PREFIX.summary.tsv`, and with `--annotation`,
  `PREFIX.cds.tsv`.
- **`compare`** reads `variants.json` (no reference needed) and writes
  `PREFIX.discordance_report.md` (start here), `PREFIX.variants.tsv` (one row per variant key:
  attributes, `cds` effect, locus and each sample's state), `PREFIX.pairwise.tsv`,
  `PREFIX.discordance.tsv`, `PREFIX.discordance_pairs.tsv`, `PREFIX.discordance_summary.tsv`, and
  with `--truth`, `PREFIX.truth_metrics.tsv` and `PREFIX.truth_detail.tsv`.
- **`run`** treats every FASTA file given (or found in a given directory: `*.fa`, `*.fasta`,
  `*.fna`, optionally `.gz`) as one set, processes the sets one after another, and writes
  `OUTDIR/<set>/{discordance_report.md, aligned.fasta, variants.json, calls.tsv, variants.vcf,
  [cds.tsv], variants.tsv, pairwise.tsv, discordance*.tsv}`, plus `OUTDIR/summary.tsv` (all sets,
  with a `set` column) and `OUTDIR/run_info.json` (tool and mafft versions, command line, inputs,
  reference SHA-256, thresholds, fixed behaviours, skipped sets). Before aligning, `run` skips,
  with a warning, any set in which fewer than two records contain an A, C, G or T (there is
  nothing to compare); a skipped set gets no directory and no `summary.tsv` rows. `run` does not
  compute truth metrics; run `compare --truth` on a
  set's `variants.json`.

All commands accept `--version`, `-v/--verbose` and `-q/--quiet`, and write deterministic,
sorted output. The reference must currently have a single record (see
[Multi-segment references](#multi-segment-references)).

## Thresholds

| Flag | Default | Meaning |
|---|---|---|
| `--hp-min N` (≥ 2) | 5 | an event is `homopolymer` if its `hp_len ≥ N` |
| `--locus-join-gap N` (≥ 0) | 2 | slide spans separated by ≤ N reference bases join one locus |
| `--end-dist N` (≥ 0) | 10 | a locus within N bp of a sample's first/last covered base gets `near_end` |

Defaults live in `src/aubergine/params.py`. `call` records all three in `variants.json`
(`params`). `compare` uses, for each threshold, **the value given on its command line, otherwise
the value recorded in `variants.json`**; it warns when the two differ, and a JSON without a
recorded value is an error. An `--hp-min` override at `compare` re-derives `homopolymer` from the
stored `hp_len`; the `flags` column in `calls.tsv` keeps the call-time value.

Fixed behaviours (recorded in `run_info.json`): mafft `--6merpair --addfragments`; non-N IUPAC
codes are ambiguous, N and terminal gaps are missing; an insertion touching a deletion is two
events and SNPs are always single-base; indels are always left-normalised; NCBI genetic code 1.

## Conventions

### Coordinates and variant keys

Coordinates are 1-based on the ungapped reference. Events are written VCF-style:

| Event | POS / REF / ALT |
|---|---|
| SNP at `p` | `p`, `ref[p]`, `alt` |
| deletion of `s..e` | `s-1`, `ref[s-1..e]`, `ref[s-1]`; when `s = 1`: `1`, `ref[1..e+1]`, `ref[e+1]` |
| insertion of `I` after `p` | `p`, `ref[p]`, `ref[p]+I`; when `p = 0`: `1`, `ref[1]`, `I+ref[1]` |

The **variant key** is `type:contig:POS:REF:ALT` with type `snp`, `ins` or `del`, e.g.
`del:MN908947.3:10:GT:G`. (The legacy pipeline used positions without a contig.) Keys are parsed
from both ends, so contig IDs containing `:` are fine.

In TSV/JSON, `start`/`end` describe the event itself: the SNP position, the deleted bases, or for
an insertion `start = p` and `end = p + 1` (the bases it sits between).

### Calling

Each alignment column is classified per sample:

| reference | sample | class |
|---|---|---|
| `-` | `-` | ignored (does not split runs) |
| base | same base | match |
| A/C/G/T | different A/C/G/T | SNP |
| base | IUPAC code other than N | ambiguous call (reported, not a variant) |
| base | `N`, or a gap before the sample's first / after its last base | missing |
| base | internal gap | deletion |
| `-` | any character, inside the sample's covered span | insertion |
| `-` | any character, outside the covered span | overhang (counted only) |
| N/IUPAC | base | uncallable (no SNP) |

Consecutive deletion (or insertion) columns form one event; an event can start in the column
right after another ends. An insertion made only of `N` is a *missing insertion*, not a variant.
An insertion containing N or IUPAC codes stays a variant: ALT has `N` in their place, the original
sequence is kept in `variants.json` (`calls[].ins_seq_raw`) and the call is flagged
`contains_ambiguity`.

**Left-normalisation** shifts a deletion left while `ref[s-1] == ref[e]`, and an insertion left
while `ref[p] == I[-1]`, but each step only if the sample itself is plain reference at the newly
included base (no SNP, N, ambiguity, deletion or other insertion of that sample there). So keys do
not depend on where mafft put a gap, and normalisation never changes the sample's sequence. Two
calls are the same variant if and only if their keys are identical.

Per-call annotation (all in `variants.json`; `seq_start` and `flags` also in `calls.tsv`):
`length` (signed), `len_mod3`, `hp_len`/`hp_run` (longest reference
homopolymer containing or directly flanking the event; for an insertion, the run touching `p` or
`p+1` whose base occurs in the inserted sequence), `seq_start` (1-based position in the sample's
ungapped sequence) and `flags`: `homopolymer`, `contains_ambiguity`, `near_missing` (missing data
directly next to the event) and `cds_boundary`.

### States (`variants.tsv`, VCF `FORMAT/ST`)

For each variant V and sample S, the first matching rule:

1. `alt`: S has V's key.
2. `missing`: missing data in S overlaps V's callable span (SNP: its base; deletion: the deleted
   bases; insertion: both flanking bases), or S has a missing insertion at V's point.
3. `ambiguous`: S has an ambiguity code in the callable span.
4. `deleted`: a different deletion of S covers part of the callable span.
5. `other_allele`: S has a different SNP in the span, or a different insertion at the same point
   (or strictly inside a deleted span).
6. `ref`.

`pairwise.tsv` counts, per pair and type (`snp`, `insertion`, `deletion`, `all`), keys that both,
only one, or neither sample carries, among keys where neither sample is `missing`/`ambiguous`
(`comparable`; the rest are `excluded`). `jaccard = shared/(shared+only_a+only_b)`, `NA` when 0/0.

## Disagreement analysis

### Loci and haplotypes

Each indel has a **slide span**: every reference position it could occupy in its repeat or
homopolymer, computed from the reference alone. Per contig, slide spans (SNPs: their position)
separated by at most `--locus-join-gap` bases are merged into a **locus** `L<contig
index>:<start>-<end>`. So two indels in the same homopolymer always share a locus, whatever their
length or alignment.

For each locus, each sample's events are applied to the reference window to give its **local
haplotype**, which is labelled `H0` (reference), `H1`, `H2`, … (most common first, ties by
sequence), `missing` (missing data or a missing insertion in the window) or `ambiguous` (an
ambiguity code from the sample in the window). A locus is **discordant** if at least two samples
carry different sets of keys there. `haplotype_concordant` means all callable samples share one
label; `callable_discordant` means at least two callable samples have different labels (the
"real disagreement" filter).

### Reasons

For every pair of samples with different keys at a locus (`Da`/`Db` = keys only in a/b), the
first matching reason:

| # | Reason | Condition |
|---|---|---|
| 1 | `missing_data` | either haplotype is `missing` |
| 2 | `ambiguity` | either haplotype is `ambiguous` |
| 3 | `equivalent_representation` | identical haplotypes (different decomposition or placement) |
| 4 | `paired_indel` | equal haplotype lengths, and `Da` or `Db` contains both an insertion and a deletion |
| 5 | `indel_length` | only indels differ, both sides have some, haplotype lengths differ |
| 6 | `indel_sequence` | only indels differ, both sides have some, same haplotype length |
| 7 | `indel_presence` | only indels differ, one side has none |
| 8 | `different_allele` | only SNPs differ, and both sides have an SNP at some position |
| 9 | `snp_presence` | only SNPs differ, no shared position |
| 10 | `mixed` | anything else (SNP and indel differences) |

A locus lists all its pair reasons (`reasons`) and a `primary_reason`, which is the first one
present in this ranking, so that real differences come before artefacts and missing data:
`paired_indel, indel_length, indel_sequence, indel_presence, different_allele, snp_presence,
mixed, equivalent_representation, ambiguity, missing_data`.

Flags (per pair, and their union per locus): `near_end` (the locus is within `--end-dist` of
either sample's first or last covered base, where `--addfragments` placement is least
reliable), `homopolymer` (a differing event has `hp_len ≥ --hp-min`; "homopolymer indel length"
is `indel_length` + `homopolymer`), `frameshift` (a differing call shifts a CDS frame).

Outputs:

- `discordance_report.md`: the first place to look. A headline count, then per sample (coverage,
  genome-wide ambiguity codes, how often it is the odd one out, missing or ambiguous at
  discordant loci), a matrix of real differences per sample pair (missing data, ambiguity and
  equivalent representations excluded; `+N equiv` counts the last; above 12 samples, the 20 pairs
  with most differences), and the discordant loci, most substantive first (up to 50), with who
  carries which sequence and the coding effect.
- `discordance.tsv`: one row per discordant locus with its keys, flags, overlapping genes,
  haplotype sequences (`H0=…;H1=…`, long ones abbreviated), `odd_samples` (callable samples
  differing from the unique majority haplotype; `.` on a tie), `cds_effect` and each sample's
  label. `cds_effect` lists `gene:fs(samples)` / `gene:fs_compensated(samples)` for frameshifting
  calls at the locus and `gene:premature_stop(samples)` when a sample's call at the locus overlaps
  any base of its first stop codon (a stop caused further downstream by a frameshift is not
  attributed; the locus shows `fs`).
- `discordance_pairs.tsv`: one row per discordant locus × disagreeing pair, with the reason and
  the keys only in a / only in b (coordinates: join `locus` to `discordance.tsv`).
- `discordance_summary.tsv`: per reason, the number of loci where it is primary (and of those,
  callable-discordant), the number of loci where it occurs at all, and the number of pairs.

### Worked example

Reference `ACGTCAAAAAAGTCCGATGCATGC` (A6–A11 is a homopolymer of 6). s1: one A deleted from
A6–A11, G16 deleted and G inserted after A17, T22C. s2: two As deleted, G16A and A17G, T22C. s3:
reference except N at 22–24. (`tests/fixtures/discord_example`)

| Locus | s1 | s2 | s3 | Pair reasons | primary | callable_discordant |
|---|---|---|---|---|---|---|
| `L1:6-11` | H2 `AAAAA` | H1 `AAAA` | H0 | s1–s2 `indel_length`; s1–s3, s2–s3 `indel_presence` | `indel_length` | yes |
| `L1:16-18` | H1 `AGT` | H1 `AGT` | H0 `GAT` | s1–s2 `equivalent_representation`; s1–s3 `paired_indel`; s2–s3 `snp_presence` | `paired_indel` | yes |
| `L1:22-22` | H1 `C` | H1 `C` | missing | s1–s3, s2–s3 `missing_data` | `missing_data` | no |

All three loci are `near_end` (the example is short); `L1:6-11` is also `homopolymer`. Key-based
`pairwise.tsv` gives s1 vs s2 a Jaccard of 0.1429, but `discordance_pairs.tsv` shows the four
keys at `L1:16-18` are the same sequence, and the report's pair matrix gives s1–s2 as
`1 (+1 equiv)`. s3 is the `odd_samples` entry at `L1:16-18`; `L1:6-11` (three haplotypes, one
sample each) has no majority.

## Truth metrics

`compare --truth SAMPLE_ID` scores every other sample in the set against that sample. A truth
variant and a query call match only if their keys are identical.

| Class | Definition |
|---|---|
| TP | truth variant whose state in the query is `alt` |
| FN | truth variant whose state in the query is `ref`, `deleted` or `other_allele` |
| uncalled | truth variant whose state in the query is `missing`/`ambiguous` (not counted as FN) |
| FP | query call not in the truth, whose state in the truth is `ref`, `deleted` or `other_allele` |
| truth_uncallable | query call not in the truth, whose state in the truth is `missing`/`ambiguous` (not counted as FP) |

`precision = TP/(TP+FP)`, `recall = TP/(TP+FN)`, `F1 = 2TP/(2TP+FP+FN)`, `NA` when the
denominator is 0. `truth_detail.tsv` gives every item's locus and the pair reason between truth
and query there; `equivalent` in `truth_metrics.tsv` counts FP/FN items whose reason is
`equivalent_representation` (same sequence, different keys; not credited).

To use a **truth VCF**, turn it into a genome and add it to the set:
`bcftools consensus -f REF.fasta truth.vcf.gz > truth.fasta`, then pass its record ID to
`--truth`.

Worked example (`tests/fixtures/truth_example`): reference `ACGTAACCGGTTTACGATCG`; truth has
G3A, A5G, a T deleted from TTT (`del:10:GT:G`) and TT inserted after C15; the query has A5G, the
same T deletion placed at a different column, N at 15–16 and T18C. Result: SNP TP 1 / FP 1 / FN 1,
insertion uncalled 1, deletion TP 1; overall precision = recall = F1 = 0.6667. The FP T18C shares
locus `L1:11-18` with the insertion the query cannot see, so its detail reason is `missing_data`.

## VCF

`variants.vcf` is VCF 4.3, haploid, multi-sample. Variants sharing a POS form one record; REF is
the longest REF among them, and each sample's allele is its sequence over that span (an SNP and
an insertion at the same POS in one sample give one allele with `TYPE=complex`). ALT alleles are
ordered by (length, sequence); ambiguity codes never become ALT alleles. INFO: `TYPE`, `LEN`,
`LEN_MOD3`, `HP_LEN`, `FS` (CDSs in which that allele alone shifts the frame), `AC`, `AN`.
FORMAT: `GT` (allele index, `0`, or `.`) and `ST` (the state). `deleted` is written as `GT=.`
with `ST=deleted`, not as a `*` allele, because many tools miscount `*`. The file is plain text;
`bgzip` and `tabix` it if needed.

## CDS analysis (`--annotation`)

GFF3 only. Rows of type `CDS` sharing an `ID` (or `Parent`) form one CDS; segments are ordered in
transcription order, so multi-segment CDSs and ribosomal-slippage overlaps (SARS-CoV-2 ORF1ab)
are handled. The label is `gene`, else `Name`, else `ID`; labels used by more than one CDS
become `label[ID]`. Every `seqid` must be a reference record.

**Frameshifts.** A deletion changes the CDS length by minus the number of CDS coordinates it
removes (a slippage base counts twice). An insertion is coding only if the bases on either side
are consecutive in the CDS, so insertions just before the start or after the stop do not count.
Per call (`frameshift` column, `gene:value`): `inframe`, `fs`, or `fs_compensated` (the sample's
net change over that CDS is a multiple of 3). Per sample × CDS (`cds.tsv`): `net_delta`, `status`
(`intact`, `inframe`, `restored` when compensating indels restore the frame, `frameshift`),
`out_of_frame` (reference spans where the frame is shifted) and `incomplete` (missing data
overlaps the CDS). A deletion over the first or last codon is flagged `cds_boundary`.

**Premature stops.** The reference CDS is checked first: its length must be a multiple of 3,
it must end with a stop and have no internal stop, and `transl_table` must be 1. Otherwise
`stop_status = not_evaluated` with a `note`. The sample's CDS sequence is then rebuilt from the
reference and its calls (missing positions become N; only coding insertions are included) and
read codon by codon. A codon containing N or any IUPAC code is unknown, never a stop.

| `stop_status` | Meaning |
|---|---|
| `premature_stop` | the first stop starts before the reference stop |
| `normal` | the first stop starts at (or after) the reference stop position |
| `stop_lost` | no stop in the CDS, and no unknown codons |
| `unknown` | no stop found, but unknown codons could hide one |
| `not_evaluated` | the reference CDS failed validation |

Also reported: `stop_aa_pos`, `protein_len`, `ref_protein_len`, `stop_ref_pos`,
`stop_incomplete` (an unknown codon before the stop, or a missing insertion in the CDS) and
`stop_frameshift` (whether the premature stop lies in an out-of-frame region; for `stop_lost`,
whether the frame stays shifted to the CDS end). Frameshift status and stop status are reported
independently: a `restored` CDS can still carry a `premature_stop` inside its out-of-frame
region, so "restored" alone does not mean the protein is intact.

## Multi-segment references

Only single-record references are accepted for now; the check is in one place (`cli.py`). The
data model is contig-aware throughout (keys, intervals, JSON, VCF contig lines, GFF seqids, the
per-contig alignment and calling units), so adding segmented references mostly needs a step that
assigns sample records to segments.

## Known limitations

- Left-normalisation unifies placements of the same indel, but not different decompositions of
  the same change (e.g. `del+SNP` vs `SNP+del`). Those keep different keys in `variants.tsv`,
  pairwise and truth tables, and are reported as `equivalent_representation` by the
  disagreement analysis.
- A large deletion makes a large locus, which can absorb nearby variants.
- Sample bases beyond the reference ends (overhang) are counted, not called.
- Stop detection does not read past the CDS end to find a read-through stop, and does not expand
  IUPAC codes.

## Mapping from the legacy pipeline

`bin/msa_call_variants.py` TSV → `calls.tsv` joined on `key` to `variants.tsv`:

| Legacy column | New column(s) | Notes |
|---|---|---|
| `sample` | `sample` | every non-reference sequence is now reported, including those with no calls |
| `position` | `pos` (VCF-anchored; `start`/`end` in `variants.json`) | legacy insertions used the next reference base, deletions the first deleted base |
| `ref_allele`, `alt_allele` | `ref`, `alt` | anchored alleles; legacy used `.` for the empty side |
| `variant_type` | `type` (`snp`/`ins`/`del`) | `N` is no longer an SNP (see `missing`); IUPAC codes are `ambiguous:…` rows in `calls.tsv` |
| `seq_position` | `seq_start` | |
| `context.homopolymer` | `hp_len`, `homopolymer` flag | |
| `context.frameshift` | `cds`, `frameshift` | CDS-aware; needs `--annotation` |

`bin/process_variants.py` `present`/`absent`/`low_coverage` → `variants.tsv` `alt` /
`ref`, `deleted`, `other_allele` / `missing`, `ambiguous`.

Bugs in the legacy caller that the new caller avoids: a deletion directly followed by an
insertion (or the reverse) gave a bogus `-` SNP and a truncated or lost indel; an SNP right after
an insertion overwrote it; a trailing indel was dropped; a deletion was split by another sample's
insertion column; terminal gaps were called as deletions; insertions starting with `N` and every
`N` became SNPs; the reference row was assumed to be first; the homopolymer check was unreliable;
"frameshift" ignored coding sequence.

## License

MIT
