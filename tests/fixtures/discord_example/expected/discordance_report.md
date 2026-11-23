# Discordance report: example

3 samples, 7 variant keys, **3 discordant loci** (2 where callable samples carry different haplotypes; 1 explained only by missing data or ambiguity codes).

Files: `discordance_summary.tsv` (counts per reason), `discordance.tsv` (one row per locus), `discordance_pairs.tsv` (one row per locus and sample pair), `variants.tsv` (per-variant states), `calls.tsv` (per-sample calls).

## Samples

`odd one out`: discordant loci where the sample's haplotype differs from the majority haplotype. `missing`/`ambiguous`: discordant loci where the sample cannot be called. `covered` and `ambiguity codes` are genome-wide.

| sample | covered | ambiguity codes (genome) | odd one out | missing | ambiguous |
|---|--:|--:|--:|--:|--:|
| s1 | 100.0% | 0 | 0 | 0 | 0 |
| s2 | 100.0% | 0 | 0 | 0 | 0 |
| s3 | 87.5% | 0 | 1 | 1 | 0 |

## Sample pairs

Loci where the two samples carry different sequences (missing data, ambiguity and equivalent representations excluded); `+N equiv` = loci where they differ only in how the same sequence is written as variant keys.

| | S1 | S2 | S3 |
|---|--:|--:|--:|
| S1 |  | 1 (+1 equiv) | 2 |
| S2 |  |  | 2 |
| S3 |  |  |  |

S1 = s1; S2 = s2; S3 = s3

## Discordant loci

Most substantive first.

Who carries what: each haplotype's sequence over the locus (`ref` = reference) and its samples; the largest group is given as a count.

| locus | genes | primary reason | flags | who carries what | CDS effect |
|---|---|---|---|---|---|
| L1:16-18 | . | paired_indel | near_end | AGT: 2 samples / GAT (ref): s3 | . |
| L1:6-11 | . | indel_length | near_end;homopolymer | AAAAAA (ref): s3 / AAAA: s2 / AAAAA: s1 | . |
| L1:22-22 | . | missing_data | near_end | C: 2 samples / missing: s3 | . |
