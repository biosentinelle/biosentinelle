# BioSentinelle strain-identity database

Organism-scoped reference data used to identify the strain behind a set of
submitted reads. Everything that is species-specific lives under `db/<organism>/`,
so adding a new model organism means adding a directory, not editing the pipeline.

```
db/
  yeast/                      # Saccharomyces cerevisiae, proof of concept
    organism.tsv              # organism-level metadata and reference genome id
    tools/
      tools.fasta             # one record per biomolecular tool element
      tools.tsv               # id, name, class, length, source, provenance
    strains/
      strains.tsv             # strain id, name, assembly used, notes
    panels/
      <strain>.snp.tsv        # discriminant SNPs + flanking k-mers (generated)
      panel_index.tsv         # which panels exist, k, number of sites, checksum
```

## Two kinds of evidence

| Evidence | What it answers | How it is detected |
| :--- | :--- | :--- |
| **Tool elements** (`tools/`) | Which markers, tags, resistance cassettes, reporters are present in the submitted sample | Read alignment to each element, scored by breadth of coverage and identity |
| **Element origin** (`tools/tools.tsv`, `origin` column) | Whether a hit is evidence of engineering at all | Fraction of the element's k-mers found in the host reference genome |
| **Strain SNP panels** (`panels/`) | Which laboratory background the sample most resembles | Exact k-mer matching of the reference and alternative allele k-mers in the raw reads |

## Native elements are not evidence

Cassettes are assembled partly from the host's own sequence. In the yeast set,
`TRP1`, the `TRP1` and `GAL1` promoters and the `ADH1` and `CYC1` terminators
are all *S. cerevisiae* sequence; `kanMX`, `NrsR`, `HygR`, `ble`, `GFP` and the
*K. lactis* `URA3` and `TRP1` are not. An unmodified strain covers every native
element end to end at full identity, simply because it carries the native
locus.

On the real SRR26617475 run, the only elements scoring above 90 were exactly
those five native ones, and the strain carries no cassette at all. Reported
without the distinction, that reads as "TRP1 present, GAL1 promoter present",
and anyone would conclude the strain is tagged.

So `origin` splits them, `score_strains.py` counts only foreign elements, and a
native hit raises a warning saying what it does and does not mean. Calling a
native element genuinely present needs reads spanning the junction with an
insertion site, or a copy number above one; neither is implemented yet.

The two kinds of evidence are deliberately independent. Tool elements are foreign sequence, so
presence/absence is a coverage question. Strain identity is a question about
which allele is carried at thousands of informative positions, which exact
k-mer counting answers without any alignment.

## Provenance rule

No sequence in this database is written by hand. Every record is derived by a
script in `bin/` from a cited source file, and the source is recorded in the
`source` and `provenance` columns of the manifests. If you cannot point at where
a sequence came from, it does not belong here.

## Building

```bash
# tool element database, from the curated feature package in dataset/
python3 bin/build_tool_db.py \
    --manifest dataset/features_test/manifest_features_only.tsv \
    --fasta    dataset/features_test/integrated_features_feature_orientation_69_records.fasta \
    --outdir   db/yeast/tools

# then mark which of those elements are the host's own sequence.
# build_tool_db.py rewrites tools.tsv, so this must be run after it, every time.
python3 bin/flag_native_elements.py \
    --tools     db/yeast/tools/tools.fasta \
    --manifest  db/yeast/tools/tools.tsv \
    --reference dataset/genome_references/test/S288C_reference_sequence_R64-1-1_20110203.fsa

# strain SNP panels, from whole-genome assemblies aligned to the reference
nextflow run build_db.nf -c nextflow.config
```

The yeast panels are committed so the pipeline runs out of the box, and they
are fully reproducible: `build_db.nf` regenerates them from the assemblies
listed in `strains/strains.tsv`. Rebuild them whenever the strain list, the
reference or `panel_kmer_size` changes, since any of the three invalidates
every panel.

## Choosing k

The binding constraint is read length, not specificity. A read of length L
yields L - k + 1 k-mers, so a long k on short reads leaves almost nothing, and
a single sequencing error destroys what is left. Measured on SRR26617475,
whose reads are 35 bp:

| k | sites called, of 36694 |
| :--- | :--- |
| 31 | 14 |
| 21 | 34003 |

Specificity is not what limits k for a 12 Mb genome, and `build_snp_panel.py`
checks it empirically anyway by requiring each probe to occur exactly once in
the reference. Keep k well below the read length.

## Current yeast panel

Built 2026-10-23 at k=21: 36694 informative sites over BY4742, SK1, W303, Y55
and S288C. Differences from the reference: W303 8611, SK1 68529, Y55 59853,
BY4742 241. The W303 figure agrees with the published count of roughly 9500,
which is the main external check that the alignment settings are right.

Sites that separate each pair:

| pair | sites |
| :--- | :--- |
| SK1 vs Y55 | 25514 |
| SK1 vs W303 | 21789 |
| W303 vs Y55 | 21473 |
| SK1 vs S288C | 19277 |
| Y55 vs S288C | 19154 |
| W303 vs S288C | 5180 |
| **BY4742 vs S288C** | **42** |

The last row is the honest limit of this panel. BY4742 is an S288C derivative,
so the two are separated by a few dozen sites and no amount of sequencing depth
changes that. `score_strains.py` reports the number of discriminating sites
alongside the ranking and warns below 200, because a probability computed over
36694 sites will read as certainty even when the decision rests on 40 of them.
