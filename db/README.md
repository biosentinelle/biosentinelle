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
| **Strain SNP panels** (`panels/`) | Which laboratory background the sample most resembles | Exact k-mer matching of the reference and alternative allele k-mers in the raw reads |

The two are deliberately independent. Tool elements are foreign sequence, so
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

# strain SNP panels, from whole-genome assemblies aligned to the reference
nextflow run build_db.nf -c nextflow.config
```

The SNP panels are not committed as sequence; they are generated from the
assemblies listed in `strains/strains.tsv` so that they can be rebuilt and
versioned against a named reference.
