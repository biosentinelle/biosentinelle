Yeast common molecular-biology sequence FASTA set (exact plasmid-derived records)
===============================================================================

Scope
-----
This package does NOT claim to contain every possible yeast molecular biology sequence.
Names such as GFP, kanMX, GAL1 promoter, ADH1 terminator, 13Myc, FLAG, GST, AmpR, and ori
have multiple variants and/or variable boundaries. The FASTA files here are exact sequences
extracted from the downloaded GenBank records listed in source_genbank/.

Contents
--------
- full_plasmids_fasta/: one FASTA per full plasmid.
- annotated_features_fasta/: one FASTA per annotated GenBank feature; complement(...) features
  are reverse-complemented so the FASTA is in the feature's annotated functional orientation.
- all_exact_records.fasta: all full plasmids + extracted features in one multi-FASTA.
- all_annotated_features_only.fasta: extracted features only.
- manifest.tsv: coordinates, strand, length, source notes, and caveats.
- source_genbank/: original downloaded GenBank files used for extraction.

Important caveat
----------------
pFA6a-13Myc-kanMX6 has an explicit GenBank feature for only one Myc repeat (66..95).
I also included an INFERRED full 13Myc ORF (66..590, no stop codon) because the plasmid name
and upstream ORF contain 13 Myc epitope repeats. Treat this one as inferred, not source-annotated.

Exclusions
----------
I did not include arbitrary FASTA records for broad sequence families from the previous answer
(native promoters such as TDH3p/TEF1p, native terminators, auxotrophic markers, tet systems,
CRISPR guide spacers, YKO barcodes, etc.). For those, an exact plasmid/accession or genomic
coordinate boundary is required before a reliable FASTA can be made.
