Yeast integrated individual feature FASTA package v1

Contents
- individual_feature_fastas/: 69 one-record FASTA files for individual integration-relevant features only.
- integrated_features_feature_orientation_69_records.fasta: the same 69 feature records combined into one multi-FASTA.
- manifest_features_only.tsv: source plasmid, feature name/type, coordinates, strand, length, and source URLs.
- checksums.sha256: SHA-256 checksums for files in this cleaned package.

What was removed from the previous larger package
- Full integrated cassette spans.
- Plasmid backbone/source GenBank files.
- Bacterial propagation features such as AmpR, bacterial origin, T7, SP6.

Orientation
- FASTA sequences are in FEATURE orientation.
- For features annotated on the minus strand in the source record, the sequence was reverse-complemented.
- The coordinates and strand in the FASTA header and manifest still refer to the original source plasmid record.

Important limitation
This is not an official exhaustive list of all yeast genome-modification modules. It is the cleaned subset of the 69 individual features from the previously generated classical yeast integration cassette package. Target-locus homology arms are not included.
