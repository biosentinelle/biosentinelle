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
    markers/
      marker_genes.tsv        # the native marker genes, by systematic name (curated)
      marker_loci.tsv|fasta|bed  # ORF +- 1 kb, coordinates from Ensembl (generated)
      marker_repeats.tsv      # locus bases also found elsewhere in the reference (generated)
      marker_alleles.tsv      # allele catalogue, one row per assembly and gene (generated)
      marker_alleles_report.txt  # what the build found, and the documented-strain check
      documented_genotypes.tsv   # published genotypes of lab strains, each with its citation
```

## Two kinds of evidence

| Evidence | What it answers | How it is detected |
| :--- | :--- | :--- |
| **Tool elements** (`tools/`) | Which markers, tags, resistance cassettes, reporters are present in the submitted sample | Read alignment to each element, scored by breadth of coverage and identity |
| **Element origin** (`tools/tools.tsv`, `origin` column) | Whether a hit is evidence of engineering at all | Fraction of the element's k-mers found in the host reference genome |
| **Native marker genes** (`markers/`) | Whether the strain's own URA3, HIS3, LEU2, TRP1, LYS2, ADE2, MET15 and CAN1 are active, inactivated in place, deleted or interrupted, and which known allele each one is | Reads aligned to the reference genome; pileup, depth and junction reads over each locus, judged by the protein the allele encodes |
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

## Native marker genes

The selection markers that matter most in a lab strain are its own genes, and
what has been done to them. They are native, so a coverage test says nothing:
a wild-type strain covers them end to end. Worse, three different situations
look alike to it:

| situation | example | coverage | what decides |
| :--- | :--- | :--- | :--- |
| active, with the strain's own polymorphisms | SK1 LYS2, four missense changes | full | the same protein in active strains |
| inactivated in place | W303 `ura3-1` (G234E), `ade2-1` (E64\*), `leu2-3,112` (two +1 G) | full | the protein: a stop, a frameshift, or a documented form |
| deleted, or interrupted | BY4742 `ura3D0`, `his3D1` (187 bp); YPH499 `ura3-52` (Ty1) | a gap, or junction reads | depth and soft-clipped reads |

Every allele, from an assembly or from reads, is written as events on the
locus and passed through one classifier (`bin/marker_effects.py`), which
applies them and translates the result. Nonsense, frameshift, start loss, a
deleted ORF, an insertion, or an in-frame deletion of a fifth of the protein is
loss of function from the sequence itself. A missense form is never judged on
its own. It is called active when the same protein is carried by catalogued
strains with no documented mutation there, inactive when it is a documented
laboratory allele, and uncertain when no catalogued strain carries it.

### The allele catalogue

`build_marker_alleles.py` aligns every locus to every ScRAPdb assembly with
minimap2 and reads each strain's allele off the alignment: small events from
the cs tag, and deletions, insertions and replacements from the way the
alignment breaks up. Built 2026-10-23: 290 assemblies, 2320 alleles. 33 more
assemblies miss four or more loci outright. They are other Saccharomyces
species in ScRAPdb, too divergent for asm5, and are excluded rather than read
as deletions.

The documented strains are the external check, and every documented allele
is recovered from its assembly, with the expected lesion: BY4742 `his3D1 leu2D0
lys2D0 ura3D0`; W303 `ura3-1` G234E, `his3-11,15`, `leu2-3,112` (c.244_245insG
and c.788_789insG), `trp1-1` E83\*, `ade2-1` E64\*, `can1-100` K47\*; YPH499
`ura3-52` (a 5.9 kb insertion at c.116), `lys2-801` W655\*, `ade2-101` E64\*,
`trp1-D63`, `his3-D200`, `leu2-D1` (161 aa in frame); CEN.PK2-1 `ura3-52 his3D1
leu2-3,112 trp1-289` (Q135\*). A missense allele such as `ura3-1` cannot be
confirmed from sequence; the catalogue is what names it.

The check also caught a mislabelled assembly. ScRAPdb's `CENPK1137D` carries
`leu2-3,112`, as do `CENPKL` and `IMX2600`, although CEN.PK113-7D is a
documented prototroph, so it is left out of `documented_genotypes.tsv`. The
build reports a disagreement; it never resolves one silently.

Two loci have a flank repeated in the genome, a delta element upstream of
LEU2 and of MET15. `marker_repeats.tsv` records them, found from the
reference's own hits outside the locus. They are kept out of deciding where a
gene sits in an assembly, and out of the calls from reads, since the reads of
a repeat are shared out between all its copies.

Foreign sequence inside a locus is matched against the foreign tool elements:
two isolates carry `ura3::natMX6`, named from the inserted sequence.

### Calling from reads

Reads are aligned to the whole reference genome (`MarkerAlign`, bowtie2 in
local mode), not to the loci alone, so the reads of a paralog get a low mapping
quality instead of landing in a marker. `genotype_markers.py` then combines:

* the pileup, for substitutions and small indels. An indel is counted only over
  reads that span it with margin: local alignment clips a read ending near an
  indel instead of giving it a gap, which with 35 bp reads put the real G
  insertions of `leu2-3,112` at 40 to 70 per cent.
* depth, for deletions. Depth places a breakpoint to within a few bases, so a
  catalogued deletion whose two ends both lie within 15 bp is taken as the
  exact event, and its protein with it.
* reads soft-clipped at one place from both sides, for insertions. A clip
  whose bases just continue the reference with a few mismatches comes from a
  dense SNP cluster, not a junction, and is not counted.
* positions with too few good reads, which are unknown, never the reference.
  Where a catalogued allele differs from the reference at such positions, the
  call depends on it. It is made when every allele that fits agrees, and is
  uncertain otherwise.

`test/validate_markers.py` simulates reads from the strain assemblies and
checks the calls against the documented genotypes, or for undocumented
strains against their own assembly. Eight strains, eight
markers each:

| reads | depth | correct | what is missed |
| :--- | :--- | :--- | :--- |
| 150 bp | 40x | 64 of 64 | nothing |
| 100 bp | 25x | 64 of 64 | nothing, including both Ty1 insertions (`ura3-52`) |
| 35 bp | 12x | 62 of 64 | the two `ura3-52` Ty1 insertions |

At 35 bp a Ty1 leaves the locus sequence intact and too few reads this short
cross its junction. In YPH499 the call is uncertain, because 17 bp around the
catalogued insertion site cannot be read and the reads cannot exclude it. In
CEN.PK2-1 nothing showed and URA3 is called active. The report states this
limit whenever reads are shorter than 50 bp. A one-copy deletion in a diploid,
which halves the depth, is not called either.

On SRR26617475 the marker genotype is `his3D1 leu2D0 met15D0 ura3D0`, with
TRP1, LYS2, ADE2 and CAN1 wild type. That is the documented genotype of BY4741
at all eight markers, consistent with the SNP panel placing the sample in the
S288C lineage.

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

# native marker loci, with coordinates fetched from Ensembl and checked against
# the reference: each ORF must start with ATG and end with its only stop
python3 bin/build_marker_loci.py \
    --genes     db/yeast/markers/marker_genes.tsv \
    --reference dataset/genome_references/test/S288C_reference_sequence_R64-1-1_20110203.fsa \
    --outdir    db/yeast/markers

# native marker allele catalogue, from every assembly in the directory
python3 bin/build_marker_alleles.py \
    --loci db/yeast/markers/marker_loci.tsv --loci-fasta db/yeast/markers/marker_loci.fasta \
    --assembly-dir dataset/genome_references/ScRAPdb_20240610_nuclear_genome \
    --reference dataset/genome_references/test/S288C_reference_sequence_R64-1-1_20110203.fsa \
    --documented db/yeast/markers/documented_genotypes.tsv \
    --tools db/yeast/tools/tools.fasta --tools-manifest db/yeast/tools/tools.tsv \
    --paf-dir <scratch dir> --minimap2 'apptainer exec -B /mnt/c <minimap2 image> minimap2' \
    --repeats-out db/yeast/markers/marker_repeats.tsv \
    --out db/yeast/markers/marker_alleles.tsv --report db/yeast/markers/marker_alleles_report.txt
```

Rebuild the marker catalogue after the loci, and after any change to the tool
elements, since inserted sequence is named against them.

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
