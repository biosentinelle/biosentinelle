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
      strains.tsv             # strain catalogue: one row per strain and the source of its SNPs (generated)
      aliases.tsv             # strain names that normalisation cannot match (curated)
      reads/                  # per-strain consensus VCF and callable BED of the read strains (generated)
    sources/
      read_runs.tsv           # SRA runs of the S288C family, built from or held out (curated)
    panels/
      panel_sites.tsv.gz      # sites, with the strains carrying alt, het or unknown (generated)
      panel_kmers.tsv.gz      # every k-mer of every site, with the allele at its centre (generated)
      panel_strains.tsv       # strain index used by panel_sites, and sites carried (generated)
      panel_report.txt        # what the build kept and dropped, and the strains it cannot separate
    mating/
      mating_kmers.tsv        # MATa- and MATalpha-specific k-mers (generated)
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
| **Strain SNP panel** (`panels/`) | Which strain, among about 1,200 laboratory strains and isolates, the sample is | Exact matching, in the raw reads, of each strain's own k-mers around each site |
| **Mating type** (`mating/`) | MATa, MATalpha or a/alpha | Ratio of Ya- to Yalpha-specific k-mer counts |

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

# sources of the strain panel (downloaded once, into dataset/, not committed)
#   dataset/1011genomes/1011Matrix.gvcf.gz   http://1002genomes.u-strasbg.fr/files/ (md5 in md5.txt)
#   dataset/1011genomes/supp/table_s1.tsv    Table S1 of Peter et al. 2018, from Europe PMC PMC6784862
#   dataset/sgd_strains/                      SGD strain assemblies, s3 sgd-archive.yeastgenome.org sequence/strains/
#   dataset/sra_lab/<run>_[12].fastq.gz       ENA, runs listed in db/yeast/sources/read_runs.tsv
# the population matrix as a sparse table: per record, the samples not hom-ref
tabix -p vcf 1011Matrix.gvcf.gz
for c in $(tabix -l 1011Matrix.gvcf.gz); do
    bcftools query -r $c -i 'GT!="RR"' -f '%CHROM\t%POS\t%REF\t%ALT\t[%SAMPLE:%GT:%AD,]\n' \
        1011Matrix.gvcf.gz | gzip -1
done > dataset/1011genomes/sparse_gt.tsv.gz
bcftools query -l 1011Matrix.gvcf.gz > dataset/1011genomes/samples.txt

# strain catalogue: every strain once, with its best source
python3 bin/build_strain_catalog.py \
    --matrix-table dataset/1011genomes/supp/table_s1.tsv \
    --matrix-samples dataset/1011genomes/samples.txt \
    --reads db/yeast/sources/read_runs.tsv --reads-dir dataset/sra_lab \
    --scrapdb-dir dataset/genome_references/ScRAPdb_20240610_nuclear_genome \
    --sgd-dir dataset/sgd_strains \
    --exclude "<the 33 other-species assemblies of marker_alleles_report.txt>,scrapdb:BY,scrapdb:S288C,scrapdb:SGDref,scrapdb:CENPK1137D" \
    --aliases db/yeast/strains/aliases.tsv \
    --path-map "dataset=/mnt/c/Users/gmillot/Documents/Git_projects/biosentinelle/dataset" \
    --out db/yeast/strains/strains.tsv

# strain panel: assemblies aligned, read runs called, then every strain at once
nextflow run build_db.nf -c nextflow.config

# mating-type k-mers, with the cassette coordinates fetched from Ensembl
python3 bin/build_mating_kmers.py \
    --reference dataset/genome_references/test/S288C_reference_sequence_R64-1-1_20110203.fsa \
    --out db/yeast/mating/mating_kmers.tsv

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

The yeast panel is committed so the pipeline runs out of the box, and it is
fully reproducible: `build_db.nf` regenerates it from the sources listed in
`strains/strains.tsv`. Rebuild it whenever the strain list, the reference or
`panel_kmer_size` changes, since any of the three invalidates it.

## The strain panel

### Where each strain's SNPs come from

A strain does not need an assembly, only its alleles against the reference
and where they are known. Three sources give that, and the catalogue
(`bin/build_strain_catalog.py`) takes each strain once, from the best of them:

| source | strains | why this one |
| :--- | :--- | :--- |
| reads | the S288C family: S288C, FY4, FY1679, BY4741, BY4742, X2180-1A, YPH499 | they differ from each other by tens of SNPs, as many as the base errors of a short-read assembly |
| matrix | the 1,011 isolates of the 1011 Yeast Genomes project (Peter et al. 2018), W303 and SK1 among them | one caller, one reference, explicit no-calls |
| assembly | everything else in ScRAPdb (long reads) and SGD: CEN.PK, Sigma1278b, D273-10B, FL100, JK9-3d, SEY6210, 74D694, industrial and wild isolates | the only source there is |

Names are normalised before comparison (`RM11-1a`, `RM11-1A` and the matrix
isolate `RM11-1a` are one strain), and `aliases.tsv` holds the matches
normalisation cannot make. Two matrix samples sharing a name are both kept:
the 1011 set holds some isolates twice, as the original and as a spore of it.
Left out: the 33 ScRAPdb assemblies of other species, the ScRAPdb `BY`,
`S288C` and `SGDref` assemblies (the read strains replace them) and ScRAPdb
`CENPK1137D`, which carries `leu2-3,112` although CEN.PK113-7D is a
prototroph (SGD's Delft assembly stands for it).

**Read strains.** Each run is aligned end to end with bowtie2, called with
bcftools, and made callable where at least 10 reads of mapping quality 20 cover
a position and depth stays under 2.5 times the median. Where a strain has
runs from independent stocks, `bin/consensus_strain_calls.py` keeps a variant
only when every run calls it identically and every run can see it; anywhere
else the strain is unknown, never reference. A heterozygous call in a haploid
strain is a paralog collapsed onto one copy, and the position is made unknown
too. The runs and the stocks they come from are in `sources/read_runs.tsv`;
the runs marked `heldout` are never used to build, only to test.

**Matrix strains.** The matrix lists every position variable in the
population, with each isolate's genotype and read counts there. It writes many
no-calls over reads that all carry the reference base (`./.:14,0`); a no-call
with at least two reference reads and no other is read as reference. A
position the matrix does not list is taken as reference in every isolate.

### Each strain's own k-mers

A probe k-mer is valid for a strain only if the k-1 bases around the SNP are
the ones that strain carries. With a thousand strains, carrying together a
variant every few bases, almost every site has another strain's variant in
its window, so a single reference/alternative pair per site would leave
almost nothing. `bin/build_population_panel.py` rebuilds the window of every
strain from that strain's own alleles and keeps one k-mer per distinct
haplotype, labelled with the allele at its centre; `scan_kmers.py` sums them
per allele. A strain whose window cannot be rebuilt (an indel in it, a
no-call, more than two heterozygous neighbours) is unknown at that site.
Every k-mer must be absent from the reference genome except the reference
window itself, which must occur once.

### Indels

Short indels, up to `--max-indel` (20) bases gained or lost, are sites too.
Each is left-aligned first, so the matrix and an assembly writing one event at
two positions give one site. The reference and alternative haplotypes are
written out with the allele between 10 bases on each side, and each probe is
the k-mer of its haplotype centred on the allele: the junction of a deletion,
the inserted bases of an insertion. k stays 21, so 35 bp reads still carry the
probes. An indel that only shifts a tandem repeat leaves the reference probe
readable in the alternative allele and is dropped.

A short deletion's reference probe is centred on a deleted base, which makes
it the reference probe of a SNV at that base as well. When a k-mer is shared
that way, the SNV is kept and the indel dropped; any other shared k-mer drops
every site holding it.

`--trace-strain` writes, for every variant one strain carries, the step that
kept or dropped it. For BY4742's 47 variants: 4 kept, 20 carried by more than
20 strains and thinned as common sites, 13 unknown in more than a tenth of
the strains, 7 at positions where other strains carry a different
alternative allele, 2 low complexity. None is lost to indels or collisions:
what limits BY4742 is how little separates it from the rest of its family.

### Which sites are kept

Two clonal relatives differ only at variants rare in the population, often
private to one of them, so for every strain its rarest sites are kept first,
up to `--per-strain-sites`. Common sites, which place a sample in its clade,
are added evenly along the genome. The site table lists, per site, the strains
carrying the alternative allele, the heterozygous ones and the unknown ones;
every other strain carries the reference allele.

The scorer holds each strain as bit sets over the sites, so ranking a sample
against every strain is a few popcounts per strain, not a loop over sites.

### Telling the S288C family apart

The SNP panel places a sample in the S288C family, and inside it the strains
differ by few SNPs or none. Two more measurements decide, and they are the
ones the laboratory derivatives were built to differ by:

* **mating type**, from `mating/mating_kmers.tsv`. Every genome carries a
  silent a cassette (HMR) and a silent alpha cassette (HML), so the presence of
  either says nothing; the ratio does. Per k-mer, a / (a + alpha) is about 2/3
  in a MATa haploid, 1/3 in a MATalpha one and 1/2 in an a/alpha diploid. The
  k-mers are those of the a1, a2, alpha1 and alpha2 genes found in no other
  cassette. On SRR26617475 (BY4741, MATa) the ratio is 0.76.
* **native marker genotype**, as described above: BY4741 `his3D1 leu2D0
  met15D0 ura3D0`, FY1679 `ura3-52` and heterozygous `trp1D63 leu2D1 his3D200`,
  YPH499 `ura3-52 lys2-801 ade2-101 trp1-D63 his3-D200 leu2-D1`, FY4 and
  X2180-1A prototrophs.

`score_strains.py` collects every strain separated from the best one by fewer
than `--min-discriminating-sites` called sites, drops those whose documented
mating type or markers the sample contradicts, and names what is left. A strain
with no documented genotype is compared with the marker alleles of its own
assembly in `markers/marker_alleles.tsv` instead. FY4 (MATa) and S288C
(MATalpha) differ at no marker, so mating type alone parts them; S288C and a
MATalpha prototroph derivative such as FY5 cannot be parted, and the report
says so.

Two more rules decide the call:

* Among the strains nothing contradicts, the one whose genotype was checked
  most fully comes first: a documented mating type and documented markers
  above markers read from an assembly, above nothing. SY14, built from BY4742,
  ties with it on SNPs and markers, but only BY4742's are documented.
* A strain also stands as the parent of the sample when the only differences
  are the two ways laboratory strains are made from one another: more genes
  knocked out (a marker inactive in the sample, active in the strain), and an
  a/alpha diploid of a haploid strain's background (X2180 from X2180-1A and
  -1B, BY4743 from BY4741 and BY4742). The call is then "BY4741-derived (trp1,
  lys2, ade2, can1 inactivated as well)" or "X2180-1A-derived (as an a/alpha
  diploid)". Parents are ranked with the fitting strains by how much of their
  genotype agrees, then by how few changes they need, since a prototroph would
  otherwise fit as the parent of anything.

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

Built 2026-10-28 at k=21: **1,172 strains, 264,652 sites (251,660 SNVs and
12,992 indels), 1,485,095 k-mers** (234,652 rare sites, kept per strain, and
30,000 common ones). 1,011 strains
come from the 1011 matrix, 154 from assemblies (ScRAPdb and SGD), 7 from reads
(the S288C family). One pair of strains has identical expected alleles at every
kept site, IMDO050523 and LW2591Y, two ScRAPdb assemblies.

What the build dropped, of 1,896,949 candidate SNVs and indels: 211,528
common sites thinned by distance, 126,068 unknown in more than a tenth of the
strains, 61,797 low complexity, 21,672 whose k-mers are repeated in the
reference, 25,149 uninformative, 5,934 with a k-mer shared with another
site, 1,934 indels that only shift a repeat, 28 whose strain haplotypes give
one k-mer to both alleles.

Against the SNV-only build of 2026-10-27 (261,469 sites), the indels add
specific sites to 338 strains and take some from 115; S288C, BY4741 and
FY1679 gain one each, BY4742 stays at none. Since the per-strain quota and the
common sites are shared with indels, most strains carry a few dozen fewer
sites in all (median 56 fewer).

Read strains, after the agreement between runs (`strains/reads/`):

| strain | runs | variants kept | dropped, runs disagree | panel sites carried |
| :--- | :--- | :--- | :--- | :--- |
| S288C | 2 | 155 | 1,198 | 15 |
| FY4 | 2 | 135 | 146 | 10 |
| BY4741 | 2 | 146 | 284 | 7 |
| BY4742 | 2 | 46 | 721 | 4 |
| X2180-1A | 1 | 198 | - | 19 |
| FY1679 | 1 | 423 | - | 18 |
| YPH499 | 2 | 4,686 | 897 | 636 |

### Checked against genomes it was not built from

`test/validate_panel.py` simulates 100 bp reads at 20x from a long-read
ScRAPdb assembly of a strain the panel takes from another source (the 1011
matrix, or the reads), and scores them on SNPs alone, without markers. Of 22
strains, 15 rank first and 21 are first or among the strains the report says
it cannot separate from the first:

| expected | ranked first | tied with it |
| :--- | :--- | :--- |
| W303, RM11-1a, YPH499, L1528 and 11 isolates | itself | - |
| SK1 | S799 | SK1 |
| YPS128 | YPS606, a YPS128 derivative | YPS128 and 12 more oak isolates |
| UWOPS03-433.3 | BMC_2a | it and 5 more Malaysian isolates |
| NC_02_b | AGY731 | it and 5 more |
| BY4742 | SY14, built from BY4742 | BY4742 and 16 more of the S288C family |
| S288C | BY4741 | S288C and 16 more of the S288C family |
| DBVPG6044 | Y5 | not DBVPG6044 |

The last row is a disagreement between sources, not a pipeline error. The
ScRAPdb DBVPG6044 assembly matches the West African isolates (Y5, YJM1248,
CBS403) at 99.96 per cent, where DBVPG6044 is classically placed; the 1011
sample labelled DBVPG6044 (AKI, "Mosaic region 3" in Table S1) differs from
that assembly at 3,623 sites. One of the two stocks is not the classic
DBVPG6044. The panel keeps the matrix sample, as for every isolate in the
matrix, and this note is the record of the conflict.

### The honest limit: stocks drift

Two stocks of one strain differ by tens of mutations: the two S288C runs
disagree at 1,198 positions. Inside the S288C family, and between a strain and
a derivative built from it (YPS606 from YPS128, SY14 from BY4742), that is as
many differences as separate two strains. So the SNP panel names the family,
and the strain inside it is named by mating type and markers. The report lists
every strain separated from the best one by fewer than
`--min-discriminating-sites` (200) called sites, drops those whose documented
genotype, or whose own assembly's marker alleles, the sample contradicts, and
says which are left.

### Held-out runs

Five real runs never used to build the panel (`sources/read_runs.tsv`, role
`heldout`, and SRR26617475), through the whole pipeline. In each, the SNP panel
ties 18 strains of the S288C family; mating type and markers make the call.

| run, as labelled | mating type (a fraction) | markers | call |
| :--- | :--- | :--- | :--- |
| SRR26617475, BY4741 | MATa (0.76) | ura3D0 his3D1 leu2D0 met15D0 | BY4741 |
| SRR3279672, BY4741 (a ctf18D replication study) | MATa (0.83) | BY4741's, plus trp1 lys2 ade2 can1 | BY4741-derived (trp1, lys2, ade2, can1 inactivated as well) |
| SRR5431143, BY4742 | MATalpha (0.38) | ura3D0 his3D1 leu2D0 lys2D0 | BY4742 |
| SRR5312396, FY4 | MATa | all eight active | FY4; S288C, the same markers in MATalpha, is excluded |
| DRR505219, X2180 | a/alpha (0.49) | all eight active | X2180-1A-derived (as an a/alpha diploid), which is what X2180 is |

The strains each report still lists as not excluded are those nothing measured
contradicts but whose genotype is less documented: CBS7765, a 1011 isolate only
100 SNPs from the reference and almost certainly an S288C-derived laboratory
strain deposited in a collection, with no marker genotype; and, for the
prototrophs, X2180-1A and the DGY lines of ScRAPdb.
