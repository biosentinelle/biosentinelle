<img src="logo.png" align="center" width="160" alt="lab_strain_watch logo">

# lab_strain_watch

| Usage | Requirement |
| :--- | :--- |
| [![Nextflow](https://img.shields.io/badge/code-Nextflow-blue?style=plastic)](https://www.nextflow.io/) | [![Dependencies: Nextflow Version](https://img.shields.io/badge/Nextflow-v24.10.4-blue?style=plastic)](https://github.com/nextflow-io/nextflow) |
| [![License: GPL-3.0](https://img.shields.io/badge/licence-GPL%20(%3E%3D3)-green?style=plastic)](https://www.gnu.org/licenses) | [![Dependencies: Apptainer Version](https://img.shields.io/badge/Apptainer-v1.3.5-blue?style=plastic)](https://github.com/apptainer/apptainer) |
| | [![Dependencies: Graphviz Version](https://img.shields.io/badge/Graphviz-v2.42.3-blue?style=plastic)](https://www.graphviz.org/download/) |

<br><br>
## TABLE OF CONTENTS


   - [AIM](#aim)
   - [WARNINGS](#warnings)
   - [CONTENT](#content)
   - [HOW IT WORKS](#how-it-works)
   - [INPUT](#input)
   - [STRAIN IDENTITY DATABASE](#strain-identity-database)
   - [HOW TO RUN](#how-to-run)
   - [OUTPUT](#output)
   - [VERSIONS](#versions)
   - [LICENCE](#licence)
   - [CITATION](#citation)
   - [CREDITS](#credits)
   - [ACKNOWLEDGEMENTS](#acknowledgements)
   - [WHAT'S NEW IN](#whats-new-in)

<br><br>
## AIM

Check the identity of a sequenced laboratory strain, i.e., answer "is this tube really the strain written on its label?". From the sequences of a single sample (reads or assembly), lab_strain_watch reports:
- **Yeast**
    - the most probable strain among the strains of a database (1,172 *Saccharomyces cerevisiae* strains for the yeast database), or "unassigned" if none fits well enough,
    - the mating type (MATa, MATalpha or MATa/alpha),
    - the genotype of the native auxotrophic and counter-selection marker genes (URA3, HIS3, LEU2, TRP1, LYS2, ADE2, MET15, CAN1, ARG4, LEU5),
    - the status of the mitochondrial genome (rho+, rho- or rho0), hence whether the strain is expected to grow on glycerol,
    - the laboratory tool elements present (plasmids, cassettes, drug resistance markers, tags, promoters, ...).

<br><br>
## WARNINGS

- Currently dedicated to *Saccharomyces cerevisiae*: only the yeast database (`db/yeast` of [yeast_db](https://github.com/lab-strain-watch/yeast_db)) is provided. Another organism needs its own database, built with the *build_db.nf* file of yeast_db (see [STRAIN IDENTITY DATABASE](#strain-identity-database)).
- The strain probability is a posterior over the strains of the database, which are assumed to be exhaustive. It is not the probability that the sample *is* that strain. The concordance (fraction of the expected alleles found) is the measure of fit. See the [OUTPUT](#output) section.
- A single sample per run. All the files of the `fastq_path` folder are analyzed together as one sample.
- Laboratory stocks drift: a few mismatching sites between a sample and its documented strain are expected (see [here](https://github.com/lab-strain-watch/yeast_db/tree/agent/db)).

<br><br>
## CONTENT

| Files and folder | Description |
| :--- | :--- |
| **main.nf** | File that can be executed using a linux terminal, a MacOS terminal or Windows 10 WSL2. Analysis of a sample. |
| **nextflow.config** | Parameter settings for the *main.nf* file. Users have to open this file, set the desired settings and save these modifications before execution. Of note, this configuration file is systematically saved in the *reports* folder (see [below](#output)) during each execution, to save the parameter settings. |
| **bin folder** | Contains the scripts (Python, R, Rmarkdown) required by the *main.nf* file. |
| **modules folder** | Contains the Nextflow processes required by the *main.nf* file. |
| **conf folder** | Contains the parameter checks of the *main.nf* file. |
| **test folder** | Validation scripts of the database and of the strain identity branch (*smoke_test_strain_id.py*, *validate_markers.py*, *validate_panel.py*). |
| **galaxy folder** | Galaxy tool running the whole pipeline in a single Galaxy job, whose form recapitulates *nextflow.config*, and planemo scripts to run it in a local Galaxy (see [galaxy/README.md](galaxy/README.md)). |
| **.shed.yml** | Description of the Tool Shed repository of the Galaxy tool. |
| **web folder** | Public web page (Shiny form generated from *nextflow.config*) and submission daemon to Galaxy, with its Dockerfile, built by *.gitlab-ci.yml* (see [web/README.md](web/README.md)). |
| **logo.png** | Logo of the project. |
| **Licence.txt** | Licence of the release. |


<br><br>
## HOW IT WORKS

No assembly is made. The sample sequences are first compared with a strain identity database (see [here](https://github.com/orgs/lab-strain-watch/repositories)), through three categories:

| Evidence | How | Reference genome used at run time? |
| :--- | :--- | :--- |
| **Strain SNP and indel panel** | Alignment-free. For each panel site, the database holds the k-mers of each allele. The sample sequences are streamed once and the k-mers of each allele are counted, giving one allele per site where the depth is at least `kmer_min_depth`. | No |
| **Sex** | Alignment-free. Same as SNP and indels. | No |
| **Native marker genes** (URA3, HIS3, LEU2, TRP1, LYS2, ADE2, MET15, CAN1, ARG4, LEU5) | The reads are aligned with bowtie2 against a whole reference genome (`panel_reference` parameter of the *nextflow.config file*), so that reads from paralogs and repeats get a low mapping quality. The allele of each marker is read from the alignment (substitutions and small indels from the pileup, deletions from the depth, insertions from soft-clipped reads) and compared with the alleles catalogued in the database, judged by the protein it encodes. | Yes |


The panel strains tested are ranked by likelihood. If the best strain matches less than `min_concordance` (parameter of the *nextflow.config file*) of its expected alleles, the sample is reported as unassigned rather than forced onto the closest strain, which is what flags a strain absent from the database. 

Sex (if relevant, e.g., BY4741 / BY4742) and Native marker genes status are then compared to the most probable strain computed.

What is also checked:

| Evidence | How | Reference genome used at run time? |
| **Mitochondrial genome, growth on glycerol** | From the same alignment, the depth over the exons of the mitochondrial genes needed for respiration (database *mito/mito_genes.tsv*), against the nuclear depth: rho+ (every gene covered, growth on glycerol expected), rho- (whole genes missing), rho0 (no mtDNA), uncertain (genes covered in part). Optional introns are not counted, point mutations (mit-) and nuclear respiratory mutations are not seen. In *marker_genotype.json* under `mitochondria`, in the Decision of the report (Glycerol) and in its own section. | Yes |
| **Lab tool elements** (plasmids, cassettes, drug resistance markers, ...) | The reads are aligned with bowtie2 against the tool sequences of the database. Each element is scored by breadth of coverage and identity. The elements present are then placed in the genome from the alignment against the whole reference genome made for the native marker genes: reads crossing the junction between the genome and a construct (soft-clipped where the construct begins, the clipped bases being found in the element), read pairs with one read in the element and its mate on the genome, and, for an element the reference itself carries (a yeast promoter or terminator), its native locus. | Only the tool sequences for the scores; the whole reference genome for the positions |


<br><br>
## INPUT

| Parameter | Description |
| :--- | :--- |
| `fastq_path` | The sample sequences, despite the parameter name not only FASTQ. Either a run accession of SRA (SRR...), ENA (ERR...) or DDBJ (DRR...), downloaded with fastq-dump, or a file, or a folder of files. All the files of a folder are analyzed together as a single sample. |
| `seq_kind` | `"reads"` (default) for sequencing reads, short or long, or `"assembly"` for contigs, scaffolds or a complete genome sequence. |
| `seq_molecule` | `"DNA"` (default) or `"RNA"` (RNA-seq or direct RNA reads). The analysis is the same, but the report warns that depths then follow gene expression, and that a strain site or a marker gene that is not expressed is reported as not covered. |
| `cram_reference` | Only for a CRAM sample: the FASTA reference it was compressed against. If "NULL", samtools finds the reference from its checksum, which needs network access. |
| `db_path`, `organism_dir` | Strain identity database used: *\<db_path\>/\<organism_dir\>* (*../yeast_db/db/yeast* by default, i.e., yeast_db cloned next to lab_strain_watch). |
| `panel_reference` | Reference genome of the database (S288C R64-1-1 for yeast), used to genotype the native marker genes. |

See the *nextflow.config* file for the other parameters (`kmer_min_depth`, `marker_min_depth`, `strain_error_rate`, `min_concordance`, `kmer_max_bases`, `shred_max_length`, `shred_length`, ...), each one described there.

Accepted sample formats (recognized from the file content, except BAM, SAM and CRAM, recognized from the extension):

| Format | Extensions | Note |
| :--- | :--- | :--- |
| FASTQ | .fastq, .fq | Paired files are recognized by their _1/_2 or _R1/_R2 suffix (example: sample_R1.fastq.gz and sample_R2.fastq.gz). |
| FASTA | .fasta, .fa, .fna, .fas, .seq, .txt | Reads or assembly contigs. A uniform base quality of 40 is given. |
| GenBank | .gb, .gbk, .gbff | The sequence of the ORIGIN section is used. |
| EMBL | .embl, .dat | The sequence of the SQ section is used. |
| BAM, SAM, CRAM | .bam, .sam, .cram | Converted back to reads with samtools (pairs kept). The original alignment is not reused. |

Each of these files can be compressed (.gz, .bz2) or put into a .zip archive (except BAM, SAM and CRAM). U is read as T (RNA) and IUPAC ambiguity letters as N.

Long reads (Oxford Nanopore, PacBio) and assemblies: the k-mer scan reads each sequence whole. For the bowtie2 alignments, which are designed for short reads, the sequences longer than `shred_max_length` (1000 bp by default) are cut into pieces of `shred_length` bases (250 by default), without overlap for reads (the depth stays the sequencing depth), and with half overlapping pieces for an assembly. With `seq_kind = "assembly"`, the depth thresholds `kmer_min_depth` and `marker_min_depth` are set to 1, since an assembly has a single copy of each sequence, and the depths of the tool elements and of the native marker genes no longer measure copy number. The mating type of an assembly is read from the number of copies of the a and alpha cassettes (2 : 1 for MATa, 1 : 2 for MATalpha). A short-read assembly often collapses the three nearly identical cassettes (MAT, HML, HMR) into one a and one alpha contig: the mating type is then reported as unknown rather than as a MATa/alpha diploid. The *tsv/prepare_reads.tsv* output file describes what was done with each input file.

Not accepted: raw Nanopore signal (FAST5, POD5), which must be basecalled first (with dorado, for instance), and Illumina BCL run folders, which must be demultiplexed first (with bcl-convert, for instance).

<br><br>
## STRAIN IDENTITY DATABASE

The database is in the separate repository [yeast_db](https://github.com/lab-strain-watch/yeast_db), which also holds *build_db.nf*, the workflow that builds it. The yeast database is provided in its *db/yeast* folder and is only read by *main.nf*: nothing has to be built to analyze a sample. On galaxy.pasteur.fr, it is installed on the server by the Galaxy administrators (see [galaxy/README.md](galaxy/README.md)). It contains:

| Folder | Content |
| :--- | :--- |
| *organism.tsv* | Organism, reference genome and summary of the database. |
| *panels* | Strain SNP and indel panel: *panel_sites.tsv.gz* (sites and expected allele of each strain), *panel_kmers.tsv.gz* (k-mers of each allele), *panel_strains.tsv* (strains of the panel), *panel_report.txt* (build log), *panel_tree.nwk* (tree of the panel strains, opened from the report) and *panel_tree_report.txt* (its build log). |
| *strains* | Strain catalogue (*strains.tsv*) and strain name aliases (*aliases.tsv*). |
| *mating* | MATa and MATalpha k-mers (*mating_kmers.tsv*). |
| *markers* | Native marker gene loci (*marker_loci.tsv*, *.bed*, *.fasta*), allele catalogue (*marker_alleles.tsv*), documented genotypes of laboratory strains (*documented_genotypes.tsv*) and repeats masked in the flanks (*marker_repeats.tsv*). |
| *tools* | Laboratory tool element sequences (*tools.fasta*) and their description (*tools.tsv*). |

The way each part was made, the sources, the choice of k and the validation of the yeast panel are described in the *db/README.md* file of yeast_db, together with the commands to rebuild them. Warning: `panel_reference` must be the reference genome the database was built on.

<br><br>
## HOW TO RUN

### 1. Connection

- Open [https://galaxy.pasteur.fr](https://galaxy.pasteur.fr) in a web browser and log in with your Institut Pasteur account.
- Type *lab_strain_watch* in the search field of the *Tools* panel (left), then click on the tool: its form opens in the central panel.

<br><br>
### 2. Parameters

- Click on *Upload* (left) and add the sequence files of the sample (see [INPUT](#input)): *Choose local file*, or *Paste/Fetch data* for URLs, then *Start*. They appear in the history (right panel). Check their datatype (e.g., *fastqsanger.gz*), to change with the *Edit attributes* icon (pencil) if wrong.
- The form recapitulates the *nextflow.config* file: same sections, one field per parameter, with the name of the parameter in parentheses and its description below.
- *Sample sequences (fastq_path)*:
    - *Sequence file(s)*: one or several files of the sample, in any accepted format, the mates being paired by the _1/_2 or _R1/_R2 suffix of their names.
    - *Paired-end reads*: forward and reverse files, or a paired collection.
    - *Run accession*: SRR..., ERR... or DRR..., downloaded during the run. This requires an internet access of the Galaxy computing nodes: otherwise, import the reads into the history first (e.g., with the *Faster Download and Extract Reads in FASTQ* tool, if available) and use *Sequence file(s)*.
    - *Sample name*: optional, name of the sample in the reports (default: name of the first file).
- *Kind of sequences (seq_kind)*: reads or assembly. *Molecule (seq_molecule)*: DNA or RNA.
- *Sample handling* (optional): FASTA reference of a CRAM sample, and the cutting of the long sequences (long reads, contigs) into pieces before the bowtie2 alignments.
- *Strain identity* (optional): thresholds of the strain call. The strain identity database (yeast, 1,172 strains) is installed on the server.
- *Outputs*: optional additional datasets (see [Results](#4-results) below).
- *Parameter file* (optional): JSON file `{"parameter": "value"}` whose values replace those of the form, e.g., to reuse a set of parameters.
- The execution parameters are set by Galaxy.
- Optional: *Email notification* (bottom of the form, if available), to receive an e-mail at the end of the run.
- Click on *Run Tool*.

<br><br>
### 3. Run follow-up

- The datasets of the run appear at the top of the history: grey (queued), yellow (running), green (finished) or red (error).
- From a few minutes for an assembly to about an hour for a whole-genome sequencing run, plus the waiting time of the job on the Maestro cluster.
- Click on a dataset, then on *Dataset Details* (i icon): parameters of the run and nextflow messages (*Tool Standard Output* and *Tool Standard Error*).
- Red dataset: click on it to read the error message. *Error* icon (bug): report to the Galaxy administrators (or post an issue [here](https://github.com/lab-strain-watch/lab_strain_watch/issues)). *Run Job Again* icon: form with the same parameters, to rerun.
- Warnings of the run at the end of *report.html* (e.g., NO SNP PANEL or NO NATIVE MARKER DATABASE: database of the server incomplete, to report to the Galaxy administrators).

<br><br>
### 4. Results

- *strain identity report* dataset, *View* icon (eye): main result, *strain_identity_report.html* (decision, warnings, strain ranking, tool elements, native marker genes, mitochondrial genome, see [OUTPUT](#output)); its *tree* button, next to *Strains in the panel*, opens the tree of all the panel strains, the strain called highlighted.
- *report and result folder* dataset:
    - *View* icon (eye): *report.html* (database used, warnings) in the central panel, with working links to the tables and reports.
    - *Download* icon (floppy disk, after clicking on the dataset): zip file of the whole result folder (see [OUTPUT](#output)), saved where the web browser proposes.
- Optional datasets (*Outputs* section of the form): *tables* (.tsv files as a collection, for other Galaxy tools) and *result folder (zip)*.
- The datasets stay in the history: *Delete* icon (bin) when no longer needed.

<br><br>
## OUTPUT

By default, all the results are returned in a *results* folder where the *main.nf* executed file is located (created if does not exist). This can be changed using the `out_path_ini` parameter of the *nextflow.config* file. By default, each execution produces a new folder named *lab_strain_watch_\<ID\>*, created inside the *results* folder and containing all the outputs of the execution. The name of the folder can be changed using the `result_folder_name` parameter of the *nextflow.config* file. The new name file will be followed by an \<ID\> (number of seconds since 1970-01-01) in all cases.
<br><br>
Files of the strain identity branch (*strain_identity_report.html*, *site_calls.tsv*, *strain_ranking.tsv*, *strain_id.json*, *tool_coverage.tsv*, *tool_positions.tsv*, *marker_genotype.tsv*, *marker_genotype.json*) are absent if the strain panel is missing from the database. The *marker_genotype.\** and *tool_positions.tsv* files are also absent if the native marker database or the `panel_reference` file is missing.
<br><br>
<table>
    <tr>
        <th align="left">lab_strain_watch_&lt;ID&gt; folders and files</th>
        <th align="left">Description</th>
    </tr>
    <tr>
        <th align="left">strain_identity_report.html</th>
        <td>
            Main result. HTML report of the strain identity of the sample, made from the <i>tsv/strain_id.json</i> file, with a table of contents of its sections on the left side:
            <ul>
                <li><b>Header</b>: organism, reference genome, date, numbers of panel sites and strains. The <b>tree</b> button next to <i>Strains in the panel</i> opens, in a new page, the tree of all the panel strains (BIONJ over the panel sites, <i>panels/panel_tree.nwk</i> of the database), with the strain called highlighted and the 10 best ranked strains marked: rectangular or circular layout, zoom, part of the tree around a strain, search, details and links of each strain, Newick and SVG download. Its branch lengths are distances over the panel sites, not genome-wide divergences (see <i>db/README.md</i> of yeast_db). The button is absent, with a warning, if the database has no tree.</li>
                <li><b>Decision</b>: bullet points only: the strain called, then whether the mating type, each native marker and each foreign tool element present agree or disagree with the documented genotype of that strain; or an <b>UNASSIGNED</b> banner when no strain of the panel matches the sample well enough (with the reason). The strain is called in two steps, shown in the <b>Decided by</b> point. First the SNP panel: the first strain of the ranking is called when no other strain is close to it. Second, when other strains are too close to the first one for the SNPs alone (a log-likelihood within ln 1000 of the best, or fewer than 200 separating sites, i.e., called sites where the genotypes of the two strains are known and differ, whatever allele the sample carries there: close relatives such as BY4741 and BY4742, which differ by their mating type and by whole-gene deletions, met15D0 against lys2D0, rather than by SNPs), each of them is compared with the mating type and the native markers measured in the sample, from its documented genotype or, without one, from the marker alleles of its assembly. A strain contradicted at any point is ruled out; among the others, the one agreeing at the most points is called, the SNP ranking breaking the ties. Each of these strains is listed as called, ruled out (with what contradicts it) or not excluded, with its SNP rank, its PP, its log-likelihood below the best, its separating sites and the number of them where the sample carries the allele of the best strain (and its own): the strain called first, the others in a drop-down list under it. A strain close by its few separating sites can have a PP of 0 and rank below the 10 strains shown in <b>Strain ranking</b>, even when the sample carries the allele of the best strain at all of them: a lead built on so few sites is not trusted alone, as panel errors, a misassembled region or the drift of a laboratory stock can make it. When every one of them is contradicted, or none can be compared, the first of the SNP ranking is kept and a warning says so.</li>
                <li><b>Warnings</b>: every point that limits the call (low coverage, many heterozygous sites, native markers differing from those documented for the best strain, tool elements native to the host genome, ...).</li>
                <li><b>Strain ranking</b>: the top strains of the ranking (see <i>strain_ranking.tsv</i> below), with, for each, the detail of the panel sites (every mismatching and uncovered site, the matching sites being thinned).</li>
                <li><b>Tool elements</b>: score of each element of the tool database (see <i>tool_coverage.tsv</i> below), as green (score ≥ 99, present), orange (90 ≤ score &lt; 99, present but partial or diverged) or red (score &lt; 90, absent), with score = 100 x breadth x identity. Short elements are shown for information only. Elements that are native to the host genome (a yeast promoter, terminator or gene used in a cassette) are not counted as evidence of a construct. The <b>position in the genome</b> column gives, for each element present, where the reads place it on the reference genome (chromosome:first-last base, see <i>tool_positions.tsv</i> below): insertion (both junctions with the genome seen, the construct replacing the bases given, or between the two bases given), one junction, approximate (~, from read pairs only, within about a fragment length), native locus (the reference genome itself carries this sequence there), or not found (no read joins the element to the genome: a plasmid, or a construct whose ends are not in the database). Mouse over a position, or click the ? of the element, for the evidence.</li>
                <li><b>Native marker genes</b>: genotype of each marker (see <i>marker_genotype.tsv</i> below), the documented laboratory genotypes in a drop-down list (closed, the genotype agreeing at the most markers first, each with the documented allele of every marker and whether the sample agrees), and the limits of the analysis.</li>
            </ul>
        </td>
    </tr>
    <tr>
        <th align="left">report.html</th>
        <td>
            HTML report of the run: link to the parameters used (<i>nextflow.config</i>), summary of the strain identity database used (strains, native markers, tool elements) and the warnings of the run.
        </td>
    </tr>
    <tr>
        <th align="left">reports</th>
        <td>
            Folder containing:
            <ul>
                <li><b>*.log</b>: the log of each process of the pipeline (<i>prepare_reads_report.log</i>, <i>kmer_scan_report.log</i>, <i>bowtie2_tools_report.log</i>, <i>tool_depth_report.log</i>, <i>tool_score_report.log</i>, <i>marker_align_report.log</i>, <i>tool_position_report.log</i>, <i>marker_genotype_report.log</i>, <i>strain_score_report.log</i>, <i>fastq_dump_report.log</i> if an accession was downloaded, <i>aln_to_fastq_report.log</i> if a BAM, SAM or CRAM was given, <i>print_report.log</i>).</li>
                <li><b>nextflow.config</b>: copy of the <i>nextflow.config</i> file used, and <b>nextflow.config.html</b>, its display linked from <i>report.html</i>.</li>
                <li><b>Run_info.txt</b>: git repository and commit of the code used, command line, execution mode, loaded modules (<code>--modules</code> option), result path and nextflow version.</li>
                <li><b>warnings.txt</b>: all the warnings of the run.</li>
                <li><b>.nextflow.log</b> and <b>Log_info.txt</b>: copy of the nextflow log (incomplete, the complete one is where the run was launched, as indicated in <i>Log_info.txt</i>).</li>
                <li><b>report.html</b>, <b>timeline.html</b>, <b>trace.txt</b> and <b>nf_dag.png</b>: the nextflow execution report, timeline, trace of the processes and map of the processes execution. Not to be confused with the <i>report.html</i> file of the result folder.</li>
            </ul>
        </td>
    </tr>
    <tr>
        <th align="left">tsv</th>
        <td>
            Folder containing all the <b>tsv files</b> described below (and two <b>json</b> files). NA means not applicable.
        </td>
    </tr>
    <tr>
        <th align="left">- prepare_reads.tsv</th>
        <td>
            What was done with each input file of the sample, before the analyses. One row per file.
            <ul>
                <li>file: name of the input file (BAM, SAM and CRAM files appear as the FASTQ files they were converted into).</li>
                <li>format: format recognized from the file content (fastq, fasta, genbank or embl).</li>
                <li>role: R1 or R2 (first or second mate of a pair), U (unpaired sequences) or S (sequences cut into pieces for the bowtie2 alignments, see <code>shred_max_length</code>).</li>
                <li>action: "used as is" (FASTQ file used without any change), "converted" (rewritten as FASTQ: other format, U to T, IUPAC letters to N) or "cut into &lt;N&gt; pieces of &lt;L&gt; bp".</li>
                <li>sequences: number of sequences in the file (NA if used as is).</li>
                <li>bases: number of bases in the file (NA if used as is).</li>
                <li>longest: length of the longest sequence, in bp (NA if used as is).</li>
            </ul>
        </td>
    </tr>
    <tr>
        <th align="left">- site_calls.tsv</th>
        <td>
            Allele observed at each site of the strain panel (<i>db/&lt;organism&gt;/panels/panel_sites.tsv.gz</i>), from the alignment-free k-mer scan of the sample. One row per panel site. The counts are the number of occurrences of the k-mers of each allele (both strands).
            <ul>
                <li>site_id: identifier of the site in the panel (see <i>panel_sites.tsv.gz</i> for its position and alleles).</li>
                <li>ref_count: count of the reference allele k-mers.</li>
                <li>alt_count: count of the alternative allele k-mers.</li>
                <li>depth: ref_count + alt_count.</li>
                <li>alt_fraction: alt_count / depth (0 if the site is not called).</li>
                <li>call: allele called. R (reference, alt_fraction ≤ 0.2), A (alternative, alt_fraction ≥ 0.8), H (both alleles, heterozygous diploid, mixed culture or cross) or ? (not covered: depth below the <code>kmer_min_depth</code> parameter of the <i>nextflow.config</i> file). A site not covered is never considered as reference.</li>
            </ul>
        </td>
    </tr>
    <tr>
        <th align="left">- strain_ranking.tsv</th>
        <td>
            All the strains of the panel, ranked by likelihood of the observed site calls (<i>site_calls.tsv</i>) given the alleles expected for each strain. One row per strain.
            <ul>
                <li>rank: rank of the strain (1 is the most probable).</li>
                <li>strain_id: strain identifier in the database (see <i>db/&lt;organism&gt;/panels/panel_strains.tsv</i>).</li>
                <li>display_name: strain name displayed in the report.</li>
                <li>clade: clade or family of the strain (e.g., "S288C family").</li>
                <li>source_type: source of the expected alleles of the strain in the database: matrix (1011 Yeast Genomes population matrix, Peter et al. 2018), assembly (genome assembly) or reads (sequencing runs called by <i>build_db.nf</i>).</li>
                <li>probability: posterior probability of the strain among the strains of the panel (uniform prior, per site error rate set by the <code>strain_error_rate</code> parameter of the <i>nextflow.config</i> file). Warning: the panel is assumed to be exhaustive, so a sample of a strain absent from the panel still gives a high probability to its closest strain. Use the concordance to judge the fit.</li>
                <li>concordance: sites_matching / sites_called. If below the <code>min_concordance</code> parameter of the <i>nextflow.config</i> file for the first strain, the sample is reported as unassigned.</li>
                <li>sites_called: number of panel sites called in the sample and with a known allele in the strain (sites_matching + sites_mismatching). A heterozygous site (H) matches a strain heterozygous at this site, and mismatches the others.</li>
                <li>sites_matching: number of called sites carrying the allele expected for the strain.</li>
                <li>sites_mismatching: number of called sites carrying another allele than the one expected for the strain (including a homozygous call at a site where the strain is heterozygous).</li>
                <li>sites_uncovered: number of panel sites not usable for the strain: site not called in the sample (?), or allele unknown in the strain. A high fraction of uncovered sites with an adequate depth is the sign of a strain absent from the panel.</li>
            </ul>
        </td>
    </tr>
    <tr>
        <th align="left">- strain_id.json</th>
        <td>
            All the results of the strain identity branch, from which <i>strain_identity_report.html</i> is made. Main keys:
            <ul>
                <li>sample, generated (date), organism, reference (database folder), panel (number of sites and strains).</li>
                <li>call: strain called. resolution: method (snp: the first of the SNP ranking; snp+markers: chosen by mating type and native markers among the strains too close to the best for the SNPs alone), tie_log_ratio and min_separating_sites (the two thresholds that make a strain close), summary, and candidates: the best strain and each strain too close to it (SNP rank, probability, close_by: the rules that make it close, likelihood and/or separating sites, separating sites, separating_with_best and separating_with_strain: the number of them where the sample carries the allele of the best strain and of this strain, mismatches, log-likelihood gap, documented mating type, source of its expected markers, the points that agree and those that contradict it, outcome: called, ruled out or not excluded). mating_check: the measured mating type against the one documented for the strain called (agree, disagree or not compared). tool_check: the foreign tool elements present, and those inserted by the documented genotype (after ::), each agreeing or not.</li>
                <li>unassigned (true / false) and unassigned_reason.</li>
                <li>mating_type: mean count of the MATa and MATalpha k-mers, a_fraction (MATa / (MATa + MATalpha)) and call: MATa (a_fraction ≥ 0.58), MATalpha (a_fraction ≤ 0.42), MATa/alpha (in between) or unknown (too few k-mers seen).</li>
                <li>coverage: median depth of the called sites and fraction of uncovered sites.</li>
                <li>strains: the top strains of the ranking, with the columns of <i>strain_ranking.tsv</i> plus source (the sequencing runs of a reads strain, the code of a matrix isolate, or the collection of an assembly strain: ScRAPdb or SGD, never the path of the assembly file), log_likelihood, sites_het (number of H sites), the site detail, and the number of sites separating each strain from the next one.</li>
                <li>tools and tool_summary: score of each tool element (columns of <i>tool_coverage.tsv</i> plus name, class, score, colour, present and origin: native or foreign to the host genome, and, for the elements present, positions and position_note: their positions in the genome as in <i>tool_positions.tsv</i>) and number of green, orange and red elements. tool_positions: whether the positions were searched, the read layout and length, notes, and the insertion sites with their junctions.</li>
                <li>markers: content of <i>marker_genotype.json</i>. markers_expected: comparison of the native markers with those expected for the first strain of the ranking.</li>
                <li>warnings: all the warnings shown in the report.</li>
            </ul>
        </td>
    </tr>
    <tr>
        <th align="left">- tool_coverage.tsv</th>
        <td>
            Coverage of each laboratory tool element of the database (<i>db/&lt;organism&gt;/tools/tools.fasta</i>, described in <i>tools.tsv</i>) by the reads (bowtie2 alignment against the tool sequences only). One row per element.
            <ul>
                <li>element_id: element identifier, &lt;NAME&gt;__&lt;HASH&gt; (the hash being the first characters of the sha256 of the sequence).</li>
                <li>length_bp: length of the element, in bp.</li>
                <li>covered_bases: number of bases of the element covered by at least one read.</li>
                <li>breadth: covered_bases / length_bp. Is the whole element there? A partial breadth is the signature of a truncated cassette, a shared promoter or terminator matching only its conserved part, or a chance match.</li>
                <li>identity: fraction of the aligned bases identical to the element (from the NM tag of the alignments). Is it the same sequence? A diverged homologue (the yeast's own TRP1 on the TRP1 of a cassette) can cover an element end to end with a low identity.</li>
                <li>mean_depth: mean number of reads per base of the element.</li>
                <li>reads: number of reads aligned on the element.</li>
            </ul>
            The score shown in <i>strain_identity_report.html</i> is 100 x breadth x identity.
        </td>
    </tr>
    <tr>
        <th align="left">- tool_positions.tsv</th>
        <td>
            Positions in the genome of the tool elements present (score &ge; 99), from the alignment of the reads against the whole reference genome (<code>panel_reference</code>) made for the native marker genes. A read crossing the junction between the genome and a construct aligns on its genomic part and is soft-clipped where the construct begins: the clipped bases (15 or more, mapping quality of the genomic part 20 or more) are looked up in the elements (15-mers), which gives the junction to the base, give or take the few bases the genome and the element may share there; it needs reads of about 50 bases or more. With paired reads, a read in an element whose mate is on the genome places the junction within about a fragment length, which is what places a construct with shorter reads. For long sequences cut into pieces (an assembly, long reads), the pieces of the same sequence near a junction are also used, which tells apart cassettes that share their ends (kanMX, natMX6, hphMX6 and HIS3MX6 all start with the same TEF promoter and end with the same TEF terminator). Junctions within 20 kb of each other on one chromosome make one insertion site. One row per position of each element present, or one row "not found".
            <ul>
                <li>element_id, name, origin: as in <i>tools.tsv</i> (origin native: the yeast genome carries this sequence too).</li>
                <li>position_kind: insertion (both junctions seen to the base), one junction (one side seen), approximate (a side seen by read pairs only), native (native locus: the reference genome carries the sequence there, so an unmodified strain has it there too), or not found.</li>
                <li>chromosome, reference_sequence: chromosome name (chrI to chrXVI, chrM) and name of the sequence in the reference FASTA.</li>
                <li>start, end: first and last base (1-based) of the reference replaced by the construct; for an insertion between two bases, these two bases; for one junction, the last base before the construct or the first base after it.</li>
                <li>strand: orientation of the element on the reference (+ or -).</li>
                <li>junction_reads, read_pairs: reads crossing the junctions of the site, read pairs with one read in the construct and its mate on the genome.</li>
                <li>position: chromosome:start-end, ~ when approximate.</li>
                <li>or_instead: other elements that may be at this site in place of this one: the reads seen there are in sequence they share (two cassettes built on the same TEF promoter and terminator, with reads too short to reach their resistance gene).</li>
                <li>detail: the same in words, with the element that holds this one when it is part of a larger construct (the TEF promoter of a kanMX cassette), or the element a native one is joined to (the ADH1 terminator between a GFP tag and its marker).</li>
            </ul>
        </td>
    </tr>
    <tr>
        <th align="left">- marker_genotype.tsv</th>
        <td>
            Genotype of each native marker gene (URA3, HIS3, LEU2, TRP1, LYS2, ADE2, MET15, CAN1, ARG4, LEU5), from the reads aligned on the whole reference genome (<code>panel_reference</code>). One row per gene. For a heterozygous gene, only the first allele is shown here (both are in <i>marker_genotype.json</i>).
            <ul>
                <li>gene: marker gene.</li>
                <li>status: active (functional protein), inactive (loss of function), uncertain (protein form not recognized, never called active or inactive on a guess), mixed (two alleles of different status, heterozygous) or no_data (the marker loci are covered below the <code>marker_min_depth</code> parameter of the <i>nextflow.config</i> file).</li>
                <li>effect: effect of the allele on the protein. reference or synonymous (active); absent, orf_deleted, insertion, large_deletion, start_lost, frameshift or nonsense (loss of function); missense, inframe_indel, late_truncation, stop_lost or partial (part of the ORF without data) (uncertain unless the protein form is in the allele catalogue).</li>
                <li>protein_change: description of the change (e.g., "p.F104fs, protein of 136 instead of 220 aa", "804 of 804 bp of the ORF deleted"). Empty for the reference protein.</li>
                <li>known_as: documented laboratory allele(s) encoding the same protein, with the strain documented, separated by ; (e.g., "ura3D0 (BY4742)"). Empty if none.</li>
                <li>reason: why this status was given (e.g., "loss of function from the sequence itself (frameshift)", "protein form of the documented allele ...", "protein form carried by documented wild-type ...").</li>
                <li>structural: structural evidence, separated by " | ": regions covered by no read (deletion, genomic coordinates) and foreign sequence joined to the locus (insertion, from soft-clipped reads).</li>
            </ul>
        </td>
    </tr>
    <tr>
        <th align="left">- marker_genotype.json</th>
        <td>
            Complete native marker results: baseline_depth and read_length of the sample, limits of the analysis, genotype (e.g., "ura3D0 his3D1 leu2D0 TRP1 lys2D0 ADE2 MET15 CAN1"), markers (per gene: locus, depths, status, all the alleles with their coding, upstream and downstream events, protein form, documented names and carriers in the catalogue, structural evidence, heterozygous calls), documented_genotype_matches (documented laboratory strains with the number of genes agreeing and the differences) and catalogue_assemblies (number of assemblies in the allele catalogue).
        </td>
    </tr>
</table>


<br><br>
## VERSIONS


The different releases are [here](https://github.com/lab-strain-watch/lab_strain_watch/releases).

The different tagged versions are [here](https://github.com/lab-strain-watch/lab_strain_watch/tags).


<br><br>
## LICENCE


This package of scripts can be redistributed and/or modified under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
Distributed in the hope that it will be useful, but without any warranty; without even the implied warranty of merchantability or fitness for a particular purpose.
See the GNU General Public License for more details at https://www.gnu.org/licenses or in the Licence.txt attached file.


<br><br>
## CITATION

Not yet published.

The yeast database uses the population matrix of the 1011 Yeast Genomes project: Peter J. et al. Genome evolution across 1,011 *Saccharomyces cerevisiae* isolates. Nature 556, 339–344 (2018).


<br><br>
## CREDITS

[Gael A. Millot](https://gitlab.pasteur.fr/gmillot), Institut Pasteur, Université Paris Cité, Bioinformatics and Biostatistics Hub, 75015 Paris, France

<br><br>
## ACKNOWLEDGEMENTS


The developers & maintainers of the mentioned softwares and packages, including:

- [Nextflow](https://www.nextflow.io/)
- [Apptainer](https://apptainer.org/)
- [Docker](https://www.docker.com/)
- [Bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/)
- [Samtools](https://www.htslib.org/)
- [SRA Toolkit](https://github.com/ncbi/sra-tools)
- [Python](https://www.python.org/)
- [R](https://www.r-project.org/)
- [The 1011 Yeast Genomes project](http://1002genomes.u-strasbg.fr/)
- [Saccharomyces Genome Database](https://www.yeastgenome.org/)
- [Github](https://github.com/)
- [Bash](https://www.gnu.org/software/bash/)
- [Ubuntu](https://ubuntu.com/)


<br><br>
## WHAT'S NEW IN

#### v1.0

Everything
