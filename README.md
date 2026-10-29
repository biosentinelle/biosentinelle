<img src="logo.png" align="center" width="300" alt="Biosentinelle Logo">

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

Check the identity of a sequenced laboratory strain, i.e., answer "is this tube really the strain written on its label?". From the sequences of a single sample (reads or assembly), biosentinelle reports:

- the most probable strain among the strains of a database (1,172 *Saccharomyces cerevisiae* strains for the yeast database), or "unassigned" if none fits well enough,
- the mating type (MATa, MATalpha or MATa/alpha),
- the genotype of the native auxotrophic and counter-selection marker genes (URA3, HIS3, LEU2, TRP1, LYS2, ADE2, MET15, CAN1),
- the laboratory tool elements present (plasmids, cassettes, drug resistance markers, tags, promoters, ...),
- coverage and variant counts against user-provided genomes, features and mating type sequences.

<br><br>
## WARNINGS

- Currently dedicated to *Saccharomyces cerevisiae*: only the yeast database (`db/yeast`) is provided. Another organism needs its own database, built with *build_db.nf* (see [STRAIN IDENTITY DATABASE](#strain-identity-database)).
- The strain probability is a posterior over the strains of the database, which are assumed to be exhaustive. It is not the probability that the sample *is* that strain. The concordance (fraction of the expected alleles found) is the measure of fit. See the [OUTPUT](#output) section.
- A single sample per run. All the files of the `fastq_path` folder are analyzed together as one sample.
- Laboratory stocks drift: a few mismatching sites between a sample and its documented strain are expected (see *db/README.md*).

<br><br>
## CONTENT

| Files and folder | Description |
| :--- | :--- |
| **main.nf** | File that can be executed using a linux terminal, a MacOS terminal or Windows 10 WSL2. Analysis of a sample. |
| **nextflow.config** | Parameter settings for the *main.nf* and *build_db.nf* files. Users have to open this file, set the desired settings and save these modifications before execution. Of note, this configuration file is systematically saved in the *reports* folder (see [below](#output)) during each execution, to save the parameter settings. |
| **build_db.nf** | Builds the strain panel of the strain identity database (*db/\<organism\>/panels*), once per organism. Not required to analyze a sample with the provided yeast database. See [STRAIN IDENTITY DATABASE](#strain-identity-database). |
| **db folder** | Strain identity database, one subfolder per organism (*db/yeast*). See *db/README.md* for its full description. |
| **bin folder** | Contains the scripts (Python, R, Rmarkdown) required by the *main.nf* and *build_db.nf* files. |
| **modules folder** | Contains the Nextflow processes required by the *main.nf* and *build_db.nf* files. |
| **conf folder** | Contains the parameter checks of the *main.nf* file. |
| **test folder** | Validation scripts of the database and of the strain identity branch (*smoke_test_strain_id.py*, *validate_markers.py*, *validate_panel.py*). |
| **logo.png** | Logo of the project. |
| **Licence.txt** | Licence of the release. |


<br><br>
## HOW IT WORKS

No assembly is made. The sample sequences are compared with a strain identity database (*db/\<organism\>*, built once by *build_db.nf* from a reference genome and the sequences of 1,172 *S. cerevisiae* strains), through three independent kinds of evidence:

| Evidence | How | Reference genome used at run time? |
| :--- | :--- | :--- |
| **Strain SNP and indel panel** | Alignment-free. For each panel site, the database holds the k-mers of each allele. The sample sequences are streamed once and the k-mers of each allele are counted, giving one allele per site where the depth is at least `kmer_min_depth`. The same is done for the mating type (MATa / MATalpha) k-mers. | No |
| **Native marker genes** (URA3, HIS3, LEU2, TRP1, LYS2, ADE2, MET15, CAN1) | The reads are aligned with bowtie2 against the whole reference genome (`panel_reference`), so that reads from paralogs and repeats get a low mapping quality. The allele of each marker is read from the alignment (substitutions and small indels from the pileup, deletions from the depth, insertions from soft-clipped reads) and compared with the alleles catalogued in the database, judged by the protein it encodes. | Yes |
| **Lab tool elements** (plasmids, cassettes, drug resistance markers, ...) | The reads are aligned with bowtie2 against the tool sequences of the database. Each element is scored by breadth of coverage and identity. | No, only the tool sequences |

The panel strains are then ranked by likelihood. If the best strain matches less than `min_concordance` of its expected alleles, the sample is reported as unassigned rather than forced onto the closest strain, which is what flags a strain absent from the database. The result is the *strain_identity_report.html* file.

In addition, the historical part of the pipeline aligns the reads (bowtie2) against the genomes of `ref_path`, the features of `features_path` and the mating type sequences of `sex_path`, for coverage plots, read counts after mapping quality 20 filtering, and variant counts (bcftools), summarized in *report.html*. This part does not use the database.

<br><br>
## INPUT

| Parameter | Description |
| :--- | :--- |
| `fastq_path` | The sample sequences, despite the parameter name not only FASTQ. Either a run accession of SRA (SRR...), ENA (ERR...) or DDBJ (DRR...), downloaded with fastq-dump, or a file, or a folder of files. All the files of a folder are analyzed together as a single sample. |
| `seq_kind` | `"reads"` (default) for sequencing reads, short or long, or `"assembly"` for contigs, scaffolds or a complete genome sequence. |
| `seq_molecule` | `"DNA"` (default) or `"RNA"` (RNA-seq or direct RNA reads). The analysis is the same, but the report warns that depths then follow gene expression, and that a strain site or a marker gene that is not expressed is reported as not covered. |
| `cram_reference` | Only for a CRAM sample: the FASTA reference it was compressed against. If "NULL", samtools finds the reference from its checksum, which needs network access. |
| `ref_path` | Genome(s) used for the variant counts (*variant_counts.tsv*). A file or a folder of files, or a .zip archive. A single multi-sequence FASTA file is split into one file per sequence. |
| `features_path` | Feature sequences (cassettes, genes, ...) used for the coverage plots and the *coverage_feat.tsv* file. Same format as `ref_path`. |
| `sex_path` | Mating type sequences (MATa, MATalpha) used for the coverage plots and the *coverage_sex.tsv* file. Same format as `ref_path`. |
| `db_path`, `organism_dir` | Strain identity database used: *\<db_path\>/\<organism_dir\>* (*db/yeast* by default). |
| `panel_reference` | Reference genome of the database (S288C R64-1-1 for yeast), used to genotype the native marker genes. |

See the *nextflow.config* file for the other parameters (`strain_id_on`, `kmer_min_depth`, `marker_min_depth`, `strain_error_rate`, `min_concordance`, `kmer_max_bases`, `shred_max_length`, `shred_length`, `ylab`, ...), each one described there.

Accepted sample formats (recognized from the file content, except BAM, SAM and CRAM, recognized from the extension):

| Format | Extensions | Note |
| :--- | :--- | :--- |
| FASTQ | .fastq, .fq | Paired files are recognized by their _1/_2 or _R1/_R2 suffix (example: sample_R1.fastq.gz and sample_R2.fastq.gz). |
| FASTA | .fasta, .fa, .fna, .fas, .seq, .txt | Reads or assembly contigs. A uniform base quality of 40 is given. |
| GenBank | .gb, .gbk, .gbff | The sequence of the ORIGIN section is used. |
| EMBL | .embl, .dat | The sequence of the SQ section is used. |
| BAM, SAM, CRAM | .bam, .sam, .cram | Converted back to reads with samtools (pairs kept). The original alignment is not reused. |

Each of these files can be compressed (.gz, .bz2) or put into a .zip archive (except BAM, SAM and CRAM). U is read as T (RNA) and IUPAC ambiguity letters as N.

Long reads (Oxford Nanopore, PacBio) and assemblies: the k-mer scan reads each sequence whole. For the bowtie2 alignments, which are designed for short reads, the sequences longer than `shred_max_length` (1000 bp by default) are cut into pieces of `shred_length` bases (250 by default), without overlap for reads (the depth stays the sequencing depth), and with half overlapping pieces for an assembly. With `seq_kind = "assembly"`, the depth thresholds `kmer_min_depth` and `marker_min_depth` are set to 1, since an assembly has a single copy of each sequence, and the coverage plots no longer measure copy number. The *tsv/prepare_reads.tsv* output file describes what was done with each input file.

Not accepted: raw Nanopore signal (FAST5, POD5), which must be basecalled first (with dorado, for instance), and Illumina BCL run folders, which must be demultiplexed first (with bcl-convert, for instance).

The reference, feature and mating type sequences (`ref_path`, `features_path`, `sex_path`) can be FASTA, GenBank or EMBL files, compressed with gzip or not, or a .zip archive of them.

<br><br>
## STRAIN IDENTITY DATABASE

The yeast database is provided in *db/yeast* and is only read by *main.nf*: nothing has to be built to analyze a sample. It contains:

| Folder | Content |
| :--- | :--- |
| *organism.tsv* | Organism, reference genome and summary of the database. |
| *panels* | Strain SNP and indel panel: *panel_sites.tsv.gz* (sites and expected allele of each strain), *panel_kmers.tsv.gz* (k-mers of each allele), *panel_strains.tsv* (strains of the panel), *panel_report.txt* (build log). |
| *strains* | Strain catalogue (*strains.tsv*) and strain name aliases (*aliases.tsv*). |
| *mating* | MATa and MATalpha k-mers (*mating_kmers.tsv*). |
| *markers* | Native marker gene loci (*marker_loci.tsv*, *.bed*, *.fasta*), allele catalogue (*marker_alleles.tsv*), documented genotypes of laboratory strains (*documented_genotypes.tsv*) and repeats masked in the flanks (*marker_repeats.tsv*). |
| *tools* | Laboratory tool element sequences (*tools.fasta*) and their description (*tools.tsv*). |

The way each part was made, the sources, the choice of k and the validation of the yeast panel are described in [db/README.md](db/README.md), together with the commands to rebuild them. The strain panel is rebuilt with:

<pre>
nextflow run build_db.nf -c nextflow.config
</pre>

using the "Database construction" parameters of the *nextflow.config* file (`strain_table`, `population_matrix`, `reads_dir`, `reads_min_depth`, `panel_reference`, `panel_kmer_size`, `panel_options`). Warning: changing `panel_reference` or `panel_kmer_size` invalidates the whole database.

<br><br>
## HOW TO RUN

### 1. Prerequisite

Installation of:<br>
[nextflow DSL2](https://gael-millot.github.io/protocols/docs/Protocol%20152-rev0%20DSL2.html#_Toc208504071). Please, use the version indicated above.<br>
[Graphviz](https://www.graphviz.org/download/), `sudo apt install graphviz` for Linux ubuntu.<br>
[Apptainer](https://gael-millot.github.io/protocols/docs/Protocol%20135-rev0%20APPTAINER.html#_Toc160091693).<br>
<br>

Set the `apptainer_local_images` (local run) or `apptainer_cluster_images` (cluster run) parameter of the *nextflow.config* file to a folder where the apptainer images are stored. Images absent from this folder are pulled from Docker Hub at the first run.

<br>

### 2. Local running (personal computer)

#### 2.1. *main.nf* file in the personal computer

- Mount a server if required:

<pre>
DRIVE="Z" # change the letter to fit the correct drive
sudo mkdir /mnt/share
sudo mount -t drvfs $DRIVE: /mnt/share
</pre>

Warning: if no mounting, it is possible that nextflow does nothing, or displays a message like:
<pre>
Launching `main.nf` [loving_morse] - revision: d5aabe528b
/mnt/share/Users
</pre>

- Run the following command from where the *main.nf* and *nextflow.config* files are (example: \\wsl$\Ubuntu-20.04\home\gael):

<pre>
nextflow run main.nf -c nextflow.config # or nextflow run main.nf
</pre>

with -c to specify the name of the config file used.

<br><br>
#### 2.2. *main.nf* file in the public github repository

Run the following command from where you want the results:

<pre>
nextflow run biosentinelle/biosentinelle -r main
</pre>

<br><br>

### 3. Distant running (example with the Pasteur cluster)

#### 3.1. Pre-execution

Set `system_exec = 'slurm'` (and `slurm_queue`, `slurm_qos` if required) in the *nextflow.config* file.

Go into the directory where the main.nf and nextflow.config files are.
Copy-paste this code:

<pre>
EXEC_PATH=$(pwd) # where the bin folder of biosentinelle is located (by default, the same path as for the main.nf file)
export CONF_BEFORE=/opt/gensoft/exe # on maestro

export JAVA_CONF=java/13.0.2
export JAVA_CONF_AFTER=bin/java # on maestro
export APP_CONF=apptainer/1.3.5
export APP_CONF_AFTER=bin/apptainer # on maestro
export GIT_CONF=git/2.39.1
export GIT_CONF_AFTER=bin/git # on maestro
export GRAPHVIZ_CONF=graphviz/2.42.3
export GRAPHVIZ_CONF_AFTER=bin/graphviz # on maestro
export GRAALVM_CONF=graalvm/ce-java23-23.0.1 # required for nextflow
export GRAALVM_CONF_AFTER=bin/graalvm # on maestro
export NEXTFLOW_CONF=nextflow/24.10.3
export NEXTFLOW_CONF_AFTER=bin/nextflow # on maestro

MODULES="${CONF_BEFORE}/${JAVA_CONF}/${JAVA_CONF_AFTER},${CONF_BEFORE}/${APP_CONF}/${APP_CONF_AFTER},${CONF_BEFORE}/${GIT_CONF}/${GIT_CONF_AFTER},${CONF_BEFORE}/${GRAPHVIZ_CONF}/${GRAPHVIZ_CONF_AFTER},${CONF_BEFORE}/${GRAALVM_CONF}/${GRAALVM_CONF_AFTER},${CONF_BEFORE}/${NEXTFLOW_CONF}/${NEXTFLOW_CONF_AFTER}"
cd ${EXEC_PATH}
chmod -R 755 ${EXEC_PATH}/bin/*.* # nextflow needs x authorization for all files and subfolder files
module load ${JAVA_CONF} ${APP_CONF} ${GIT_CONF} ${GRAPHVIZ_CONF} ${GRAALVM_CONF}
module load ${NEXTFLOW_CONF}
</pre>

<br><br>

#### 3.2. *main.nf* file in a cluster folder

Run from where the *main.nf* and *nextflow.config* files are (which has been set thanks to the EXEC_PATH variable above):

<pre>
HOME_INI=$HOME
HOME="${HELIXHOME}/biosentinelle/" # $HOME changed to allow the creation of .nextflow into /$HELIXHOME/biosentinelle/, for instance. See NFX_HOME in the nextflow software script
nextflow run main.nf -c nextflow.config -resume --modules ${MODULES} # --modules ${MODULES} in order to have all the used module versions recorded into the report file, -resume to continue a potential job that has already been submitted
HOME=$HOME_INI
</pre>

<br><br>

#### 3.3. *main.nf* file in the public github repository

Modify the first line of the code below, in order to call the desired version of biosentinelle (branch or tag, see [here](https://github.com/biosentinelle/biosentinelle/tags)), and run (results will be where the EXEC_PATH variable has been set above):

<pre>
VERSION="main"
HOME_INI=$HOME
HOME="${HELIXHOME}/biosentinelle/" # $HOME changed to allow the creation of .nextflow into /$HELIXHOME/biosentinelle/, for instance. See NFX_HOME in the nextflow software script
nextflow run biosentinelle/biosentinelle -r $VERSION -c $HOME/nextflow.config --modules ${MODULES}
HOME=$HOME_INI
</pre>

<br><br>
### 4. Error messages and solutions

#### Message 1
```
Unknown error accessing project `biosentinelle/biosentinelle` -- Repository may be corrupted: /pasteur/sonic/homes/gmillot/.nextflow/assets/biosentinelle/biosentinelle
```

Purge using:
<pre>
rm -rf /pasteur/sonic/homes/gmillot/.nextflow/assets/biosentinelle*
</pre>

#### Message 2

```
permission denied
```

Use chmod to change the user rights. Example linked to files in the bin folder: 
```
chmod 755 bin/*.*
```

#### Message 3

```
WARNING:
NO SNP PANEL IN .../db/yeast/panels.
STRAIN RANKING IS SKIPPED. RUN build_db.nf FIRST.
```

The *panel_sites.tsv.gz*, *panel_kmers.tsv.gz* or *panel_strains.tsv* file is missing in *\<db_path\>/\<organism_dir\>/panels*. Check the `db_path` and `organism_dir` parameters of the *nextflow.config* file, or build the panel (see [STRAIN IDENTITY DATABASE](#strain-identity-database)). No *strain_identity_report.html* file is produced in that case.

#### Message 4

```
WARNING:
NO NATIVE MARKER DATABASE IN .../db/yeast/markers, OR panel_reference NOT FOUND.
NATIVE MARKER GENOTYPING IS SKIPPED.
```

Check that the `panel_reference` parameter of the *nextflow.config* file points to the reference genome of the database (S288C R64-1-1 for yeast).

#### Message 5

```
ERROR ~ Error executing process > 'Print_report (1)'
...
Command error:
  Quitting from lines ... (biosentinelle_template.rmd)
```

If obtained using `nextflow run main.nf -resume`, then rerun the same command once. Sometimes, nextflow shows difficulties to build the report.html file with `-resume`. If the problem persists, try `nextflow run main.nf`. If the problem still occurs, please send an issue [here](https://github.com/biosentinelle/biosentinelle/issues).

#### Error using -resume

Do not hesitate to rerun once the same command to reobtain the same error message. Indeed, sometimes nextflow has difficulties to reassemble all the cached files and reruning solves the problem.


<br><br>
## OUTPUT

By default, all the results are returned in a *results* folder where the *main.nf* executed file is located (created if does not exist). This can be changed using the `out_path_ini` parameter of the *nextflow.config* file. By default, each execution produces a new folder named *biosentinelle_\<ID\>*, created inside the *results* folder and containing all the outputs of the execution. The name of the folder can be changed using the `result_folder_name` parameter of the *nextflow.config* file. The new name file will be followed by an \<ID\> (number of seconds since 1970-01-01) in all cases.
<br><br>
Files of the strain identity branch (*strain_identity_report.html*, *site_calls.tsv*, *strain_ranking.tsv*, *strain_id.json*, *tool_coverage.tsv*, *marker_genotype.tsv*, *marker_genotype.json*) are absent if the `strain_id_on` parameter of the *nextflow.config* file is "FALSE", or if the strain panel is missing from the database. The *marker_genotype.\** files are also absent if the native marker database or the `panel_reference` file is missing.
<br><br>
<table>
    <tr>
        <th align="left">biosentinelle_&lt;ID&gt; folders and files</th>
        <th align="left">Description</th>
    </tr>
    <tr>
        <th align="left">strain_identity_report.html</th>
        <td>
            Main result. HTML report of the strain identity of the sample, made from the <i>tsv/strain_id.json</i> file:
            <ul>
                <li><b>Call</b>: the strain called, the mating type and the coverage, or an <b>UNASSIGNED</b> banner when no strain of the panel matches the sample well enough (with the reason). When several strains cannot be separated by the panel sites (close relatives, a parent and its derivatives), the mating type and the native markers are used to choose between them, and the strains left tied are listed.</li>
                <li><b>Warnings</b>: every point that limits the call (low coverage, many heterozygous sites, native markers differing from those documented for the best strain, tool elements native to the host genome, ...).</li>
                <li><b>Strain ranking</b>: the top strains of the ranking (see <i>strain_ranking.tsv</i> below), with, for each, the detail of the panel sites (every mismatching and uncovered site, the matching sites being thinned).</li>
                <li><b>Tool elements</b>: score of each element of the tool database (see <i>tool_coverage.tsv</i> below), as green (score ≥ 99, present), orange (90 ≤ score &lt; 99, present but partial or diverged) or red (score &lt; 90, absent), with score = 100 x breadth x identity. Short elements are shown for information only. Elements that are native to the host genome (a yeast promoter, terminator or gene used in a cassette) are not counted as evidence of a construct.</li>
                <li><b>Native marker genes</b>: genotype of each marker (see <i>marker_genotype.tsv</i> below), the documented laboratory genotypes it agrees with, and the limits of the analysis.</li>
            </ul>
        </td>
    </tr>
    <tr>
        <th align="left">report.html</th>
        <td>
            HTML report of the historical part of the pipeline: variant counts against the <code>ref_path</code> genomes (<i>variant_counts.tsv</i>), coverage of the <code>features_path</code> features (<i>coverage_feat.tsv</i>), and a summary of the strain identity database used (strains, native markers, tool elements). Also contains the warnings of the run.
        </td>
    </tr>
    <tr>
        <th align="left">figures</th>
        <td>
            Folder containing the coverage plots (<i>.png</i>), one per feature of <code>features_path</code> and per mating type sequence of <code>sex_path</code>, named <i>plot_&lt;SEQUENCE_NAME&gt;_mini.png</i>. Each plot shows the number of reads (mapping quality ≥ 20) covering each position of the sequence (<code>bedtools genomecov -bga</code>), the y-axis range being set by the <code>ylab</code> parameter of the <i>nextflow.config</i> file. The title indicates the number of reads after the mapping quality 20 filtering. A plot shows "EMPTY ... NO PLOT DRAWN" if too few positions are covered.
        </td>
    </tr>
    <tr>
        <th align="left">reports</th>
        <td>
            Folder containing:
            <ul>
                <li><b>*.log</b>: the log of each process of the pipeline (<i>prepare_reads_report.log</i>, <i>kmer_scan_report.log</i>, <i>bowtie2_*_report.log</i>, <i>q20_*_report.log</i>, <i>coverage_*_report.log</i>, <i>variantCall_*_report.log</i>, <i>tool_depth_report.log</i>, <i>tool_score_report.log</i>, <i>marker_align_report.log</i>, <i>marker_genotype_report.log</i>, <i>strain_score_report.log</i>, <i>fastq_dump_report.log</i> if an accession was downloaded, <i>aln_to_fastq_report.log</i> if a BAM, SAM or CRAM was given, <i>split_fasta.log</i>, <i>print_report.log</i>).</li>
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
                <li>call: strain called. resolution: how ties between strains not separated by the panel were resolved (tied strains, evidence from mating type and native markers).</li>
                <li>unassigned (true / false) and unassigned_reason.</li>
                <li>mating_type: mean count of the MATa and MATalpha k-mers, a_fraction (MATa / (MATa + MATalpha)) and call: MATa (a_fraction ≥ 0.58), MATalpha (a_fraction ≤ 0.42), MATa/alpha (in between) or unknown (too few k-mers seen).</li>
                <li>coverage: median depth of the called sites and fraction of uncovered sites.</li>
                <li>strains: the top strains of the ranking, with the columns of <i>strain_ranking.tsv</i> plus log_likelihood, sites_het (number of H sites), the site detail, and the number of sites separating each strain from the next one.</li>
                <li>tools and tool_summary: score of each tool element (columns of <i>tool_coverage.tsv</i> plus name, class, score, colour, present and origin: native or foreign to the host genome) and number of green, orange and red elements.</li>
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
        <th align="left">- marker_genotype.tsv</th>
        <td>
            Genotype of each native marker gene (URA3, HIS3, LEU2, TRP1, LYS2, ADE2, MET15, CAN1), from the reads aligned on the whole reference genome (<code>panel_reference</code>). One row per gene. For a heterozygous gene, only the first allele is shown here (both are in <i>marker_genotype.json</i>).
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
    <tr>
        <th align="left">- variant_counts.tsv</th>
        <td>
            Number of variants of the sample against each genome of <code>ref_path</code>, sorted by increasing number. The genome with the fewest variants is the most probable genetic background of the sample. One row per genome.
            <ul>
                <li>Name: name of the genome file (without extension).</li>
                <li>Nb_of_variants: number of variants (SNPs and indels) called by <code>bcftools mpileup -q 20 -Q 20 | bcftools call -mv</code>, minus 1.</li>
            </ul>
        </td>
    </tr>
    <tr>
        <th align="left">- coverage_feat.tsv</th>
        <td>
            Coverage of each feature of <code>features_path</code> by the reads of mapping quality ≥ 20 (bowtie2). One row per feature (features with too few covered positions to be plotted are absent).
            <ul>
                <li>REF: name of the feature (file name without extension).</li>
                <li>MAX: maximum number of reads covering a position of the feature.</li>
                <li>MEAN: mean of the depth values of the <code>bedtools genomecov -bga</code> intervals of the feature (intervals of constant depth, not weighted by their length).</li>
            </ul>
        </td>
    </tr>
    <tr>
        <th align="left">- coverage_sex.tsv</th>
        <td>
            As <i>coverage_feat.tsv</i>, but for the mating type sequences of <code>sex_path</code> (MATA, MATALPHA). The mating type call of the strain identity branch is in <i>strain_id.json</i>.
        </td>
    </tr>
    <tr>
        <th align="left">- variant_calling</th>
        <td>
            Folder containing, for each genome of <code>ref_path</code> and each mating type sequence of <code>sex_path</code>, the variants of the sample (<i>&lt;NAME&gt;_variants.tsv</i>, <code>bcftools query</code> of the VCF). One row per variant. The header line starts with #.
            <ul>
                <li>[1]CHROM: chromosome or sequence name.</li>
                <li>[2]POS: position of the variant (1-based).</li>
                <li>[3]REF: reference allele.</li>
                <li>[4]ALT: alternative allele(s), separated by commas.</li>
                <li>[5]QUAL: quality of the variant (phred scale).</li>
                <li>[6]&lt;SAMPLE&gt;:GT: genotype of the sample (0/1 heterozygous, 1/1 homozygous alternative, 1/2 two different alternative alleles).</li>
                <li>[7]&lt;SAMPLE&gt;:TYPE: SNP, MNP, INDEL or OTHER.</li>
            </ul>
        </td>
    </tr>
</table>


<br><br>
## VERSIONS


The different releases are [here](https://github.com/biosentinelle/biosentinelle/releases).

The different tagged versions are [here](https://github.com/biosentinelle/biosentinelle/tags).


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
- [Samtools and bcftools](https://www.htslib.org/)
- [Bedtools](https://bedtools.readthedocs.io/)
- [SRA Toolkit](https://github.com/ncbi/sra-tools)
- [Python](https://www.python.org/)
- [R](https://www.r-project.org/)
- [ggplot2](https://ggplot2.tidyverse.org/)
- [The 1011 Yeast Genomes project](http://1002genomes.u-strasbg.fr/)
- [Saccharomyces Genome Database](https://www.yeastgenome.org/)
- [Github](https://github.com/)
- [Bash](https://www.gnu.org/software/bash/)
- [Ubuntu](https://ubuntu.com/)


<br><br>
## WHAT'S NEW IN

#### v1.0

Everything
