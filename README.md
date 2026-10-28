<img src="logo.png" align="center" width="300" alt="Biosentinelle Logo">
<img src="images/repertoire_profiler.gif" align="right" width="300" alt="Demo GIF">

| Usage | Requirement |
| :--- | :--- |
| [![Nextflow](https://img.shields.io/badge/code-Nextflow-blue?style=plastic)](https://www.nextflow.io/) | [![Dependencies: Nextflow Version](https://img.shields.io/badge/Nextflow-v24.10.4-blue?style=plastic)](https://github.com/nextflow-io/nextflow) |
| [![License: GPL-3.0](https://img.shields.io/badge/licence-GPL%20(%3E%3D3)-green?style=plastic)](https://www.gnu.org/licenses) | [![Dependencies: Apptainer Version](https://img.shields.io/badge/Apptainer-v1.3.5-blue?style=plastic)](https://github.com/apptainer/apptainer) |
| | [![Dependencies: Graphviz Version](https://img.shields.io/badge/Graphviz-v2.42.3-blue?style=plastic)](https://www.graphviz.org/download/) |

<br><br>
## TABLE OF CONTENTS


   - [AIM](#aim)
   - [WARNING](#warning)
   - [CONTENT](#content)
   - [HOW IT WORKS](#how-it-works)
   - [INPUT](#input)
   - [HOW TO RUN](#how-to-run)
   - [OUTPUT](#output)
   - [VERSIONS](#versions)
   - [LICENCE](#licence)
   - [CITATION](#citation)
   - [CREDITS](#credits)
   - [ACKNOWLEDGEMENTS](#Acknowledgements)
   - [WHAT'S NEW IN](#what's-new-in)

<br><br>
## AIM

- Check of the lab model strain sequenced.

<br><br>
## WARNINGS

- Right now, only dedicated to Saccharomyces cerevisiae.

<br><br>
## CONTENT

| Files and folder | Description |
| :--- | :--- |
| **main.nf** | File that can be executed using a linux terminal, a MacOS terminal or Windows 10 WSL2. |
| **nextflow.config** | Parameter settings for the *main.nf* file. Users have to open this file, set the desired settings and save these modifications before execution. Of note, this configuration file is systematically saved in the reports folder (see [below](#output)) during each execution, to save the parameter settings. |
| **bin folder** | Contains files required by the *main.nf* file. |
| **modules folder** | Contains files required by the *main.nf* file. |
| **conf folder** | Contains files required by the *main.nf* file. |
| **Licence.txt** | Licence of the release. |


<br><br>
## HOW IT WORKS

No assembly is made. The sample sequences are compared with a strain identity database (`db/<organism>`, built once by *build_db.nf* from a reference genome and the sequences of 1,172 *S. cerevisiae* strains), through three independent kinds of evidence:

| Evidence | How | Reference genome used at run time? |
| :--- | :--- | :--- |
| **Strain SNP and indel panel** | Alignment-free. For each panel site, the database holds the k-mers of each allele. The sample sequences are streamed once and the k-mers of each allele are counted, giving one allele per site where the depth is at least `kmer_min_depth`. The same is done for the mating type (MATa / MATalpha) k-mers. | No |
| **Native marker genes** (URA3, HIS3, LEU2, TRP1, LYS2, ADE2, MET15, CAN1) | The reads are aligned with bowtie2 against the whole reference genome (`panel_reference`), so that reads from paralogs and repeats get a low mapping quality. The allele of each marker is read from the alignment and compared with the alleles catalogued in the database, judged by the protein it encodes. | Yes |
| **Lab tool elements** (plasmids, cassettes, drug resistance markers, ...) | The reads are aligned with bowtie2 against the tool sequences of the database. Each element is scored by breadth of coverage and identity. | No, only the tool sequences |

The panel strains are then ranked by likelihood. If the best strain matches less than `min_concordance` of its expected alleles, the sample is reported as unassigned rather than forced onto the closest strain, which is what flags a strain absent from the database. The result is the *strain_identity_report.html* file.

In addition, the historical part of the pipeline aligns the reads against the genomes of `ref_path`, the features of `features_path` and the mating type sequences of `sex_path`, for coverage plots, Q20 read counts and variant counts (bcftools), summarized in *report.html*.

<br>

## INPUT

| Parameter | Description |
| :--- | :--- |
| `fastq_path` | The sample sequences, despite the parameter name not only FASTQ. Either a run accession of SRA (SRR...), ENA (ERR...) or DDBJ (DRR...), downloaded with fastq-dump, or a file, or a folder of files. All the files of a folder are analyzed together as a single sample. |
| `seq_kind` | `"reads"` (default) for sequencing reads, short or long, or `"assembly"` for contigs, scaffolds or a complete genome sequence. |
| `seq_molecule` | `"DNA"` (default) or `"RNA"` (RNA-seq or direct RNA reads). The analysis is the same, but the report warns that depths then follow gene expression, and that a strain site or a marker gene that is not expressed is reported as not covered. |
| `cram_reference` | Only for a CRAM sample: the FASTA reference it was compressed against. If "NULL", samtools finds the reference from its checksum, which needs network access. |

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

<br>

<br><br>
## HOW TO RUN

### 1. Prerequisite

Installation of:<br>
[nextflow DSL2](file:///C:/Users/gmillot/Documents/Git_projects/protocols/docs/Protocol%20152-rev0%20DSL2.htm#_Toc208504071). Please, use the version indicated above.<br>
[Graphviz](https://www.graphviz.org/download/), `sudo apt install graphviz` for Linux ubuntu.<br>
[Apptainer](https://gael-millot.github.io/protocols/docs/Protocol%20135-rev0%20APPTAINER.html#_Toc160091693).<br>
<br>

Optional installation (to avoid reccurent message) of:<br>
[Gocryptfs](https://github.com/rfjakob/gocryptfs), `sudo apt install gocryptfs` for Linux ubuntu.<br> 

Itol key:<br>
If you need sequence phylogenic trees in the output, you can freely register at https://itol.embl.de/itol_account.cgi to get your own itol key. Once registered, go to https://itol.embl.de/userInfo.cgi and click on the Toggle API access button. Then, add the key in the `phylo_tree_itolkey` parameter of the *nextflow.config* file, and set the `phylo_tree_itol_subscription` parameter to `TRUE`. If you experience problem with registration, set the `phylo_tree_itol_subscription` parameter to `FALSE`. The html output file explains how to see the trees without ITOL key.

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
#### 2.2. *main.nf* file in a public github / gitlab repository

Run the following command from where you want the results:

<pre>
nextflow run -hub pasteur gmillot/repertoire_profiler -r v1.0.0
</pre>

<br><br>

### 3. Distant running (example with the Pasteur cluster)

#### 3.1. Pre-execution

Go into the directory where the main.nf and nextflow.config files are.
Copy-paste this code:

<pre>
EXEC_PATH=$(pwd) # where the bin folder of 19583_loot is located (by default, the same path as for the main.nf file)
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
HOME="${HELIXHOME}/repertoire_profiler/" # $HOME changed to allow the creation of .nextflow into /$HELIXHOME/repertoire_profiler/, for instance. See NFX_HOME in the nextflow software script
nextflow run main.nf -c nextflow.config -resume --modules ${MODULES} # --modules ${MODULES} in order to have all the used module versions recorded into the report file, -resume to continue a potential job that has already been submitted
HOME=$HOME_INI
</pre>

<br><br>

#### 3.3. *main.nf* file in the public gitlab repository

Modify the first line of the code below, in order to call the desired version of repertoire_profiler indicated [here](https://github.com/gael-millot/repertoire_profiler/tags), and run (results will be where the EXEC_PATH variable has been set above):

<pre>
VERSION="v1.0"
HOME_INI=$HOME
HOME="${HELIXHOME}/repertoire_profiler/" # $HOME changed to allow the creation of .nextflow into /$HELIXHOME/repertoire_profiler/, for instance. See NFX_HOME in the nextflow software script
nextflow run -hub pasteur gmillot/repertoire_profiler -r $VERSION -c $HOME/nextflow.config
HOME=$HOME_INI
</pre>

<br><br>
### 4. Error messages and solutions

#### Message 1
```
Unknown error accessing project `gmillot/repertoire_profiler` -- Repository may be corrupted: /pasteur/sonic/homes/gmillot/.nextflow/assets/gmillot/repertoire_profiler
```

Purge using:
<pre>
rm -rf /pasteur/sonic/homes/gmillot/.nextflow/assets/gmillot*
</pre>

#### Message 2
```
WARN: Cannot read project manifest -- Cause: Remote resource not found: https://gitlab.pasteur.fr/api/v4/projects/gmillot%2Frepertoire_profiler
```

Contact Gael Millot (distant repository is not public).

#### Message 3

```
permission denied
```

Use chmod to change the user rights. Example linked to files in the bin folder: 
```
chmod 755 bin/*.*
```

#### Message 4

```
ERROR ~ Error executing process > 'print_report (1)'

Caused by:
  Process `print_report (1)` terminated with an error exit status (1)

...

Command error:
  Quitting from lines 280-301 (report_file.rmd)
  Error in if (file.info(file_path)$size == 0) { :
    missing value where TRUE/FALSE needed
  Calls: <Anonymous> ... eval_with_user_handlers -> eval -> eval -> read_tsv_with_dummy
  Execution halted
```

If obtained using `nextflow run main.nf -resume`, then rerun the same command once. Sometimes, nextflow shows difficulties to build the report.html file with `-resume`. If the problem persists, try `nextflow run main.nf`. If the problem still occurs, please send an issue [here](https://github.com/gael-millot/repertoire_profiler/issues).

#### Message 5

```
ERROR ~ Error executing process > 'ITOL (2)'

Caused by:
  Process `ITOL (2)` terminated with an error exit status (1)
  
INFO:    underlay of /etc/localtime required more than 50 (83) bind mounts
```

Register at Itol as explained in the [Prerequisite](#1-prerequisite) section above, or set the `phylo_tree_itol_subscription` parameter of the *nextflow.config* file to `FALSE` and rerun.

#### Error using -resume

Do not hesitate to rerun once the same command to reobtain the same error message. Indeed, sometimes nextflow has difficulties to reassemble all the cached files and reruning solves the problem.




<br><br>
## OUTPUT

By default, all the results are returned in a *result* folder where the *main.nf* executed file is located (created if does not exist). This can be changed using the `out_path_ini` parameter of the *nextflow.config* file. By default, each execution produces a new folder named *repertoire_profiler_\<ID\>*, created inside the *result* folder and containing all the outputs of the execution. The name of the folder can be changed using the `result_folder_name` parameter of the *nextflow.config* file. The new name file will be followed by an \<ID\> in all cases.
<br><br>
An example of results obtained with the dataset is present at this address: https://zenodo.org/records/18633133/files/repertoire_profiler_1770995657.zip.


<br><br>
## VERSIONS


The different releases are [here](https://github.com/gael-millot/repertoire_profiler/releases).

The different tagged versions are [here](https://github.com/gael-millot/repertoire_profiler/tags).


<br><br>
## LICENCE


This package of scripts can be redistributed and/or modified under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
Distributed in the hope that it will be useful, but without any warranty; without even the implied warranty of merchantability or fitness for a particular purpose.
See the GNU General Public License for more details at https://www.gnu.org/licenses or in the Licence.txt attached file.


<br><br>
## CITATION


<br><br>
## CREDITS

[Gael A. Millot](https://gitlab.pasteur.fr/gmillot), Institut Pasteur, Université Paris Cité, Bioinformatics and Biostatistics Hub, 75015 Paris, France

<br><br>
## ACKNOWLEDGEMENTS


The developers & maintainers of the mentioned softwares and packages, including:

- [R](https://www.r-project.org/)
- [ggplot2](https://ggplot2.tidyverse.org/)
- [immcantation](https://immcantation.readthedocs.io/en/stable/)
- [ggtree](https://yulab-smu.top/treedata-book/)
- [Python](https://www.python.org/)
- [Abalign](http://cao.labshare.cn/abalign/)
- [Mafft](https://mafft.cbrc.jp/alignment/server/index.html)
- [Nextflow](https://www.nextflow.io/)
- [Apptainer](https://apptainer.org/)
- [Docker](https://www.docker.com/)
- [Gitlab](https://about.gitlab.com/)
- [Github](https://github.com/)
- [Bash](https://www.gnu.org/software/bash/)
- [Ubuntu](https://ubuntu.com/)

Special acknowledgement to:


<br><br>
## WHAT'S NEW IN

#### v1.0

Everything


