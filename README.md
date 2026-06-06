
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
## INPUT

| Required files |
| :--- |
| A single file (zipped or not) containing multiple nucleotide fasta sequences. Allowed extensions for this file are: fasta, fa, fas, fna, txt, seq or faa.<br>It can also be a folder (zipped or not) containing nucleotide fasta files, each containing a single sequence. Allowed extensions are the same.<br>Use table2fasta.R ([https://github.com/gael-millot/table2fasta](https://github.com/gael-millot/table2fasta)) if sequences are in a .table file.<br>In fasta files, the sequence can be split into several lines (\n and or \r separated). In addition, spaces and tabs can be present in the header (they will be replaced by an underscore).<br>Sequences do not have to be Immunonoglobulin, BCR or TCR sequences: filtering is performed by the pipeline. |
| A metadata file (optional) for adding informations in the results. |

<br>

The dataset used in the *nextflow.config* file, as an example, is available at https://zenodo.org/records/18631337/files/human_IGH.zip.

<br>

The metadata file used in the *nextflow.config* file, as an example, is available at https://zenodo.org/records/18294088/files/human_IGH_metadata.tsv.

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


