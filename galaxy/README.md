# Galaxy tool (local Galaxy with planemo, galaxy.pasteur.fr)

Developer notes. To run the tool on galaxy.pasteur.fr: [HOW TO RUN](../README.md#how-to-run) section of the README.

The *galaxy/biosentinelle.xml* Galaxy tool runs the whole pipeline inside a single Galaxy job (nextflow local executor, processes in apptainer containers). Its form, i.e., the web page of the tool in Galaxy, recapitulates the *nextflow.config* file: same sections, one field per parameter, with the name of the parameter in parentheses and its comment as help. The *Run tool* button runs the pipeline. The strain identity database and the execution parameters (*Local / Cluster execution* section) are set by the Galaxy server.

**Results**: the *strain identity report* dataset of the Galaxy history is *strain_identity_report.html*, the main result. The *report and result folder* dataset is *report.html*, which also carries the whole result folder (see [OUTPUT](../README.md#output)): the eye icon shows the report, whose links to the tables and reports work in Galaxy, and **the download icon saves the whole result folder as a zip file**, where the web browser proposes (enable "Ask where to save each file" in the browser settings to choose the folder each time). Optionally, the tool also returns the .tsv tables as a collection (for other Galaxy tools) and the result folder as a separate zip dataset. On a Galaxy server, the datasets can also be exported to a storage configured by the administrators (Export datasets).

**Local Galaxy with planemo** ([planemo](https://planemo.readthedocs.io), e.g., in ~/.venv: `python3 -m venv ~/.venv && ~/.venv/bin/pip install planemo`). The strain identity database must be cloned next to this repository (`git clone https://github.com/biosentinelle/biosentinelle_yeast_db.git` from the parent folder), or its *db* folder given with the BS_DB_PATH variable:

<pre>
bash galaxy/local/serve.sh    # local Galaxy with the tool, at http://127.0.0.1:9091 (Ctrl+C to stop)
bash galaxy/local/test.sh     # lint (tool and Tool Shed repository) and tests of the tool (each test runs the whole pipeline)
python3 galaxy/local/watch.py # in a second terminal: progress of the nextflow processes of the running job (serve.sh or test.sh)
</pre>

Then, in the web browser: upload the sample files (reads or assembly) and the optional reference files (*Upload* button), open the *biosentinelle* tool, set the parameters and click on *Run tool*. The histories are kept between restarts (planemo profile). The first start installs Galaxy (release_26.1, in ~/.planemo, shared with the local Galaxy of other pipelines), the first run installs nextflow, apptainer, graphviz and python with conda (~/miniforge3), then pulls the container images of the pipeline into ~/apptainer (about 4 GB, long). Settings (cores, memory, folders) are environment variables, described in *galaxy/local/settings.sh*, e.g., `BS_CPUS=8 BS_MEMORY_MB=12000 bash galaxy/local/serve.sh`.

**Tool Shed and galaxy.pasteur.fr**: the *.shed.yml* file at the root of the repository describes the Tool Shed repository: the tool with its test data and the pipeline files (main.nf, nextflow.config, bin, conf, modules and references folders), in the same layout as this repository. After `bash galaxy/local/test.sh lint`, from the root of the repository (API key of the Tool Shed account):

<pre>
planemo shed_create --shed_target https://toolshed.pasteur.fr --shed_key API_KEY .                 # first upload
planemo shed_update --shed_target https://toolshed.pasteur.fr --shed_key API_KEY --check_diff .    # next versions
</pre>

The Galaxy administrators then install the repository from the Tool Shed and set the job environment variables described at the top of *galaxy/biosentinelle.xml*: BS_APPTAINER_CACHE (folder of the apptainer images, shared by the jobs) and BS_DB_PATH (*db* folder of a clone of [biosentinelle_yeast_db](https://github.com/biosentinelle/biosentinelle_yeast_db), the strain identity database). If the computing nodes have no internet access, the images are pulled beforehand into the image folder with `bash galaxy/pull_images.sh <folder>`, and the *Run accession* option of the form, which downloads the reads during the job, does not work.

- *galaxy/biosentinelle.xml*: the Galaxy tool. Inputs: the sample files (or forward and reverse reads, a paired collection, a run accession), optional reference, feature and mating type files, the parameters of *nextflow.config*, and an optional JSON parameter file (whose values replace those of the form, e.g., to reuse a set of parameters).
- *galaxy/fill_config.py*: writes the *nextflow.config* of a job from the form (parameters are in env {} blocks and cannot be changed with --param on the command line). User values are strictly checked (no quote, backslash or $), so that no Groovy code can be injected into the config.
- *galaxy/locked_params.txt*: parameters that users cannot set (input paths, database, execution, output location).
- *galaxy/pull_images.sh*: pulls the container images of *nextflow.config* into a folder, with the file names used by nextflow.
- *galaxy/test-data*: files of the tests of the tool. The sequence files are made from the references of the repository by *galaxy/local/make_test_data.py*.
- *galaxy/local/*: planemo scripts (not in the Tool Shed repository). *check_tool.py*, run by *test.sh*, checks that the form still matches *nextflow.config*: run it after any change of the parameters of *nextflow.config*.
