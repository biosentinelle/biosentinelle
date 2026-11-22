#!/usr/bin/env bash
#########################################################################
##                                                                     ##
##     serve.sh                                                        ##
##                                                                     ##
##     Local Galaxy web page to run lab_strain_watch (planemo)         ##
##     Bioinformatics and Biostatistics Hub                            ##
##     Institut Pasteur Paris                                          ##
##                                                                     ##
#########################################################################

# Starts a local Galaxy with the lab_strain_watch tool (planemo serve), then:
#   1) open http://127.0.0.1:9091 in a web browser (port: LSW_PORT variable)
#   2) upload the sample files (reads or assembly) and the optional reference files (Upload button, top left)
#   3) open the lab_strain_watch tool: its form recapitulates nextflow.config, one field per parameter
#   4) click on "Run tool"
#   5) while the job runs (yellow), follow the nextflow processes:
#      - in the Galaxy page: click on the running dataset, then on its i icon (Dataset details): the Tool Standard Output
#        shows the lines of nextflow as they come (one line per started task, refreshed every 3 seconds)
#      - or in a second terminal, the progress table of the nextflow console: python3 galaxy/local/watch.py
#   6) when the jobs are green, open the "strain identity report" dataset (eye icon), and the "report and result folder"
#      dataset: the eye icon shows report.html, the download icon saves the whole result folder as a zip, where the web
#      browser proposes.
# The histories and results are kept between restarts, in the planemo profile (LSW_PROFILE variable).
# Stop the server with Ctrl+C.
#
# Usage, from any folder:
#   bash galaxy/local/serve.sh
# Settings: see settings.sh, plus
#   LSW_PORT     port of the web server (default: 9091, so that the local Galaxy of repertoire_profiler can use 9090)
#   LSW_HOST     address of the web server (default: 127.0.0.1, i.e., this computer only)
#   LSW_PROFILE  planemo profile keeping the Galaxy database and datasets (default: lab_strain_watch)
# Requirements: planemo (python3 -m venv ~/.venv && ~/.venv/bin/pip install planemo), git, internet access, and the strain
# identity database cloned next to this repository (see LSW_DB_PATH in settings.sh).
# The first start installs Galaxy (about 15 minutes), the first run of the tool installs nextflow and apptainer (conda),
# then pulls the container images of the pipeline (several GB, long).

source "$(dirname "${BASH_SOURCE[0]}")/settings.sh"

LSW_PORT="${LSW_PORT:-9091}"
LSW_HOST="${LSW_HOST:-127.0.0.1}"
LSW_PROFILE="${LSW_PROFILE:-lab_strain_watch}"

# profile: sqlite database and datasets kept in ~/.planemo/profiles/<LSW_PROFILE>
if ! "${PLANEMO}" profile_list 2>/dev/null | grep -q -x -F "${LSW_PROFILE}"; then
    "${PLANEMO}" profile_create "${LSW_PROFILE}"
fi

# the HTML reports keep their scripts and styles: the tool is in the sanitization allowlist of this Galaxy
# (on a Galaxy server, the administrators do the same with Admin > Allowlist)
printf 'lab_strain_watch\n' > "${LSW_WORK_DIR}/sanitize_allowlist.txt"
export GALAXY_CONFIG_OVERRIDE_SANITIZE_ALLOWLIST_FILE="${LSW_WORK_DIR}/sanitize_allowlist.txt"

echo "Galaxy will be at http://${LSW_HOST}:${LSW_PORT} (profile ${LSW_PROFILE})"
exec "${PLANEMO}" serve \
    "${PLANEMO_GALAXY_OPTIONS[@]}" \
    --profile "${LSW_PROFILE}" \
    --host "${LSW_HOST}" \
    --port "${LSW_PORT}" \
    "${LSW_TOOL_XML}"
