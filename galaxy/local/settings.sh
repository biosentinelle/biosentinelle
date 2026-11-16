#!/usr/bin/env bash
#########################################################################
##                                                                     ##
##     settings.sh                                                     ##
##                                                                     ##
##     Common settings of serve.sh and test.sh (sourced, not run)      ##
##     Bioinformatics and Biostatistics Hub                            ##
##     Institut Pasteur Paris                                          ##
##                                                                     ##
#########################################################################

# Settings of the local Galaxy used by planemo to run the biosentinelle tool (galaxy/biosentinelle.xml).
# All of them can be changed with environment variables, e.g., BS_CPUS=8 bash galaxy/local/serve.sh
#   PLANEMO               planemo executable (default: planemo of the PATH, else ~/.venv/bin/planemo)
#   BS_GALAXY_BRANCH      Galaxy release used by planemo (default: release_26.1)
#   BS_GALAXY_ROOT        folder of this Galaxy, cloned at the first use (default: ~/.planemo/galaxy_<BS_GALAXY_BRANCH>)
#   BS_CONDA_PREFIX       conda used by Galaxy to install the requirements of the tool: nextflow, apptainer, graphviz and python
#                         (default: ~/miniforge3, installed at the first use)
#   BS_APPTAINER_CACHE    folder of the apptainer images of the pipeline (default: ~/apptainer). The first run pulls them (long)
#   BS_APPTAINER_OPTIONS  apptainer run options of the pipeline processes (default: --no-home -B /run/shm:/run/shm, as the local
#                         runOptions of nextflow.config)
#   BS_DB_PATH            strain identity database: db folder of a clone of https://github.com/biosentinelle/biosentinelle_yeast_db
#                         (default: biosentinelle_yeast_db/db next to this repository, as the db_path default of nextflow.config)
#   BS_CPUS               cores of a biosentinelle job (default: all the cores of the computer)
#   BS_MEMORY_MB          memory of a biosentinelle job in MB (default: 80% of the memory of the computer)
#   BS_WORK_DIR           folder of the generated job configuration and of the test reports (default: ~/.planemo/biosentinelle)
#   BS_TMPDIR             temporary folder of planemo and of the Galaxy jobs, on disk (default: BS_WORK_DIR/tmp)

set -euo pipefail

BS_TOOL_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)" # galaxy folder of the repository
BS_TOOL_XML="${BS_TOOL_DIR}/biosentinelle.xml"

PLANEMO="${PLANEMO:-$(command -v planemo || echo "${HOME}/.venv/bin/planemo")}"
BS_GALAXY_BRANCH="${BS_GALAXY_BRANCH:-release_26.1}"
BS_GALAXY_ROOT="${BS_GALAXY_ROOT:-${HOME}/.planemo/galaxy_${BS_GALAXY_BRANCH}}"
BS_CONDA_PREFIX="${BS_CONDA_PREFIX:-${HOME}/miniforge3}"
BS_APPTAINER_CACHE="${BS_APPTAINER_CACHE:-${HOME}/apptainer}"
BS_APPTAINER_OPTIONS="${BS_APPTAINER_OPTIONS:---no-home -B /run/shm:/run/shm}"
BS_DB_PATH="${BS_DB_PATH:-$(cd "${BS_TOOL_DIR}/../.." && pwd)/biosentinelle_yeast_db/db}"
BS_CPUS="${BS_CPUS:-$(nproc)}"
BS_MEMORY_MB="${BS_MEMORY_MB:-$(awk '/^MemTotal:/ {printf "%d", $2 / 1024 * 0.8}' /proc/meminfo)}"
BS_WORK_DIR="${BS_WORK_DIR:-${HOME}/.planemo/biosentinelle}"

if [[ ! -x "${PLANEMO}" ]]; then
    echo "ERROR: planemo not found. Install it (python3 -m venv ~/.venv && ~/.venv/bin/pip install planemo) or set the PLANEMO variable." >&2
    exit 1
fi
for i1 in BS_CPUS BS_MEMORY_MB; do
    if [[ ! "${!i1}" =~ ^[1-9][0-9]*$ ]]; then
        echo "ERROR: ${i1} must be a positive integer: ${!i1}" >&2
        exit 1
    fi
done
for i1 in BS_APPTAINER_CACHE BS_APPTAINER_OPTIONS BS_DB_PATH BS_WORK_DIR BS_CONDA_PREFIX; do
    if [[ "${!i1}" == *'"'* || "${!i1}" == *'\'* ]]; then
        echo "ERROR: ${i1} cannot contain a double quote or a backslash: ${!i1}" >&2
        exit 1
    fi
done
if [[ ! -d "${BS_DB_PATH}/yeast" ]]; then
    echo "ERROR: no strain identity database in ${BS_DB_PATH} (BS_DB_PATH variable). Clone it next to this repository:" >&2
    echo "    git clone https://github.com/biosentinelle/biosentinelle_yeast_db.git $(cd "${BS_TOOL_DIR}/../.." && pwd)/biosentinelle_yeast_db" >&2
    exit 1
fi

mkdir -p "${BS_APPTAINER_CACHE}" "${BS_WORK_DIR}"

# temporary folder of planemo, hence of the Galaxy jobs (work folder of nextflow, builds of the apptainer images):
# on disk, since /tmp can be a small folder in memory (tmpfs, e.g., in WSL)
export TMPDIR="${BS_TMPDIR:-${BS_WORK_DIR}/tmp}"
mkdir -p "${TMPDIR}"

# Galaxy used by planemo: a fixed release, kept between runs (planemo installs its virtual environment and web client at the first use)
if [[ ! -d "${BS_GALAXY_ROOT}" ]]; then
    echo "Cloning Galaxy ${BS_GALAXY_BRANCH} into ${BS_GALAXY_ROOT}"
    # core.autocrlf=false: a global core.autocrlf=true would give CRLF line endings, which break the Galaxy shell scripts
    git -c core.autocrlf=false clone --depth 1 --branch "${BS_GALAXY_BRANCH}" https://github.com/galaxyproject/galaxy.git "${BS_GALAXY_ROOT}"
fi

# Galaxy job configuration: other tools run with 1 core, the biosentinelle tool runs with BS_CPUS cores
# (GALAXY_SLOTS) and BS_MEMORY_MB MB (GALAXY_MEMORY_MB), with the apptainer and database settings read by the tool command.
# live_tool_output_reporting: the Galaxy page shows the standard output of the running job (one line per nextflow task),
# in the dataset details (i icon of the running dataset), Tool Standard Output. Progress table in a terminal: watch.py
BS_JOB_CONF="${BS_WORK_DIR}/job_conf.yml"
cat > "${BS_JOB_CONF}" <<EOF
# Generated by ${BS_TOOL_DIR}/local/settings.sh: changes are overwritten at the next run
runners:
  local:
    load: galaxy.jobs.runners.local:LocalJobRunner
    workers: 4
execution:
  default: local
  environments:
    local:
      runner: local
      tmp_dir: true
    biosentinelle:
      runner: local
      local_slots: ${BS_CPUS}
      tmp_dir: true
      live_tool_output_reporting: true
      env:
        - name: GALAXY_MEMORY_MB
          value: "${BS_MEMORY_MB}"
        - name: BS_APPTAINER_CACHE
          value: "${BS_APPTAINER_CACHE}"
        - name: BS_APPTAINER_OPTIONS
          value: "${BS_APPTAINER_OPTIONS}"
        - name: BS_DB_PATH
          value: "${BS_DB_PATH}"
tools:
  - class: local
    environment: local
  - id: biosentinelle
    environment: biosentinelle
EOF

# options common to planemo serve and planemo test
PLANEMO_GALAXY_OPTIONS=(
    --galaxy_root "${BS_GALAXY_ROOT}"
    --job_config_file "${BS_JOB_CONF}"
    --conda_prefix "${BS_CONDA_PREFIX}"
    --conda_auto_init
    --conda_auto_install
)

echo "planemo:          ${PLANEMO}"
echo "Galaxy:           ${BS_GALAXY_ROOT}"
echo "conda:            ${BS_CONDA_PREFIX}"
echo "apptainer images: ${BS_APPTAINER_CACHE} (options: ${BS_APPTAINER_OPTIONS})"
echo "database:         ${BS_DB_PATH}"
echo "job resources:    ${BS_CPUS} cores, ${BS_MEMORY_MB} MB"
