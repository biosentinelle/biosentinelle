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

# Settings of the local Galaxy used by planemo to run the lab_strain_watch tool (galaxy/lab_strain_watch.xml).
# All of them can be changed with environment variables, e.g., LSW_CPUS=8 bash galaxy/local/serve.sh
#   PLANEMO               planemo executable (default: planemo of the PATH, else ~/.venv/bin/planemo)
#   LSW_GALAXY_BRANCH     Galaxy release used by planemo (default: release_26.1)
#   LSW_GALAXY_ROOT       folder of this Galaxy, cloned at the first use (default: ~/.planemo/galaxy_<LSW_GALAXY_BRANCH>)
#   LSW_CONDA_PREFIX      conda used by Galaxy to install the requirements of the tool: nextflow, apptainer, graphviz and python
#                         (default: ~/miniforge3, installed at the first use)
#   LSW_APPTAINER_CACHE   folder of the apptainer images of the pipeline (default: ~/apptainer). The first run pulls them (long)
#   LSW_APPTAINER_OPTIONS apptainer run options of the pipeline processes (default: --no-home -B /run/shm:/run/shm, as the local
#                         runOptions of nextflow.config)
#   LSW_DB_PATH           strain identity database: db folder of a clone of https://github.com/lab-strain-watch/yeast_db
#                         (default: yeast_db/db next to this repository, as the db_path default of nextflow.config)
#   LSW_CPUS              cores of a lab_strain_watch job (default: all the cores of the computer)
#   LSW_MEMORY_MB         memory of a lab_strain_watch job in MB (default: 80% of the memory of the computer)
#   LSW_WORK_DIR          folder of the generated job configuration and of the test reports (default: ~/.planemo/lab_strain_watch)
#   LSW_TMPDIR            temporary folder of planemo and of the Galaxy jobs, on disk (default: LSW_WORK_DIR/tmp)

set -euo pipefail

LSW_TOOL_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)" # galaxy folder of the repository
LSW_TOOL_XML="${LSW_TOOL_DIR}/lab_strain_watch.xml"

PLANEMO="${PLANEMO:-$(command -v planemo || echo "${HOME}/.venv/bin/planemo")}"
LSW_GALAXY_BRANCH="${LSW_GALAXY_BRANCH:-release_26.1}"
LSW_GALAXY_ROOT="${LSW_GALAXY_ROOT:-${HOME}/.planemo/galaxy_${LSW_GALAXY_BRANCH}}"
LSW_CONDA_PREFIX="${LSW_CONDA_PREFIX:-${HOME}/miniforge3}"
LSW_APPTAINER_CACHE="${LSW_APPTAINER_CACHE:-${HOME}/apptainer}"
LSW_APPTAINER_OPTIONS="${LSW_APPTAINER_OPTIONS:---no-home -B /run/shm:/run/shm}"
LSW_DB_PATH="${LSW_DB_PATH:-$(cd "${LSW_TOOL_DIR}/../.." && pwd)/yeast_db/db}"
LSW_CPUS="${LSW_CPUS:-$(nproc)}"
LSW_MEMORY_MB="${LSW_MEMORY_MB:-$(awk '/^MemTotal:/ {printf "%d", $2 / 1024 * 0.8}' /proc/meminfo)}"
LSW_WORK_DIR="${LSW_WORK_DIR:-${HOME}/.planemo/lab_strain_watch}"

if [[ ! -x "${PLANEMO}" ]]; then
    echo "ERROR: planemo not found. Install it (python3 -m venv ~/.venv && ~/.venv/bin/pip install planemo) or set the PLANEMO variable." >&2
    exit 1
fi
for i1 in LSW_CPUS LSW_MEMORY_MB; do
    if [[ ! "${!i1}" =~ ^[1-9][0-9]*$ ]]; then
        echo "ERROR: ${i1} must be a positive integer: ${!i1}" >&2
        exit 1
    fi
done
for i1 in LSW_APPTAINER_CACHE LSW_APPTAINER_OPTIONS LSW_DB_PATH LSW_WORK_DIR LSW_CONDA_PREFIX; do
    if [[ "${!i1}" == *'"'* || "${!i1}" == *'\'* ]]; then
        echo "ERROR: ${i1} cannot contain a double quote or a backslash: ${!i1}" >&2
        exit 1
    fi
done
if [[ ! -d "${LSW_DB_PATH}/yeast" ]]; then
    echo "ERROR: no strain identity database in ${LSW_DB_PATH} (LSW_DB_PATH variable). Clone it next to this repository:" >&2
    echo "    git clone https://github.com/lab-strain-watch/yeast_db.git $(cd "${LSW_TOOL_DIR}/../.." && pwd)/yeast_db" >&2
    exit 1
fi

mkdir -p "${LSW_APPTAINER_CACHE}" "${LSW_WORK_DIR}"

# temporary folder of planemo, hence of the Galaxy jobs (work folder of nextflow, builds of the apptainer images):
# on disk, since /tmp can be a small folder in memory (tmpfs, e.g., in WSL)
export TMPDIR="${LSW_TMPDIR:-${LSW_WORK_DIR}/tmp}"
mkdir -p "${TMPDIR}"

# Galaxy used by planemo: a fixed release, kept between runs (planemo installs its virtual environment and web client at the first use)
if [[ ! -d "${LSW_GALAXY_ROOT}" ]]; then
    echo "Cloning Galaxy ${LSW_GALAXY_BRANCH} into ${LSW_GALAXY_ROOT}"
    # core.autocrlf=false: a global core.autocrlf=true would give CRLF line endings, which break the Galaxy shell scripts
    git -c core.autocrlf=false clone --depth 1 --branch "${LSW_GALAXY_BRANCH}" https://github.com/galaxyproject/galaxy.git "${LSW_GALAXY_ROOT}"
fi

# Galaxy job configuration: other tools run with 1 core, the lab_strain_watch tool runs with LSW_CPUS cores
# (GALAXY_SLOTS) and LSW_MEMORY_MB MB (GALAXY_MEMORY_MB), with the apptainer and database settings read by the tool command.
# live_tool_output_reporting: the Galaxy page shows the standard output of the running job (one line per nextflow task),
# in the dataset details (i icon of the running dataset), Tool Standard Output. Progress table in a terminal: watch.py
LSW_JOB_CONF="${LSW_WORK_DIR}/job_conf.yml"
cat > "${LSW_JOB_CONF}" <<EOF
# Generated by ${LSW_TOOL_DIR}/local/settings.sh: changes are overwritten at the next run
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
    lab_strain_watch:
      runner: local
      local_slots: ${LSW_CPUS}
      tmp_dir: true
      live_tool_output_reporting: true
      env:
        - name: GALAXY_MEMORY_MB
          value: "${LSW_MEMORY_MB}"
        - name: LSW_APPTAINER_CACHE
          value: "${LSW_APPTAINER_CACHE}"
        - name: LSW_APPTAINER_OPTIONS
          value: "${LSW_APPTAINER_OPTIONS}"
        - name: LSW_DB_PATH
          value: "${LSW_DB_PATH}"
tools:
  - class: local
    environment: local
  - id: lab_strain_watch
    environment: lab_strain_watch
EOF

# options common to planemo serve and planemo test
PLANEMO_GALAXY_OPTIONS=(
    --galaxy_root "${LSW_GALAXY_ROOT}"
    --job_config_file "${LSW_JOB_CONF}"
    --conda_prefix "${LSW_CONDA_PREFIX}"
    --conda_auto_init
    --conda_auto_install
)

echo "planemo:          ${PLANEMO}"
echo "Galaxy:           ${LSW_GALAXY_ROOT}"
echo "conda:            ${LSW_CONDA_PREFIX}"
echo "apptainer images: ${LSW_APPTAINER_CACHE} (options: ${LSW_APPTAINER_OPTIONS})"
echo "database:         ${LSW_DB_PATH}"
echo "job resources:    ${LSW_CPUS} cores, ${LSW_MEMORY_MB} MB"
