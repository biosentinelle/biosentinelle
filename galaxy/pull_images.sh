#!/usr/bin/env bash
#########################################################################
##                                                                     ##
##     pull_images.sh                                                  ##
##                                                                     ##
##     Pre-pulls the container images of lab_strain_watch              ##
##     Bioinformatics and Biostatistics Hub                            ##
##     Institut Pasteur Paris                                          ##
##                                                                     ##
#########################################################################

# Pulls the images of the container directives of nextflow.config into an apptainer image folder, with the file names
# used by nextflow (image name with / and : replaced by -, plus .img), so that nextflow finds them and does not pull.
# Use it to fill the LSW_APPTAINER_CACHE folder of the Galaxy tool (or the apptainer_local_images / apptainer_cluster_images
# folders of nextflow.config) from a computer with internet access, when the computing nodes have none.
# Images already present are skipped.
#
# Usage:
#   bash galaxy/pull_images.sh <image_folder> [number of parallel pulls, default 3]
# Requirements: apptainer (or singularity) in the PATH, internet access, about 4 GB in the image folder.
# Temporary files: APPTAINER_TMPDIR variable (default: <image_folder>/.pull_tmp).

set -euo pipefail

if [[ $# -lt 1 ]]; then
    echo "Usage: bash $0 <image_folder> [parallel pulls]" >&2
    exit 1
fi
CONFIG="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/nextflow.config"
export IMAGE_DIR="$(mkdir -p "$1" && cd "$1" && pwd)"
PARALLEL="${2:-3}"
export APPTAINER_TMPDIR="${APPTAINER_TMPDIR:-${IMAGE_DIR}/.pull_tmp}"
mkdir -p "${APPTAINER_TMPDIR}"
export APPTAINER_BIN="$(command -v apptainer || command -v singularity || true)"
if [[ -z "${APPTAINER_BIN}" ]]; then
    echo "ERROR: apptainer (or singularity) not found in the PATH" >&2
    exit 1
fi

pull() {
    local image="$1"
    local name
    name="$(echo "${image}" | sed 's#[:/]#-#g').img"
    if [[ -s "${IMAGE_DIR}/${name}" ]]; then
        echo "present: ${name}"
        return 0
    fi
    echo "pulling: docker://${image} -> ${name}"
    if "${APPTAINER_BIN}" pull --disable-cache --name "${IMAGE_DIR}/${name}.pulling" "docker://${image}" > "${APPTAINER_TMPDIR}/${name}.log" 2>&1; then
        mv "${IMAGE_DIR}/${name}.pulling" "${IMAGE_DIR}/${name}"
        rm -f "${APPTAINER_TMPDIR}/${name}.log"
        echo "pulled:  ${name}"
    else
        rm -f "${IMAGE_DIR}/${name}.pulling"
        echo "FAILED:  ${name} (see ${APPTAINER_TMPDIR}/${name}.log)"
        return 1
    fi
}
export -f pull

grep -o -E "container *= *['\"][^'\"]+['\"]" "${CONFIG}" | sed -E "s/container *= *['\"]//; s/['\"]$//" | sort -u \
    | xargs -P "${PARALLEL}" -I{} bash -c 'pull "$@"' _ {}
echo "Images in ${IMAGE_DIR}:"
ls -1sh "${IMAGE_DIR}"/*.img
