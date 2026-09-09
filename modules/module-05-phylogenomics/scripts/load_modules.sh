#!/bin/bash
###############################################################################
# Module 5 — Phylogenomics: environment-module loading.
#
# SOURCE this file from each script -- do not execute it:
#   SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"
#   source "${SCRIPT_DIR}/load_modules.sh"
#
# See lib/module_loader.sh for the try-pinned / discover-by-name /
# refuse-if-ambiguous fallback and how to override with e.g.
# IQTREE_MODULE=<exact-name>.
#
# Scripts 2, 4, 5, and 6 are plain Python and have no module system of their
# own -- they only need Biopython/ete3, which are Python packages, not
# environment modules. Only script 1 (bcftools/samtools, via Module 4's
# load_modules.sh) and script 3 (IQ-TREE, below) need this file.
###############################################################################

if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    echo "ERROR: load_modules.sh must be sourced, not executed." >&2
    exit 1
fi

_M5_LIB_DIR="${PIPELINE_SCRIPT_DIR:-${SLURM_SUBMIT_DIR:+${SLURM_SUBMIT_DIR}/modules/module-05-phylogenomics/scripts}}"
_M5_LIB_DIR="${_M5_LIB_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
source "${_M5_LIB_DIR}/../../../lib/module_loader.sh"

load_iqtree()   { load_module_capability iqtree   "IQ-TREE/2.3.6-gompi-2024a" "IQ-TREE"; }
load_samtools() { load_module_capability samtools "SAMtools/1.21-GCC-13.3.0" "SAMtools"; }
load_bcftools() { load_module_capability bcftools "BCFtools/1.21-GCC-13.3.0" "BCFtools"; }
