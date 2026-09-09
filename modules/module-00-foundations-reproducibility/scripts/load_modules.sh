#!/bin/bash
###############################################################################
# Module 0 — Foundations & Reproducibility: environment-module loading.
#
# SOURCE this file from a script -- do not execute it:
#   SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"
#   source "${SCRIPT_DIR}/load_modules.sh"
#
# Scripts 1-9 in this module intentionally call bare `module load
# <exact-version>` -- that IS the lesson at that point in the course (see
# script 11, "Robust HPC Module Loading", for why a bare pinned `module
# load` breaks on a different cluster, and how this file's helper functions
# fix that). Only script 10 (the capstone, which integrates everything) and
# script 11 (which teaches this pattern directly) source this file.
#
# See lib/module_loader.sh for the underlying try-pinned / discover-by-name
# / refuse-if-ambiguous fallback and how to override with e.g.
# SAMTOOLS_MODULE=<exact-name>.
###############################################################################

if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    echo "ERROR: load_modules.sh must be sourced, not executed." >&2
    exit 1
fi

_M0_LIB_DIR="${PIPELINE_SCRIPT_DIR:-${SLURM_SUBMIT_DIR:+${SLURM_SUBMIT_DIR}/modules/module-00-foundations-reproducibility/scripts}}"
_M0_LIB_DIR="${_M0_LIB_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
source "${_M0_LIB_DIR}/../../../lib/module_loader.sh"

load_sra_toolkit() { load_module_capability sra_toolkit "SRA-Toolkit/3.0.3-gompi-2022a" "SRA-Toolkit"; }
load_fastqc()      { load_module_capability fastqc      "FastQC/0.11.9-Java-11"         "FastQC"; }
load_fastp()       { load_module_capability fastp       "fastp/0.23.4-GCC-13.2.0"       "fastp"; }
load_multiqc()     { load_module_capability multiqc     "MultiQC/1.28-foss-2024a"        "MultiQC"; }
load_samtools()    { load_module_capability samtools    "SAMtools/1.16.1-GCC-11.3.0"     "SAMtools"; }
