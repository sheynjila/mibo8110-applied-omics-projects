#!/bin/bash
###############################################################################
# Module 1 — Raw Reads & QC: environment-module loading.
#
# SOURCE this file from each script -- do not execute it:
#   SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"
#   source "${SCRIPT_DIR}/load_modules.sh"
#
# Then, per phase (matching this module's existing `module purge` habit):
#   module purge
#   load_sra_toolkit || exit 1
#
# Each pinned version below is what the scripts in this folder were authored
# and tested against; on a different cluster the loader falls back to
# discovering the tool by name (see lib/module_loader.sh for the exact
# try-pinned / discover-by-name / refuse-if-ambiguous rule, and how to force
# a specific module with e.g. FASTQC_MODULE=<name>).
###############################################################################

if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    echo "ERROR: load_modules.sh must be sourced, not executed." >&2
    exit 1
fi

_M1_LIB_DIR="${PIPELINE_SCRIPT_DIR:-${SLURM_SUBMIT_DIR:+${SLURM_SUBMIT_DIR}/modules/module-01-raw-reads-qc/scripts}}"
_M1_LIB_DIR="${_M1_LIB_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
source "${_M1_LIB_DIR}/../../../lib/module_loader.sh"

load_sra_toolkit() { load_module_capability sra_toolkit "SRA-Toolkit/3.0.3-gompi-2022a" "SRA-Toolkit"; }
load_fastqc()      { load_module_capability fastqc      "FastQC/0.11.9-Java-11"        "FastQC"; }
load_fastp()       { load_module_capability fastp       "fastp/0.23.4-GCC-13.2.0"      "fastp"; }
load_multiqc()     { load_module_capability multiqc     "MultiQC/1.28-foss-2024a"      "MultiQC"; }
