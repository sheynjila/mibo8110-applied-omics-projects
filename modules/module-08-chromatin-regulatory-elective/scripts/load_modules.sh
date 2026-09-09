#!/bin/bash
###############################################################################
# Module 8 — Chromatin & Regulatory Genomics (Elective): environment-module
# loading.
#
# SOURCE this file from each script -- do not execute it:
#   SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"
#   source "${SCRIPT_DIR}/load_modules.sh"
#
# See lib/module_loader.sh for the try-pinned / discover-by-name /
# refuse-if-ambiguous fallback and how to override with e.g.
# MACS2_MODULE=<exact-name>.
###############################################################################

if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    echo "ERROR: load_modules.sh must be sourced, not executed." >&2
    exit 1
fi

_M8_LIB_DIR="${PIPELINE_SCRIPT_DIR:-${SLURM_SUBMIT_DIR:+${SLURM_SUBMIT_DIR}/modules/module-08-chromatin-regulatory-elective/scripts}}"
_M8_LIB_DIR="${_M8_LIB_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
source "${_M8_LIB_DIR}/../../../lib/module_loader.sh"

load_sra_toolkit() { load_module_capability sra_toolkit "SRA-Toolkit/3.2.0-gompi-2024a" "SRA-Toolkit"; }
load_fastp()       { load_module_capability fastp       "fastp/0.23.4-GCC-13.2.0"       "fastp"; }
load_bowtie2()     { load_module_capability bowtie2     "Bowtie2/2.5.4-GCC-13.2.0"       "Bowtie2"; }
load_samtools()    { load_module_capability samtools    "SAMtools/1.18-GCC-12.3.0"       "SAMtools"; }
load_macs2()       { load_module_capability macs2       "MACS2/2.2.9.1-foss-2022a"       "MACS2"; }
load_bedtools()    { load_module_capability bedtools    "BEDTools/2.31.0-GCC-12.3.0"     "BEDTools"; }
load_multiqc()     { load_module_capability multiqc     "MultiQC/1.28-foss-2024a"         "MultiQC"; }
# MultiQC's `env python` shebang can resolve to the wrong python3 if another
# module (R, QIIME2, ...) is loaded first -- load this immediately before
# load_multiqc, in the calling script, to force the right interpreter.
load_python()      { load_module_capability python      "Python/3.11.3-GCCcore-12.3.0"    "Python"; }
