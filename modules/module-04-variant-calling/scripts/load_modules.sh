#!/bin/bash
###############################################################################
# Module 4 — Variant Calling: environment-module loading.
#
# SOURCE this file from each script -- do not execute it:
#   SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"
#   source "${SCRIPT_DIR}/load_modules.sh"
#
# Pinned versions below match what master_snp_pipeline.sh (Tier A) already
# used; see lib/module_loader.sh for the try-pinned / discover-by-name /
# refuse-if-ambiguous fallback and how to override with e.g.
# SAMTOOLS_MODULE=<exact-name>.
###############################################################################

if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    echo "ERROR: load_modules.sh must be sourced, not executed." >&2
    exit 1
fi

_M4_LIB_DIR="${PIPELINE_SCRIPT_DIR:-${SLURM_SUBMIT_DIR:+${SLURM_SUBMIT_DIR}/modules/module-04-variant-calling/scripts}}"
_M4_LIB_DIR="${_M4_LIB_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
source "${_M4_LIB_DIR}/../../../lib/module_loader.sh"

load_sra_toolkit() { load_module_capability sra_toolkit "SRA-Toolkit/3.2.0-gompi-2024a" "SRA-Toolkit"; }
load_fastp()       { load_module_capability fastp       "fastp/0.23.4-GCC-13.2.0"       "fastp"; }
load_bwa()         { load_module_capability bwa         "BWA/0.7.18-GCC-13.3.0"         "BWA"; }
load_samtools()    { load_module_capability samtools    "SAMtools/1.21-GCC-13.3.0"       "SAMtools"; }
load_bcftools()    { load_module_capability bcftools    "BCFtools/1.21-GCC-13.3.0"       "BCFtools"; }
load_multiqc()     { load_module_capability multiqc     "MultiQC/1.28-foss-2024a"        "MultiQC"; }
# MultiQC's `env python` shebang can resolve to the wrong python3 if another
# module (R, QIIME2, ...) is loaded first -- load this immediately before
# load_multiqc, in the calling script, to force the right interpreter.
load_python()      { load_module_capability python      "Python/3.11.3-GCCcore-12.3.0"   "Python"; }
