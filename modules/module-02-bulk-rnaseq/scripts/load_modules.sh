#!/bin/bash
###############################################################################
# Module 2 — Bulk RNA-seq: environment-module loading.
#
# SOURCE this file from each script -- do not execute it:
#   SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"
#   source "${SCRIPT_DIR}/load_modules.sh"
#
# Only master_rnaseq_pipeline_consolidated.sh (Tier A) needs this -- the
# seven Tier B R scripts (1-7) are run interactively in R/RStudio and manage
# their own package availability via `library()`, not environment modules.
#
# See lib/module_loader.sh for the try-pinned / discover-by-name /
# refuse-if-ambiguous fallback and how to override with e.g.
# HISAT2_MODULE=<exact-name>.
###############################################################################

if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    echo "ERROR: load_modules.sh must be sourced, not executed." >&2
    exit 1
fi

_M2_LIB_DIR="${PIPELINE_SCRIPT_DIR:-${SLURM_SUBMIT_DIR:+${SLURM_SUBMIT_DIR}/modules/module-02-bulk-rnaseq/scripts}}"
_M2_LIB_DIR="${_M2_LIB_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
source "${_M2_LIB_DIR}/../../../lib/module_loader.sh"

load_sra_toolkit() { load_module_capability sra_toolkit "SRA-Toolkit/3.0.3-gompi-2022a" "SRA-Toolkit"; }
load_fastqc()      { load_module_capability fastqc      "FastQC/0.11.9-Java-11"         "FastQC"; }
load_fastp()       { load_module_capability fastp       "fastp/0.23.4-GCC-13.2.0"       "fastp"; }
load_hisat2()      { load_module_capability hisat2      "HISAT2/2.2.1-gompi-2022a"      "HISAT2"; }
load_samtools()    { load_module_capability samtools    "SAMtools/1.16.1-GCC-11.3.0"    "SAMtools"; }
load_subread()     { load_module_capability subread     "Subread/2.0.6-GCC-12.3.0"      "Subread"; }
load_multiqc()     { load_module_capability multiqc     "MultiQC/1.28-foss-2024a"        "MultiQC"; }
# MultiQC's `env python` shebang can resolve to the wrong python3 if another
# module (R, QIIME2, ...) is loaded first -- load this immediately before
# load_multiqc, in the calling script, to force the right interpreter. Also
# needed before HISAT2's own Python dependency (hisat2-build/hisat2 scripts).
load_python()      { load_module_capability python      "Python/3.11.3-GCCcore-12.3.0"   "Python"; }
