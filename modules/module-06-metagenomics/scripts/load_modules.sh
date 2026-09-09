#!/bin/bash
###############################################################################
# Module 6 — Metagenomics (shotgun + amplicon): environment-module loading.
#
# SOURCE this file from each script -- do not execute it:
#   SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"
#   source "${SCRIPT_DIR}/load_modules.sh"
#
# See lib/module_loader.sh for the try-pinned / discover-by-name /
# refuse-if-ambiguous fallback and how to override with e.g.
# KRAKEN2_MODULE=<exact-name>.
###############################################################################

if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    echo "ERROR: load_modules.sh must be sourced, not executed." >&2
    exit 1
fi

_M6_LIB_DIR="${PIPELINE_SCRIPT_DIR:-${SLURM_SUBMIT_DIR:+${SLURM_SUBMIT_DIR}/modules/module-06-metagenomics/scripts}}"
_M6_LIB_DIR="${_M6_LIB_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
source "${_M6_LIB_DIR}/../../../lib/module_loader.sh"

load_sra_toolkit() { load_module_capability sra_toolkit "SRA-Toolkit/3.2.0-gompi-2024a"   "SRA-Toolkit"; }
load_fastqc()      { load_module_capability fastqc      "FastQC/0.11.9-Java-11"           "FastQC"; }
load_fastp()       { load_module_capability fastp       "fastp/0.23.4-GCC-13.2.0"         "fastp"; }
load_cutadapt()    { load_module_capability cutadapt    "cutadapt/4.9-GCCcore-13.3.0"     "cutadapt"; }
load_bowtie2()     { load_module_capability bowtie2     "Bowtie2/2.5.4-GCC-13.2.0"        "Bowtie2"; }
load_kraken2()     { load_module_capability kraken2     "Kraken2/2.1.2-gompi-2022a"       "Kraken2"; }
load_bracken()     { load_module_capability bracken     "Bracken/2.8-GCC-11.3.0"          "Bracken"; }
load_multiqc()     { load_module_capability multiqc     "MultiQC/1.28-foss-2024a"          "MultiQC"; }
# MultiQC's `env python` shebang can resolve to the wrong python3 if another
# module (R, QIIME2, ...) is loaded first -- load this immediately before
# load_multiqc, in the calling script, to force the right interpreter.
load_python()      { load_module_capability python      "Python/3.11.3-GCCcore-12.3.0"     "Python"; }
