#!/bin/bash
###############################################################################
# Module 7 — Single-Cell (Elective): environment-module loading.
#
# SOURCE this file from each script -- do not execute it:
#   SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"
#   source "${SCRIPT_DIR}/load_modules.sh"
#
# See lib/module_loader.sh for the try-pinned / discover-by-name /
# refuse-if-ambiguous fallback and how to override with e.g.
# STAR_MODULE=<exact-name>. The original scripts never pinned an exact STAR
# build (they said "verify exact version with module spider STAR") -- the
# pin below is an unverified best guess, not a tested value like this
# module's other tools; expect load_star() to fall back to discovery on
# most clusters, and treat a pinned-build hit as a coincidence, not
# confirmation this exact build was ever tested here.
###############################################################################

if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    echo "ERROR: load_modules.sh must be sourced, not executed." >&2
    exit 1
fi

_M7_LIB_DIR="${PIPELINE_SCRIPT_DIR:-${SLURM_SUBMIT_DIR:+${SLURM_SUBMIT_DIR}/modules/module-07-single-cell-elective/scripts}}"
_M7_LIB_DIR="${_M7_LIB_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
source "${_M7_LIB_DIR}/../../../lib/module_loader.sh"

load_sra_toolkit() { load_module_capability sra_toolkit "SRA-Toolkit/3.2.0-gompi-2024a" "SRA-Toolkit"; }
load_fastqc()      { load_module_capability fastqc      "FastQC/0.11.9-Java-11"         "FastQC"; }
load_multiqc()     { load_module_capability multiqc     "MultiQC/1.28-foss-2024a"        "MultiQC"; }
load_star()        { load_module_capability star        "STAR/2.7.11b-GCC-13.2.0"        "STAR"; }
# MultiQC's `env python` shebang can resolve to the wrong python3 if another
# module (R, QIIME2, ...) is loaded first -- load this immediately before
# load_multiqc, in the calling script, to force the right interpreter.
load_python()      { load_module_capability python      "Python/3.11.3-GCCcore-12.3.0"   "Python"; }
