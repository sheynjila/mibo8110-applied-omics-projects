#!/bin/bash
###############################################################################
# Module 3 — Genome Assembly: environment-module loading.
#
# SOURCE this file from each script -- do not execute it:
#   SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"
#   source "${SCRIPT_DIR}/load_modules.sh"
#
# NOTE: this module's scripts previously called `spades.py` / `flye` /
# `minimap2` / `samtools` with no `module load` at all -- they only worked
# if those tools happened to already be on PATH. The pins below are
# reasonable defaults (not yet verified against a real cluster's module
# tree) -- run `module spider <name>` and override via e.g.
# SPADES_MODULE=<exact-name> if the pinned build below isn't available.
#
# Scripts 4-7 are plain Python (no #SBATCH/module system of their own) --
# `module load` in a shell cannot reach into a Python subprocess, so if you
# run them directly, source this file and call the relevant load_* function
# in your *shell* first, so spades.py/samtools/minimap2/busco/checkm are on
# PATH before Python starts.
###############################################################################

if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    echo "ERROR: load_modules.sh must be sourced, not executed." >&2
    exit 1
fi

_M3_LIB_DIR="${PIPELINE_SCRIPT_DIR:-${SLURM_SUBMIT_DIR:+${SLURM_SUBMIT_DIR}/modules/module-03-genome-assembly/scripts}}"
_M3_LIB_DIR="${_M3_LIB_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
source "${_M3_LIB_DIR}/../../../lib/module_loader.sh"

load_spades()   { load_module_capability spades   "SPAdes/3.15.5-GCC-13.2.0"     "SPAdes"; }
load_flye()     { load_module_capability flye     "Flye/2.9.6-GCC-13.2.0"        "Flye"; }
load_minimap2() { load_module_capability minimap2 "minimap2/2.28-GCCcore-13.2.0" "minimap2"; }
load_samtools() { load_module_capability samtools "SAMtools/1.21-GCC-13.3.0"     "SAMtools"; }
load_busco()    { load_module_capability busco    "BUSCO/5.7.1-foss-2023a"       "BUSCO"; }
load_checkm()   { load_module_capability checkm   "CheckM/1.2.2-foss-2022a"      "CheckM"; }
