#!/bin/bash
###############################################################################
# Supplemental Case Studies (AMR + MTB): environment-module loading.
#
# SOURCE this file from each script -- do not execute it:
#   SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"
#   source "${SCRIPT_DIR}/load_modules.sh"
#
# See lib/module_loader.sh for the try-pinned / discover-by-name /
# refuse-if-ambiguous fallback and how to override with e.g.
# AMRFINDERPLUS_MODULE=<exact-name>.
#
# load_sra_toolkit() takes an OPTIONAL pinned-version argument because the
# two case studies here were authored against different SRA-Toolkit builds
# (Final_master_end_to_end_amr.sh: 3.0.3-gompi-2022a;
#  PRJNA1425489_mtb_single_end.sh: 3.2.0-gompi-2024a) -- pass the version
# each script already expects rather than picking one arbitrarily:
#   load_sra_toolkit "SRA-Toolkit/3.0.3-gompi-2022a"
###############################################################################

if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    echo "ERROR: load_modules.sh must be sourced, not executed." >&2
    exit 1
fi

_MSUPP_LIB_DIR="${PIPELINE_SCRIPT_DIR:-${SLURM_SUBMIT_DIR:+${SLURM_SUBMIT_DIR}/modules/supplemental-case-studies/scripts}}"
_MSUPP_LIB_DIR="${_MSUPP_LIB_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
source "${_MSUPP_LIB_DIR}/../../../lib/module_loader.sh"

load_sra_toolkit() {
    local pinned="${1:-SRA-Toolkit/3.2.0-gompi-2024a}"
    load_module_capability sra_toolkit "${pinned}" "SRA-Toolkit"
}
load_fastp()          { load_module_capability fastp          "fastp/0.23.4-GCC-13.2.0"              "fastp"; }
load_fastqc()         { load_module_capability fastqc         "FastQC/0.11.9-Java-11"                "FastQC"; }
load_spades()         { load_module_capability spades         "SPAdes/3.15.5-GCC-13.2.0"             "SPAdes"; }
load_quast()          { load_module_capability quast          "QUAST/5.2.0-foss-2022a"               "QUAST"; }
load_amrfinderplus()  { load_module_capability amrfinderplus  "ncbi-amrfinderplus/3.11.11-foss-2022a" "amrfinderplus"; }
load_hisat2()         { load_module_capability hisat2         "HISAT2/2.2.1-gompi-2022a"             "HISAT2"; }
load_samtools()       { load_module_capability samtools       "SAMtools/1.18-GCC-12.3.0"             "SAMtools"; }
load_subread()        { load_module_capability subread        "Subread/2.0.6-GCC-12.3.0"             "Subread"; }
load_multiqc()        { load_module_capability multiqc        "MultiQC/1.28-foss-2024a"               "MultiQC"; }
load_python()         { load_module_capability python         "Python/3.11.3-GCCcore-12.3.0"          "Python"; }
