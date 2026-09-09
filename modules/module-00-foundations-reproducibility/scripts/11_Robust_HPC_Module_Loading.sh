#!/bin/bash
###############################################################################
# 11_Robust_HPC_Module_Loading.sh
#
# NEW IDEA on top of scripts 1-10: every earlier script in this module
# called `module load <exact-build-string>` directly, e.g.
# `module load FastQC/0.11.9-Java-11`. That is fine on the cluster this
# course was authored on. It is NOT fine the moment a student runs the same
# script on a different HPC site, because module trees are not
# standardized across institutions -- the same tool is commonly available
# as a different exact build string, or under a different capitalization,
# or not at all under that name. A bare `module load` then fails the whole
# job with a message that says nothing about what to do next.
#
# THIS SCRIPT teaches the fix this project uses everywhere from here on:
# lib/module_loader.sh, sourced by a small per-module load_modules.sh, is
# a *hybrid* loader:
#   1. Try the exact pinned build first (what the course was tested on).
#   2. If that's missing, search `module avail <loose-name>` for this
#      cluster's equivalent.
#   3. If exactly ONE candidate turns up, load it -- but WARN loudly, since
#      a different build can (rarely) give different results, and that
#      needs to be visible in your run log / report, not silently hidden.
#   4. If ZERO or MULTIPLE candidates turn up, REFUSE to guess. Ambiguity
#      here is exactly the kind of silent-wrong-answer risk Module 1's
#      "check inputs before trusting outputs" habit exists to prevent --
#      picking the wrong module version is a garbage-in problem one layer
#      below the data.
#
# This script does NOT re-teach that logic from scratch (it already lives,
# fully commented, in ../../../lib/module_loader.sh). Instead it runs it
# against this module's own five tools so you can see the three outcomes
# (pinned hit / single-candidate fallback / refusal) for real, on your own
# cluster, before you rely on it in scripts 10 and every module after this
# one.
#
# USAGE:
#   bash 11_Robust_HPC_Module_Loading.sh
#
# Run this on the login node (it only calls `module`, nothing heavy).
###############################################################################

SCRIPT_DIR="${PIPELINE_SCRIPT_DIR:-${SLURM_SUBMIT_DIR:+${SLURM_SUBMIT_DIR}/modules/module-00-foundations-reproducibility/scripts}}"
SCRIPT_DIR="${SCRIPT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)}"
source "${SCRIPT_DIR}/load_modules.sh"

echo "=================================================================="
echo " MODULE 0 / SCRIPT 11 — Robust HPC module loading, demonstrated"
echo "=================================================================="
echo "Sourced: ${SCRIPT_DIR}/load_modules.sh"
echo "Which itself sourced: ${SCRIPT_DIR}/../../../lib/module_loader.sh"
echo ""
echo "Each line below is one call to load_module_capability(). Watch for"
echo "three different outcomes as this runs on YOUR cluster:"
echo "  [load_modules] <tool>: loaded pinned module '...'          <- exact match"
echo "  [load_modules] WARNING: ... using discovered module '...'  <- fallback (1 candidate)"
echo "  [load_modules] ERROR: ...                                  <- refused (0 or 2+ candidates)"
echo ""

module purge

FAILED_TOOLS=()
for tool_loader in load_sra_toolkit load_fastqc load_fastp load_multiqc load_samtools; do
    echo "------------------------------------------------------------------"
    echo " ${tool_loader}"
    echo "------------------------------------------------------------------"
    if ! "${tool_loader}"; then
        FAILED_TOOLS+=("${tool_loader}")
    fi
done

echo ""
echo "=================================================================="
if [ ${#FAILED_TOOLS[@]} -eq 0 ]; then
    echo " All five tools resolved to a loadable module on this cluster."
    echo " If any of them printed a WARNING above (fallback, not exact"
    echo " match), note that in your Module 0 reflection -- it is exactly"
    echo " the kind of environment difference script 7's reproducibility"
    echo " record exists to surface."
else
    echo " ${#FAILED_TOOLS[@]} tool(s) could not be resolved automatically:"
    printf '   %s\n' "${FAILED_TOOLS[@]}"
    echo ""
    echo " This is not a bug in this script -- it is the loader correctly"
    echo " refusing to guess. Run 'module spider <tool-name>' yourself to"
    echo " find the exact name on this cluster, then re-run with e.g.:"
    echo "   SAMTOOLS_MODULE=<exact-name-you-found> bash 11_Robust_HPC_Module_Loading.sh"
fi
echo "=================================================================="

###############################################################################
# WHAT TO DO WITH THIS, GOING FORWARD:
#   - Every module from here on (1, 3, 4, 5, 6, 7, 8, and the supplemental
#     case studies) ships its own scripts/load_modules.sh that follows this
#     exact pattern, tuned to that module's own tools and pinned versions.
#   - When you write your OWN pipeline script that needs a module, do not
#     write a bare `module load <build-string>`. Add one load_<tool>()
#     function to that module's load_modules.sh (or lib/module_loader.sh
#     directly, for a brand-new tool) and call it instead.
###############################################################################
