#!/bin/bash
#SBATCH --job-name=iqtree_cohort
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=iqtree_%j.out
#SBATCH --error=iqtree_%j.err

# ==============================================================================
# MODULE 5 · SCRIPT 3 of 6 — TREE INFERENCE WITH IQ-TREE (ASCERTAINMENT-BIAS
# CORRECTED) AND ULTRAFAST BOOTSTRAP SUPPORT
# ==============================================================================
# ONE NEW IDEA on top of Script 2: Script 2 confirmed the alignment is worth
# building a tree from. This script builds the tree — and, critically, does
# NOT run IQ-TREE with a plain substitution model, because Script 1's
# alignment is SNP-ONLY (invariant sites already excluded). A plain model
# implicitly assumes every site in the alignment, including hypothetically
# invariant ones, had an equal chance to be observed; that assumption is
# false here by construction, and silently biases branch-length and model
# estimates (Lewis 2001's "ascertainment bias" problem). IQ-TREE's `+ASC`
# model correction exists specifically for this situation and is used here
# by default, not as an optional flag.
#
# TEACHING NOTE — model selection (MFP) + ultrafast bootstrap (-B) together:
#   `-m MFP+ASC` runs IQ-TREE's ModelFinder to pick the best-fitting
#   substitution model FOR THIS ALIGNMENT (never assume GTR or any other
#   model is automatically right), with ascertainment-bias correction
#   applied to every candidate model, not bolted on afterward. `-B 1000`
#   runs 1000 ultrafast bootstrap (UFBoot) replicates, producing a support
#   value at every internal branch. IQ-TREE's own guidance (Hoang et al.
#   2018) is that UFBoot support should be read against a ~95% cutoff for
#   "well supported" — NOT the classic Felsenstein bootstrap's traditional
#   70% cutoff, because UFBoot is a less conservative, bias-reduced
#   estimator by design. Script 4 uses the 95% threshold explicitly, and
#   states which kind of bootstrap it is reading, rather than reporting a
#   bare percentage with no interpretation attached.
#
# TEACHING NOTE — the curriculum's noted alternative (RAxML-NG), documented
# not run, mirrors Module 2's apeglm-fallback and Module 3's BUSCO/CheckM
# documentation pattern:
#   raxml-ng --all --msa core_snp_alignment.fasta --model GTR+ASC_LEWIS \
#     --bs-trees 1000 --threads 8 --prefix cohort_tree_raxml
#   RAxML-NG's `--all` runs tree search + bootstrapping in one invocation;
#   `--model GTR+ASC_LEWIS` is its equivalent ascertainment-bias correction.
#   Prefer this over IQ-TREE only if your curriculum/instructor specifically
#   requires RAxML-NG output format.
#
# Input:  core_snp_alignment.fasta (Script 1).
# Output: cohort_tree.treefile (Newick, UFBoot support at internal nodes),
#         cohort_tree.iqtree (full report incl. selected model),
#         cohort_tree.log.
# ==============================================================================

set -e
set -o pipefail

SCRIPT_DIR="${PIPELINE_SCRIPT_DIR:-${SLURM_SUBMIT_DIR:+${SLURM_SUBMIT_DIR}/modules/module-05-phylogenomics/scripts}}"
SCRIPT_DIR="${SCRIPT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)}"
source "${SCRIPT_DIR}/load_modules.sh"

ALIGNMENT="${1:-core_snp_alignment.fasta}"
PREFIX="${2:-cohort_tree}"
THREADS="${3:-8}"

echo "=================================================================="
echo " MODULE 5 / SCRIPT 3 — IQ-TREE inference (ModelFinder + ASC + UFBoot)"
echo "=================================================================="
echo "Alignment: ${ALIGNMENT}"
echo "Prefix:    ${PREFIX}"
echo ""

if [ ! -s "${ALIGNMENT}" ]; then
    echo "[FAIL] alignment not found or empty: ${ALIGNMENT}"
    exit 1
fi

module purge
load_iqtree || exit 1

iqtree2 -s "${ALIGNMENT}" \
        -m MFP+ASC \
        -B 1000 \
        -T "${THREADS}" \
        --prefix "${PREFIX}" \
        -redo

echo ""
echo "=================================================================="
echo " IQ-TREE COMPLETE"
echo " Tree:       ${PREFIX}.treefile"
echo " Full report (selected model, log-likelihood, etc.): ${PREFIX}.iqtree"
echo " Run Script 4 next to interpret the UFBoot support values -"
echo " do not read ${PREFIX}.treefile topology alone as if every branch"
echo " were equally trustworthy."
echo "=================================================================="
