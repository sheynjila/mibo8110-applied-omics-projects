#!/bin/bash
#SBATCH --job-name=scrnaseq_smoke_test
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=24:00:00
#SBATCH --output=scrnaseq_smoke_%j.out
#SBATCH --error=scrnaseq_smoke_%j.err

# ==============================================================================
# SMOKE TEST WRAPPER for master_scrnaseq_pipeline.sh
# ==============================================================================
# WHY THIS SCRIPT EXISTS:
#   master_scrnaseq_pipeline.sh's own built-in fallback (used only when
#   srr_list.txt is missing) hardcodes SRR11092056/SRR11092057 - the same
#   "default test list" already verified elsewhere in this repository, via
#   NCBI eutils, as SARS-CoV-2 metagenomic RNA-seq, not single-cell 10x data
#   of any kind. Per this repository's "conserve what's already delivered"
#   policy, that fallback is NOT edited here; this wrapper supplies a real,
#   verified srr_list.txt instead, the same pattern used for every other
#   module's smoke test in this repository.
#
# A DISCLOSED FINDING ABOUT THIS MODULE'S SOURCE MATERIAL - the worked
# example this wrapper does NOT use, and why:
#   The source materials handed off for this module included two draft
#   scripts (PRJNA450892_GSE114725_scRNA_processing.sh and a duplicate
#   misfiled as PRJNA482620_TCGA_BRCA_alignment.sh - its own filename does
#   not match its content, which is entirely about PRJNA450892/GSE114725)
#   both hardcoding SRR7227261 and claiming it is BioProject PRJNA450892 /
#   GSE114725, "Single-cell map of immune phenotypes in breast cancer."
#   Verified via NCBI eutils, BOTH claims are wrong, independently:
#     1. SRR7227261 actually belongs to BioProject PRJNA438778 ("TCM
#        visualizes trajectories and cell populations from single cell
#        data"), a small Illumina MiSeq dataset (~129K spots - far too
#        small to be a real 10x single-cell run) with no connection to
#        GSE114725 at all.
#     2. GSE114725 itself (verified via its own GEO record) is really
#        SRA study SRP148597 / BioProject PRJNA472383 - and its actual
#        sequencing protocol is inDrop v2, NOT 10x Genomics. Running
#        master_scrnaseq_pipeline.sh's 10x-specific STARsolo
#        --soloCBwhitelist/CB_UMI_Simple logic against inDrop reads would
#        not error - per that script's own "METHODOLOGY NOTES" #2, a
#        chemistry mismatch is a SILENT failure (near-zero valid cells),
#        exactly the trap this smoke test exists to avoid walking into.
#   Neither draft script's claimed worked example is usable here, for two
#   independent reasons, not one - so this wrapper uses a different,
#   independently verified real 10x dataset instead (below).
#
# WORKED-EXAMPLE COHORT (verified real, 10x Chromium V3, Illumina NovaSeq
# 6000, via NCBI eutils esearch/esummary against BioProject 593571):
#   BioProject PRJNA593571 - "Benchmarking Single-Cell RNA Sequencing
#   Protocols for Cell Atlas Projects." Both runs are a deliberate
#   species-mixing ("barnyard") design: human PBMC (60%) + mouse colon
#   (30%) + HEK293T/NIH3T3/MDCK cell lines (10% combined), R1=28bp
#   (16bp cell barcode + 12bp UMI, confirming 10x v3 chemistry), R2=89bp
#   (cDNA).
#     SRR10587809  "10XV3_AN4491" - Chromium V3 WITH viability sorting
#     SRR10587810  "10XV3_AN4492" - Chromium V3 WITHOUT viability sorting
#   A genuine two-condition comparison (viability sorting on/off), not two
#   replicates of one condition. TEACHING NOTE: this dataset's barnyard
#   design is a real bonus for Module 7's own doublet-detection learning
#   outcome - a droplet containing both a human and a mouse cell's RNA is
#   an unambiguous, ground-truth doublet, which is the classic experimental
#   validation method for any doublet-detection algorithm. This module's
#   Tier B pipeline (see scripts/1-5) aligns only against the human
#   reference (Tier A's own default) and detects doublets computationally
#   (scDblFinder), NOT via species-mixing ground truth - extending Tier A
#   to a combined human+mouse reference to actually cross-check
#   scDblFinder's calls against this dataset's real barnyard ground truth
#   is a genuine, not-yet-built extension, noted in this module's README.
#
# Usage: sbatch smoke_test_prjna593571.sh
#        (run BEFORE committing srr_list.txt to a full cohort submission)
# ==============================================================================

set -e
set -o pipefail

echo "=========================================================="
echo " SINGLE-CELL PIPELINE SMOKE TEST — PRJNA593571"
echo " (2 real 10x Chromium V3 species-mixing benchmark runs)"
echo "=========================================================="

export CHEMISTRY="10x_v3"
WORKDIR="/scratch/$(whoami)/master_scrnaseq_pipeline"
mkdir -p ${WORKDIR}
cd ${WORKDIR}

cat > srr_list.txt <<'EOF'
SRR10587809
SRR10587810
EOF

echo "Wrote srr_list.txt:"
cat srr_list.txt

if [ ! -s "${WHITELIST_PATH:-/path/to/10x_whitelists/${CHEMISTRY}_barcodes.txt}" ]; then
    echo ""
    echo "NOTE: WHITELIST_PATH is not yet set to a real 10x v3 barcode whitelist file."
    echo "master_scrnaseq_pipeline.sh will refuse to run (by design, Phase 0) until"
    echo "WHITELIST_PATH points at one - see that script's header for where to obtain it."
fi

echo -e "\nLaunching master_scrnaseq_pipeline.sh (CHEMISTRY=${CHEMISTRY}) against the verified 2-sample cohort...\n"
bash "$(dirname "$0")/master_scrnaseq_pipeline.sh"

echo "=========================================================="
echo " SMOKE TEST COMPLETE — inspect ${WORKDIR}/alignment/<SRR>/Solo.out/Gene/filtered/"
echo " before submitting a larger cohort."
echo "=========================================================="
