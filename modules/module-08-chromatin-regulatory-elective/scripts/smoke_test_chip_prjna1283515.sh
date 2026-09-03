#!/bin/bash
#SBATCH --job-name=chipseq_smoke_test
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=24:00:00
#SBATCH --output=chipseq_smoke_%j.out
#SBATCH --error=chipseq_smoke_%j.err

# ==============================================================================
# SMOKE TEST WRAPPER for master_chipseq_pipeline.sh (ChIP mode)
# ==============================================================================
# WHY THIS SCRIPT EXISTS:
#   master_chipseq_pipeline.sh's own built-in fallback (used only when
#   srr_list.txt is missing) hardcodes SRR11092056/SRR11092057 - the same
#   "default test list" already verified elsewhere in this repository, via
#   NCBI eutils, as SARS-CoV-2 metagenomic RNA-seq, not ChIP-seq of any
#   kind. This wrapper supplies a real, verified srr_list.txt instead, the
#   same pattern used for every other module's smoke test in this repository.
#
# WORKED-EXAMPLE COHORT (verified real, paired-end, Illumina NextSeq 2000,
# via NCBI eutils esearch/esummary against BioProject 1283515):
#   BioProject PRJNA1283515 - "KLF5 controls subtype-independent highly
#   interactive enhancers in pancreatic cancer to regulate cell survival."
#   CTCF ChIP-seq, T3M4 pancreatic cancer cell line, 3 real biological
#   replicates (not technical replicates of one run):
#     SRR35978255  T3M4 CTCF Rep1
#     SRR35978254  T3M4 CTCF Rep2
#     SRR35978253  T3M4 CTCF Rep3
#
# A REAL, DISCLOSED LIMITATION OF THIS WORKED EXAMPLE — read before running:
#   This BioProject's 149 deposited runs (checked via NCBI eutils: no hit
#   for "IgG", "Input", or "control" anywhere in the project) include NO
#   input/IgG control sample for ANY of its ChIP-seq experiments (CTCF,
#   dNp63, H3K27ac, across three cell lines). This is not a search failure
#   on this wrapper's part - it is a real, published study that did not
#   deposit a matched input control. CONTROL_SRR is deliberately left UNSET
#   below as a result - see master_chipseq_pipeline.sh's own header
#   TEACHING NOTE for what running without one actually means, and this
#   module's README "Worked-example dataset" for the full account. Do NOT
#   treat this smoke test's peak calls as equivalent to a properly
#   input-controlled ChIP-seq run; Tier B Script 1 flags this condition
#   explicitly rather than silently proceeding as if it were fine.
#
# Usage: sbatch smoke_test_chip_prjna1283515.sh
#        (run BEFORE committing srr_list.txt to a full cohort submission)
# ==============================================================================

set -e
set -o pipefail

echo "=========================================================="
echo " CHROMATIN PIPELINE SMOKE TEST (ChIP, CTCF) — PRJNA1283515"
echo " (3 real biological replicates, T3M4 cell line, NO input control)"
echo "=========================================================="

export ASSAY_TYPE="ChIP"
# CONTROL_SRR intentionally left unset - see header above.
WORKDIR="/scratch/$(whoami)/master_chipseq_pipeline_${ASSAY_TYPE}"
mkdir -p ${WORKDIR}
cd ${WORKDIR}

cat > srr_list.txt <<'EOF'
SRR35978255
SRR35978254
SRR35978253
EOF

echo "Wrote srr_list.txt:"
cat srr_list.txt

echo -e "\nLaunching master_chipseq_pipeline.sh (ASSAY_TYPE=${ASSAY_TYPE}, no input control) against the verified 3-replicate cohort...\n"
bash "$(dirname "$0")/master_chipseq_pipeline.sh"

echo "=========================================================="
echo " SMOKE TEST COMPLETE — inspect ${WORKDIR}/macs2_output/*/*_peaks.narrowPeak"
echo " before submitting a larger cohort."
echo "=========================================================="
