#!/bin/bash
#SBATCH --job-name=atacseq_smoke_test
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=24:00:00
#SBATCH --output=atacseq_smoke_%j.out
#SBATCH --error=atacseq_smoke_%j.err

# ==============================================================================
# SMOKE TEST WRAPPER for master_chipseq_pipeline.sh (ATAC mode)
# ==============================================================================
# WHY THIS SCRIPT EXISTS: same reason as smoke_test_chip_prjna1283515.sh -
# master_chipseq_pipeline.sh's own hardcoded SRR11092056/SRR11092057
# fallback is SARS-CoV-2 RNA-seq, not ATAC-seq. This wrapper supplies a
# real, verified srr_list.txt instead.
#
# WORKED-EXAMPLE COHORT (verified real, paired-end, Illumina NextSeq 2000,
# via NCBI eutils esearch/esummary against BioProject 1283515):
#   BioProject PRJNA1283515 - the SAME study and, deliberately, the SAME
#   T3M4 cell line as smoke_test_chip_prjna1283515.sh's CTCF ChIP-seq
#   cohort, so this module's two worked examples can be read side by side
#   (accessible chromatin vs. CTCF binding sites, same cells) rather than
#   being two unrelated datasets:
#     SRR35978262  T3M4 ATACseq Rep1
#     SRR35978261  T3M4 ATACseq Rep2
#     SRR35978259  T3M4 ATACseq Rep3
#
#   ATAC-seq has no input-control concept (see master_chipseq_pipeline.sh's
#   header) - CONTROL_SRR is not applicable in ATAC mode regardless of
#   whether it is set.
#
# Usage: sbatch smoke_test_atac_prjna1283515.sh
#        (run BEFORE committing srr_list.txt to a full cohort submission)
# ==============================================================================

set -e
set -o pipefail

echo "=========================================================="
echo " CHROMATIN PIPELINE SMOKE TEST (ATAC) — PRJNA1283515"
echo " (3 real biological replicates, T3M4 cell line)"
echo "=========================================================="

export ASSAY_TYPE="ATAC"
WORKDIR="/scratch/$(whoami)/master_chipseq_pipeline_${ASSAY_TYPE}"
mkdir -p ${WORKDIR}
cd ${WORKDIR}

cat > srr_list.txt <<'EOF'
SRR35978262
SRR35978261
SRR35978259
EOF

echo "Wrote srr_list.txt:"
cat srr_list.txt

echo -e "\nLaunching master_chipseq_pipeline.sh (ASSAY_TYPE=${ASSAY_TYPE}) against the verified 3-replicate cohort...\n"
bash "$(dirname "$0")/master_chipseq_pipeline.sh"

echo "=========================================================="
echo " SMOKE TEST COMPLETE — inspect ${WORKDIR}/macs2_output/*/*_peaks.narrowPeak"
echo " before submitting a larger cohort."
echo "=========================================================="
