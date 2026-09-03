#!/bin/bash
#SBATCH --job-name=amplicon_smoke_test
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --output=amplicon_smoke_%j.out
#SBATCH --error=amplicon_smoke_%j.err

# ==============================================================================
# SMOKE TEST WRAPPER for master_amplicon_pipeline.sh (16S route)
# ==============================================================================
# WHY THIS SCRIPT EXISTS:
#   master_amplicon_pipeline.sh has no built-in default sample list at all
#   (see that script's own PHASE 1 error message) - so unlike the shotgun
#   smoke test above, this wrapper is not replacing a wrong hardcoded
#   default, it is supplying the only real worked example this pipeline has.
#
# WORKED-EXAMPLE COHORT (verified real, paired-end, Illumina MiSeq, via
# NCBI eutils esearch/esummary against BioProject 917645):
#   BioProject PRJNA917645 - "Shallow shotgun sequencing reduces technical
#   variation in microbiome analysis" (human gut/stool metagenome). This
#   study sequenced the SAME underlying stool samples by both shallow
#   shotgun AND 16S amplicon - the shotgun arm of this exact BioProject is
#   this module's shotgun-route worked example too (see
#   smoke_test_shotgun_prjna1520859.sh's header for why a DIFFERENT
#   BioProject was used there instead: this study's own shotgun runs are
#   single-end, and master_microbiome_pipeline.sh's Tier A pipeline is
#   paired-end throughout). The 16S arm used here IS paired-end and runs
#   cleanly through this pipeline as delivered:
#     SRR22959853  human gut 16S rRNA amplicon (V4), MiSeq, paired-end
#     SRR22959854  human gut 16S rRNA amplicon (V4), MiSeq, paired-end
#     SRR22959855  human gut 16S rRNA amplicon (V4), MiSeq, paired-end
#
#   These are human (animal-associated, not plant-associated) stool samples.
#   amplicon_dada2_analysis.R's chloroplast/mitochondrial ASV filter (Phase F)
#   is 16S-generic and will run correctly against this cohort, but a human
#   stool sample is not expected to yield chloroplast hits the way a plant
#   tissue or rhizosphere 16S sample would - do not read "few or no
#   chloroplast ASVs removed" here as evidence the filter itself is broken;
#   it is evidence this particular worked example has little plant material
#   in it. Point AMPLICON_TYPE=16S / this pipeline at a real plant- or
#   soil-associated 16S dataset if you want to see that filter do more work.
#
# Usage: sbatch smoke_test_amplicon_prjna917645.sh
#        (run BEFORE committing srr_list.txt to a full cohort submission)
# ==============================================================================

set -e
set -o pipefail

echo "=========================================================="
echo " AMPLICON PIPELINE SMOKE TEST (16S) — PRJNA917645"
echo " (3 real human-gut 16S V4 amplicon runs)"
echo "=========================================================="

export AMPLICON_TYPE="16S"
WORKDIR="/scratch/$(whoami)/master_amplicon_pipeline_${AMPLICON_TYPE}"
mkdir -p ${WORKDIR}
cd ${WORKDIR}

cat > srr_list.txt <<'EOF'
SRR22959853
SRR22959854
SRR22959855
EOF

echo "Wrote srr_list.txt:"
cat srr_list.txt

echo -e "\nLaunching master_amplicon_pipeline.sh (AMPLICON_TYPE=${AMPLICON_TYPE}) against the verified 3-sample cohort...\n"
bash "$(dirname "$0")/master_amplicon_pipeline.sh"

echo "=========================================================="
echo " SMOKE TEST COMPLETE — trimmed reads in ${WORKDIR}/trimmed_reads/"
echo " Next: run amplicon_dada2_analysis.R (AMPLICON_TYPE=\"16S\") against this directory."
echo "=========================================================="
