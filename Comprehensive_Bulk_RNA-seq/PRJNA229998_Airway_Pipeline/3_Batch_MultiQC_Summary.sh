#!/bin/bash
#SBATCH --job-name=qc03_multiqc
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=01:30:00

###############################################################################
# 3_Batch_MultiQC_Summary.sh
#
# NARRATIVE:
# This script aggregates individual FastQC outputs into a single interactive 
# MultiQC dashboard. In a real-world multi-sample study, manually reviewing 
# dozens of individual QC HTML files is prone to human error. MultiQC compiles 
# these metrics into parallel tracks, instantly highlighting systemic batch 
# effects or isolated sample failures prior to preprocessing.
###############################################################################
set -e
set -o pipefail

# V9 FIX: Purge stale modules before loading toolchains
module purge
module load MultiQC/1.28-foss-2024a

WORKDIR="/scratch/$(whoami)/PRJNA229998_airway_pipeline_stepbystep"
cd "${WORKDIR}"

echo "Aggregating qc_before/ reports with MultiQC..."
multiqc qc_before/ -n multiqc_raw_batch.html -o . --force

echo "OK: Combined report generated as multiqc_raw_batch.html"