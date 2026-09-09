#!/bin/bash
#SBATCH --job-name=qc_single
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=12G
#SBATCH --time=02:00:00
#SBATCH --output=%x_%j.out

###############################################################################
# SINGLE-SAMPLE WORKFLOW
#
# This version of the pipeline is designed for a single accession and
# serves as a simple introductory example before moving to automated
# batch processing workflows.
#
# You should create the script file using the specified filename
# and follow the setup instructions provided in the course materials.
###############################################################################

# ==============================================================================
# The Hardcoded Script (Single Run)
# ==============================================================================
# This is the safest script for testing a pipeline. 
# It runs one specific sample from start to finish.
# ==============================================================================

# Set up directory structure
cd /scratch/$(whoami)
mkdir -p hardcoded_single_run/{raw_reads,qc_before,trimmed_reads,qc_after}
cd hardcoded_single_run

# Hardcoded single run
SRR="SRR11092056"
echo "Processing ${SRR}..."

# 1. Retrieve the reads (Isolated Environment + Prefetch Fix)
module purge
module load SRA-Toolkit/3.2.0-gompi-2024a
prefetch ${SRR} -O raw_reads/
fasterq-dump raw_reads/${SRR} --split-files --outdir raw_reads --threads 4

# 2. Initial FastQC (Isolated Environment)
module purge
module load FastQC/0.12.1-Java-11
fastqc raw_reads/${SRR}_1.fastq raw_reads/${SRR}_2.fastq -o qc_before -t 4

# 3. Evidence-based trimming (Isolated Environment)
module purge
module load fastp/0.23.4-GCC-13.2.0
fastp -i raw_reads/${SRR}_1.fastq -I raw_reads/${SRR}_2.fastq -o trimmed_reads/${SRR}_1_clean.fastq -O trimmed_reads/${SRR}_2_clean.fastq --thread 4 --html ${SRR}_fastp.html

# 4. Post-trimming FastQC (Isolated Environment)
module purge
module load FastQC/0.12.1-Java-11
fastqc trimmed_reads/${SRR}_1_clean.fastq trimmed_reads/${SRR}_2_clean.fastq -o qc_after -t 4

# 5. MultiQC Reporting (Isolated Environment)
module purge
module load MultiQC/1.28-foss-2024a
multiqc qc_before/ -n multiqc_raw.html
multiqc qc_after/ -n multiqc_trimmed.html

echo "Done!"