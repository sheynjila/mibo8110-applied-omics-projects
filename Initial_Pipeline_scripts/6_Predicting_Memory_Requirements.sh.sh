#!/bin/bash
#SBATCH --job-name=master_qc_array
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=12G
#SBATCH --time=02:00:00
#SBATCH --array=1-2     # Adjust array bounds to match srr_list.txt line count
#SBATCH --output=logs/qc_array_%A_%a.out

###############################################################################
# PART 1: HIGH-PERFORMANCE PARALLELIZATION
###############################################################################
# This module introduces HPC parallel processing with SLURM Job Arrays.
#
# In a SLURM Job Array workflow, each task processes an individual sample
# independently and generates its own quality-control outputs.
###############################################################################

# Setup Directories
cd /scratch/$(whoami)
mkdir -p master_array_qc/{raw_reads,qc_before,trimmed_reads,qc_after,logs}
cd master_array_qc

# Array Logic: Grab the exact SRR for this specific task ID
# This reads line 'n' from srr_list.txt where 'n' is the array ID
SRR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" srr_list.txt)

if [ -z "$SRR" ]; then
    echo "No SRR found for ID $SLURM_ARRAY_TASK_ID"
    exit 0
fi

echo "Array Task ${SLURM_ARRAY_TASK_ID}: Processing ${SRR}..."

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

echo "Task ${SLURM_ARRAY_TASK_ID} Complete!"

###############################################################################
# PART 2: MULTIQC REPORT GENERATION (POST-ARRAY)
###############################################################################
# Because MultiQC aggregates results across all samples, it is typically run
# only after every array task has completed successfully.
#
# Once the entire job array has finished, you must generate the consolidated
# QC report manually by running:
#
#     module load MultiQC/1.28-foss-2024a
#     multiqc qc_before/ -n multiqc_raw_array.html
#     multiqc qc_after/ -n multiqc_trimmed_array.html
###############################################################################

###############################################################################
# PART 3: ADVANCED VARIANTS - LONG-READ SEQUENCING (ONT)
###############################################################################
# Previous examples assume Illumina short-read sequencing data. While highly
# accurate, short reads may struggle to fully resolve complex genomic regions,
# including mobile genetic elements and antimicrobial resistance (AMR)
# plasmids.
#
# Oxford Nanopore Technologies (ONT) produces long reads that can span
# thousands of bases, making them valuable for genome assembly and structural
# analysis.
#
# Key Differences:
#   - Illumina: shorter reads, lower error rates.
#   - ONT: longer reads, higher error rates.
#
# Workflow Modifications for ONT:
#   - FastQC   --> NanoPlot
#   - fastp    --> Filtlong
###############################################################################