#!/bin/bash
#SBATCH --job-name=qc_20G_limit
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=12G
#SBATCH --time=10:00:00
#SBATCH --output=%x_%j.out

# ==============================================================================
# SAFETY CHECK: AUTOMATED STORAGE GUARD
# ==============================================================================
# This script checks the total size of your working directory before 
# downloading the next file. 
#
# If the folder has reached or exceeded 20GB, it breaks out of the loop, 
# skips the rest of the list, and generates the final reports for whatever 
# it successfully processed.
# ==============================================================================

# Set up directory structure
cd /scratch/$(whoami)
mkdir -p automated_hardcoded_20GB_limit/{raw_reads,qc_before,trimmed_reads,qc_after}
cd automated_hardcoded_20GB_limit

LIMIT_GB=20

for SRR in $(cat srr_list.txt); do
    # Check current directory size in Gigabytes
    CURRENT_GB=$(du -s -BG . | awk '{print $1}' | tr -d 'G')
    
    # Stop if we hit or exceed the limit
    if [ "$CURRENT_GB" -ge "$LIMIT_GB" ]; then
        echo "Storage reached ${CURRENT_GB}GB (Limit is ${LIMIT_GB}GB). Stopping downloads."
        break
    fi

    echo "Current space: ${CURRENT_GB}GB. Processing ${SRR}..."
    
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
    
done

# Run MultiQC once at the end for all samples
module purge
module load MultiQC/1.28-foss-2024a
multiqc qc_before/ -n multiqc_raw_limited.html
multiqc qc_after/ -n multiqc_trimmed_limited.html
echo "Limited run complete!"

###############################################################################
# INPUT FILE REQUIREMENT
#
# The workflow expects srr_list.txt to exist in the configured working
# directory before execution. If the file is missing, the accession loop
# cannot start and the job will exit immediately.
#
# Verify:
#   - The file is named exactly: srr_list.txt
#   - The file is located in the expected working directory
#   - One accession is listed per line
#
# Test accessions may be used to confirm that the workflow, directory
# structure, and batch-processing logic are functioning correctly before
# running larger datasets.

#mkdir -p /scratch/$(whoami)/automated_hardcoded_20GB_limit
#  cd /scratch/$(whoami)/automated_hardcoded_20GB_limit
#  echo -e "SRR11092056\nSRR11092057" > srr_list.txt
#  dos2unix srr_list.txt
#  cd ~
#  mkdir -p /scratch/$(whoami)/automated_hardcoded_20GB_limit
#  cd /scratch/$(whoami)/automated_hardcoded_20GB_limit
#  echo -e "SRR11092056\nSRR11092057" > srr_list.txt
#  dos2unix srr_list.txt
#  cd ~
#  sbatch 3_Automated_Script_Hardcoded_20GB_Limit.sh
###############################################################################