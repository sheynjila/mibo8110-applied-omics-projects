#!/bin/bash
#SBATCH --job-name=qc_all
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=12G
#SBATCH --time=10:00:00
#SBATCH --output=%x_%j.out

# ==============================================================================
# BATCH PROCESSING CONFIGURATION
# ==============================================================================
# This script uses a text file named srr_list.txt (which you must create 
# in the same directory, containing one SRR per line) and processes 
# every single run in that file.
# ==============================================================================

# Set up directory structure
cd /scratch/$(whoami)
mkdir -p automated_allruns/{raw_reads,qc_before,trimmed_reads,qc_after}
cd automated_allruns

# Loop through every run in the list
for SRR in $(cat srr_list.txt); do
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
    
done

# Run MultiQC once at the end for all samples
module purge
module load MultiQC/1.28-foss-2024a
multiqc qc_before/ -n multiqc_raw_all.html
multiqc qc_after/ -n multiqc_trimmed_all.html

echo "All samples processed!"


###############################################################################
# INPUT FILE SETUP
#
# Create srr_list.txt in the directory expected by the workflow.
# The pipeline reads accession IDs from this file at runtime and will
# fail if the file is missing or located elsewhere.
#
# Add one accession per line and verify the file exists before
# submitting the job.
###############################################################################
# 1. Create the new working directory
#       mkdir -p /scratch/$(whoami)/automated_allruns

# 2. Move into that directory
#       cd /scratch/$(whoami)/automated_allruns

# 3. Create and edit the list
#       nano srr_list.txt

# Type in your accessions, like SRR11092056 and SRR11092057, one on each line, then save and exit

# 4. Clean the file of hidden formatting
#       dos2unix srr_list.txt

# Step 3: Resubmit!
# Now that the script is clean and the list exists in the correct folder, you can go back to your home directory and submit it safely:

#           cd ~
#            sbatch 2_Automated_Script_Processes_ALLRuns.sh

###############################################################################

###############################################################################
# INPUT FILE REQUIREMENTS
#
# The accession list must be named exactly:
#
#   srr_list.txt
#
# The script searches for this file in the configured working directory.
# If the filename or location is incorrect, the accession-processing
# loop will not start.
#
# Before submitting the job:
#   - Verify the file exists in the expected directory.
#   - Confirm the filename matches exactly (case-sensitive).
#   - Ensure the file contains one accession per line.
#
# Keeping the input file in the correct location allows the workflow
# to automatically locate and process all listed accessions.
###############################################################################
