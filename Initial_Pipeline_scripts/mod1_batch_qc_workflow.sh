#!/bin/bash
#SBATCH --job-name=batch_sra_qc
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=10:00:00 
#SBATCH --output=%x_%j.out

echo "=========================================================="
echo " INITIATING AUTOMATED BATCH PIPELINE"
echo "=========================================================="

# Set up directory structure
cd /scratch/$(whoami)
mkdir -p capstone_automated/module1/{raw_reads,qc_before,trimmed_reads,qc_after}
cd capstone_automated/module1

# Note: Accessions are no longer hardcoded here. 
# The script will read them dynamically from srr_list.txt in the loop below.

# THE LOOP: Read each line from srr_list.txt
for SRR in $(cat srr_list.txt); do
    
    echo "========================================"
    echo "STARTING PIPELINE FOR: ${SRR}"
    echo "========================================"
    
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
    fastp -i raw_reads/${SRR}_1.fastq -I raw_reads/${SRR}_2.fastq -o trimmed_reads/${SRR}_1_clean.fastq -O trimmed_reads/${SRR}_2_clean.fastq --thread 4 --html ${SRR}_fastp_report.html

    # 4. Post-trimming FastQC (Isolated Environment)
    module purge
    module load FastQC/0.12.1-Java-11
    fastqc trimmed_reads/${SRR}_1_clean.fastq trimmed_reads/${SRR}_2_clean.fastq -o qc_after -t 4
    
    echo "FINISHED: ${SRR}"
    
done

# 5. Run MultiQC ONCE at the very end to summarize ALL samples together
echo "========================================"
echo "Compiling MultiQC reports..."
module purge
module load MultiQC/1.28-foss-2024a
multiqc qc_before/ -n multiqc_report_raw_all_samples.html
multiqc qc_after/ -n multiqc_report_trimmed_all_samples.html

echo "Batch workflow complete!"

###############################################################################
# BEFORE SUBMITTING
#
# This pipeline reads accession IDs from srr_list.txt. Ensure the file
# has Unix line endings before running to avoid loop and parsing issues.
#
# Convert the file if it was created or edited on Windows:
#
#   dos2unix srr_list.txt
#
# After conversion, submit the job with sbatch.
#
# This workflow is designed for scalable batch processing, allowing
# multiple SRA accessions to be analyzed from a single input list.
###############################################################################

###############################################################################
# SRR LIST SETUP
#
# This pipeline expects srr_list.txt to exist in the working directory
# on scratch storage. If the file is missing, the accession loop cannot
# start and the job will fail.
#
# Create srr_list.txt in the module1 directory and place one accession
# per line. Example:
#
#   SRR11092056
#   SRR11092057
#
# After creating the file, convert it to Unix format if it was edited
# on Windows to avoid line-ending issues during processing.
###############################################################################

###############################################################################
# FINAL PREPARATION
#
# Add one accession per line in srr_list.txt. Example:
#
#   SRR11092056
#   SRR11092057
#
# After saving the file, convert it to Unix format if needed to prevent
# line-ending issues during the accession-processing loop.
#           dos2unix srr_list.txt
# Once the file has been validated, submit the workflow to SLURM for
# automated batch processing.
#      sbatch ~/batch_qc_workflow.sh
###############################################################################

