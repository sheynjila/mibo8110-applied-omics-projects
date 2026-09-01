#!/bin/bash
#SBATCH --job-name=sra_qc_fastp
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=12G
#SBATCH --time=02:00:00
#SBATCH --output=%x_%j.out

echo "=========================================================="
echo " INITIATING MODULE 1: PATHOGEN QC & PREPROCESSING "
echo "=========================================================="

# Set up directory structure
cd /scratch/$(whoami)
mkdir -p capstone/module1/{raw_reads,qc_before,trimmed_reads,qc_after}
cd capstone/module1
ACCESSION="SRR11092056"

# Step 1: Retrieve the reads (Isolated Environment)
echo ">> Fetching raw reads for ${ACCESSION}..."
module purge
module load SRA-Toolkit/3.2.0-gompi-2024a

# First, safely download the compressed data
prefetch ${ACCESSION} -O raw_reads/

# Then, extract the fastq files from the downloaded data
fasterq-dump raw_reads/${ACCESSION} --split-files --outdir raw_reads --threads 4

# Step 2: Inspect raw reads (Isolated Environment)
echo ">> Running initial FastQC..."
module purge
module load FastQC/0.12.1-Java-11
module load MultiQC/1.28-foss-2024a
fastqc raw_reads/${ACCESSION}_1.fastq raw_reads/${ACCESSION}_2.fastq -o qc_before -t 4
multiqc qc_before/ -n multiqc_report_raw.html

# Step 3: Evidence-based trimming (Isolated Environment)
echo ">> Trimming reads with fastp..."
module purge
module load fastp/0.23.4-GCC-13.2.0
fastp -i raw_reads/${ACCESSION}_1.fastq -I raw_reads/${ACCESSION}_2.fastq -o trimmed_reads/${ACCESSION}_1_clean.fastq -O trimmed_reads/${ACCESSION}_2_clean.fastq --thread 4 --html fastp_report.html

# Step 4: Repeat QC to compare (Isolated Environment)
echo ">> Running post-trimming FastQC..."
module purge
module load FastQC/0.12.1-Java-11
module load MultiQC/1.28-foss-2024a
fastqc trimmed_reads/${ACCESSION}_1_clean.fastq trimmed_reads/${ACCESSION}_2_clean.fastq -o qc_after -t 4
multiqc qc_after/ -n multiqc_report_trimmed.html

echo "=========================================================="
echo " Workflow Complete! "
echo "=========================================================="


###############################################################################
# SUCCESS VALIDATION NOTES
#
# You have successfully built a fully functioning, end-to-end HPC pipeline.
#
# Why this run is considered successful:
#
# GOLDEN METRIC: Job Runtime = 00:04:52 (Job ID: 29834)
#
# Previous Job (29833):
#   - Failed after approximately 2 minutes.
#   - NCBI connection timed out during data retrieval.
#
# Current Job (29834):
#   - Completed in approximately 5 minutes.
#   - Runtime is consistent with the expected workflow:
#       1. prefetch         -> Download SRA data safely from NCBI
#       2. fasterq-dump     -> Convert SRA archive to FASTQ files
#       3. FastQC           -> Generate pre-trim quality report
#       4. fastp            -> Trim adapters and low-quality reads
#       5. FastQC           -> Generate post-trim quality report
#       6. MultiQC          -> Aggregate all QC reports
#
# Job Status:
#   State     = COMPLETED
#   ExitCode  = 0:0
#
# Interpretation:
#   - The pipeline executed from start to finish without fatal errors.
#   - The two-step prefetch + fasterq-dump approach successfully bypassed
#     the NCBI timeout issues encountered previously.
#   - The workflow is stable and suitable for production-scale testing
#     and student training exercises.
#
# POST-RUN VALIDATION
#
# Verify that output files were created successfully:
#
#   ls -lh *.fastq*
#
# Review generated reports:
#
#   ls -lh *_fastqc.html
#   ls -lh multiqc_report.html
#
# Check job accounting information:
#
#   sacct -j 29834 --format=JobID,State,Elapsed,ExitCode
#
# Recommended next step:
#   Open the FastQC and MultiQC reports to confirm sequence quality,
#   trimming performance, and overall dataset health before proceeding
#   to downstream analyses.
###############################################################################
#   ls -lh /scratch/osnjila/capstone/module1/raw_reads/
#   ls -lh /scratch/osnjila/capstone/module1/trimmed_reads/

# Look at the final, clean log:
#   cat sra_qc_fastp_29834.out



###############################################################################
# DOWNLOADING AND REVIEWING QUALITY CONTROL REPORTS
#
# This is one of the most rewarding steps in the bioinformatics workflow.
# The pipeline has generated interactive visual quality reports that allow
# you to inspect sequencing quality before and after read trimming.
#
# VIEWING THE REPORTS
#
# 1. Open WinSCP and connect to the HPC system.
#
# 2. Navigate to the project directory:
#
#      /scratch/osnjila/capstone/module1/
#
# 3. Locate the quality control output folders:
#
#      qc_before/
#      qc_after/
#
# 4. Download the generated HTML reports to your local computer.
#    Common report files include:
#
#      multiqc_report_raw.html
#      multiqc_report_trimmed.html
#
# 5. Open the downloaded HTML files in a modern web browser such as:
#
#      - Google Chrome
#      - Microsoft Edge
#
#
# INTERPRETING THE RESULTS
#
# Compare the reports generated before and after trimming.
#
# Focus on the following MultiQC sections:
#
#      - Per Base Sequence Quality
#      - Sequence Quality Histograms
#      - Adapter Content
#      - Sequence Length Distribution
#
# Expected Outcome:
#
#      The trimmed dataset should show improved overall quality scores
#      compared to the raw dataset.
#
# Typical improvements include:
#
#      - Higher average Phred quality scores
#      - Reduced low-quality tail regions
#      - Fewer adapter-contaminated reads
#      - Cleaner sequence quality distributions
#
# These improvements confirm that the fastp trimming step successfully
# removed low-quality bases and sequencing artifacts.
#
#
# TEACHING NOTE
#
# For student demonstrations, the "Sequence Quality Histograms" and
# "Per Base Sequence Quality" plots usually provide the clearest visual
# example of how raw sequencing data is improved through quality control.
#
# Encourage students to compare:
#
#      Raw Reads  --->  Trimmed Reads
#
# and identify how the quality metrics change throughout the workflow.
#
###############################################################################

###############################################################################
# SLURM LOG FILE LOCATION
#
# SLURM creates the .out log file in the directory where the sbatch
# command was submitted, not necessarily where the script later changes
# directories during execution.
#
# Common locations to check:
#
#   Home directory:
#       cd ~
#       ls -l *.out
#
#   Testing directory:
#       cd ~/testing
#      ls -l *.out
#
# If the log file cannot be found, search your account:
#
#       find ~ -name "sra_qc_fastp_29834.out"
#
# To view the contents of the log file:
#
#       cat sra_qc_fastp_29834.out
#
###############################################################################