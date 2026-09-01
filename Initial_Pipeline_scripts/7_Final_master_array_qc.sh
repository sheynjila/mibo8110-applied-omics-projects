#!/bin/bash
#SBATCH --job-name=module1_post_array_summary
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --output=%x_%j.out

###############################################################################
# POST-ARRAY AGGREGATION: MULTIQC
###############################################################################
# In a SLURM Job Array workflow, each sample is processed by an independent
# array task running on the cluster. Because these tasks execute separately,
# it is not practical to generate a project-wide MultiQC report from within
# the array script itself.
#
# Instead, each job writes its own QC outputs (e.g., FastQC reports, fastp
# summaries) to the project directory. After all array tasks have completed,
# the results can be aggregated into a single report.
#
# WHY RUN MULTIQC SEPARATELY?
#   - Ensures all sample analyses have completed.
#   - Prevents partial or incomplete project summaries.
#   - Produces a consolidated view of quality metrics across all samples.
#   - Simplifies identification of failed runs or outlier samples.
#
# EDUCATIONAL NOTE:
#   This separation of per-sample processing (Job Array) and project-level
#   aggregation (MultiQC) reflects a common workflow pattern used in
#   production HPC bioinformatics pipelines.
###############################################################################

# 1. Navigate to the master array directory
cd /scratch/$(whoami)/master_array_qc

# 2. Run MultiQC on both pre- and post-trimming directories using modern modules
module purge
module load MultiQC/1.28-foss-2024a

echo "Generating aggregate MultiQC report for raw reads..."
multiqc qc_before/ -n multiqc_raw_array.html

echo "Generating aggregate MultiQC report for trimmed reads..."
multiqc qc_after/ -n multiqc_trimmed_array.html

###############################################################################
# DATA HOUSEKEEPING AND STORAGE TRIAGE
###############################################################################
# Before leaving the project directory, students should practice good data
# stewardship and storage management.
#
# Since GACRC scratch storage is subject to automatic purging after a defined
# retention period, important results should be preserved before cleanup.
#
# RECOMMENDED ACTIONS:
#   - Archive final reports.
#   - Save filtered FASTQ files required for downstream analyses.
#   - Remove temporary or intermediate files that can be regenerated.
#   - Verify that critical outputs have been copied to long-term storage.
###############################################################################

# 3. Preserve critical lightweight summary outputs into home directory
mkdir -p ~/capstone_results ~/capstone_scripts

echo "Archiving lightweight summaries and scripts to home directory..."
cp multiqc_raw_array.html ~/capstone_results/
cp multiqc_trimmed_array.html ~/capstone_results/
cp srr_list.txt ~/capstone_results/
cp 6_master_array_qc.sh ~/capstone_scripts/

# NOTE: Large raw and trimmed FASTQ files remain safely in /scratch for downstream assembly.

echo "Post-aggregation housekeeping complete!"

###############################################################################
# FINAL DELIVERABLE: MODULE 1 AUDIT LOG
###############################################################################
# The primary learning objective of Module 1 is not simply to run quality
# control software, but to evaluate whether the sequencing reads are suitable
# for downstream biological analysis and to justify any preprocessing steps
# that were performed.
#
# GUIDING QUESTION:
#   "Are the reads suitable for the intended biological analysis, and what
#   preprocessing is justified?"
#
# After MultiQC has generated the consolidated reports, students should transfer
# the HTML files to their local computer, open them in a web browser, and complete
# a short Audit Log summarizing their findings.
#
# ---------------------------------------------------------------------------
# AUDIT LOG QUESTIONS FOR STUDENTS
# ---------------------------------------------------------------------------
#
# 1. Quality Assessment
#    - Did the average Phred quality scores decrease toward the 3' ends of
#      the reads?
#    - Is this pattern expected for Illumina sequencing data?
#
# 2. Adapter Content and Filtering
#    - Did fastp successfully remove adapter contamination?
#    - What percentage of reads survived the filtering process?
#    - Was any substantial loss of sequencing data observed?
#
# 3. Biological Verdict
#    - Based on the quality control results, are the reads suitable for
#      downstream genome assembly?
#    - Is the dataset of sufficient quality to support assembly of a
#      Salmonella genome with confidence?
#
# ---------------------------------------------------------------------------
# TRANSITION TO MODULE 2: PATHOGEN GENOME ASSEMBLY
# ---------------------------------------------------------------------------
# Primary Tool: SPAdes (or Unicycler)
# Core Objective: Reconstruct the pathogen genome by identifying overlaps among
# sequencing reads and assembling them into longer contiguous sequences (contigs).
#
# Key Questions for Module 2:
#   - How can millions of short sequencing reads be combined into a genome?
#   - What is a contig, scaffold, and assembly graph?
#   - How do we evaluate assembly quality?
#   - Is the assembled genome suitable for downstream analyses such as
#     antimicrobial resistance detection, serotyping, and comparative genomics?
###############################################################################