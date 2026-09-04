# 🧬 Automated RNA-Seq Processing Pipeline (HPC/Slurm)

**Project:** PRJNA229998 Airway Transcriptomics  
**Environment:** Slurm High-Performance Computing (HPC)  
**Input:** Raw SRA Accessions  
**Output:** Gene expression count matrices ready for DESeq2  

## 📌 Architecture Overview

This repository contains a robust, end-to-end RNA-sequencing pipeline optimized for Slurm-managed HPC clusters. It features a two-tiered architecture designed for both development safety and production speed:

1. **Modular Suite (Scripts 1–7, 9):** A step-by-step pipeline for isolation, testing, debugging, and validation.
2. **V9 Monolithic Suite (Script 8):** A production-grade, automated bash loop (`8_Automated_Pipeline_AllSamples.sh`) that processes an entire sample cohort continuously.

## 🛠️ Crucial HPC Safeguards (The "V9" Fixes)

During validation on the cluster, three critical stability patches were implemented globally to prevent silent failures and memory crashes:

* **Memory Allocation Limits:** STAR indexing requires massive RAM footprints. All alignment and indexing scripts explicitly request `#SBATCH --mem=36G` to prevent `std::bad_alloc` crashes.
* **Read-Pair Synchronization Parity:** Fastq dumping utilizes `fasterq-dump --split-3` (replacing `--split-files`). This ensures that dropped or empty reads do not break Read 1 / Read 2 parity, preventing catastrophic downstream STAR mapping failures.
* **Environment Isolation:** Every script issues `module purge` prior to loading toolchains to prevent hidden cluster library conflicts.
* **Dynamic Rsubread Bootstrapping:** Bypasses administrative cluster blocks by automatically injecting and compiling Bioconductor/Rsubread into a local `R_libs` scratch directory.

## 📁 Directory Structure & Data Flow

```text
/scratch/$(whoami)/
  ├── PRJNA229998_airway_pipeline_ref/        # Shared STAR Genome Index & Ensembl GTF
  ├── PRJNA229998_airway_pipeline_stepbystep/ # Modular validation workspace
  └── PRJNA229998_airway_pipeline_automated/  # V9 Cohort production workspace
      ├── raw_reads/                          # Fetched via SRA Toolkit
      ├── trimmed_reads/                      # Cleaned via fastp
      ├── alignments/                         # BAM files mapped via STAR
      ├── counts/                             # Final matrices via featureCounts
      └── R_libs/                             # Transient user-space R packages



Execution Procedures
A. Running the Automated Cohort (Production)
For processing the complete multi-sample cohort without intervention:

Bash
sbatch 8_Automated_Pipeline_AllSamples.sh
B. Running Modular Validation (Debugging/Testing)
Execute in sequence, ensuring the previous Slurm job finishes before submitting the next:

sbatch 1_Fetch_SRA.sh

sbatch 2_FastQC_Raw.sh

sbatch 3_Trim_fastp.sh

sbatch 4_FastQC_Clean.sh

sbatch 5_STAR_Index.sh (Run once per reference genome)

sbatch 6_STAR_Align.sh

sbatch 7_Quantification_Rsubread.sh

bash 9_QC_Report_Generator.sh (Generates markdown validation report)

🔎 Monitoring & Log Tracking
Check active jobs: squeue -u $(whoami)

Cancel a running job: scancel <JOB_ID>

Stream active logs: tail -f $(ls -t ~/slurm-*.out | head -n 1)

Audit for errors: grep -iE "error|fatal|cannot|execution halted" ~/slurm-*.out

Verify mapping stats: cat counts/*_counts.txt.summary

📊 Downstream Handoff
The ultimate output is the raw count matrix located at counts/airway_raw_counts.txt. This tab-delimited file is structured with Ensembl Gene IDs as rows and Sample IDs as columns, completely formatted and ready for direct import into DESeq2 in R for differential gene expression modeling.


<FollowUp label="Commit README" query="Now that you have the raw text, are you ready to i