# MIBO 8110: Applied Omics Projects 🧬

**Author:** Oliver Shey Njila 
**Environment:** Linux HPC (SLURM Workload Manager)

## 📌 Project Overview
This repository contains automated, reproducible, and scalable bioinformatics pipelines developed for the analysis of whole-genome sequencing (WGS) and metagenomic data. Designed for a high-performance computing (HPC) environment, these workflows process raw SRA reads through quality control, assembly, and functional annotation.

The initial pipeline scripts provided here focus on foundational data retrieval, quality control (QC), and progressive HPC scaling, serving as the required preprocessing engine for downstream metagenomic profiling and antimicrobial resistance (AMR) characterization.

---

## 📂 Repository Structure
```text
mibo8110-applied-omics-projects/
├── Initial_Pipeline_scripts/
│   ├── 1_the_hardcoded_script_single_run.sh           # Introductory single-sample workflow
│   ├── 2_Automated_Script_Processes_ALLRuns.sh        # Automated batch processing
│   ├── 3_Automated_Script_Hardcoded_20GB_Limit.sh     # Hardcoded directory size tracking
│   ├── 4_Automated_Script_User-DefinedMaximumSpace.sh # Configurable maximum storage limit
│   ├── 5_Automated_Script_System-BoundDynamicLimit.sh # GACRC dynamic filesystem monitoring
│   ├── 6_master_array_qc_ok.sh                        # High-performance SLURM Job Arrays
│   ├── 7_Final_master_array_qc.sh                     # Post-array MultiQC aggregation
│   ├── mod1_batch_qc_workflow.sh                      # Scalable batch execution script
│   └── mod1_qc_workflow.sh                            # Standard QC and trimming workflow
├── .gitignore                                         # HPC and genomic data ignore rules
└── README.md                                          # Project documentation
Phased Pipeline DevelopmentThis project incorporates progressive module development to teach and implement robust HPC data management. Key features include:Robust Data Retrieval: Utilizes a two-step prefetch and fasterq-dump method to safely download SRA data and successfully bypass NCBI connection timeout issues.  Batch Processing: Reads dynamic SRA accessions at runtime from a specifically formatted srr_list.txt file.  Dynamic Storage Guards: Queries the HPC file system directly to determine remaining free space and halts the processing loop before storage drops below a safety threshold, preventing file corruption and disk exhaustion.  High-Performance Parallelization: Implements SLURM Job Arrays where each task processes an individual sample independently to maximize computational efficiency.  Separation of Aggregation: Runs MultiQC as a separate, post-array script to guarantee all independent tasks have completed successfully before generating a consolidate project report.

Assembly & AMR Workflow Architecture
Assembly & AMR Workflow Architecture
The initial QC and preprocessing modules above lay the groundwork for the downstream end-to-end genomic analysis visualized below:

flowchart TD
    %% Input Data
    Input[srr_list.txt] --> Phase1

    %% Phase 1


    Out1 --> RScript[R Visualization & Analysis]



Usage1. Setup Input FilesThe batch automated pipelines require an input text file named exactly srr_list.txt located in the expected working directory.  The file must contain exactly one accession per line.  If you create or edit this file on a Windows machine, you must clean the hidden formatting by running the dos2unix srr_list.txt command before job submission to prevent loop and parsing issues.  

# Example list setup
echo -e "SRR11092056\nSRR11092057" > srr_list.txt
dos2unix srr_list.txt
2. Execution
Submit the desired script version to the SLURM scheduler:

# Run a simple single-sample test
sbatch Initial_Pipeline_scripts/mod1_qc_workflow.sh

# Run the parallelized Array workflow
sbatch --array=1-2 Initial_Pipeline_scripts/6_master_array_qc_ok.sh

# Run the post-array MultiQC aggregation
sbatch Initial_Pipeline_scripts/7_Final_master_array_qc.sh

3. Reviewing ResultsOnce processing is completed, interactive visual quality reports will be generated.  Use WinSCP (or a similar tool) to connect to the HPC system and download the HTML reports to your local computer.  Open the downloaded multiqc_report_raw.html and multiqc_report_trimmed.html files (or their batch equivalents) in a modern web browser.  Compare the reports to visually confirm that the fastp trimming step successfully removed low-quality tail regions and adapter-contaminated reads.  
