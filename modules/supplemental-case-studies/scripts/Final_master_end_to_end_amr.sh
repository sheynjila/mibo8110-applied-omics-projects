###############################################################################
# WORKFLOW OVERVIEW
#
# Final production version of the end-to-end assembly and AMR pipeline.
#
# Optimized for bacterial isolate whole-genome sequencing (WGS) data
# and configured with increased HPC resources to support memory-intensive
# genome assembly workloads.
#
# Key updates:
#   - Increased SLURM resources for assembly-scale processing
#   - Extended runtime for long-running jobs
#   - SPAdes configured with isolate mode for high-coverage isolates
#   - Designed to support assembly, quality assessment, and AMR profiling
#
# This configuration is intended for production-scale isolate genome
# analysis rather than lightweight alignment-based workflows.
#
# ------------------------------------------------------------------------
# CORRECTED VERSION — see pipeline_appraisal_2026-08-31.md for the full
# review. Fixes applied:
#   1. The per-sample loop is now fault-tolerant: a failed download, assembly,
#      QUAST run, or AMR scan is logged and the pipeline moves to the next
#      sample instead of calling `exit 1` and killing the whole 16-CPU/64GB
#      job (which would previously discard all already-completed samples'
#      work too).
#   2. AMR_ORGANISM is pulled into a variable at the top instead of being
#      hardcoded inline, so reusing this script for a non-Salmonella dataset
#      is a one-line change instead of an edit buried in Phase 2.
#   3. Aggregation (Phase 3) now only includes samples that actually produced
#      an AMR report, and a summary of any failed samples is printed at the
#      end.
###############################################################################


#!/bin/bash
#SBATCH --job-name=end_to_end_amr
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16       # SPAdes requires heavy multi-threading
#SBATCH --mem=64G                # Assembly is extremely memory-intensive
#SBATCH --time=48:00:00          # Assembly takes significantly longer than alignment
#SBATCH --output=master_amr_%j.out
#SBATCH --error=master_amr_%j.err

# Enable strict error handling
set -e
set -o pipefail

# FIX #2: organism now configurable in one place.
AMR_ORGANISM="${AMR_ORGANISM:-Salmonella}"

echo "=========================================================="
echo " INITIATING END-TO-END ASSEMBLY & AMR PIPELINE (organism=${AMR_ORGANISM}) "
echo "=========================================================="

# ==========================================
# PHASE 1: DIRECTORY SETUP & AMR DB UPDATE
# ==========================================
WORKDIR="/scratch/$(whoami)/end_to_end_amr_pipeline"
mkdir -p ${WORKDIR}/{raw_reads,qc,trimmed_reads,spades_output,final_assemblies,quast_qc,amr_reports,logs}
cd ${WORKDIR}

if [ ! -f "srr_list.txt" ]; then
    echo "WARNING: srr_list.txt not found. Generating default test list..."
    echo -e "SRR11092056\nSRR11092057" > srr_list.txt
fi
sed -i 's/\r$//' srr_list.txt

echo -e "\n[PHASE 1] Updating AMRFinder Database..."
module purge
module load ncbi-amrfinderplus/3.11.11-foss-2022a
amrfinder -u

# ==========================================
# PHASE 2: PROCESSING, ASSEMBLY & ANNOTATION LOOP
# ==========================================
echo -e "\n[PHASE 2] EXECUTING ASSEMBLY & AMR SCANNING..."

# FIX #1: track failures instead of a hard abort on the first bad sample.
FAILED_SAMPLES=()

for SRR in $(cat srr_list.txt); do
    echo "----------------------------------------------------------"
    echo " PROCESSING SAMPLE: ${SRR}"
    echo "----------------------------------------------------------"

    set +e
    (
        set -e
        set -o pipefail

        # Step A: Download
        echo ">> Downloading from SRA..."
        module purge
        module load SRA-Toolkit/3.2.0-gompi-2024a
        prefetch ${SRR} -O raw_reads/
        fasterq-dump raw_reads/${SRR} --split-files --outdir raw_reads --threads 4

        # Step B: Quality Control & Trimming
        echo ">> Trimming with fastp..."
        module purge
        module load fastp/0.23.4-GCC-13.2.0
        fastp -i raw_reads/${SRR}_1.fastq -I raw_reads/${SRR}_2.fastq \
              -o trimmed_reads/${SRR}_1_clean.fastq -O trimmed_reads/${SRR}_2_clean.fastq \
              --thread 8 --html qc/${SRR}_fastp.html 2> qc/${SRR}_fastp.log

        # Step C: De Novo Assembly (SPAdes)
        echo ">> Running SPAdes Assembler (This may take several hours)..."
        module purge
        module load SPAdes/3.15.5-GCC-13.2.0

        spades.py -1 trimmed_reads/${SRR}_1_clean.fastq \
                  -2 trimmed_reads/${SRR}_2_clean.fastq \
                  -o spades_output/${SRR} \
                  --isolate \
                  -t 16 -m 64

        if [ -f "spades_output/${SRR}/contigs.fasta" ]; then
            cp spades_output/${SRR}/contigs.fasta final_assemblies/${SRR}_contigs.fasta
        else
            echo ">> ERROR: Assembly failed for ${SRR} (no contigs.fasta produced)."
            exit 1
        fi

        # Step D: Quality Assessment (QUAST)
        echo ">> Assessing assembly quality with QUAST..."
        module purge
        module load QUAST/5.2.0-foss-2022a
        quast.py final_assemblies/${SRR}_contigs.fasta \
                 -o quast_qc/${SRR} \
                 --threads 8

        # Step E: AMR Annotation (AMRFinderPlus)
        echo ">> Scanning assembled contigs for AMR genes..."
        module purge
        module load ncbi-amrfinderplus/3.11.11-foss-2022a
        amrfinder -n final_assemblies/${SRR}_contigs.fasta \
                  -O ${AMR_ORGANISM} \
                  --plus \
                  --threads 16 \
                  > amr_reports/${SRR}_amr_results.tsv
    )
    SAMPLE_STATUS=$?
    set -e

    if [ ${SAMPLE_STATUS} -ne 0 ]; then
        echo ">> ERROR: Sample ${SRR} failed (exit code ${SAMPLE_STATUS}) - skipping to next sample."
        FAILED_SAMPLES+=("${SRR}")
        continue
    fi

    echo ">> Completed ${SRR}"
done

# ==========================================
# PHASE 3: AGGREGATE RESULTS FOR R & MULTIQC
# ==========================================
echo -e "\n[PHASE 3] AGGREGATING FINAL REPORTS..."

# FIX #3: only aggregate samples that actually produced an AMR report.
cd amr_reports
AMR_TSVS=$(ls *_amr_results.tsv 2>/dev/null || true)

if [ -z "${AMR_TSVS}" ]; then
    echo "WARNING: No AMR results were produced - skipping aggregation."
else
    echo -n -e "Sample\t" > ALL_SAMPLES_AMR_SUMMARY.tsv
    head -n 1 $(echo "${AMR_TSVS}" | head -n 1) >> ALL_SAMPLES_AMR_SUMMARY.tsv

    for TSV in ${AMR_TSVS}; do
        BASENAME=$(basename ${TSV} _amr_results.tsv)
        awk -v sample="${BASENAME}" 'NR>1 {print sample "\t" $0}' ${TSV} >> ALL_SAMPLES_AMR_SUMMARY.tsv
    done
fi
cd ${WORKDIR}

# 2. Compile standard MultiQC report
module purge
module load MultiQC/1.28-foss-2024a
multiqc qc/ quast_qc/ -n final_assembly_report.html

echo "=========================================================="
if [ ${#FAILED_SAMPLES[@]} -gt 0 ]; then
    echo " PIPELINE COMPLETED WITH ${#FAILED_SAMPLES[@]} FAILED SAMPLE(S):"
    printf '   %s\n' "${FAILED_SAMPLES[@]}"
    echo " Re-run those specific SRR IDs individually if you need complete data."
else
    echo " PIPELINE SUCCESSFULLY COMPLETED! "
fi
echo " Assembly QA: final_assembly_report.html"
echo " R-ready AMR Data: amr_reports/ALL_SAMPLES_AMR_SUMMARY.tsv"
echo "=========================================================="

###############################################################################
# CAPSTONE METHODOLOGY SUMMARY (carried over)
#
# 1. De Novo Assembly Advantages
# ---------------------------------------------------------------------------
# De novo assembly identifies shared sequence patterns among reads to build
# long contigs without a known reference sequence. This is especially
# useful for discovering genomic content absent from references, such as
# novel antimicrobial-resistance plasmids and mobile genetic elements.
#
# 2. Interpreting QUAST Quality Metrics
# ---------------------------------------------------------------------------
# Total Length: For Salmonella, expect roughly 4.8-5.0 Mb. If the length is
# around 9.5 Mb, this strongly indicates contamination or a mixed-species
# sample. N50: a higher N50 indicates a more contiguous assembly.
#
# 3. AMR Guardrails
# ---------------------------------------------------------------------------
# Assembly alone does not prove that a contig is a plasmid or carries an
# AMR determinant — the downstream AMRFinderPlus module provides the
# functional annotation required to make those conclusions.
###############################################################################
