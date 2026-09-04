#!/bin/bash
#SBATCH --job-name=qc02_fastqc
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=6G
#SBATCH --time=00:45:00

###############################################################################
# 2_Interpret_FastQC_Diagnostics.sh
# 
# NARRATIVE:
# This script performs the initial quality control assessment on raw FASTQ files.
# Beyond simply generating HTML reports, it uses the `--extract` flag to expose 
# the underlying `summary.txt` data. It then programmatically parses this text 
# to synthesize the PASS/WARN/FAIL metrics. This extraction is critical, as it 
# converts visual QC plots into machine-readable text that will drive automated, 
# evidence-based trimming decisions in Script 4.
###############################################################################
set -e
set -o pipefail

module load FastQC/0.11.9-Java-11

WORKDIR="/scratch/$(whoami)/PRJNA229998_airway_pipeline_stepbystep"
mkdir -p "${WORKDIR}/qc_before"
cd "${WORKDIR}"

SRR="SRR1039508"

echo "Running FastQC on ${SRR}..."
fastqc "raw_reads/${SRR}_1.fastq" "raw_reads/${SRR}_2.fastq" -o qc_before -t 4 --extract

for MATE in 1 2; do
    SUMMARY="qc_before/${SRR}_${MATE}_fastqc/summary.txt"
    echo "---- ${SRR} mate ${MATE}: FastQC verdicts ----"
    cat "${SUMMARY}"
done

FAIL_COUNT=$(cat qc_before/${SRR}_1_fastqc/summary.txt qc_before/${SRR}_2_fastqc/summary.txt | grep -c "^FAIL" || true)
WARN_COUNT=$(cat qc_before/${SRR}_1_fastqc/summary.txt qc_before/${SRR}_2_fastqc/summary.txt | grep -c "^WARN" || true)

echo "Synthesis: ${FAIL_COUNT} FAIL, ${WARN_COUNT} WARN (across both mates)."