#!/bin/bash
#SBATCH --job-name=qc04_fastp
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=01:00:00

###############################################################################
# 4_Evidence_Based_Fastp_Trimming.sh
#
# NARRATIVE:
# This script executes read trimming as an evidence-based decision rather than 
# a blind default[cite: 8]. It reads the FastQC `summary.txt` generated in Script 2[cite: 8]. 
# If (and only if) adapter contamination or 3' quality decay is flagged, it 
# dynamically appends the appropriate `fastp` parameters[cite: 8]. It then documents 
# this programmatic logic into a plain-text rationale file to ensure strict 
# methodological reproducibility[cite: 8].
###############################################################################
set -e
set -o pipefail

# V9 FIX: Purge stale modules before loading toolchains
module purge
module load FastQC/0.11.9-Java-11
module load fastp/0.23.4-GCC-13.2.0

WORKDIR="/scratch/$(whoami)/PRJNA229998_airway_pipeline_stepbystep"
mkdir -p "${WORKDIR}"/{trimmed_reads,fastp_reports,qc_after}
cd "${WORKDIR}"

SRR="SRR1039508"
SUMMARY_1="qc_before/${SRR}_1_fastqc/summary.txt"
SUMMARY_2="qc_before/${SRR}_2_fastqc/summary.txt"

FASTP_FLAGS=("--length_required" "36")
RATIONALE_LINES=("- Minimum length floor of 36bp applied unconditionally.")

if grep -qE "^(WARN|FAIL)\s+Adapter Content" "${SUMMARY_1}" "${SUMMARY_2}"; then
    FASTP_FLAGS+=("--detect_adapter_for_pe")
    RATIONALE_LINES+=("- Adapter trimming ENABLED: FastQC flagged Adapter Content.")
fi

if grep -qE "^(WARN|FAIL)\s+Per base sequence quality" "${SUMMARY_1}" "${SUMMARY_2}"; then
    FASTP_FLAGS+=("--cut_tail" "--cut_tail_mean_quality" "20")
    RATIONALE_LINES+=("- 3' quality trimming ENABLED: FastQC flagged quality decay.")
fi

echo "Running fastp with flags: ${FASTP_FLAGS[*]}"
fastp -i "raw_reads/${SRR}_1.fastq" -I "raw_reads/${SRR}_2.fastq" \
      -o "trimmed_reads/${SRR}_1_clean.fastq" -O "trimmed_reads/${SRR}_2_clean.fastq" \
      "${FASTP_FLAGS[@]}" --thread 4 \
      --json "fastp_reports/${SRR}_fastp.json" --html "fastp_reports/${SRR}_fastp.html"

{
    echo "# Trimming rationale for ${SRR}"
    printf '%s\n' "${RATIONALE_LINES[@]}"
} > "fastp_reports/${SRR}_THRESHOLD_RATIONALE.md"

fastqc "trimmed_reads/${SRR}_1_clean.fastq" "trimmed_reads/${SRR}_2_clean.fastq" -o qc_after -t 4 --extract