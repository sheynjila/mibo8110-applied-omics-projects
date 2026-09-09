#!/bin/bash
#SBATCH --job-name=qc06_pipeline
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=03:00:00

###############################################################################
# 6_Automated_QC_Pipeline_AllSamples.sh
# Module 1 (Raw Sequencing Reads & QC) — Script 6 of 7
#
# ONE NEW IDEA on top of scripts 1-5: put every habit from this module —
# verified retrieval, diagnostic FastQC, evidence-based fastp trimming, and
# paired-end integrity checking — into ONE fault-tolerant loop that runs
# unattended over every sample in the study. This IS the reusable
# `qc_pipeline.sh` this module's README asked for as its "suggested next
# step": a standalone, dedicated QC module script, not QC logic buried
# inside a larger domain pipeline (as it is in modules 2/4/6/7 today).
###############################################################################

set -e
set -o pipefail

SCRIPT_DIR="${PIPELINE_SCRIPT_DIR:-${SLURM_SUBMIT_DIR:+${SLURM_SUBMIT_DIR}/modules/module-01-raw-reads-qc/scripts}}"
SCRIPT_DIR="${SCRIPT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)}"
source "${SCRIPT_DIR}/load_modules.sh"

module purge
load_sra_toolkit || exit 1
load_fastqc      || exit 1
load_fastp       || exit 1
load_multiqc     || exit 1

cd /scratch/$(whoami)
mkdir -p module1_qc/{raw_reads,checksums,qc_before,trimmed_reads,fastp_reports,qc_after,logs,pairing_reports}
cd module1_qc

SRR_LIST=("SRR40359383" "SRR40359384")
FAILED_SRRS=()
PASSED_SRRS=()

check_pair() {
    local MATE1="$1" MATE2="$2"
    [ -s "${MATE1}" ] && [ -s "${MATE2}" ] || return 1
    local C1 C2
    C1=$(( $(wc -l < "${MATE1}") / 4 ))
    C2=$(( $(wc -l < "${MATE2}") / 4 ))
    [ "${C1}" -eq "${C2}" ] || return 1
    local IDS1 IDS2
    IDS1=$(mktemp); IDS2=$(mktemp)
    awk 'NR%4==1' "${MATE1}" | sed -E 's#[ /][12]?(:.*)?$##' > "${IDS1}"
    awk 'NR%4==1' "${MATE2}" | sed -E 's#[ /][12]?(:.*)?$##' > "${IDS2}"
    diff -q "${IDS1}" "${IDS2}" > /dev/null
    local RC=$?
    rm -f "${IDS1}" "${IDS2}"
    return ${RC}
}

for SRR in "${SRR_LIST[@]}"; do
    echo "=========================================================="
    echo " ${SRR}: starting full QC pipeline"
    echo "=========================================================="
    LOG="logs/${SRR}_pipeline.log"

    # ------------------------------------------------------------------------
    # Fault-tolerant per-sample subshell (same pattern as script 3): a
    # failure at ANY stage below aborts only THIS sample's subshell, is
    # logged, and the outer loop continues to the next accession.
    # ------------------------------------------------------------------------
    (
        set -e
        set -o pipefail

        # -- Stage 1: fetch + verify (script 1's logic) --
        if [ ! -s "raw_reads/${SRR}_1.fastq" ] || [ ! -s "raw_reads/${SRR}_2.fastq" ]; then
            prefetch "${SRR}" --output-directory sra_cache >> "${LOG}" 2>&1
            vdb-validate "sra_cache/${SRR}/${SRR}.sra" 2>&1 | tee -a "${LOG}" | grep -q "is consistent"
            fasterq-dump "sra_cache/${SRR}/${SRR}.sra" --split-files --outdir raw_reads --threads 4 >> "${LOG}" 2>&1
        fi
        [ -s "raw_reads/${SRR}_1.fastq" ] && [ -s "raw_reads/${SRR}_2.fastq" ]
        md5sum "raw_reads/${SRR}_1.fastq" "raw_reads/${SRR}_2.fastq" > "checksums/${SRR}_raw.md5"
        echo "  [1/5] fetch+verify OK" | tee -a "${LOG}"

        # -- Stage 2: diagnostic FastQC on raw reads (script 2's logic) --
        fastqc "raw_reads/${SRR}_1.fastq" "raw_reads/${SRR}_2.fastq" -o qc_before -t 4 --extract >> "${LOG}" 2>&1
        echo "  [2/5] raw FastQC OK" | tee -a "${LOG}"

        # -- Stage 3: evidence-based fastp trimming (script 4's logic) --
        FASTP_FLAGS=("--length_required" "36")
        if grep -qE "^(WARN|FAIL)\s+Adapter Content" "qc_before/${SRR}_1_fastqc/summary.txt" "qc_before/${SRR}_2_fastqc/summary.txt"; then
            FASTP_FLAGS+=("--detect_adapter_for_pe")
        fi
        if grep -qE "^(WARN|FAIL)\s+Per base sequence quality" "qc_before/${SRR}_1_fastqc/summary.txt" "qc_before/${SRR}_2_fastqc/summary.txt"; then
            FASTP_FLAGS+=("--cut_tail" "--cut_tail_mean_quality" "20")
        fi
        fastp -i "raw_reads/${SRR}_1.fastq" -I "raw_reads/${SRR}_2.fastq" \
              -o "trimmed_reads/${SRR}_1_clean.fastq" -O "trimmed_reads/${SRR}_2_clean.fastq" \
              "${FASTP_FLAGS[@]}" --thread 4 \
              --json "fastp_reports/${SRR}_fastp.json" --html "fastp_reports/${SRR}_fastp.html" >> "${LOG}" 2>&1
        echo "  [3/5] fastp trim OK (flags: ${FASTP_FLAGS[*]})" | tee -a "${LOG}"

        # -- Stage 4: pairing integrity, raw AND trimmed (script 5's logic) --
        check_pair "raw_reads/${SRR}_1.fastq" "raw_reads/${SRR}_2.fastq" \
            && echo "RAW: pairing intact" > "pairing_reports/${SRR}.txt" \
            || { echo "RAW: PAIRING BROKEN"  > "pairing_reports/${SRR}.txt"; exit 1; }
        check_pair "trimmed_reads/${SRR}_1_clean.fastq" "trimmed_reads/${SRR}_2_clean.fastq" \
            && echo "TRIMMED: pairing intact" >> "pairing_reports/${SRR}.txt" \
            || { echo "TRIMMED: PAIRING BROKEN" >> "pairing_reports/${SRR}.txt"; exit 1; }
        echo "  [4/5] pairing checks OK" | tee -a "${LOG}"

        # -- Stage 5: diagnostic FastQC on trimmed reads (before/after pair) --
        fastqc "trimmed_reads/${SRR}_1_clean.fastq" "trimmed_reads/${SRR}_2_clean.fastq" -o qc_after -t 4 --extract >> "${LOG}" 2>&1
        echo "  [5/5] trimmed FastQC OK" | tee -a "${LOG}"
    )

    if [ $? -eq 0 ]; then
        PASSED_SRRS+=("${SRR}")
        echo "  RESULT: ${SRR} PASSED the full pipeline."
    else
        FAILED_SRRS+=("${SRR}")
        echo "  RESULT: ${SRR} FAILED — see logs/${SRR}_pipeline.log. Continuing."
    fi
done

echo ""
echo "Aggregating before/after reports with MultiQC..."
multiqc qc_before/ -n multiqc_raw_batch.html -o . --force
multiqc qc_after/  -n multiqc_trimmed_batch.html -o . --force

echo ""
echo "=========================================================="
echo " PIPELINE SUMMARY: ${#PASSED_SRRS[@]} passed, ${#FAILED_SRRS[@]} failed"
echo " Passed: ${PASSED_SRRS[*]:-none}"
echo " Failed: ${FAILED_SRRS[*]:-none}"
echo "=========================================================="


# ==============================================================================
# Nuances & Pitfalls
# ==============================================================================
# - This script deliberately re-derives the fastp flag logic from script 4
#   inline rather than sourcing script 4 as a library — that keeps this file
#   runnable standalone via sbatch (Slurm does not "source" companion scripts
#   across job boundaries cleanly), at the cost of duplicating a few lines.
#   If you extend this module, consider factoring shared logic into a
#   `lib/qc_functions.sh` that BOTH the teaching scripts and this integrated
#   pipeline source — the same lib-extraction idea Module 0's capstone notes
#   apply here too.
# - A sample with PAIRING BROKEN after fastp almost always indicates fastp
#   was invoked incorrectly (single-mate mode, or mismatched -i/-o pairing)
#   rather than a fastp bug — fastp preserves pairing correctly when both
#   mates are passed together, as this script does.
# - Like script 3, `--force` on the MultiQC calls means these two combined
#   reports are overwritten on every re-run — that is intentional for a
#   pipeline meant to be re-run as data changes, but back them up first if
#   you need to preserve a specific run's exact report.
# ==============================================================================
