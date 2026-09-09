#!/bin/bash
#SBATCH --job-name=qc03_multiqc
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=01:30:00

###############################################################################
# 3_Batch_MultiQC_Summary.sh
# Module 1 (Raw Sequencing Reads & QC) — Script 3 of 7
#
# ONE NEW IDEA on top of script 2: scripts 1-2 only ever looked at ONE run
# (SRR40359383). Real studies have multiple runs, and comparing them side by
# side is exactly what catches a single bad sample before it contaminates a
# downstream analysis. This script loops the fetch+FastQC steps over BOTH
# real runs in this study, fault-tolerantly, then aggregates everything with
# MultiQC into one comparison report.
###############################################################################

set -e
set -o pipefail

SCRIPT_DIR="${PIPELINE_SCRIPT_DIR:-${SLURM_SUBMIT_DIR:+${SLURM_SUBMIT_DIR}/modules/module-01-raw-reads-qc/scripts}}"
SCRIPT_DIR="${SCRIPT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)}"
source "${SCRIPT_DIR}/load_modules.sh"

module purge
load_sra_toolkit || exit 1
load_fastqc      || exit 1
load_multiqc     || exit 1

cd /scratch/$(whoami)/module1_qc
mkdir -p raw_reads checksums qc_before logs

SRR_LIST=("SRR40359383" "SRR40359384")
FAILED_SRRS=()

for SRR in "${SRR_LIST[@]}"; do
    echo "=========================================================="
    echo " Processing ${SRR}"
    echo "=========================================================="

    # ------------------------------------------------------------------------
    # TEACHING NOTE (fault-tolerant loop, same pattern used across this repo):
    # run each sample's fetch+QC in its own subshell with its OWN `set -e`.
    # If SRR40359384 fails, the loop logs it and moves on to check any
    # remaining samples, instead of one bad accession stopping the entire
    # batch — you get a full report on everything that DID work, plus a clear
    # list of what needs a second look.
    # ------------------------------------------------------------------------
    (
        set -e
        set -o pipefail

        if [ ! -s "raw_reads/${SRR}_1.fastq" ] || [ ! -s "raw_reads/${SRR}_2.fastq" ]; then
            prefetch "${SRR}" --output-directory sra_cache >> "logs/${SRR}_fetch.log" 2>&1
            vdb-validate "sra_cache/${SRR}/${SRR}.sra" 2>&1 | grep -q "is consistent"
            fasterq-dump "sra_cache/${SRR}/${SRR}.sra" --split-files --outdir raw_reads --threads 4 >> "logs/${SRR}_fetch.log" 2>&1
        else
            echo "  raw_reads/${SRR}_{1,2}.fastq already present — skipping re-fetch."
        fi

        [ -s "raw_reads/${SRR}_1.fastq" ] && [ -s "raw_reads/${SRR}_2.fastq" ]
        md5sum "raw_reads/${SRR}_1.fastq" "raw_reads/${SRR}_2.fastq" > "checksums/${SRR}_raw.md5"

        fastqc "raw_reads/${SRR}_1.fastq" "raw_reads/${SRR}_2.fastq" -o qc_before -t 4
        echo "  OK: ${SRR} fetched, verified, and FastQC'd."
    )

    if [ $? -ne 0 ]; then
        echo "  FAILED: ${SRR} — see logs/${SRR}_fetch.log. Continuing with remaining samples."
        FAILED_SRRS+=("${SRR}")
    fi
done

echo ""
echo "Aggregating all qc_before/ reports with MultiQC..."
multiqc qc_before/ -n multiqc_raw_batch.html -o . --force

echo ""
echo "=========================================================="
echo " Batch summary: ${#SRR_LIST[@]} samples attempted, ${#FAILED_SRRS[@]} failed"
if [ ${#FAILED_SRRS[@]} -gt 0 ]; then
    echo " Failed accessions: ${FAILED_SRRS[*]}"
fi
echo " Combined report: multiqc_raw_batch.html"
echo "=========================================================="


# ==============================================================================
# Why MultiQC instead of reading each FastQC report separately
# ==============================================================================
# MultiQC does not run any new QC — it parses the FastQC outputs already sitting
# in qc_before/ and renders ONE report with every sample's metrics as parallel
# bars/lines you can compare directly. This is what turns "does SRR40359384
# look worse than SRR40359383?" from a manual side-by-side comparison of two
# HTML files into a single glance at one chart. At 2 samples this is a
# convenience; at 20+ samples it is the only practical way to spot the one
# outlier sample before it silently drags down a whole analysis.
#
# Nuances & Pitfalls:
#   - MultiQC only sees what's IN the directory you point it at. If a sample's
#     FastQC step failed (as tracked by FAILED_SRRS above), it simply will not
#     appear in the combined report — silence, not a visible warning. Always
#     cross-check the batch summary line count against your FAILED_SRRS list.
#   - `--force` overwrites a previous multiqc_raw_batch.html. That is
#     intentional here (this script is meant to be re-run as more samples are
#     added) but means the report does not accumulate history on its own —
#     rename or archive it first if you need to keep an old version.
# ==============================================================================
