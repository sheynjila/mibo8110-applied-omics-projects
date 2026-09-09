#!/bin/bash
#SBATCH --job-name=qc02_fastqc
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=6G
#SBATCH --time=00:45:00

###############################################################################
# 2_Interpret_FastQC_Diagnostics.sh
# Module 1 (Raw Sequencing Reads & QC) — Script 2 of 7
#
# ONE NEW IDEA on top of script 1: it is not enough to RUN FastQC — you have
# to read its verdict programmatically and decide what it means. This script
# runs FastQC on the mates fetched in script 1, then parses FastQC's own
# summary.txt (bundled inside its output zip) into a plain PASS/WARN/FAIL
# synthesis, and prints a short interpretation note for the modules students
# most often misread: per-base quality, adapter content, and duplication.
###############################################################################

set -e
set -o pipefail

SCRIPT_DIR="${PIPELINE_SCRIPT_DIR:-${SLURM_SUBMIT_DIR:+${SLURM_SUBMIT_DIR}/modules/module-01-raw-reads-qc/scripts}}"
SCRIPT_DIR="${SCRIPT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)}"
source "${SCRIPT_DIR}/load_modules.sh"

module purge
load_fastqc || exit 1

cd /scratch/$(whoami)/module1_qc
mkdir -p qc_before

SRR="SRR40359383"

if [ ! -s "raw_reads/${SRR}_1.fastq" ] || [ ! -s "raw_reads/${SRR}_2.fastq" ]; then
    echo "ERROR: raw_reads/${SRR}_{1,2}.fastq not found. Run script 1 first." >&2
    exit 1
fi

echo "Running FastQC on ${SRR} (both mates)..."
fastqc "raw_reads/${SRR}_1.fastq" "raw_reads/${SRR}_2.fastq" -o qc_before -t 4 --extract

# ------------------------------------------------------------------------------
# TEACHING NOTE: FastQC writes a plain-text summary.txt inside each
# <sample>_fastqc/ extraction directory with one PASS/WARN/FAIL line per
# module (Per base sequence quality, Adapter Content, Sequence Duplication
# Levels, etc.). Parsing that file — instead of only eyeballing the HTML
# report — is what lets this become a repeatable, scriptable decision later
# in script 4, rather than a one-off human judgment call.
# ------------------------------------------------------------------------------
for MATE in 1 2; do
    SUMMARY="qc_before/${SRR}_${MATE}_fastqc/summary.txt"
    if [ ! -f "${SUMMARY}" ]; then
        echo "ERROR: expected summary at ${SUMMARY} was not produced." >&2
        exit 1
    fi
    echo ""
    echo "---- ${SRR} mate ${MATE}: FastQC verdicts ----"
    cat "${SUMMARY}"
done

# ------------------------------------------------------------------------------
# Interpretation synthesis: count FAIL/WARN lines across both mates and print
# a plain-language note on the three metrics students most often misread.
# ------------------------------------------------------------------------------
FAIL_COUNT=$(cat qc_before/${SRR}_1_fastqc/summary.txt qc_before/${SRR}_2_fastqc/summary.txt | grep -c "^FAIL" || true)
WARN_COUNT=$(cat qc_before/${SRR}_1_fastqc/summary.txt qc_before/${SRR}_2_fastqc/summary.txt | grep -c "^WARN" || true)

echo ""
echo "=========================================================="
echo " Synthesis for ${SRR}: ${FAIL_COUNT} FAIL, ${WARN_COUNT} WARN (across both mates)"
echo "=========================================================="
echo ""
echo "How to read the three most commonly misread modules:"
echo ""
echo "1. Per base sequence quality — a WARN/FAIL near the READ TAIL is normal"
echo "   for most Illumina chemistries (quality decays toward the 3' end) and"
echo "   is exactly what evidence-based trimming in script 4 exists to fix."
echo "   A FAIL near the READ START is more concerning and worth inspecting"
echo "   the HTML plot directly before deciding anything."
echo ""
echo "2. Adapter Content — a WARN/FAIL here means FastQC detected adapter"
echo "   sequence contaminating reads (common when insert size < read length)."
echo "   This is diagnostic evidence to justify adapter trimming in script 4 —"
echo "   FastQC never trims anything itself, it only tells you whether to."
echo ""
echo "3. Sequence Duplication Levels — a WARN/FAIL here is AMBIGUOUS on its"
echo "   own: high duplication can mean PCR over-amplification (a technical"
echo "   artifact worth flagging) OR genuinely highly-expressed transcripts"
echo "   in RNA-seq data (real biology). Do not auto-deduplicate reads based"
echo "   on this flag alone — that is exactly the kind of \"passing QC report"
echo "   does not prove biological validity\" misconception this curriculum"
echo "   calls out explicitly."

echo ""
echo "Done: FastQC diagnostics for ${SRR} parsed and interpreted."


# ==============================================================================
# Nuances & Pitfalls
# ==============================================================================
# - FastQC's PASS/WARN/FAIL thresholds are generic defaults tuned for
#   general-purpose sequencing, not your specific assay. A WARN on GC content
#   is expected and NOT concerning for organisms with unusual genome-wide GC
#   composition — always check what "normal" looks like for your organism
#   before treating a WARN as a problem.
# - Reading only the HTML report (opening it in a browser) does not scale
#   past a handful of samples. Parsing summary.txt, as this script does, is
#   the bridge to script 3's MultiQC batch aggregation and script 4's
#   evidence-based automated trimming decisions.
# - `--extract` is required to get the unzipped summary.txt directly; without
#   it FastQC only writes the .zip and you would need to unzip it yourself.
# ==============================================================================
