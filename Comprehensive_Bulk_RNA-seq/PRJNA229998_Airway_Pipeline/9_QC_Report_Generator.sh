#!/bin/bash
#SBATCH --job-name=qc09_report
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:30:00

###############################################################################
# 9_QC_Report_Generator.sh
#
# NARRATIVE:
# This script consolidates the audit trail. By gathering disparate outputs—
# checksums, trimming rationales, pairing verdicts, and STAR alignment metrics—
# it dynamically constructs a comprehensive, reproducible Markdown portfolio 
# artifact. This final QC report serves as the formal "go/no-go" boundary 
# before proceeding to differential expression analysis.
###############################################################################
set -e

WORKDIR="/scratch/$(whoami)/PRJNA229998_airway_pipeline_automated"
cd "${WORKDIR}"

SRR_LIST=("SRR1039508" "SRR1039509" "SRR1039512" "SRR1039513")
REPORT="QC_REPORT.md"

{
echo "# PRJNA229998 (Airway) Automated Pipeline Report"
echo "Generated on $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo ""
echo "## 1. Checksums"
for SRR in "${SRR_LIST[@]}"; do
    [ -f "checksums/${SRR}_raw.md5" ] && cat "checksums/${SRR}_raw.md5" || echo "${SRR}: MISSING"
done
echo ""
echo "## 2. Preprocessing & Integrity"
for SRR in "${SRR_LIST[@]}"; do
    if [ -f "fastp_reports/${SRR}_THRESHOLD_RATIONALE.md" ]; then
        echo -n "**${SRR} Flags:** "
        cat "fastp_reports/${SRR}_THRESHOLD_RATIONALE.md"
        echo -n "**Pairing:** "
        cat "pairing_reports/${SRR}.txt"
        echo ""
    fi
done
echo "## 3. Alignment & Quantification Metrics"
echo "| Sample | STAR Unique Mapping | featureCounts Assigned |"
echo "|---|---|---|"
for SRR in "${SRR_LIST[@]}"; do
    if [ -f "alignments/${SRR}_Log.final.out" ]; then
        STAR_RATE=$(grep "Uniquely mapped reads %" "alignments/${SRR}_Log.final.out" | awk -F'|\t+' '{print $2}' | tr -d ' \t')
        echo "| ${SRR} | ${STAR_RATE} | See counts/airway_raw_counts.txt.summary |"
    else
        echo "| ${SRR} | Failed/Missing | Failed/Missing |"
    fi
done
} > "${REPORT}"

echo "OK: Wrote ${REPORT} in ${WORKDIR}"

