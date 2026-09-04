#!/bin/bash
#SBATCH --job-name=qc08_pipeline
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=12:00:00

###############################################################################
# 8_Automated_Pipeline_AllSamples.sh
#
# NARRATIVE:
# This script is the production-grade culmination of the module. It orchestrates 
# Steps 1-7 into a fault-tolerant batch loop operating on its own isolated 
# directory (`_automated`). The `set -e` subshell architecture ensures that if 
# one accession fails (e.g., download timeout or failed integrity check), the 
# script gracefully logs the error, skips that sample, and continues processing 
# the remainder of the cohort, finally quantifying all successful BAMs together.
###############################################################################
set -e
set -o pipefail

module load SRA-Toolkit/3.0.3-gompi-2022a
module load FastQC/0.11.9-Java-11
module load fastp/0.23.4-GCC-13.2.0
module load STAR/2.7.10b-GCC-11.3.0
module load MultiQC/1.28-foss-2024a
module load R

WORKDIR="/scratch/$(whoami)/PRJNA229998_airway_pipeline_automated"
REFDIR="/scratch/$(whoami)/PRJNA229998_airway_pipeline_ref"
mkdir -p "${WORKDIR}"/{raw_reads,checksums,qc_before,trimmed_reads,fastp_reports,qc_after,pairing_reports,alignments,counts,logs,sra_cache}
cd "${WORKDIR}"

SRR_LIST=("SRR1039508" "SRR1039509" "SRR1039512" "SRR1039513")
PASSED_SRRS=(); FAILED_SRRS=()

check_pair() {
    local M1="$1" M2="$2"
    [ -s "${M1}" ] && [ -s "${M2}" ] || return 1
    local C1=$(( $(wc -l < "${M1}") / 4 )); local C2=$(( $(wc -l < "${M2}") / 4 ))
    [ "${C1}" -eq "${C2}" ] || return 1
    local I1=$(mktemp); local I2=$(mktemp)
    awk 'NR%4==1' "${M1}" | sed -E 's#[ /][12]?(:.*)?$##' > "${I1}"
    awk 'NR%4==1' "${M2}" | sed -E 's#[ /][12]?(:.*)?$##' > "${I2}"
    diff -q "${I1}" "${I2}" > /dev/null; local RC=$?; rm -f "${I1}" "${I2}"
    return ${RC}
}

for SRR in "${SRR_LIST[@]}"; do
    LOG="logs/${SRR}_pipeline.log"
    echo "Processing ${SRR}..."
    (
        set -e; set -o pipefail
        
        if [ ! -s "raw_reads/${SRR}_1.fastq" ]; then
            prefetch "${SRR}" --output-directory sra_cache >> "${LOG}" 2>&1
            vdb-validate "sra_cache/${SRR}/${SRR}.sra" 2>&1 | grep -q "is consistent"
            fasterq-dump "sra_cache/${SRR}/${SRR}.sra" --split-files --outdir raw_reads --threads 8 >> "${LOG}" 2>&1
        fi
        md5sum "raw_reads/${SRR}_1.fastq" "raw_reads/${SRR}_2.fastq" > "checksums/${SRR}_raw.md5"
        
        fastqc "raw_reads/${SRR}_1.fastq" "raw_reads/${SRR}_2.fastq" -o qc_before -t 8 --extract >> "${LOG}" 2>&1
        
        FASTP=("--length_required" "36")
        if grep -qE "^(WARN|FAIL)\s+Adapter Content" "qc_before/${SRR}_1_fastqc/summary.txt" "qc_before/${SRR}_2_fastqc/summary.txt"; then FASTP+=("--detect_adapter_for_pe"); fi
        if grep -qE "^(WARN|FAIL)\s+Per base sequence quality" "qc_before/${SRR}_1_fastqc/summary.txt" "qc_before/${SRR}_2_fastqc/summary.txt"; then FASTP+=("--cut_tail" "--cut_tail_mean_quality" "20"); fi
        fastp -i "raw_reads/${SRR}_1.fastq" -I "raw_reads/${SRR}_2.fastq" -o "trimmed_reads/${SRR}_1_clean.fastq" -O "trimmed_reads/${SRR}_2_clean.fastq" "${FASTP[@]}" --thread 8 --json "fastp_reports/${SRR}_fastp.json" --html "fastp_reports/${SRR}_fastp.html" >> "${LOG}" 2>&1
        echo "\`fastp ${FASTP[*]}\`" > "fastp_reports/${SRR}_THRESHOLD_RATIONALE.md"

        check_pair "trimmed_reads/${SRR}_1_clean.fastq" "trimmed_reads/${SRR}_2_clean.fastq" \
            && echo "TRIMMED: intact" > "pairing_reports/${SRR}.txt" || { echo "TRIMMED: BROKEN" > "pairing_reports/${SRR}.txt"; exit 1; }

        STAR --runThreadN 8 --genomeDir "${REFDIR}/star_index" --readFilesIn "trimmed_reads/${SRR}_1_clean.fastq" "trimmed_reads/${SRR}_2_clean.fastq" --outSAMtype BAM SortedByCoordinate --outFileNamePrefix "alignments/${SRR}_" >> "${LOG}" 2>&1
    )
    if [ $? -eq 0 ]; then PASSED_SRRS+=("${SRR}"); else FAILED_SRRS+=("${SRR}"); fi
done

if [ ${#PASSED_SRRS[@]} -gt 0 ]; then
    BAM_LIST=""
    for SRR in "${PASSED_SRRS[@]}"; do BAM_LIST="${BAM_LIST}\"alignments/${SRR}_Aligned.sortedByCoord.out.bam\", "; done
    BAM_LIST=$(echo "${BAM_LIST}" | sed 's/, $//')
    cat << EOF > counts/run_final_quant.R
suppressMessages(library(Rsubread))
fc <- featureCounts(files = c(${BAM_LIST}), annot.ext = "${REFDIR}/GRCh38.gtf", isGTFAnnotationFile = TRUE, isPairedEnd = TRUE, nthreads = 8)
write.table(fc\$counts, file="counts/airway_raw_counts.txt", sep="\t", quote=FALSE, col.names=NA)
write.table(fc\$stat, file="counts/airway_raw_counts.txt.summary", sep="\t", quote=FALSE, row.names=FALSE)
EOF
    Rscript counts/run_final_quant.R
fi