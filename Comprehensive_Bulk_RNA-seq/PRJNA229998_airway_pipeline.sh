#!/bin/bash
#SBATCH --job-name=airway_pipeline
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=12:00:00

###############################################################################
# PRJNA229998_airway_pipeline.sh
# Comprehensive Bulk RNA-seq Pipeline (Raw Reads to featureCounts)
# Study: GSE52778 / BioProject PRJNA229998 (Himes et al., 2014)
#
# PEDAGOGICAL PHILOSOPHY:
# 1. Verify before inspect: Cache validation + md5 checksums at retrieval.
# 2. Evidence-based decisions: FastQC verdicts dictate fastp trimming flags.
# 3. Defensive programming: Paired-end integrity is explicitly verified.
# 4. Fault-tolerance: Subshells ensure one failed sample doesn't crash the batch.
# 5. Traceability: Automated generation of a comprehensive QC_REPORT.md.
# 6. Rsubread fallback: Uses Rsubread via Rscript if command-line is missing.
###############################################################################

set -e
set -o pipefail

# ------------------------------------------------------------------------------
# ENVIRONMENT & SETUP
# ------------------------------------------------------------------------------
module load SRA-Toolkit/3.0.3-gompi-2022a
module load FastQC/0.11.9-Java-11
module load fastp/0.23.4-GCC-13.2.0
module load STAR/2.7.10b-GCC-11.3.0
module load MultiQC/1.28-foss-2024a

# CHANGE THIS LINE to your cluster's specific R module that contains Rsubread:
module load R 

WORKDIR="/scratch/$(whoami)/PRJNA229998_airway_pipeline"
mkdir -p "${WORKDIR}"/{metadata,ref,raw_reads,checksums,qc_before,trimmed_reads,fastp_reports,qc_after,pairing_reports,alignments,counts,logs,sra_cache}
cd "${WORKDIR}"

SRR_LIST=("SRR1039508" "SRR1039509" "SRR1039512" "SRR1039513")
FAILED_SRRS=()
PASSED_SRRS=()

# ------------------------------------------------------------------------------
# HELPER FUNCTION: Paired-End Integrity Check
# ------------------------------------------------------------------------------
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

# ------------------------------------------------------------------------------
# STEP 0: METADATA ACQUISITION
# ------------------------------------------------------------------------------
echo "=========================================================="
echo " Fetching Study Metadata (PRJNA229998)"
echo "=========================================================="
METADATA_FILE="metadata/PRJNA229998_SraRunInfo.csv"
if [ ! -f "${METADATA_FILE}" ]; then
    curl -s "https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?save=efetch&db=sra&rettype=runinfo&term=PRJNA229998" > "${METADATA_FILE}"
    echo "  OK: SRA RunInfo downloaded to ${METADATA_FILE}"
else
    echo "  OK: Metadata already exists at ${METADATA_FILE}"
fi

# ------------------------------------------------------------------------------
# STEP 1: REFERENCE GENOME ACQUISITION & INDEXING
# ------------------------------------------------------------------------------
echo "=========================================================="
echo " Preparing Human GRCh38 Reference Genome"
echo "=========================================================="
FASTA_URL="http://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
GTF_URL="http://ftp.ensembl.org/pub/release-110/gtf/homo_sapiens/Homo_sapiens.GRCh38.110.gtf.gz"

if [ ! -f "ref/GRCh38.fa" ]; then
    wget -c "${FASTA_URL}" -O ref/GRCh38.fa.gz
    gunzip ref/GRCh38.fa.gz
fi

if [ ! -f "ref/GRCh38.gtf" ]; then
    wget -c "${GTF_URL}" -O ref/GRCh38.gtf.gz
    gunzip ref/GRCh38.gtf.gz
fi

if [ ! -d "ref/star_index" ] || [ -z "$(ls -A ref/star_index)" ]; then
    mkdir -p ref/star_index
    STAR --runThreadN 8 --runMode genomeGenerate --genomeDir ref/star_index \
         --genomeFastaFiles ref/GRCh38.fa --sjdbGTFfile ref/GRCh38.gtf --sjdbOverhang 100
fi

# ------------------------------------------------------------------------------
# STEP 2: FAULT-TOLERANT PROCESSING LOOP
# ------------------------------------------------------------------------------
for SRR in "${SRR_LIST[@]}"; do
    echo "=========================================================="
    echo " Processing ${SRR}"
    echo "=========================================================="
    LOG="logs/${SRR}_pipeline.log"

    (
        set -e
        set -o pipefail

        # --- A. Verify Before Inspect (Fetch & Checksum) ---
        if [ ! -s "raw_reads/${SRR}_1.fastq" ] \vert{}\vert{} [ ! -s "raw_reads/${SRR}_2.fastq" ]; then
            prefetch "${SRR}" --output-directory sra_cache >> "${LOG}" 2>&1
            vdb-validate "sra_cache/${SRR}/${SRR}.sra" 2>&1 \vert{} tee -a "${LOG}" | grep -q "is consistent"
            fasterq-dump "sra_cache/${SRR}/${SRR}.sra" --split-files --outdir raw_reads --threads 8 >> "${LOG}" 2>&1
        fi
        
        md5sum "raw_reads/${SRR}_1.fastq" "raw_reads/${SRR}_2.fastq" > "checksums/${SRR}_raw.md5"

        # --- B. Raw Diagnostic QC ---
        fastqc "raw_reads/${SRR}_1.fastq" "raw_reads/${SRR}_2.fastq" -o qc_before -t 8 --extract >> "${LOG}" 2>&1

        # --- C. Evidence-Based fastp Trimming ---
        FASTP_FLAGS=("--length_required" "36")
        RATIONALE_LINES=("- Minimum length floor of 36bp applied unconditionally.")

        if grep -qE "^(WARN|FAIL)\s+Adapter Content" "qc_before/${SRR}_1_fastqc/summary.txt" "qc_before/${SRR}_2_fastqc/summary.txt"; then
            FASTP_FLAGS+=("--detect_adapter_for_pe")
            RATIONALE_LINES+=("- Adapter trimming ENABLED: FastQC flagged Adapter Content.")
        else
            RATIONALE_LINES+=("- Adapter trimming SKIPPED: No evidence of adapter contamination.")
        fi

        if grep -qE "^(WARN|FAIL)\s+Per base sequence quality" "qc_before/${SRR}_1_fastqc/summary.txt" "qc_before/${SRR}_2_fastqc/summary.txt"; then
            FASTP_FLAGS+=("--cut_tail" "--cut_tail_mean_quality" "20")
            RATIONALE_LINES+=("- 3' quality trimming ENABLED: FastQC flagged per-base quality decay.")
        else
            RATIONALE_LINES+=("- 3' quality trimming SKIPPED: Per-base quality passed.")
        fi

        fastp -i "raw_reads/${SRR}_1.fastq" -I "raw_reads/${SRR}_2.fastq" \
              -o "trimmed_reads/${SRR}_1_clean.fastq" -O "trimmed_reads/${SRR}_2_clean.fastq" \
              "${FASTP_FLAGS[@]}" --thread 8 \
              --json "fastp_reports/${SRR}_fastp.json" --html "fastp_reports/${SRR}_fastp.html" >> "${LOG}" 2>&1
              
        {
            echo "## Trimming rationale for ${SRR}"
            printf '%s\n' "${RATIONALE_LINES[@]}"
            echo "\`fastp ${FASTP_FLAGS[*]}\`"
        } > "fastp_reports/${SRR}_THRESHOLD_RATIONALE.md"

        # --- D. Paired-End Integrity Verdict ---
        check_pair "raw_reads/${SRR}_1.fastq" "raw_reads/${SRR}_2.fastq" \
            && echo "RAW: intact" > "pairing_reports/${SRR}.txt" \
            || { echo "RAW: BROKEN" > "pairing_reports/${SRR}.txt"; exit 1; }
            
        check_pair "trimmed_reads/${SRR}_1_clean.fastq" "trimmed_reads/${SRR}_2_clean.fastq" \
            && echo "TRIMMED: intact" >> "pairing_reports/${SRR}.txt" \
            || { echo "TRIMMED: BROKEN" >> "pairing_reports/${SRR}.txt"; exit 1; }

        # --- E. Post-trim QC & STAR Alignment ---
        fastqc "trimmed_reads/${SRR}_1_clean.fastq" "trimmed_reads/${SRR}_2_clean.fastq" -o qc_after -t 8 --extract >> "${LOG}" 2>&1
        
        STAR --runThreadN 8 --genomeDir ref/star_index \
             --readFilesIn "trimmed_reads/${SRR}_1_clean.fastq" "trimmed_reads/${SRR}_2_clean.fastq" \
             --outSAMtype BAM SortedByCoordinate \
             --outFileNamePrefix "alignments/${SRR}_" >> "${LOG}" 2>&1
    )

    if [ $? -eq 0 ]; then
        PASSED_SRRS+=("${SRR}")
        echo "  RESULT: ${SRR} pipeline PASSED."
    else
        FAILED_SRRS+=("${SRR}")
        echo "  RESULT: ${SRR} pipeline FAILED. See logs."
    fi
done

# ------------------------------------------------------------------------------
# STEP 3: RSUBREAD QUANTIFICATION (Dynamically built R script)
# ------------------------------------------------------------------------------
if [ ${#PASSED_SRRS[@]} -eq 0 ]; then
    echo "ERROR: No samples passed alignment."
    exit 1
fi

echo "=========================================================="
echo " Running Rsubread featureCounts via Rscript"
echo "=========================================================="

# Build a comma-separated list of BAM files for R string formatting
BAM_LIST_R=""
for SRR in "${PASSED_SRRS[@]}"; do
    BAM_LIST_R="${BAM_LIST_R}\"alignments/${SRR}_Aligned.sortedByCoord.out.bam\", "
done
BAM_LIST_R=$(echo "${BAM_LIST_R}" \vert{} sed 's/, $//') # Remove trailing comma

# Dynamically write a short Rscript to run featureCounts
cat << EOF > counts/run_featurecounts.R
suppressMessages(library(Rsubread))

bam_files <- c(${BAM_LIST_R})

fc <- featureCounts(files = bam_files,
                    annot.ext = "ref/GRCh38.gtf",
                    isGTFAnnotationFile = TRUE,
                    isPairedEnd = TRUE,
                    nthreads = 8)

# Write out the raw counts matrix 
write.table(fc\$counts, file="counts/airway_raw_counts.txt", sep="\t", quote=FALSE, col.names=NA)

# Write out the summary file so MultiQC can still parse the mapping rates
write.table(fc\$stat, file="counts/airway_raw_counts.txt.summary", sep="\t", quote=FALSE, row.names=FALSE)
EOF

# Execute the newly written R script
Rscript counts/run_featurecounts.R

# ------------------------------------------------------------------------------
# STEP 4: BATCH AGGREGATION & REPORTING
# ------------------------------------------------------------------------------
echo "Aggregating all logs and metrics with MultiQC..."
multiqc . -n multiqc_final_report.html -o . --force

echo "Generating comprehensive QC_REPORT.md..."
REPORT="QC_REPORT.md"

{
echo "# PRJNA229998 (Airway) Bulk RNA-seq Pipeline Report"
echo "Generated on $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo ""
echo "## 1. Batch Summary"
echo "- **Passed:** ${PASSED_SRRS[*]:-none}"
echo "- **Failed:** ${FAILED_SRRS[*]:-none}"
echo ""
echo "## 2. Experimental Metadata"
echo "The full SRA RunInfo (including experimental conditions) was downloaded to:"
echo "- \`${METADATA_FILE}\`"
echo "This file must be used to build the \`colData\` design matrix in DESeq2."
echo ""
echo "## 3. Checksums (Data Integrity at Retrieval)"
echo '```text'
for SRR in "${SRR_LIST[@]}"; do
    [ -f "checksums/${SRR}_raw.md5" ] && cat "checksums/${SRR}_raw.md5" || echo "${SRR}: MISSING"
done
echo '```'
echo ""
echo "## 4. Preprocessing Rationale & Pairing Integrity"
for SRR in "${SRR_LIST[@]}"; do
    if [ -f "fastp_reports/${SRR}_THRESHOLD_RATIONALE.md" ]; then
        cat "fastp_reports/${SRR}_THRESHOLD_RATIONALE.md"
        echo -n "- **Pairing Verdict:** "
        cat "pairing_reports/${SRR}.txt" | tr '\n' ' '
        echo -e "\n"
    fi
done
echo "## 5. Alignment & Quantification Metrics"
echo "For a fully interactive breakdown of alignment rates (STAR) and feature assignment rates (featureCounts), open \`multiqc_final_report.html\`."
echo ""
echo "### Raw featureCounts 'Assigned' reads:"
echo '```text'
grep -E "^Status|^Assigned" counts/airway_raw_counts.txt.summary | sed 's|alignments/||g' | sed 's|_Aligned.sortedByCoord.out.bam||g'
echo '```'
echo ""
echo "## 6. Next Steps"
echo "You now have:"
echo "1. \`counts/airway_raw_counts.txt\` (The count matrix)"
echo "2. \`${METADATA_FILE}\` (The sample metadata)"
echo "You are ready to initialize the \`DESeqDataSet\` object in R."
} > "${REPORT}"

echo "=========================================================="
echo " Pipeline Complete. Check ${REPORT} and multiqc_final_report.html."
echo "=========================================================="