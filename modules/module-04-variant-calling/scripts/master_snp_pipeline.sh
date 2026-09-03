#!/bin/bash
#SBATCH --job-name=master_snp_caller
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --output=master_snp_%j.out
#SBATCH --error=master_snp_%j.err

# ==============================================================================
# CORRECTED VERSION of master_snp_pipeline.sh
# See pipeline_appraisal_2026-08-31.md for the full review. Fixes applied:
#
#   1. The multi-sample VCF merge (Phase 5) is now automated inside the
#      script instead of being a manual copy-paste block in a trailing
#      comment. Any *_filtered.vcf produced this run is indexed and merged
#      into merged_cohort.vcf automatically.
#   2. Reference FASTA download verifies non-empty output before treating a
#      rerun as "already downloaded".
#   3. The per-sample loop is fault-tolerant: a failed sample (download,
#      alignment, or calling failure) is logged and skipped instead of
#      aborting the whole run, and failed samples are excluded from the
#      cohort merge.
# ==============================================================================

set -e
set -o pipefail

echo "=========================================================="
echo " INITIATING MASTER BACTERIAL VARIANT CALLING PIPELINE "
echo "=========================================================="

# ==========================================
# PHASE 1: DIRECTORY & ENVIRONMENT SETUP
# ==========================================
WORKDIR="/scratch/$(whoami)/master_snp_pipeline"
mkdir -p ${WORKDIR}/{ref,raw_reads,qc,trimmed_reads,alignment,variants,logs}
cd ${WORKDIR}

if [ ! -f "srr_list.txt" ]; then
    echo "WARNING: srr_list.txt not found. Generating default test list..."
    echo -e "SRR11092056\nSRR11092057" > srr_list.txt
fi
sed -i 's/\r$//' srr_list.txt

# ==========================================
# PHASE 2: ACQUIRE & PREPARE DNA REFERENCE
# ==========================================
echo -e "\n[PHASE 2] PREPARING BACTERIAL DNA REFERENCE..."
cd ref

FASTA_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/006/945/GCF_000006945.2_ASM694v2/GCF_000006945.2_ASM694v2_genomic.fna.gz"
FASTA_FILE="genome.fna"

# FIX #2: -s checks non-empty, so a truncated download gets re-fetched.
if [ ! -s "${FASTA_FILE}" ]; then
    echo "Downloading reference FASTA..."
    wget -qO- ${FASTA_URL} | gunzip -c > ${FASTA_FILE}
    if [ ! -s "${FASTA_FILE}" ]; then
        echo "CRITICAL ERROR: ${FASTA_FILE} downloaded empty/corrupt. Aborting."
        rm -f "${FASTA_FILE}"
        exit 1
    fi
fi

if [ ! -f "${FASTA_FILE}.bwt" ]; then
    echo "Building BWA reference index..."
    module purge
    module load BWA/0.7.18-GCC-13.3.0
    bwa index ${FASTA_FILE}
fi

if [ ! -f "${FASTA_FILE}.fai" ]; then
    echo "Building SAMtools FASTA index..."
    module purge
    module load SAMtools/1.21-GCC-13.3.0
    samtools faidx ${FASTA_FILE}
fi
cd ${WORKDIR}

# ==========================================
# PHASE 3: SAMPLE PROCESSING LOOP
# ==========================================
echo -e "\n[PHASE 3] EXECUTING WGS SAMPLE PROCESSING LOOP..."

# FIX #3: track failures instead of a hard abort on the first bad sample.
FAILED_SAMPLES=()

for SRR in $(cat srr_list.txt); do
    echo "----------------------------------------------------------"
    echo " PROCESSING SAMPLE: ${SRR}"
    echo "----------------------------------------------------------"

    set +e
    (
        set -e
        set -o pipefail

        # Step A: Download (Two-step Prefetch Fix)
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

        # Step C: BWA-MEM DNA Alignment & Direct Sorting
        echo ">> Aligning DNA reads to reference..."
        module purge
        module load BWA/0.7.18-GCC-13.3.0
        module load SAMtools/1.21-GCC-13.3.0

        bwa mem -t 8 -R "@RG\tID:${SRR}\tSM:${SRR}\tPL:ILLUMINA" ref/${FASTA_FILE} \
                trimmed_reads/${SRR}_1_clean.fastq \
                trimmed_reads/${SRR}_2_clean.fastq | \
        samtools sort -@ 4 -o alignment/${SRR}_unsorted.bam -

        # Step D: Mark and Remove PCR/Optical Duplicates
        echo ">> Marking and filtering PCR duplicates with SAMtools..."
        samtools collate -@ 4 -O -u alignment/${SRR}_unsorted.bam | \
        samtools fixmate -@ 4 -m -u - - | \
        samtools sort -@ 4 -u - | \
        samtools markdup -@ 4 -r - alignment/${SRR}_dedup.bam

        samtools index alignment/${SRR}_dedup.bam
        rm alignment/${SRR}_unsorted.bam # Clean up intermediate file

        # Step E: Variant Calling (Generating Raw VCF)
        echo ">> Calling variants (SNPs and INDELs)..."
        module purge
        module load BCFtools/1.21-GCC-13.3.0

        bcftools mpileup -O b -f ref/${FASTA_FILE} alignment/${SRR}_dedup.bam | \
        bcftools call -mv -O v -o variants/${SRR}_raw.vcf

        # Step F: Multi-Metric Robust Variant Filtering
        echo ">> Applying robust multi-metric filtering (QUAL >= 20, DP >= 10, MQ >= 30)..."
        bcftools filter -s LowQual -e 'QUAL<20 || DP<10 || MQ<30' variants/${SRR}_raw.vcf > variants/${SRR}_filtered.vcf
    )
    SAMPLE_STATUS=$?
    set -e

    if [ ${SAMPLE_STATUS} -ne 0 ]; then
        echo ">> ERROR: Sample ${SRR} failed (exit code ${SAMPLE_STATUS}) - skipping to next sample."
        FAILED_SAMPLES+=("${SRR}")
        continue
    fi

    echo ">> Completed variant calling for ${SRR}"
done

# ==========================================
# PHASE 4: GLOBAL REPORTING
# ==========================================
echo -e "\n[PHASE 4] COMPILING MULTIQC REPORT..."
module purge
module load MultiQC/1.28-foss-2024a
multiqc qc/ -n final_multiqc_report.html

# ==========================================
# PHASE 5: MULTI-SAMPLE VCF MERGE (FIX #1 — now automated)
# ==========================================
echo -e "\n[PHASE 5] MERGING SAMPLE VCFs INTO A COHORT-LEVEL VCF..."
cd variants
module purge
module load BCFtools/1.21-GCC-13.3.0

FILTERED_VCFS=$(ls *_filtered.vcf 2>/dev/null || true)

if [ -z "${FILTERED_VCFS}" ]; then
    echo "WARNING: No filtered VCFs available - skipping cohort merge."
elif [ $(echo "${FILTERED_VCFS}" | wc -l) -eq 1 ]; then
    echo "Only one sample succeeded - copying its filtered VCF as the 'cohort' VCF (no merge needed)."
    cp ${FILTERED_VCFS} merged_cohort.vcf
else
    # bcftools merge needs bgzipped + tabix-indexed inputs
    for vcf in ${FILTERED_VCFS}; do
        bgzip -f -k ${vcf}
        bcftools index -t -f ${vcf}.gz
    done
    bcftools merge $(for vcf in ${FILTERED_VCFS}; do echo "${vcf}.gz"; done) -O v -o merged_cohort.vcf
    echo "Cohort VCF written to ${WORKDIR}/variants/merged_cohort.vcf"
fi
cd ${WORKDIR}

echo "=========================================================="
if [ ${#FAILED_SAMPLES[@]} -gt 0 ]; then
    echo " PIPELINE COMPLETED WITH ${#FAILED_SAMPLES[@]} FAILED SAMPLE(S):"
    printf '   %s\n' "${FAILED_SAMPLES[@]}"
    echo " These samples were excluded from merged_cohort.vcf."
else
    echo " PIPELINE SUCCESSFULLY COMPLETED! "
fi
echo " Filtered VCF files stored in: ${WORKDIR}/variants/"
echo " Cohort VCF (if >=1 sample succeeded): ${WORKDIR}/variants/merged_cohort.vcf"
echo "=========================================================="

###############################################################################
# COMPREHENSIVE EDUCATIONAL COMMENTARY & METHODOLOGY (carried over)
#
# 1. PCR Duplicate Removal (samtools markdup)
# ---------------------------------------------------------------------------
# During PCR amplification, identical DNA fragments can be over-represented.
# If an amplification error occurs early in PCR, it can mimic a biological
# mutation. Removing duplicate reads (-r flag in markdup) ensures that each
# candidate variant is supported by independent biological fragments.
#
# 2. Robust Multi-Metric Variant Filtering
# ---------------------------------------------------------------------------
# This pipeline enforces three strict criteria:
#   - QUAL < 20 : Removes low-confidence variant calls.
#   - DP < 10   : Eliminates calls backed by fewer than 10 reads (low depth).
#   - MQ < 30   : Excludes reads that map ambiguously to repetitive regions.
#
# 3. Read Group Metadata (-R flag in BWA-MEM)
# ---------------------------------------------------------------------------
# Adding explicit Read Group tags (@RG) embeds sample identity metadata
# directly inside the BAM header — essential for multi-sample analysis,
# IGV visualization, and workflow provenance.
#
# 4. Epidemiological Outbreak Tracing Application
# ---------------------------------------------------------------------------
# The merged cohort VCF records SNPs relative to the reference strain.
# Downstream tools (e.g., snp-dists or IQ-TREE) can compare these SNP
# profiles across isolates to construct phylogenetic trees. A shared SNP
# count alone does not confirm an outbreak link — interpretation must also
# consider reference selection, filtering criteria, missing data,
# recombination, organism-specific thresholds, sampling context, and
# epidemiological evidence.
###############################################################################
