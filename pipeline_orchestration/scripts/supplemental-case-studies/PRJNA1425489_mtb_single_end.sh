#!/bin/bash
#SBATCH --job-name=PRJNA1425489_mtb_SE
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=PRJNA1425489_SE_%j.out
#SBATCH --error=PRJNA1425489_SE_%j.err

# ==============================================================================
# CORRECTED VERSION of PRJNA1425489_mtb_single_end.sh
# See pipeline_appraisal_2026-08-31.md for the full review. Fixes applied:
#
#   1. [CRITICAL] Phase 4 featureCounts no longer uses -p/-B/-C. Those flags
#      require paired-end fragments; this script aligns single-end reads
#      (`hisat2 -U`), so the BAMs are single-end. Keeping -p/-B/-C causes
#      featureCounts to error out ("no paired-end reads detected") after all
#      12 samples have already been downloaded, trimmed, and aligned.
#   2. WORKDIR changed to .../PRJNA1425489_rnaseq_SE (was identical to the
#      paired-end PRJNA1425489_mtb_pipeline.sh's WORKDIR). Running both
#      scripts against the same directory lets one overwrite the other's
#      reference index, BAMs, and counts.
#   3. Reference downloads now verify the file is non-empty before treating
#      it as "already present" on reruns, so a truncated wget/gunzip doesn't
#      silently pass the existence check forever.
#   4. The per-sample loop no longer aborts the whole batch on one failure:
#      each sample runs in its own strict subshell, failures are logged and
#      skipped, and a summary of any failed SRRs is printed at the end.
# ==============================================================================

# ORCHESTRATION COPY -- see pipeline_orchestration/README.md. Differs from
# the original module script only by RUN_TAG-suffixed WORKDIR and
# PIPELINE_SCRIPT_DIR (so it reuses that module's load_modules.sh).

set -e
set -o pipefail

RUN_TAG="${RUN_TAG:-$(date -u +%Y%m%dT%H%M%SZ)}"
export PIPELINE_SCRIPT_DIR="${PIPELINE_SCRIPT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")/../../../modules/supplemental-case-studies/scripts" && pwd)}"
source "${PIPELINE_SCRIPT_DIR}/load_modules.sh"

echo "=========================================================="
echo " INITIATING SINGLE-END RNA-SEQ PIPELINE FOR PRJNA1425489 "
echo "=========================================================="

# ==========================================
# PHASE 1: DIRECTORY SETUP
# ==========================================
# FIX #2: renamed with _SE suffix so this never shares a WORKDIR with the
# paired-end PRJNA1425489_mtb_pipeline.sh.
WORKDIR="/scratch/$(whoami)/PRJNA1425489_rnaseq_SE_${RUN_TAG}"
mkdir -p ${WORKDIR}/{ref,raw_reads,qc,trimmed_reads,alignment,counts,logs}
cd ${WORKDIR}

# Create sample list for the 12 PRJNA1425489 samples
if [ ! -f "srr_list.txt" ]; then
    echo "Creating srr_list.txt for PRJNA1425489..."
    echo -e "SRR37277506\nSRR37277507\nSRR37277508\nSRR37277509\nSRR37277510\nSRR37277511\nSRR37277512\nSRR37277513\nSRR37277514\nSRR37277515\nSRR37277516\nSRR37277517" > srr_list.txt
fi
sed -i 's/\r$//' srr_list.txt

# ==========================================
# PHASE 2: REFERENCE GENOME PREPARATION
# ==========================================
echo -e "\n[PHASE 2] PREPARING M. TUBERCULOSIS REFERENCE GENOME..."
cd ref

FASTA_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/195/955/GCF_000195955.2_ASM19595v2/GCF_000195955.2_ASM19595v2_genomic.fna.gz"
GFF_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/195/955/GCF_000195955.2_ASM19595v2/GCF_000195955.2_ASM19595v2_genomic.gff.gz"

FASTA_FILE="genome.fna"
GFF_FILE="annotation.gff"
INDEX_PREFIX="genome_index"

# FIX #3: -s checks the file exists AND is non-empty, so a prior truncated
# download gets re-fetched instead of being treated as complete.
if [ ! -s "${FASTA_FILE}" ]; then
    echo "Downloading FASTA..."
    wget -qO- ${FASTA_URL} | gunzip -c > ${FASTA_FILE}
    if [ ! -s "${FASTA_FILE}" ]; then
        echo "CRITICAL ERROR: ${FASTA_FILE} downloaded empty/corrupt. Aborting."
        rm -f "${FASTA_FILE}"
        exit 1
    fi
fi

if [ ! -s "${GFF_FILE}" ]; then
    echo "Downloading GFF..."
    wget -qO- ${GFF_URL} | gunzip -c > ${GFF_FILE}
    if [ ! -s "${GFF_FILE}" ]; then
        echo "CRITICAL ERROR: ${GFF_FILE} downloaded empty/corrupt. Aborting."
        rm -f "${GFF_FILE}"
        exit 1
    fi
fi

if [ ! -f "${INDEX_PREFIX}.1.ht2" ]; then
    echo "Building HISAT2 index..."
    module purge
    load_python || exit 1
    load_hisat2 || exit 1
    hisat2-build -p 8 ${FASTA_FILE} ${INDEX_PREFIX}
fi
cd ${WORKDIR}

# ==========================================
# PHASE 3: SAMPLE PROCESSING LOOP
# ==========================================
echo -e "\n[PHASE 3] EXECUTING SAMPLE PROCESSING LOOP..."

# FIX #4: track failures instead of letting one bad sample kill the batch.
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
        load_sra_toolkit "SRA-Toolkit/3.0.3-gompi-2022a" || exit 1
        prefetch ${SRR} -O raw_reads/
        fasterq-dump raw_reads/${SRR} --force --outdir raw_reads --threads 4

        # Step B: Pre-trimming QC (SINGLE END)
        echo ">> Running Initial FastQC..."
        module purge
        load_fastqc || exit 1
        fastqc raw_reads/${SRR}.fastq -o qc -t 8 -q

        # Step C: Quality & Adapter Trimming (SINGLE END)
        echo ">> Trimming with fastp..."
        module purge
        load_fastp || exit 1
        fastp -i raw_reads/${SRR}.fastq \
              -o trimmed_reads/${SRR}_clean.fastq \
              --thread 8 --html qc/${SRR}_fastp.html 2> qc/${SRR}_fastp.log

        # Step D: Post-trimming QC (SINGLE END)
        echo ">> Running Post-trim QC..."
        module purge
        load_fastqc || exit 1
        fastqc trimmed_reads/${SRR}_clean.fastq -o qc -t 8 -q

        # Step E: Alignment & BAM Sorting (SINGLE END)
        # Note: HISAT2 uses -U for single-end instead of -1 and -2
        echo ">> Aligning and sorting BAM..."
        module purge
        load_python   || exit 1
        load_hisat2   || exit 1
        load_samtools || exit 1

        hisat2 -p 8 --no-spliced-alignment -x ref/${INDEX_PREFIX} \
               -U trimmed_reads/${SRR}_clean.fastq | \
        samtools sort -@ 4 -o alignment/${SRR}_sorted.bam -

        # Step F: Indexing BAM
        echo ">> Indexing BAM..."
        samtools index alignment/${SRR}_sorted.bam
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
# PHASE 4: GLOBAL FEATURE QUANTIFICATION
# ==========================================
echo -e "\n[PHASE 4] GENERATING GENE EXPRESSION MATRIX..."
module purge
load_subread || exit 1

BAM_FILES=$(ls alignment/*_sorted.bam 2>/dev/null || true)

if [ -z "${BAM_FILES}" ]; then
    echo "CRITICAL ERROR: No aligned BAM files were produced. Aborting quantification."
    exit 1
fi

# FIX #1: removed -p -B -C. Those flags assume paired-end fragments; this
# pipeline aligns with `hisat2 -U` (single-end), so featureCounts would
# reject these BAMs with "no paired-end reads detected" if -p were kept.
featureCounts -T 8 \
              -t gene -g locus_tag \
              -a ref/${GFF_FILE} \
              -o counts/gene_counts_raw.txt \
              ${BAM_FILES}

# ==========================================
# PHASE 5: MATRIX CLEANUP
# ==========================================
echo -e "\n[PHASE 5] CLEANING MATRIX FOR DESeq2..."
cd counts
cut -f1,7- gene_counts_raw.txt | sed '1d' > gene_counts_clean.tsv
sed -i 's/alignment\///g; s/_sorted.bam//g' gene_counts_clean.tsv
cd ${WORKDIR}

# ==========================================
# PHASE 6: GLOBAL REPORTING
# ==========================================
echo -e "\n[PHASE 6] COMPILING MULTIQC REPORT..."
module purge
load_python  || exit 1   # must be loaded before MultiQC -- see load_modules.sh
load_multiqc || exit 1
multiqc qc/ -n final_multiqc_report.html

echo "=========================================================="
if [ ${#FAILED_SAMPLES[@]} -gt 0 ]; then
    echo " PIPELINE COMPLETED WITH ${#FAILED_SAMPLES[@]} FAILED SAMPLE(S):"
    printf '   %s\n' "${FAILED_SAMPLES[@]}"
    echo " Re-run those specific SRR IDs individually if you need complete data."
else
    echo " PIPELINE SUCCESSFULLY COMPLETED! "
fi
echo "=========================================================="
