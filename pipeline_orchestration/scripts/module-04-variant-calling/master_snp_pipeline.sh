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
# ORCHESTRATION COPY of modules/module-04-variant-calling/scripts/master_snp_pipeline.sh
# See pipeline_orchestration/README.md for what's different and why (short
# version: RUN_TAG-suffixed WORKDIR so concurrent/repeated runs don't
# collide, PIPELINE_SCRIPT_DIR so this copy reuses the original module's
# load_modules.sh). Everything below this block is otherwise identical to
# the original -- if you fix a bug in one, apply the same fix to the other.
# ==============================================================================

set -e
set -o pipefail

RUN_TAG="${RUN_TAG:-$(date -u +%Y%m%dT%H%M%SZ)}"
export PIPELINE_SCRIPT_DIR="${PIPELINE_SCRIPT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")/../../../modules/module-04-variant-calling/scripts" && pwd)}"
source "${PIPELINE_SCRIPT_DIR}/load_modules.sh"

echo "=========================================================="
echo " INITIATING MASTER BACTERIAL VARIANT CALLING PIPELINE "
echo " RUN_TAG: ${RUN_TAG}"
echo "=========================================================="

# ==========================================
# PHASE 1: DIRECTORY & ENVIRONMENT SETUP
# ==========================================
WORKDIR="/scratch/$(whoami)/master_snp_pipeline_${RUN_TAG}"
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
    load_bwa || exit 1
    bwa index ${FASTA_FILE}
fi

if [ ! -f "${FASTA_FILE}.fai" ]; then
    echo "Building SAMtools FASTA index..."
    module purge
    load_samtools || exit 1
    samtools faidx ${FASTA_FILE}
fi
cd ${WORKDIR}

# ==========================================
# PHASE 3: SAMPLE PROCESSING LOOP
# ==========================================
echo -e "\n[PHASE 3] EXECUTING WGS SAMPLE PROCESSING LOOP..."

FAILED_SAMPLES=()

for SRR in $(cat srr_list.txt); do
    echo "----------------------------------------------------------"
    echo " PROCESSING SAMPLE: ${SRR}"
    echo "----------------------------------------------------------"

    set +e
    (
        set -e
        set -o pipefail

        echo ">> Downloading from SRA..."
        module purge
        load_sra_toolkit || exit 1
        prefetch ${SRR} -O raw_reads/
        fasterq-dump raw_reads/${SRR} --split-files --outdir raw_reads --threads 4

        echo ">> Trimming with fastp..."
        module purge
        load_fastp || exit 1
        fastp -i raw_reads/${SRR}_1.fastq -I raw_reads/${SRR}_2.fastq \
              -o trimmed_reads/${SRR}_1_clean.fastq -O trimmed_reads/${SRR}_2_clean.fastq \
              --thread 8 --html qc/${SRR}_fastp.html 2> qc/${SRR}_fastp.log

        echo ">> Aligning DNA reads to reference..."
        module purge
        load_bwa      || exit 1
        load_samtools || exit 1

        bwa mem -t 8 -R "@RG\tID:${SRR}\tSM:${SRR}\tPL:ILLUMINA" ref/${FASTA_FILE} \
                trimmed_reads/${SRR}_1_clean.fastq \
                trimmed_reads/${SRR}_2_clean.fastq | \
        samtools sort -@ 4 -o alignment/${SRR}_unsorted.bam -

        echo ">> Marking and filtering PCR duplicates with SAMtools..."
        samtools collate -@ 4 -O -u alignment/${SRR}_unsorted.bam | \
        samtools fixmate -@ 4 -m -u - - | \
        samtools sort -@ 4 -u - | \
        samtools markdup -@ 4 -r - alignment/${SRR}_dedup.bam

        samtools index alignment/${SRR}_dedup.bam
        rm alignment/${SRR}_unsorted.bam

        echo ">> Calling variants (SNPs and INDELs)..."
        module purge
        load_bcftools || exit 1

        bcftools mpileup -O b -f ref/${FASTA_FILE} alignment/${SRR}_dedup.bam | \
        bcftools call -mv -O v -o variants/${SRR}_raw.vcf

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
load_python  || exit 1
load_multiqc || exit 1
multiqc qc/ -n final_multiqc_report.html

# ==========================================
# PHASE 5: MULTI-SAMPLE VCF MERGE
# ==========================================
echo -e "\n[PHASE 5] MERGING SAMPLE VCFs INTO A COHORT-LEVEL VCF..."
cd variants
module purge
load_bcftools || exit 1

FILTERED_VCFS=$(ls *_filtered.vcf 2>/dev/null || true)

if [ -z "${FILTERED_VCFS}" ]; then
    echo "WARNING: No filtered VCFs available - skipping cohort merge."
elif [ $(echo "${FILTERED_VCFS}" | wc -l) -eq 1 ]; then
    echo "Only one sample succeeded - copying its filtered VCF as the 'cohort' VCF (no merge needed)."
    cp ${FILTERED_VCFS} merged_cohort.vcf
else
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
echo " RUN_TAG: ${RUN_TAG}"
echo " Filtered VCF files stored in: ${WORKDIR}/variants/"
echo " Cohort VCF (if >=1 sample succeeded): ${WORKDIR}/variants/merged_cohort.vcf"
echo "=========================================================="
