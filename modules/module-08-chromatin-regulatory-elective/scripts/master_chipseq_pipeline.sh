#!/bin/bash
#SBATCH --job-name=master_chipseq
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=24:00:00
#SBATCH --output=master_chipseq_%j.out
#SBATCH --error=master_chipseq_%j.err

# ==============================================================================
# NEW PIPELINE — CHIP-SEQ / ATAC-SEQ PEAK CALLING
# See this module's README "Relationship to existing repo content" - no
# pipeline in this repository touched chromatin/regulatory assays before
# this one. This fills that gap the same way this repository's single-cell
# RNA-seq and amplicon metagenomics pipelines were built: a brand-new
# pipeline, not a correction of an existing script.
#
# ONE SCRIPT, TWO MODES (ASSAY_TYPE):
#   ChIP — protein-DNA binding (a transcription factor or histone mark).
#          MACS2 is run WITH an input/IgG control when CONTROL_SRR is set
#          (the field-standard way to model ChIP-seq's own background:
#          open/accessible chromatin and sequencing bias pull down some
#          background DNA even with a "clean" antibody, and a matched
#          control is what lets MACS2 tell a real binding peak apart from
#          that background rather than guessing from the treatment track
#          alone).
#   ATAC — chromatin accessibility (Tn5 transposase tagmentation). No input
#          control exists for ATAC-seq by design (there is no antibody to
#          have a "no-antibody" background counterpart for) - background is
#          instead handled by MACS2's own ATAC-specific parameters
#          (--nomodel --shift -75 --extsize 150), which correct for Tn5's
#          own well-characterized insertion-site bias instead of trying to
#          model fragment length the way ChIP-seq's default does.
#
# TEACHING NOTE — what happens when CONTROL_SRR is left unset for ChIP mode:
#   MACS2 will still run and call peaks from the treatment track alone. This
#   is NOT equivalent to running with a real input control - see this
#   module's README "Worked-example dataset" for a real, disclosed case
#   where a published CTCF ChIP-seq study deposited no input control at all,
#   and Tier B Script 1's explicit flag for exactly this condition rather
#   than silently treating a no-control run as equivalent to a controlled one.
# ==============================================================================

set -e
set -o pipefail

# ==========================================
# CONFIGURATION
# ==========================================
ASSAY_TYPE="${ASSAY_TYPE:-ChIP}"   # ChIP | ATAC
CONTROL_SRR="${CONTROL_SRR:-}"     # ChIP mode only; leave unset to run without an input control (see TEACHING NOTE above)

case "${ASSAY_TYPE}" in
    ChIP|ATAC) ;;
    *) echo "CRITICAL ERROR: ASSAY_TYPE must be 'ChIP' or 'ATAC', got '${ASSAY_TYPE}'."; exit 1 ;;
esac

# Default reference: human GRCh38 (Ensembl release 110) - same source this
# repository's single-cell pipeline uses, for consistency.
FASTA_URL="${FASTA_URL:-https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz}"
GTF_URL="${GTF_URL:-https://ftp.ensembl.org/pub/release-110/gtf/homo_sapiens/Homo_sapiens.GRCh38.110.gtf.gz}"

echo "=========================================================="
echo " INITIATING CHROMATIN PIPELINE (ASSAY_TYPE=${ASSAY_TYPE}) "
echo " Input control: ${CONTROL_SRR:-<none set>}"
echo "=========================================================="

# ==========================================
# PHASE 1: DIRECTORY SETUP
# ==========================================
WORKDIR="/scratch/$(whoami)/master_chipseq_pipeline_${ASSAY_TYPE}"
mkdir -p ${WORKDIR}/{ref,raw_reads,trimmed_reads,alignment,macs2_output,annotation,qc,logs}
cd ${WORKDIR}

if [ ! -f "srr_list.txt" ]; then
    echo "WARNING: srr_list.txt not found. Generating default test list..."
    echo -e "SRR11092056\nSRR11092057" > srr_list.txt
fi
sed -i 's/\r$//' srr_list.txt

# ==========================================
# PHASE 2: REFERENCE GENOME PREPARATION
# ==========================================
echo -e "\n[PHASE 2] PREPARING REFERENCE GENOME, BOWTIE2 INDEX, AND GENE-FEATURES BED..."
cd ref

FASTA_FILE="genome.fna"
GTF_FILE="annotation.gtf"

if [ ! -s "${FASTA_FILE}" ]; then
    echo "Downloading FASTA..."
    wget -qO- ${FASTA_URL} | gunzip -c > ${FASTA_FILE}
    if [ ! -s "${FASTA_FILE}" ]; then
        echo "CRITICAL ERROR: ${FASTA_FILE} downloaded empty/corrupt. Aborting."
        rm -f "${FASTA_FILE}"
        exit 1
    fi
fi

if [ ! -s "${GTF_FILE}" ]; then
    echo "Downloading GTF..."
    wget -qO- ${GTF_URL} | gunzip -c > ${GTF_FILE}
    if [ ! -s "${GTF_FILE}" ]; then
        echo "CRITICAL ERROR: ${GTF_FILE} downloaded empty/corrupt. Aborting."
        rm -f "${GTF_FILE}"
        exit 1
    fi
fi

if [ ! -f "bt2_index.1.bt2" ]; then
    echo "Building Bowtie2 index (this takes a while for a full human genome)..."
    module purge
    module load Bowtie2/2.5.4-GCC-13.2.0
    bowtie2-build --threads 8 ${FASTA_FILE} bt2_index
fi

# Gene-features BED for Tier B's peak-annotation step: one row per gene,
# (chrom, start, end, strand, gene_name), 0-based half-open per BED
# convention (GTF's start is 1-based, hence the -1).
if [ ! -s "gene_features.bed" ]; then
    echo "Building gene_features.bed from ${GTF_FILE}..."
    awk -F'\t' '$3 == "gene" {
        match($9, /gene_name "([^"]+)"/, arr);
        name = (arr[1] != "" ? arr[1] : "NA");
        print $1"\t"($4-1)"\t"$5"\t"name"\t.\t"$7
    }' ${GTF_FILE} > gene_features.bed
fi
cd ${WORKDIR}

# ==========================================
# PHASE 3: OPTIONAL INPUT-CONTROL ALIGNMENT (ChIP mode only, run ONCE)
# ==========================================
CONTROL_BAM=""
if [ "${ASSAY_TYPE}" == "ChIP" ] && [ -n "${CONTROL_SRR}" ]; then
    echo -e "\n[PHASE 3] ALIGNING INPUT CONTROL (${CONTROL_SRR})..."
    if [ ! -f "alignment/${CONTROL_SRR}_dedup.bam" ]; then
        module purge
        module load SRA-Toolkit/3.2.0-gompi-2024a
        prefetch ${CONTROL_SRR} -O raw_reads/
        fasterq-dump raw_reads/${CONTROL_SRR} --split-files --outdir raw_reads --threads 4

        module purge
        module load fastp/0.23.4-GCC-13.2.0
        fastp -i raw_reads/${CONTROL_SRR}_1.fastq -I raw_reads/${CONTROL_SRR}_2.fastq \
              -o trimmed_reads/${CONTROL_SRR}_1_clean.fastq -O trimmed_reads/${CONTROL_SRR}_2_clean.fastq \
              --thread 8 --html qc/${CONTROL_SRR}_fastp.html 2> qc/${CONTROL_SRR}_fastp.log

        module purge
        module load Bowtie2/2.5.4-GCC-13.2.0 SAMtools/1.18-GCC-12.3.0
        bowtie2 -p 16 -x ref/bt2_index \
                -1 trimmed_reads/${CONTROL_SRR}_1_clean.fastq -2 trimmed_reads/${CONTROL_SRR}_2_clean.fastq 2> qc/${CONTROL_SRR}_bowtie2.log | \
        samtools sort -@ 4 -o alignment/${CONTROL_SRR}_sorted.bam -
        samtools collate -@ 4 -o alignment/${CONTROL_SRR}_collate.bam alignment/${CONTROL_SRR}_sorted.bam
        samtools fixmate -@ 4 -m alignment/${CONTROL_SRR}_collate.bam alignment/${CONTROL_SRR}_fixmate.bam
        samtools sort -@ 4 -o alignment/${CONTROL_SRR}_possort.bam alignment/${CONTROL_SRR}_fixmate.bam
        samtools markdup -@ 4 alignment/${CONTROL_SRR}_possort.bam alignment/${CONTROL_SRR}_dedup.bam
        samtools index alignment/${CONTROL_SRR}_dedup.bam
    fi
    CONTROL_BAM="alignment/${CONTROL_SRR}_dedup.bam"
    echo "Input control ready: ${CONTROL_BAM}"
elif [ "${ASSAY_TYPE}" == "ChIP" ]; then
    echo -e "\n[PHASE 3] No CONTROL_SRR set - proceeding WITHOUT an input control for ChIP mode."
    echo "This is a real limitation, not a formality - see this pipeline's header TEACHING NOTE."
fi

# ==========================================
# PHASE 4: PER-SAMPLE ALIGNMENT + PEAK CALLING LOOP
# ==========================================
echo -e "\n[PHASE 4] EXECUTING PER-SAMPLE ALIGNMENT AND PEAK CALLING..."

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

        # Step C: Alignment (Bowtie2) + duplicate marking
        # Duplicate marking matters more here than in a bulk RNA-seq context:
        # peak calling directly counts read pileup height, and PCR
        # duplicates from library amplification inflate that pileup without
        # adding real signal - the same samtools collate/fixmate/sort/markdup
        # chain this repository's Module 4 SNP pipeline already established.
        echo ">> Aligning with Bowtie2..."
        module purge
        module load Bowtie2/2.5.4-GCC-13.2.0 SAMtools/1.18-GCC-12.3.0
        bowtie2 -p 16 -x ref/bt2_index \
                -1 trimmed_reads/${SRR}_1_clean.fastq -2 trimmed_reads/${SRR}_2_clean.fastq 2> qc/${SRR}_bowtie2.log | \
        samtools sort -@ 4 -o alignment/${SRR}_sorted.bam -
        samtools collate -@ 4 -o alignment/${SRR}_collate.bam alignment/${SRR}_sorted.bam
        samtools fixmate -@ 4 -m alignment/${SRR}_collate.bam alignment/${SRR}_fixmate.bam
        samtools sort -@ 4 -o alignment/${SRR}_possort.bam alignment/${SRR}_fixmate.bam
        samtools markdup -@ 4 alignment/${SRR}_possort.bam alignment/${SRR}_dedup.bam
        samtools index alignment/${SRR}_dedup.bam
        rm -f alignment/${SRR}_collate.bam alignment/${SRR}_fixmate.bam alignment/${SRR}_possort.bam

        # Step D: Peak calling (MACS2) - parameters differ by ASSAY_TYPE
        echo ">> Calling peaks with MACS2 (${ASSAY_TYPE} mode)..."
        module purge
        module load MACS2/2.2.9.1-foss-2022a
        if [ "${ASSAY_TYPE}" == "ChIP" ]; then
            CONTROL_FLAG=""
            [ -n "${CONTROL_BAM}" ] && CONTROL_FLAG="-c ${CONTROL_BAM}"
            macs2 callpeak -t alignment/${SRR}_dedup.bam ${CONTROL_FLAG} \
                  -f BAMPE -g hs -n ${SRR} --outdir macs2_output/${SRR} \
                  2> qc/${SRR}_macs2.log
        else
            # ATAC-specific: --nomodel/--shift/--extsize correct for Tn5
            # insertion bias instead of MACS2's default fragment-length
            # model, which assumes a ChIP-seq-style sonication size
            # distribution that does not apply to Tn5 tagmentation.
            macs2 callpeak -t alignment/${SRR}_dedup.bam \
                  -f BAMPE -g hs -n ${SRR} --outdir macs2_output/${SRR} \
                  --nomodel --shift -75 --extsize 150 \
                  2> qc/${SRR}_macs2.log
        fi

        if [ ! -s "macs2_output/${SRR}/${SRR}_peaks.narrowPeak" ]; then
            echo ">> ERROR: MACS2 produced no peaks file for ${SRR}."
            exit 1
        fi
        N_PEAKS=$(wc -l < macs2_output/${SRR}/${SRR}_peaks.narrowPeak)
        echo ">> ${N_PEAKS} peaks called for ${SRR}."

        # Step E: Annotate peaks against nearest gene (BEDTools)
        echo ">> Annotating peaks against nearest gene feature (BEDTools)..."
        module purge
        module load BEDTools/2.31.0-GCC-12.3.0
        sort -k1,1 -k2,2n macs2_output/${SRR}/${SRR}_peaks.narrowPeak > macs2_output/${SRR}/${SRR}_peaks.sorted.bed
        sort -k1,1 -k2,2n ref/gene_features.bed > ref/gene_features.sorted.bed
        bedtools closest -d \
                 -a macs2_output/${SRR}/${SRR}_peaks.sorted.bed \
                 -b ref/gene_features.sorted.bed \
                 > annotation/${SRR}_peaks_annotated.bed
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
# PHASE 5: GLOBAL REPORTING
# ==========================================
echo -e "\n[PHASE 5] COMPILING MULTIQC REPORT..."
module purge
module load MultiQC/1.28-foss-2024a
multiqc qc/ -n final_chromatin_${ASSAY_TYPE}_report.html

echo "=========================================================="
if [ ${#FAILED_SAMPLES[@]} -gt 0 ]; then
    echo " PIPELINE COMPLETED WITH ${#FAILED_SAMPLES[@]} FAILED SAMPLE(S):"
    printf '   %s\n' "${FAILED_SAMPLES[@]}"
else
    echo " PIPELINE SUCCESSFULLY COMPLETED! "
fi
echo " Peaks: macs2_output/<SRR>/<SRR>_peaks.narrowPeak"
echo " Peak-to-gene annotation: annotation/<SRR>_peaks_annotated.bed"
echo "=========================================================="

###############################################################################
# METHODOLOGY NOTES
# ---------------------------------------------------------------------------
# 1. Why Bowtie2, not BWA or HISAT2: Bowtie2 is a standard end-to-end short-
#    read aligner with no splice-awareness assumption, matching ChIP-seq/
#    ATAC-seq reads (genomic DNA fragments, not spliced transcripts - HISAT2's
#    splice model would be actively wrong here) and no need for BWA-MEM's
#    long-read/split-alignment features this data doesn't require.
# 2. Why duplicate marking happens on genomic coordinates, not UMIs: this
#    dataset (like most ChIP-seq/ATAC-seq) has no UMI, so samtools markdup's
#    position-based duplicate definition (identical 5' mapping coordinates
#    for both mates) is what's available - a real, standard limitation
#    (two independent fragments that truly started at the same position by
#    chance are indistinguishable from PCR duplicates without a UMI), not a
#    shortcut specific to this script.
# 3. -f BAMPE (not BAM): MACS2 is told these are paired-end fragments so it
#    uses the actual insert size from each read pair instead of extending
#    single reads by an estimated fragment length - more accurate whenever
#    paired-end data is actually available, as it is for this pipeline's
#    worked example.
# 4. bedtools closest -d: reports the nearest gene feature AND the distance
#    to it (0 if the peak overlaps the gene) for every peak, not just peaks
#    that directly overlap a gene - most regulatory peaks (enhancers) are
#    NOT inside the gene they regulate, so an overlap-only annotation would
#    silently drop most of the biology this pipeline is trying to capture.
###############################################################################
