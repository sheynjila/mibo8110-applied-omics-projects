#!/bin/bash
#SBATCH --job-name=master_scrnaseq
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=24:00:00
#SBATCH --output=master_scrnaseq_%j.out
#SBATCH --error=master_scrnaseq_%j.err

# ==============================================================================
# NEW PIPELINE — SINGLE-CELL RNA-SEQ (10x Genomics-style, STARsolo)
# See domain_coverage_audit_2026-08-31.md — this fills the single-cell RNA-seq
# gap that was entirely missing from the original 22-script collection.
#
# WHY THIS IS A SEPARATE PIPELINE FROM master_rnaseq_pipeline_consolidated.sh:
#   Bulk RNA-seq (HISAT2 + featureCounts) treats each sample as one pooled
#   read population. Single-cell data instead has one cDNA read (R2) paired
#   with a cell-barcode + UMI read (R1) per read pair, and needs a
#   barcode-aware aligner/quantifier (here: STARsolo) that demultiplexes
#   reads into per-cell counts as part of alignment.
#
# SCOPE (per your answers):
#   - Generic 10x Chromium-style teaching example (v3 chemistry: 16bp cell
#     barcode + 12bp UMI). Swap CHEMISTRY/CB_LEN/UMI_LEN below for v2 data
#     (16bp CB + 10bp UMI) or a different platform.
#   - Stops at the gene-by-cell count matrix. STARsolo's output
#     (matrix.mtx / features.tsv / barcodes.tsv) is already in standard 10x
#     format, directly loadable in R via Seurat::Read10X() or in Python via
#     scanpy.read_10x_mtx() — no separate export step needed.
#
# REQUIRED MANUAL SETUP BEFORE FIRST RUN (cannot be scripted/downloaded
# automatically — see Phase 0 below):
#   1. A 10x cell-barcode whitelist file matching your chemistry version.
#      These ship with Cell Ranger's reference bundles and are also
#      available from 10x Genomics support; place the path in
#      WHITELIST_PATH.
#   2. Confirm the STAR module name/version available on your cluster
#      (`module spider STAR`) — the version below is a placeholder.
# ==============================================================================

# ORCHESTRATION COPY -- see pipeline_orchestration/README.md. Differs from
# the original module script only by RUN_TAG-suffixed WORKDIR and
# PIPELINE_SCRIPT_DIR (so it reuses that module's load_modules.sh).

set -e
set -o pipefail

RUN_TAG="${RUN_TAG:-$(date -u +%Y%m%dT%H%M%SZ)}"
export PIPELINE_SCRIPT_DIR="${PIPELINE_SCRIPT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")/../../../modules/module-07-single-cell-elective/scripts" && pwd)}"
source "${PIPELINE_SCRIPT_DIR}/load_modules.sh"

# ==========================================
# CONFIGURATION
# ==========================================
CHEMISTRY="${CHEMISTRY:-10x_v3}"   # 10x_v3 (16bp CB + 12bp UMI) | 10x_v2 (16bp CB + 10bp UMI)
case "${CHEMISTRY}" in
  10x_v3) CB_LEN=16; UMI_LEN=12 ;;
  10x_v2) CB_LEN=16; UMI_LEN=10 ;;
  *) echo "CRITICAL ERROR: Unknown CHEMISTRY '${CHEMISTRY}'. Use 10x_v3 or 10x_v2."; exit 1 ;;
esac

# Path to the 10x barcode whitelist matching CHEMISTRY (see Phase 0 notice).
WHITELIST_PATH="${WHITELIST_PATH:-/path/to/10x_whitelists/${CHEMISTRY}_barcodes.txt}"

# Default reference: human GRCh38.p14 (swap for mouse/other via these two URLs).
FASTA_URL="${FASTA_URL:-https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.fna.gz}"
GTF_URL="${GTF_URL:-https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.gtf.gz}"

# cDNA (R2) read length minus 1, per STAR's --sjdbOverhang convention.
# 10x v3 R2 is typically 91bp -> overhang 90. Adjust for your run's actual length.
SJDB_OVERHANG="${SJDB_OVERHANG:-90}"

echo "=========================================================="
echo " INITIATING SINGLE-CELL RNA-SEQ PIPELINE (STARsolo, ${CHEMISTRY}) "
echo "=========================================================="

# ==========================================
# PHASE 0: PRE-FLIGHT CHECKS
# ==========================================
if [ ! -s "${WHITELIST_PATH}" ]; then
    echo "CRITICAL ERROR: Barcode whitelist not found at ${WHITELIST_PATH}."
    echo "Set WHITELIST_PATH to a valid 10x barcode whitelist file for ${CHEMISTRY}"
    echo "before running this pipeline (see header comment for where to obtain one)."
    exit 1
fi

# ==========================================
# PHASE 1: DIRECTORY SETUP
# ==========================================
WORKDIR="/scratch/$(whoami)/master_scrnaseq_pipeline_${RUN_TAG}"
mkdir -p ${WORKDIR}/{ref,raw_reads,qc,alignment,logs}
cd ${WORKDIR}

if [ ! -f "srr_list.txt" ]; then
    echo "WARNING: srr_list.txt not found. Generating default test list..."
    echo -e "SRR11092056\nSRR11092057" > srr_list.txt
fi
sed -i 's/\r$//' srr_list.txt

# ==========================================
# PHASE 2: REFERENCE GENOME PREPARATION
# ==========================================
echo -e "\n[PHASE 2] PREPARING REFERENCE GENOME & STAR INDEX..."
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

if [ ! -f "star_index/SAindex" ]; then
    echo "Building STAR genome index (this takes a while for large genomes)..."
    module purge
    load_star || exit 1
    mkdir -p star_index
    STAR --runMode genomeGenerate \
         --runThreadN 8 \
         --genomeDir star_index \
         --genomeFastaFiles ${FASTA_FILE} \
         --sjdbGTFfile ${GTF_FILE} \
         --sjdbOverhang ${SJDB_OVERHANG}
fi
cd ${WORKDIR}

# ==========================================
# PHASE 3: SAMPLE PROCESSING LOOP
# ==========================================
echo -e "\n[PHASE 3] EXECUTING PER-SAMPLE STARsolo QUANTIFICATION..."

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
        # SRA convention for 10x runs: _1 = cell barcode + UMI (short read),
        # _2 = cDNA (biological read). Confirm this matches your dataset's
        # actual read layout before trusting it blindly (check read lengths).
        echo ">> Downloading from SRA..."
        module purge
        load_sra_toolkit || exit 1
        prefetch ${SRR} -O raw_reads/
        fasterq-dump raw_reads/${SRR} --split-files --outdir raw_reads --threads 4

        # Step B: Raw QC (informational only - no trimming before STARsolo;
        # STARsolo handles barcode/UMI extraction itself and trimming cDNA
        # reads pre-alignment is not standard practice for single-cell data).
        echo ">> Running FastQC..."
        module purge
        load_fastqc || exit 1
        fastqc raw_reads/${SRR}_1.fastq raw_reads/${SRR}_2.fastq -o qc -t 8 -q

        # Step C: STARsolo alignment + per-cell quantification
        echo ">> Running STARsolo..."
        module purge
        load_star || exit 1

        mkdir -p alignment/${SRR}
        STAR --runThreadN 8 \
             --genomeDir ref/star_index \
             --readFilesIn raw_reads/${SRR}_2.fastq raw_reads/${SRR}_1.fastq \
             --soloType CB_UMI_Simple \
             --soloCBwhitelist ${WHITELIST_PATH} \
             --soloCBstart 1 --soloCBlen ${CB_LEN} \
             --soloUMIstart $((CB_LEN + 1)) --soloUMIlen ${UMI_LEN} \
             --soloFeatures Gene \
             --outSAMtype BAM SortedByCoordinate \
             --outFileNamePrefix alignment/${SRR}/

        # Sanity check: confirm a count matrix was actually produced
        if [ ! -f "alignment/${SRR}/Solo.out/Gene/raw/matrix.mtx" ]; then
            echo ">> ERROR: STARsolo did not produce a count matrix for ${SRR}."
            exit 1
        fi
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
# PHASE 4: GLOBAL REPORTING
# ==========================================
echo -e "\n[PHASE 4] COMPILING MULTIQC REPORT..."
module purge
load_python  || exit 1   # must be loaded before MultiQC -- see load_modules.sh
load_multiqc || exit 1
multiqc qc/ alignment/ -n final_scrnaseq_report.html

echo "=========================================================="
if [ ${#FAILED_SAMPLES[@]} -gt 0 ]; then
    echo " PIPELINE COMPLETED WITH ${#FAILED_SAMPLES[@]} FAILED SAMPLE(S):"
    printf '   %s\n' "${FAILED_SAMPLES[@]}"
else
    echo " PIPELINE SUCCESSFULLY COMPLETED! "
fi
echo " Per-sample count matrices: alignment/<SRR>/Solo.out/Gene/{raw,filtered}/"
echo " Load into Seurat with:  Seurat::Read10X(\"alignment/<SRR>/Solo.out/Gene/filtered/\")"
echo " Load into Scanpy with:  scanpy.read_10x_mtx(\"alignment/<SRR>/Solo.out/Gene/filtered/\")"
echo "=========================================================="

###############################################################################
# METHODOLOGY NOTES
# ---------------------------------------------------------------------------
# 1. Why "filtered" vs "raw" matrices: STARsolo emits both a raw matrix
#    (every barcode seen, mostly empty droplets/noise) and a filtered matrix
#    (STARsolo's own cell-calling, analogous to Cell Ranger's). Use the
#    filtered matrix for standard downstream analysis; keep raw around only
#    if you plan to re-run cell calling with a different algorithm (e.g.
#    EmptyDrops) yourself.
# 2. Chemistry mismatches are silent failures: if CB_LEN/UMI_LEN don't match
#    how the reads were actually generated, STARsolo will still run but
#    produce garbage barcodes/near-zero valid cells. Always confirm chemistry
#    version from the study's methods/SRA metadata before trusting defaults.
# 3. This pipeline does not perform ambient RNA correction, doublet removal,
#    or clustering — those are downstream analysis steps (Seurat/Scanpy),
#    intentionally out of scope here per the requested stopping point.
###############################################################################
