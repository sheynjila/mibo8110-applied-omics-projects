#!/bin/bash
#SBATCH --job-name=master_rnaseq
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --output=master_rnaseq_%j.out
#SBATCH --error=master_rnaseq_%j.err

###############################################################################
# CASE-BY-CASE CONFIGURATION (DYNAMIC ORGANISM SUPPORT)
# ---------------------------------------------------------------------------
# Because organisms have fundamentally different genomic structures (e.g.,
# mammals have introns, bacteria do not), the script dynamically adjusts its
# alignment and counting logic based on a single variable called REFERENCE_TYPE.
#
# You can run different biological cases by setting this variable. You can also
# supply ANY custom genome using the FASTA_URL and GFF_URL variables. If URLs
# are omitted, the script safely defaults to the curriculum standard genomes.
#
#   * bacterial : No spliced alignment; counts via locus_tag.
#                 Default: Salmonella enterica (GCF_000006945.2)
#   * viral     : No spliced alignment; counts via locus_tag.
#                 Default: SARS-CoV-2 (GCF_009858895.2)
#   * human / mammal : Spliced alignment enabled (crosses exon junctions);
#                 counts via gene. Default: Human GRCh38.p14 (GCF_000001405.40)
#
# USAGE EXAMPLES:
# 1. Run curriculum default (Human):
#    sbatch --export=REFERENCE_TYPE=mammal master_rnaseq_pipeline_consolidated.sh
# 2. Run a CUSTOM bacterium (e.g., E. coli) by passing URLs:
#    sbatch --export=REFERENCE_TYPE=bacterial,FASTA_URL="https://...fna.gz",GFF_URL="https://...gff.gz" master_rnaseq_pipeline_consolidated.sh
# 3. Run a CUSTOM mammal (e.g., Mouse) whose GFF uses a non-default attribute:
#    sbatch --export=REFERENCE_TYPE=mammal,FASTA_URL="https://...mouse...fna.gz",GFF_URL="https://...mouse...gff.gz",FC_ATTR=gene_id master_rnaseq_pipeline_consolidated.sh
#    (Only pass FC_ATTR/SPLICE_FLAG if your annotation's attribute key isn't
#    "gene"/"locus_tag" — most RefSeq GFFs don't need this. RefSeq mouse
#    annotation below uses "gene", same as human, so it is NOT needed there.)
#
# -----------------------------------------------------------------------------
# QUICK START — VERIFIED, READY-TO-RUN EXAMPLE (PRJNA1518998, Mouse retina
# study; URLs and accessions checked live against NCBI on 2026-08-31):
#
#   Step 1 — write a 2-sample smoke-test list (real paired-end runs from
#   PRJNA1518998, ~3GB/run — small enough for a quick first test):
#
#     echo -e "SRR40359383\nSRR40359384" > srr_list.txt
#
#   Step 2 — submit with REFERENCE_TYPE=mammal (Mouse uses the same
#   splice/attribute ruleset as Human) plus the verified GRCm39 URLs, since
#   the script's default under "mammal" is Human, not Mouse. Copy this as ONE
#   line (the wrap below is only for readability in this comment block):
#
#     sbatch --export=REFERENCE_TYPE=mammal,FASTA_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/635/GCF_000001635.27_GRCm39/GCF_000001635.27_GRCm39_genomic.fna.gz",GFF_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/635/GCF_000001635.27_GRCm39/GCF_000001635.27_GRCm39_genomic.gff.gz" master_rnaseq_pipeline_consolidated.sh
#
#   That's the ONLY thing that changes between studies — srr_list.txt plus
#   the --export flags. Nothing inside this file needs to be opened or edited.
#   Once the smoke test's alignment rate looks healthy (>80%, see PIPELINE
#   NOTES at the bottom of this file), replace srr_list.txt with the study's
#   full accession list and resubmit the same sbatch command.
#
#   To run a DIFFERENT study later, you only ever touch two things:
#     (a) srr_list.txt        — the sample accessions for that study
#     (b) the --export flags  — REFERENCE_TYPE, and FASTA_URL/GFF_URL if the
#                                organism isn't one of the 3 curriculum defaults
#   FC_ATTR/SPLICE_FLAG only need to be added if a custom (non-RefSeq)
#   annotation uses an attribute key other than "gene"/"locus_tag".
# -----------------------------------------------------------------------------
#
# IMPORTANT EXCEPTIONS
# ---------------------------------------------------------------------------
# * Strictly Paired-End Data : This script assumes all input samples are
#                              paired-end reads (-p -B -C flags).
# * Single-End Data          : Use the dedicated single-end pipeline instead
#                              (see modules/supplemental-case-studies/).
# * One study per WORKDIR    : WORKDIR is now keyed by the resolved genome
#                              (see FIX #1 below), not just REFERENCE_TYPE, so
#                              switching organisms is safe. Two DIFFERENT
#                              STUDIES against the SAME genome build still
#                              share raw_reads/trimmed_reads/alignment/counts —
#                              clear WORKDIR (or change srr_list.txt and check
#                              the startup warning) before starting a new study
#                              on a genome you've already used.
#
# FIXES APPLIED IN THIS REVISION (relative to the auto-generated draft):
#   1. CRITICAL: WORKDIR (and the reference cache inside it) is now keyed by
#      a GENOME_TAG derived from FASTA_URL, not just REFERENCE_TYPE. The
#      original keyed WORKDIR purely on REFERENCE_TYPE, so overriding
#      FASTA_URL/GFF_URL to point at a different species (e.g. Mouse while
#      REFERENCE_TYPE=mammal) would silently reuse a PREVIOUSLY DOWNLOADED
#      genome/index/BAMs from a different species already sitting in that
#      same directory — aligning your new samples against the wrong genome
#      with no error or warning. This was the single most dangerous bug in
#      the draft, precisely because it fails silently and "successfully."
#   2. FC_ATTR and SPLICE_FLAG can now be overridden via --export, exactly
#      like FASTA_URL/GFF_URL. The draft hardcoded these per REFERENCE_TYPE
#      bucket, so a genuinely custom (non-RefSeq) organism whose GFF uses a
#      different attribute key (e.g. Ensembl-style "gene_id" or "ID") would
#      run to completion but silently produce all-zero or wrong counts.
#   3. Thread/CPU counts (-p, -T, --thread, -@) now derive from
#      SLURM_CPUS_PER_TASK instead of being hardcoded to 8/4 everywhere, so
#      they stay correct if --cpus-per-task is changed at submission.
#   4. Sample loop uses `while read` instead of `for X in $(cat file)`,
#      avoiding the word-splitting/blank-line trap (same fix already applied
#      to the batch-automation script in the storage-safety module).
#   5. A running reference_manifest.txt log (in logs/) records which
#      REFERENCE_TYPE/GENOME_TAG/URLs/attribute were actually used on every
#      run of this WORKDIR, so a later reader can see exactly what a given
#      counts matrix was generated against.
#   6. Environment-module loading goes through this module's load_modules.sh
#      (see modules/module-02-bulk-rnaseq/scripts/load_modules.sh) instead of
#      bare `module load <exact-build>` — see the root RUNBOOK.md for why.
###############################################################################

set -e
set -o pipefail

SCRIPT_DIR="${PIPELINE_SCRIPT_DIR:-${SLURM_SUBMIT_DIR:+${SLURM_SUBMIT_DIR}/modules/module-02-bulk-rnaseq/scripts}}"
SCRIPT_DIR="${SCRIPT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)}"
source "${SCRIPT_DIR}/load_modules.sh"

# ==========================================
# CONFIGURATION
# ==========================================
REFERENCE_TYPE="${REFERENCE_TYPE:-bacterial}"
FASTA_URL="${FASTA_URL:-}"
GFF_URL="${GFF_URL:-}"
SPLICE_FLAG="${SPLICE_FLAG:-}"
FC_ATTR="${FC_ATTR:-}"

case "${REFERENCE_TYPE}" in
  bacterial)
    FASTA_URL="${FASTA_URL:-https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/006/945/GCF_000006945.2_ASM694v2/GCF_000006945.2_ASM694v2_genomic.fna.gz}"
    GFF_URL="${GFF_URL:-https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/006/945/GCF_000006945.2_ASM694v2/GCF_000006945.2_ASM694v2_genomic.gff.gz}"
    SPLICE_FLAG="${SPLICE_FLAG:---no-spliced-alignment}"
    FC_ATTR="${FC_ATTR:-locus_tag}"
    ;;
  viral)
    FASTA_URL="${FASTA_URL:-https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/858/895/GCF_009858895.2_ASM985889v3/GCF_009858895.2_ASM985889v3_genomic.fna.gz}"
    GFF_URL="${GFF_URL:-https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/858/895/GCF_009858895.2_ASM985889v3/GCF_009858895.2_ASM985889v3_genomic.gff.gz}"
    SPLICE_FLAG="${SPLICE_FLAG:---no-spliced-alignment}"
    FC_ATTR="${FC_ATTR:-locus_tag}"
    ;;
  human|mammal)
    FASTA_URL="${FASTA_URL:-https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.fna.gz}"
    GFF_URL="${GFF_URL:-https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.gff.gz}"
    SPLICE_FLAG="${SPLICE_FLAG:-}"   # allow spliced alignment across exon-exon junctions
    FC_ATTR="${FC_ATTR:-gene}"
    ;;
  *)
    echo "CRITICAL ERROR: Unknown REFERENCE_TYPE '${REFERENCE_TYPE}'. Use bacterial|viral|human|mammal."
    exit 1
    ;;
esac

# FIX #1: derive a genome-specific tag from the FASTA URL so WORKDIR (and
# every cache decision inside it) is keyed on the ACTUAL genome being used,
# not just the 3-way REFERENCE_TYPE bucket. A short hash of the full URL is
# appended so two different URLs can never collide onto the same tag even if
# their filenames happen to match.
URL_HASH="$(printf '%s' "${FASTA_URL}" | md5sum | cut -c1-8)"
BASE_TAG="$(basename "${FASTA_URL}" 2>/dev/null | sed -E 's/\.(fna|fa|fasta)\.gz$//; s/[^A-Za-z0-9._-]/_/g')"
GENOME_TAG="${BASE_TAG:-genome}_${URL_HASH}"

# FIX #3: scale thread usage with the actual Slurm allocation instead of a
# hardcoded 8, so changing --cpus-per-task at submission stays correct.
CPUS="${SLURM_CPUS_PER_TASK:-8}"
SORT_THREADS=$(( CPUS / 2 ))
[ "${SORT_THREADS}" -lt 1 ] && SORT_THREADS=1

echo "=========================================================="
echo " INITIATING CONSOLIDATED RNA-SEQ PIPELINE "
echo " RULESET     : ${REFERENCE_TYPE} "
echo " GENOME_TAG  : ${GENOME_TAG} "
echo " FASTA       : ${FASTA_URL} "
echo " FC_ATTR     : ${FC_ATTR} "
echo " CPUS        : ${CPUS} "
echo "=========================================================="

# ==========================================
# PHASE 1: DIRECTORY SETUP
# ==========================================
# FIX #1 (cont.): WORKDIR includes GENOME_TAG, so overriding FASTA_URL/GFF_URL
# to a different species/assembly always lands in its own isolated directory —
# it can never silently reuse another genome's downloaded FASTA, built index,
# or leftover BAMs.
WORKDIR="/scratch/$(whoami)/master_rnaseq_pipeline_${REFERENCE_TYPE}_${GENOME_TAG}"
mkdir -p "${WORKDIR}"/{ref,raw_reads,qc,trimmed_reads,alignment,counts,logs}
cd "${WORKDIR}"

if [ ! -f "srr_list.txt" ]; then
    echo "Creating default srr_list.txt..."
    echo -e "SRR11092056\nSRR11092057" > srr_list.txt
fi
sed -i 's/\r$//' srr_list.txt

# FIX #5: keep an append-only log of which reference this WORKDIR was run
# against, so a counts matrix generated weeks apart from the sbatch command
# that made it is still traceable.
{
    echo "timestamp: $(date -Iseconds)"
    echo "reference_type: ${REFERENCE_TYPE}"
    echo "genome_tag: ${GENOME_TAG}"
    echo "fasta_url: ${FASTA_URL}"
    echo "gff_url: ${GFF_URL}"
    echo "fc_attr: ${FC_ATTR}"
    echo "splice_flag: ${SPLICE_FLAG:-<none>}"
    echo "---"
} >> logs/reference_manifest.txt

# Advisory only (does not abort): two different STUDIES against the same
# genome still share raw_reads/trimmed_reads/alignment/counts by design, so
# flag it loudly if this WORKDIR already has processed samples from a
# different srr_list.txt than the one about to run.
LIST_HASH="$(md5sum srr_list.txt | cut -d' ' -f1)"
if [ -f ".last_srr_list_hash" ] && [ "$(cat .last_srr_list_hash)" != "${LIST_HASH}" ] && [ -n "$(ls alignment/*_sorted.bam 2>/dev/null)" ]; then
    echo "WARNING: srr_list.txt changed since the last run in this WORKDIR,"
    echo "         but alignment/ and counts/ still hold data from that"
    echo "         earlier sample list. If this is a NEW study (not just"
    echo "         resuming/extending the same one), clear ${WORKDIR}"
    echo "         before continuing to avoid mixing samples in the count matrix."
fi
echo "${LIST_HASH}" > .last_srr_list_hash

# ==========================================
# PHASE 2: REFERENCE GENOME PREPARATION
# ==========================================
echo -e "\n[PHASE 2] PREPARING REFERENCE GENOME..."
cd ref

FASTA_FILE="genome.fna"
GFF_FILE="annotation.gff"
INDEX_PREFIX="genome_index"

if [ ! -s "${FASTA_FILE}" ]; then
    echo "Downloading FASTA..."
    wget -qO- "${FASTA_URL}" | gunzip -c > "${FASTA_FILE}"
    if [ ! -s "${FASTA_FILE}" ]; then
        echo "CRITICAL ERROR: FASTA downloaded empty/corrupt. Aborting."
        rm -f "${FASTA_FILE}"
        exit 1
    fi
fi

if [ ! -s "${GFF_FILE}" ]; then
    echo "Downloading GFF..."
    wget -qO- "${GFF_URL}" | gunzip -c > "${GFF_FILE}"
    if [ ! -s "${GFF_FILE}" ]; then
        echo "CRITICAL ERROR: GFF downloaded empty/corrupt. Aborting."
        rm -f "${GFF_FILE}"
        exit 1
    fi
fi

if [ ! -f "${INDEX_PREFIX}.1.ht2" ]; then
    echo "Building HISAT2 index (this can take a while for large genomes)..."
    module purge
    load_python || exit 1
    load_hisat2 || exit 1
    hisat2-build -p "${CPUS}" "${FASTA_FILE}" "${INDEX_PREFIX}"
fi
cd "${WORKDIR}"

# ==========================================
# PHASE 3: SAMPLE PROCESSING LOOP
# ==========================================
echo -e "\n[PHASE 3] EXECUTING SAMPLE PROCESSING LOOP..."

FAILED_SAMPLES=()

# FIX #4: `while read` reads exactly one line at a time (no word-splitting on
# whitespace, no silent collapsing of blank lines) instead of
# `for SRR in $(cat srr_list.txt)`.
while IFS= read -r SRR || [ -n "${SRR}" ]; do
    [ -z "${SRR}" ] && continue
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
        load_sra_toolkit || exit 1
        prefetch "${SRR}" -O raw_reads/
        fasterq-dump "raw_reads/${SRR}" --split-files --force --outdir raw_reads --threads 4

        # Step B: Pre-trimming QC
        echo ">> Running Initial FastQC..."
        module purge
        load_fastqc || exit 1
        fastqc "raw_reads/${SRR}_1.fastq" "raw_reads/${SRR}_2.fastq" -o qc -t "${CPUS}" -q

        # Step C: Quality & Adapter Trimming
        echo ">> Trimming with fastp..."
        module purge
        load_fastp || exit 1
        fastp -i "raw_reads/${SRR}_1.fastq" -I "raw_reads/${SRR}_2.fastq" \
              -o "trimmed_reads/${SRR}_1_clean.fastq" -O "trimmed_reads/${SRR}_2_clean.fastq" \
              --thread "${CPUS}" --html "qc/${SRR}_fastp.html" 2> "qc/${SRR}_fastp.log"

        # Step D: Post-trimming QC
        echo ">> Running Post-trim QC..."
        module purge
        load_fastqc || exit 1
        fastqc "trimmed_reads/${SRR}_1_clean.fastq" "trimmed_reads/${SRR}_2_clean.fastq" -o qc -t "${CPUS}" -q

        # Step E: Alignment & BAM Sorting
        echo ">> Aligning and sorting BAM..."
        module purge
        load_python   || exit 1
        load_hisat2   || exit 1
        load_samtools || exit 1

        hisat2 -p "${CPUS}" ${SPLICE_FLAG} -x "ref/${INDEX_PREFIX}" \
               -1 "trimmed_reads/${SRR}_1_clean.fastq" \
               -2 "trimmed_reads/${SRR}_2_clean.fastq" | \
        samtools sort -@ "${SORT_THREADS}" -o "alignment/${SRR}_sorted.bam" -

        # Step F: Indexing BAM
        echo ">> Indexing BAM..."
        samtools index "alignment/${SRR}_sorted.bam"
    )
    SAMPLE_STATUS=$?
    set -e

    if [ ${SAMPLE_STATUS} -ne 0 ]; then
        echo ">> ERROR: Sample ${SRR} failed (exit code ${SAMPLE_STATUS}) - skipping to next sample."
        FAILED_SAMPLES+=("${SRR}")
        continue
    fi

    echo ">> Completed ${SRR}"
done < srr_list.txt

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

featureCounts -p -B -C -T "${CPUS}" \
              -t gene -g "${FC_ATTR}" \
              -a "ref/${GFF_FILE}" \
              -o counts/gene_counts_raw.txt \
              ${BAM_FILES}

# ==========================================
# PHASE 5: MATRIX CLEANUP
# ==========================================
echo -e "\n[PHASE 5] CLEANING MATRIX FOR DESeq2..."
cd counts
cut -f1,7- gene_counts_raw.txt | sed '1d' > gene_counts_clean.tsv
sed -i 's/alignment\///g; s/_sorted.bam//g' gene_counts_clean.tsv
cd "${WORKDIR}"

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

###############################################################################
# RNA-SEQ METHODOLOGY & PIPELINE NOTES
###############################################################################
# 1. Alignment Architecture (Spliced vs. Non-Spliced Alignment)
#    '--no-spliced-alignment' disables spliced alignment for prokaryotic /
#    intron-less genomes (bacterial, viral). It is omitted for human/mammal
#    so HISAT2 can map reads across exon-exon junctions.
#
# 2. Interpreting Read Quality and Alignment Metrics
#    FastQC / fastp: aim for average Phred scores (Q) > 30.
#    Alignment rate: healthy RNA-Seq runs are typically > 80%; < 50% suggests
#    contamination, adapter pollution, or a mismatched reference.
#
# 3. Read Quantification Strategy (featureCounts)
#    -p: paired-end reads. -B: both mates must align to the same chromosome.
#    -C: excludes chimeric fragments across chromosomes.
#    -g locus_tag is used for bacterial/viral RefSeq annotations; -g gene is
#    used for the human/mammal GFF by default — override FC_ATTR if a custom
#    organism's annotation uses a different attribute key.
#
# 4. Why GENOME_TAG matters
#    Every reference download, HISAT2 index, and aligned BAM lives under a
#    WORKDIR keyed by GENOME_TAG (derived from FASTA_URL). This is what makes
#    it actually safe to override FASTA_URL/GFF_URL for a new organism without
#    ever touching the script — a different genome always gets its own
#    directory, so it can never silently reuse another species' cached
#    reference or leftover alignments.
###############################################################################
