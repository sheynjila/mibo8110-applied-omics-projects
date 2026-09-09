#!/bin/bash
###############################################################################
# smoke_test_prjna1518998.sh
# ---------------------------------------------------------------------------
# ONE-COMMAND SMOKE TEST for master_rnaseq_pipeline_consolidated.sh
#
# Wraps the "Quick Start" example from that script's header into a single
# runnable helper: submits a 2-sample test run against REAL PRJNA1518998
# (Mouse retina) data using the verified GRCm39 RefSeq FASTA/GFF, waits for
# the Slurm job to leave the queue, then automatically checks the three
# things that matter before you trust the pipeline on a full study:
#
#   1. Did it actually use the Mouse genome (not a stale cached one from a
#      previous species run under the same REFERENCE_TYPE)?
#   2. Did both samples align well (>=80% overall alignment rate)?
#   3. Does the resulting count matrix have real, non-zero numbers?
#
# WHY A SEPARATE WRAPPER (teaching note):
#   Keeping this logic OUT of master_rnaseq_pipeline_consolidated.sh matters.
#   That script is a general-purpose pipeline meant to run unattended inside
#   an sbatch job; it should not know about any specific test dataset. This
#   wrapper is disposable, study-specific glue that lives on the login node —
#   it submits the job, then polls and inspects results from OUTSIDE the job.
#   That separation is exactly what lets you swap in a different smoke test
#   (different organism, different accessions) later without touching the
#   pipeline itself.
#
# USAGE:
#   ./smoke_test_prjna1518998.sh
#   (run from the login node, in the same directory as
#    master_rnaseq_pipeline_consolidated.sh — or set PIPELINE_SCRIPT below)
#
# WHAT THIS DOES NOT DO:
#   It is a fast PASS/FAIL gate, not a full report. If something fails, go
#   read the actual .out/.err logs it points you to.
###############################################################################

set -e
set -o pipefail

# ---- Configuration (only these three lines change for a different smoke test) ----
PIPELINE_SCRIPT="${PIPELINE_SCRIPT:-./master_rnaseq_pipeline_consolidated.sh}"
REFERENCE_TYPE="mammal"
FASTA_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/635/GCF_000001635.27_GRCm39/GCF_000001635.27_GRCm39_genomic.fna.gz"
GFF_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/635/GCF_000001635.27_GRCm39/GCF_000001635.27_GRCm39_genomic.gff.gz"
SRR_LIST=("SRR40359383" "SRR40359384")   # real paired-end runs from PRJNA1518998

if [ ! -f "${PIPELINE_SCRIPT}" ]; then
    echo "ERROR: Can't find ${PIPELINE_SCRIPT}."
    echo "       Run this wrapper from the directory containing"
    echo "       master_rnaseq_pipeline_consolidated.sh, or set:"
    echo "       PIPELINE_SCRIPT=/full/path/to/it ./smoke_test_prjna1518998.sh"
    exit 1
fi

# Replicate the EXACT GENOME_TAG derivation used inside the main pipeline so
# this wrapper always inspects the SAME WORKDIR the sbatch job will use, even
# if FASTA_URL above is later changed to a different organism/assembly.
URL_HASH="$(printf '%s' "${FASTA_URL}" | md5sum | cut -c1-8)"
BASE_TAG="$(basename "${FASTA_URL}" 2>/dev/null | sed -E 's/\.(fna|fa|fasta)\.gz$//; s/[^A-Za-z0-9._-]/_/g')"
GENOME_TAG="${BASE_TAG:-genome}_${URL_HASH}"
WORKDIR="/scratch/$(whoami)/master_rnaseq_pipeline_${REFERENCE_TYPE}_${GENOME_TAG}"

echo "=========================================================="
echo " SMOKE TEST: PRJNA1518998 (Mouse retina), ${#SRR_LIST[@]} sample(s)"
echo " WORKDIR: ${WORKDIR}"
echo "=========================================================="

# ---- Step 1: pre-create WORKDIR and drop the sample list in it ----
# The pipeline cd's into WORKDIR before it looks for srr_list.txt, so it must
# be written there in advance, not in the submission directory.
mkdir -p "${WORKDIR}"
printf '%s\n' "${SRR_LIST[@]}" > "${WORKDIR}/srr_list.txt"
echo "Wrote srr_list.txt (${#SRR_LIST[@]} accession[s]) to ${WORKDIR}"

# ---- Step 2: submit and capture the job ID ----
SUBMIT_OUT=$(sbatch --export="REFERENCE_TYPE=${REFERENCE_TYPE},FASTA_URL=${FASTA_URL},GFF_URL=${GFF_URL}" "${PIPELINE_SCRIPT}")
echo "${SUBMIT_OUT}"
JOB_ID=$(echo "${SUBMIT_OUT}" | grep -oE '[0-9]+$')

if [ -z "${JOB_ID}" ]; then
    echo "ERROR: Could not parse a job ID from sbatch's output above."
    echo "       Check 'squeue -u $(whoami)' manually - the job may still"
    echo "       have been submitted even though this wrapper can't track it."
    exit 1
fi

OUT_LOG="master_rnaseq_${JOB_ID}.out"
ERR_LOG="master_rnaseq_${JOB_ID}.err"

# ---- Step 3: wait for the job to leave the queue ----
echo "Submitted job ${JOB_ID}. Polling squeue every 30s..."
echo "(Ctrl+C any time - the Slurm job keeps running; rerun 'squeue -u $(whoami)'"
echo " or tail ${OUT_LOG}/${ERR_LOG} yourself if you stop watching.)"
while squeue -j "${JOB_ID}" -h 2>/dev/null | grep -q "${JOB_ID}"; do
    sleep 30
done
echo -e "\nJob ${JOB_ID} has left the queue. Running PASS/FAIL checks...\n"

# ---- Step 4: automated PASS/FAIL checks ----
PASS=true

echo "---- CHECK 1: Reference manifest (correct genome was used) ----"
if [ -f "${WORKDIR}/logs/reference_manifest.txt" ]; then
    tail -8 "${WORKDIR}/logs/reference_manifest.txt"
    if tail -8 "${WORKDIR}/logs/reference_manifest.txt" | grep -q "GRCm39"; then
        echo "PASS: manifest confirms the GRCm39 Mouse genome was used."
    else
        echo "FAIL: manifest does not mention GRCm39 - possible stale WORKDIR."
        PASS=false
    fi
else
    echo "FAIL: no reference_manifest.txt found - job may not have reached Phase 1."
    PASS=false
fi

echo -e "\n---- CHECK 2: Alignment rate (>=80% expected) ----"
if [ -f "${ERR_LOG}" ] && grep -q "overall alignment rate" "${ERR_LOG}"; then
    grep -B4 "overall alignment rate" "${ERR_LOG}"
    LOW_RATE=$(grep "overall alignment rate" "${ERR_LOG}" | grep -oE '^[0-9]+\.[0-9]+' | awk '$1<80{print;exit}')
    if [ -n "${LOW_RATE}" ]; then
        echo "FAIL: at least one sample aligned below 80% (${LOW_RATE}%) - possible reference/species mismatch."
        PASS=false
    else
        echo "PASS: all reported alignment rates are >=80%."
    fi
else
    echo "FAIL: no alignment-rate line found in ${ERR_LOG} - inspect it manually."
    PASS=false
fi

echo -e "\n---- CHECK 3: Count matrix has real, non-zero data ----"
COUNTS_FILE="${WORKDIR}/counts/gene_counts_clean.tsv"
if [ -s "${COUNTS_FILE}" ]; then
    NONZERO_ROWS=$(awk -F'\t' 'NR>1 { s=0; for(i=2;i<=NF;i++) s+=$i; if (s>0) c++ } END{print c+0}' "${COUNTS_FILE}")
    echo "Genes with non-zero total counts: ${NONZERO_ROWS}"
    if [ "${NONZERO_ROWS}" -gt 0 ]; then
        echo "PASS: count matrix has real expression data."
    else
        echo "FAIL: every gene has zero counts across both samples - check FC_ATTR/annotation match."
        PASS=false
    fi
else
    echo "FAIL: ${COUNTS_FILE} not found or empty."
    PASS=false
fi

echo -e "\n=========================================================="
if [ "${PASS}" = true ]; then
    echo " SMOKE TEST PASSED - safe to scale up to the full sample list."
    echo " Next: replace ${WORKDIR}/srr_list.txt with the full PRJNA1518998"
    echo " accession list, then resubmit with the SAME sbatch command:"
    echo ""
    echo "   sbatch --export=REFERENCE_TYPE=${REFERENCE_TYPE},FASTA_URL=\"${FASTA_URL}\",GFF_URL=\"${GFF_URL}\" ${PIPELINE_SCRIPT}"
else
    echo " SMOKE TEST FAILED - see the FAIL line(s) above before scaling up."
    echo " Logs to inspect: ${OUT_LOG}  ${ERR_LOG}"
fi
echo "=========================================================="
