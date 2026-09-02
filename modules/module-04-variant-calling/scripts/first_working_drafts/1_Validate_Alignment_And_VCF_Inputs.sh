#!/bin/bash
# ==============================================================================
# MODULE 4 · SCRIPT 1 of 7 — VALIDATE ALIGNMENT AND VCF INPUTS
# ==============================================================================
# ONE NEW IDEA on top of master_snp_pipeline.sh: that script produces BAMs
# and VCFs but never checks that they actually belong to each other. This
# script validates the pipeline's own OUTPUT before any downstream script
# trusts it — the same "check inputs before trusting outputs" habit Module 2
# Script 1 applied to count matrices, one level down the stack at the
# alignment/calling layer.
#
# TEACHING NOTE — why this has to be its own step, not folded into Script 3:
#   samtools and bcftools will both happily operate on a BAM whose @RG SM:
#   tag does not match the sample name in its own filename, or on a VCF
#   whose lone sample column has been silently renamed by a copy/rename
#   mistake. Neither tool errors — they just produce output labeled with
#   whatever name is actually in the file header, which may not be the
#   sample you think you are looking at. In a cohort merge (Phase 5 of
#   master_snp_pipeline.sh) a name collision or silent mismatch is exactly
#   how two different isolates end up merged under one sample column, or
#   how an outbreak-cluster call in variant_analysis.R gets attributed to
#   the wrong isolate. Catching this here is cheap; catching it after a
#   cluster call has already been reported is not.
#
# Input:  WORKDIR — the master_snp_pipeline.sh working directory, containing
#         ref/genome.fna[.fai], alignment/<SRR>_dedup.bam[.bai], and
#         variants/<SRR>_filtered.vcf for one or more samples.
# Output: prints a PASS/FAIL/WARN validation report to the console per
#         sample; stops with a clear error on the first hard failure.
# ==============================================================================

set -o pipefail

WORKDIR="${1:-/scratch/$(whoami)/master_snp_pipeline}"

echo "=================================================================="
echo " MODULE 4 / SCRIPT 1 — Alignment + VCF input validation"
echo "=================================================================="
echo "Working directory: ${WORKDIR}"
echo ""

fail() { echo "[FAIL] $1"; exit 1; }
warn_flag() { echo "[WARN] $1"; }
pass() { echo "[PASS] $1"; }

for tool in samtools bcftools; do
    command -v "${tool}" >/dev/null 2>&1 || fail "required tool '${tool}' not found on PATH — load its module before running this script"
done
pass "samtools and bcftools are both available"

REF="${WORKDIR}/ref/genome.fna"
[ -s "${REF}" ]         || fail "reference FASTA not found or empty: ${REF}"
[ -s "${REF}.fai" ]     || fail "reference FASTA index missing: ${REF}.fai (run samtools faidx)"
[ -s "${REF}.bwt" ]     || fail "BWA index missing: ${REF}.bwt (run bwa index)"
pass "reference FASTA and its samtools/BWA indexes are present and non-empty"

shopt -s nullglob
VCFS=("${WORKDIR}"/variants/*_filtered.vcf)
shopt -u nullglob

[ ${#VCFS[@]} -gt 0 ] || fail "no *_filtered.vcf files found in ${WORKDIR}/variants — run master_snp_pipeline.sh first"

HARD_FAILS=0
for VCF in "${VCFS[@]}"; do
    SRR=$(basename "${VCF}" | sed 's/_filtered\.vcf$//')
    BAM="${WORKDIR}/alignment/${SRR}_dedup.bam"
    echo "------------------------------------------------------------"
    echo " Sample: ${SRR}"
    echo "------------------------------------------------------------"

    # --- BAM must exist, be indexed, and pass samtools' own structural check ---
    if [ ! -s "${BAM}" ]; then
        echo "[FAIL] BAM not found or empty: ${BAM}"
        HARD_FAILS=$((HARD_FAILS + 1))
        continue
    fi
    if [ ! -s "${BAM}.bai" ] && [ ! -s "${BAM}.csi" ]; then
        echo "[FAIL] BAM index missing for ${SRR} (expected ${BAM}.bai)"
        HARD_FAILS=$((HARD_FAILS + 1))
        continue
    fi
    if ! samtools quickcheck "${BAM}" 2>/dev/null; then
        echo "[FAIL] samtools quickcheck failed on ${BAM} — truncated or corrupt BAM"
        HARD_FAILS=$((HARD_FAILS + 1))
        continue
    fi
    echo "[PASS] BAM exists, is indexed, and passes samtools quickcheck"

    # --- BAM read-group SM: tag must match the sample name, BY NAME not by
    #     filename-only assumption ---
    RG_SM=$(samtools view -H "${BAM}" | grep '^@RG' | head -n1 | \
            grep -oP '(?<=\bSM:)[^\t]+' || true)
    if [ -z "${RG_SM}" ]; then
        echo "[WARN] no @RG SM: tag found in BAM header for ${SRR} — cannot cross-check sample identity"
    elif [ "${RG_SM}" != "${SRR}" ]; then
        echo "[FAIL] @RG SM:${RG_SM} in BAM header does not match filename-derived sample ${SRR} — likely a rename/copy mistake"
        HARD_FAILS=$((HARD_FAILS + 1))
        continue
    else
        echo "[PASS] BAM @RG SM: tag matches sample name (${SRR})"
    fi

    # --- VCF must exist, be non-empty, and its sample column must match too ---
    if [ ! -s "${VCF}" ]; then
        echo "[FAIL] VCF not found or empty: ${VCF}"
        HARD_FAILS=$((HARD_FAILS + 1))
        continue
    fi
    VCF_SAMPLES=$(bcftools query -l "${VCF}" 2>/dev/null || true)
    N_VCF_SAMPLES=$(echo "${VCF_SAMPLES}" | grep -c . || true)
    if [ "${N_VCF_SAMPLES}" -ne 1 ]; then
        echo "[FAIL] expected exactly 1 sample column in ${VCF}, found ${N_VCF_SAMPLES}"
        HARD_FAILS=$((HARD_FAILS + 1))
        continue
    fi
    if [ "${VCF_SAMPLES}" != "${SRR}" ]; then
        echo "[FAIL] VCF sample column '${VCF_SAMPLES}' does not match filename-derived sample ${SRR}"
        HARD_FAILS=$((HARD_FAILS + 1))
        continue
    fi
    echo "[PASS] VCF exists, is non-empty, and its sample column matches ${SRR} by name"

    N_RECORDS=$(bcftools view -H "${VCF}" 2>/dev/null | wc -l)
    echo "       ${N_RECORDS} variant records (PASS + LowQual-tagged) in ${VCF}"
done

echo ""
echo "=================================================================="
if [ ${HARD_FAILS} -gt 0 ]; then
    fail "${HARD_FAILS} sample(s) failed validation — fix before proceeding to Script 2"
fi
echo " VALIDATION COMPLETE for ${#VCFS[@]} sample(s) — safe to proceed to Script 2"
echo "=================================================================="
