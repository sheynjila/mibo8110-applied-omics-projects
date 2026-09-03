#!/bin/bash
# ==============================================================================
# MODULE 8 · SCRIPT 1 of 5 — VALIDATE PEAK-CALLING INPUTS
# ==============================================================================
# ONE NEW IDEA on top of master_chipseq_pipeline.sh: that script produces a
# deduplicated BAM and a MACS2 peaks file per sample but never checks that
# they actually belong to each other, or that a ChIP-mode run's peaks were
# actually called with an input control - the same "check the pipeline's own
# output before trusting it" habit Module 4 Script 1 applied to BAM/VCF
# pairs and Module 6 Script 1 applied to Kraken2/Bracken pairs.
#
# TEACHING NOTE — why the missing-input-control check matters more here than
# a generic file-existence check would suggest:
#   A MACS2 run without `-c` does not fail, warn loudly in its own log by
#   default, or produce an obviously different-looking output format from a
#   properly controlled run - narrowPeak files from both look identical in
#   shape. The only way to know a given peak set was called without a
#   background model is to check the pipeline's own MACS2 log for the `-c`
#   flag actually being present, which is exactly what this script does
#   instead of assuming every peak set in a batch was called the same way.
#
# Input:  WORKDIR - the master_chipseq_pipeline.sh working directory,
#         containing alignment/<SRR>_dedup.bam[.bai], macs2_output/<SRR>/
#         (with <SRR>_peaks.narrowPeak and MACS2's own log), and
#         annotation/<SRR>_peaks_annotated.bed.
# Output: prints a PASS/FAIL/WARN validation report to the console per
#         sample; stops with a clear error on the first hard failure.
# ==============================================================================

set -o pipefail

WORKDIR="${1:-/scratch/$(whoami)/master_chipseq_pipeline_ChIP}"

echo "=================================================================="
echo " MODULE 8 / SCRIPT 1 — Peak-calling input validation"
echo "=================================================================="
echo "Working directory: ${WORKDIR}"
echo ""

fail() { echo "[FAIL] $1"; exit 1; }
pass() { echo "[PASS] $1"; }

shopt -s nullglob
PEAK_FILES=("${WORKDIR}"/macs2_output/*/*_peaks.narrowPeak)
shopt -u nullglob

[ ${#PEAK_FILES[@]} -gt 0 ] || fail "no *_peaks.narrowPeak files found under ${WORKDIR}/macs2_output — run master_chipseq_pipeline.sh first"

HARD_FAILS=0
N_NO_CONTROL=0
for PEAKS in "${PEAK_FILES[@]}"; do
    SRR=$(basename "$(dirname "${PEAKS}")")
    BAM="${WORKDIR}/alignment/${SRR}_dedup.bam"
    MACS2_LOG="${WORKDIR}/qc/${SRR}_macs2.log"
    ANNOT="${WORKDIR}/annotation/${SRR}_peaks_annotated.bed"
    echo "------------------------------------------------------------"
    echo " Sample: ${SRR}"
    echo "------------------------------------------------------------"

    if [ ! -s "${BAM}" ]; then
        echo "[FAIL] deduplicated BAM not found or empty: ${BAM}"
        HARD_FAILS=$((HARD_FAILS + 1))
        continue
    fi
    if [ ! -s "${BAM}.bai" ]; then
        echo "[FAIL] BAM index missing for ${SRR} (expected ${BAM}.bai)"
        HARD_FAILS=$((HARD_FAILS + 1))
        continue
    fi
    if ! samtools quickcheck "${BAM}" 2>/dev/null; then
        echo "[FAIL] samtools quickcheck failed on ${BAM} — truncated or corrupt BAM"
        HARD_FAILS=$((HARD_FAILS + 1))
        continue
    fi
    echo "[PASS] deduplicated BAM exists, is indexed, and passes samtools quickcheck"

    N_PEAKS=$(wc -l < "${PEAKS}")
    if [ "${N_PEAKS}" -eq 0 ]; then
        echo "[FAIL] ${PEAKS} exists but has zero peaks called"
        HARD_FAILS=$((HARD_FAILS + 1))
        continue
    fi
    echo "[PASS] ${N_PEAKS} peak(s) called for ${SRR}"

    if [ ! -s "${ANNOT}" ]; then
        echo "[WARN] no gene annotation found for ${SRR} at ${ANNOT} — Script 4 will have nothing to summarize for this sample"
    else
        echo "[PASS] peak-to-gene annotation present (${ANNOT})"
    fi

    # Detect whether this sample's MACS2 run was actually given an input
    # control, from MACS2's own log rather than assuming from CONTROL_SRR
    # alone (a batch could mix controlled and uncontrolled runs across
    # different sbatch submissions).
    if [ -s "${MACS2_LOG}" ]; then
        if grep -q -- "-c " "${MACS2_LOG}" || grep -qi "control file" "${MACS2_LOG}"; then
            echo "[PASS] MACS2 log indicates an input control was used for ${SRR}"
        else
            N_NO_CONTROL=$((N_NO_CONTROL + 1))
            echo "[WARN] MACS2 log shows NO input control was used for ${SRR} — these peaks were called"
            echo "       from the treatment track alone. This is not automatically wrong (see this"
            echo "       module's README for a real published dataset with no deposited input control),"
            echo "       but every downstream interpretation of these peaks must say so explicitly."
        fi
    else
        echo "[WARN] MACS2 log not found for ${SRR} (${MACS2_LOG}) — cannot confirm input-control status"
    fi
done

echo ""
echo "=================================================================="
if [ ${HARD_FAILS} -gt 0 ]; then
    fail "${HARD_FAILS} sample(s) failed validation — fix before proceeding to Script 2"
fi
if [ ${N_NO_CONTROL} -gt 0 ]; then
    echo " VALIDATION COMPLETE for ${#PEAK_FILES[@]} sample(s) — ${N_NO_CONTROL} of them called WITHOUT an input control."
    echo " Carry this forward explicitly into Script 5's report; do not let it get lost."
else
    echo " VALIDATION COMPLETE for ${#PEAK_FILES[@]} sample(s), all with a confirmed input control — safe to proceed to Script 2"
fi
echo "=================================================================="
