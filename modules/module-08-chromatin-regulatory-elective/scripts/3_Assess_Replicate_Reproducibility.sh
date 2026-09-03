#!/bin/bash
# ==============================================================================
# MODULE 8 · SCRIPT 3 of 5 — ASSESS REPLICATE REPRODUCIBILITY
# ==============================================================================
# ONE NEW IDEA on top of Scripts 1-2: those scripts characterized ONE
# sample's peaks in isolation. Neither can answer this module's own named
# learning outcome - "interpret peak scores and reproducibility across
# replicates" - because reproducibility is not a property of a single peak
# set at all; it only exists as a comparison between two or more. A single
# replicate calling 20,000 confident-looking peaks tells you nothing about
# whether a SECOND replicate of the same experiment would call the same
# 20,000 regions or a mostly-different set - and a peak set that cannot be
# reproduced in a second replicate is exactly the "looks real, isn't" case
# this module's guiding question is about.
#
# TEACHING NOTE — why bedtools jaccard, not a raw peak-count comparison:
#   Two replicates could report an identical PEAK COUNT while sharing almost
#   no actual genomic positions - count alone cannot distinguish that from
#   genuine reproducibility. bedtools jaccard instead measures the actual
#   overlap: (bases covered by BOTH peak sets) / (bases covered by EITHER),
#   giving a single 0-1 score per pair that is 0 only when the two sets
#   share no genomic positions at all and 1 only when they cover exactly the
#   same bases - a direct, position-aware reproducibility measure.
#
# SCOPE NOTE: this script computes pairwise Jaccard similarity for EVERY
# pair of peak sets found under WORKDIR/macs2_output/, not just pairs that
# are actually biological replicates of each other - srr_list.txt carries no
# replicate-group metadata for this script to use. For this module's own
# worked example (three real biological replicates of the same condition -
# see this module's README), all pairs ARE replicate comparisons; for a
# batch mixing multiple conditions, the analyst must judge which pairwise
# scores below are meaningful replicate comparisons and which are
# expected-to-differ cross-condition comparisons - this script does not
# and cannot make that distinction on its own.
#
# Input:  WORKDIR - macs2_output/<SRR>/<SRR>_peaks.narrowPeak (Script 1's inputs).
# Output: prints a pairwise Jaccard reproducibility report; writes
#         replicate_reproducibility.tsv
# ==============================================================================

set -o pipefail

WORKDIR="${1:-/scratch/$(whoami)/master_chipseq_pipeline_ChIP}"
OUTDIR="${2:-.}"
LOW_JACCARD=0.20   # below this, treat the pair as poorly reproduced (see TEACHING NOTE)

mkdir -p "${OUTDIR}"
OUT_TSV="${OUTDIR}/replicate_reproducibility.tsv"

echo "=================================================================="
echo " MODULE 8 / SCRIPT 3 — Replicate reproducibility (pairwise Jaccard)"
echo "=================================================================="
echo "Working directory: ${WORKDIR}"
echo ""

shopt -s nullglob
PEAK_FILES=("${WORKDIR}"/macs2_output/*/*_peaks.narrowPeak)
shopt -u nullglob

if [ ${#PEAK_FILES[@]} -lt 2 ]; then
    echo "[FAIL] fewer than 2 peak sets found under ${WORKDIR}/macs2_output — reproducibility needs at least a pair."
    exit 1
fi

echo -e "sample_a\tsample_b\tjaccard\tn_intersections\tflag" > "${OUT_TSV}"

N_LOW=0
N_PAIRS=0
for ((i = 0; i < ${#PEAK_FILES[@]}; i++)); do
    for ((j = i + 1; j < ${#PEAK_FILES[@]}; j++)); do
        A="${PEAK_FILES[$i]}"
        B="${PEAK_FILES[$j]}"
        SRR_A=$(basename "$(dirname "${A}")")
        SRR_B=$(basename "$(dirname "${B}")")

        SORTED_A=$(mktemp)
        SORTED_B=$(mktemp)
        sort -k1,1 -k2,2n "${A}" > "${SORTED_A}"
        sort -k1,1 -k2,2n "${B}" > "${SORTED_B}"

        JACCARD_LINE=$(bedtools jaccard -a "${SORTED_A}" -b "${SORTED_B}" | tail -n1)
        JACCARD=$(echo "${JACCARD_LINE}" | cut -f3)
        N_INTERSECT=$(echo "${JACCARD_LINE}" | cut -f4)
        rm -f "${SORTED_A}" "${SORTED_B}"

        FLAG="ok"
        if awk -v j="${JACCARD}" -v t="${LOW_JACCARD}" 'BEGIN { exit !(j < t) }'; then
            FLAG="LOW_REPRODUCIBILITY"
            N_LOW=$((N_LOW + 1))
        fi
        N_PAIRS=$((N_PAIRS + 1))

        echo "  ${SRR_A} vs ${SRR_B}: jaccard=${JACCARD}  n_intersections=${N_INTERSECT}   ${FLAG}"
        echo -e "${SRR_A}\t${SRR_B}\t${JACCARD}\t${N_INTERSECT}\t${FLAG}" >> "${OUT_TSV}"
    done
done

echo ""
if [ ${N_LOW} -gt 0 ]; then
    echo "[WARN] ${N_LOW}/${N_PAIRS} pair(s) fall below Jaccard ${LOW_JACCARD} — if these ARE supposed to be"
    echo "       replicates of the same condition (confirm from srr_list.txt / sample metadata, this script"
    echo "       cannot tell), that is a real reproducibility concern worth investigating before trusting"
    echo "       either peak set's calls in isolation."
else
    echo "[PASS] all ${N_PAIRS} pair(s) meet the ${LOW_JACCARD} Jaccard threshold."
fi
echo ""
echo "Wrote ${OUT_TSV}"
echo "=================================================================="
