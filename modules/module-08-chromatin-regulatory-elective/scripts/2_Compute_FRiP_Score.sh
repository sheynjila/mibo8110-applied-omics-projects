#!/bin/bash
# ==============================================================================
# MODULE 8 · SCRIPT 2 of 5 — COMPUTE FRiP (FRACTION OF READS IN PEAKS)
# ==============================================================================
# ONE NEW IDEA on top of Script 1: Script 1 confirmed each sample's BAM and
# peaks file are structurally trustworthy; it did not ask whether the peaks
# actually explain where this sample's signal IS. FRiP (Fraction of Reads in
# Peaks) is the ENCODE-standard answer to this module's own guiding question
# - "how do I distinguish a real regulatory signal from assay background" -
# made into one number per sample: what fraction of this sample's mapped
# reads fall inside the peaks MACS2 called, versus scattered as background
# across the rest of the genome. A high peak count with a low FRiP is a real,
# common failure mode this script exists to catch: MACS2 can call thousands
# of "peaks" out of noise given a low enough significance threshold, and peak
# COUNT alone (Script 1's own metric) cannot tell that apart from a real,
# concentrated signal the way FRiP does.
#
# TEACHING NOTE — why FRiP thresholds are guidance, not a pass/fail law:
#   ENCODE's own published guidance treats FRiP > 1% as a minimum viability
#   bar and > 5% as a marker of a well-behaved transcription-factor ChIP-seq
#   experiment - but a broad histone mark (H3K27ac, H3K4me1) or an ATAC-seq
#   sample legitimately spreads signal across more of the genome than a
#   sharp TF binding profile does, and can have a lower FRiP while still
#   being a perfectly good experiment. This script flags low FRiP as
#   something to look at, not something to discard on sight.
#
# Input:  WORKDIR - alignment/<SRR>_dedup.bam and
#         macs2_output/<SRR>/<SRR>_peaks.narrowPeak (Script 1's inputs).
# Output: prints a per-sample FRiP report; writes frip_summary.tsv
# ==============================================================================

set -o pipefail

WORKDIR="${1:-/scratch/$(whoami)/master_chipseq_pipeline_ChIP}"
OUTDIR="${2:-.}"
LOW_FRIP_PCT=1.0   # ENCODE's own minimum-viability guidance (see header)

mkdir -p "${OUTDIR}"
OUT_TSV="${OUTDIR}/frip_summary.tsv"

echo "=================================================================="
echo " MODULE 8 / SCRIPT 2 — FRiP (Fraction of Reads in Peaks)"
echo "=================================================================="
echo "Working directory: ${WORKDIR}"
echo ""

shopt -s nullglob
PEAK_FILES=("${WORKDIR}"/macs2_output/*/*_peaks.narrowPeak)
shopt -u nullglob

if [ ${#PEAK_FILES[@]} -eq 0 ]; then
    echo "[FAIL] no *_peaks.narrowPeak files found under ${WORKDIR}/macs2_output — run Script 1 first"
    exit 1
fi

echo -e "sample\ttotal_mapped_reads\treads_in_peaks\tfrip_pct\tflag" > "${OUT_TSV}"

N_LOW=0
for PEAKS in "${PEAK_FILES[@]}"; do
    SRR=$(basename "$(dirname "${PEAKS}")")
    BAM="${WORKDIR}/alignment/${SRR}_dedup.bam"

    if [ ! -s "${BAM}" ]; then
        echo "  ${SRR}: [FAIL] BAM not found — skipping"
        continue
    fi

    TOTAL=$(samtools view -c -F 4 "${BAM}")
    IN_PEAKS=$(samtools view -c -F 4 -L "${PEAKS}" "${BAM}")

    if [ "${TOTAL}" -eq 0 ]; then
        echo "  ${SRR}: [FAIL] zero mapped reads in ${BAM} — cannot compute FRiP"
        continue
    fi

    FRIP=$(awk -v n="${IN_PEAKS}" -v d="${TOTAL}" 'BEGIN { printf "%.2f", (100.0*n/d) }')
    FLAG="ok"
    if awk -v f="${FRIP}" -v t="${LOW_FRIP_PCT}" 'BEGIN { exit !(f < t) }'; then
        FLAG="LOW_FRIP"
        N_LOW=$((N_LOW + 1))
    fi

    echo "  ${SRR}: total=${TOTAL}  in_peaks=${IN_PEAKS}  FRiP=${FRIP}%   ${FLAG}"
    echo -e "${SRR}\t${TOTAL}\t${IN_PEAKS}\t${FRIP}\t${FLAG}" >> "${OUT_TSV}"
done

echo ""
if [ ${N_LOW} -gt 0 ]; then
    echo "[WARN] ${N_LOW} sample(s) have FRiP below ${LOW_FRIP_PCT}% (ENCODE's minimum-viability guidance)."
    echo "       Check whether this reflects the assay type (broad marks/ATAC legitimately run lower -"
    echo "       see this script's header) before concluding the experiment itself failed."
else
    echo "[PASS] every sample meets ENCODE's ${LOW_FRIP_PCT}% minimum-viability FRiP guidance."
fi
echo ""
echo "Wrote ${OUT_TSV}"
echo "=================================================================="
