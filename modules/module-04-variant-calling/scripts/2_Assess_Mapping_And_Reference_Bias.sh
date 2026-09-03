#!/bin/bash
# ==============================================================================
# MODULE 4 · SCRIPT 2 of 7 — ASSESS MAPPING AND REFERENCE BIAS
# ==============================================================================
# ONE NEW IDEA on top of Script 1: Script 1 confirmed the files are
# structurally sound and correctly labeled. This script asks a different
# question — even a perfectly valid, correctly labeled BAM can still be a
# POOR alignment, and every variant call downstream inherits that weakness
# silently. This is a first look at alignment quality before trusting any
# variant filter built on top of it (Script 3).
#
# TEACHING NOTE — what "reference bias" means operationally here:
#   Reference bias is the tendency for reads carrying the reference allele
#   (or from a strain close to the reference) to map more easily than reads
#   from a divergent strain, which under-recruits into low-identity or
#   absent regions of the reference. Two directly measurable symptoms:
#     1. Overall mapping rate (samtools flagstat) — a low rate can mean the
#        sample is genuinely divergent from the reference, or contaminated,
#        or simply mislabeled (see Script 1).
#     2. Breadth of coverage (samtools depth) — the fraction of the
#        reference genome with ZERO read depth. A region with no coverage
#        cannot produce a variant call there at all — not "no variant",
#        but "no information". A high zero-coverage fraction is a silent
#        blind spot in the call set, not a clean negative result.
#   Neither symptom alone proves reference bias (a genuinely small/complete
#   genome, or a sample with lower sequencing depth, can look similar) —
#   but a sample flagged on BOTH counts is a real reference-bias risk that
#   the eventual limitations statement (Script 6) needs to name explicitly.
#
# Input:  WORKDIR — the master_snp_pipeline.sh working directory.
# Output: prints a per-sample mapping/coverage table; writes
#         mapping_bias_summary.tsv into WORKDIR; flags samples for Script 6.
# ==============================================================================

set -o pipefail

WORKDIR="${1:-/scratch/$(whoami)/master_snp_pipeline}"
REF="${WORKDIR}/ref/genome.fna"
OUT_TSV="${WORKDIR}/mapping_bias_summary.tsv"

echo "=================================================================="
echo " MODULE 4 / SCRIPT 2 — Mapping and reference-bias assessment"
echo "=================================================================="
echo "Working directory: ${WORKDIR}"
echo ""

command -v samtools >/dev/null 2>&1 || { echo "[FAIL] samtools not found on PATH"; exit 1; }

GENOME_LEN=$(awk '/^>/{next}{n+=length($0)}END{print n}' "${REF}" 2>/dev/null)
if [ -z "${GENOME_LEN}" ] || [ "${GENOME_LEN}" -eq 0 ]; then
    echo "[FAIL] could not determine reference genome length from ${REF}"
    exit 1
fi
echo "Reference genome length: ${GENOME_LEN} bp"
echo ""

# Mapping rate below this is flagged as a mapping-bias risk. Zero-coverage
# fraction above this is flagged as a coverage-gap risk. Both thresholds are
# documented starting points, not universal cutoffs — Script 6 restates them
# next to the samples they flagged so a reader does not have to guess why.
MIN_MAPPED_PCT=95.0
MAX_ZERO_COV_PCT=5.0

printf "sample\tpct_mapped\tpct_properly_paired\tmean_depth\tpct_zero_coverage\tflag\n" > "${OUT_TSV}"

shopt -s nullglob
BAMS=("${WORKDIR}"/alignment/*_dedup.bam)
shopt -u nullglob
[ ${#BAMS[@]} -gt 0 ] || { echo "[FAIL] no *_dedup.bam files found in ${WORKDIR}/alignment"; exit 1; }

printf "%-14s %10s %10s %10s %10s   %s\n" "sample" "pct_mapped" "pct_pair" "mean_dp" "pct_0cov" "flag"

for BAM in "${BAMS[@]}"; do
    SRR=$(basename "${BAM}" _dedup.bam)

    FLAGSTAT=$(samtools flagstat "${BAM}")
    PCT_MAPPED=$(echo "${FLAGSTAT}" | grep 'mapped (' | head -n1 | grep -oP '(?<=\()[0-9.]+(?=%)')
    PCT_PAIRED=$(echo "${FLAGSTAT}" | grep 'properly paired' | grep -oP '(?<=\()[0-9.]+(?=%)')
    PCT_MAPPED=${PCT_MAPPED:-0}
    PCT_PAIRED=${PCT_PAIRED:-0}

    # samtools depth -a reports every reference position, including 0-depth
    # ones (without -a, zero-depth positions are silently omitted, which
    # would make the genome look fully covered no matter what).
    read -r MEAN_DEPTH ZERO_COV_PCT < <(
        samtools depth -a "${BAM}" | \
        awk -v glen="${GENOME_LEN}" '
            { sum += $3; if ($3 == 0) zero++ ; n++ }
            END {
                mean = (n > 0) ? sum / n : 0
                zpct = (glen > 0) ? 100.0 * zero / glen : 0
                printf "%.2f %.2f", mean, zpct
            }'
    )

    FLAG="ok"
    if awk -v p="${PCT_MAPPED}" -v m="${MIN_MAPPED_PCT}" 'BEGIN{exit !(p<m)}'; then
        FLAG="LOW_MAPPING_RATE"
    fi
    if awk -v p="${ZERO_COV_PCT}" -v m="${MAX_ZERO_COV_PCT}" 'BEGIN{exit !(p>m)}'; then
        FLAG="${FLAG}+COVERAGE_GAPS"
        FLAG="${FLAG#ok+}"
    fi

    printf "%-14s %9s%% %9s%% %10s %9s%%   %s\n" "${SRR}" "${PCT_MAPPED}" "${PCT_PAIRED}" "${MEAN_DEPTH}" "${ZERO_COV_PCT}" "${FLAG}"
    printf "%s\t%s\t%s\t%s\t%s\t%s\n" "${SRR}" "${PCT_MAPPED}" "${PCT_PAIRED}" "${MEAN_DEPTH}" "${ZERO_COV_PCT}" "${FLAG}" >> "${OUT_TSV}"
done

echo ""
echo "Wrote ${OUT_TSV}"
echo "Thresholds used: flag LOW_MAPPING_RATE if % mapped < ${MIN_MAPPED_PCT}%,"
echo "                 flag COVERAGE_GAPS   if % zero-coverage bases > ${MAX_ZERO_COV_PCT}%."
echo "A sample flagged on BOTH is the strongest reference-bias signal — carried"
echo "forward into Script 6's limitations statement."
echo "=================================================================="
