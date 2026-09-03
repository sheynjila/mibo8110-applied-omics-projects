#!/bin/bash
# ==============================================================================
# MODULE 4 · SCRIPT 3 of 7 — DESIGN A DEFENSIBLE VARIANT FILTER
# ==============================================================================
# ONE NEW IDEA on top of master_snp_pipeline.sh: that script's Phase 3 Step F
# applies one fixed, unexamined filter to every sample, every time:
#     QUAL >= 20, DP >= 10, MQ >= 30
# with no strand-bias check at all, despite this module's own learning
# outcomes naming depth, quality, AND strand bias as the three pillars of a
# defensible filter. This script replaces "trust the fixed numbers" with
# "look at THIS cohort's actual QUAL/DP/MQ distributions, add the missing
# strand-bias check, and show exactly which calls the naive fixed filter
# got wrong" — the same naive-vs-improved, side-by-side pattern Module 2
# Script 3 used for its low-count gene filter, applied here to variant QC.
#
# TEACHING NOTE — the strand-bias gap:
#   A true variant's supporting reads should come from both strands roughly
#   proportionally to the library's own strand balance. A candidate variant
#   whose ALT-supporting reads come almost entirely from ONE strand despite
#   adequate total depth is a classic PCR/mapping-artifact signature (the
#   same signal GATK's StrandOddsRatio/FS filters target) — and the fixed
#   QUAL/DP/MQ filter in master_snp_pipeline.sh cannot catch it, because
#   none of those three numbers encode strand information. This script
#   reads it directly from bcftools' own DP4 tag (ref-fwd, ref-rev,
#   alt-fwd, alt-rev read counts) instead of adding a new dependency.
#
# TEACHING NOTE — why the thresholds are recomputed per cohort, not reused:
#   QUAL/DP/MQ are call-quality metrics whose useful range depends on this
#   run's sequencing depth and organism — a fixed "DP >= 10" is defensible
#   for a 30x cohort and needlessly permissive for a 200x one. Percentile
#   thresholds anchor the filter to what this data actually looks like, and
#   are printed explicitly so the choice is auditable, not hidden.
#
# Input:  WORKDIR — the master_snp_pipeline.sh working directory, with
#         variants/<SRR>_raw.vcf (UNFILTERED — this script needs the full
#         candidate population, not what Phase 3 Step F already discarded).
# Output: variants/<SRR>_filtered_v2.vcf per sample; prints a side-by-side
#         naive-vs-data-driven comparison table.
# ==============================================================================

set -o pipefail

WORKDIR="${1:-/scratch/$(whoami)/master_snp_pipeline}"

echo "=================================================================="
echo " MODULE 4 / SCRIPT 3 — Designing a defensible variant filter"
echo "=================================================================="
echo "Working directory: ${WORKDIR}"
echo ""

command -v bcftools >/dev/null 2>&1 || { echo "[FAIL] bcftools not found on PATH"; exit 1; }

# Naive fixed thresholds, unchanged from master_snp_pipeline.sh Phase 3 Step F —
# kept here ONLY as the comparison baseline, never edited.
NAIVE_QUAL=20
NAIVE_DP=10
NAIVE_MQ=30

# Minimum ALT read count before a strand-bias check is even meaningful — at
# very low depth, "all N alt reads came from one strand" is expected by
# chance alone, not evidence of an artifact.
MIN_ALT_FOR_STRAND_CHECK=4
# An ALT allele supported by less than this fraction of reads from its
# weaker strand, despite having enough total ALT depth, is flagged.
MIN_STRAND_FRACTION=0.10

shopt -s nullglob
RAW_VCFS=("${WORKDIR}"/variants/*_raw.vcf)
shopt -u nullglob
[ ${#RAW_VCFS[@]} -gt 0 ] || { echo "[FAIL] no *_raw.vcf files found in ${WORKDIR}/variants — run master_snp_pipeline.sh first"; exit 1; }

printf "%-14s %8s %8s %8s %8s   %s\n" "sample" "p10_QUAL" "p10_DP" "p10_MQ" "strand_bias" "naive_vs_v2"

for RAW in "${RAW_VCFS[@]}"; do
    SRR=$(basename "${RAW}" _raw.vcf)
    V2="${WORKDIR}/variants/${SRR}_filtered_v2.vcf"

    # --- Pull QUAL/DP/MQ/DP4 for every candidate record in one pass ---------
    METRICS_TSV="${WORKDIR}/variants/${SRR}_raw_metrics.tsv"
    bcftools query -f '%QUAL\t%INFO/DP\t%INFO/MQ\t%INFO/DP4\n' "${RAW}" > "${METRICS_TSV}" 2>/dev/null

    # --- Data-driven thresholds: 10th percentile of each metric across this
    #     sample's own candidate calls -----------------------------------
    P10_QUAL=$(awk -F'\t' '$1!="."{print $1}' "${METRICS_TSV}" | sort -n | \
               awk '{a[NR]=$1} END{idx=int(0.10*NR); if(idx<1)idx=1; print (NR>0)?a[idx]:0}')
    P10_DP=$(awk -F'\t' '$2!="."{print $2}' "${METRICS_TSV}" | sort -n | \
             awk '{a[NR]=$1} END{idx=int(0.10*NR); if(idx<1)idx=1; print (NR>0)?a[idx]:0}')
    P10_MQ=$(awk -F'\t' '$3!="."{print $3}' "${METRICS_TSV}" | sort -n | \
             awk '{a[NR]=$1} END{idx=int(0.10*NR); if(idx<1)idx=1; print (NR>0)?a[idx]:0}')
    P10_QUAL=${P10_QUAL:-0}; P10_DP=${P10_DP:-0}; P10_MQ=${P10_MQ:-0}

    # Never go BELOW the naive threshold just because the data-driven
    # percentile is looser — the redesigned filter should be at least as
    # strict as the original, only smarter about strand bias on top.
    V2_QUAL=$(awk -v a="${P10_QUAL}" -v b="${NAIVE_QUAL}" 'BEGIN{print (a>b)?a:b}')
    V2_DP=$(awk -v a="${P10_DP}" -v b="${NAIVE_DP}" 'BEGIN{print (a>b)?a:b}')
    V2_MQ=$(awk -v a="${P10_MQ}" -v b="${NAIVE_MQ}" 'BEGIN{print (a>b)?a:b}')

    # --- Count strand-biased sites (for the report below) directly from DP4,
    #     independent of bcftools' own expression engine, so the printed
    #     count is trustworthy even if a given bcftools build parses vector
    #     subscripts differently ------------------------------------------
    STRAND_SITES="${WORKDIR}/variants/${SRR}_strand_biased_sites.tsv"
    bcftools query -f '%CHROM\t%POS\t%INFO/DP4\n' "${RAW}" 2>/dev/null | \
    awk -F'\t' -v min_alt="${MIN_ALT_FOR_STRAND_CHECK}" -v min_frac="${MIN_STRAND_FRACTION}" '
        {
            split($3, dp4, ",")
            alt_fwd = dp4[3]; alt_rev = dp4[4]
            alt_total = alt_fwd + alt_rev
            if (alt_total >= min_alt) {
                min_strand = (alt_fwd < alt_rev) ? alt_fwd : alt_rev
                if (min_strand < min_frac * alt_total) print $1"\t"$2
            }
        }' > "${STRAND_SITES}"
    N_STRAND_BIASED=$(wc -l < "${STRAND_SITES}")

    # --- Apply the v2 filter in one pass: naive-or-stricter QUAL/DP/MQ, plus
    #     strand bias read directly from INFO/DP4's own vector elements
    #     (DP4 = ref-fwd,ref-rev,alt-fwd,alt-rev; 0-based subscripts) --------
    STRAND_EXPR="((INFO/DP4[2]+INFO/DP4[3])>=${MIN_ALT_FOR_STRAND_CHECK} && (INFO/DP4[2]<${MIN_STRAND_FRACTION}*(INFO/DP4[2]+INFO/DP4[3]) || INFO/DP4[3]<${MIN_STRAND_FRACTION}*(INFO/DP4[2]+INFO/DP4[3])))"
    bcftools filter -s LowQual+StrandBias \
        -e "QUAL<${V2_QUAL} || INFO/DP<${V2_DP} || INFO/MQ<${V2_MQ} || ${STRAND_EXPR}" \
        "${RAW}" > "${V2}" 2>/dev/null

    # --- Side-by-side comparison: what did v2 change relative to Phase 3's
    #     naive filter? -------------------------------------------------
    NAIVE_VCF="${WORKDIR}/variants/${SRR}_filtered.vcf"
    NAIVE_PASS=$(bcftools view -H -f PASS "${NAIVE_VCF}" 2>/dev/null | wc -l)
    V2_PASS=$(bcftools view -H -f PASS "${V2}" 2>/dev/null | wc -l)
    DROPPED_BY_V2=$((NAIVE_PASS - V2_PASS))

    printf "%-14s %8s %8s %8s %8d   naive_PASS=%d -> v2_PASS=%d (%+d)\n" \
        "${SRR}" "${V2_QUAL}" "${V2_DP}" "${V2_MQ}" "${N_STRAND_BIASED}" \
        "${NAIVE_PASS}" "${V2_PASS}" $((-DROPPED_BY_V2))
done

echo ""
echo "Wrote <sample>_filtered_v2.vcf, <sample>_strand_biased_sites.tsv per sample."
echo "v2 thresholds are max(naive fixed value, this cohort's own 10th percentile) —"
echo "never looser than Phase 3's original filter, and strand-bias-aware on top."
echo "=================================================================="
