#!/bin/bash
# ==============================================================================
# MODULE 5 · SCRIPT 1 of 6 — EXTRACT A CORE-SNP ALIGNMENT FROM A COHORT VCF
# ==============================================================================
# ONE NEW IDEA on top of Module 4: Module 4 produced a cohort VCF and a
# pairwise SNP-distance matrix (variant_analysis.R), but neither of those is
# an ALIGNMENT — the input format tree-inference tools like IQ-TREE actually
# need. This script is the missing link this module's own README named
# explicitly: turn Module 4's cohort VCF into a FASTA alignment suitable for
# tree inference, using nothing beyond bcftools (already a dependency of
# this repo's variant-calling module).
#
# TEACHING NOTE — what a "core-SNP alignment" is, and is not:
#   Each output sequence is built ONLY from the reference-relative genotype
#   at every site the cohort's variant caller found polymorphic — invariant
#   (monomorphic) positions are never included. This is standard practice in
#   bacterial outbreak phylogenomics (the same approach tools like
#   snp-sites/Lyve-SET/CFSAN SNP Pipeline use) and keeps the alignment small
#   enough to build by hand from a VCF with no external aligner — but it is
#   NOT a whole-genome alignment, and its branch lengths do NOT represent
#   ordinary per-site substitution rates, because invariant sites have been
#   deliberately excluded (ascertainment bias). Script 3 corrects for this
#   explicitly with IQ-TREE's `+ASC` model correction rather than silently
#   treating the tree as if it were built from a full genome alignment.
#
# TEACHING NOTE — heterozygous calls become 'N', not an ambiguity code:
#   This reference is haploid (Module 4 Script 4). A het call surviving into
#   this cohort VCF is evidence of the diploid-default artifact Module 4
#   Script 4 exists to catch, not real biological heterozygosity — treating
#   it as a trustworthy IUPAC ambiguity code would smuggle that artifact
#   into the tree. It is coded 'N' (missing/uncertain) instead, and Script 2
#   reports how often this happens per sample.
#
# Input:  $1 = cohort VCF (use Module 4's ploidy-corrected
#         merged_cohort_haploid.vcf — see this module's README for why).
#         $2 = output directory (default: current directory).
# Output: core_snp_alignment.fasta — one sequence per cohort sample, plus a
#         "Reference" sequence (the REF allele at every extracted site, for
#         rooting) — and snp_sites.tsv, the intermediate per-site table.
# ==============================================================================

set -o pipefail

VCF="${1:?Usage: $0 <cohort.vcf> [output_dir]}"
OUTDIR="${2:-.}"
mkdir -p "${OUTDIR}"

echo "=================================================================="
echo " MODULE 5 / SCRIPT 1 — Extracting a core-SNP alignment"
echo "=================================================================="
echo "Cohort VCF: ${VCF}"
echo "Output dir: ${OUTDIR}"
echo ""

command -v bcftools >/dev/null 2>&1 || { echo "[FAIL] bcftools not found on PATH"; exit 1; }
[ -s "${VCF}" ] || { echo "[FAIL] cohort VCF not found or empty: ${VCF}"; exit 1; }

SAMPLE_NAMES=$(bcftools query -l "${VCF}")
N_SAMPLES=$(echo "${SAMPLE_NAMES}" | grep -c .)
[ "${N_SAMPLES}" -ge 3 ] || echo "[WARN] only ${N_SAMPLES} sample(s) in ${VCF} — a tree needs at least 3 taxa to show any topology at all, and at least 4 to have an internal branch worth a support value."

SITES_TSV="${OUTDIR}/snp_sites.tsv"
bcftools view -m2 -M2 -v snps "${VCF}" 2>/dev/null | \
bcftools query -f '%REF\t%ALT[\t%GT]\n' - > "${SITES_TSV}"

N_SITES=$(wc -l < "${SITES_TSV}")
if [ "${N_SITES}" -eq 0 ]; then
    echo "[FAIL] no biallelic SNP sites found in ${VCF} (after excluding indels/multiallelic sites) — nothing to align"
    exit 1
fi
echo "[PASS] ${N_SITES} biallelic SNP sites found across ${N_SAMPLES} samples"

OUT_FASTA="${OUTDIR}/core_snp_alignment.fasta"
awk -v n="${N_SAMPLES}" -v names="${SAMPLE_NAMES}" '
BEGIN {
    split(names, name_arr, "\n")
}
{
    ref = $1; alt = $2
    for (i = 1; i <= n; i++) {
        gt = $(i + 2)
        gsub(/\|/, "/", gt)
        n_alleles = split(gt, alleles, "/")
        base = "N"
        if (n_alleles == 1) {
            if (alleles[1] == "0") base = ref
            else if (alleles[1] == "1") base = alt
        } else if (n_alleles == 2) {
            if (alleles[1] == alleles[2]) {
                if (alleles[1] == "0") base = ref
                else if (alleles[1] == "1") base = alt
            } else {
                base = "N"   # heterozygous under a haploid model - not trusted, see header note
                het_count[i]++
            }
        }
        seq[i] = seq[i] base
    }
    seq["REF"] = seq["REF"] ref
    n_sites++
}
END {
    for (i = 1; i <= n; i++) {
        printf(">%s\n%s\n", name_arr[i], seq[i])
    }
    printf(">Reference\n%s\n", seq["REF"])
}
' "${SITES_TSV}" > "${OUT_FASTA}"

echo "Wrote ${OUT_FASTA} (${N_SITES} sites x $((N_SAMPLES + 1)) sequences, including the Reference row)"
echo "Wrote ${SITES_TSV}"
echo ""
echo "This is a SNP-ONLY alignment (ascertainment-biased — see header note)."
echo "Run Script 2 before trusting it, then Script 3 with IQ-TREE's +ASC"
echo "model correction, not a plain substitution model."
echo "=================================================================="
