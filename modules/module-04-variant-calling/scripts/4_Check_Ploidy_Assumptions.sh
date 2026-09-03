#!/bin/bash
# ==============================================================================
# MODULE 4 · SCRIPT 4 of 7 — CHECK PLOIDY ASSUMPTIONS
# ==============================================================================
# ONE NEW IDEA on top of Script 3: Script 3 redesigned WHICH calls survive
# filtering. This script questions something upstream of filtering
# entirely — whether the calls were even generated under the right genetic
# model in the first place.
#
# TEACHING NOTE — the bug, stated precisely:
#   master_snp_pipeline.sh's Phase 3 Step E calls variants with
#   `bcftools call -mv -O v` and no --ploidy flag. bcftools call's default
#   ploidy, absent an explicit --ploidy/-p argument or --ploidy-file, is
#   DIPLOID (2) genome-wide — a human/eukaryote-shaped assumption. The
#   reference this pipeline downloads and indexes (GCF_000006945.2,
#   Salmonella enterica serovar Typhimurium LT2, ASM694v2) is a single
#   circular bacterial chromosome: genuinely HAPLOID. Calling a haploid
#   organism under a diploid model does not crash — bcftools will happily
#   emit heterozygous genotypes (0/1, 1/2, ...) at positions where read
#   noise, a misalignment, or a genuinely mixed/contaminated culture
#   produces two different alleles at meaningful frequency. For a true
#   haploid genome, every real GT call should be homozygous (0/0 or 1/1);
#   any 0/1-style call is not "heterozygosity" in the biological sense —
#   it is either a calling-model artifact or evidence worth investigating
#   on its own (contamination, mixed infection), not something a diploid
#   genotype model was ever designed to represent correctly.
#
# TEACHING NOTE — why this matters for the OTHER learning outcomes:
#   Script 3's filter design and Script 5's cohort Ti/Tv sanity check are
#   both computed from genotype calls. If those calls were made under the
#   wrong ploidy, their apparent "quality" can look fine while their
#   genetic interpretation is wrong — this is exactly the kind of problem
#   that produces a plausible-looking but incorrect result, the same
#   failure mode Module 2 Script 1's header warns about for count/metadata
#   mismatches.
#
# Input:  WORKDIR — the master_snp_pipeline.sh working directory, with
#         variants/<SRR>_raw.vcf (diploid-called) and the indexed reference.
# Output: variants/<SRR>_raw_haploid.vcf (re-called with --ploidy 1); prints
#         a before/after heterozygous-call count per sample. Downstream
#         scripts should be re-run against the *_raw_haploid.vcf files for a
#         ploidy-correct result — see this module's README/RUNBOOK.
# ==============================================================================

set -o pipefail

WORKDIR="${1:-/scratch/$(whoami)/master_snp_pipeline}"
REF="${WORKDIR}/ref/genome.fna"

echo "=================================================================="
echo " MODULE 4 / SCRIPT 4 — Ploidy assumption check and correction"
echo "=================================================================="
echo "Working directory: ${WORKDIR}"
echo ""
echo "Reference organism (per this pipeline's own download URL) is a single"
echo "bacterial chromosome (HAPLOID). master_snp_pipeline.sh's Phase 3 Step E"
echo "calls variants with 'bcftools call -mv' and no --ploidy flag, which"
echo "defaults to DIPLOID. Checking for genotype calls that assumption cannot"
echo "actually justify..."
echo ""

command -v bcftools >/dev/null 2>&1 || { echo "[FAIL] bcftools not found on PATH"; exit 1; }

shopt -s nullglob
RAW_VCFS=("${WORKDIR}"/variants/*_raw.vcf)
shopt -u nullglob
[ ${#RAW_VCFS[@]} -gt 0 ] || { echo "[FAIL] no *_raw.vcf files found in ${WORKDIR}/variants — run master_snp_pipeline.sh first"; exit 1; }

printf "%-14s %14s %20s %s\n" "sample" "het_calls(dip)" "het_calls(hap,expect 0)" "status"

TOTAL_HET_DIPLOID=0
for RAW in "${RAW_VCFS[@]}"; do
    SRR=$(basename "${RAW}" _raw.vcf)
    BAM="${WORKDIR}/alignment/${SRR}_dedup.bam"
    HAP_VCF="${WORKDIR}/variants/${SRR}_raw_haploid.vcf"

    # --- Count heterozygous GT calls under the pipeline's actual (diploid)
    #     output — GT field where the two alleles differ, e.g. 0/1, 1/2 ----
    HET_DIPLOID=$(bcftools query -f '[%GT]\n' "${RAW}" 2>/dev/null | \
        awk -F'[/|]' '{if ($1!="." && $2!="" && $1!=$2) c++} END{print c+0}')
    TOTAL_HET_DIPLOID=$((TOTAL_HET_DIPLOID + HET_DIPLOID))

    # --- Re-call from scratch with the correct haploid model. mpileup is
    #     re-run rather than reused because Phase 3 Step E piped mpileup
    #     directly into call without saving the intermediate BCF ------------
    if [ -s "${BAM}" ] && [ -s "${REF}" ]; then
        bcftools mpileup -O b -f "${REF}" "${BAM}" 2>/dev/null | \
        bcftools call -mv --ploidy 1 -O v -o "${HAP_VCF}" 2>/dev/null
        HET_HAPLOID=$(bcftools query -f '[%GT]\n' "${HAP_VCF}" 2>/dev/null | \
            awk -F'[/|]' '{if ($1!="." && $2!="" && $1!=$2) c++} END{print c+0}')
    else
        HET_HAPLOID="n/a (BAM/ref missing)"
    fi

    if [ "${HET_DIPLOID}" -gt 0 ]; then
        STATUS="diploid default produced ${HET_DIPLOID} impossible het call(s) for a haploid genome"
    else
        STATUS="no het calls even under the wrong default (does not mean ploidy was right — just lucky here)"
    fi
    printf "%-14s %14s %20s   %s\n" "${SRR}" "${HET_DIPLOID}" "${HET_HAPLOID}" "${STATUS}"
done

echo ""
echo "Wrote <sample>_raw_haploid.vcf (bcftools call --ploidy 1) per sample."
if [ "${TOTAL_HET_DIPLOID}" -gt 0 ]; then
    echo "[WARN] ${TOTAL_HET_DIPLOID} total heterozygous call(s) found under the pipeline's"
    echo "       default diploid model, across a genome this pipeline itself treats as"
    echo "       haploid. Re-run Script 3 (filter design) and Script 5 (cohort QC)"
    echo "       against *_raw_haploid.vcf, not *_raw.vcf, before drawing conclusions."
else
    echo "No heterozygous calls found in this cohort's diploid-default output — but the"
    echo "*_raw_haploid.vcf files are still the ploidy-correct calls to build on, since"
    echo "absence of het calls this run is not a guarantee for every cohort."
fi
echo "=================================================================="
