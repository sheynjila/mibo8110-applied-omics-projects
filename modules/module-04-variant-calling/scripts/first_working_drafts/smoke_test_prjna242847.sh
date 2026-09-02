#!/bin/bash
#SBATCH --job-name=snp_smoke_test
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --output=snp_smoke_%j.out
#SBATCH --error=snp_smoke_%j.err

# ==============================================================================
# SMOKE TEST WRAPPER for master_snp_pipeline.sh
# ==============================================================================
# WHY THIS SCRIPT EXISTS:
#   master_snp_pipeline.sh's own built-in fallback (used only when srr_list.txt
#   is missing) hardcodes SRR11092056/SRR11092057 as a "default test list."
#   Verified via NCBI eutils (esearch/esummary against the SRA and Assembly
#   databases, same method Module 3 used to verify its own accessions):
#
#     SRR11092056 = SARS-CoV-2 metagenomic RNA-Seq, bronchoalveolar lavage
#                   fluid, Wuhan Institute of Virology (BioProject
#                   PRJNA605983 / study SRP249613) - NOT bacterial WGS.
#     GCF_000006945.2 = Salmonella enterica subsp. enterica serovar
#                   Typhimurium str. LT2 (ASM694v2) - a bacterial reference.
#
#   Aligning viral metagenomic RNA-seq reads against a bacterial DNA
#   reference with BWA-MEM will not error - it will just produce near-zero,
#   scattered coverage and a variant call set that is pure alignment noise.
#   That is exactly the kind of silent, plausible-looking-but-wrong result
#   this module's Tier B scripts (esp. Script 2, mapping/reference-bias
#   assessment) are built to catch — but it is far better to never run the
#   smoke test against data that cannot possibly succeed in the first place.
#
#   Per this repository's "conserve what's already delivered" policy,
#   master_snp_pipeline.sh's internal fallback is NOT edited here (it is a
#   convenience default, not something any real run should depend on). This
#   wrapper supplies a real, verified, matching bacterial WGS srr_list.txt
#   instead, exactly the way Module 2's smoke-test wrapper supplied real
#   PRJNA1518998 accessions rather than relying on unverified defaults.
#
# WORKED-EXAMPLE COHORT (verified real, paired-end, Illumina MiSeq):
#   BioProject PRJNA242847 (NCBI Pathogen Detection / GenomeTrakr routine
#   public-health WGS surveillance), organism Salmonella enterica subsp.
#   enterica serovar Typhimurium - the same serovar as the LT2 reference
#   master_snp_pipeline.sh already downloads and indexes.
#     SRR40474618  (BioSample SAMN62852532)
#     SRR40472300  (BioSample SAMN62840582)
#     SRR40472280  (BioSample SAMN62840173)
#   Three isolates (not two) so the cohort merge and variant_analysis.R's
#   SNP-distance clustering have more than one pair to actually compare.
#
# Usage: sbatch smoke_test_prjna242847.sh
#        (run BEFORE committing srr_list.txt to a full cohort submission)
# ==============================================================================

set -e
set -o pipefail

echo "=========================================================="
echo " SNP PIPELINE SMOKE TEST — PRJNA242847 (3 real S. Typhimurium isolates)"
echo "=========================================================="

WORKDIR="/scratch/$(whoami)/master_snp_pipeline"
mkdir -p ${WORKDIR}
cd ${WORKDIR}

cat > srr_list.txt <<'EOF'
SRR40474618
SRR40472300
SRR40472280
EOF

echo "Wrote srr_list.txt:"
cat srr_list.txt

echo -e "\nLaunching master_snp_pipeline.sh against the verified 3-sample cohort...\n"
bash "$(dirname "$0")/master_snp_pipeline.sh"

echo "=========================================================="
echo " SMOKE TEST COMPLETE — inspect ${WORKDIR}/variants/merged_cohort.vcf"
echo " before submitting a larger cohort."
echo "=========================================================="
