#!/bin/bash
#SBATCH --job-name=microbiome_smoke_test
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=microbiome_smoke_%j.out
#SBATCH --error=microbiome_smoke_%j.err

# ==============================================================================
# SMOKE TEST WRAPPER for master_microbiome_pipeline.sh (shotgun route)
# ==============================================================================
# WHY THIS SCRIPT EXISTS:
#   master_microbiome_pipeline.sh's own built-in fallback (used only when
#   srr_list.txt is missing) hardcodes SRR11092056/SRR11092057 as a "default
#   test list" - the exact same two accessions Module 4's smoke test already
#   verified (via NCBI eutils) are SARS-CoV-2 metagenomic RNA-seq from a
#   bronchoalveolar lavage sample (BioProject PRJNA605983), not a bacterial
#   or environmental shotgun metagenome of any kind. Kraken2 will still
#   "run" against them - it will just classify almost nothing at the
#   bacterial/archaeal level and produce a report that looks like an empty
#   or failed community profile rather than the obviously-wrong-input
#   problem it actually is. Per this repository's "conserve what's already
#   delivered" policy, that fallback is NOT edited here; this wrapper
#   supplies a real, verified srr_list.txt instead, the same pattern Module
#   4's smoke test used for the identical underlying default.
#
# WORKED-EXAMPLE COHORT (verified real, paired-end, Illumina NovaSeq 6000,
# via NCBI eutils esearch/esummary against BioProject 1520859):
#   BioProject PRJNA1520859 - "Shotgun metagenomic analysis of rhizosphere
#   microbiomes from healthy and Fusarium-diseased tomato plants" (rhizosphere
#   metagenome; Nanjing Agricultural University).
#     SRR40412948  Fusarium-diseased tomato rhizosphere, replicate 3
#     SRR40412949  Fusarium-diseased tomato rhizosphere, replicate 2
#     SRR40412951  Healthy tomato rhizosphere, replicate 3
#   Two conditions (not three replicates of one) so the cohort abundance
#   table this module's Tier B scripts consume has a real healthy-vs-diseased
#   contrast to report on, not just three near-identical samples.
#
#   This cohort is plant-associated (rhizosphere = soil directly surrounding
#   living tomato roots), which is exactly the sample type
#   master_microbiome_pipeline.sh's optional HOST_REFERENCE host-depletion
#   step targets. To exercise that feature, set HOST_REFERENCE to a tomato
#   genome FASTA (e.g. Solanum lycopersicum ITAG / RefSeq GCF_000188115.5)
#   before running this smoke test:
#     export HOST_REFERENCE=/path/to/Slycopersicum_genome.fna
#   Left unset, host depletion is skipped and every trimmed read goes
#   straight to Kraken2 - a legitimate run, just not one that demonstrates
#   the host-depletion feature. Rhizosphere soil (unlike root tissue itself)
#   is typically a modest, not dominant, host-DNA fraction - do not expect
#   the dramatic host percentages a root- or gut-tissue sample would show.
#
# Usage: sbatch smoke_test_shotgun_prjna1520859.sh
#        (run BEFORE committing srr_list.txt to a full cohort submission)
# ==============================================================================

set -e
set -o pipefail

echo "=========================================================="
echo " MICROBIOME PIPELINE SMOKE TEST (shotgun) — PRJNA1520859"
echo " (2 Fusarium-diseased + 1 healthy tomato rhizosphere sample)"
echo "=========================================================="

WORKDIR="/scratch/$(whoami)/master_microbiome_pipeline"
mkdir -p ${WORKDIR}
cd ${WORKDIR}

cat > srr_list.txt <<'EOF'
SRR40412948
SRR40412949
SRR40412951
EOF

echo "Wrote srr_list.txt:"
cat srr_list.txt

echo -e "\nLaunching master_microbiome_pipeline.sh against the verified 3-sample cohort..."
echo "HOST_REFERENCE=${HOST_REFERENCE:-<unset — host depletion will be skipped>}"
echo ""
bash "$(dirname "$0")/master_microbiome_pipeline.sh"

echo "=========================================================="
echo " SMOKE TEST COMPLETE — inspect ${WORKDIR}/bracken_output/ALL_SAMPLES_MICROBIOME.tsv"
echo " before submitting a larger cohort."
echo "=========================================================="
