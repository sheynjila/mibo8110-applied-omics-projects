#!/bin/bash
#SBATCH --job-name=amr_smoke_test
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=amr_smoke_%j.out
#SBATCH --error=amr_smoke_%j.err

# ==============================================================================
# SMOKE TEST WRAPPER for Final_master_end_to_end_amr.sh
# ==============================================================================
# WHY THIS SCRIPT EXISTS:
#   Final_master_end_to_end_amr.sh's own built-in fallback (used only when
#   srr_list.txt is missing) hardcodes SRR11092056/SRR11092057 - the same
#   "default test list" already verified elsewhere in this repository, via
#   NCBI eutils, as SARS-CoV-2 metagenomic RNA-seq, not bacterial isolate WGS
#   of any kind. SPAdes will still "assemble" viral RNA-seq reads into
#   something - just not a ~4.8-5.0 Mb Salmonella genome, and AMRFinderPlus
#   would be scanning a meaningless assembly for the wrong organism entirely.
#   Per this repository's "conserve what's already delivered" policy, that
#   fallback is NOT edited here; this wrapper supplies a real, verified
#   srr_list.txt instead, the same pattern used for every other module's
#   smoke test in this repository.
#
# WORKED-EXAMPLE COHORT: this wrapper deliberately reuses Module 4 and
# Module 5's own verified 5-isolate PRJNA242847 cohort (GenomeTrakr/USDA-FSIS
# real-world Salmonella enterica subsp. enterica serovar Typhimurium
# surveillance WGS) rather than sourcing a new dataset - the same organism
# AMR_ORGANISM defaults to here, real paired-end Illumina MiSeq data already
# verified via NCBI eutils, and it lets this case study's own AMR findings be
# read directly alongside Module 4's variant calls and Module 5's tree for
# the SAME five isolates, not a disconnected fourth dataset:
#     SRR40474618  (BioSample SAMN62852532)
#     SRR40472300  (BioSample SAMN62840582)
#     SRR40472280  (BioSample SAMN62840173)
#     SRR40426677  (BioSample SAMN62799163)
#     SRR40376221  (BioSample SAMN62746372)
#
# Usage: sbatch smoke_test_amr_prjna242847.sh
#        (run BEFORE committing srr_list.txt to a full cohort submission -
#        SPAdes assembly is the most expensive step in this repository's
#        entire script collection; do not learn about a bad input on the
#        full cohort.)
# ==============================================================================

set -e
set -o pipefail

echo "=========================================================="
echo " AMR PIPELINE SMOKE TEST — PRJNA242847"
echo " (5 real S. Typhimurium isolates, same cohort as Modules 4-5)"
echo "=========================================================="

export AMR_ORGANISM="Salmonella"
WORKDIR="/scratch/$(whoami)/end_to_end_amr_pipeline"
mkdir -p ${WORKDIR}
cd ${WORKDIR}

cat > srr_list.txt <<'EOF'
SRR40474618
SRR40472300
SRR40472280
SRR40426677
SRR40376221
EOF

echo "Wrote srr_list.txt:"
cat srr_list.txt

echo -e "\nLaunching Final_master_end_to_end_amr.sh (AMR_ORGANISM=${AMR_ORGANISM}) against the verified 5-sample cohort...\n"
bash "$(dirname "$0")/Final_master_end_to_end_amr.sh"

echo "=========================================================="
echo " SMOKE TEST COMPLETE — inspect ${WORKDIR}/amr_reports/ALL_SAMPLES_AMR_SUMMARY.tsv"
echo " and ${WORKDIR}/quast_qc/ (expect ~4.8-5.0 Mb assemblies - see Tier A's own"
echo " 'CAPSTONE METHODOLOGY SUMMARY' footer for what a ~9.5 Mb assembly would mean)"
echo " before submitting a larger cohort."
echo "=========================================================="
