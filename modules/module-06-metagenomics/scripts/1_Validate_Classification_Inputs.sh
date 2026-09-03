#!/bin/bash
# ==============================================================================
# MODULE 6 · SCRIPT 1 of 5 — VALIDATE CLASSIFICATION INPUTS
# ==============================================================================
# ONE NEW IDEA on top of master_microbiome_pipeline.sh: that script produces
# a Kraken2 report and a Bracken abundance table per sample but never checks
# that they actually belong to each other, or that Bracken's re-estimation
# ran against the report Kraken2 actually produced for THAT sample - the
# same "check the pipeline's own output before trusting it" habit Module 4
# Script 1 applied to BAM/VCF pairs.
#
# TEACHING NOTE — the strategy choice this script's header documents rather
# than automates (Learning outcome 1: "choose the right strategy for the
# data"):
#   This module ships TWO Tier A pipelines - master_microbiome_pipeline.sh
#   (shotgun, Kraken2/Bracken) and master_amplicon_pipeline.sh +
#   amplicon_dada2_analysis.R (amplicon, DADA2) - deliberately, because
#   there is no single "right" pipeline for a mixed-community sample, only a
#   right pipeline FOR A GIVEN QUESTION AND BUDGET:
#     - Shotgun sees every gene in every read (functional profiling is
#       possible downstream, not just taxonomy) and is not limited to
#       whatever a "universal" marker-gene primer pair happens to amplify -
#       but it costs more per sample for equivalent taxonomic depth, and its
#       taxonomic call is only as good as Kraken2's reference database
#       (Script 4 below quantifies that limit directly).
#     - Amplicon (16S/ITS) is cheaper per sample and gives fine within-marker
#       resolution for whatever the primers DO amplify, but sees only that
#       one marker gene - no functional profiling - and is blind to whatever
#       the primer pair itself amplifies poorly (a real, silent gap, not a
#       hypothetical one; there is no script in this module that can detect
#       a taxon the primers never touched).
#   This module's own worked example (see RUNBOOK) ran BOTH routes against
#   comparable real cohorts specifically so this tradeoff is something you
#   can see in this module's own output, not just read about here. Scripts
#   2-5 below walk the SHOTGUN route's output in depth, because it produces
#   the concrete per-read classification-yield numbers Script 4's
#   database-completeness point needs; the amplicon route is validated by
#   its own smoke test (smoke_test_amplicon_prjna917645.sh) and DADA2's
#   built-in per-step diagnostics (filterAndTrim's retained-read counts,
#   removeBimeraDenovo's kept-read fraction), not by a separate Tier B
#   walkthrough.
#
# Input:  WORKDIR — the master_microbiome_pipeline.sh working directory,
#         containing kraken_output/<SRR>_report.txt and
#         bracken_output/<SRR>_bracken.tsv for one or more samples.
# Output: prints a PASS/FAIL/WARN validation report to the console per
#         sample; stops with a clear error on the first hard failure.
# ==============================================================================

set -o pipefail

WORKDIR="${1:-/scratch/$(whoami)/master_microbiome_pipeline}"

echo "=================================================================="
echo " MODULE 6 / SCRIPT 1 — Classification input validation"
echo "=================================================================="
echo "Working directory: ${WORKDIR}"
echo ""

fail() { echo "[FAIL] $1"; exit 1; }
pass() { echo "[PASS] $1"; }

shopt -s nullglob
REPORTS=("${WORKDIR}"/kraken_output/*_report.txt)
shopt -u nullglob

[ ${#REPORTS[@]} -gt 0 ] || fail "no *_report.txt files found in ${WORKDIR}/kraken_output — run master_microbiome_pipeline.sh first"

HARD_FAILS=0
HOST_DEPLETED_SAMPLES=0
for REPORT in "${REPORTS[@]}"; do
    SRR=$(basename "${REPORT}" | sed 's/_report\.txt$//')
    BRACKEN="${WORKDIR}/bracken_output/${SRR}_bracken.tsv"
    HOST_LOG="${WORKDIR}/qc/${SRR}_host_depletion.log"
    echo "------------------------------------------------------------"
    echo " Sample: ${SRR}"
    echo "------------------------------------------------------------"

    if [ ! -s "${REPORT}" ]; then
        echo "[FAIL] Kraken2 report not found or empty: ${REPORT}"
        HARD_FAILS=$((HARD_FAILS + 1))
        continue
    fi
    N_REPORT_LINES=$(wc -l < "${REPORT}")
    if [ "${N_REPORT_LINES}" -lt 2 ]; then
        echo "[FAIL] Kraken2 report has ${N_REPORT_LINES} line(s) — too few to contain even the standard 'unclassified' + 'root' rows"
        HARD_FAILS=$((HARD_FAILS + 1))
        continue
    fi
    N_COLS=$(head -n1 "${REPORT}" | awk -F'\t' '{print NF}')
    if [ "${N_COLS}" -ne 6 ]; then
        echo "[FAIL] Kraken2 report does not have the expected 6 tab-separated columns (found ${N_COLS}) — is this really a --report output?"
        HARD_FAILS=$((HARD_FAILS + 1))
        continue
    fi
    echo "[PASS] Kraken2 report exists, non-empty, and structurally well-formed (${N_REPORT_LINES} taxa, 6 columns)"

    if [ ! -s "${BRACKEN}" ]; then
        echo "[FAIL] Bracken output not found or empty: ${BRACKEN} — Kraken2 succeeded for ${SRR} but Bracken did not (or has not yet run)"
        HARD_FAILS=$((HARD_FAILS + 1))
        continue
    fi
    BRACKEN_HEADER=$(head -n1 "${BRACKEN}")
    case "${BRACKEN_HEADER}" in
        *fraction_total_reads*) ;;
        *) echo "[FAIL] ${BRACKEN} does not look like a Bracken abundance table (missing 'fraction_total_reads' column)"
           HARD_FAILS=$((HARD_FAILS + 1))
           continue
           ;;
    esac
    N_TAXA=$(($(wc -l < "${BRACKEN}") - 1))
    echo "[PASS] Bracken output exists and is well-formed (${N_TAXA} species-level taxon/taxa re-estimated)"

    if [ -s "${HOST_LOG}" ]; then
        HOST_DEPLETED_SAMPLES=$((HOST_DEPLETED_SAMPLES + 1))
        echo "[PASS] host-depletion log present for ${SRR} — HOST_REFERENCE was set for this run (see Script 2 for the actual host fraction)"
    else
        echo "[WARN] no host-depletion log for ${SRR} — HOST_REFERENCE was unset for this run; every trimmed read went straight to Kraken2"
    fi
done

echo ""
echo "=================================================================="
if [ ${HARD_FAILS} -gt 0 ]; then
    fail "${HARD_FAILS} sample(s) failed validation — fix before proceeding to Script 2"
fi
echo " VALIDATION COMPLETE for ${#REPORTS[@]} sample(s) (${HOST_DEPLETED_SAMPLES} host-depleted) — safe to proceed to Script 2"
echo "=================================================================="
