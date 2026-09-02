#!/bin/bash
#SBATCH --job-name=mod3_choose_strategy
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:10:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
###############################################################################
# Module 3, Script 1 — Choose an Assembly Strategy From the Reads Themselves
#
# NEW IDEA on top of the earlier modules: nothing here assembles anything yet.
# The curriculum's first Module 3 learning outcome is "choose among short-read,
# long-read, and hybrid assembly strategies" — that is a decision, and this
# script makes the decision explicit and justified from evidence in the FASTQ
# files, instead of hardcoding "run SPAdes" the way a copy-pasted pipeline
# script typically does.
#
# TEACHING NOTE: the existing repo pipeline
# (modules/supplemental-case-studies/scripts/Final_master_end_to_end_amr.sh)
# always calls SPAdes with --isolate. That is a *correct* default for
# short-read Illumina isolate data, but it is a default the author chose once
# for one dataset — it is not a rule that generalizes to long-read or hybrid
# data. This script exists so a student learns to check *which* situation
# they are in before picking the tool, rather than inheriting someone else's
# choice unexamined.
#
# WHAT THIS SCRIPT DOES:
#   1. Reads read-length statistics directly out of a FASTQ file (first N
#      records), instead of trusting a filename or a platform label.
#   2. Classifies the run as SHORT_READ, LONG_READ, or UNCERTAIN using a
#      simple, documented length threshold (short reads cluster tightly
#      around ~50-300bp; long reads have a long tail into the thousands).
#   3. If both a short-read and a long-read FASTQ are supplied, recommends
#      HYBRID and explains why (short reads correct long-read errors within
#      an assembly; long reads resolve repeats short reads cannot span).
#   4. Prints the exact next script to run for the chosen strategy — it does
#      NOT invoke SPAdes or Flye itself. Decision and execution are kept in
#      separate scripts so the decision is inspectable on its own.
#
# USAGE:
#   ./1_Choose_Assembly_Strategy.sh --short reads_1.fastq [--long long.fastq]
#
# REAL WORKED-EXAMPLE DATASETS FOR THIS MODULE (see README.md for citations):
#   Short-read: SRR1770413  (E. coli K-12, Illumina MiSeq, paired-end, PRJNA272917)
#   Long-read:  SRR39619343 (E. coli K-12 MG1655, Oxford Nanopore, PRJNA1490942)
###############################################################################

set -euo pipefail
mkdir -p logs

usage() {
    echo "Usage: $0 --short <R1.fastq> [--long <long_reads.fastq>]" >&2
    exit 1
}

SHORT_FASTQ=""
LONG_FASTQ=""
SAMPLE_N=2000   # how many reads to inspect — enough to be representative, cheap to scan

while [[ $# -gt 0 ]]; do
    case "$1" in
        --short) SHORT_FASTQ="$2"; shift 2 ;;
        --long)  LONG_FASTQ="$2";  shift 2 ;;
        *) usage ;;
    esac
done

if [[ -z "$SHORT_FASTQ" && -z "$LONG_FASTQ" ]]; then
    usage
fi

# --- Read-length summary from a FASTQ file (sequence line is line 2 of every
#     4-line record). Prints: count, mean, min, max.
summarize_lengths() {
    local fastq="$1"
    local n="$2"
    awk -v maxrec="$n" '
        NR % 4 == 2 {
            len = length($0)
            sum += len
            if (NR == 2 || len < min) min = len
            if (len > max) max = len
            count++
            if (count >= maxrec) exit
        }
        END {
            if (count == 0) { print "0 0 0 0"; exit }
            printf "%d %.1f %d %d\n", count, sum/count, min, max
        }
    ' "$fastq"
}

classify() {
    # Argument: mean read length. Threshold rationale: Illumina short reads
    # are essentially always <= 300bp (even 2x300 MiSeq kits); anything with
    # a mean length solidly above that is a long-read platform (ONT/PacBio),
    # where individual reads commonly run 1,000-50,000+ bp.
    local mean_len="$1"
    if awk -v m="$mean_len" 'BEGIN { exit !(m <= 350) }'; then
        echo "SHORT_READ"
    elif awk -v m="$mean_len" 'BEGIN { exit !(m > 350) }'; then
        echo "LONG_READ"
    else
        echo "UNCERTAIN"
    fi
}

echo "=========================================================="
echo " MODULE 3 — ASSEMBLY STRATEGY SELECTION"
echo "=========================================================="

SHORT_CLASS=""
LONG_CLASS=""

if [[ -n "$SHORT_FASTQ" ]]; then
    if [[ ! -f "$SHORT_FASTQ" ]]; then
        echo "ERROR: --short file not found: $SHORT_FASTQ" >&2
        exit 1
    fi
    read -r S_COUNT S_MEAN S_MIN S_MAX <<< "$(summarize_lengths "$SHORT_FASTQ" "$SAMPLE_N")"
    SHORT_CLASS=$(classify "$S_MEAN")
    echo ">> $SHORT_FASTQ: n=$S_COUNT (sampled) mean=${S_MEAN}bp min=${S_MIN}bp max=${S_MAX}bp -> classified as $SHORT_CLASS"
fi

if [[ -n "$LONG_FASTQ" ]]; then
    if [[ ! -f "$LONG_FASTQ" ]]; then
        echo "ERROR: --long file not found: $LONG_FASTQ" >&2
        exit 1
    fi
    read -r L_COUNT L_MEAN L_MIN L_MAX <<< "$(summarize_lengths "$LONG_FASTQ" "$SAMPLE_N")"
    LONG_CLASS=$(classify "$L_MEAN")
    echo ">> $LONG_FASTQ: n=$L_COUNT (sampled) mean=${L_MEAN}bp min=${L_MIN}bp max=${L_MAX}bp -> classified as $LONG_CLASS"
fi

echo ""
echo "---------------------------------------------------------"
echo " RECOMMENDATION"
echo "---------------------------------------------------------"

if [[ -n "$SHORT_FASTQ" && -n "$LONG_FASTQ" ]]; then
    echo "Both a short-read and a long-read FASTQ were supplied -> recommend HYBRID."
    echo "Rationale: short reads (low error rate, ~0.1-1%) correct the higher"
    echo "per-base error rate typical of long-read platforms, while the long"
    echo "reads span repeats that short reads alone cannot resolve into a single"
    echo "contig. Run BOTH 2_Run_SPAdes_Short_Read_Assembly.sh (with --hybrid"
    echo "pointing at the long reads) and inspect the resulting graph with"
    echo "4_Interpret_Assembly_Graph.sh before deciding a pure-long-read Flye"
    echo "assembly (3_Run_Flye_Long_Read_Assembly.sh) is unnecessary."
elif [[ "$SHORT_CLASS" == "SHORT_READ" ]]; then
    echo "Short-read data only -> recommend SPAdes."
    echo "Next step: 2_Run_SPAdes_Short_Read_Assembly.sh"
elif [[ "$LONG_CLASS" == "LONG_READ" ]]; then
    echo "Long-read data only -> recommend Flye."
    echo "Next step: 3_Run_Flye_Long_Read_Assembly.sh"
else
    echo "WARNING: could not confidently classify the supplied read set(s) as"
    echo "short- or long-read from length statistics alone. Do NOT default to"
    echo "SPAdes here just because it is the more familiar tool — check the"
    echo "sequencing platform metadata (SRA run's 'Platform'/'Model' fields, or"
    echo "the FASTQ header if the run was basecalled locally) before choosing."
fi

echo ""
echo "This script made a recommendation only — it did not assemble anything."
echo "=========================================================="
