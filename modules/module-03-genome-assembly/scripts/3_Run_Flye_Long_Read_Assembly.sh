#!/bin/bash
#SBATCH --job-name=mod3_flye
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
###############################################################################
# Module 3, Script 3 — Run Flye on Long Reads, With the Read-Type Flag Justified
#
# NEW IDEA on top of script 2: script 2 taught documented-parameter SPAdes
# calls for short reads. This script does the same for the long-read branch,
# and the one new idea it adds is that Flye's most consequential flag isn't
# a performance knob — it's a statement about your data's error PROFILE, and
# guessing it wrong silently degrades the assembly rather than erroring out.
#
# WHY THE --nano-raw / --nano-hq / --pacbio-raw DISTINCTION MATTERS:
# Flye's read-type flag tells its error-correction model what kind of noise
# to expect. Older ONT chemistries and default (non-"super-accuracy")
# basecalling produce ~5-15% error reads -> --nano-raw. Reads basecalled with
# ONT's higher-accuracy (Q20+) models are cleaner -> --nano-hq, which uses a
# lighter correction pass appropriate for ~1-5% error reads. Using --nano-raw
# on --nano-hq-quality reads over-corrects and can erase real variants;
# using --nano-hq on truly noisy --nano-raw-quality reads under-corrects and
# leaves more errors in the consensus. This script requires the read type as
# an explicit argument for exactly that reason — it will NOT guess for you.
#
# USAGE:
#   ./3_Run_Flye_Long_Read_Assembly.sh <long_reads.fastq> <output_dir> <read_type> [--threads N] [--low-mem]
#     <read_type> must be one of: nano-raw | nano-hq | nano-corr | pacbio-raw | pacbio-hifi
#
# MEMORY NOTE (a real, verified finding, not a hypothetical): Flye's k-mer
# counting step allocates a flat array sized 4^k / 2 bytes, and that size
# depends ONLY on the k-mer size k, NOT on your genome size. For --nano-raw's
# default k=17, that is 4^17 / 2 = 8,589,934,592 bytes (~8.6GB) -- allocated
# up front even to assemble a small bacterial genome. If your assembly node
# is memory-constrained (a laptop, a small VM, a shared login node), this
# step can fail with std::bad_alloc despite the input data itself being tiny.
# This module's own development sandbox (7.8GB total RAM) hit exactly this
# failure on Flye's own official bundled test dataset. The fix is --low-mem
# below, which drops k to 15 (4^15/2 = ~537MB) via Flye's --extra-params
# passthrough. On a real HPC node with the 64GB this script's SBATCH header
# requests, you do NOT need --low-mem -- use Flye's tuned default (k=17) for
# real assemblies, since a larger k gives better repeat resolution.
#
# REAL WORKED EXAMPLE: SRR39619343 (E. coli K-12 MG1655, Oxford Nanopore
# PromethION, PRJNA1490942) — see README.md for the full citation.
###############################################################################

set -euo pipefail
mkdir -p logs

SCRIPT_DIR="${PIPELINE_SCRIPT_DIR:-${SLURM_SUBMIT_DIR:+${SLURM_SUBMIT_DIR}/modules/module-03-genome-assembly/scripts}}"
SCRIPT_DIR="${SCRIPT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)}"
source "${SCRIPT_DIR}/load_modules.sh"

module purge
load_flye || exit 1

VALID_TYPES=("nano-raw" "nano-hq" "nano-corr" "pacbio-raw" "pacbio-hifi")

usage() {
    echo "Usage: $0 <long_reads.fastq> <output_dir> <read_type> [--threads N] [--low-mem]" >&2
    echo "  <read_type> must be one of: ${VALID_TYPES[*]}" >&2
    exit 1
}

if [[ $# -lt 3 ]]; then
    usage
fi

LONG_FASTQ="$1"; OUTDIR="$2"; READ_TYPE="$3"
shift 3
THREADS=16
LOW_MEM=0
while [[ $# -gt 0 ]]; do
    case "$1" in
        --threads) THREADS="$2"; shift 2 ;;
        --low-mem) LOW_MEM=1; shift ;;
        *) echo "ERROR: unrecognized argument: $1" >&2; exit 1 ;;
    esac
done

if [[ ! -f "$LONG_FASTQ" ]]; then
    echo "ERROR: input read file not found: $LONG_FASTQ" >&2
    exit 1
fi

VALID=0
for t in "${VALID_TYPES[@]}"; do
    [[ "$t" == "$READ_TYPE" ]] && VALID=1
done
if [[ "$VALID" -ne 1 ]]; then
    echo "ERROR: '$READ_TYPE' is not a recognized Flye read type." >&2
    usage
fi

mkdir -p "$OUTDIR"

echo "=========================================================="
echo " MODULE 3 — Flye LONG-READ ASSEMBLY"
echo "=========================================================="
echo "Reads: $LONG_FASTQ"
echo "Output: $OUTDIR"
echo "Read type: --$READ_TYPE"
echo "Threads: $THREADS"
[[ "$LOW_MEM" -eq 1 ]] && echo "Low-memory mode: ON (k-mer size forced to 15, see MEMORY NOTE above)"
echo ""

# --scaffold and --meta are deliberately NOT used here: --scaffold pads gaps
# with N's for draft scaffolding (premature before you've evaluated the raw
# assembly), and --meta assumes multi-organism/uneven-coverage input, which
# contradicts the single-isolate assumption this module teaches. Both are
# legitimate flags for other situations — just not this one.
FLYE_ARGS=(--"$READ_TYPE" "$LONG_FASTQ" --out-dir "$OUTDIR" --threads "$THREADS")
if [[ "$LOW_MEM" -eq 1 ]]; then
    FLYE_ARGS+=(--extra-params kmer_size=15)
fi

echo ">> Running: flye ${FLYE_ARGS[*]}"
flye "${FLYE_ARGS[@]}" 2>&1 | tee "logs/flye_$(basename "$OUTDIR").log"

ASSEMBLY="$OUTDIR/assembly.fasta"
GRAPH="$OUTDIR/assembly_graph.gfa"

if [[ ! -s "$ASSEMBLY" ]]; then
    echo "ERROR: expected $ASSEMBLY was not produced (or is empty)." >&2
    exit 1
fi
if [[ ! -s "$GRAPH" ]]; then
    echo "WARNING: expected graph file $GRAPH was not produced." >&2
fi

N_CONTIGS=$(grep -c "^>" "$ASSEMBLY")
echo ""
echo "---------------------------------------------------------"
echo " RESULT: $N_CONTIGS contig(s)/scaffold(s) written to $ASSEMBLY"
echo " Graph: $GRAPH"
echo " Flye's own summary: $OUTDIR/assembly_info.txt (coverage + circularity per contig)"
echo " Next: 4_Interpret_Assembly_Graph.sh and 5_Compute_Assembly_Statistics.sh"
echo "---------------------------------------------------------"
