#!/bin/bash
#SBATCH --job-name=mod3_spades
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16      # SPAdes benefits heavily from multi-threading
#SBATCH --mem=64G               # de novo assembly is memory-intensive at real isolate scale
#SBATCH --time=12:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
###############################################################################
# Module 3, Script 2 — Run SPAdes on Short Reads, With Every Parameter Named
#
# NEW IDEA on top of script 1: script 1 only recommended a strategy. This
# script actually assembles, but the point is not "SPAdes ran" — it is that
# every non-default flag is written down with a one-line reason, so a later
# reader (including you, in six months) can tell which choices were dataset-
# specific decisions versus copy-pasted defaults.
#
# TEACHING NOTE (relationship to existing repo content): this project's
# existing modules/supplemental-case-studies/scripts/Final_master_end_to_end_amr.sh
# already calls `spades.py --isolate -t 16 -m 64` inline as one step of a much
# larger AMR-calling pipeline. That call is correct but undocumented — this
# script extracts the same core operation as its own standalone lesson, adds
# --hybrid support for long reads, and explains --isolate rather than just
# using it.
#
# WHY --isolate: SPAdes's --isolate mode is tuned for high-coverage
# (>50x is common), low-diversity bacterial isolate data — it disables some
# of SPAdes's error-correction heuristics that exist mainly to help with
# uneven-coverage or multi-cell/metagenomic input, because a clean isolate
# doesn't need them and running without --isolate on isolate data commonly
# produces MORE fragmented assemblies, not fewer. If your input is NOT a
# clean single-isolate culture (e.g. it's a metagenome or a low-coverage
# sample), --isolate is the wrong flag — use SPAdes's default mode instead.
#
# USAGE:
#   ./2_Run_SPAdes_Short_Read_Assembly.sh <R1.fastq> <R2.fastq> <output_dir> \
#       [--threads N] [--mem GB] [--hybrid long_reads.fastq]
#
# REAL WORKED EXAMPLE: SRR1770413 (E. coli K-12, Illumina MiSeq, paired-end,
# PRJNA272917, BioSample SAMN03287706) — see README.md for the full citation.
###############################################################################

set -euo pipefail
mkdir -p logs

SCRIPT_DIR="${PIPELINE_SCRIPT_DIR:-${SLURM_SUBMIT_DIR:+${SLURM_SUBMIT_DIR}/modules/module-03-genome-assembly/scripts}}"
SCRIPT_DIR="${SCRIPT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)}"
source "${SCRIPT_DIR}/load_modules.sh"

module purge
load_spades || exit 1

if [[ $# -lt 3 ]]; then
    echo "Usage: $0 <R1.fastq> <R2.fastq> <output_dir> [--threads N] [--mem GB] [--hybrid long_reads.fastq]" >&2
    exit 1
fi

R1="$1"; R2="$2"; OUTDIR="$3"
shift 3

THREADS=16
MEM_GB=64
HYBRID_LONG=""

# Named flags only from here on — no positional ambiguity between an optional
# thread count, an optional memory cap, and an optional hybrid long-read file.
while [[ $# -gt 0 ]]; do
    case "$1" in
        --threads) THREADS="$2"; shift 2 ;;
        --mem)     MEM_GB="$2";  shift 2 ;;
        --hybrid)  HYBRID_LONG="$2"; shift 2 ;;
        *) echo "ERROR: unrecognized argument: $1" >&2; exit 1 ;;
    esac
done

for f in "$R1" "$R2"; do
    if [[ ! -f "$f" ]]; then
        echo "ERROR: input read file not found: $f" >&2
        exit 1
    fi
done

mkdir -p "$OUTDIR"

echo "=========================================================="
echo " MODULE 3 — SPAdes SHORT-READ ASSEMBLY"
echo "=========================================================="
echo "R1: $R1"
echo "R2: $R2"
echo "Output: $OUTDIR"
echo "Threads: $THREADS   Memory cap: ${MEM_GB}G"
[[ -n "$HYBRID_LONG" ]] && echo "Hybrid long-read input: $HYBRID_LONG"
echo ""

SPADES_ARGS=(-1 "$R1" -2 "$R2" -o "$OUTDIR" --isolate -t "$THREADS" -m "$MEM_GB")

if [[ -n "$HYBRID_LONG" ]]; then
    if [[ ! -f "$HYBRID_LONG" ]]; then
        echo "ERROR: --hybrid long-read file not found: $HYBRID_LONG" >&2
        exit 1
    fi
    # --nanopore tells SPAdes to use these reads for hybrid gap-filling/repeat
    # resolution alongside the short reads, NOT as the primary assembly input.
    SPADES_ARGS+=(--nanopore "$HYBRID_LONG")
    echo "Mode: HYBRID (short reads primary, long reads for repeat resolution)"
else
    echo "Mode: SHORT-READ ONLY (--isolate)"
fi

echo ""
echo ">> Running: spades.py ${SPADES_ARGS[*]}"
spades.py "${SPADES_ARGS[@]}" 2>&1 | tee "logs/spades_$(basename "$OUTDIR").log"

# --- Verify the two artifacts this module's later scripts depend on both exist
#     and are non-empty, rather than assuming a zero exit code means everything
#     downstream will work.
CONTIGS="$OUTDIR/contigs.fasta"
GRAPH="$OUTDIR/assembly_graph_with_scaffolds.gfa"

if [[ ! -s "$CONTIGS" ]]; then
    echo "ERROR: expected $CONTIGS was not produced (or is empty)." >&2
    exit 1
fi
if [[ ! -s "$GRAPH" ]]; then
    echo "WARNING: expected graph file $GRAPH was not produced — script 4" >&2
    echo "(graph interpretation) will have nothing to read for this assembly." >&2
fi

N_CONTIGS=$(grep -c "^>" "$CONTIGS")
echo ""
echo "---------------------------------------------------------"
echo " RESULT: $N_CONTIGS contig(s) written to $CONTIGS"
echo " Graph (if present): $GRAPH"
echo " Next: 4_Interpret_Assembly_Graph.sh and 5_Compute_Assembly_Statistics.sh"
echo "---------------------------------------------------------"
