#!/usr/bin/env python3
"""
Module 3, Script 5 — Compute Assembly Statistics (and Why N50 Alone Lies)
==============================================================================
NEW IDEA on top of script 4: script 4 exposed the graph's TOPOLOGY. This
script computes the numbers most people mean when they say "assembly
quality": total length, contig count, N50/L50, GC%, longest contig. But the
one new idea this script insists on is printed as a caveat every single time
it runs, not buried in a comment: N50 is a CONTIGUITY statistic, not a
correctness, completeness, or contamination statistic. A curriculum-mandated
framing (see this module's README, "Quality framework: 5 dimensions") is that
treating N50 as a stand-in for "good assembly" is one of the most common
mistakes in genome assembly teaching -- an assembly can have an excellent N50
and still be wrong, incomplete, or contaminated. This script computes N50
alongside everything needed to build the OTHER dimensions in scripts 6-7, and
never lets it stand alone.

WHAT THIS SCRIPT COMPUTES (pure-Python FASTA parser, no dependencies):
  - Total assembly length (sum of all contig lengths)
  - Number of contigs
  - N50 / L50 (and N90 / L90 for a fuller contiguity picture)
  - Longest and shortest contig length
  - Overall GC%
  - Per-contig length + GC% table (feeds directly into script 6's
    per-contig outlier detection)

DEFINITIONS (spelled out because these are commonly misremembered):
  N50: sort contigs longest to shortest; walk down the list summing lengths;
       N50 is the length of the contig at which the running sum first
       reaches >= 50% of the total assembly length. A HIGH N50 means the
       genome is represented in a few long pieces (fewer, bigger contigs) --
       an assembly evidence for CONTIGUITY, not correctness.
  L50: the number of contigs needed to reach that 50% point (i.e. the RANK of
       the N50 contig in the sorted list). LOW L50 = more contiguous.

USAGE:
    python3 5_Compute_Assembly_Statistics.py <assembly.fasta> [--label NAME]
"""

import sys
import argparse


def read_fasta(path):
    """Parse a multi-record FASTA file into {header: sequence} preserving
    input order. Pure Python -- no biopython dependency."""
    records = {}
    header = None
    seq_chunks = []
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records[header] = "".join(seq_chunks)
                header = line[1:].split()[0]  # first whitespace-delimited token, like samtools faidx
                seq_chunks = []
            else:
                seq_chunks.append(line)
        if header is not None:
            records[header] = "".join(seq_chunks)
    return records


def gc_percent(seq):
    if not seq:
        return 0.0
    gc = sum(1 for b in seq.upper() if b in "GC")
    return 100.0 * gc / len(seq)


def compute_nx_lx(lengths, x):
    """Generic N50/L50 (and N90/L90 etc.) calculator.
    lengths: list of contig lengths (any order).
    x: percentage threshold, e.g. 50 for N50/L50.
    Returns (Nx_length, Lx_count).
    """
    sorted_lengths = sorted(lengths, reverse=True)
    total = sum(sorted_lengths)
    threshold = total * (x / 100.0)
    running = 0
    for i, length in enumerate(sorted_lengths, start=1):
        running += length
        if running >= threshold:
            return length, i
    return 0, 0  # only reached if lengths is empty


def main():
    parser = argparse.ArgumentParser(
        description="Compute contiguity/composition statistics for a genome assembly FASTA."
    )
    parser.add_argument("fasta_path", help="Path to an assembly FASTA (SPAdes contigs.fasta or Flye assembly.fasta)")
    parser.add_argument("--label", default=None, help="Optional label for this assembly (e.g. 'SPAdes' or 'Flye')")
    args = parser.parse_args()

    records = read_fasta(args.fasta_path)
    if not records:
        print(f"ERROR: no FASTA records found in {args.fasta_path}", file=sys.stderr)
        sys.exit(1)

    label = args.label or args.fasta_path
    lengths = [len(seq) for seq in records.values()]
    total_len = sum(lengths)
    n_contigs = len(lengths)
    longest = max(lengths)
    shortest = min(lengths)
    overall_gc = gc_percent("".join(records.values()))

    n50, l50 = compute_nx_lx(lengths, 50)
    n90, l90 = compute_nx_lx(lengths, 90)

    print("==========================================================")
    print(f" MODULE 3 — ASSEMBLY STATISTICS: {label}")
    print("==========================================================")
    print(f"Total assembly length: {total_len:,} bp")
    print(f"Number of contigs: {n_contigs}")
    print(f"Longest contig: {longest:,} bp")
    print(f"Shortest contig: {shortest:,} bp")
    print(f"Overall GC content: {overall_gc:.2f}%")
    print(f"N50: {n50:,} bp   (L50 = {l50} contig(s))")
    print(f"N90: {n90:,} bp   (L90 = {l90} contig(s))")
    print("")
    print("---------------------------------------------------------")
    print(" PER-CONTIG TABLE (feeds into script 6's outlier detection)")
    print("---------------------------------------------------------")
    print(f"{'Contig':<25}{'Length (bp)':>12}{'GC%':>8}")
    for name, seq in sorted(records.items(), key=lambda kv: len(kv[1]), reverse=True):
        print(f"{name:<25}{len(seq):>12,}{gc_percent(seq):>8.2f}")

    print("")
    print("---------------------------------------------------------")
    print(" MANDATORY CAVEAT — N50 IS A CONTIGUITY STATISTIC, NOT A")
    print(" QUALITY GUARANTEE")
    print("---------------------------------------------------------")
    print("A high N50 tells you the genome is represented in a small number of")
    print("long pieces. It says NOTHING about whether those pieces are:")
    print("  - CORRECT       (misassembled repeats can inflate a contig's")
    print("                    apparent length while silently joining two")
    print("                    genomic regions that don't actually belong")
    print("                    together)")
    print("  - COMPLETE      (a high-N50 assembly can still be missing genes")
    print("                    entirely if coverage was too low in places)")
    print("  - UNCONTAMINATED (a contaminant genome sequenced at high enough")
    print("                    coverage can itself assemble into one long,")
    print("                    high-N50-boosting contig -- see script 6)")
    print("A short, fragmented assembly with a low N50 can still be 100%")
    print("correct and complete; it is just harder to work with. Judge")
    print("contiguity (this script), completeness/contamination (script 6),")
    print("and fitness for purpose (script 7) as SEPARATE questions.")
    print("==========================================================")


if __name__ == "__main__":
    main()
