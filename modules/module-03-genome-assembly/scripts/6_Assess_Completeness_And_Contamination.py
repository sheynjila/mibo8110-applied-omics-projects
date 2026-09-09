#!/usr/bin/env python3
"""
Module 3, Script 6 — Lightweight Completeness/Contamination Proxy Checks
==============================================================================
NEW IDEA on top of script 5: script 5 computed per-contig length and GC%, and
explicitly said that number alone does not tell you about completeness or
contamination. This script is where those two quality dimensions actually get
addressed -- with cheap, dependency-light PROXY checks that run in seconds,
not the real tools.

WHAT THIS SCRIPT ACTUALLY DOES (two independent proxy signals):
  1. GC-CONTENT OUTLIER DETECTION: computes each contig's GC% (reusing the
     same math as script 5) and flags any contig whose GC% deviates from the
     assembly-wide mean by more than a configurable threshold (default: 10
     percentage points). Rationale: within one bacterial genome, GC content
     is remarkably uniform across the chromosome (a real biological property
     called "compositional homogeneity"); a contig with wildly different GC%
     is a classic, well-established signature of contamination from a
     different organism, or of a plasmid/mobile element with different
     compositional signature that deserves a closer look either way.
  2. READ-DEPTH OUTLIER DETECTION: maps the ORIGINAL reads back to the
     assembled contigs with minimap2 + samtools, computes mean per-contig
     read depth, and flags contigs whose depth deviates sharply from the
     assembly-wide median. Rationale: if 95% of your reads came from your
     target organism and 5% from a contaminant, the contaminant's contig(s)
     will show a depth roughly proportional to that minority fraction --
     depth is an independent line of evidence from GC%, and agreement
     between the two is stronger evidence than either alone.

CRITICAL, CURRICULUM-MANDATED CAVEAT (mirrors this project's Module 2
apeglm-fallback pattern: document the exact real tool and command, then
explain the simplified stand-in actually used here):
    The REAL tools for completeness and contamination assessment are BUSCO
    (Benchmarking Universal Single-Copy Orthologs) and CheckM, which compare
    your assembly against curated marker-gene sets for your organism's
    taxonomic lineage. Those require multi-gigabyte lineage databases that
    are impractical to bundle with a teaching module and were NOT run here.
    Real commands you would run on an HPC node with those databases installed:

        busco -i assembly.fasta -o busco_out -m genome -l bacteria_odb10 --cpu 8
        checkm lineage_wf -x fasta assembly_dir/ checkm_out/ --threads 8

    This script's GC/depth outlier checks are a fast, dependency-light PROXY
    for contamination only -- they say nothing about gene-level completeness
    (a genome can have uniform GC and depth and still be missing 20% of its
    genes from low coverage in specific regions). Use BUSCO for that question.

USAGE:
    python3 6_Assess_Completeness_And_Contamination.py <assembly.fasta> \\
        --reads1 <R1.fastq> [--reads2 <R2.fastq>] [--gc-threshold 10] [--depth-fold 3]

    (--reads2 omitted => treats --reads1 as single-end/long reads for minimap2 -ax map-ont)
"""

import sys
import argparse
import subprocess
import shutil
import tempfile
import os


def read_fasta(path):
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
                header = line[1:].split()[0]
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


def length_weighted_median_depth(depth_by_contig, records):
    """A plain statistics.median() over PER-CONTIG depth values is dominated
    by how many contigs exist, not by how much of the genome they represent.
    A real genome fragmented into many small low-depth pieces plus one large
    high-depth contig would incorrectly put the 'median' near the small
    fragments' depth (since they outnumber the one big contig), which then
    flags the correct, dominant contig as the 'outlier' and misses the real
    outliers. Weighting by contig length (so the value is effectively 'the
    depth of the base-pair sitting in the middle of the assembly, if you
    laid all contigs end to end sorted by depth') fixes this: a few short,
    low-depth contigs cannot out-vote one contig that makes up most of the
    assembly's total length, which is what we actually want as the baseline
    for 'typical, likely-correct' coverage.
    """
    pairs = sorted(((depth_by_contig.get(name, 0.0), len(seq)) for name, seq in records.items()))
    total_len = sum(length for _depth, length in pairs)
    if total_len == 0:
        return 0.0
    half = total_len / 2.0
    running = 0
    for depth, length in pairs:
        running += length
        if running >= half:
            return depth
    return pairs[-1][0]


def check_tool(name):
    if shutil.which(name) is None:
        print(f"ERROR: required tool '{name}' not found on PATH.", file=sys.stderr)
        sys.exit(1)


def compute_depth_per_contig(assembly_path, reads1, reads2, workdir):
    """Map reads back to the assembly and compute mean per-contig depth.
    Uses minimap2 (works for both short paired-end and long reads with the
    right preset) + samtools sort/index/depth. Returns {contig_name: mean_depth}.
    """
    check_tool("minimap2")
    check_tool("samtools")

    bam_path = os.path.join(workdir, "mapped.sorted.bam")
    preset = "sr" if reads2 else "map-ont"

    if reads2:
        minimap_cmd = ["minimap2", "-ax", preset, "-t", "2", assembly_path, reads1, reads2]
    else:
        minimap_cmd = ["minimap2", "-ax", preset, "-t", "2", assembly_path, reads1]

    sort_cmd = ["samtools", "sort", "-o", bam_path, "-"]

    p1 = subprocess.Popen(minimap_cmd, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    p2 = subprocess.Popen(sort_cmd, stdin=p1.stdout, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    p1.stdout.close()
    p2.communicate()
    p1.wait()

    subprocess.run(["samtools", "index", bam_path], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    depth_result = subprocess.run(
        ["samtools", "depth", "-a", bam_path],
        check=True, capture_output=True, text=True
    )

    per_contig_sum = {}
    per_contig_count = {}
    for line in depth_result.stdout.splitlines():
        contig, _pos, depth = line.split("\t")
        depth = int(depth)
        per_contig_sum[contig] = per_contig_sum.get(contig, 0) + depth
        per_contig_count[contig] = per_contig_count.get(contig, 0) + 1

    return {
        c: (per_contig_sum[c] / per_contig_count[c] if per_contig_count[c] else 0.0)
        for c in per_contig_sum
    }


def main():
    parser = argparse.ArgumentParser(
        description="Proxy checks for contamination via per-contig GC-content and read-depth outliers."
    )
    parser.add_argument("assembly_path", help="Assembly FASTA (SPAdes contigs.fasta or Flye assembly.fasta)")
    parser.add_argument("--reads1", required=True, help="R1 (paired-end) or the only read file (long-read)")
    parser.add_argument("--reads2", default=None, help="R2 for paired-end short reads (omit for long reads)")
    parser.add_argument("--gc-threshold", type=float, default=10.0,
                         help="Flag a contig if its GC%% differs from the assembly mean by more than this many percentage points (default: 10)")
    parser.add_argument("--depth-fold", type=float, default=3.0,
                         help="Flag a contig if its mean depth differs from the assembly median by more than this fold-change (default: 3x)")
    args = parser.parse_args()

    records = read_fasta(args.assembly_path)
    if not records:
        print(f"ERROR: no FASTA records found in {args.assembly_path}", file=sys.stderr)
        sys.exit(1)

    print("==========================================================")
    print(" MODULE 3 — COMPLETENESS/CONTAMINATION PROXY CHECKS")
    print("==========================================================")
    print(f"Assembly: {args.assembly_path} ({len(records)} contig(s))")
    print("")
    print("REMINDER: these are fast PROXY checks, not BUSCO/CheckM. See the")
    print("script's header comment for the real commands to run when lineage")
    print("databases are available.")
    print("")

    # --- GC outlier check ---
    # IMPORTANT: use the LENGTH-WEIGHTED (whole-assembly) GC%% as the
    # reference point, not a plain average of each contig's GC%%. With only
    # a handful of contigs, an unweighted average is easily dragged toward a
    # short outlier contig -- e.g. one 49,400 bp contig at 50%% GC and one
    # 2,400 bp contig at 83%% GC unweighted-average to ~67%%, which is far
    # from BOTH real values and would flag the dominant, almost certainly
    # correct, majority contig as if it were the outlier. Weighting by
    # length (equivalent to computing GC%% over the whole concatenated
    # assembly, same math as script 5's "Overall GC content") correctly
    # reflects that the genome IS mostly the long contig, so short
    # compositionally different contigs are properly seen as the anomaly.
    gc_by_contig = {name: gc_percent(seq) for name, seq in records.items()}
    mean_gc = gc_percent("".join(records.values()))

    print("---------------------------------------------------------")
    print(f" GC-CONTENT OUTLIER CHECK (threshold: +/-{args.gc_threshold} percentage points from mean)")
    print("---------------------------------------------------------")
    print(f"Assembly-wide (length-weighted) GC%%: {mean_gc:.2f}%%".replace("%%", "%"))
    gc_flagged = set()
    for name, gc in sorted(gc_by_contig.items(), key=lambda kv: -len(records[kv[0]])):
        delta = abs(gc - mean_gc)
        flag = delta > args.gc_threshold
        marker = "  <-- FLAGGED (GC outlier)" if flag else ""
        print(f"  {name:<35} GC={gc:6.2f}%  (delta from mean: {delta:5.2f}pp){marker}")
        if flag:
            gc_flagged.add(name)

    # --- Depth outlier check ---
    depth_flagged = set()
    depth_by_contig = {}
    if shutil.which("minimap2") and shutil.which("samtools"):
        with tempfile.TemporaryDirectory() as workdir:
            try:
                depth_by_contig = compute_depth_per_contig(args.assembly_path, args.reads1, args.reads2, workdir)
            except subprocess.CalledProcessError as e:
                print(f"WARNING: read-mapping-based depth check failed to run ({e}); skipping depth check.", file=sys.stderr)

        if depth_by_contig:
            median_depth = length_weighted_median_depth(depth_by_contig, records)
            print("")
            print("---------------------------------------------------------")
            print(f" READ-DEPTH OUTLIER CHECK (threshold: >{args.depth_fold}x deviation from median)")
            print("---------------------------------------------------------")
            print(f"Assembly length-weighted median depth: {median_depth:.1f}x")
            for name in sorted(records, key=lambda n: -len(records[n])):
                depth = depth_by_contig.get(name, 0.0)
                if median_depth > 0:
                    fold = max(depth, 1e-9) / median_depth if depth < median_depth else depth / median_depth
                else:
                    fold = 1.0
                is_outlier = (depth < median_depth / args.depth_fold) or (depth > median_depth * args.depth_fold)
                marker = "  <-- FLAGGED (depth outlier)" if is_outlier else ""
                print(f"  {name:<35} depth={depth:7.1f}x  (fold vs median: {fold:5.2f}x){marker}")
                if is_outlier:
                    depth_flagged.add(name)
    else:
        print("")
        print("NOTE: minimap2/samtools not found on PATH -- skipping the read-depth")
        print("outlier check. GC-content check above still ran.")

    # --- Combined verdict ---
    print("")
    print("---------------------------------------------------------")
    print(" COMBINED VERDICT")
    print("---------------------------------------------------------")
    both_flagged = gc_flagged & depth_flagged
    either_flagged = gc_flagged | depth_flagged
    if both_flagged:
        print(f"Contig(s) flagged by BOTH GC and depth checks (strongest evidence of")
        print(f"contamination or a distinct mobile element): {', '.join(sorted(both_flagged))}")
    if either_flagged - both_flagged:
        print(f"Contig(s) flagged by only ONE check (weaker, worth a manual look --")
        print(f"e.g. BLAST the contig against NCBI nt, or run BUSCO/CheckM):")
        print(f"  {', '.join(sorted(either_flagged - both_flagged))}")
    if not either_flagged:
        print("No contigs flagged by either proxy check. This is REASSURING but")
        print("NOT PROOF of a clean, complete assembly -- these proxies can miss")
        print("low-level contamination or completeness gaps that only BUSCO/CheckM")
        print("would catch. Report this as 'no proxy evidence of contamination'")
        print("in your portfolio artifact, not as 'confirmed contamination-free'.")

    print("==========================================================")


if __name__ == "__main__":
    main()
