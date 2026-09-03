#!/usr/bin/env python3
"""
Module 5, Script 2 of 6 — Assess Alignment Quality
==============================================================================
NEW IDEA on top of Script 1: Script 1 produced an alignment; it did not
check whether that alignment is actually trustworthy to build a tree from.
This is the same "look at the data before trusting the model" habit Module 2
Script 4 applied to a count matrix (PCA, outlier detection) and Module 4
Script 2 applied to a BAM (mapping rate, coverage), applied here one level
further downstream: a SNP alignment that is mostly 'N' for one sample will
still produce A tree — IQ-TREE does not refuse to run on it — but that
sample's branch length and placement will be driven by missingness, not by
real relatedness, and nothing about the tree image itself will say so.

TEACHING NOTE — why per-sample AND per-site missingness both matter:
  High missingness in one SAMPLE (row) means that sample's placement in the
  tree is unreliable, specifically. High missingness at one SITE (column)
  means that position's contribution to every sample's relationships is
  unreliable, generally. A sample can look "clean" in the per-site view and
  still be untrustworthy if its own missingness is concentrated non-randomly
  (e.g. an entire low-coverage genomic region for that isolate) -- this
  script reports both, deliberately, rather than collapsing them into a
  single alignment-wide number.

Input:  core_snp_alignment.fasta (Script 1)
Output: prints a PASS/WARN per-sample and per-site missingness report;
        writes alignment_qc_summary.tsv
"""

import sys
from pathlib import Path

MAX_SAMPLE_MISSING_PCT = 10.0   # flag a sample with > 10% N across all sites
MAX_SITE_MISSING_PCT = 20.0     # flag a site with > 20% N across all samples


def read_fasta(path):
    sequences = {}
    name = None
    chunks = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    sequences[name] = "".join(chunks)
                name = line[1:].strip()
                chunks = []
            else:
                chunks.append(line.strip())
        if name is not None:
            sequences[name] = "".join(chunks)
    return sequences


def main():
    if len(sys.argv) < 2:
        print("Usage: 2_Assess_Alignment_Quality.py <core_snp_alignment.fasta> [output_dir]")
        sys.exit(1)

    fasta_path = Path(sys.argv[1])
    outdir = Path(sys.argv[2]) if len(sys.argv) >= 3 else Path(".")
    outdir.mkdir(parents=True, exist_ok=True)

    print("==================================================================")
    print(" MODULE 5 / SCRIPT 2 — Alignment quality assessment")
    print("==================================================================")
    print(f"Alignment: {fasta_path}\n")

    if not fasta_path.exists() or fasta_path.stat().st_size == 0:
        print(f"[FAIL] alignment not found or empty: {fasta_path}")
        sys.exit(1)

    seqs = read_fasta(fasta_path)
    if len(seqs) == 0:
        print("[FAIL] no sequences parsed from FASTA — is Script 1's output well-formed?")
        sys.exit(1)

    lengths = {name: len(seq) for name, seq in seqs.items()}
    if len(set(lengths.values())) != 1:
        print("[FAIL] sequences are not all the same length — this is not a valid alignment:")
        for name, ln in lengths.items():
            print(f"       {name}: {ln} bp")
        sys.exit(1)
    n_sites = next(iter(lengths.values()))
    n_samples = len(seqs)
    print(f"[PASS] {n_samples} sequences, all {n_sites} bp — structurally a valid alignment\n")

    # --- Per-sample missingness --------------------------------------------
    print("Per-sample missingness:")
    flagged_samples = []
    sample_rows = []
    for name, seq in sorted(seqs.items()):
        n_count = seq.upper().count("N")
        pct = 100.0 * n_count / n_sites if n_sites else 0.0
        flag = "FLAGGED" if pct > MAX_SAMPLE_MISSING_PCT else "ok"
        if flag == "FLAGGED":
            flagged_samples.append(name)
        print(f"  {name:20s} {n_count:6d}/{n_sites} N ({pct:5.1f}%)   {flag}")
        sample_rows.append((name, n_count, round(pct, 1), flag))

    if flagged_samples:
        print(f"\n[WARN] {len(flagged_samples)} sample(s) exceed {MAX_SAMPLE_MISSING_PCT}% missingness: "
              f"{', '.join(flagged_samples)}")
        print("       Their placement/branch length in Script 3's tree should be treated with caution.")
    else:
        print(f"\n[PASS] no sample exceeds {MAX_SAMPLE_MISSING_PCT}% missingness")

    # --- Per-site missingness ------------------------------------------------
    names_in_order = sorted(seqs.keys())
    site_missing_pct = []
    for i in range(n_sites):
        n_count = sum(1 for name in names_in_order if seqs[name][i].upper() == "N")
        site_missing_pct.append(100.0 * n_count / n_samples if n_samples else 0.0)

    n_flagged_sites = sum(1 for p in site_missing_pct if p > MAX_SITE_MISSING_PCT)
    print(f"\nPer-site missingness: {n_flagged_sites}/{n_sites} site(s) exceed {MAX_SITE_MISSING_PCT}% missing "
          f"across samples")
    if n_flagged_sites > 0:
        print(f"[WARN] consider whether these sites should be excluded before tree inference "
              f"(not automated here — a site-exclusion decision changes what the alignment IS, "
              f"and belongs to the analyst, the same reasoning Module 2 used for the gene-filter rule).")
    else:
        print("[PASS] no site exceeds the per-site missingness threshold")

    out_tsv = outdir / "alignment_qc_summary.tsv"
    with open(out_tsv, "w") as fh:
        fh.write("sample\tn_count\tpct_missing\tflag\n")
        for row in sample_rows:
            fh.write("\t".join(str(x) for x in row) + "\n")
    print(f"\nWrote {out_tsv}")
    print("==================================================================")


if __name__ == "__main__":
    main()
