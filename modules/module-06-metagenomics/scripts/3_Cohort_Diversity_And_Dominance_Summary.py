#!/usr/bin/env python3
"""
Module 6, Script 3 of 5 — Cohort Diversity and Dominance Summary
==============================================================================
NEW IDEA on top of Scripts 1-2: those scripts asked "can we trust this
sample's classification at all." This script asks the first actual biology
question about a trusted sample: how diverse is its community, and is any
one taxon suspiciously dominant? Bracken's per-sample TSV and the pipeline's
own ALL_SAMPLES_MICROBIOME.tsv report relative abundance per taxon; neither
computes a single per-sample diversity number, and neither flags dominance —
both are left as something a human would eyeball in a spreadsheet, which
does not scale past a handful of samples and does not happen consistently
across a whole cohort report.

TEACHING NOTE — why Shannon diversity, and why dominance is checked
separately rather than inferred from a low diversity score:
  Shannon entropy (H' = -sum(p_i * ln(p_i))) rewards both richness (how many
  taxa) and evenness (how balanced their abundances are) in one number, which
  is exactly why it should NOT be over-interpreted alone: a sample with 40
  taxa where one holds 95% of reads and a sample with 3 nearly-even taxa can
  land at a similar H' despite being biologically very different pictures —
  one is a single dominant organism with a long tail of low-abundance noise
  (contamination, a spillover taxon, or PCR/database artifact are all
  plausible), the other is a genuinely low-richness but balanced community.
  This script reports the DOMINANCE flag (single taxon exceeding
  DOMINANCE_THRESHOLD) as its own independent check for exactly that reason
  — a diversity number alone would not have caught the first case.

Input:  bracken_output/ALL_SAMPLES_MICROBIOME.tsv (columns: Sample, Taxon,
        Estimated_Reads, Fraction_Total) — written by
        master_microbiome_pipeline.sh's own PHASE 3 aggregation step.
Output: prints a per-sample diversity + dominance report;
        writes cohort_diversity_summary.tsv
"""

import math
import sys
from collections import defaultdict
from pathlib import Path

DOMINANCE_THRESHOLD = 50.0   # flag a sample where one taxon exceeds this % of reads


def read_cohort_table(path):
    """Return {sample: [(taxon, fraction_total), ...]} from ALL_SAMPLES_MICROBIOME.tsv."""
    by_sample = defaultdict(list)
    lines = path.read_text().splitlines()
    if not lines:
        return by_sample
    header = lines[0].split("\t")
    idx = {name: i for i, name in enumerate(header)}
    for col in ("Sample", "Taxon", "Fraction_Total"):
        if col not in idx:
            raise ValueError(f"expected column '{col}' not found in {path} header: {header}")
    for line in lines[1:]:
        if not line.strip():
            continue
        cols = line.split("\t")
        sample = cols[idx["Sample"]]
        taxon = cols[idx["Taxon"]]
        fraction = float(cols[idx["Fraction_Total"]])
        by_sample[sample].append((taxon, fraction))
    return by_sample


def shannon_diversity(fractions):
    """Shannon entropy H' over a set of relative abundances that sum to ~1."""
    total = sum(fractions)
    if total <= 0:
        return 0.0
    h = 0.0
    for f in fractions:
        p = f / total
        if p > 0:
            h -= p * math.log(p)
    return h


def main():
    if len(sys.argv) < 2:
        print("Usage: 3_Cohort_Diversity_And_Dominance_Summary.py <ALL_SAMPLES_MICROBIOME.tsv> [output_dir]")
        sys.exit(1)

    cohort_path = Path(sys.argv[1])
    outdir = Path(sys.argv[2]) if len(sys.argv) >= 3 else Path(".")
    outdir.mkdir(parents=True, exist_ok=True)

    print("==================================================================")
    print(" MODULE 6 / SCRIPT 3 — Cohort diversity and dominance summary")
    print("==================================================================")
    print(f"Cohort table: {cohort_path}\n")

    if not cohort_path.exists() or cohort_path.stat().st_size == 0:
        print(f"[FAIL] cohort table not found or empty: {cohort_path}")
        sys.exit(1)

    by_sample = read_cohort_table(cohort_path)
    if not by_sample:
        print("[FAIL] no sample rows parsed — is this really master_microbiome_pipeline.sh's PHASE 3 output?")
        sys.exit(1)

    rows = []
    flagged = []
    for sample in sorted(by_sample):
        taxa = by_sample[sample]
        fractions = [f for _name, f in taxa]
        h = shannon_diversity(fractions)
        top_taxon, top_fraction = max(taxa, key=lambda t: t[1])
        top_pct = 100.0 * top_fraction / sum(fractions) if sum(fractions) > 0 else 0.0

        flag = "DOMINANT" if top_pct > DOMINANCE_THRESHOLD else "ok"
        if flag == "DOMINANT":
            flagged.append((sample, top_taxon, top_pct))

        print(f"  {sample:20s} n_taxa={len(taxa):3d}  Shannon_H={h:5.2f}  "
              f"top_taxon='{top_taxon}' ({top_pct:5.1f}%)   {flag}")
        rows.append((sample, len(taxa), round(h, 2), top_taxon, round(top_pct, 1), flag))

    print()
    if flagged:
        print(f"[WARN] {len(flagged)} sample(s) have a single taxon exceeding {DOMINANCE_THRESHOLD:.0f}% of reads:")
        for sample, taxon, pct in flagged:
            print(f"       {sample}: '{taxon}' at {pct:.1f}%")
        print("       A dominant taxon is not automatically wrong — a real infection, bloom, or a genuinely")
        print("       low-complexity sample can all look like this. It IS a reason to check that taxon's")
        print("       identity against the sample's expected biology before reporting the diversity number")
        print("       alone as 'this community's diversity'.")
    else:
        print(f"[PASS] no sample has a single taxon exceeding {DOMINANCE_THRESHOLD:.0f}% of reads")

    out_tsv = outdir / "cohort_diversity_summary.tsv"
    with open(out_tsv, "w") as fh:
        fh.write("sample\tn_taxa\tshannon_h\ttop_taxon\ttop_taxon_pct\tflag\n")
        for row in rows:
            fh.write("\t".join(str(x) for x in row) + "\n")
    print(f"\nWrote {out_tsv}")
    print("==================================================================")


if __name__ == "__main__":
    main()
