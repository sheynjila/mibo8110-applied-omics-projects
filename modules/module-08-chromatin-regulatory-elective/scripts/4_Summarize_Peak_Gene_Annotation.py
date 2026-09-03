#!/usr/bin/env python3
"""
Module 8, Script 4 of 5 — Summarize Peak-to-Gene Annotation
==============================================================================
NEW IDEA on top of Scripts 1-3: those scripts asked whether each sample's
peaks are trustworthy (Script 1), concentrated rather than scattered
background (Script 2), and reproducible across replicates (Script 3). None
of them answer this module's guiding question's other half: "where in the
genome" a real signal actually sits, relative to genes. master_chipseq_pipeline.sh's
own bedtools closest step already computed the nearest gene and distance for
every peak; this script is the first to actually summarize that per-sample,
turning a large per-peak BED file into the promoter-proximal-vs-distal
picture a portfolio report needs.

TEACHING NOTE — why "promoter-proximal" uses a distance threshold, not
overlap alone:
  bedtools closest -d reports 0 only when a peak genuinely overlaps a gene
  feature - but a peak 200bp upstream of a transcription start site is
  still very plausibly regulating that gene's promoter, and an overlap-only
  definition would misclassify it as "distal" alongside a peak 500kb away
  with no more evidence either way. This script uses a fixed distance
  threshold (PROMOTER_DISTANCE_BP) as the promoter-proximal/distal split -
  a common, defensible convention, not the only one; the threshold is a
  parameter precisely so it can be changed for a different regulatory
  question, not treated as a fixed fact about the genome.

Input:  WORKDIR (default /scratch/$USER/master_chipseq_pipeline_ChIP),
        containing annotation/<SRR>_peaks_annotated.bed (master_chipseq_pipeline.sh's
        bedtools closest -d output: 10 narrowPeak columns + 6 gene-feature
        columns + 1 distance column, 17 total).
Output: prints a per-sample promoter-proximal/distal breakdown;
        writes peak_annotation_summary.tsv
"""

import sys
from pathlib import Path

PROMOTER_DISTANCE_BP = 1000   # peaks within this distance of a gene feature are called promoter-proximal (see TEACHING NOTE)


def parse_annotated_bed(path):
    """Return a list of distances (int) - 0 means the peak overlaps the gene feature directly."""
    distances = []
    for line in path.read_text().splitlines():
        if not line.strip():
            continue
        cols = line.split("\t")
        if len(cols) < 17:
            continue   # malformed row - bedtools closest with no match ('.' gene columns) still has 17 cols normally; skip anything shorter
        try:
            d = int(cols[16])
        except ValueError:
            continue
        if d == -1:
            continue   # bedtools closest's "no feature found on this chromosome" sentinel
        distances.append(d)
    return distances


def main():
    if len(sys.argv) < 2:
        print("Usage: 4_Summarize_Peak_Gene_Annotation.py <WORKDIR> [output_dir]")
        sys.exit(1)
    workdir = Path(sys.argv[1])
    outdir = Path(sys.argv[2]) if len(sys.argv) >= 3 else Path(".")
    outdir.mkdir(parents=True, exist_ok=True)

    print("==================================================================")
    print(" MODULE 8 / SCRIPT 4 — Peak-to-gene annotation summary")
    print("==================================================================")
    print(f"Working directory: {workdir}\n")

    annot_dir = workdir / "annotation"
    files = sorted(annot_dir.glob("*_peaks_annotated.bed"))
    if not files:
        print(f"[FAIL] no *_peaks_annotated.bed files found in {annot_dir} — run master_chipseq_pipeline.sh (and Script 1) first")
        sys.exit(1)

    rows = []
    for f in files:
        srr = f.name.removesuffix("_peaks_annotated.bed")
        distances = parse_annotated_bed(f)
        n_total = len(distances)
        if n_total == 0:
            print(f"  {srr:20s} [WARN] no usable annotated peaks found")
            rows.append((srr, 0, 0, 0, "NO_DATA"))
            continue

        n_promoter = sum(1 for d in distances if d <= PROMOTER_DISTANCE_BP)
        n_distal = n_total - n_promoter
        pct_promoter = 100.0 * n_promoter / n_total
        median_distance = sorted(distances)[n_total // 2]

        flag = "ok"
        print(f"  {srr:20s} n_peaks={n_total:6d}  promoter-proximal(<={PROMOTER_DISTANCE_BP}bp)={n_promoter:6d} "
              f"({pct_promoter:5.1f}%)  distal={n_distal:6d}  median_distance={median_distance}bp   {flag}")
        rows.append((srr, n_total, n_promoter, n_distal, round(pct_promoter, 1)))

    print()
    print(f"Promoter-proximal threshold: {PROMOTER_DISTANCE_BP}bp (see this script's header for why this is a")
    print("parameter, not a fixed rule) — distal peaks are not automatically less real, they are simply")
    print("further from an annotated gene feature and need different evidence (e.g. enhancer marks, Hi-C/HiChIP")
    print("contact data) before being attributed to a specific target gene.")

    out_tsv = outdir / "peak_annotation_summary.tsv"
    with open(out_tsv, "w") as fh:
        fh.write("sample\tn_peaks\tn_promoter_proximal\tn_distal\tpct_promoter_proximal\n")
        for row in rows:
            fh.write("\t".join(str(x) for x in row) + "\n")
    print(f"\nWrote {out_tsv}")
    print("==================================================================")


if __name__ == "__main__":
    main()
