#!/usr/bin/env python3
"""
Module 7, Script 5 of 5 — Assemble the Single-Cell Report
==============================================================================
NEW IDEA on top of Scripts 1-4: each prior script computed ONE fact about
this dataset in isolation (QC filtering, doublet rate, cluster sizes, marker
genes). This script assembles those facts into the curriculum's required
portfolio artifact -- "a QC'd, clustered single-cell dataset with marker-gene
identification and an interpretation of what the clusters represent
biologically" -- the same auto-fill-the-objective-sections/leave-
interpretation-blank pattern every other module's report script in this
repository uses (DE_REPORT.md, ASSEMBLY_REPORT.md, VARIANT_REPORT.md,
PHYLO_REPORT.md, MICROBIOME_REPORT.md). The one thing this script refuses to
automate is cell-type identity itself -- Script 4 already stops at marker
GENES, deliberately, and this script does not go further than Script 4 did.

Input:  --qc-tsv, --doublet-tsv, --cluster-sizes-tsv, --top-markers-tsv
        (all optional; missing inputs are noted explicitly, not silently
        skipped -- see Module 6 Script 5 for the same pattern).
Output: SINGLECELL_REPORT.md
"""

import argparse
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path


def read_tsv_as_md_table(path, missing_note):
    if path is None or not Path(path).exists():
        return missing_note
    lines = Path(path).read_text().splitlines()
    rows = [line.split("\t") for line in lines if line.strip() != ""]
    if len(rows) < 2:
        return missing_note
    header, body = rows[0], rows[1:]
    md = ["| " + " | ".join(header) + " |", "|" + "|".join(["---"] * len(header)) + "|"]
    for r in body:
        md.append("| " + " | ".join(r) + " |")
    return "\n".join(md)


def read_tsv_as_dict(path):
    if path is None or not Path(path).exists():
        return {}
    lines = Path(path).read_text().splitlines()
    out = {}
    for line in lines[1:]:
        if not line.strip():
            continue
        cols = line.split("\t")
        if len(cols) >= 2:
            out[cols[0]] = cols[1]
    return out


def markers_by_cluster_md(path, missing_note):
    if path is None or not Path(path).exists():
        return missing_note
    lines = Path(path).read_text().splitlines()
    if len(lines) < 2:
        return missing_note
    header = lines[0].split("\t")
    idx = {name: i for i, name in enumerate(header)}
    by_cluster = defaultdict(list)
    for line in lines[1:]:
        if not line.strip():
            continue
        cols = line.split("\t")
        by_cluster[cols[idx["cluster"]]].append(cols[idx["gene"]])
    parts = []
    for cluster in sorted(by_cluster, key=lambda c: (len(c), c)):
        genes = ", ".join(by_cluster[cluster])
        parts.append(f"- **Cluster {cluster}:** {genes}")
    return "\n".join(parts)


def main():
    parser = argparse.ArgumentParser(description="Assemble Scripts 1-4's output into SINGLECELL_REPORT.md.")
    parser.add_argument("--qc-tsv", default=None, help="qc_summary.tsv (Script 1)")
    parser.add_argument("--doublet-tsv", default=None, help="doublet_summary.tsv (Script 2)")
    parser.add_argument("--cluster-sizes-tsv", default=None, help="cluster_sizes.tsv (Script 3)")
    parser.add_argument("--top-markers-tsv", default=None, help="top_markers_per_cluster.tsv (Script 4)")
    parser.add_argument("--out", default="SINGLECELL_REPORT.md")
    args = parser.parse_args()

    print("==================================================================")
    print(" MODULE 7 / SCRIPT 5 — Assembling SINGLECELL_REPORT.md")
    print("==================================================================")

    timestamp = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
    qc = read_tsv_as_dict(Path(args.qc_tsv)) if args.qc_tsv else {}
    doublet = read_tsv_as_dict(Path(args.doublet_tsv)) if args.doublet_tsv else {}
    cluster_md = read_tsv_as_md_table(
        args.cluster_sizes_tsv,
        "_Not available — run Script 3 (`3_Cluster_Cells.R`)._")
    markers_md = markers_by_cluster_md(
        args.top_markers_tsv,
        "_Not available — run Script 4 (`4_Identify_Cluster_Markers.R`)._")

    qc_line = (f"{qc.get('cells_after_qc', '?')}/{qc.get('cells_before_qc', '?')} cells retained "
               f"({qc.get('pct_retained', '?')}%)") if qc else "_Not available — run Script 1 (`1_Load_And_QC_Filter_Cells.R`)._"
    doublet_line = (f"{doublet.get('cells_flagged_doublet', '?')}/{doublet.get('cells_in', '?')} cells "
                    f"flagged as doublets ({doublet.get('pct_doublet', '?')}%), "
                    f"{doublet.get('cells_out', '?')} singlets retained") if doublet else \
                   "_Not available — run Script 2 (`2_Detect_And_Remove_Doublets.R`)._"

    report = f"""# Single-Cell Report

Generated: {timestamp}

Auto-assembled by `5_Assemble_SingleCell_Report.py` from Scripts 1-4's own
output. Section 5 is intentionally left blank — no script in this module
assigns cell-type identity on its own; Script 4 stops at marker GENES,
deliberately.

---

## 1. QC filtering (Script 1)

{qc_line}

Thresholds applied: `nFeature_RNA >= 200`, `nCount_RNA >= 500`,
`percent.mt <= 15%` (Script 1's own defaults — confirm these suited this
specific sample via `qc_metrics_violin.pdf` before trusting them unchanged
on different tissue).

## 2. Doublet detection (Script 2)

{doublet_line}

Detected computationally (scDblFinder), not via this dataset's own real
species-mixing ground truth — see Script 2's header for why, and for what
extending this check to that ground truth would require.

## 3. Clustering (Script 3)

{cluster_md}

Clustering ran in PCA space (20 components) at resolution=0.8 — a first-pass
default, not a claim that this is the correct number of clusters for this
tissue (see Script 3's header).

## 4. Cluster markers (Script 4)

Top marker genes per cluster (`p_val_adj < 0.05`, sorted by adjusted
p-value):

{markers_md}

## 5. Interpretation **[FILL IN]**

**Before writing anything below, this module's core caution applies
unconditionally:**

> A cluster is a statement about transcriptional similarity, not a
> confirmed cell type. Marker genes above are evidence for a cell-type
> assignment, not an assignment themselves — cross-reference them against
> known marker genes for the tissue actually sampled before naming any
> cluster. A cluster can also reflect a technical axis (batch, cell-cycle
> phase, ambient RNA contamination) rather than a real, distinct cell
> population — nothing in Scripts 1-4 rules this out on its own, and this
> dataset's own species-mixing design (Script 2's header) is a reminder that
> what looks like a clean, well-separated cluster is not automatically free
> of pipeline-specific artifacts.

_Given Sections 1-4 above: what cell type(s) does each cluster's marker list
suggest, and how confident is that call given the cluster's size (Section 3)
and whether it had any significant markers at all (Script 4's console
output)? Does the doublet rate (Section 2) change how much you trust any
particular cluster's boundary? What would you need beyond this pipeline's
output — a reference atlas for automated label transfer, additional marker
genes from the literature, or (for this specific dataset) a combined
human+mouse alignment to check scDblFinder's calls against real ground
truth — before treating any cluster's identity as established rather than a
working hypothesis?_

(blank — to be completed by the analyst)
"""

    out_path = Path(args.out)
    out_path.write_text(report)
    print(f"Wrote {out_path} ({len(report.splitlines())} lines, {len(report)} bytes)")
    print("Section 5 is intentionally blank, and its own caveat text is NOT optional filler —")
    print("fill in an interpretation only alongside that caveat, not instead of it.")
    print("==================================================================")


if __name__ == "__main__":
    main()
