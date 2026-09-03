#!/usr/bin/env python3
"""
Module 5, Script 6 of 6 — Assemble the Phylogenomics Report
==============================================================================
NEW IDEA on top of Scripts 1-5: each prior script computed ONE fact about
this cohort in isolation (alignment quality, tree, branch support, tree-vs-
SNP-cluster consistency). This script assembles those facts into the
curriculum's required portfolio artifact -- "a phylogenetic tree +
branch-support summary + a written outbreak interpretation that respects the
relatedness-vs-transmission distinction" -- the same auto-fill-the-objective-
sections/leave-interpretation-blank pattern every other module's report
script in this repository uses (DE_REPORT.md, ASSEMBLY_REPORT.md,
VARIANT_REPORT.md). The one thing this script refuses to automate is the
actual outbreak interpretation -- and it prints this module's core caution
as unconditional boilerplate, not just as a prompt the analyst could skip
past, because the curriculum names this specific distinction as the module's
central learning outcome, not an optional caveat.

Input:  --tree, --support-tsv, --alignment-qc-tsv, --crossref-tsv, --iqtree-report
        (all optional except --tree; missing inputs are noted explicitly,
        not silently skipped -- see Module 4 Script 7 for the same pattern).
Output: PHYLO_REPORT.md
"""

import argparse
import re
from datetime import datetime, timezone
from pathlib import Path


def read_tsv_as_md_table(path, missing_note):
    if path is None or not Path(path).exists():
        return missing_note
    lines = Path(path).read_text().splitlines()
    if len(lines) < 1:
        return missing_note
    rows = [line.split("\t") for line in lines if line.strip() != ""]
    if len(rows) < 2:
        return missing_note
    header, body = rows[0], rows[1:]
    md = ["| " + " | ".join(header) + " |", "|" + "|".join(["---"] * len(header)) + "|"]
    for r in body:
        md.append("| " + " | ".join(r) + " |")
    return "\n".join(md)


def extract_best_model(iqtree_report_path):
    if iqtree_report_path is None or not Path(iqtree_report_path).exists():
        return "_Not available — pass --iqtree-report cohort_tree.iqtree (Script 3's full report)._"
    text = Path(iqtree_report_path).read_text()
    m = re.search(r"Best-fit model according to.*?:\s*(\S+)", text)
    if m:
        return f"`{m.group(1)}` (selected by ModelFinder, ascertainment-bias corrected — see Script 3)"
    return "_Model line not found in the given .iqtree report — check the file is Script 3's output._"


def main():
    parser = argparse.ArgumentParser(description="Assemble Scripts 1-5's output into PHYLO_REPORT.md.")
    parser.add_argument("--tree", required=True, help="cohort_tree.treefile (Script 3)")
    parser.add_argument("--iqtree-report", default=None, help="cohort_tree.iqtree (Script 3)")
    parser.add_argument("--alignment-qc-tsv", default=None, help="alignment_qc_summary.tsv (Script 2)")
    parser.add_argument("--support-tsv", default=None, help="branch_support_summary.tsv (Script 4)")
    parser.add_argument("--crossref-tsv", default=None, help="tree_snp_cluster_crossref.tsv (Script 5)")
    parser.add_argument("--out", default="PHYLO_REPORT.md")
    args = parser.parse_args()

    print("==================================================================")
    print(" MODULE 5 / SCRIPT 6 — Assembling PHYLO_REPORT.md")
    print("==================================================================")

    if not Path(args.tree).exists():
        print(f"[FAIL] tree file not found: {args.tree}")
        raise SystemExit(1)

    timestamp = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
    model_line = extract_best_model(args.iqtree_report)
    qc_md = read_tsv_as_md_table(
        args.alignment_qc_tsv,
        "_Not available — run Script 2 (`2_Assess_Alignment_Quality.py`)._")
    support_md = read_tsv_as_md_table(
        args.support_tsv,
        "_Not available — run Script 4 (`4_Interpret_Branch_Support.py`)._")
    crossref_md = read_tsv_as_md_table(
        args.crossref_tsv,
        "_Not available — run Script 5 (`5_Cross_Reference_Tree_And_SNP_Clusters.py`). "
        "Note: no `outbreak_clusters.tsv` from Module 4 is itself a legitimate outcome, not an error._")

    report = f"""# Phylogenomics Report

Generated: {timestamp}
Tree file: `{args.tree}`

Auto-assembled by `6_Assemble_Phylogenomics_Report.py` from Scripts 1-5's own
output. Section 5 is intentionally left blank — no script in this module
generates an outbreak/transmission conclusion on its own.

---

## 1. Alignment (Scripts 1-2)

- Alignment: core-SNP, ascertainment-biased by construction (invariant sites
  excluded) — see Script 1's header note.

**Alignment quality (Script 2):**

{qc_md}

## 2. Tree inference (Script 3)

- Model selected: {model_line}
- Bootstrap: 1000 ultrafast bootstrap (UFBoot) replicates.

## 3. Branch support (Script 4)

UFBoot support >= 95% is read as well-supported (IQ-TREE's own recommended
cutoff for this bootstrap method — NOT the classic Felsenstein 70% cutoff;
see Script 4's header note for why the two are not interchangeable).

{support_md}

## 4. Tree vs. SNP-distance cluster consistency (Script 5)

Cross-references Module 4's `outbreak_clusters.tsv` (pairwise SNP-distance
threshold calls) against this tree's actual topology and support.

{crossref_md}

## 5. Outbreak interpretation **[FILL IN]**

**Before writing anything below, this module's core caution applies
unconditionally:**

> A close phylogenetic distance — even an exclusive, well-supported
> (UFBoot >= 95%) sister-clade relationship — is **necessary but not
> sufficient** evidence of a direct transmission event. Shared ancestry is
> not the same claim as "isolate A infected isolate B" or "both isolates
> came from the same source at the same time." A real epidemiological
> conclusion requires evidence this module does not have access to:
> sampling dates and locations, exposure history, and independent outbreak
> investigation data. Do not write a transmission claim below without that
> evidence, regardless of how strong the numbers above look.

_Given Sections 1-4 above: which isolates, if any, show both SNP-distance
proximity AND exclusive, well-supported clade membership? Where do the two
signals (SNP distance, tree topology+support) disagree, and what does that
disagreement tell you about which signal to trust less for that pair? What
would you need — beyond what this pipeline produced — before treating any
grouping here as an actual transmission cluster rather than a relatedness
hypothesis?_

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
