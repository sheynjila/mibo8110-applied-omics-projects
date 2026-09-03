#!/usr/bin/env python3
"""
Module 8, Script 5 of 5 — Assemble the Chromatin Report
==============================================================================
NEW IDEA on top of Scripts 1-4: each prior script computed ONE fact about
this cohort in isolation (input validity/control status, signal
concentration, replicate reproducibility, gene proximity). This script
assembles those facts into the curriculum's required portfolio artifact --
"a peak/methylation call set + QC summary + an interpretation of the
regulatory signal, with explicit discussion of assay-specific background
handling" -- the same auto-fill-the-objective-sections/leave-interpretation-
blank pattern every other module's report script in this repository uses.
The one thing this script refuses to automate is the actual regulatory
interpretation -- and it prints this module's core background-handling
caution as unconditional boilerplate, not just a prompt the analyst could
skip past, because "distinguish a real regulatory signal from assay
background" is this module's own guiding question, not an optional caveat.

Input:  --control-status (free text, from Script 1's console output --
        Script 1 does not write a machine-readable file, see its own
        header), --frip-tsv, --reproducibility-tsv, --annotation-tsv (all
        optional except --control-status; missing inputs are noted
        explicitly, not silently skipped).
Output: CHROMATIN_REPORT.md
"""

import argparse
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


def main():
    parser = argparse.ArgumentParser(description="Assemble Scripts 1-4's output into CHROMATIN_REPORT.md.")
    parser.add_argument("--assay-type", default="ChIP", help="'ChIP' or 'ATAC' (must match the actual pipeline run)")
    parser.add_argument("--control-status", required=True,
                         help="Free-text summary of Script 1's input-control finding, e.g. "
                              "'3/3 samples called WITHOUT an input control (none deposited for this study)'")
    parser.add_argument("--frip-tsv", default=None, help="frip_summary.tsv (Script 2)")
    parser.add_argument("--reproducibility-tsv", default=None, help="replicate_reproducibility.tsv (Script 3)")
    parser.add_argument("--annotation-tsv", default=None, help="peak_annotation_summary.tsv (Script 4)")
    parser.add_argument("--out", default="CHROMATIN_REPORT.md")
    args = parser.parse_args()

    print("==================================================================")
    print(" MODULE 8 / SCRIPT 5 — Assembling CHROMATIN_REPORT.md")
    print("==================================================================")

    timestamp = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
    frip_md = read_tsv_as_md_table(
        args.frip_tsv,
        "_Not available — run Script 2 (`2_Compute_FRiP_Score.sh`)._")
    repro_md = read_tsv_as_md_table(
        args.reproducibility_tsv,
        "_Not available — run Script 3 (`3_Assess_Replicate_Reproducibility.sh`)._")
    annot_md = read_tsv_as_md_table(
        args.annotation_tsv,
        "_Not available — run Script 4 (`4_Summarize_Peak_Gene_Annotation.py`)._")

    report = f"""# Chromatin Report

Generated: {timestamp}
Assay type: {args.assay_type}

Auto-assembled by `5_Assemble_Chromatin_Report.py` from Scripts 1-4's own
output. Section 5 is intentionally left blank — no script in this module
generates a regulatory interpretation on its own.

---

## 1. Input validation and background handling (Script 1)

**Input-control status:** {args.control_status}

{"An input/IgG control models a ChIP-seq experiment's own assay background directly; peaks called without one rest entirely on the treatment track's own signal, which is a real limitation to carry into Section 5, not a footnote." if args.assay_type == "ChIP" else "ATAC-seq has no input-control concept — background is instead handled by the ATAC-specific MACS2 parameters (--nomodel --shift -75 --extsize 150) used by this pipeline's Tier A run, correcting for Tn5 insertion bias rather than modeling a ChIP-seq-style fragment-length distribution."}

## 2. Signal concentration — FRiP (Script 2)

Fraction of Reads in Peaks per sample — ENCODE's own minimum-viability
guidance is >1%, with >5% expected for a well-behaved transcription-factor
ChIP-seq experiment (broad marks and ATAC-seq legitimately run lower; see
Script 2's header):

{frip_md}

## 3. Replicate reproducibility (Script 3)

Pairwise Jaccard similarity (bases shared / bases covered by either peak
set) between every pair of peak sets in this cohort:

{repro_md}

A pair flagged `LOW_REPRODUCIBILITY` is only a real concern if that pair is
actually a replicate comparison — confirm from your own sample metadata
before treating a low score between two genuinely different conditions as a
problem.

## 4. Peak-to-gene annotation (Script 4)

Promoter-proximal vs. distal breakdown per sample (threshold and rationale
in Script 4's own header):

{annot_md}

## 5. Regulatory interpretation **[FILL IN]**

**Before writing anything below, this module's core caution applies
unconditionally:**

> A peak with a strong score, a high FRiP contribution, and a reproducible
> position across replicates is still not, on its own, proof of a
> functional regulatory event — it is evidence of a real, reproducible
> binding or accessibility signal at that position. Whether that signal
> actually regulates the nearest gene (Section 4), a more distant gene, or
> nothing at all requires evidence this pipeline does not have access to:
> functional validation (reporter assays, CRISPR perturbation), chromatin
> conformation data (Hi-C/HiChIP) linking the peak to a specific promoter,
> or expression data showing the putative target gene actually responds.
> Do not write a regulatory-function claim below without naming what
> additional evidence would be needed, regardless of how clean the numbers
> above look.

_Given Sections 1-4 above: which peaks have both a strong FRiP contribution
and reproducible replicate support (Sections 2-3) — and does the
input-control status (Section 1) change how much weight you put on any of
them? Which promoter-proximal peaks (Section 4) plausibly regulate the gene
they're nearest to, and which distal peaks would need additional evidence
(e.g. this study's own HiChIP data, if available) before attributing them to
any specific target?_

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
