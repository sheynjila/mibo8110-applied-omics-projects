#!/usr/bin/env python3
"""
Module 6, Script 5 of 5 — Assemble the Microbiome Report
==============================================================================
NEW IDEA on top of Scripts 1-4: each prior script computed ONE fact about
this cohort in isolation (input validity, classification yield/host
fraction, diversity/dominance, database-detection caveats). This script
assembles those facts into the curriculum's required portfolio artifact --
"a microbiome profile report: method chosen and why, QC/contamination
handling, abundance results, and an interpretation that acknowledges
database and detection-limit caveats" -- the same auto-fill-the-objective-
sections/leave-interpretation-blank pattern every other module's report
script in this repository uses (DE_REPORT.md, ASSEMBLY_REPORT.md,
VARIANT_REPORT.md, PHYLO_REPORT.md). The one thing this script refuses to
automate is the actual interpretation, and it prints this module's core
database/detection-limit caveat as unconditional boilerplate, not just a
prompt the analyst could skip past, because Script 4 already established
that no script in this module can resolve that ambiguity on its own.

Input:  --validation-log, --yield-tsv, --diversity-tsv, --caveats-tsv
        (all optional; missing inputs are noted explicitly, not silently
        skipped -- see Module 4 Script 7 / Module 5 Script 6 for the same
        pattern).
Output: MICROBIOME_REPORT.md
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
    parser = argparse.ArgumentParser(description="Assemble Scripts 1-4's output into MICROBIOME_REPORT.md.")
    parser.add_argument("--method", default="shotgun (Kraken2/Bracken)",
                         help="Which route this report covers, e.g. 'shotgun (Kraken2/Bracken)' or '16S amplicon (DADA2)'.")
    parser.add_argument("--yield-tsv", default=None, help="classification_yield_summary.tsv (Script 2)")
    parser.add_argument("--diversity-tsv", default=None, help="cohort_diversity_summary.tsv (Script 3)")
    parser.add_argument("--caveats-tsv", default=None, help="database_detection_caveats.tsv (Script 4)")
    parser.add_argument("--out", default="MICROBIOME_REPORT.md")
    args = parser.parse_args()

    print("==================================================================")
    print(" MODULE 6 / SCRIPT 5 — Assembling MICROBIOME_REPORT.md")
    print("==================================================================")

    timestamp = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
    yield_md = read_tsv_as_md_table(
        args.yield_tsv,
        "_Not available — run Script 2 (`2_Assess_Classification_Yield_And_Host_Fraction.py`)._")
    diversity_md = read_tsv_as_md_table(
        args.diversity_tsv,
        "_Not available — run Script 3 (`3_Cohort_Diversity_And_Dominance_Summary.py`)._")
    caveats_md = read_tsv_as_md_table(
        args.caveats_tsv,
        "_Not available — run Script 4 (`4_Detection_Limit_And_Database_Caveats.py`). "
        "A taxon present in every sample is a legitimate outcome with nothing to report here, not an error._")

    report = f"""# Microbiome Profile Report

Generated: {timestamp}
Method: {args.method}

Auto-assembled by `5_Assemble_Microbiome_Report.py` from Scripts 1-4's own
output. Section 5 is intentionally left blank — no script in this module
generates a biological interpretation on its own.

---

## 1. Method chosen and why

This report covers the **{args.method}** route. This module ships two
worked pipelines — shotgun (`master_microbiome_pipeline.sh`, Kraken2/Bracken)
and amplicon (`master_amplicon_pipeline.sh` + `amplicon_dada2_analysis.R`,
DADA2) — because there is no single right choice for every mixed-community
sample, only a right choice for a given question and budget. See Script 1's
header for the full tradeoff (functional-profiling breadth and no primer
bias vs. lower per-sample cost and finer within-marker resolution). Confirm
this matches the actual analytical goal for this sample set before treating
anything below as the community's true composition rather than this
method's view of it.

## 2. QC and contamination handling (Scripts 1-2)

Script 1 confirmed every sample's Kraken2 report and Bracken table are
structurally valid and belong to each other before anything below was
trusted. Classification yield and host-read fraction (if `HOST_REFERENCE`
was set):

{yield_md}

A sample flagged `LOW_YIELD` above classified less than half its reads —
treat that sample's abundance fractions in Section 3 as computed over a
small, possibly biased subset of its actual community.

## 3. Abundance results (Script 3)

Per-sample diversity (Shannon H') and single-taxon dominance check:

{diversity_md}

A sample flagged `DOMINANT` above has one taxon exceeding the dominance
threshold — not automatically wrong, but worth confirming that taxon's
identity matches the sample's expected biology.

## 4. Database and detection-limit caveats (Script 4)

Taxa present in some cohort samples but absent from others, cross-referenced
against those samples' own classification yield:

{caveats_md}

**Unconditional caveat — applies regardless of what the table above shows:**

> A missing row for a taxon in this cohort's Bracken output is consistent
> with three different realities this pipeline cannot distinguish on its
> own: the organism was genuinely absent, it was present below Bracken's
> minimum-read re-estimation threshold, or it was present but is not
> represented in the Kraken2 reference database used for classification at
> all. Treating "not in the table" as "confirmed absent" overstates what
> this pipeline actually measured.

## 5. Interpretation **[FILL IN]**

**Before writing anything below, this module's core caution applies
unconditionally (restated from Section 4 — it is not optional filler):**

> Do not write a claim of true absence, or a claim that one sample's
> community is definitively less diverse/rich than another's, without first
> checking that both samples being compared had comparable classification
> yield (Section 2) and that neither report is being pushed past what Bracken
> and Kraken2's reference database can actually support.

_Given Sections 1-4 above: what does this cohort's composition and diversity
actually say about the biological question that motivated collecting these
samples? Which flagged samples (LOW_YIELD, DOMINANT, or WEAK_ABSENCE_EVIDENCE
taxa) change how much weight you put on which comparisons? What would you
need beyond this pipeline's output — deeper sequencing, a different or larger
reference database, functional profiling, or a validated amplicon
cross-check — before this report's numbers could support a stronger claim
than the one you are about to write?_

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
