#!/usr/bin/env python3
"""
Module 6, Script 4 of 5 — Detection Limit and Database Caveats
==============================================================================
NEW IDEA on top of Scripts 1-3: those scripts characterized what IS in each
sample's Bracken table. This script characterizes what is NOT — specifically,
what it does and does not mean when a taxon that shows up in one cohort
sample has no row at all in another. Kraken2/Bracken's abundance table is not
a presence/absence assay: a missing row can mean "genuinely absent," "present
but below Bracken's minimum-read re-estimation threshold," or "present but
not represented in Kraken2's reference database at all" (a k-mer database
cannot classify a genome it has never seen, however good the read is). This
module's guiding question — "how much of that signal is sequencing/database
artifact versus real biology" — is exactly this ambiguity, and no script
upstream of this one resolves it, because it cannot be resolved from this
pipeline's own output alone. This script's job is to make the ambiguity
explicit and quantified per taxon, not to pretend it can be settled here.

TEACHING NOTE — why a low-classification-yield sample (Script 2) weakens an
"absent" call more than a high-yield one:
  A taxon's row being missing from a sample that classified 95% of its reads
  is much stronger (though still not proof-positive) evidence of true
  absence than the same missing row in a sample that only classified 40% of
  its reads — the low-yield sample simply explored less of its own read
  space before Bracken ever got a chance to see that taxon. This script
  cross-references each "missing in sample X" case against Script 2's own
  classification-yield number for sample X, rather than treating every
  absence in the cohort table as equally informative.

Input:  ALL_SAMPLES_MICROBIOME.tsv (Script 3's input) and
        classification_yield_summary.tsv (Script 2's output).
Output: prints a per-taxon presence/absence caveat report;
        writes database_detection_caveats.tsv
"""

import sys
from collections import defaultdict
from pathlib import Path

LOW_CLASSIFIED_PCT = 50.0   # must match Script 2's threshold


def read_cohort_taxa(path):
    """Return {sample: {taxon}} and the sorted list of all samples/taxa seen."""
    by_sample = defaultdict(set)
    lines = path.read_text().splitlines()
    header = lines[0].split("\t")
    idx = {name: i for i, name in enumerate(header)}
    for line in lines[1:]:
        if not line.strip():
            continue
        cols = line.split("\t")
        by_sample[cols[idx["Sample"]]].add(cols[idx["Taxon"]])
    return by_sample


def read_yield_summary(path):
    """Return {sample: pct_classified} from Script 2's output, or {} if absent."""
    if not path.exists():
        return {}
    lines = path.read_text().splitlines()
    header = lines[0].split("\t")
    idx = {name: i for i, name in enumerate(header)}
    out = {}
    for line in lines[1:]:
        if not line.strip():
            continue
        cols = line.split("\t")
        out[cols[idx["sample"]]] = float(cols[idx["pct_classified"]])
    return out


def main():
    if len(sys.argv) < 2:
        print("Usage: 4_Detection_Limit_And_Database_Caveats.py <ALL_SAMPLES_MICROBIOME.tsv> "
              "[classification_yield_summary.tsv] [output_dir]")
        sys.exit(1)

    cohort_path = Path(sys.argv[1])
    yield_path = Path(sys.argv[2]) if len(sys.argv) >= 3 else Path("classification_yield_summary.tsv")
    outdir = Path(sys.argv[3]) if len(sys.argv) >= 4 else Path(".")
    outdir.mkdir(parents=True, exist_ok=True)

    print("==================================================================")
    print(" MODULE 6 / SCRIPT 4 — Detection limit and database caveats")
    print("==================================================================")
    print(f"Cohort table: {cohort_path}")
    print(f"Yield summary: {yield_path}\n")

    if not cohort_path.exists() or cohort_path.stat().st_size == 0:
        print(f"[FAIL] cohort table not found or empty: {cohort_path}")
        sys.exit(1)

    by_sample = read_cohort_taxa(cohort_path)
    yields = read_yield_summary(yield_path)
    if not yields:
        print("[WARN] classification_yield_summary.tsv not found — run Script 2 first for the")
        print("       yield cross-reference below. Continuing with presence/absence only.\n")

    all_samples = sorted(by_sample)
    all_taxa = sorted({t for taxa in by_sample.values() for t in taxa})

    print(f"{len(all_taxa)} distinct taxon/taxa across {len(all_samples)} sample(s).\n")

    rows = []
    n_ambiguous = 0
    for taxon in all_taxa:
        present_in = [s for s in all_samples if taxon in by_sample[s]]
        absent_in = [s for s in all_samples if taxon not in by_sample[s]]
        if not absent_in:
            continue   # present in every sample — nothing ambiguous to report

        low_yield_absences = [s for s in absent_in if yields.get(s, 100.0) < LOW_CLASSIFIED_PCT]
        weak = len(low_yield_absences) > 0
        if weak:
            n_ambiguous += 1

        note = (f"absent from {len(absent_in)}/{len(all_samples)} sample(s); "
                f"{len(low_yield_absences)} of those had <{LOW_CLASSIFIED_PCT:.0f}% classification yield "
                f"(Script 2) — weak evidence of true absence there")
        print(f"  {taxon}")
        print(f"    present in: {', '.join(present_in) if present_in else '(none)'}")
        print(f"    {note}")
        rows.append((taxon, len(present_in), len(absent_in), len(low_yield_absences),
                     "WEAK_ABSENCE_EVIDENCE" if weak else "ok"))

    print()
    if n_ambiguous > 0:
        print(f"[WARN] {n_ambiguous} taxon/taxa have at least one 'absent' call in a low-yield sample —")
        print("       do not report those as confirmed absences without checking Script 2's yield number.")
    else:
        print("[PASS] no taxon's absence calls are undermined by a low-classification-yield sample")

    print("\nUnconditional caveat — applies to every taxon in this report, not only the flagged ones:")
    print("  A missing row for a taxon in Bracken's output is consistent with three different realities")
    print("  this pipeline cannot distinguish on its own: (1) the organism was genuinely absent from that")
    print("  sample, (2) it was present below Bracken's minimum-read re-estimation threshold, or (3) it was")
    print("  present in the actual sample but is not represented in the Kraken2 database used for")
    print("  classification at all — a database cannot classify a genome it has never seen, however good")
    print("  the underlying read is. None of these three are distinguishable from ALL_SAMPLES_MICROBIOME.tsv")
    print("  alone; treating 'not in the table' as 'confirmed absent' overstates what this pipeline measured.")

    out_tsv = outdir / "database_detection_caveats.tsv"
    with open(out_tsv, "w") as fh:
        fh.write("taxon\tn_present\tn_absent\tn_low_yield_absences\tflag\n")
        for row in rows:
            fh.write("\t".join(str(x) for x in row) + "\n")
    print(f"\nWrote {out_tsv}")
    print("==================================================================")


if __name__ == "__main__":
    main()
