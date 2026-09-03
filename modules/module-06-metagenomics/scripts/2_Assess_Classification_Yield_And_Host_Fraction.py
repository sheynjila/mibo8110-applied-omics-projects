#!/usr/bin/env python3
"""
Module 6, Script 2 of 5 — Assess Classification Yield and Host Fraction
==============================================================================
NEW IDEA on top of Script 1: Script 1 confirmed the Kraken2/Bracken files
exist and are well-formed; it did not look at what fraction of each sample's
reads those files actually represent. A Kraken2 report can be structurally
perfect and still mean "we could confidently place almost none of this
sample" -- and Bracken re-estimates SPECIES-LEVEL abundance only among the
reads Kraken2 already classified, so a cohort abundance table built on a
30%-classified sample and a 90%-classified sample are not comparable on
their face, even though both produce a row in ALL_SAMPLES_MICROBIOME.tsv
that looks equally confident. This is the module's guiding-question split
made concrete: "unclassified" is either real signal (an organism truly
absent from Kraken2's reference database -- see Script 4) or an artifact of
sequencing depth/quality, and this script's job is only to report the
number, not to decide which explanation applies to any one sample.

TEACHING NOTE — why host fraction is reported here, not folded into the
classification-yield number:
  A host-depleted sample's classification % is computed AFTER host removal
  (Kraken2 never saw the host reads Bowtie2 already pulled out), so a sample
  that was 60% host DNA and 95%-classified-of-the-remainder is answering a
  different question than a sample with no host depletion that is also
  95%-classified. Reporting both numbers side by side, instead of one
  blended figure, is what lets you tell those two cases apart.

Input:  WORKDIR (default /scratch/$USER/master_microbiome_pipeline),
        containing kraken_output/<SRR>_report.txt and, optionally,
        qc/<SRR>_host_depletion.log (Script 1's inputs).
Output: prints a per-sample classification-yield + host-fraction report;
        writes classification_yield_summary.tsv
"""

import re
import sys
from pathlib import Path

LOW_CLASSIFIED_PCT = 50.0   # flag a sample where less than half of reads classified


def parse_kraken_report(path):
    """Return (pct_unclassified, pct_classified) from a Kraken2 --report file."""
    lines = path.read_text().splitlines()
    pct_unclassified = 0.0
    for line in lines:
        cols = line.split("\t")
        if len(cols) != 6:
            continue
        pct, _clade, _direct, rank_code, _taxid, name = cols
        if rank_code == "U" or name.strip().lower() == "unclassified":
            pct_unclassified = float(pct)
            break
    return pct_unclassified, 100.0 - pct_unclassified


def parse_host_depletion_log(path):
    """Return the Bowtie2 'overall alignment rate' (host %) if present, else None."""
    if not path.exists():
        return None
    text = path.read_text()
    m = re.search(r"([\d.]+)%\s*overall alignment rate", text)
    return float(m.group(1)) if m else None


def main():
    if len(sys.argv) < 2:
        print("Usage: 2_Assess_Classification_Yield_And_Host_Fraction.py <WORKDIR> [output_dir]")
        sys.exit(1)
    workdir = Path(sys.argv[1])
    outdir = Path(sys.argv[2]) if len(sys.argv) >= 3 else Path(".")
    outdir.mkdir(parents=True, exist_ok=True)

    kraken_dir = workdir / "kraken_output"
    qc_dir = workdir / "qc"

    print("==================================================================")
    print(" MODULE 6 / SCRIPT 2 — Classification yield and host fraction")
    print("==================================================================")
    print(f"Working directory: {workdir}\n")

    reports = sorted(kraken_dir.glob("*_report.txt"))
    if not reports:
        print(f"[FAIL] no *_report.txt files found in {kraken_dir} — run Script 1 first")
        sys.exit(1)

    rows = []
    flagged = []
    for report in reports:
        srr = report.name.removesuffix("_report.txt")
        pct_unclassified, pct_classified = parse_kraken_report(report)
        host_pct = parse_host_depletion_log(qc_dir / f"{srr}_host_depletion.log")
        host_str = f"{host_pct:.1f}" if host_pct is not None else "NA"

        flag = "LOW_YIELD" if pct_classified < LOW_CLASSIFIED_PCT else "ok"
        if flag == "LOW_YIELD":
            flagged.append(srr)

        print(f"  {srr:20s} classified={pct_classified:5.1f}%  unclassified={pct_unclassified:5.1f}%  "
              f"host_removed={host_str:>5s}%   {flag}")
        rows.append((srr, round(pct_classified, 1), round(pct_unclassified, 1), host_str, flag))

    print()
    if flagged:
        print(f"[WARN] {len(flagged)} sample(s) classified less than {LOW_CLASSIFIED_PCT:.0f}% of reads: "
              f"{', '.join(flagged)}")
        print("       Low classification yield does not mean 'few organisms present' — it means Kraken2's")
        print("       database could not confidently place most of this sample's reads. Treat this sample's")
        print("       Bracken abundance fractions as computed over a small, possibly biased subset of its")
        print("       actual community, not the whole community.")
    else:
        print(f"[PASS] every sample classified at least {LOW_CLASSIFIED_PCT:.0f}% of reads")

    out_tsv = outdir / "classification_yield_summary.tsv"
    with open(out_tsv, "w") as fh:
        fh.write("sample\tpct_classified\tpct_unclassified\tpct_host_removed\tflag\n")
        for row in rows:
            fh.write("\t".join(str(x) for x in row) + "\n")
    print(f"\nWrote {out_tsv}")
    print("==================================================================")


if __name__ == "__main__":
    main()
