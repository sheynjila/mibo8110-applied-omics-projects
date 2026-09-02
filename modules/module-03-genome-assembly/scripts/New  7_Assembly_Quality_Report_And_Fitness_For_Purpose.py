#!/usr/bin/env python3
"""
Module 3, Script 7 — Assembly Quality Report & Fitness-for-Purpose
==============================================================================
NEW IDEA on top of script 6: scripts 4-6 each computed ONE dimension of
assembly quality in isolation (topology, contiguity, completeness/
contamination). This script is where those pieces get assembled into a
single, portfolio-ready ASSEMBLY_REPORT.md -- and, critically, where the
module forces the one judgment call that no script can make for you:
whether this assembly is FIT FOR THE INTENDED PURPOSE. A report with
excellent numbers can still be unfit for a purpose that needs single-base
accuracy (e.g. calling a specific antibiotic-resistance SNP); a report with
mediocre numbers can be perfectly fit for a purpose that only needs rough
gene content (e.g. species identification). This script therefore:

  1. Re-runs (by importing, not shelling out) scripts 4-6's logic to gather
     every quantitative fact into one place.
  2. Auto-writes all of the OBJECTIVE sections (topology, statistics,
     completeness/contamination proxy results) directly into the report.
  3. Leaves the SUBJECTIVE sections -- interpretation and fitness-for-purpose
     conclusion -- as explicit blank prompts for a human to fill in. No
     script should auto-generate a "yes this assembly is good" verdict; that
     is a judgment call tied to what the assembly will be USED for, which no
     amount of automation can decide on the analyst's behalf.
  4. Explicitly distinguishes ASSEMBLY EVALUATION (the 5-dimension quality
     framework this module covers: contiguity/completeness/correctness/
     contamination/fitness-for-purpose) from GENOME ANNOTATION (identifying
     genes, operons, and function within the assembled sequence -- a
     DIFFERENT downstream task this module does NOT cover; see the README's
     "What this module does not do" section).

USAGE:
    python3 7_Assembly_Quality_Report_And_Fitness_For_Purpose.py \\
        <assembly.fasta> <graph.gfa> --reads1 <R1> [--reads2 <R2>] \\
        --assembler "SPAdes" --out ASSEMBLY_REPORT.md
"""

import sys
import argparse
import subprocess
import importlib.util
from pathlib import Path
from datetime import datetime, timezone


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def capture_script_output(script_path, args_list):
    """Run one of scripts 4/5/6 as a subprocess and capture its stdout
    verbatim, so the report contains exactly what a human would see running
    it directly -- no duplicated/out-of-sync report-formatting logic."""
    result = subprocess.run(
        [sys.executable, str(script_path)] + args_list,
        capture_output=True, text=True
    )
    if result.returncode != 0:
        return f"[ERROR running {script_path.name}: {result.stderr.strip()}]"
    return result.stdout


def main():
    parser = argparse.ArgumentParser(
        description="Assemble scripts 4-6's output into one portfolio-ready ASSEMBLY_REPORT.md."
    )
    parser.add_argument("assembly_fasta", help="Assembly FASTA (contigs.fasta or assembly.fasta)")
    parser.add_argument("graph_gfa", help="Assembly graph GFA file")
    parser.add_argument("--reads1", required=True)
    parser.add_argument("--reads2", default=None)
    parser.add_argument("--assembler", default="unspecified", help="e.g. 'SPAdes' or 'Flye'")
    parser.add_argument("--out", default="ASSEMBLY_REPORT.md")
    args = parser.parse_args()

    script_dir = Path(__file__).resolve().parent
    graph_output = capture_script_output(script_dir / "4_Interpret_Assembly_Graph.py", [args.graph_gfa])
    stats_output = capture_script_output(
        script_dir / "5_Compute_Assembly_Statistics.py",
        [args.assembly_fasta, "--label", args.assembler]
    )
    contam_args = [args.assembly_fasta, "--reads1", args.reads1]
    if args.reads2:
        contam_args += ["--reads2", args.reads2]
    contam_output = capture_script_output(
        script_dir / "6_Assess_Completeness_And_Contamination.py", contam_args
    )

    timestamp = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")

    report = f"""# Assembly Quality Report

Generated: {timestamp}
Assembler: {args.assembler}
Assembly file: `{args.assembly_fasta}`
Graph file: `{args.graph_gfa}`
Reads: `{args.reads1}`{f" + `{args.reads2}`" if args.reads2 else ""}

This report was auto-assembled by `7_Assembly_Quality_Report_And_Fitness_For_Purpose.py`
from scripts 4 (graph topology), 5 (contiguity statistics), and 6
(completeness/contamination proxy checks). Sections marked **[FILL IN]** are
intentionally left for a human analyst -- no script in this module
auto-generates a quality verdict.

---

## 1. Topology (from script 4)

```
{graph_output.strip()}
```

## 2. Contiguity statistics (from script 5)

```
{stats_output.strip()}
```

## 3. Completeness / contamination proxy checks (from script 6)

```
{contam_output.strip()}
```

## 4. Interpretation **[FILL IN]**

_Summarize, in your own words, what sections 1-3 tell you about this
assembly. Does the topology show unresolved repeats? Is the N50 high or low
relative to the expected genome size, and does that matter for your use
case? Were any contigs flagged as likely contamination, and if so, did BOTH
the GC and depth checks agree (stronger evidence) or only one (weaker,
needing a follow-up such as BLAST against NCBI nt, or a real BUSCO/CheckM
run)?_

(blank -- to be completed by the analyst)

## 5. Fitness-for-purpose conclusion **[FILL IN]**

_State explicitly what this assembly WILL and WILL NOT be used for, and
whether the evidence above supports that specific use. The same assembly can
be "fit for purpose" for one downstream use (e.g. species-level
identification) and "not fit for purpose" for another (e.g. calling a
specific point mutation) -- there is no single universal "good assembly"
verdict independent of intended use._

(blank -- to be completed by the analyst)

## 6. Scope note: assembly evaluation vs. genome annotation

This report evaluates the ASSEMBLY -- did we correctly and completely
reconstruct the genome's sequence, and how confident are we in that
reconstruction? It does **not** perform ANNOTATION -- identifying which
genes, operons, regulatory elements, or functional features exist within
that sequence. Annotation is a separate downstream task (typical tools:
Prokka, Bakta, NCBI PGAP) that assumes you already trust the assembly
underneath it. Do not skip straight to annotation-driven biological
conclusions ("this strain carries gene X") without first completing sections
4-5 above -- a gene call on a misassembled or contaminated contig is not
trustworthy regardless of how confident the annotation tool's per-gene score
looks.
"""

    out_path = Path(args.out)
    out_path.write_text(report)
    print(f"Wrote assembly quality report to {out_path}")
    print(f"({len(report.splitlines())} lines, {len(report)} bytes)")
    print("Sections 4 and 5 are intentionally blank -- fill them in before")
    print("treating this report as a finished portfolio artifact.")


if __name__ == "__main__":
    main()
