#!/usr/bin/env python3
"""
Module 5, Script 4 of 6 — Interpret Branch Support
==============================================================================
NEW IDEA on top of Script 3: Script 3 produced a tree TOPOLOGY. This module's
own learning outcomes explicitly separate "infer a tree" from "read
branch-support values (not just topology) as a measure of confidence" -- a
topology with no support values attached is a hypothesis about relationships,
not evidence for them. This script is that second, distinct step: it reads
ONLY the UFBoot support values IQ-TREE attached to internal branches, and
classifies each one, rather than treating the tree image as self-evidently
trustworthy everywhere it draws a line.

TEACHING NOTE — why 95%, not the classic 70%, and why this distinction
matters enough to state explicitly:
  Classic (Felsenstein) bootstrap's traditional "well supported" cutoff is
  ~70%. IQ-TREE's ULTRAFAST bootstrap (UFBoot, what Script 3 actually ran
  with `-B 1000`) is a different, less conservative estimator by
  construction -- its own documentation (Hoang et al. 2018, the UFBoot2
  paper) recommends reading UFBoot support against a ~95% cutoff instead.
  Applying the classic 70% threshold to a UFBoot value overstates
  confidence in every branch between 70-94% support. This script uses 95%,
  and labels every number it reports as "UFBoot support", not the generic
  and easily-conflated term "bootstrap support".

Pure-Python, dependency-light Newick parsing (no ete3/dendropy/Biopython) --
same "no third-party packages, deliberately" choice Module 3 made for its
GFA graph parsing, for the same reason: this module should not gain a
dependency chain beyond the tree-inference tool itself.

Input:  cohort_tree.treefile (Script 3's IQ-TREE output; Newick with UFBoot
        support labels at internal nodes)
Output: prints a per-branch UFBoot support report; writes
        branch_support_summary.tsv
"""

import sys
from pathlib import Path

UFBOOT_WELL_SUPPORTED = 95.0


class Node:
    __slots__ = ("children", "name", "support", "length")

    def __init__(self):
        self.children = []
        self.name = None
        self.support = None
        self.length = None

    def is_leaf(self):
        return len(self.children) == 0


def parse_newick(text):
    s = text.strip()
    if s.endswith(";"):
        s = s[:-1]
    pos = [0]

    def parse_label():
        start = pos[0]
        while pos[0] < len(s) and s[pos[0]] not in ",():;":
            pos[0] += 1
        return s[start:pos[0]]

    def parse_number():
        label = parse_label()
        try:
            return float(label)
        except ValueError:
            return None

    def parse_node():
        node = Node()
        if pos[0] < len(s) and s[pos[0]] == "(":
            pos[0] += 1
            node.children.append(parse_node())
            while pos[0] < len(s) and s[pos[0]] == ",":
                pos[0] += 1
                node.children.append(parse_node())
            if pos[0] < len(s) and s[pos[0]] == ")":
                pos[0] += 1
            label = parse_label()
            node.support = label if label else None
        else:
            node.name = parse_label()
        if pos[0] < len(s) and s[pos[0]] == ":":
            pos[0] += 1
            node.length = parse_number()
        return node

    return parse_node()


def get_leaves(node):
    if node.is_leaf():
        return [node.name]
    leaves = []
    for c in node.children:
        leaves.extend(get_leaves(c))
    return leaves


def collect_internal_nodes(node, acc):
    if not node.is_leaf():
        acc.append(node)
        for c in node.children:
            collect_internal_nodes(c, acc)
    return acc


def main():
    if len(sys.argv) < 2:
        print("Usage: 4_Interpret_Branch_Support.py <cohort_tree.treefile> [output_dir]")
        sys.exit(1)

    tree_path = Path(sys.argv[1])
    outdir = Path(sys.argv[2]) if len(sys.argv) >= 3 else Path(".")
    outdir.mkdir(parents=True, exist_ok=True)

    print("==================================================================")
    print(" MODULE 5 / SCRIPT 4 — Interpreting UFBoot branch support")
    print("==================================================================")
    print(f"Tree file: {tree_path}\n")

    if not tree_path.exists() or tree_path.stat().st_size == 0:
        print(f"[FAIL] tree file not found or empty: {tree_path}")
        sys.exit(1)

    text = tree_path.read_text().strip()
    root = parse_newick(text)

    internal_nodes = collect_internal_nodes(root, [])
    if len(internal_nodes) == 0:
        print("[WARN] no internal branches found — a tree this small (<=3 taxa, one unrooted topology) "
              "has no branch support to interpret. Add more isolates for a meaningful support analysis.")
        sys.exit(0)

    rows = []
    n_well_supported = 0
    n_weak = 0
    n_unlabeled = 0
    print(f"{'clade (leaf set)':60s} {'UFBoot':>8s}  status")
    for i, node in enumerate(internal_nodes):
        leaves = sorted(get_leaves(node))
        clade_label = ",".join(leaves)
        if node.support in (None, ""):
            status = "unlabeled (root or no support value written here)"
            support_val = ""
            n_unlabeled += 1
        else:
            try:
                support_val = float(node.support)
            except ValueError:
                support_val = None
            if support_val is None:
                status = f"unparsed label '{node.support}'"
                n_unlabeled += 1
            elif support_val >= UFBOOT_WELL_SUPPORTED:
                status = f"WELL SUPPORTED (UFBoot >= {UFBOOT_WELL_SUPPORTED:.0f}%)"
                n_well_supported += 1
            else:
                status = f"WEAK (UFBoot < {UFBOOT_WELL_SUPPORTED:.0f}%) — do not over-interpret this grouping"
                n_weak += 1

        display = clade_label if len(clade_label) <= 58 else clade_label[:55] + "..."
        support_display = f"{support_val:.0f}" if isinstance(support_val, float) else str(support_val)
        print(f"{display:60s} {support_display:>8s}  {status}")
        rows.append((clade_label, support_val, status))

    print(f"\n{n_well_supported} well-supported, {n_weak} weak, {n_unlabeled} unlabeled internal branch(es) "
          f"out of {len(internal_nodes)} total.")
    if n_weak > 0:
        print(f"[WARN] {n_weak} branch(es) fall below the UFBoot {UFBOOT_WELL_SUPPORTED:.0f}% threshold — "
              "the corresponding groupings should not be treated as confirmed relationships in Script 6's "
              "written interpretation.")

    out_tsv = outdir / "branch_support_summary.tsv"
    with open(out_tsv, "w") as fh:
        fh.write("clade_leaves\tufboot_support\tstatus\n")
        for clade_label, support_val, status in rows:
            sv = f"{support_val:.1f}" if isinstance(support_val, float) else ""
            fh.write(f"{clade_label}\t{sv}\t{status}\n")
    print(f"\nWrote {out_tsv}")
    print("==================================================================")


if __name__ == "__main__":
    main()
