#!/usr/bin/env python3
"""
Module 5, Script 5 of 6 — Cross-Reference the Tree Against Module 4's
SNP-Distance Clusters
==============================================================================
NEW IDEA on top of Script 4: Module 4's variant_analysis.R already flags
isolate PAIRS within a SNP-distance threshold as a probable
outbreak/transmission cluster (outbreak_clusters.tsv), independent of any
tree. A SNP-distance pair and a phylogenetic clade are related but NOT the
same claim: a pairwise distance says nothing about whether OTHER isolates in
the cohort fall between that pair in shared ancestry, and a distance
threshold has no branch-support concept at all. This script checks, for each
Module 4 cluster pair, whether the tree actually places them as an exclusive
(sister) clade, and if so, at what UFBoot support -- turning "these two are
SNP-close" into "these two are SNP-close AND phylogenetically exclusive AND
that grouping has X% support", or explicitly flagging when any one of those
three does not hold.

TEACHING NOTE — this is the module's core caution, made structural, not just
stated:
  A close phylogenetic distance (or a tight SNP distance) is NECESSARY but
  NOT SUFFICIENT evidence of a direct transmission event. This script cannot
  and does not decide whether two isolates transmitted between hosts/sources
  -- it only checks internal CONSISTENCY between two computational signals
  (SNP distance and tree topology+support) that a real epidemiological
  investigation would need to look considerably beyond (sampling context,
  timing, exposure data) before drawing that conclusion. Disagreement
  between the two signals here is itself informative -- it means at least
  one of them should not be trusted alone.

Pure-Python Newick parsing (no ete3/dendropy/Biopython), matching Script 4.

Input:  cohort_tree.treefile (Script 3), outbreak_clusters.tsv (Module 4's
        variant_analysis.R — pass its path; if it does not exist, this
        script reports that plainly and exits without error, since "no
        SNP-close pairs were found" is a legitimate outcome, not a failure).
Output: prints a per-pair consistency report; writes
        tree_snp_cluster_crossref.tsv
"""

import csv
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


def find_mrca(root, a, b):
    result = [None]

    def visit(node):
        if node.is_leaf():
            return {node.name} & {a, b}
        found = set()
        for c in node.children:
            found |= visit(c)
        if a in found and b in found and result[0] is None:
            result[0] = node
        return found

    visit(root)
    return result[0]


def main():
    if len(sys.argv) < 2:
        print("Usage: 5_Cross_Reference_Tree_And_SNP_Clusters.py <cohort_tree.treefile> "
              "[outbreak_clusters.tsv] [output_dir]")
        sys.exit(1)

    tree_path = Path(sys.argv[1])
    clusters_path = Path(sys.argv[2]) if len(sys.argv) >= 3 else Path("outbreak_clusters.tsv")
    outdir = Path(sys.argv[3]) if len(sys.argv) >= 4 else Path(".")
    outdir.mkdir(parents=True, exist_ok=True)

    print("==================================================================")
    print(" MODULE 5 / SCRIPT 5 — Cross-referencing tree vs. SNP-distance clusters")
    print("==================================================================")
    print(f"Tree file:      {tree_path}")
    print(f"Cluster file:   {clusters_path}\n")

    if not tree_path.exists() or tree_path.stat().st_size == 0:
        print(f"[FAIL] tree file not found or empty: {tree_path}")
        sys.exit(1)
    root = parse_newick(tree_path.read_text().strip())
    all_leaves = set(get_leaves(root))

    if not clusters_path.exists():
        print(f"[INFO] {clusters_path} not found — Module 4's variant_analysis.R found no isolate pair "
              "within its SNP threshold (or was not yet run). Nothing to cross-reference; this is a "
              "legitimate outcome, not a script failure.")
        sys.exit(0)

    rows_out = []
    with open(clusters_path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        pairs = list(reader)

    if not pairs:
        print(f"[INFO] {clusters_path} exists but has no rows — no SNP-close pairs to cross-reference.")
        sys.exit(0)

    print(f"{'pair':30s} {'SNP dist':>8s}  {'MRCA leaves':>11s}  {'UFBoot':>8s}  verdict")
    for row in pairs:
        a, b = row["Sample_A"], row["Sample_B"]
        snp_dist = row.get("SNP_Distance", "?")
        if a not in all_leaves or b not in all_leaves:
            pair_label = f"{a}/{b}"
            print(f"{pair_label:30s} {snp_dist:>8s}  {'n/a':>11s}  {'n/a':>8s}  "
                  f"[WARN] one or both names not found among tree leaves — sample-name mismatch "
                  f"between Module 4's cluster file and this module's alignment/tree; check naming")
            rows_out.append((a, b, snp_dist, "", "", "name_mismatch"))
            continue

        mrca = find_mrca(root, a, b)
        mrca_leaves = sorted(get_leaves(mrca))
        extra = [x for x in mrca_leaves if x not in (a, b)]
        support = mrca.support if mrca.support not in (None, "") else None
        try:
            support_val = float(support) if support is not None else None
        except ValueError:
            support_val = None

        if len(mrca_leaves) == 2:
            if support_val is not None and support_val >= UFBOOT_WELL_SUPPORTED:
                verdict = "CONSISTENT: exclusive sister pair, well-supported (UFBoot >= 95%)"
            elif support_val is not None:
                verdict = "PARTIALLY CONSISTENT: exclusive sister pair, but WEAK support (UFBoot < 95%)"
            else:
                verdict = "exclusive sister pair, support unlabeled (root-adjacent branch?)"
        else:
            verdict = (f"INCONSISTENT: SNP-close pair is NOT phylogenetically exclusive — "
                       f"{', '.join(extra)} fall within the same clade")

        support_display = f"{support_val:.0f}" if support_val is not None else "n/a"
        pair_label = f"{a}/{b}"
        print(f"{pair_label:30s} {snp_dist:>8s}  {len(mrca_leaves):>11d}  {support_display:>8s}  {verdict}")
        rows_out.append((a, b, snp_dist, len(mrca_leaves), support_display, verdict))

    out_tsv = outdir / "tree_snp_cluster_crossref.tsv"
    with open(out_tsv, "w") as fh:
        fh.write("sample_a\tsample_b\tsnp_distance\tmrca_leaf_count\tufboot_support\tverdict\n")
        for r in rows_out:
            fh.write("\t".join(str(x) for x in r) + "\n")
    print(f"\nWrote {out_tsv}")
    print("\nReminder: even a CONSISTENT, well-supported result above is evidence of shared ancestry,")
    print("not proof of direct transmission — see this module's README and Script 6's report caveat.")
    print("==================================================================")


if __name__ == "__main__":
    main()
