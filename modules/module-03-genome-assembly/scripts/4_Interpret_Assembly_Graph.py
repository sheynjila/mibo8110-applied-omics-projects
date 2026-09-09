#!/usr/bin/env python3
"""
Module 3, Script 4 — Interpret the Assembly Graph, Not Just the Final FASTA
==============================================================================
NEW IDEA on top of scripts 2-3: those scripts produced a final contigs.fasta
(SPAdes) or assembly.fasta (Flye). Most beginners stop there and treat the
FASTA as ground truth. This script exposes the assembly GRAPH those FASTA
files were collapsed from -- because the graph is where "this looks like a
single clean genome" and "this genome has an unresolved repeat that could be
one contig or three" become visibly different things.

WHY THIS MATTERS: both SPAdes and Flye report their true internal structure
as a graph (GFA format) with nodes (contig-like segments) and edges (adjacency
between them, from overlapping k-mers or read alignments). A perfectly clean,
single-chromosome bacterial genome collapses to ONE long path with no branches
-- a straight line. A genome with unresolved repeats (rRNA operons, transposons,
duplicated genes -- all common in real bacteria) produces BRANCH POINTS: nodes
with more than one incoming or outgoing edge, because the assembler could not
determine which path is the "real" one without more information (e.g. longer
reads, mate-pair distance, or manual curation). The final FASTA hides this --
it just picks one path (or reports fragments) and moves on. This script does
not hide it.

WHAT THIS SCRIPT DOES (pure-Python GFA parser, no dependencies):
  1. Parses a GFA file's S (segment/node) and L (link/edge) lines.
  2. Reports node count, edge count, and total sequence length.
  3. Detects branch points: any node where more than one distinct edge
     touches the same end (its "+" or "-" orientation), which means the
     assembler had more than one way to continue the path there.
  4. Flags likely-circular contigs: a node with a self-referencing edge
     (its own end links back to its own start) -- the classic signature of
     a fully closed circular chromosome or plasmid.
  5. Prints a plain-English summary a student can act on, not just numbers.

USAGE:
    python3 4_Interpret_Assembly_Graph.py <assembly_graph.gfa>

WORKS ON BOTH ASSEMBLER OUTPUTS:
    SPAdes: assembly_graph_with_scaffolds.gfa
    Flye:   assembly_graph.gfa
"""

import sys
import argparse
from collections import defaultdict


def parse_gfa(path):
    """Read S (segment) and L (link) lines from a GFA file.

    Returns:
        segments: dict {name: length_in_bp}
        links: list of (from_name, from_orient, to_name, to_orient)
    """
    segments = {}
    links = []
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            fields = line.split("\t")
            record_type = fields[0]

            if record_type == "S":
                # S <name> <sequence-or-*> [tags...]
                name = fields[1]
                seq = fields[2]
                if seq == "*":
                    # Sequence omitted (common for large graphs) -- look for
                    # an LN:i:<length> tag instead of assuming length 0.
                    length = 0
                    for tag in fields[3:]:
                        if tag.startswith("LN:i:"):
                            length = int(tag.split(":")[-1])
                            break
                else:
                    length = len(seq)
                segments[name] = length

            elif record_type == "L":
                # L <from> <from_orient> <to> <to_orient> <overlap> [tags...]
                from_name, from_orient, to_name, to_orient = fields[1:5]
                links.append((from_name, from_orient, to_name, to_orient))

    return segments, links


def find_branch_points(segments, links):
    """A branch point is a (node, end) occupied by more than one DISTINCT
    edge, where 'end' distinguishes the node's + (head/3') side from its -
    (tail/5') side. More than one distinct edge on the same end means the
    assembler had multiple candidate continuations there and could not
    resolve a single path -- almost always caused by a repeated sequence
    shared by >1 genomic locus.

    IMPORTANT: a self-loop edge (e.g. "L edge_1 + edge_1 +", the standard way
    assemblers mark a fully closed circular contig) touches the SAME end from
    both its 'from' and 'to' roles. That is one edge, not two competing
    edges, so it must be de-duplicated by edge identity -- otherwise every
    circular single-contig genome would be falsely reported as a branch
    point when it is actually the simplest, cleanest topology possible.
    """
    end_occupants = defaultdict(set)   # end -> set of edge indices touching it
    end_neighbors = defaultdict(dict)  # end -> {edge index: (neighbor_node, neighbor_orient)}
    for edge_idx, (f, fo, t, to) in enumerate(links):
        end_a = (f, fo)
        end_b = (t, to)
        end_occupants[end_a].add(edge_idx)
        end_occupants[end_b].add(edge_idx)
        end_neighbors[end_a][edge_idx] = (t, to)
        end_neighbors[end_b][edge_idx] = (f, fo)

    branch_points = {}
    for end, edge_ids in end_occupants.items():
        if len(edge_ids) > 1:
            branch_points[end] = [end_neighbors[end][eid] for eid in sorted(edge_ids)]
    return branch_points


def find_self_loops(links):
    """A node whose own end links back to itself is a strong circularity
    signal -- the assembler found the sequence forms a closed loop, which is
    the expected topology for a fully assembled circular chromosome or
    plasmid (common in bacteria), NOT a defect.
    """
    self_loops = []
    for (f, fo, t, to) in links:
        if f == t:
            self_loops.append((f, fo, to))
    return self_loops


def main():
    parser = argparse.ArgumentParser(
        description="Interpret an assembly graph (GFA) for branch points and circularity."
    )
    parser.add_argument("gfa_path", help="Path to a GFA file (SPAdes or Flye output)")
    args = parser.parse_args()

    segments, links = parse_gfa(args.gfa_path)

    if not segments:
        print(f"ERROR: no S (segment) lines found in {args.gfa_path} -- is this a valid GFA file?", file=sys.stderr)
        sys.exit(1)

    total_len = sum(segments.values())
    branch_points = find_branch_points(segments, links)
    self_loops = find_self_loops(links)

    print("==========================================================")
    print(" MODULE 3 — ASSEMBLY GRAPH INTERPRETATION")
    print("==========================================================")
    print(f"Graph file: {args.gfa_path}")
    print(f"Nodes (segments): {len(segments)}")
    print(f"Edges (links): {len(links)}")
    print(f"Total node sequence length: {total_len:,} bp")
    print("")

    print("---------------------------------------------------------")
    print(" TOPOLOGY")
    print("---------------------------------------------------------")
    if len(links) == 0 and len(segments) == 1:
        print("Single node, zero edges: this is the simplest possible graph --")
        print("one linear sequence with no adjacency information recorded.")
        print("(This is typical of Flye's final assembly_graph.gfa for a")
        print("single-contig genome, since Flye resolves the graph internally")
        print("before writing the final assembly.)")
    elif not branch_points:
        print("No branch points detected: every node's ends connect to at most")
        print("one neighbor. If the graph forms a single connected path, that")
        print("is the strongest topology-based evidence this assembly reflects")
        print("one full, unambiguous genome sequence.")
    else:
        print(f"{len(branch_points)} branch point(s) detected. Each one is a node-end")
        print("touched by more than one edge -- the assembler found more than one")
        print("way to continue the path here and could not choose. This is the")
        print("signature of an UNRESOLVED REPEAT (e.g. an rRNA operon, a")
        print("transposon, or a duplicated gene present at >1 genomic location).")
        print("A repeat this long relative to your read length cannot be spanned")
        print("by any assembler without additional evidence (longer reads,")
        print("mate-pair/linked-read distance, or manual scaffolding).")
        print("")
        shown = 0
        for (node, orient), neighbors in branch_points.items():
            if shown >= 10:
                print(f"  ... and {len(branch_points) - shown} more branch point(s) (truncated)")
                break
            neighbor_desc = ", ".join(f"{n}{o}" for n, o in neighbors)
            print(f"  Node {node}{orient} connects to {len(neighbors)} neighbors: {neighbor_desc}")
            shown += 1

    print("")
    print("---------------------------------------------------------")
    print(" CIRCULARITY")
    print("---------------------------------------------------------")
    if self_loops:
        print(f"{len(self_loops)} self-referencing edge(s) detected -- likely CIRCULAR contig(s):")
        for (node, orient, to_orient) in self_loops:
            print(f"  Node {node}: end {orient} links back to its own end {to_orient}")
        print("")
        print("A circular signature is a GOOD sign for a bacterial chromosome or")
        print("plasmid (both are naturally circular molecules) -- it means the")
        print("assembler found the sequence closes on itself with no gap.")
    else:
        print("No self-referencing edges detected. This assembly's contigs are")
        print("linear (not reported as closed circles) -- which is expected for")
        print("a fragmented or incomplete assembly, but would be UNEXPECTED if")
        print("you know the target organism's chromosome is circular and your")
        print("assembly is a single contig at close to the expected genome size")
        print("(worth double-checking assembly_info.txt / contigs.paths for a")
        print("circularity flag your assembler may report separately from the GFA).")

    print("")
    print("This script describes the graph topology only. It does NOT tell you")
    print("whether the assembly is CORRECT, COMPLETE, or free of CONTAMINATION --")
    print("see 5_Compute_Assembly_Statistics.py and")
    print("6_Assess_Completeness_And_Contamination.py for those dimensions.")
    print("==========================================================")


if __name__ == "__main__":
    main()
