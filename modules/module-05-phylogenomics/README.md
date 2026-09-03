# Module 5 — Phylogenomics and Outbreak Context

**Status: built.**

**Curriculum module:** [Module 5 — Phylogenomics and Outbreak Context](../../docs/curriculum/modular_bioinformatics_curriculum.md).

## Guiding question

Given a set of related genomes or variant calls, how closely related are they, and what can (and can't) that relatedness tell me about a real-world outbreak?

## Learning outcomes

- Build a multiple sequence alignment suitable for tree inference — directly from a cohort VCF, with the ascertainment-bias consequence of doing so made explicit, not hidden.
- Infer a tree with IQ-TREE and read branch-support values (not just topology) as a measure of confidence — using the correct support threshold for the specific bootstrap method actually run (UFBoot), not a generically-remembered number.
- Apply the curriculum's core caution explicitly: genomic relatedness is not the same claim as a transmission event — a close phylogenetic distance is necessary but not sufficient evidence of direct transmission. This module makes that caution a structural, always-printed part of its output, not just a sentence in the README.

## Portfolio artifact (per curriculum)

A phylogenetic tree + branch-support summary + a written outbreak interpretation that respects the relatedness-vs-transmission distinction above. Script 6 below auto-assembles this as `PHYLO_REPORT.md` — sections 1-4 are filled in automatically from Scripts 1-5's own output; section 5 (the actual outbreak interpretation) is intentionally left blank, printed alongside the module's core caution as unconditional boilerplate rather than optional filler, mirroring Module 2's `DE_REPORT.md`, Module 3's `ASSEMBLY_REPORT.md`, and Module 4's `VARIANT_REPORT.md`.

## Relationship to existing repo content

This module has no pre-existing "Tier A" pipeline to build on (unlike Modules 2-4) — this draft's own status was "planned, not yet built. No dedicated tree-building script exists yet." It is, however, deliberately **not** an isolated module: as this draft's own "Relationship to existing repo content" note anticipated, [Module 4's `variant_analysis.R`](../module-04-variant-calling/README.md) already produces a merged cohort VCF and a SNP-distance/`outbreak_clusters.tsv` — this module's Script 1 turns that same cohort VCF into an alignment, and Script 5 cross-references its tree directly against Module 4's own cluster calls, rather than duplicating pairwise-distance logic that already exists.

## Worked-example dataset (real, verified via NCBI eutils — same method Modules 3-4 used)

Module 4's own 3-isolate *Salmonella enterica* Typhimurium cohort (BioProject PRJNA242847, GenomeTrakr/USDA-FSIS) has only one possible unrooted topology at 3 taxa — not enough structure to teach "some branches are well-supported, others are not." This module extends that exact cohort with 2 additional isolates from the same verified BioProject, so the tree has multiple internal branches with genuinely heterogeneous support:

| Run | BioSample | Note |
|---|---|---|
| `SRR40474618` | SAMN62852532 | From Module 4's original 3-isolate cohort |
| `SRR40472300` | SAMN62840582 | From Module 4's original 3-isolate cohort |
| `SRR40472280` | SAMN62840173 | From Module 4's original 3-isolate cohort |
| `SRR40426677` | SAMN62799163 | New — added for this module's tree topology |
| `SRR40376221` | SAMN62746372 | New — added for this module's tree topology |

All five are real, verified, paired-end Illumina MiSeq *S.* Typhimurium WGS runs under PRJNA242847. To build this module's worked example, re-run Module 4's pipeline (Tier A `master_snp_pipeline.sh` + Script 4 ploidy correction) against a 5-sample `srr_list.txt` containing all five accessions above, producing a 5-sample `merged_cohort_haploid.vcf` — see this module's [`RUNBOOK.md`](RUNBOOK.md) Step 0 for the exact commands. This module does not re-implement alignment or variant calling; it strictly reuses Module 4's.

## Scripts in this module

No Tier A/Tier B split — this module has one pipeline, six scripts, each adding exactly one new idea, in the same style used throughout this repository:

| # | File | New idea introduced |
|---|---|---|
| 1 | [`scripts/1_Extract_Core_SNP_Alignment_From_VCF.sh`](scripts/1_Extract_Core_SNP_Alignment_From_VCF.sh) | Turn a cohort VCF into a FASTA alignment using only `bcftools` — the missing link this module's own draft explicitly named as not yet existing. Heterozygous calls (an artifact against this haploid reference — see Module 4 Script 4) are coded `N`, not a trusted IUPAC ambiguity code. |
| 2 | [`scripts/2_Assess_Alignment_Quality.py`](scripts/2_Assess_Alignment_Quality.py) | Report per-sample AND per-site missingness on the alignment before trusting a tree built from it — a high-missingness sample will still produce a tree, with a placement driven by missing data, not relatedness, and nothing in the tree image alone would say so. |
| 3 | [`scripts/3_Run_Tree_Inference_IQTREE.sh`](scripts/3_Run_Tree_Inference_IQTREE.sh) | Run IQ-TREE with ModelFinder (`-m MFP+ASC`) and 1000 ultrafast bootstrap replicates (`-B 1000`) — critically, with ascertainment-bias correction (`+ASC`) applied, because Script 1's alignment is SNP-only (invariant sites excluded) and a plain substitution model would silently mis-estimate branch lengths and the model itself as a result. Documents RAxML-NG's equivalent invocation (the curriculum's noted alternative) rather than running it. |
| 4 | [`scripts/4_Interpret_Branch_Support.py`](scripts/4_Interpret_Branch_Support.py) | Parse the tree's UFBoot support values and classify each internal branch against a 95% threshold — IQ-TREE's own recommended cutoff for ultrafast bootstrap, explicitly **not** the classic Felsenstein bootstrap's traditional 70% cutoff, a distinction this script states rather than silently picking one. Pure-Python Newick parsing, no third-party tree library. |
| 5 | [`scripts/5_Cross_Reference_Tree_And_SNP_Clusters.py`](scripts/5_Cross_Reference_Tree_And_SNP_Clusters.py) | Cross-reference Module 4's `outbreak_clusters.tsv` (pairwise SNP-distance calls) against this module's own tree: for each SNP-close pair, is it also an exclusive, well-supported clade — or does the tree place other isolates between them? Disagreement between the two signals is itself informative, not an error. |
| 6 | [`scripts/6_Assemble_Phylogenomics_Report.py`](scripts/6_Assemble_Phylogenomics_Report.py) | Auto-assemble `PHYLO_REPORT.md` from Scripts 1-5, with the objective sections filled in automatically and an interpretation section left blank — printed alongside this module's relatedness-vs-transmission caution as unconditional boilerplate, not optional filler, since the curriculum names that distinction as this module's central point. |

### Verification performed

All shell scripts (`1_Extract_Core_SNP_Alignment_From_VCF.sh`, `3_Run_Tree_Inference_IQTREE.sh`) pass `bash -n`. This session's environment has no `bcftools`, `iqtree2`, or Python interpreter installed, so the Python scripts (2, 4, 5, 6) could not be syntax-checked with `python3 -m py_compile` or functionally run against real or synthetic fixtures the way Modules 2 and 3 were — they received careful manual line-by-line review only. This is stated plainly, matching Module 4's disclosure, rather than claimed away. Completing that functional testing (a small synthetic Newick tree and cohort VCF, run through all six scripts end-to-end) on a machine with the actual toolchain is the top item in "still to develop" below.

## How to run

0. **Prerequisite (not part of this module's own scripts):** re-run Module 4's pipeline against the 5-sample cohort above to produce `merged_cohort_haploid.vcf` — see [`RUNBOOK.md`](RUNBOOK.md) Step 0 for exact commands.
1. `bash 1_Extract_Core_SNP_Alignment_From_VCF.sh merged_cohort_haploid.vcf .` → `core_snp_alignment.fasta`
2. `python3 2_Assess_Alignment_Quality.py core_snp_alignment.fasta .` → `alignment_qc_summary.tsv`
3. `bash 3_Run_Tree_Inference_IQTREE.sh core_snp_alignment.fasta cohort_tree` → `cohort_tree.treefile`, `cohort_tree.iqtree`
4. `python3 4_Interpret_Branch_Support.py cohort_tree.treefile .` → `branch_support_summary.tsv`
5. `python3 5_Cross_Reference_Tree_And_SNP_Clusters.py cohort_tree.treefile outbreak_clusters.tsv .` → `tree_snp_cluster_crossref.tsv`
6. `python3 6_Assemble_Phylogenomics_Report.py --tree cohort_tree.treefile --iqtree-report cohort_tree.iqtree --alignment-qc-tsv alignment_qc_summary.tsv --support-tsv branch_support_summary.tsv --crossref-tsv tree_snp_cluster_crossref.tsv --out PHYLO_REPORT.md`
7. Open `PHYLO_REPORT.md` and fill in Section 5 (outbreak interpretation) by hand, alongside its printed caveat — the intentional human-judgment step the curriculum requires.

See [`RUNBOOK.md`](RUNBOOK.md) for the full step-by-step execution order, required tools, expected outputs, and troubleshooting.

## Teaching materials

- Lecture deck: not yet provided for this module — still to develop.
- Runbook: [`RUNBOOK.md`](RUNBOOK.md).
- Still to develop per [Section 5 of the curriculum](../../docs/curriculum/modular_bioinformatics_curriculum.md): the outstanding functional (not just syntax/manual-review) testing noted above, a guided student workbook, a dataset card/provenance record for the 5-isolate PRJNA242847 cohort, and a transfer-task assessment.
