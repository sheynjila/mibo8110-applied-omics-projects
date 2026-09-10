# Module 5 — Phylogenomics and Outbreak Context: Runbook

This is an **operational** document — step-by-step execution instructions,
expected outputs, and troubleshooting. For the conceptual "why," see this
module's [`README.md`](README.md). For the full curriculum context, see
[`docs/curriculum/modular_bioinformatics_curriculum.md`](../../docs/curriculum/modular_bioinformatics_curriculum.md).

## Prerequisites

- Module 4 completed — this module has no wet-lab or variant-calling
  pipeline of its own; it starts from a cohort VCF Module 4 produces.
- `bcftools` on `PATH` — used by Script 1 (plain shell, checks via
  `command -v`; if missing, load it via
  [`scripts/load_modules.sh`](scripts/load_modules.sh)'s `load_bcftools`
  in your shell first).
- IQ-TREE2 — Script 3 loads it automatically via
  [`scripts/load_modules.sh`](scripts/load_modules.sh) (pinned to
  `IQ-TREE/2.3.6-gompi-2024a`, with an automatic fallback if your cluster
  doesn't have that exact build — see the root
  [`RUNBOOK.md`](../../RUNBOOK.md) §4).
- Python 3, standard library only — Scripts 2, 4, 5, and 6 use no
  third-party packages, deliberately (same choice Module 3 made for its
  graph-parsing scripts), so this module has no dependency chain beyond
  `bcftools` and IQ-TREE2 themselves.

## Step 0 — build the 5-isolate cohort VCF (in Module 4, not this module)

This module's tree needs more than 3 taxa to show heterogeneous branch
support (see README). Re-run Module 4's pipeline against the extended
cohort:

```bash
cd path/to/modules/module-04-variant-calling/scripts
cat > srr_list.txt <<'EOF'
SRR40474618
SRR40472300
SRR40472280
SRR40426677
SRR40376221
EOF
bash master_snp_pipeline.sh          # or sbatch, per Module 4's own RUNBOOK
bash 4_Check_Ploidy_Assumptions.sh   # writes <SRR>_raw_haploid.vcf x5

cd "$WORKDIR/variants"   # WORKDIR defaults to /scratch/$(whoami)/master_snp_pipeline
for f in *_raw_haploid.vcf; do
    bcftools view -Oz -o "${f%.vcf}.vcf.gz" "$f"
    bcftools index -t "${f%.vcf}.vcf.gz"
done
bcftools merge *_raw_haploid.vcf.gz -Ov -o merged_cohort_haploid.vcf
```

Copy or symlink `merged_cohort_haploid.vcf` into this module's own working
directory before Step 1 below.

## Execution order

| Step | Command | Expected result |
|---|---|---|
| 1 | `bash 1_Extract_Core_SNP_Alignment_From_VCF.sh merged_cohort_haploid.vcf .` | `core_snp_alignment.fasta` (6 sequences: 5 isolates + `Reference`), `snp_sites.tsv` |
| 2 | `python3 2_Assess_Alignment_Quality.py core_snp_alignment.fasta .` | Console per-sample/per-site missingness report; `alignment_qc_summary.tsv` |
| 3 | `bash 3_Run_Tree_Inference_IQTREE.sh core_snp_alignment.fasta cohort_tree` | `cohort_tree.treefile`, `cohort_tree.iqtree`, `cohort_tree.log` |
| 4 | `python3 4_Interpret_Branch_Support.py cohort_tree.treefile .` | Console UFBoot support table; `branch_support_summary.tsv` |
| 5 | `python3 5_Cross_Reference_Tree_And_SNP_Clusters.py cohort_tree.treefile outbreak_clusters.tsv .` | Console consistency report; `tree_snp_cluster_crossref.tsv` |
| 6 | `python3 6_Assemble_Phylogenomics_Report.py --tree cohort_tree.treefile --iqtree-report cohort_tree.iqtree --alignment-qc-tsv alignment_qc_summary.tsv --support-tsv branch_support_summary.tsv --crossref-tsv tree_snp_cluster_crossref.tsv` | `PHYLO_REPORT.md`, sections 1-4 auto-filled, section 5 blank |
| 7 | Fill in Section 5 (outbreak interpretation) of `PHYLO_REPORT.md` by hand, alongside its printed caveat | Completed portfolio artifact |

Step 5's `outbreak_clusters.tsv` comes from Module 4's `variant_analysis.R`,
run against the same 5-sample `merged_cohort_haploid.vcf` from Step 0 — run
it (or re-run it) before Step 5 if you have not already. If no isolate pair
falls within `variant_analysis.R`'s `SNP_THRESHOLD`, no file is written;
Script 5 reports this as a legitimate outcome, not an error, and exits
cleanly.

## Expected outputs, by step

- **Step 1:** `core_snp_alignment.fasta`, `snp_sites.tsv`.
- **Step 2:** `alignment_qc_summary.tsv`.
- **Step 3:** `cohort_tree.treefile`, `cohort_tree.iqtree`, `cohort_tree.log`.
- **Step 4:** `branch_support_summary.tsv`.
- **Step 5:** `tree_snp_cluster_crossref.tsv` (or nothing, if Module 4 found no SNP-close pairs).
- **Step 6:** `PHYLO_REPORT.md`.

## Troubleshooting

| Symptom | Likely cause | Fix |
|---|---|---|
| Script 1 fails with `no biallelic SNP sites found` | The input VCF has only indels/multiallelic sites, or is genuinely empty/invariant across the cohort | Confirm Step 0's `merged_cohort_haploid.vcf` actually contains variant records: `bcftools view -H merged_cohort_haploid.vcf \| wc -l` |
| Script 1's output has more/fewer sequences than expected | `bcftools query -l` sample count does not match your `srr_list.txt` — one sample failed in Module 4's fault-tolerant loop and was excluded from the merge | Check Module 4's `master_snp_pipeline.sh` console output for `FAILED_SAMPLES`; re-run that sample or accept the smaller cohort and note it in `PHYLO_REPORT.md` Section 5 |
| Script 2 flags a sample with high missingness | That isolate had low mapping rate or coverage gaps in Module 4 (see Module 4's `mapping_bias_summary.tsv`) | Cross-check against Module 4's Script 2 output for that sample before deciding whether to exclude it and re-run from Step 1 |
| Script 3 (IQ-TREE) errors on `+ASC` with a message about invariant sites | The alignment is not actually SNP-only — an invariant site slipped through, most likely because the input VCF was not restricted to variant sites (e.g., you accidentally passed a VCF containing reference-only/non-variant records) | Confirm the input to Script 1 is a true multi-sample **variant** VCF from `bcftools call -mv`, not a genotyped-all-sites VCF |
| Script 3 runs but the log shows very few or one candidate model tested | Alignment has too few sites (very few SNPs across the cohort) for ModelFinder to distinguish models meaningfully | Expected for a small, closely related cohort; note this explicitly in `PHYLO_REPORT.md` Section 5 rather than treating the reported model as strongly evidenced |
| Script 4 reports "no internal branches found" | Fewer than 4 taxa in the alignment (an unrooted tree of 3 taxa has exactly one, unlabeled internal branch) | Use the full 5-isolate cohort from Step 0, not Module 4's original 3-isolate one, if you want a support analysis at all |
| Script 4/5 report a support value that looks unexpectedly low for a clade you expected to be solid | UFBoot support reflects the ALIGNMENT's actual signal, including any ascertainment/missingness issues Script 2 already flagged — a "solid-looking" clade in a coverage-poor alignment is exactly the case this module's pipeline is built to catch, not a false negative to explain away | Re-check Script 2's per-sample missingness for every taxon in that clade before assuming the number is simply wrong |
| Script 5 reports `INCONSISTENT` for a pair Module 4 flagged as SNP-close | The tree places other isolates between that pair — i.e., they share SNP-level similarity to the reference but are not each other's closest relative in the full alignment | This is the module's central teaching point, not a bug — report it as exactly that disagreement in `PHYLO_REPORT.md` Section 5, not as a discarded result |
| Script 6's `PHYLO_REPORT.md` shows "Not available" for a section | The corresponding earlier script's output file was not passed via its `--*-tsv`/`--iqtree-report` flag, or was not generated for a different reason (e.g., Script 4 above) | Re-run the missing step, then re-run Script 6 with the correct file path |
| Script 6's Section 5 is blank | Expected — intentionally left for the analyst, alongside the printed caveat | Fill it in referencing Sections 1-4's actual numbers, and do not remove the caveat text above it |

## Completion checklist (portfolio artifact)

The curriculum's Module 5 assessment asks for: *"A phylogenetic tree +
branch-support summary + a written outbreak interpretation that respects
the relatedness-vs-transmission distinction."* Confirm each is present
before considering the module complete:

- [ ] **Phylogenetic tree** — `cohort_tree.treefile` (Script 3) is retained
      alongside the report, built with ascertainment-bias correction
      (`+ASC`), not a plain substitution model.
- [ ] **Branch-support summary** — Section 3 of `PHYLO_REPORT.md` reflects
      a real Script 4 run, using the 95% UFBoot threshold (not the classic
      70% bootstrap cutoff), against your actual tree.
- [ ] **Tree-vs-SNP-cluster consistency** — Section 4 reflects a real
      Script 5 run against your actual Module 4 `outbreak_clusters.tsv`,
      including any `INCONSISTENT` results, not just agreeing ones.
- [ ] **Written outbreak interpretation** — Section 5 is filled in by
      hand, explicitly distinguishing what the tree/SNP evidence supports
      (shared ancestry) from what it does NOT establish on its own (a
      confirmed transmission event) — the caveat text above it is not
      removed or paraphrased away.
