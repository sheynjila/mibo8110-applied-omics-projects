# Module 2 — Bulk RNA-seq: Runbook

This is an **operational** document — step-by-step execution instructions,
expected outputs, and troubleshooting. For the conceptual "why," see this
module's [`README.md`](README.md). For the full curriculum context, see
[`docs/curriculum/modular_bioinformatics_curriculum.md`](../../docs/curriculum/modular_bioinformatics_curriculum.md).

## Prerequisites

- Slurm environment for Tier A (`scripts/master_rnaseq_pipeline_consolidated.sh`) —
  environment modules are loaded automatically via
  [`scripts/load_modules.sh`](scripts/load_modules.sh); see the root
  [`RUNBOOK.md`](../../RUNBOOK.md) §4 for how the pinned-build /
  discover-by-name fallback works and how to override a module name for
  your own cluster.
- R (interactive session, RStudio, or `Rscript` on a login node) with
  `DESeq2`, `apeglm` (optional — Script 6 falls back gracefully without
  it), and `ggplot2` installed, for Tier B (scripts 1-7).
- Module 1 completed or reviewed — Script 1's "validate before you trust
  it" habit is the same one taught there, applied one level up the
  pipeline (a count matrix instead of raw reads).

## Execution order

| Step | Command | Expected result |
|---|---|---|
| 1 | `bash smoke_test_prjna1518998.sh` (or `sbatch --export=REFERENCE_TYPE=... master_rnaseq_pipeline_consolidated.sh` directly, with your own `srr_list.txt`) | `counts/gene_counts_clean.tsv`, `logs/reference_manifest.txt`, `final_multiqc_report.html` |
| 2 | Write `metadata.csv` (`sample,condition[,batch,...]`) for your samples | — |
| 3 | `Rscript 1_Validate_Counts_And_Metadata.R counts.txt metadata.csv` | Console PASS/FAIL report; stops on the first hard failure |
| 4 | `Rscript 2_Assess_Design_And_Confounding.R metadata.csv` | Console cross-tab; hard failure if any column is perfectly confounded with `condition` |
| 5 | `Rscript 3_Filter_Low_Count_Genes.R counts.txt metadata.csv` | `counts_filtered.tsv`; console before/after gene counts |
| 6 | `Rscript 4_Exploratory_Analysis_VST_PCA.R counts_filtered.tsv metadata.csv` | `pca_plot.pdf`, `sample_distance_heatmap.pdf`; console outlier-sample flags |
| 7 | `Rscript 5_Specify_Design_Formula_And_Contrasts.R counts_filtered.tsv metadata.csv <reference_level>` | `dds_prepared.rds`; console design formula + named contrast |
| 8 | `Rscript 6_Run_DE_Test_With_Shrinkage_And_FDR.R dds_prepared.rds` | `deseq2_results.csv`, `deseq2_results_shrunk.csv`; console effect-size/FDR summary |
| 9 | `Rscript 7_Accessible_Figures_And_Report.R deseq2_results_shrunk.csv` | `volcano_plot.pdf`, `ma_plot.pdf`, `top_genes.csv`, `DE_REPORT.md` (sections 1-7 auto-filled, 8-9 blank) |
| 10 | Fill in Sections 8-9 (interpretation, limitations) of `DE_REPORT.md` by hand | Completed portfolio artifact |

## Expected outputs, by step

- **Step 1 (Tier A):** `raw_reads/`, `qc/`, `trimmed_reads/`, `alignment/<SRR>_sorted.bam[.bai]`, `counts/gene_counts_raw.txt`, `counts/gene_counts_clean.tsv`, `logs/reference_manifest.txt`, `final_multiqc_report.html`.
- **Step 3:** stdout only — validation report.
- **Step 4:** stdout only — confounding cross-tab.
- **Step 5:** `counts_filtered.tsv`.
- **Step 6:** `pca_plot.pdf`, `sample_distance_heatmap.pdf`.
- **Step 7:** `dds_prepared.rds`.
- **Step 8:** `deseq2_results.csv`, `deseq2_results_shrunk.csv`.
- **Step 9:** `volcano_plot.pdf`, `ma_plot.pdf`, `top_genes.csv`, `DE_REPORT.md`.

## Troubleshooting

| Symptom | Likely cause | Fix |
|---|---|---|
| `smoke_test_prjna1518998.sh` reports "manifest does not mention GRCm39" | A previous run used a different `REFERENCE_TYPE`/`FASTA_URL` and this WORKDIR is stale, or `GENOME_TAG` derivation changed | Clear the reported `WORKDIR` and re-run — never reuse a WORKDIR across different genomes |
| Smoke test alignment rate `<80%` | Reference/species mismatch, or a genuinely poor-quality run | Confirm `FASTA_URL`/`GFF_URL` actually match your `srr_list.txt` organism before treating this as a data-quality issue |
| Script 1 fails with "count column does not match any metadata row by name" | `metadata.csv`'s `sample` column doesn't exactly match the BAM-derived column headers in `counts.txt` (case, extra whitespace, `_sorted` suffix not stripped) | Compare the printed column names directly; Tier A's Phase 5 matrix cleanup already strips `alignment/` and `_sorted.bam` — don't re-add them in `metadata.csv` |
| Script 2 aborts: "column X is perfectly confounded with condition" | A real design flaw — every sample in one condition group also shares one value of another metadata column | This cannot be fixed downstream; the experiment (or the specific contrast you can test) is limited by its own design — document this rather than proceeding |
| Script 4 flags a sample as a PCA outlier | Either a real technical issue (check that sample's FastQC/fastp/MultiQC output from Step 1) or genuine biological variability | Don't drop the sample without a documented, non-circular reason (i.e. not "because it's the outlier") |
| Script 5 error: "reference level not found" | The `<reference_level>` argument doesn't exactly match a value in `metadata.csv`'s `condition` column | Check the exact string with `Rscript -e 'levels(factor(read.csv("metadata.csv")$condition))'` |
| Script 6 warns `apeglm` not installed | Optional package missing | Install it (`BiocManager::install("apeglm")`) for shrunk LFCs, or proceed with the unshrunk `deseq2_results.csv` and note this in `DE_REPORT.md` |
| Script 6/7 report far more or fewer significant genes than expected | Usually traced back to Script 5's contrast direction/reference level, not the test itself | Re-check Script 5's printed contrast before re-running the test |
| `DE_REPORT.md` shows "Not available" for a section | The corresponding earlier script wasn't run, or was run in a different directory | Re-run the missing step in the same working directory, then re-run Script 7 |

## Completion checklist (portfolio artifact)

The curriculum's Module 2 assessment asks for: *"A reproducible report
containing metadata validation, exploratory plots, model specification,
differential-expression results, effect-size visualization,
interpretation, and limitations."* Confirm each is present in `DE_REPORT.md`
before considering the module complete:

- [ ] **Metadata validation** — Section 1 reflects a real Script 1 run against your actual `counts.txt`/`metadata.csv`, not a placeholder.
- [ ] **Design/confounding check** — Section 2 states plainly whether any column was flagged as perfectly confounded, and what that means for which contrasts are testable.
- [ ] **Exploratory plots** — Section 3/4 include the real `pca_plot.pdf`/`sample_distance_heatmap.pdf`, with any flagged outlier samples addressed.
- [ ] **Model specification** — Section 5 states the exact design formula, reference level, and named contrast Script 5 actually used.
- [ ] **DE results** — Section 6 reports the real Script 6 effect size, uncertainty, raw p-value, and FDR-adjusted p-value as four distinct numbers, not conflated.
- [ ] **Effect-size visualization** — Section 7 includes the real `volcano_plot.pdf`/`ma_plot.pdf`.
- [ ] **Interpretation and limitations** — Sections 8-9 are filled in by hand, referencing the actual design, confounding, and DE findings above — not left as the auto-generated prompt text.
