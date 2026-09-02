# Module 2 — Bulk RNA-seq and Gene Expression

**Status: built.**

**Curriculum module:** [Module 2 — Bulk RNA-seq and Gene Expression](../../docs/curriculum/modular_bioinformatics_curriculum.md).

## Guiding question

Given raw reads from two or more conditions, which genes are reliably differentially expressed, and how confident should I be in that call?

## Learning outcomes

- Validate integer count data and sample metadata before testing anything.
- Explain biological replication, confounding, covariates, and batch effects — and detect a perfectly confounded design programmatically, not just by inspection.
- Filter low-information features using a documented, group-size-aware rule instead of an unexamined magic number.
- Perform exploratory analysis using transformations, PCA, distances, and sample-level plots.
- Specify designs and contrasts in DESeq2 explicitly — an explicit reference level and named contrast, never an alphabetical or positional default.
- Interpret effect sizes, uncertainty, p-values, and false-discovery-rate adjustment as four separate, distinct numbers.
- Produce accessible figures (colorblind-safe, directly labeled) and a biologically cautious report with stated limitations.

## Portfolio artifact (per curriculum)

A reproducible report containing metadata validation, exploratory plots, model specification, differential-expression results, effect-size visualization, interpretation, and limitations. Script 7 below auto-assembles this as `DE_REPORT.md` — sections 1-7 are filled in automatically from the pipeline's own output; sections 8-9 (interpretation, limitations) are intentionally left for the analyst to complete by hand, since no script can respond to biological meaning on its own.

## Scripts in this module

Two tiers, in the same "add exactly one idea" style used throughout this repository:

**Tier A — the wet-lab pipeline (pre-existing, unchanged):** turns raw reads into a gene x sample count matrix. See the table below.

**Tier B — the statistical analysis pipeline (new this session):** seven scripts that take that count matrix from "raw numbers" to "a defensible, documented differential-expression report." Each script adds exactly one new idea on top of the previous one, built around the same real paired-end study used elsewhere in this repository for continuity: **PRJNA1518998** ("Effect of adipocyte depletion on the retina," *Mus musculus*, Vanderbilt University Medical Center; GEO accession GSE345194).

### Tier A — wet-lab pipeline (existing)

| File | Purpose |
|---|---|
| [`scripts/master_rnaseq_pipeline_consolidated.sh`](scripts/master_rnaseq_pipeline_consolidated.sh) | End-to-end Slurm pipeline: SRA download → FastQC → fastp → HISAT2 alignment → featureCounts. One script parameterized by `REFERENCE_TYPE` (bacterial/viral/human/dynamic-custom) instead of 8 separate hardcoded variants. Ready-to-run Quick Start block includes verified GRCm39 (mouse) URLs and real `PRJNA1518998` SRR accessions. |
| [`scripts/smoke_test_prjna1518998.sh`](scripts/smoke_test_prjna1518998.sh) | Two-sample smoke-test wrapper — runs the full pipeline end-to-end on 2 real accessions before committing to a full cohort submission. |

### Tier B — statistical analysis pipeline (new)

| # | File | New idea introduced |
|---|---|---|
| 1 | [`scripts/1_Validate_Counts_And_Metadata.R`](scripts/1_Validate_Counts_And_Metadata.R) | Strip featureCounts annotation columns automatically, confirm counts are non-negative integers (not TPM/FPKM), and match every count column to its metadata row **by sample name**, never by position — the exact bug the original `rnaseq_deseq2_analysis.R` had silently. |
| 2 | [`scripts/2_Assess_Design_And_Confounding.R`](scripts/2_Assess_Design_And_Confounding.R) | Cross-tabulate condition against every other metadata column (batch, sex, etc.) and programmatically flag any column that is **perfectly confounded** with condition — a design flaw no amount of downstream modeling can fix. |
| 3 | [`scripts/3_Filter_Low_Count_Genes.R`](scripts/3_Filter_Low_Count_Genes.R) | Replace the naive `rowSums(counts) >= 10` rule with a group-size-aware rule (`>= 10 reads in >= size of the smallest condition group`), and print a side-by-side comparison showing exactly which genes the naive rule would have wrongly kept. |
| 4 | [`scripts/4_Exploratory_Analysis_VST_PCA.R`](scripts/4_Exploratory_Analysis_VST_PCA.R) | Variance-stabilize the filtered counts, plot PCA and a sample-distance heatmap, and programmatically flag any sample sitting >3 SD from its group's PCA centroid — a first look at the data before any model is fit. |
| 5 | [`scripts/5_Specify_Design_Formula_And_Contrasts.R`](scripts/5_Specify_Design_Formula_And_Contrasts.R) | Set an **explicit** reference level with `relevel()` (never let it default alphabetically) and build a **named** contrast vector — the design decision that determines the sign of every fold change in the results. |
| 6 | [`scripts/6_Run_DE_Test_With_Shrinkage_And_FDR.R`](scripts/6_Run_DE_Test_With_Shrinkage_And_FDR.R) | Fit DESeq2 with the named contrast, apply `apeglm` log-fold-change shrinkage (falls back gracefully if `apeglm` isn't installed), and explicitly separate and label all four result numbers: effect size, uncertainty, raw p-value, and FDR-adjusted p-value — with a hard warning against reporting raw p-value counts. |
| 7 | [`scripts/7_Accessible_Figures_And_Report.R`](scripts/7_Accessible_Figures_And_Report.R) | Generate a colorblind-safe (blue/orange), directly-labeled volcano plot and MA plot, then auto-assemble the curriculum's required portfolio report (`DE_REPORT.md`) from every prior script's output. |

All seven Tier B scripts were syntax-checked with `Rscript`'s `parse()` (R's equivalent of `bash -n` — this repository's standard pre-delivery validation). All seven were also **functionally run end-to-end twice** against synthetic count-matrix fixtures: once with pure noise (correctly found 0 significant genes) and once with 30 seeded true positives (correctly recovered them, correctly detected a perfectly-confounded batch design in a regression test, and correctly matched the constructed DESeq2 coefficient name against `resultsNames()`). One real bug found during that functional testing — `vst()`'s default 1000-gene subsampling floor erroring out on smaller filtered gene sets — was fixed rather than worked around; a second cosmetic bug (a PCA legend rendering literal R code as its title instead of "batch") was also fixed and re-verified visually.

## How to run

1. Run the Tier A wet-lab pipeline first (see its own Quick Start block) to produce `gene_counts_clean.tsv`/`counts.txt` and write your own `metadata.csv` (`sample,condition[,batch,...]`).
2. Run Tier B scripts 1→7 in order, each taking the previous script's output:
   - `Rscript 1_Validate_Counts_And_Metadata.R counts.txt metadata.csv`
   - `Rscript 2_Assess_Design_And_Confounding.R metadata.csv`
   - `Rscript 3_Filter_Low_Count_Genes.R counts.txt metadata.csv` → writes `counts_filtered.tsv`
   - `Rscript 4_Exploratory_Analysis_VST_PCA.R counts_filtered.tsv metadata.csv` → writes `pca_plot.pdf`, `sample_distance_heatmap.pdf`
   - `Rscript 5_Specify_Design_Formula_And_Contrasts.R counts_filtered.tsv metadata.csv <reference_level>` → writes `dds_prepared.rds`
   - `Rscript 6_Run_DE_Test_With_Shrinkage_And_FDR.R dds_prepared.rds` → writes `deseq2_results.csv`, `deseq2_results_shrunk.csv`
   - `Rscript 7_Accessible_Figures_And_Report.R deseq2_results_shrunk.csv` → writes `volcano_plot.pdf`, `ma_plot.pdf`, `top_genes.csv`, `DE_REPORT.md`
3. Open `DE_REPORT.md` and fill in Sections 8-9 (interpretation and limitations) by hand — this is the intentional human-judgment step the curriculum requires and no script should try to automate.

See [`RUNBOOK.md`](RUNBOOK.md) for the full step-by-step execution order, required R packages, expected outputs, and troubleshooting.

## Corrections applied (see [changelog](../../docs/corrected_scripts_changelog.md) for full detail)

`rnaseq_deseq2_analysis.R`'s original by-name metadata matching fix (documented in that file's own header) is now a reusable, first-class, standalone validation step (Script 1) instead of logic buried inside the analysis script — and Script 2 adds a check that script never had: programmatic confounding detection.

## Teaching materials

- Lecture deck: [`lecture/Bulk_RNAseq_Lecture.pptx`](lecture/Bulk_RNAseq_Lecture.pptx).
- Runbook: [`RUNBOOK.md`](RUNBOOK.md).
- Still to develop per [Section 5 of the curriculum](../../docs/curriculum/modular_bioinformatics_curriculum.md): guided student workbook, dataset card/provenance record for PRJNA1518998, and a transfer-task assessment.
