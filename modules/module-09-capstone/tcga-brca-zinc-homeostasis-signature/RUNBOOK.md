# Capstone — Zinc-Homeostasis BC Signature: Runbook

This is an **operational** document — step-by-step execution instructions,
expected outputs, and troubleshooting. For the guiding question, data
provenance, and the required submission-package narrative, see this
project's [`README.md`](README.md).

## Prerequisites

- R (paper's own version: 4.4.0) with `TCGAbiolinks`, `GEOquery`, `DESeq2`,
  `GSVA`, `WGCNA`, `glmnet`, `survival`, `timeROC`, `rms`, `ggvenn`,
  `clusterProfiler`, `Seurat`, `Monocle`, `immunedeconv`, `oncoPredict` — see
  README "Environment" for exact versions.
- Network access to the GDC API (`GDCquery`/`GDCdownload`) and NCBI GEO
  (`getGEOSuppFiles`) — this repository's other modules only ever need SRA;
  this project is the first to need GDC directly.
- Substantial local storage and time for Script 1: TCGA-BRCA's full
  gene-expression quantification set is large (over 1,000 samples); GDC
  downloads and `GDCprepare()` are not fast operations. Consider a scratch
  or project-scoped `GDCdata/` location, not a home directory.
- The real 51-gene zinc-homeostasis list from Guo et al. 2024 (README "Known,
  disclosed gaps" #1) — Tier A's `zh_genes` has only 6 of them as delivered.

## Execution order

| Step | Command | Expected result |
|---|---|---|
| 1 | `Rscript 1_Acquire_TCGA_BRCA_Bulk_Data.R` | `TCGA_BRCA_raw_counts.tsv`, `TCGA_BRCA_tpm.tsv`, `TCGA_BRCA_clinical.csv` |
| 2 | `Rscript 2_Acquire_GSE114725_SingleCell_Data.R` | `geo_supp/GSE114725/<downloaded files>` |
| 2.5 | Inspect Step 2's downloaded files; if not the 10x trio, adjust Tier A's Step 6 loader (see Script 2's header and README gap #2) | A confirmed, working loader for GSE114725's actual format |
| 2.6 | Replace Tier A's `zh_genes` list with the real 51-gene list (README gap #1) | `zh_genes` reflects Guo et al. 2024, not a 6-gene placeholder |
| 3 | `Rscript PRJNA482620_TCGA_BRCA_downstream.R` (from the same working directory as Steps 1-2's output) | `variant_quality_qc.pdf`-equivalent plots from this script's own `pdf()` calls, `outbreak_clusters.tsv`-equivalent intermediate objects in the R session, a fitted `cv_fit` LASSO Cox object, a Kaplan-Meier plot, and a single-cell UMAP |
| 3.5 | Export `cv_fit`'s coefficients (see `3_Validate_Against_Published_Coefficients.py`'s header for the exact snippet) | `derived_coefficients.tsv` |
| 4 | `python3 3_Validate_Against_Published_Coefficients.py --derived-tsv derived_coefficients.tsv` | Console sign/magnitude QA report; `coefficient_qa_summary.tsv` |
| 5 | Fill in README's Interpretation and Reflection sections by hand | Completed submission package |

## Troubleshooting

| Symptom | Likely cause | Fix |
|---|---|---|
| Script 1 fails during `GDCquery`/`GDCdownload` | No network access to the GDC API from this environment, or a `TCGAbiolinks` version mismatch with GDC's current API | Confirm network access first; check `TCGAbiolinks` is a recent version (GDC's API has broken older client versions before) |
| Script 1's sample count differs noticeably from the paper's 1,066 tumor + 99 normal | Expected — TCGA-BRCA on GDC has drifted since the paper's own download date (Script 1's header) | Do not force the count to match; note the actual count and date in this project's "Data provenance" |
| Script 1 reports many patients with missing survival time | Real, expected gap in TCGA's clinical annotation (Script 1's own `[WARN]`) | These patients are dropped by `Surv()`-based analysis downstream; do not impute |
| Script 2 downloads a file that isn't obviously 10x-formatted | Expected — see README gap #2; this dataset is InDrop-derived | Inspect the file, then adjust Tier A's `Read10X()` call per Script 2's header guidance — a disclosed correction, not a silent edit |
| Tier A's `Read10X("GSE114725_matrix/")` call fails outright | The directory name doesn't match where Script 2 actually put its output (`geo_supp/GSE114725/`) | Either symlink/copy to `GSE114725_matrix/`, or edit Tier A's path — record which you did |
| Tier A's ssGSEA step runs but the ZH score looks unlike the paper's Fig. 1A | `zh_genes` still has only 6 of the real 51 genes (gap #1) | Retrieve the full list from Guo et al. 2024 before trusting this step's output |
| Script 4 (QA) reports several `SIGN_DISAGREEMENT` genes | Could be real cohort drift, a sample-filtering difference from the paper, or a genuine bug upstream | Investigate before trusting the risk stratification; do not treat disagreement as "the QA script is too strict" without checking the upstream steps first |
| Script 4 fails with "file not found" for `derived_coefficients.tsv` | That file is never written by Tier A as delivered (a real, disclosed gap) | Add the `coef(cv_fit, ...)` export snippet from Script 4's header to a real Tier A run, then re-run Script 4 |

## Completion checklist (submission package)

- [ ] **Data provenance recorded** — the actual GDC/GEO access date and
      resulting sample counts from Steps 1-2 are written into the README's
      "Data provenance" row, not left as the paper's own numbers.
- [ ] **Quality evidence generated** — `coefficient_qa_summary.tsv` reflects
      a real Step 4 run, including any `SIGN_DISAGREEMENT` flags, not just
      "PASS."
- [ ] **Known gaps resolved or explicitly still open** — the 51-gene ZH list
      and the GSE114725 format question (README gaps #1-2) are either fixed,
      with the fix disclosed, or still flagged as open — never silently
      left unresolved without a note.
- [ ] **Interpretation filled in** — referencing Script 4's actual QA output
      and Tier A's actual Kaplan-Meier/UMAP results, not the source paper's
      own reported numbers.
- [ ] **Reflection filled in** — a genuine assessment of what this project's
      own process revealed, not a restatement of the source paper's
      conclusions.
