# Capstone Project — A Zinc-Homeostasis Gene Signature for Breast Cancer Prognosis (Bulk + Single-Cell Integration)

**Status: built as a scaffold — data acquisition, Tier A analysis, and a QA script are real and ready to run; not executed in this environment (no R, no network access to GDC/GEO in this session — see "Verification performed").**

**Parent module:** [Module 9 — Integrated Capstone and Professional Portfolio](../README.md). This is the first concrete `modules/module-09-capstone/<project-slug>/` instance the parent module's own README asked for.

## Guiding question

Do zinc-homeostasis-related genes carry independent prognostic information in breast cancer, and does a compact gene signature built from bulk RNA-seq — cross-checked against single-cell data for which cell types actually express those genes — stratify patients into meaningfully different survival groups?

This project draws on **Module 2** (bulk RNA-seq / DESeq2, the differential-expression half of this analysis) and **Module 7** (single-cell QC/clustering, the cell-type-identification half) — exactly the "combine whatever modules you've completed" synthesis the parent module's README describes, plus survival modeling (LASSO Cox, Kaplan-Meier) that no other module in this repository covers on its own.

## Data source (real, with full provenance)

This project reproduces a real, published analysis: **Xu, W., Zheng, Y. & Xu, W. "Construction of a novel zinc homeostasis-related gene signature to predict breast cancer prognosis by integrating single-cell and bulk RNA-seq data." *Sci Rep* (2026). https://doi.org/10.1038/s41598-026-68769-z** — found among this project's source materials as `s41598-026-68769-z_reference.pdf`, initially unread; reading it in full is what turned a set of draft scripts assuming pre-processed input files into a project with an actual, citable data-acquisition plan. See Scripts 1-2 for how each dataset is actually pulled, not assumed to already exist:

| Dataset | Accession / portal | Role | Real N (per the paper) |
|---|---|---|---|
| TCGA-BRCA | [GDC portal](https://portal.gdc.cancer.gov/), project `TCGA-BRCA` | Bulk RNA-seq training cohort | 1,066 tumor + 99 normal |
| GSE20685 | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE20685) | Validation cohort (not yet scripted — see "Still to develop") | 327 samples |
| GSE96058 | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE96058) | External cohort (not yet scripted — see "Still to develop") | 2,969 samples |
| GSE114725 | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE114725) / SRA study SRP148597 | Single-cell immune profiling (Azizi et al. 2018, *Cell*) | 8 tumor + 8 normal breast tissue samples |

This project's own scope is the **training cohort (TCGA-BRCA) + GSE114725** only — GSE20685/GSE96058 validation is real, citable, and a natural extension, but out of scope for this instance (see "Still to develop").

## The 8 phases, applied to this project

1. **Proposal** — this README's Guiding question.
2. **Data and compute plan** — Scripts 1-2 (below), and the "Environment" section.
3. **Analysis plan** — Script 3 of `PRJNA482620_TCGA_BRCA_downstream.R` (Tier A): ssGSEA zinc-homeostasis scoring → WGCNA module identification → DESeq2 DEGs → LASSO Cox 9-gene signature → Kaplan-Meier → single-cell UMAP.
4. **Pilot run** — not yet performed in this environment (see "Verification performed").
5. **Full execution** — not yet performed.
6. **Quality-assurance audit** — Script 3, `3_Validate_Against_Published_Coefficients.py`: compares a real re-run's LASSO coefficients against the paper's own published 9-gene signature, gene by gene, by sign and magnitude.
7. **Interpretive report** — this README's "Interpretation" section, intentionally left for the analyst (see below).
8. **Presentation / peer review** — not yet performed.

## Required submission package

| Component | This project's content |
|---|---|
| **README** | This file. |
| **Environment** | R packages, with the exact versions the source paper itself reports using (see its Methods): `TCGAbiolinks`, `GEOquery`, `DESeq2` (1.44.0), `GSVA` (2.0.2, for ssGSEA), `WGCNA` (1.72-5), `glmnet` (4.1-8), `survival` (3.5-8), `timeROC` (0.4), `rms` (8.1-1), `ggvenn` (0.1.10), `clusterProfiler` (4.14.3), `Seurat` (4.1.0), `Monocle`/Monocle2 (2.22.0), `immunedeconv` (2.1.0), `oncoPredict` (1.2). R version 4.4.0 (per the paper). No container/lockfile is provided yet — see "Still to develop". |
| **Workflow** | `scripts/1_Acquire_TCGA_BRCA_Bulk_Data.R` → `scripts/2_Acquire_GSE114725_SingleCell_Data.R` → `scripts/PRJNA482620_TCGA_BRCA_downstream.R` → `scripts/3_Validate_Against_Published_Coefficients.py`. |
| **Data provenance** | Table above. GDC/GEO access dates are whenever Scripts 1-2 are actually run — record that date when you run them; TCGA-BRCA and GEO are living databases and will have drifted since the source paper's own download (see Script 1's header). |
| **Quality evidence** | Script 3's `coefficient_qa_summary.tsv` (sign/magnitude agreement with the published model) — this is the "Quality-assurance audit" phase's actual artifact. |
| **Results** | Not yet generated — this environment cannot run R/DESeq2/WGCNA/Seurat (see "Verification performed"). Running Scripts 1→2→Tier A→3 on a machine with the real toolchain produces `TCGA_BRCA_raw_counts.tsv`/`TCGA_BRCA_tpm.tsv`/`TCGA_BRCA_clinical.csv` (Script 1), `geo_supp/GSE114725/` (Script 2), `variant_quality_qc.pdf`-equivalent outputs from Tier A's own `pdf()` calls, and `coefficient_qa_summary.tsv` (Script 3). |
| **Interpretation** | **[FILL IN]** — see below. |
| **Reflection** | **[FILL IN]** — see below. |

## Known, disclosed gaps (real limitations, not silently worked around)

1. **The zinc-homeostasis gene list is incomplete in Tier A.** `PRJNA482620_TCGA_BRCA_downstream.R`'s `zh_genes` list has 6 genes with a comment claiming "51 total genes." The source paper's own Methods name the real source: *"A total of 51 genes associated with ZH were identified in a previous study[51]"* — reference 51 is **Guo, T. et al. "Deciphering the role of zinc homeostasis in the tumor microenvironment and prognosis of prostate cancer." *Discover Oncol.* 15, 207 (2024). https://doi.org/10.1007/s12672-024-01006-z**. That paper (not in this project's source materials) is where the real 51-gene list lives. Per this repository's own standards, the other 45 gene symbols are **not fabricated here** — retrieve the real list from Guo et al. 2024 before running Tier A's ssGSEA step for real; running it against 6 genes instead of 51 would not reproduce the paper's own ZH score.
2. **GSE114725's actual file format is unconfirmed.** See `2_Acquire_GSE114725_SingleCell_Data.R`'s header — this dataset is InDrop-derived (confirmed independently in Module 7's own investigation), and `Read10X()` in Tier A's Step 6 assumes the 10x-specific barcodes/features/matrix trio. Run Script 2, inspect what actually downloads, and adjust the loader if it is not in that format — do not assume it will work unmodified.
3. **Tier A never exports its LASSO coefficients to a file.** Script 3 (`3_Validate_Against_Published_Coefficients.py`) needs them; see that script's header for the exact `coef(cv_fit, ...)` snippet to add before a real comparison is possible.
4. **GSE20685 and GSE96058 (the paper's validation/external cohorts) are not scripted here.** A real, citable, natural extension — not attempted in this instance to keep scope to one training cohort + one single-cell dataset, matching a capstone's actual size, not the full publication's.

## Scripts in this project

| File | Purpose |
|---|---|
| [`scripts/1_Acquire_TCGA_BRCA_Bulk_Data.R`](scripts/1_Acquire_TCGA_BRCA_Bulk_Data.R) | Pull TCGA-BRCA raw counts, TPM, and clinical/survival data from GDC via `TCGAbiolinks`; derive `time`/`status` survival fields Tier A's `Surv()` call expects. |
| [`scripts/2_Acquire_GSE114725_SingleCell_Data.R`](scripts/2_Acquire_GSE114725_SingleCell_Data.R) | Download GSE114725's own GEO-deposited supplementary files via `GEOquery`, with an explicit, un-silenced format check against Tier A's `Read10X()` assumption. |
| [`scripts/PRJNA482620_TCGA_BRCA_downstream.R`](scripts/PRJNA482620_TCGA_BRCA_downstream.R) | **Tier A, unchanged, pre-existing.** ssGSEA → WGCNA → DESeq2 → LASSO Cox → Kaplan-Meier → single-cell UMAP. |
| [`scripts/3_Validate_Against_Published_Coefficients.py`](scripts/3_Validate_Against_Published_Coefficients.py) | Quality-assurance audit: compares a real re-run's LASSO coefficients against the paper's own published 9-gene signature (hardcoded from its Fig. 2C), gene by gene. |

### Verification performed

`3_Validate_Against_Published_Coefficients.py` — no third-party packages, could not be functionally tested without a real `derived_coefficients.tsv` from an actual Tier A run, but was reviewed manually line by line. The two R acquisition scripts and Tier A itself could not be syntax-checked (`Rscript -e "parse(...)"`) or run — this session's environment has no R, no `TCGAbiolinks`/`GEOquery`/`Seurat`/`WGCNA`/`glmnet` installed, and (unlike this repository's other modules, which at least have `bash -n` available for their shell scripts) this project has no shell scripts to check that way either. All three R scripts received careful manual line-by-line review only, matching the disclosure every module since Module 4 has made in this same environment — with the added caveat that this project has not been run end-to-end even once, by anyone, in this repository; Modules 2-8's Tier A scripts are at minimum pre-existing, previously-used material, while Scripts 1-2 here are newly authored specifically for this capstone instance.

## Interpretation **[FILL IN]**

**Before writing anything below, note what this scaffold does and does not establish:**

> Nothing in this project has been executed yet. The published paper's own findings (Section "Results" of the PDF) are not this project's results — they are the target this project's own real re-run should be checked against (Script 3), not a substitute for running it. Do not report the paper's numbers as if this repository produced them.

_Once Scripts 1→2→Tier A→3 have actually been run: does the re-derived signature agree in sign with the published one (Script 3's output)? Does the single-cell UMAP (Tier A Step 6) show a CD8+ T cell population consistent with the paper's own finding that CD8+ T cells were the most abundant tumor-infiltrating population? What would change your confidence in the signature — a sign disagreement on a specific gene, a very different WGCNA module structure, or something in the single-cell composition?_

## Reflection **[FILL IN]**

_What would you do differently starting this project again — e.g., start from Script 1's real data-acquisition plan from day one, rather than discovering partway through that the delivered `downstream.R` assumed data that had never been fetched? What was the actual value of finding and reading the source paper (`s41598-026-68769-z_reference.pdf`) rather than treating `downstream.R` as self-contained? What's the highest-value next step — GSE20685/GSE96058 validation, the full 51-gene ZH list, or resolving the GSE114725 format question?_

## Still to develop

- Run Scripts 1→2→Tier A→3 end-to-end on a machine with the real toolchain (the actual "Pilot run" and "Full execution" phases).
- Retrieve the real 51-gene zinc-homeostasis list from Guo et al. 2024 (see "Known, disclosed gaps" #1) and update Tier A's `zh_genes` accordingly — a disclosed correction, not a silent edit, once made.
- Resolve GSE114725's actual supplementary-file format (gap #2) and, if needed, replace Tier A's `Read10X()` call with the correct loader.
- Script the GSE20685/GSE96058 validation/external cohorts the paper itself uses to demonstrate generalizability.
- A pinned environment (renv lockfile or container definition) — currently just a package/version list in this README.
- Wet-lab validation (qRT-PCR, western blot), TIDE immunotherapy-response scoring, CellChat cell-cell communication, Monocle2 pseudotime, and oncoPredict drug-sensitivity prediction are all in the source paper's own full analysis but out of scope for this capstone instance — genuine further extensions, not gaps in what was attempted here.
