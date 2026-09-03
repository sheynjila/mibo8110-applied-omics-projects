# Module 9 — Integrated Capstone and Professional Portfolio

**Status: one project instance built** ([`tcga-brca-zinc-homeostasis-signature/`](tcga-brca-zinc-homeostasis-signature/) — see that project's own README/RUNBOOK). This top-level README documents the capstone framework itself; each concrete project lives in its own `modules/module-09-capstone/<project-slug>/` subfolder with its own submission package.

**Curriculum module:** [Module 9 — Integrated Capstone and Professional Portfolio](../../docs/curriculum/modular_bioinformatics_curriculum.md).

## Purpose

The capstone is not a new analytical technique — it's the synthesis phase where a student combines whatever modules they've completed (Module 0-1 required, plus a chosen core/elective set) into one original, end-to-end, professionally documented analysis.

## The 8 phases (per curriculum)

1. Proposal
2. Data and compute plan
3. Analysis plan
4. Pilot run
5. Full execution
6. Quality-assurance audit
7. Interpretive report
8. Presentation / peer review

## Required submission package (per curriculum — applies to every capstone project)

| Component | Contents |
|---|---|
| README | Guiding question, data source, requirements, execution order, outputs |
| Environment | Pinned tool versions or container definition |
| Workflow | Scripts/notebooks with parameters and logging |
| Data provenance | Where the data came from, checksums, access date |
| Quality evidence | QC reports, acceptance criteria, what was checked |
| Results | The actual output artifacts |
| Interpretation | What the results mean, and their limitations |
| Reflection | What would be done differently, what was learned |

This is exactly the two-tier README pattern used across this whole repository (see the [root README](../../README.md)) — each capstone project's own README should read like this project's own module READMEs, scaled up to a full independent project.

## Assessment rubric (per curriculum)

| Criterion | Weight |
|---|---|
| Question / design | 15% |
| Data stewardship | 10% |
| Analytical reasoning | 20% |
| Quality evaluation | 15% |
| Reproducibility | 20% |
| Interpretation | 15% |
| Communication | 5% |

## Project instances

| Project | Draws on | Status |
|---|---|---|
| [`tcga-brca-zinc-homeostasis-signature/`](tcga-brca-zinc-homeostasis-signature/) | Module 2 (bulk RNA-seq/DESeq2) + Module 7 (single-cell QC/clustering) + survival modeling (LASSO Cox, Kaplan-Meier) not covered elsewhere in this repo | Scaffold built: real data-acquisition scripts (GDC/TCGAbiolinks, GEO/GEOquery), a pre-existing Tier A integrative analysis script, and a quality-assurance script comparing against the source publication's own reported coefficients. Not yet executed in any environment — see that project's own README "Verification performed." |

This instance came from real material found among this repository's own source folders: a delivered but non-runnable analysis script (`PRJNA482620_TCGA_BRCA_downstream.R`, which assumed pre-processed input files that were never actually produced by anything else in that source folder) plus its own uncited reference PDF, which — once read — turned out to be the real, citable, peer-reviewed publication the script was reproducing. See the project's own README for the full account, including what was corrected and what remains a disclosed, open gap.

## Suggested next step for additional instances

Per this module's own original scope: once a student (or instructor demo) selects a concrete capstone question, create `modules/module-09-capstone/<project-slug>/` with its own README following the submission-package table above, reusing scripts from whichever core modules the project draws on. Candidate material already identified elsewhere in this repository's source folders but not yet built: GSE20685/GSE96058 validation cohorts for the existing zinc-homeostasis project (a natural extension, not a new project); other combinations following the same "Module X + Module Y" synthesis pattern (e.g. Module 2 + Module 4 for an "expression changes correlate with a specific variant" style project, as the original draft README suggested).
