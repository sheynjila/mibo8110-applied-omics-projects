# mibo8110-applied-omics-projects

A competency-based, modular bioinformatics curriculum built for **MIBO8110**
— Slurm/HPC pipelines, teaching materials, and student-facing artifacts
organized around a 10-module applied-omics program (Modules 0–9), from raw
sequencing reads through capstone-level independent analysis.

## What this program is

Every module follows the same learning cycle: **frame a question → inspect
the data → choose an approach → execute it → evaluate quality → interpret
results → document limitations → make it reproducible.** The
[full curriculum document](docs/curriculum/modular_bioinformatics_curriculum.md)
is the canonical source for learning outcomes, portfolio artifacts, and the
assessment rubric — this repo is the applied, code-backed half of that
program: the actual scripts, pipelines, and per-module teaching materials
the curriculum calls for.

Not every offering of the course needs every module. Per the curriculum's
own recommendation: **Modules 0–1 are required for everyone**, then an
instructor selects a coherent core set from Modules 2–6, with Modules 7–8
as electives and Module 9 as the capstone synthesis.

## Module map — curriculum alignment and build status

| Module | Topic | Status | Folder |
|---|---|---|---|
| 0 | Computing, Data, Statistics, and Reproducibility | **Built** | [`modules/module-00-foundations-reproducibility/`](modules/module-00-foundations-reproducibility/README.md) |
| 1 | Raw Sequencing Reads and QC | **Built** | [`modules/module-01-raw-reads-qc/`](modules/module-01-raw-reads-qc/README.md) |
| 2 | Bulk RNA-seq and Gene Expression | **Built** | [`modules/module-02-bulk-rnaseq/`](modules/module-02-bulk-rnaseq/README.md) |
| 3 | Microbial Genome Assembly | **Built** | [`modules/module-03-genome-assembly/`](modules/module-03-genome-assembly/README.md) |
| 4 | Variant Calling | **Built** | [`modules/module-04-variant-calling/`](modules/module-04-variant-calling/README.md) |
| 5 | Phylogenomics and Outbreak Context | **Built** | [`modules/module-05-phylogenomics/`](modules/module-05-phylogenomics/README.md) |
| 6 | Metagenomics and Microbiome Profiles | **Built** | [`modules/module-06-metagenomics/`](modules/module-06-metagenomics/README.md) |
| 7 (elective) | Single-Cell Transcriptomics | **Built** (upstream quantification focus) | [`modules/module-07-single-cell-elective/`](modules/module-07-single-cell-elective/README.md) |
| 8 (elective) | Chromatin and Regulatory Genomics | **Built** | [`modules/module-08-chromatin-regulatory-elective/`](modules/module-08-chromatin-regulatory-elective/README.md) |
| 9 | Integrated Capstone and Professional Portfolio | **Framework built; 1 project instance** | [`modules/module-09-capstone/`](modules/module-09-capstone/README.md) |
| — | Supplemental applied case studies (AMR, *M. tuberculosis*) | **Built** | [`modules/supplemental-case-studies/`](modules/supplemental-case-studies/README.md) |
| — | Historical/exploratory scripts that predate the module structure | Kept for reference, not maintained | [`Initial_Pipeline_scripts/`](Initial_Pipeline_scripts/), [`Comprehensive_Bulk_RNA-seq/`](Comprehensive_Bulk_RNA-seq/) |

## Two ways to run every pipeline: manual and orchestrated

Every module's pipeline scripts are designed to be run **by hand first** —
`sbatch modules/module-NN-*/scripts/master_*.sh` (or the numbered
scripts, in order) against your own staged `srr_list.txt`, reading each
step as it runs. That is the primary teaching path, and it's what every
module's own `RUNBOOK.md` walks through.

For running the same pipeline **many times without collisions** — multiple
students, multiple datasets, a CI-style smoke test alongside a real
analysis — use the **orchestration copies** under
[`pipeline_orchestration/`](pipeline_orchestration/README.md). Each one is
a thin, `RUN_TAG`-parameterized copy of one master script; the manually-run
originals under `modules/*/scripts/` are never modified for this purpose.

```bash
# Manual — reads step-by-step, one fixed working directory
sbatch modules/module-04-variant-calling/scripts/master_snp_pipeline.sh

# Orchestrated — same pipeline, isolated working directory per run
RUN_TAG="section_A_$(date -u +%Y%m%dT%H%M%SZ)" \
  sbatch pipeline_orchestration/scripts/module-04-variant-calling/master_snp_pipeline.sh
```

See [`pipeline_orchestration/README.md`](pipeline_orchestration/README.md)
for the full list of orchestrated pipelines and how `RUN_TAG` isolation
works.

## Environment modules: `load_modules.sh`

Every pipeline script in this repo was authored and tested against exact
HPC environment-module builds (e.g. `FastQC/0.11.9-Java-11`). A bare
`module load <exact-build>` breaks the moment the same script runs on a
*different* cluster, because module trees are not standardized across
institutions. Instead:

- [`lib/module_loader.sh`](lib/module_loader.sh) is the shared engine: try
  the pinned build first; if it's missing, search `module avail
  <loose-name>` and load the single match with a loud warning; refuse to
  guess if there are zero or several matches.
- Each module's `scripts/load_modules.sh` declares that module's own
  tool → pinned-version map on top of the shared engine, and every
  pipeline script `source`s it instead of calling `module load` directly.
- **Override on your own cluster** with `<LABEL>_MODULE=<exact-name>`, e.g.
  `FASTQC_MODULE=FastQC/0.12.3-GCCcore-13.2.0 sbatch ...`.

[Module 0, Script 11](modules/module-00-foundations-reproducibility/scripts/11_Robust_HPC_Module_Loading.sh)
teaches this pattern directly, on this module's own five tools, before any
other module relies on it.

## Repository structure

```
mibo8110-applied-omics-projects/
├── README.md                              # this file — program-level overview
├── lib/
│   └── module_loader.sh                    # shared hybrid module-loading engine
├── pipeline_orchestration/                 # RUN_TAG-parameterized orchestration copies
│   ├── README.md
│   └── scripts/module-NN-*/...
├── docs/
│   └── curriculum/
│       └── modular_bioinformatics_curriculum.md  # canonical curriculum reference
└── modules/
    ├── module-00-foundations-reproducibility/    # scripts/ + lecture/ + README.md + RUNBOOK.md
    ├── module-01-raw-reads-qc/                    # scripts/ + docs/ + README.md + RUNBOOK.md
    ├── module-02-bulk-rnaseq/                      # scripts/ + lecture/ + README.md + RUNBOOK.md
    ├── module-03-genome-assembly/                  # scripts/ + docs/ + README.md + RUNBOOK.md
    ├── module-04-variant-calling/                  # scripts/ + docs/ + README.md + RUNBOOK.md
    ├── module-05-phylogenomics/                    # scripts/ + README.md + RUNBOOK.md
    ├── module-06-metagenomics/                     # scripts/ + README.md + RUNBOOK.md
    ├── module-07-single-cell-elective/             # scripts/ + README.md + RUNBOOK.md
    ├── module-08-chromatin-regulatory-elective/    # scripts/ + README.md + RUNBOOK.md
    ├── module-09-capstone/                         # README.md + per-project subfolders
    └── supplemental-case-studies/                  # scripts/ + README.md + RUNBOOK.md
```

## README architecture: high-level + per-module

This repo uses a **two-tier documentation structure**, mirroring how the
curriculum asks every module and every capstone project to be documented:

- **This root README** — the program-level view: what the whole curriculum
  is, the module map and build status, the folder layout, and shared
  conventions (module loading, orchestration) that apply everywhere.
- **`README.md` + `RUNBOOK.md` per module folder** — the README is the
  conceptual "why" (guiding question, learning outcomes, portfolio
  artifact); the `RUNBOOK.md` is the operational "how" (prerequisites,
  exact execution order, expected outputs, troubleshooting). Read a
  module's own README, then its RUNBOOK, before touching its scripts.

## Computing environment guidance (per curriculum §4)

- **Small/teaching datasets** (smoke tests, concept demonstrations): local WSL, a personal Linux/macOS workstation, or Galaxy ([Galaxy Training Network](https://training.galaxyproject.org/)) are all fine.
- **Full raw-read projects** (real cohort-scale runs): use the university HPC cluster or cloud compute — every script here targets Slurm.
- **Same structure, same provenance, same QC everywhere** — regardless of environment, the same acceptance criteria and documentation requirements apply.

## Conventions used throughout this repo

- **Fault-tolerant multi-sample loops.** Every multi-sample script wraps
  each sample's work in its own subshell so one failed sample is logged
  and skipped, not fatal to the whole Slurm job.
- **Robust module loading, never a bare `module load`.** See above.
- **Manual scripts stay manual; orchestration is a separate copy.** See
  [`pipeline_orchestration/README.md`](pipeline_orchestration/README.md).
- **One new idea per script/section.** Teaching progressions add exactly
  one concept at a time, with inline teaching-note comments explaining
  *why*, not just what changed.
- **Paired-end is the default assumption** for the RNA-seq pipeline
  family; single-end studies use a dedicated separate script (see
  [supplemental case studies](modules/supplemental-case-studies/README.md))
  rather than a shared flag, to avoid silently breaking the paired-end path.

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

## Reference frameworks cited in the curriculum

- [NIBLSE Core Competencies](https://qubeshub.org/community/groups/niblse/core_competencies)
- Wilson Sayres et al. 2018, [PLOS ONE](https://doi.org/10.1371/journal.pone.0196878)
- [Bioconductor DESeq2 vignette](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html)
- [Bioconductor edgeR documentation](https://bioconductor.org/packages/release/bioc/html/edgeR.html)
- [Galaxy Training Network](https://training.galaxyproject.org/)
- [Galaxy metagenomics learning pathway](https://training.galaxyproject.org/training-material/learning-pathways/metagenomics.html)
