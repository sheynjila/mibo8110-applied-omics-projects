# Pipeline Orchestration

This directory holds **orchestration copies** of this course's master
pipeline scripts (`master_*.sh` and similar all-samples driver scripts under
each module's `scripts/`). It exists so that running many students, many
datasets, or many CI-style smoke tests through the same pipeline doesn't
collide on one shared working directory, **without touching the
manually-run scripts students actually work from.**

## Why a separate copy instead of editing the original in place

The scripts under `modules/*/scripts/` are the ones a student runs by hand
(`sbatch master_snp_pipeline.sh`, etc.) and reads while learning the
pipeline. Parameterizing them in place for concurrent/automated runs would
add orchestration plumbing (`RUN_TAG`, path indirection) to scripts whose
whole point is to be readable end-to-end during manual execution. Instead,
each script here is a thin copy of one original, changed in exactly two
ways:

1. **`RUN_TAG`** — every orchestration copy's `WORKDIR` is suffixed with
   `${RUN_TAG}` (default: a UTC timestamp) instead of being a single fixed
   path under `/scratch/$(whoami)/...`. Multiple runs -- different
   students, different datasets, a CI smoke test running alongside a real
   analysis -- each get their own isolated working directory.
2. **`PIPELINE_SCRIPT_DIR`** — set (before sourcing `load_modules.sh`) to
   point back at the *original* module's `scripts/` directory, so the
   orchestration copy reuses that module's real `load_modules.sh` rather
   than needing its own duplicate. See `lib/module_loader.sh`'s
   `resolve_script_dir()` for why this matters under real `sbatch` staging
   (BASH_SOURCE can resolve to a spool path that doesn't contain sibling
   files).

Everything else -- every tool call, every flag, every comment -- is
unchanged from the original. If you fix a bug in a module's master script,
apply the same fix here; there is intentionally no shared code path between
"the script a student reads" and "the script an orchestrator calls," so
neither can silently drift out of sync with the other without a diff
showing it.

## Layout

```
pipeline_orchestration/
`-- scripts/
    |-- module-01-raw-reads-qc/
    |   `-- 6_Automated_QC_Pipeline_AllSamples.sh
    |-- module-02-bulk-rnaseq/
    |   `-- master_rnaseq_pipeline_consolidated.sh
    |-- module-04-variant-calling/
    |   `-- master_snp_pipeline.sh
    |-- module-06-metagenomics/
    |   |-- master_amplicon_pipeline.sh
    |   `-- master_microbiome_pipeline.sh
    |-- module-07-single-cell-elective/
    |   `-- master_scrnaseq_pipeline.sh
    |-- module-08-chromatin-regulatory-elective/
    |   `-- master_chipseq_pipeline.sh
    `-- supplemental-case-studies/
        |-- Final_master_end_to_end_amr.sh
        `-- PRJNA1425489_mtb_single_end.sh
```

## Running an orchestrated job

```bash
RUN_TAG="cohort_A_$(date -u +%Y%m%dT%H%M%SZ)" \
  sbatch pipeline_orchestration/scripts/module-04-variant-calling/master_snp_pipeline.sh
```

Then stage `srr_list.txt` (and any other per-run input, e.g. `HOST_REFERENCE`
for module 6) inside the resulting `WORKDIR` before the job needs it, exactly
as you would for the manual script -- the only difference is the directory
now has `_${RUN_TAG}` on the end, so it won't collide with anyone else's run
or with your own previous run.

Omit `RUN_TAG` and one is generated for you from the current UTC timestamp,
so a bare `sbatch pipeline_orchestration/scripts/.../master_*.sh` still
never collides with a prior run.

## What this is NOT

This is not a workflow-manager (Snakemake/Nextflow) integration, and it does
not schedule dependent jobs for you. It solves exactly one problem --
"multiple runs of the same master script must not write into the same
WORKDIR" -- so this course can support both "a student runs a pipeline
manually and reads every step" and "many runs happen at once" without
maintaining two different pipeline implementations.
