# Supplemental Case Studies — Applied Pathogen Pipelines

**Status: built.** These pipelines predate the formal 10-module curriculum document and don't map cleanly onto a single module — each one combines skills from several modules against a real pathogen dataset. They're kept as applied, cross-cutting case studies rather than forced into one module folder.

## Why these live outside `modules/00`–`09`

| Script | Combines skills from | Why it doesn't fit one module |
|---|---|---|
| [`scripts/Final_master_end_to_end_amr.sh`](scripts/Final_master_end_to_end_amr.sh) | [Module 1](../module-01-raw-reads-qc/README.md) (QC) + [Module 3](../module-03-genome-assembly/README.md) (assembly, via SPAdes) + AMR-specific gene annotation (AMRFinderPlus, not covered by any single module in the curriculum) | Antimicrobial resistance gene calling is its own applied domain, not a named curriculum module — it's a genuine capstone-style synthesis of assembly + annotation. |
| [`scripts/PRJNA1425489_mtb_single_end.sh`](scripts/PRJNA1425489_mtb_single_end.sh) | [Module 1](../module-01-raw-reads-qc/README.md) (QC) + [Module 2](../module-02-bulk-rnaseq/README.md) (RNA-seq quantification), single-end variant | A worked real-study example (*M. tuberculosis* RNA-seq) rather than a general-purpose teaching template — kept separate from the general Module 2 pipeline because it hardcodes single-end-specific logic that must **not** be merged with the paired-end assumption the Module 2 script relies on throughout. |

## What each script does

- **`Final_master_end_to_end_amr.sh`** — SRA download → assembly (SPAdes, `--isolate` mode) → QUAST quality assessment → AMRFinderPlus resistance-gene annotation, looped fault-tolerantly across a sample list, with a configurable `AMR_ORGANISM` variable (previously hardcoded to *Salmonella*). Aggregation only includes samples that finished successfully.
- **`PRJNA1425489_mtb_single_end.sh`** — Single-end RNA-seq quantification pipeline for the real PRJNA1425489 *M. tuberculosis* study: HISAT2 alignment with `-U` (single-end) → `featureCounts` **without** the paired-end-only `-p -B -C` flags (this was the critical bug fixed here — see the [changelog](../../docs/corrected_scripts_changelog.md), not yet present in this repo).

## Worked-example datasets (real, verified via NCBI eutils — same method Modules 3-7 used)

**`Final_master_end_to_end_amr.sh`** ships no dataset of its own — like most Tier A scripts in this repository, its only built-in default is the SRR11092056/SRR11092057 placeholder already flagged elsewhere in this repository (via NCBI eutils) as SARS-CoV-2 RNA-seq, not bacterial isolate WGS. [`scripts/smoke_test_amr_prjna242847.sh`](scripts/smoke_test_amr_prjna242847.sh) supplies a real one — and deliberately reuses **Module 4 and Module 5's own 5-isolate PRJNA242847 cohort** (GenomeTrakr/USDA-FSIS real-world *Salmonella enterica* subsp. *enterica* serovar Typhimurium surveillance WGS: `SRR40474618`, `SRR40472300`, `SRR40472280`, `SRR40426677`, `SRR40376221`) rather than sourcing a new dataset. This is deliberate, not laziness: it is the same organism `AMR_ORGANISM` defaults to, and it means this case study's AMR findings can be read directly alongside Module 4's variant calls and Module 5's tree for the **same five isolates**, not a disconnected fourth dataset.

**`PRJNA1425489_mtb_single_end.sh`** is unusual among this repository's Tier A scripts: its own `srr_list.txt` fallback is **already a real, correct, verified dataset**, not a placeholder needing replacement. Verified via NCBI eutils: all 12 hardcoded accessions (`SRR37277506`–`SRR37277517`) genuinely belong to BioProject PRJNA1425489, are single-end *Mycobacterium tuberculosis* RNA-seq (Illumina NovaSeq X Plus), and represent a real nitric-oxide + low-iron stress-response study (linked GEO series GSE319896, itself still under embargo until 2027-02-17 — the raw SRA reads are public even though the processed GEO record is not, which is expected and not a data-availability problem). **No smoke-test wrapper is provided for this script**, unlike every other Tier A pipeline in this repository — a deliberate choice, not an oversight, because there is no wrong default here to replace.

## How to run

Each script is self-contained; see its own header comment for the Slurm resource block and required `srr_list.txt` staging, following the same conventions as the numbered modules. See [`RUNBOOK.md`](RUNBOOK.md) for the operational step-by-step, required tools, expected outputs, and troubleshooting.

- **AMR pipeline:** `sbatch smoke_test_amr_prjna242847.sh` before ever submitting `Final_master_end_to_end_amr.sh` directly against a new cohort — SPAdes assembly is the single most expensive step in this entire repository's script collection (16 CPU / 64GB / 48h), and this is the only case study without Module 4/5's prior real-run history to lean on.
- **MTB pipeline:** `sbatch PRJNA1425489_mtb_single_end.sh` directly — its own built-in `srr_list.txt` default is already the real worked example (see above).

### Verification performed

Both Tier A scripts and the AMR smoke test pass `bash -n`. This session's environment has no `SPAdes`, `QUAST`, `ncbi-amrfinderplus`, `HISAT2`, `Subread`/`featureCounts`, or other bioinformatics toolchain installed, so neither script could be functionally run — they received careful manual line-by-line review only, matching the disclosure every module since Module 4 has made in this same environment.

## Suggested use in a course

Reference these from [Module 3](../module-03-genome-assembly/README.md) (as an AMR-flavored assembly example) and [Module 1](../module-01-raw-reads-qc/README.md)/[Module 2](../module-02-bulk-rnaseq/README.md) (as a single-end contrast to the paired-end default) once those modules have their own dedicated content, rather than duplicating logic into the numbered module folders. The AMR pipeline additionally cross-references cleanly with **Module 4** and **Module 5**, sharing the same 5-isolate cohort — a genuine three-way "same isolates, three analytical lenses (variants, phylogeny, resistance genes)" teaching opportunity worth calling out explicitly if this case study is ever assigned alongside those modules.

## Still to develop

- `docs/corrected_scripts_changelog.md`, linked above, does not yet exist in this repository (the same kind of pre-existing gap already disclosed in Module 2's README) — out of scope to create here.
- `amr_analysis.R` and a per-organism `srr_list_for_*_amr.sh`, both named in this repository's own domain coverage audit as part of AMR coverage alongside `Final_master_end_to_end_amr.sh`, were not found among this session's source materials and are not included here.
- Functional (not just syntax/manual-review) testing of both pipelines, on a machine with the real toolchain, following the Stage 0/Stage 1 pattern documented in this collection's own `HPC Validation Test Plan.pdf` (not itself committed to this repo).
