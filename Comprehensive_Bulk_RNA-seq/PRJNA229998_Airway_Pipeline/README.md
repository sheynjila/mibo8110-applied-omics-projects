# PRJNA229998 Airway RNA-Seq Pipeline

A reproducible, Slurm-native bulk RNA-seq pipeline that takes raw SRA reads from the
Himes et al. (2014) airway smooth-muscle glucocorticoid-response study (BioProject
[PRJNA229998](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA229998), GEO
[GSE52778](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE52778)) from
`prefetch` through to a `featureCounts` gene-level count matrix ready for `DESeq2`.

This repository is teaching material for the **MIBO8110 Applied Omics** curriculum. It
is deliberately built and documented as a *worked, debugged example* of a real HPC
RNA-seq pipeline: every fix listed below was driven by an actual failure encountered
on a Slurm cluster, not written speculatively. The full failure-and-fix history is in
[`TROUBLESHOOTING_REPORT.txt`](TROUBLESHOOTING_REPORT.txt) and
[`Detailed_Pipeline_Troubleshooting_Report.pdf`](Detailed_Pipeline_Troubleshooting_Report.pdf),
and is distilled into the [Runbook](#runbook) and the
[Known Limitations](#known-limitations-read-before-production-use) section below.

> **Status:** Pedagogical / reference pipeline. The shell logic has been reviewed
> line-by-line and validated with `bash -n` syntax checks; it has **not** been executed
> end-to-end in the environment this documentation was written in (no HPC/Slurm access
> from this workstation). Treat first production runs as a validation run, not a
> guaranteed-clean execution — see [Known Limitations](#known-limitations-read-before-production-use).

## Table of Contents

- [Study background](#study-background)
- [Repository layout](#repository-layout)
- [Pipeline architecture](#pipeline-architecture)
- [Prerequisites](#prerequisites)
- [Quickstart](#quickstart)
- [The modular suite, step by step](#the-modular-suite-step-by-step)
- [The automated cohort run (Script 8)](#the-automated-cohort-run-script-8)
- [Monitoring & troubleshooting](#monitoring--troubleshooting)
- [The "V9" HPC safeguards](#the-v9-hpc-safeguards)
- [Known limitations (read before production use)](#known-limitations-read-before-production-use)
- [Downstream handoff to DESeq2](#downstream-handoff-to-deseq2)
- [Runbook](#runbook)
- [Citation](#citation)
- [License](#license)

## Study background

| | |
|---|---|
| **Study** | Himes BE, Jiang X, Wagner P, et al. (2014). *RNA-Seq Transcriptome Profiling Identifies CRISPLD2 as a Glucocorticoid Responsive Gene that Modulates Cytokine Function in Airway Smooth Muscle Cells.* PLOS ONE. |
| **DOI** | [10.1371/journal.pone.0099625](https://doi.org/10.1371/journal.pone.0099625) |
| **GEO series** | [GSE52778](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE52778) |
| **BioProject** | [PRJNA229998](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA229998) |
| **Design** | Four human airway smooth muscle cell lines, treated with dexamethasone (a glucocorticoid) vs. untreated control — a paired design (donor + treatment) |
| **Samples used in this pipeline** | `SRR1039508`, `SRR1039509`, `SRR1039512`, `SRR1039513` (one donor pair's worth of paired-end reads; the full GSE52778 series has 8 samples / 4 donors) |
| **Reference genome** | Ensembl GRCh38, release 110 (`Homo_sapiens.GRCh38.dna.primary_assembly.fa`, `Homo_sapiens.GRCh38.110.gtf`) |

This dataset is also the canonical worked example for the Bioconductor `DESeq2` and
`airway` packages, which makes it a useful teaching dataset: students can process the
raw reads through this pipeline and then sanity-check their own count matrix against
the pre-built `airway` Bioconductor package's counts. `SRR1039513` in particular is a
useful edge case for this course — see
[Known Limitations](#known-limitations-read-before-production-use).

For a broader set of teaching datasets this pipeline's design generalizes to (batch
effects, dual RNA-seq, viral-load covariates, etc.), see
[`Host_Pathogen_Bulk_RNAseq_Guide-v2.pdf`](Host_Pathogen_Bulk_RNAseq_Guide-v2.pdf).

## Repository layout

```text
PRJNA229998_Airway_Pipeline/
├── 1_Fetch_And_Verify_Raw_Reads.sh          # Modular suite, step 1
├── 2_Interpret_FastQC_Diagnostics.sh        # Modular suite, step 2
├── 3_Batch_MultiQC_Summary.sh               # Modular suite, step 3
├── 4_Evidence_Based_Fastp_Trimming.sh       # Modular suite, step 4
├── 5_Paired_End_Integrity_Check.sh          # Modular suite, step 5
├── 6_STAR_Alignment.sh                      # Modular suite, step 6
├── 7_Quantification_Rsubread.sh             # Modular suite, step 7
├── 8_Automated_Pipeline_AllSamples.sh       # V9 production cohort run (all 4 samples)
├── 9_QC_Report_Generator.sh                 # Assembles FINAL_PIPELINE_REPORT.md
├── PRJNA229998_airway_pipeline.sh           # Historical monolithic pipeline (v6-era) — reference only, see note below
├── Final_airway_pipeline_v9_standalone.sh   # V9 validation prototype — reference only, see note below
├── README.md                                # This file
├── RNA_Seq_Master_Runbook.docx / .pdf       # Operational SOP runbook
├── TROUBLESHOOTING_REPORT.txt               # v1 → v6 failure/fix narrative
├── Detailed_Pipeline_Troubleshooting_Report.pdf  # v7 → v9 failure/fix narrative
└── Host_Pathogen_Bulk_RNAseq_Guide-v2.pdf   # Broader teaching-dataset guide (context, not executed by this pipeline)
```

**A note on the two "extra" full-pipeline `.sh` files.** `PRJNA229998_airway_pipeline.sh`
and `Final_airway_pipeline_v9_standalone.sh` are earlier, single-file monolithic
drafts of this pipeline, kept in the repository as historical/reference artifacts
documenting how the numbered, modular suite (scripts 1–9) came to be. **They are not
the recommended way to run this pipeline** — see
[Known Limitations](#known-limitations-read-before-production-use) for a specific,
verified inconsistency in `PRJNA229998_airway_pipeline.sh`. Scripts **1–9** are the
current, maintained pipeline.

## Pipeline architecture

The repository has a two-tier design:

1. **Modular suite (scripts 1, 2, 3, 4, 5, 6, 7, 9).** Runs the pipeline on a single
   worked example sample (`SRR1039508`) one Slurm job at a time, in
   `/scratch/$(whoami)/PRJNA229998_airway_pipeline_stepbystep/`. Built for teaching,
   debugging, and validating each stage's output before trusting the automated run.
2. **V9 automated cohort suite (script 8).** A single Slurm job that loops over all
   four samples in a fault-tolerant subshell, in
   `/scratch/$(whoami)/PRJNA229998_airway_pipeline_automated/`, and ends with one
   combined `featureCounts` matrix across every sample that made it through the
   loop. This is the "production" entry point once the modular walkthrough has been
   validated once against the reference genome.

Both tiers share one reference directory,
`/scratch/$(whoami)/PRJNA229998_airway_pipeline_ref/`, containing the downloaded
Ensembl FASTA/GTF and the built STAR index — **built once by Script 6**, then reused
by Script 8 (see the [dependency note](#the-automated-cohort-run-script-8) below).

```text
/scratch/$(whoami)/
├── PRJNA229998_airway_pipeline_ref/          # Shared: Ensembl GRCh38 FASTA/GTF + STAR index (built by Script 6)
├── PRJNA229998_airway_pipeline_stepbystep/   # Modular suite workspace (SRR1039508 only)
│   ├── raw_reads/  checksums/  qc_before/  trimmed_reads/  fastp_reports/
│   ├── qc_after/   alignments/  counts/     logs/           sra_cache/
│   └── FINAL_PIPELINE_REPORT.md              # Written by Script 9
└── PRJNA229998_airway_pipeline_automated/    # V9 cohort workspace (all 4 samples)
    ├── raw_reads/  checksums/  qc_before/  trimmed_reads/  fastp_reports/
    ├── pairing_reports/  alignments/  logs/  sra_cache/
    └── counts/
        ├── airway_raw_counts.txt             # Final gene x sample count matrix
        └── airway_raw_counts.txt.summary
```

## Prerequisites

- A Slurm-managed HPC cluster with internet egress from compute nodes (for
  `prefetch`/SRA and the Ensembl FTP downloads in Script 6).
- Lmod-style environment modules providing:

  | Tool | Module used in these scripts |
  |---|---|
  | SRA Toolkit ≥ 3.0 | `SRA-Toolkit/3.0.3-gompi-2022a` |
  | FastQC | `FastQC/0.11.9-Java-11` |
  | fastp | `fastp/0.23.4-GCC-13.2.0` |
  | STAR | `STAR/2.7.10b-GCC-11.3.0` |
  | MultiQC | `MultiQC/1.28-foss-2024a` |
  | R (with Bioconductor access) | `R` (cluster default) |

  Module names/versions are cluster-specific — confirm with `module avail <tool>`
  before your first run and edit the `module load` lines if your cluster's naming
  differs (see [Known Issues → Phase 2](#the-v9-hpc-safeguards)).
- **Disk/scratch space:** the GRCh38 primary-assembly FASTA + GTF + STAR index
  together need roughly 35–40 GB; budget several more GB per sample for raw/trimmed
  FASTQs and BAMs. `SRR1039513` alone extracts to ~67 million lines per mate
  (≈17M read pairs) — see the note on its zero-length reads below.
- **Memory:** STAR's `genomeGenerate` step needs ~32 GB resident memory for GRCh38 —
  never run it on a shared login node (see
  [Known Issues → Phase 3](#the-v9-hpc-safeguards)). All alignment/indexing jobs in
  this repo request `--mem=36G`.

## Quickstart

```bash
git clone <this-repo-url>
cd PRJNA229998_Airway_Pipeline

# Recommended first run: validate the pipeline end-to-end on ONE sample
sbatch 1_Fetch_And_Verify_Raw_Reads.sh
# wait for completion (squeue -u $(whoami)), then:
sbatch 2_Interpret_FastQC_Diagnostics.sh
sbatch 3_Batch_MultiQC_Summary.sh
sbatch 4_Evidence_Based_Fastp_Trimming.sh
sbatch 5_Paired_End_Integrity_Check.sh
sbatch 6_STAR_Alignment.sh          # also builds & caches the shared STAR index
sbatch 7_Quantification_Rsubread.sh
bash   9_QC_Report_Generator.sh     # local job, writes FINAL_PIPELINE_REPORT.md

# Once the modular run above has completed successfully at least once (so the
# shared STAR index in PRJNA229998_airway_pipeline_ref/ exists):
sbatch 8_Automated_Pipeline_AllSamples.sh
```

Each `sbatch` step must finish before submitting the next — the modular suite has no
internal dependency chaining (`--dependency=afterok:<jobid>`), by design, so that each
stage's output can be inspected before continuing. Full step-by-step verification
checkpoints are in the [Runbook](#runbook).

## The modular suite, step by step

| # | Script | Purpose | Key resource request |
|---|---|---|---|
| 1 | `1_Fetch_And_Verify_Raw_Reads.sh` | `prefetch` → `vdb-validate` (cache-corruption gate) → `fasterq-dump --split-3` → `md5sum` baseline | 4 CPU, 8 GB, 1 h |
| 2 | `2_Interpret_FastQC_Diagnostics.sh` | Raw-read FastQC with `--extract`, parses `summary.txt` into PASS/WARN/FAIL counts | 4 CPU, 6 GB, 45 min |
| 3 | `3_Batch_MultiQC_Summary.sh` | Aggregates `qc_before/` into one MultiQC HTML dashboard | 4 CPU, 8 GB, 1.5 h |
| 4 | `4_Evidence_Based_Fastp_Trimming.sh` | Re-reads Script 2's FastQC verdicts and *conditionally* enables adapter trimming and/or 3′ quality trimming in `fastp`, writes a rationale file, then FastQC's the cleaned reads | 4 CPU, 8 GB, 1 h |
| 5 | `5_Paired_End_Integrity_Check.sh` | Confirms mate 1 / mate 2 read counts and read IDs are still 1:1 synchronized after trimming | 2 CPU, 4 GB, 30 min |
| 6 | `6_STAR_Alignment.sh` | Downloads Ensembl GRCh38 FASTA/GTF (if absent), builds the shared STAR index (if absent), aligns the sample to a coordinate-sorted BAM | 8 CPU, 36 GB, 2 h |
| 7 | `7_Quantification_Rsubread.sh` | Gene-level counting via `Rsubread::featureCounts`, called from a dynamically generated R script; bootstraps `Rsubread` into a local `R_libs` if the module lacks it | 4 CPU, 16 GB, 1 h |
| 9 | `9_QC_Report_Generator.sh` | Assembles `FINAL_PIPELINE_REPORT.md` from the outputs above (raw/trimmed read counts, pairing verdict, STAR mapping rate, `featureCounts` assigned reads) | 1 CPU, 4 GB, 10 min — run locally with `bash`, not `sbatch` |

Script 4's evidence-based trimming is the pedagogical core of this pipeline: it does
not apply adapter or quality trimming unconditionally. It greps Script 2's FastQC
`summary.txt` for `WARN`/`FAIL` verdicts on **Adapter Content** and **Per base
sequence quality**, and only appends `--detect_adapter_for_pe` and/or
`--cut_tail --cut_tail_mean_quality 20` to the `fastp` command line if those specific
verdicts fired — then writes the decision, and why it was made, to
`fastp_reports/<SRR>_THRESHOLD_RATIONALE.md`.

## The automated cohort run (Script 8)

`8_Automated_Pipeline_AllSamples.sh` re-implements steps 1, 2, 4, 5, and 6 (fetch →
FastQC → evidence-based fastp → pairing check → STAR align) as one Slurm job, looped
over all four samples, each iteration wrapped in `( set -e; set -o pipefail; ... )` so
a single sample's failure is meant to be recorded rather than crashing the batch (see
the caveat about this in [Known Limitations](#known-limitations-read-before-production-use)).
After the loop, it runs one batch `featureCounts` call across every sample that
reported success, writing `counts/airway_raw_counts.txt`.

**Dependency you must satisfy yourself:** Script 8 references
`${REFDIR}/star_index` directly — it does **not** build the STAR index itself. Run
Script 6 at least once first (on `SRR1039508`, in the modular workspace) so the
shared index at `PRJNA229998_airway_pipeline_ref/star_index/` exists before
submitting Script 8.

Script 8 also does **not** call MultiQC or generate a report file — Script 9 only
summarizes the single-sample modular workspace. If you want a human-readable QC
summary of the 4-sample cohort run, review `logs/<SRR>_pipeline.log` per sample and
`counts/airway_raw_counts.txt.summary` directly, or extend Script 9 to point at the
`_automated` workspace.

## Monitoring & troubleshooting

```bash
# Queue status
squeue -u $(whoami)

# Cancel a job
scancel <JOB_ID>

# Tail a specific step's custom log
tail -f /scratch/$(whoami)/PRJNA229998_airway_pipeline_stepbystep/logs/SRR1039508_fetch.log

# Tail STAR's own progress log during alignment
tail -f /scratch/$(whoami)/PRJNA229998_airway_pipeline_stepbystep/alignments/SRR1039508_Log.progress.out

# Sweep Slurm stdout for common failure signatures
grep -iE "error|fatal|cannot|execution halted" ~/slurm-*.out

# Confirm STAR's mapping rate
tail -n 20 alignments/*_Log.final.out

# Confirm featureCounts assignment rate
cat counts/*_counts.txt.summary
```

For a deeper failure/fix narrative — including the exact `srun` incantations needed
to get an interactive 36 GB compute node allocation on a cluster that rejects
under-specified `srun` calls — see
[`Detailed_Pipeline_Troubleshooting_Report.pdf`](Detailed_Pipeline_Troubleshooting_Report.pdf)
and the version-by-version history in
[`TROUBLESHOOTING_REPORT.txt`](TROUBLESHOOTING_REPORT.txt).

## The "V9" HPC safeguards

Three stability patches, driven by real failures during development (documented in
full in the two troubleshooting reports), are baked into every script in the
1–9 suite:

- **`--mem=36G` on every STAR job.** Building the GRCh38 STAR index needs roughly
  ~31.5 GB resident memory; running it on a shared login node throws
  `std::bad_alloc` and gets killed by the node's own memory limits, not by Slurm.
- **`fasterq-dump --split-3`, never `--split-files`.** `SRR1039513` alone contains
  over 26.5 million zero-length reads. `--split-files` drops these asynchronously
  between mates, desynchronizing `_1.fastq`/`_2.fastq` line counts and causing
  `fastp` to fail immediately with `input files don't contain identical amount of
  reads`. `--split-3` instead routes unmated/zero-length reads to a separate
  singletons file, preserving 1:1 mate parity in the paired files.
- **`module purge` before every `module load` block.** Loading `fastp`'s GCC/13
  toolchain and STAR's GCC/11 toolchain in the same environment without purging
  between them silently drops shared C/C++ libraries (Lmod dependency clashes),
  producing cryptic runtime failures rather than a clean module-load error.

Script 7's dynamic Rsubread bootstrap (installing `Rsubread` into a local
`R_libs` scratch directory via `BiocManager` if the cluster's `R` module doesn't
already have it) is a fourth safeguard, aimed at clusters where users can't get
admin approval to add packages to the shared R library.

## Known limitations (read before production use)

These were found during a line-by-line review of the scripts against their own
documented intent, not during an actual run (no HPC access from this authoring
environment — see the [Status note](#prjna229998-airway-rna-seq-pipeline) at the
top of this file). They are disclosed here rather than silently patched, so a
future run and a future maintainer both know exactly what to check.

1. **Script 8's per-sample fault tolerance does not currently work as documented.**
   Each sample is processed inside `( set -e; set -o pipefail; ... )`, and the
   *outer* script also runs under `set -e` (line 14). Because the subshell is
   invoked as a plain top-level command — not inside an `if`/`&&`/`||` — a
   non-zero exit from any single sample's subshell will trigger the outer
   script's `set -e` and terminate the entire Slurm job immediately, before the
   `if [ $? -eq 0 ]` line that is supposed to route the failure into
   `FAILED_SRRS` ever runs. In practice this means one bad sample can still
   abort the whole cohort run, which is exactly the failure mode the subshell
   wrapping was designed to prevent. The same pattern exists in the legacy
   `PRJNA229998_airway_pipeline.sh`. **Workaround until fixed:** run the loop
   body's commands with `set +e` reinstated around the per-sample dispatch, or
   capture the subshell's exit status via `... ; RC=$?` inside a context `set -e`
   won't fire on (e.g., `if ( ... ); then …`), before relying on this script to
   silently skip a bad sample in an unattended run.
2. **Script 8's final quantification step does not bootstrap Rsubread.** Unlike
   Script 7, the R script Script 8 writes to `counts/run_final_quant.R` calls
   `suppressMessages(library(Rsubread))` directly, with no `BiocManager`/install
   fallback. If the cluster's default `R` module doesn't already have `Rsubread`
   installed, this fails only *after* all four samples have already been fetched,
   trimmed, and aligned — the most expensive part of the run. Confirm
   `Rscript -e 'library(Rsubread)'` succeeds under the `R` module you load before
   relying on Script 8 for a production cohort run, or point
   `R_LIBS_USER`/`.libPaths()` at the `R_libs` directory Script 7 already
   populated in the sibling `_stepbystep` workspace.
3. **`PRJNA229998_airway_pipeline.sh` still uses `fasterq-dump --split-files`**
   (line 125) despite carrying an appended self-assessment block that claims
   "the `--split-3` flag is hardcoded." That assessment text describes
   `Final_airway_pipeline_v9_standalone.sh`, not this file — the two were
   evidently not kept in sync. This is the specific, concrete reason
   `PRJNA229998_airway_pipeline.sh` is documented above as historical/reference
   only, not a recommended entry point: run against `SRR1039513` as-is, it would
   reproduce the exact read-desynchronization failure described in
   [`Detailed_Pipeline_Troubleshooting_Report.pdf`](Detailed_Pipeline_Troubleshooting_Report.pdf).
4. **`--sjdbOverhang` is inconsistent across the three pipeline variants in this
   repo:** `74` in `6_STAR_Alignment.sh` (the current modular suite), `99` in
   `Final_airway_pipeline_v9_standalone.sh`, `100` in
   `PRJNA229998_airway_pipeline.sh`. STAR's own guidance is `read length − 1`.
   Confirm the actual read length of your FASTQ files (`awk 'NR==2{print
   length($0); exit}' raw_reads/<SRR>_1.fastq`) before a from-scratch index
   build, rather than assuming `6_STAR_Alignment.sh`'s `74` is correct for your
   data — it is the most recently maintained value in this repo, not an
   independently re-derived one.
5. **Script 9 (`FINAL_PIPELINE_REPORT.md`) only ever summarizes `SRR1039508`** in
   the `_stepbystep` workspace. It is not wired up to the 4-sample
   `_automated` workspace Script 8 writes to — see
   [The automated cohort run](#the-automated-cohort-run-script-8) above.

None of the above are silently corrected in this repository's scripts — this README
and the Runbook are the disclosure. If you fix one, please also update this section
and cite the fix in the Runbook's version-evolution appendix.

## Downstream handoff to DESeq2

The pipeline's terminal output is a tab-delimited gene × sample count matrix:

- Modular suite (single sample): `counts/SRR1039508_counts.txt`
- Automated cohort: `counts/airway_raw_counts.txt` (rows = Ensembl gene IDs, columns
  = sample SRR accessions)

This is import-ready for a standard `DESeqDataSetFromMatrix()` call in R. Because
`PRJNA229998`/`GSE52778` is the canonical `DESeq2` and `airway`-package worked
example, you can build the accompanying `colData` design matrix (`cell_line`,
`dex` treatment/untreated) directly from the study's public SRA RunInfo/metadata, and
sanity-check your own recount against the pre-built Bioconductor `airway` object as
part of the course exercise.

## Runbook

[`RNA_Seq_Master_Runbook.docx`](RNA_Seq_Master_Runbook.docx) /
[`RNA_Seq_Master_Runbook.pdf`](RNA_Seq_Master_Runbook.pdf) is the operational
companion to this README: exact `sbatch`/`srun` commands in execution order, a
verification checklist after each step, the full troubleshooting decision history
(v1 → v9), a resource-allocation appendix, and the same known-limitations list above
in more procedural form (what to check, and when, rather than just what's wrong).

## Citation

If you use this pipeline or its worked example in coursework or derivative material,
please cite the underlying study:

> Himes BE, Jiang X, Wagner P, Hu R, Wang Q, Klanderman B, et al. (2014) RNA-Seq
> Transcriptome Profiling Identifies CRISPLD2 as a Glucocorticoid Responsive Gene
> that Modulates Cytokine Function in Airway Smooth Muscle Cells. *PLoS ONE* 9(6):
> e99625. https://doi.org/10.1371/journal.pone.0099625

## License

No license file is currently included in this repository. Before publishing this
repository publicly on GitHub, add an OSI-approved license (e.g., MIT for the
pipeline code) appropriate for your institution's and course's policies — the
underlying sequencing data and publication remain governed by their own SRA/GEO and
journal terms regardless of the license chosen for this code.

---

*Maintained as part of the MIBO8110 Applied Omics curriculum. See
[`Comprehensive_Bulk_RNA-seq`](..) for related bulk RNA-seq teaching material.*
