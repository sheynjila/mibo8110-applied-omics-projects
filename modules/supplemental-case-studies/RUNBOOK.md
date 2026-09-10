# Supplemental Case Studies — Applied Pathogen Pipelines: Runbook

This is an **operational** document — step-by-step execution instructions,
expected outputs, and troubleshooting. For why these two pipelines live
outside the numbered modules, see this folder's [`README.md`](README.md).

## Prerequisites

**AMR pipeline (`Final_master_end_to_end_amr.sh`):**
- SRA-Toolkit, fastp, SPAdes, QUAST, ncbi-amrfinderplus, MultiQC — loaded
  automatically via [`scripts/load_modules.sh`](scripts/load_modules.sh)
  (pinned to SRA-Toolkit 3.2.0, fastp 0.23.4, SPAdes 3.15.5, QUAST 5.2.0,
  ncbi-amrfinderplus 3.11.11, MultiQC 1.28, with an automatic fallback if
  your cluster doesn't have those exact builds — see the root
  [`RUNBOOK.md`](../../RUNBOOK.md) §4).
- A network-reachable AMRFinderPlus database — Phase 1 runs `amrfinder -u`
  to update it before scanning; this requires network access from the
  compute node, not just the login node, on some HPC systems.
- Substantial compute: 16 CPU / 64GB memory / 48h wall time (SPAdes
  `--isolate` mode assembly is the most resource-intensive step in this
  entire repository).

**MTB pipeline (`PRJNA1425489_mtb_single_end.sh`):**
- SRA-Toolkit, FastQC, fastp, HISAT2, SAMtools, Subread (for
  `featureCounts`), MultiQC — same
  [`scripts/load_modules.sh`](scripts/load_modules.sh), pinned to
  SRA-Toolkit 3.0.3 for this script specifically (note: the AMR pipeline
  above pins 3.2.0 — `load_sra_toolkit` takes the version as an argument
  precisely because these two scripts were authored against different
  builds), FastQC 0.11.9, fastp 0.23.4, HISAT2 2.2.1, SAMtools 1.18,
  Subread 2.0.6, MultiQC 1.28.
- No external database — the *M. tuberculosis* H37Rv reference
  (GCF_000195955.2, ASM19595v2) is downloaded directly from NCBI by Phase 2.

## Execution order — AMR pipeline

| Step | Command | Expected result |
|---|---|---|
| 0 | `sbatch smoke_test_amr_prjna242847.sh` | `$WORKDIR/amr_reports/ALL_SAMPLES_AMR_SUMMARY.tsv`, `$WORKDIR/quast_qc/<SRR>/`, `$WORKDIR/final_assembly_report.html` (MultiQC) |
| 1 | Inspect `quast_qc/<SRR>/report.txt` for each sample — total length should be ~4.8-5.0 Mb for *Salmonella*; ~9.5 Mb strongly suggests contamination/mixed-species (see Tier A's own footer commentary) | Confirmed clean assemblies before trusting the AMR calls built on them |
| 2 | Inspect `amr_reports/ALL_SAMPLES_AMR_SUMMARY.tsv` | Per-sample resistance gene calls, aggregated |

`$WORKDIR` defaults to `/scratch/$(whoami)/end_to_end_amr_pipeline` (set
inside `Final_master_end_to_end_amr.sh`). For a different organism, set
`AMR_ORGANISM` (e.g. `export AMR_ORGANISM=Escherichia`) before running
against a real `srr_list.txt` for that organism — do not reuse the
Salmonella smoke-test cohort with a different `AMR_ORGANISM`.

## Execution order — MTB pipeline

| Step | Command | Expected result |
|---|---|---|
| 0 | `sbatch PRJNA1425489_mtb_single_end.sh` (no smoke test needed — see README) | `$WORKDIR/counts/gene_counts_clean.tsv`, `$WORKDIR/final_multiqc_report.html` |
| 1 | Inspect `counts/gene_counts_clean.tsv` | A gene-by-sample raw count matrix (12 columns, one per SRR), ready for `rnaseq_deseq2_analysis.R` (Module 2) — confirm 12 sample columns are present, matching `srr_list.txt`'s 12 accessions, before assuming every sample succeeded |

`$WORKDIR` defaults to `/scratch/$(whoami)/PRJNA1425489_rnaseq_SE`.

## Troubleshooting

| Symptom | Likely cause | Fix |
|---|---|---|
| AMR Step 0 fails at `amrfinder -u` | Compute node lacks network access to the AMRFinderPlus database update server | Check your HPC's outbound-network policy for compute nodes; some clusters require this update step to run on a login/data-transfer node instead |
| AMR Step 0 reports an assembly around 9.5 Mb instead of ~4.8-5.0 Mb | Contamination or a mixed-species sample (Tier A's own footer commentary) | Do not proceed to AMRFinderPlus on that sample's contigs; investigate the raw reads (e.g. Kraken2, this repository's Module 6 shotgun pipeline) before re-assembling |
| AMR pipeline reports `FAILED_SAMPLES` for one isolate | A failed download, assembly, QUAST run, or AMR scan for that specific sample (Tier A's own fault-tolerant loop logs and continues) | Check that sample's console output; re-run just that SRR or accept the smaller cohort |
| MTB pipeline's `featureCounts` step errors with "no paired-end reads detected" | You are running the ORIGINAL (uncorrected) version of this script, not the corrected one in this folder, or you accidentally re-added `-p -B -C` | Confirm you are running `scripts/PRJNA1425489_mtb_single_end.sh` from this repository, unmodified |
| MTB pipeline's `gene_counts_clean.tsv` has fewer than 12 sample columns | One or more of the 12 SRR accessions failed in the fault-tolerant loop | Check console output for `FAILED_SAMPLES`; re-run those specific accessions individually |
| MTB pipeline's reference download step reports a "downloaded empty/corrupt" error | A truncated `wget`/`gunzip` on a prior run left a 0-byte or partial file, OR a genuinely failed download this run | The script's own `-s` (non-empty) check (Fix #3, see the script's header) already re-fetches automatically on the next run — re-run the script; if it fails again, check network access to `ftp.ncbi.nlm.nih.gov` directly |

## Cross-referencing with the numbered modules

If assigning the AMR case study alongside Modules 4-5 (same PRJNA242847
cohort): `amr_reports/ALL_SAMPLES_AMR_SUMMARY.tsv`'s resistance-gene calls,
Module 4's `outbreak_clusters.tsv` (SNP-distance clustering), and Module 5's
`cohort_tree.treefile` (phylogeny) all describe the same five isolates from
three different analytical angles — a real opportunity to ask whether
isolates that cluster together genomically also share the same resistance
gene profile, without fabricating a connection between the three outputs
that isn't actually there.
