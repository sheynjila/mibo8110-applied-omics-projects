# Module 1 — Raw Sequencing Reads & QC: Runbook

This is an **operational** document — step-by-step execution instructions,
expected outputs, and troubleshooting. For the conceptual "why," see this
module's [`README.md`](README.md). For the full curriculum context, see
[`docs/curriculum/modular_bioinformatics_curriculum.md`](../../docs/curriculum/modular_bioinformatics_curriculum.md).

## Prerequisites

- Access to an HPC cluster with Slurm (`sbatch` on `PATH`), or an equivalent
  environment.
- SRA-Toolkit, FastQC, fastp, and MultiQC — loaded automatically via
  [`scripts/load_modules.sh`](scripts/load_modules.sh) (pinned to
  `SRA-Toolkit/3.0.3-gompi-2022a`, `FastQC/0.11.9-Java-11`,
  `fastp/0.23.4-GCC-13.2.0`, `MultiQC/1.28-foss-2024a`, with an automatic
  fallback if your cluster doesn't have those exact builds — see the root
  [`RUNBOOK.md`](../../RUNBOOK.md) §4).
- A scratch/working directory with at least ~15GB free (PRJNA1518998's two
  reference SRR runs are a few GB each; raw + trimmed FASTQ + QC reports for
  both roughly doubles that; leave headroom).
- Module 0 completed or at least reviewed — this module assumes the same
  storage-safety and fault-tolerant-loop habits taught there.

## Execution order

Run every step from the module's `scripts/` directory unless noted.

| Step | Command | Where | Expected result |
|---|---|---|---|
| 1 | `sbatch 1_Fetch_And_Verify_Raw_Reads.sh` | compute | `raw_reads/SRR40359383_{1,2}.fastq` and `checksums/SRR40359383_raw.md5` created |
| 2 | `sbatch 2_Interpret_FastQC_Diagnostics.sh` | compute | `qc_before/SRR40359383_{1,2}_fastqc/` with extracted `summary.txt`; a printed PASS/WARN/FAIL synthesis |
| 3 | `sbatch 3_Batch_MultiQC_Summary.sh` | compute | Both study runs fetched + FastQC'd; `multiqc_raw_batch.html` combined report |
| 4 | `sbatch 4_Evidence_Based_Fastp_Trimming.sh` | compute | `trimmed_reads/SRR40359383_{1,2}_clean.fastq`, `fastp_reports/SRR40359383_THRESHOLD_RATIONALE.md`, and post-trim FastQC in `qc_after/` |
| 5 | `sbatch 5_Paired_End_Integrity_Check.sh` | compute | Printed PASS/FAIL for read-count and read-ID-order pairing checks, raw and trimmed |
| 6 | `sbatch 6_Automated_QC_Pipeline_AllSamples.sh` | compute | Steps 1–5's logic run end-to-end, fault-tolerantly, over **both** study runs; `multiqc_raw_batch.html` and `multiqc_trimmed_batch.html` regenerated |
| 7 | `bash 7_QC_Report_Generator.sh` | login node | `QC_REPORT.md` assembled from every prior step's output |
| 8 | Fill in Section 7 (acceptance criteria) and Section 8 (reflection questions) of `QC_REPORT.md` by hand | login node | Completed portfolio artifact |

Steps 1–5 can be run and understood independently, in order, as the
diagnose-then-decide progression this module teaches. Step 6 does not
require steps 1–5 to have been run first — it re-derives everything itself —
but running 1–5 first makes it much easier to see which stage step 6 is
doing at each point. Step 7 requires step 6 (or steps 1–5 run for both
`SRR40359383` and `SRR40359384`) to have completed, since it only reads
existing files and does not compute anything new.

## Expected outputs, by step

- **Steps 1–3:** `raw_reads/`, `checksums/`, `qc_before/`,
  `multiqc_raw_batch.html`, `logs/*_fetch.log`.
- **Step 4:** `trimmed_reads/`, `fastp_reports/` (including
  `*_THRESHOLD_RATIONALE.md`), `qc_after/`.
- **Step 5:** stdout only (no new files) — a PASS/FAIL line per pairing
  check, per read set (raw/trimmed).
- **Step 6:** everything from steps 1–5 combined, for **all** samples in
  `SRR_LIST`, under `/scratch/$(whoami)/module1_qc/`, plus
  `pairing_reports/<SRR>.txt`, `logs/<SRR>_pipeline.log`, and both combined
  MultiQC reports.
- **Step 7:** `QC_REPORT.md` in the current directory.

## Troubleshooting

| Symptom | Likely cause | Fix |
|---|---|---|
| `ERROR: vdb-validate did not confirm ... is consistent` (script 1) | Interrupted or corrupted `prefetch` download | Delete `sra_cache/<SRR>/` and re-run `prefetch`; check available disk space first |
| `ERROR: Extraction for <SRR> did not produce two non-empty mate files` | Disk quota hit mid-extraction, or the run is actually single-end | Check `df -h` on the target filesystem; confirm the accession is paired-end on the SRA run browser before assuming a bug |
| `ERROR: expected summary at qc_before/.../summary.txt was not produced` (script 2) | FastQC was run without `--extract` | Re-run with `--extract`, or `unzip` the `_fastqc.zip` manually first |
| Script 4 enables `--detect_adapter_for_pe` when you expected it not to (or vice versa) | The evidence-based logic reads `qc_before/<SRR>_*_fastqc/summary.txt` literally — check that file's `Adapter Content` line directly | This is not a bug: re-inspect the actual FastQC verdict for that sample rather than assuming the last sample's result applies |
| Script 5 reports `FAIL: read count mismatch` | One mate was truncated, or reads were filtered independently per-mate somewhere upstream | Re-fetch (script 1) or re-trim (script 4) — never proceed to alignment on a sample that fails this check |
| Script 5 reports `FAIL: read ID order does not match` with matching counts | Mates were reordered or resorted independently (e.g. by a `sort` command applied to only one file) | Re-derive both mates from the original source in the same operation; do not sort/reorder FASTQ mate files independently |
| Script 6 fails on one sample but the job continues | Expected fault-tolerant behavior — check `logs/<SRR>_pipeline.log` for that specific sample | Fix the underlying issue for that accession only; re-run script 6 (it skips re-fetching samples whose raw FASTQ already exists) |
| Script 7 shows `missing` in the diagnostic QC table for a sample | That sample's FastQC/trimming steps were not completed before running script 7 | Run script 6 (or scripts 1–5) for that sample, then re-run script 7 |
| Script 7's Section 4 shows "_No rationale file found_" | `fastp_reports/<SRR>_THRESHOLD_RATIONALE.md` doesn't exist yet for that sample | Run script 4 (or script 6, which generates it inline) for that sample |

## Completion checklist (portfolio artifact)

The curriculum's Module 1 assessment asks for: *"a concise QC report
containing the data source, checksums, diagnostic figures, preprocessing
rationale, command log, and acceptance criteria."* Confirm each is present
before considering the module complete:

- [ ] **Data source** — `QC_REPORT.md` Section 1 correctly names the
      BioProject, GEO accession, organism, and runs analyzed.
- [ ] **Checksums** — Section 2 shows a real `md5sum` line for every sample's
      raw FASTQ mates, not a placeholder/missing message.
- [ ] **Diagnostic figures** — Section 3's table has no `missing` cells; both
      `multiqc_raw_batch.html` and `multiqc_trimmed_batch.html` exist and
      open correctly.
- [ ] **Preprocessing rationale** — Section 4 shows a real, evidence-based
      rationale per sample (from script 4/6), not a "no rationale file
      found" placeholder.
- [ ] **Command log** — `logs/<SRR>_pipeline.log` exists for every sample and
      is non-empty.
- [ ] **Acceptance criteria** — Section 7's PASS/FAIL table is filled in by
      hand for every sample, based on an actual review of Sections 1–6 —
      not left blank.
- [ ] **Reflection** — Section 8's three questions are answered by hand.
- [ ] **Pairing integrity** — Section 5 shows "pairing intact" for both raw
      and trimmed reads, for every sample; any "PAIRING BROKEN" sample is
      resolved (or explicitly excluded with a documented reason) before the
      module is considered complete.
