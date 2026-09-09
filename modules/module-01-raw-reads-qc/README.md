# Module 1 — Raw Sequencing Reads and Quality Control

**Status: built.**

**Curriculum module:** [Module 1 — Raw Sequencing Reads and Quality Control](../../docs/curriculum/modular_bioinformatics_curriculum.md).

## Guiding question

Given a set of raw sequencing reads, how do I confirm they are fit for downstream analysis before I invest compute time in them?

## Learning outcomes

- Explain FASTQ record structure and Phred-quality (Q-score) encoding.
- Retrieve a selected public run and verify file/cache integrity before extracting it.
- Interpret per-base quality, adapter content, duplication, and length/composition summaries diagnostically — not just run the tool and glance at PASS/FAIL colors.
- Distinguish diagnostic QC (FastQC — reports, never modifies data) from automatic trimming (fastp — acts on the evidence).
- Process paired-end reads without breaking read pairing, and verify that pairing directly rather than assuming a tool preserved it.
- Document preprocessing decisions, compare before/after metrics, and preserve a full command/audit log.

## Portfolio artifact (per curriculum)

A QC report that documents: data source and checksums, before/after QC figures, the trimming rationale, the full command log, and explicit acceptance criteria for "reads are ready to align."

## Scripts in this module

Each script adds exactly one new idea on top of the previous one, all built around the same real paired-end study used elsewhere in this repository for continuity: **PRJNA1518998** ("Effect of adipocyte depletion on the retina," *Mus musculus*, Vanderbilt University Medical Center; GEO accession GSE345194; real runs `SRR40359383` and `SRR40359384`).

| # | File | New idea introduced |
|---|---|---|
| 1 | [`scripts/1_Fetch_And_Verify_Raw_Reads.sh`](scripts/1_Fetch_And_Verify_Raw_Reads.sh) | Retrieve one hardcoded run with `prefetch`/`fasterq-dump`, validate the SRA cache with `vdb-validate` *before* extracting, and checksum the resulting FASTQ mates. Includes a FASTQ-record and Phred+33 explainer. |
| 2 | [`scripts/2_Interpret_FastQC_Diagnostics.sh`](scripts/2_Interpret_FastQC_Diagnostics.sh) | Run FastQC on that sample, then *parse* its `summary.txt` PASS/WARN/FAIL verdicts programmatically and interpret the three most commonly misread modules (per-base quality, adapter content, duplication). |
| 3 | [`scripts/3_Batch_MultiQC_Summary.sh`](scripts/3_Batch_MultiQC_Summary.sh) | Loop the fetch+FastQC steps fault-tolerantly over **both** real runs in the study, and aggregate them into one comparison report with MultiQC. |
| 4 | [`scripts/4_Evidence_Based_Fastp_Trimming.sh`](scripts/4_Evidence_Based_Fastp_Trimming.sh) | Read the FastQC evidence and *decide* which fastp flags are justified (not apply them by default), run fastp, write a plain-language `THRESHOLD_RATIONALE.md`, and re-run FastQC on the trimmed output. |
| 5 | [`scripts/5_Paired_End_Integrity_Check.sh`](scripts/5_Paired_End_Integrity_Check.sh) | Verify paired-end integrity two ways: matching mate read **counts**, and matching mate read **ID order** — the second check catches desynchronized mates that the first one misses. Checks both raw and trimmed reads. |
| 6 | [`scripts/6_Automated_QC_Pipeline_AllSamples.sh`](scripts/6_Automated_QC_Pipeline_AllSamples.sh) | Integrate scripts 1–5 into one fault-tolerant `sbatch` pipeline across every sample — the reusable, standalone `qc_pipeline.sh` this module previously lacked. |
| 7 | [`scripts/7_QC_Report_Generator.sh`](scripts/7_QC_Report_Generator.sh) | Auto-assemble every artifact scripts 1–6 produced (checksums, before/after QC, rationale, pairing verdicts, command logs) into the single `QC_REPORT.md` portfolio artifact the curriculum requires, with a hand-completed acceptance-criteria table and reflection questions. |

All seven scripts were validated with `bash -n`. Scripts 5 and 7's core logic (pairing-integrity detection and report assembly) were also functionally smoke-tested against synthetic FASTQ/report fixtures — matched pairs, count mismatches, and read-ID-order mismatches were all correctly identified as PASS/FAIL.

## How to run

See [`RUNBOOK.md`](RUNBOOK.md) for the full step-by-step execution order, expected outputs, and troubleshooting. In short: run scripts 1–5 once on the login node/via `sbatch` to learn each stage individually, then use script 6 as the actual reusable pipeline for a full run, and finish with script 7 (login node, no Slurm needed) to generate the portfolio report.

## What this replaces

The previous version of this module was a **planned stub** noting that QC mechanics were already exercised inline inside modules 2/4/6/7's pipeline scripts, but not yet packaged as its own standalone, reusable teaching module. Those inline QC steps in other modules are untouched by this work — this module adds the missing dedicated lesson plan, diagnostic-interpretation teaching, evidence-based trimming decision process, and portfolio-artifact generator that a standalone Module 1 requires.

## Teaching materials

- Lecture deck: [`docs/Lecture.pptx`](docs/Lecture.pptx).
- Runbook: [`RUNBOOK.md`](RUNBOOK.md).
- Still to develop per [Section 5 of the curriculum](../../docs/curriculum/modular_bioinformatics_curriculum.md): a guided student workbook and a transfer-task assessment using a different accession/organism.
