# Module 8 (Elective) — Chromatin and Regulatory Genomics: Runbook

This is an **operational** document — step-by-step execution instructions,
expected outputs, and troubleshooting. For the conceptual "why," see this
module's [`README.md`](README.md). For the full curriculum context, see
[`docs/curriculum/modular_bioinformatics_curriculum.md`](../../docs/curriculum/modular_bioinformatics_curriculum.md).

## Prerequisites

- SRA-Toolkit, fastp, Bowtie2, SAMtools, MACS2, BEDTools, MultiQC — loaded
  automatically via [`scripts/load_modules.sh`](scripts/load_modules.sh)
  (pinned to SRA-Toolkit 3.2.0, fastp 0.23.4, Bowtie2 2.5.4, SAMtools 1.18,
  MACS2 2.2.9.1, BEDTools 2.31.0, MultiQC 1.28, with an automatic fallback
  if your cluster doesn't have those exact builds — see the root
  [`RUNBOOK.md`](../../RUNBOOK.md) §4).
- A GNU awk (`gawk`) providing 3-argument `match()` — used to build
  `ref/gene_features.bed` from the downloaded GTF. Most HPC systems' default
  `awk` is gawk; confirm with `awk --version` if the gene-features build
  step fails with a `match()` argument-count error.
- Bash — used by Tier B Scripts 1-3.
- Python 3, standard library only — Tier B Scripts 4-5 use no third-party
  packages, the same choice Modules 5-7 made.

## Execution order — ChIP-seq route

| Step | Command | Expected result |
|---|---|---|
| 0 | `sbatch smoke_test_chip_prjna1283515.sh` | `$WORKDIR/alignment/<SRR>_dedup.bam`, `$WORKDIR/macs2_output/<SRR>/<SRR>_peaks.narrowPeak`, `$WORKDIR/annotation/<SRR>_peaks_annotated.bed`, `$WORKDIR/final_chromatin_ChIP_report.html` (MultiQC) |
| 1 | `bash 1_Validate_Peak_Calling_Inputs.sh $WORKDIR` | Console PASS/FAIL/WARN report per sample, including input-control status |
| 2 | `bash 2_Compute_FRiP_Score.sh $WORKDIR .` | Console per-sample FRiP report; `frip_summary.tsv` |
| 3 | `bash 3_Assess_Replicate_Reproducibility.sh $WORKDIR .` | Console pairwise Jaccard report; `replicate_reproducibility.tsv` |
| 4 | `python3 4_Summarize_Peak_Gene_Annotation.py $WORKDIR .` | Console per-sample promoter/distal breakdown; `peak_annotation_summary.tsv` |
| 5 | `python3 5_Assemble_Chromatin_Report.py --assay-type ChIP --control-status "<Step 1's actual finding>" --frip-tsv frip_summary.tsv --reproducibility-tsv replicate_reproducibility.tsv --annotation-tsv peak_annotation_summary.tsv` | `CHROMATIN_REPORT.md`, sections 1-4 auto-filled, section 5 blank |
| 6 | Fill in Section 5 (interpretation) of `CHROMATIN_REPORT.md` by hand, alongside its printed caveat | Completed portfolio artifact |

`$WORKDIR` defaults to `/scratch/$(whoami)/master_chipseq_pipeline_ChIP`.

## Execution order — ATAC-seq route

Identical to the ChIP-seq route above, substituting `smoke_test_atac_prjna1283515.sh`
in Step 0 and `--assay-type ATAC --control-status "not applicable — ATAC-seq
has no input-control concept"` in Step 5. `$WORKDIR` defaults to
`/scratch/$(whoami)/master_chipseq_pipeline_ATAC`.

## Expected outputs, by step

- **Step 0:** per-sample BAM, narrowPeak, and annotated-BED files; MultiQC report.
- **Step 1:** console validation report only (no file written).
- **Step 2:** `frip_summary.tsv`.
- **Step 3:** `replicate_reproducibility.tsv`.
- **Step 4:** `peak_annotation_summary.tsv`.
- **Step 5:** `CHROMATIN_REPORT.md`.

## Troubleshooting

| Symptom | Likely cause | Fix |
|---|---|---|
| Step 0 fails building `ref/gene_features.bed` with a `match()` error | The system's default `awk` is not gawk (no 3-argument `match()` support) | Load a gawk module explicitly, or rewrite that one `awk` call with `gsub`/`substr` for a POSIX-awk-compatible extraction |
| Step 0's `bowtie2-build` step is extremely slow | Expected — a full human genome index build is a multi-hour, not multi-minute, step (same caveat Module 7's STAR index build carries) | Run it once and let `master_chipseq_pipeline.sh`'s own `[ ! -f ... ]` check skip rebuilding on subsequent runs |
| Step 0 (ChIP mode) reports MACS2 ran without `-c` even though you set `CONTROL_SRR` | `CONTROL_SRR` was set only in your current shell, not exported before `sbatch` submission (Slurm jobs do not inherit unexported shell variables) | Use `export CONTROL_SRR=...` before `sbatch`, or add it to the smoke-test wrapper directly |
| Step 1 reports a sample with NO input control for a batch you believed was controlled | Confirm what `master_chipseq_pipeline.sh`'s own Phase 3 console output said for that sbatch submission — `CONTROL_BAM` is only set if `CONTROL_SRR` was non-empty at the time that job ran | Re-run with `CONTROL_SRR` correctly exported; do not edit Script 1 to stop flagging this instead |
| Step 2 reports several samples `LOW_FRIP` | Could be a genuinely weak/failed experiment, OR an expected property of the assay (broad histone marks, ATAC-seq — see Script 2's header) | Check the assay type and mark/factor before concluding the experiment failed; do not lower `LOW_FRIP_PCT` just to silence the warning |
| Step 3 reports `LOW_REPRODUCIBILITY` for a pair that are NOT actually biological replicates of each other | Expected — Script 3 computes every pairwise comparison in the batch, including cross-condition pairs that are supposed to differ (see its own SCOPE NOTE) | Confirm from `srr_list.txt`/your own sample metadata which pairs are real replicate comparisons before treating a flag as a problem |
| Step 4 reports many peaks with distance `-1` excluded from the summary | `bedtools closest`'s sentinel for "no gene feature found on this peak's chromosome at all" — usually a contig/scaffold with no annotated genes, or a chromosome-naming mismatch between the FASTA and GTF | Confirm `ref/gene_features.bed`'s chromosome names match `macs2_output/*/*_peaks.narrowPeak`'s (e.g. `chr1` vs `1` — Ensembl FASTA/GTF use bare numbers, not a `chr` prefix) |
| Script 5's `CHROMATIN_REPORT.md` shows "Not available" for a section | The corresponding earlier script's output file was not passed via its `--*-tsv` flag | Re-run the missing step, then re-run Script 5 with the correct file path |
| Script 5's Section 5 is blank | Expected — intentionally left for the analyst, alongside the printed caveat | Fill it in referencing Sections 1-4's actual numbers, and do not remove the caveat text above it |

## Completion checklist (portfolio artifact)

The curriculum's Module 8 assessment asks for: *"A peak call set + QC
summary + an interpretation of the regulatory signal, with explicit
discussion of assay-specific background handling."* Confirm each is present
before considering the module complete:

- [ ] **Peak call set** — `macs2_output/<SRR>/<SRR>_peaks.narrowPeak` for
      every sample is retained alongside the report, not just its summary.
- [ ] **QC summary** — Sections 1-4 of `CHROMATIN_REPORT.md` reflect real
      Script 1-4 runs against your actual cohort, including any `WARN`/`LOW_*`
      flags, not just clean-looking numbers.
- [ ] **Assay-specific background handling, discussed explicitly** —
      Section 1's input-control status (ChIP) or Tn5-bias-correction note
      (ATAC) is present and specific to the actual run, not the generic
      template text.
- [ ] **Written regulatory interpretation** — Section 5 is filled in by
      hand, explicitly distinguishing what the peaks/FRiP/reproducibility
      support (a real, reproducible signal at a position) from what they do
      NOT establish on their own (that the signal has a specific regulatory
      function) — the caveat text above it is not removed or paraphrased away.
