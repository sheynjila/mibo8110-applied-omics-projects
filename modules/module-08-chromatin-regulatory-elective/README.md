# Module 8 (Elective) — Chromatin and Regulatory Genomics

**Status: built.**

**Curriculum module:** [Module 8 (Elective) — Chromatin and Regulatory Genomics](../../docs/curriculum/modular_bioinformatics_curriculum.md).

## Guiding question

Given a ChIP-seq, ATAC-seq, or bisulfite sequencing dataset, where in the genome is a protein bound, chromatin accessible, or DNA methylated — and how do I distinguish a real regulatory signal from assay background?

## Learning outcomes

- Align reads appropriate to the assay (Bowtie2, the curriculum's named aligner for this module — `master_chipseq_pipeline.sh`).
- Call peaks with MACS2 (ChIP-seq/ATAC-seq) and interpret peak scores and reproducibility across replicates — Tier B Script 2 (FRiP) and Script 3 (pairwise Jaccard reproducibility across real biological replicates).
- Use BEDTools for interval operations (peak overlap, annotation against gene features) — `master_chipseq_pipeline.sh`'s `bedtools closest` step and Tier B Script 4.
- Reason about assay-specific background and artifacts (e.g. ATAC-seq's Tn5 insertion bias, ChIP-seq's input-control necessity) — `master_chipseq_pipeline.sh`'s `ASSAY_TYPE` switch and Tier B Script 1's explicit input-control check.

## Portfolio artifact (per curriculum)

A peak call set + QC summary + an interpretation of the regulatory signal, with explicit discussion of assay-specific background handling. Script 5 below auto-assembles this as `CHROMATIN_REPORT.md` — sections 1-4 are filled in automatically from Scripts 1-4's own output; section 5 (the actual regulatory interpretation) is intentionally left blank, printed alongside this module's core "signal vs. background/function" caution as unconditional boilerplate rather than optional filler, mirroring every other built module's report script in this repository.

**Scope note:** bisulfite sequencing / DNA methylation, named in the curriculum's guiding question alongside ChIP-seq and ATAC-seq, is **not** covered by this build — the source README draft's own "Suggested next step" scoped this module to `master_chipseq_pipeline.sh` (ChIP-seq/ATAC-seq via MACS2) specifically; methylation calling uses a different toolchain entirely (e.g. Bismark, methylKit) and is a genuine, disclosed gap, not an oversight — see "Still to develop".

## Relationship to existing repo content

No pipeline in this repository touched chromatin/regulatory assays before this build — a genuine gap with no partial coverage elsewhere, exactly as this module's own draft README stated. This module has **no pre-existing Tier A** (like Module 5), built entirely fresh for this session: `master_chipseq_pipeline.sh` is a brand-new pipeline (not a correction of an existing script), the same category as this repository's single-cell RNA-seq and amplicon metagenomics pipelines. One script, two modes via `ASSAY_TYPE=ChIP|ATAC` — the same MACS2 tool, called with assay-appropriate parameters (input-control-aware for ChIP; Tn5-bias-corrected for ATAC), matching how Module 6's amplicon pipeline handles `AMPLICON_TYPE=16S|ITS` with one script. Tier B (Scripts 1-5) is new pedagogical scaffolding on top, in the numbered-script style established in Modules 2-7.

## Worked-example dataset (real, verified via NCBI eutils — same method Modules 3-7 used)

**BioProject PRJNA1283515**, "KLF5 controls subtype-independent highly interactive enhancers in pancreatic cancer to regulate cell survival" — a large (149-run) multi-assay study (ChIP-seq, ATAC-seq, HiChIP, across three pancreatic cancer cell lines) that happens to cover exactly this module's two implemented assay types **in the same cell line**, letting the ChIP and ATAC worked examples be read side by side rather than as two unrelated datasets:

| Run | Assay | Sample | Note |
|---|---|---|---|
| `SRR35978255` | CTCF ChIP-seq | T3M4, Rep1 | |
| `SRR35978254` | CTCF ChIP-seq | T3M4, Rep2 | |
| `SRR35978253` | CTCF ChIP-seq | T3M4, Rep3 | |
| `SRR35978262` | ATAC-seq | T3M4, Rep1 | |
| `SRR35978261` | ATAC-seq | T3M4, Rep2 | |
| `SRR35978259` | ATAC-seq | T3M4, Rep3 | |

All six are real biological replicates (not technical replicates of one run), paired-end, Illumina NextSeq 2000.

**A real, disclosed finding baked into this worked example, not smoothed over:** this BioProject's full 149-run deposit contains **no input/IgG control sample for any of its ChIP-seq experiments** (CTCF, dNp63, H3K27ac, across all three cell lines) — confirmed by NCBI eutils search for "IgG", "Input", and "control" across the whole project, all returning zero hits. This is a real published study, not this module's own mistake, and it lands directly on this module's own learning outcome about ChIP-seq's input-control necessity: the CTCF worked example above is deliberately run **without** an input control (`CONTROL_SRR` left unset), and Tier B Script 1 detects and flags this condition explicitly from `master_chipseq_pipeline.sh`'s own MACS2 log rather than silently treating an uncontrolled run as equivalent to a controlled one. Same pattern as Module 6's disclosed shotgun/amplicon dataset mismatch and Module 7's disclosed 10x/InDrop chemistry mismatch — a real data problem surfaced and taught from, not hidden.

## Scripts in this module

**Tier A** (new pipeline + smoke-test wrappers supplying real, verified `srr_list.txt` cohorts in place of `master_chipseq_pipeline.sh`'s own hardcoded SRR11092056/SRR11092057 default — already flagged elsewhere in this repository, via NCBI eutils, as SARS-CoV-2 RNA-seq, not ChIP-seq/ATAC-seq of any kind):

| File | Purpose |
|---|---|
| [`scripts/master_chipseq_pipeline.sh`](scripts/master_chipseq_pipeline.sh) | SRA download → fastp → Bowtie2 alignment + duplicate marking → MACS2 peak calling (`ASSAY_TYPE=ChIP` with optional input control, or `ASSAY_TYPE=ATAC` with Tn5-bias-corrected parameters) → BEDTools peak-to-gene annotation. |
| [`scripts/smoke_test_chip_prjna1283515.sh`](scripts/smoke_test_chip_prjna1283515.sh) | Runs `master_chipseq_pipeline.sh` (`ASSAY_TYPE=ChIP`, no input control — see above) against the verified 3-replicate CTCF cohort. |
| [`scripts/smoke_test_atac_prjna1283515.sh`](scripts/smoke_test_atac_prjna1283515.sh) | Runs `master_chipseq_pipeline.sh` (`ASSAY_TYPE=ATAC`) against the verified 3-replicate ATAC-seq cohort — same cell line as the ChIP smoke test. |

**Tier B** (new, one idea per script):

| # | File | New idea introduced |
|---|---|---|
| 1 | [`scripts/1_Validate_Peak_Calling_Inputs.sh`](scripts/1_Validate_Peak_Calling_Inputs.sh) | Confirm each sample's BAM and peaks file exist, are well-formed, and belong to each other — and detect, from MACS2's own log, whether that sample's peaks were actually called with an input control, rather than assuming a batch was consistent. |
| 2 | [`scripts/2_Compute_FRiP_Score.sh`](scripts/2_Compute_FRiP_Score.sh) | Compute Fraction of Reads in Peaks per sample (the ENCODE-standard metric) — this module's guiding question, "real signal vs. background," turned into one number: a high peak count with low FRiP is exactly the "looks real, isn't" case peak count alone cannot catch. |
| 3 | [`scripts/3_Assess_Replicate_Reproducibility.sh`](scripts/3_Assess_Replicate_Reproducibility.sh) | Pairwise Jaccard similarity between every pair of peak sets — a position-aware reproducibility measure, since two replicates can report the same peak COUNT while sharing almost no actual genomic positions. |
| 4 | [`scripts/4_Summarize_Peak_Gene_Annotation.py`](scripts/4_Summarize_Peak_Gene_Annotation.py) | Summarize Tier A's own `bedtools closest` output into a promoter-proximal vs. distal breakdown per sample — "where in the genome," made concrete. |
| 5 | [`scripts/5_Assemble_Chromatin_Report.py`](scripts/5_Assemble_Chromatin_Report.py) | Auto-assemble `CHROMATIN_REPORT.md` from Scripts 1-4, with the objective sections filled in automatically and the regulatory-interpretation section left blank, printed alongside this module's signal-vs-function caution as unconditional boilerplate. |

### Verification performed

All shell scripts (`master_chipseq_pipeline.sh`, Tier B Scripts 1-3, both `smoke_test_*.sh`) pass `bash -n`. This session's environment has no `Bowtie2`, `MACS2`, `BEDTools`, `samtools`, or real Python interpreter installed (only a Windows Store stub alias), so Tier B Scripts 4-5 could not be syntax-checked with `python3 -m py_compile` or functionally run against real or synthetic fixtures — they received careful manual line-by-line review only, matching the disclosure every module since Module 4 has made in this same environment. Completing that functional testing (a small synthetic BAM + narrowPeak + gene-features BED, run through Scripts 1-5 end-to-end) on a machine with the actual toolchain is the top item in "still to develop" below.

## How to run

**ChIP-seq route:**
1. `cd scripts/`
2. `sbatch smoke_test_chip_prjna1283515.sh` → `alignment/<SRR>_dedup.bam`, `macs2_output/<SRR>/<SRR>_peaks.narrowPeak`, `annotation/<SRR>_peaks_annotated.bed`
3. `bash 1_Validate_Peak_Calling_Inputs.sh $WORKDIR` (note the input-control finding for Step 8 below)
4. `bash 2_Compute_FRiP_Score.sh $WORKDIR .` → `frip_summary.tsv`
5. `bash 3_Assess_Replicate_Reproducibility.sh $WORKDIR .` → `replicate_reproducibility.tsv`
6. `python3 4_Summarize_Peak_Gene_Annotation.py $WORKDIR .` → `peak_annotation_summary.tsv`
7. `python3 5_Assemble_Chromatin_Report.py --assay-type ChIP --control-status "3/3 samples called WITHOUT an input control (none deposited for PRJNA1283515's ChIP-seq arm)" --frip-tsv frip_summary.tsv --reproducibility-tsv replicate_reproducibility.tsv --annotation-tsv peak_annotation_summary.tsv`
8. Open `CHROMATIN_REPORT.md` and fill in Section 5 by hand, alongside its printed caveat.

**ATAC-seq route:** same as above with `smoke_test_atac_prjna1283515.sh` in Step 2 and `--assay-type ATAC --control-status "not applicable — ATAC-seq has no input-control concept"` in Step 7.

See [`RUNBOOK.md`](RUNBOOK.md) for the full step-by-step execution order, required tools, expected outputs, and troubleshooting.

## Teaching materials

- Lecture deck: not yet provided for this module — still to develop.
- Runbook: [`RUNBOOK.md`](RUNBOOK.md).
- Still to develop per [Section 5 of the curriculum](../../docs/curriculum/modular_bioinformatics_curriculum.md): the outstanding functional (not just syntax/manual-review) testing noted above, a guided student workbook, a dataset card/provenance record for PRJNA1283515, and a transfer-task assessment. Bisulfite sequencing / methylation calling (the guiding question's third assay type) is out of scope for this build entirely — a genuine, separate future addition, not a small extension of `master_chipseq_pipeline.sh`.
