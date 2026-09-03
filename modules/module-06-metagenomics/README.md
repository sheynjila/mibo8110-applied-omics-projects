# Module 6 — Metagenomics and Microbiome Profiles

**Status: built.**

**Curriculum module:** [Module 6 — Metagenomics and Microbiome Profiles](../../docs/curriculum/modular_bioinformatics_curriculum.md).

## Guiding question

Given a mixed-community sample, what organisms (or genes) are present, at what relative abundance, and how much of that signal is sequencing/database artifact versus real biology?

## Learning outcomes

- Choose the right strategy for the data: amplicon (targeted marker gene) vs. shotgun (whole-community) vs. assembly-based metagenomics — and be able to state, not just intuit, why one route fits a given question and budget better than the other (Tier B Script 1).
- Run the shotgun route with Kraken2/Bracken (`master_microbiome_pipeline.sh`) — MetaPhlAn/HUMAnN are the curriculum's noted alternatives for taxonomic/functional profiling, documented but not implemented here.
- Run the amplicon route with DADA2 (`master_amplicon_pipeline.sh` + `amplicon_dada2_analysis.R`): primer trimming → denoising → ASV inference → chimera removal → taxonomy — mothur+phyloseq is the curriculum's noted alternative, documented but not implemented here.
- Recognize sample-specific contamination risks (host DNA in plant/animal-associated samples, chloroplast/mitochondrial reads in 16S data) and correct for them — and recognize what a pipeline's "unclassified"/"absent" result does and does not tell you about the real community (Tier B Scripts 2 and 4).

## Portfolio artifact (per curriculum)

A microbiome profile report: method chosen and why, QC/contamination handling, abundance results, and an interpretation that acknowledges database and detection-limit caveats. Script 5 below auto-assembles this as `MICROBIOME_REPORT.md` — sections 1-4 are filled in automatically from Scripts 1-4's own output; section 5 (the actual interpretation) is intentionally left blank, printed alongside this module's database/detection-limit caution as unconditional boilerplate rather than optional filler, mirroring Module 2's `DE_REPORT.md`, Module 3's `ASSEMBLY_REPORT.md`, Module 4's `VARIANT_REPORT.md`, and Module 5's `PHYLO_REPORT.md`.

## Relationship to existing repo content

This module's Tier A is entirely pre-existing delivered material, like Modules 2-4 — an update from this module's first build pass, which had to author the amplicon route from scratch because it could not initially be located:

- `master_microbiome_pipeline.sh` (shotgun, Kraken2/Bracken) is a pre-existing corrected script, copied in unchanged. Its own header documents three fixes: Bracken's `-r` read-length is measured per-sample instead of hardcoded to 150, the per-sample loop is fault-tolerant, and — per the [domain coverage audit](../../docs/domain_coverage_audit.md) — it now has an optional host-read depletion step (Bowtie2, `HOST_REFERENCE`, off by default).
- `master_amplicon_pipeline.sh` and `amplicon_dada2_analysis.R` (amplicon, DADA2) are also pre-existing, copied in unchanged — the same audit's later revision records both as delivered under `new_pipelines/`, addressing what its original pass had flagged as "a fundamentally different pipeline... which isn't in this collection." Not a correction of an existing script, but a brand-new pipeline, the same category as this collection's single-cell RNA-seq script. Primer trimming (cutadapt, `AMPLICON_TYPE` switch for 16S/ITS with FastQC before/after) lives in the shell script; denoising through taxonomy (DADA2, no fixed-length truncation for ITS, chloroplast/mitochondrial ASV filtering for 16S) lives in the R script.
- `04_visualize_microbiome.R` is also pre-existing, copied in unchanged — a ggplot2/pheatmap visualization companion to `master_microbiome_pipeline.sh`'s `ALL_SAMPLES_MICROBIOME.tsv` output (stacked relative-abundance bar chart + top-20-taxa heatmap). It is not part of this module's own Tier B numbering (its `04_` prefix is this file's own pre-existing name, not this module's Tier B sequence) — see the Tier A table below.

Tier B (Scripts 1-5) is new pedagogical scaffolding on top, in the numbered-script style established in Modules 2-5. It walks the **shotgun** route's output in depth (Script 1's header explains why) rather than building a second, parallel Tier B for the amplicon route; the amplicon route is validated by its own smoke test, DADA2's own built-in diagnostics, and `04_visualize_microbiome.R`'s own composition plot instead. Worth knowing: `amplicon_dada2_analysis.R`'s own header states it deliberately writes `ALL_SAMPLES_AMPLICON.tsv` in the same `Sample`/`Taxon`/`Estimated_Reads`/`Fraction_Total` long format as the shotgun route's `ALL_SAMPLES_MICROBIOME.tsv` — so Tier B Scripts 3 and 4 (which only read that generic schema, not anything Kraken2/Bracken-specific) run against the amplicon route's output too, unmodified, if you want a diversity/dominance or database-caveat pass on an amplicon cohort. Scripts 1, 2, and 5 are shotgun/Kraken2-specific as written.

**Known dataset/pipeline mismatch, disclosed rather than silently worked around:** the same BioProject (PRJNA917645) sequenced its stool samples by both shallow shotgun AND 16S amplicon — an ideal same-cohort comparison for Learning outcome 1. Its shotgun arm, however, is **single-end**, and `master_microbiome_pipeline.sh` is paired-end throughout (`fasterq-dump --split-files`, `fastp -i/-I`, `kraken2 --paired`, Bowtie2 `--un-conc`). Per this repository's "conserve what's already delivered" policy, `master_microbiome_pipeline.sh` is not edited to add single-end support. This module's shotgun worked example instead uses a different, genuinely paired-end BioProject (PRJNA1520859) — see below — while the amplicon worked example still uses PRJNA917645's own (paired-end) 16S arm. The two routes' worked examples are therefore not the same underlying samples in this module's own pipeline runs, even though a same-cohort comparison was the more pedagogically ideal option, because the actually-delivered Tier A script cannot run against it as-is.

## Worked-example datasets (real, verified via NCBI eutils — same method Modules 3-5 used)

**Shotgun route — BioProject PRJNA1520859**, "Shotgun metagenomic analysis of rhizosphere microbiomes from healthy and Fusarium-diseased tomato plants" (paired-end, Illumina NovaSeq 6000, Nanjing Agricultural University):

| Run | Condition | Note |
|---|---|---|
| `SRR40412948` | Fusarium-diseased tomato rhizosphere, replicate 3 | |
| `SRR40412949` | Fusarium-diseased tomato rhizosphere, replicate 2 | |
| `SRR40412951` | Healthy tomato rhizosphere, replicate 3 | |

Two conditions, not three replicates of one, so the cohort abundance table has a real healthy-vs-diseased contrast to report on. Plant-associated (rhizosphere), so it is a genuine (if modest — rhizosphere soil, not root tissue) candidate for the optional `HOST_REFERENCE` host-depletion step; see `smoke_test_shotgun_prjna1520859.sh`'s header for how to point it at a tomato reference genome.

**Amplicon route — BioProject PRJNA917645**, "Shallow shotgun sequencing reduces technical variation in microbiome analysis" (human gut/stool, 16S rRNA V4 amplicon arm, paired-end, Illumina MiSeq):

| Run | Note |
|---|---|
| `SRR22959853` | Human gut 16S V4 amplicon |
| `SRR22959854` | Human gut 16S V4 amplicon |
| `SRR22959855` | Human gut 16S V4 amplicon |

Human (animal-associated, not plant-associated) stool samples — `amplicon_dada2_analysis.R`'s chloroplast/mitochondrial ASV filter runs correctly against this cohort but is not expected to find much chloroplast signal in stool; see `smoke_test_amplicon_prjna917645.sh`'s header.

## Scripts in this module

**Tier A** (unchanged pipelines + smoke-test wrappers supplying real, verified `srr_list.txt` cohorts in place of both `master_microbiome_pipeline.sh`'s and `master_amplicon_pipeline.sh`'s own hardcoded SRR11092056/SRR11092057 default — already flagged elsewhere in this repository, via NCBI eutils, as SARS-CoV-2 RNA-seq, not a bacterial/environmental shotgun metagenome or a 16S amplicon sample):

| File | Purpose |
|---|---|
| [`scripts/master_microbiome_pipeline.sh`](scripts/master_microbiome_pipeline.sh) | Shotgun metagenomics: SRA download → fastp → optional Bowtie2 host depletion → Kraken2/Bracken taxonomic profiling. |
| [`scripts/04_visualize_microbiome.R`](scripts/04_visualize_microbiome.R) | Visualization companion to the shotgun route: reads `ALL_SAMPLES_MICROBIOME.tsv`, writes a top-10-taxa stacked relative-abundance bar chart and a top-20-taxa heatmap to `plots/`. |
| [`scripts/master_amplicon_pipeline.sh`](scripts/master_amplicon_pipeline.sh) | Amplicon metagenomics: SRA download → FastQC → cutadapt primer trimming (16S or ITS via `AMPLICON_TYPE`) → FastQC → MultiQC. |
| [`scripts/amplicon_dada2_analysis.R`](scripts/amplicon_dada2_analysis.R) | DADA2 core workflow (quality-profile inspection → filter/trim → denoising → ASV inference → chimera removal → taxonomy against SILVA/UNITE); for 16S, filters out chloroplast/mitochondrial ASVs and writes a composition plot. |
| [`scripts/smoke_test_shotgun_prjna1520859.sh`](scripts/smoke_test_shotgun_prjna1520859.sh) | Runs `master_microbiome_pipeline.sh` against the verified 3-sample tomato-rhizosphere cohort above. |
| [`scripts/smoke_test_amplicon_prjna917645.sh`](scripts/smoke_test_amplicon_prjna917645.sh) | Runs `master_amplicon_pipeline.sh` against the verified 3-sample human-gut 16S cohort above. |

**Tier B** (new, one idea per script, walking the shotgun route's output):

| # | File | New idea introduced |
|---|---|---|
| 1 | [`scripts/1_Validate_Classification_Inputs.sh`](scripts/1_Validate_Classification_Inputs.sh) | Confirm each sample's Kraken2 report and Bracken table exist, are well-formed, and belong to each other before anything downstream trusts them. Its header also documents this module's amplicon-vs-shotgun strategy tradeoff (Learning outcome 1) rather than leaving it implicit. |
| 2 | [`scripts/2_Assess_Classification_Yield_And_Host_Fraction.py`](scripts/2_Assess_Classification_Yield_And_Host_Fraction.py) | Report what fraction of each sample's reads Kraken2 could actually classify, and (if `HOST_REFERENCE` was used) what fraction were host and removed — the guiding question's "artifact vs. real biology" split made concrete and per-sample. |
| 3 | [`scripts/3_Cohort_Diversity_And_Dominance_Summary.py`](scripts/3_Cohort_Diversity_And_Dominance_Summary.py) | Compute per-sample Shannon diversity and flag single-taxon dominance (>50% of reads) — a dominance check kept deliberately separate from the diversity number, since the two can disagree about what "looks wrong." |
| 4 | [`scripts/4_Detection_Limit_And_Database_Caveats.py`](scripts/4_Detection_Limit_And_Database_Caveats.py) | For every taxon present in some cohort samples but absent from others, cross-reference that "absence" against Script 2's classification yield for the sample it's missing from, and state explicitly what an "absent" row does and does not prove. |
| 5 | [`scripts/5_Assemble_Microbiome_Report.py`](scripts/5_Assemble_Microbiome_Report.py) | Auto-assemble `MICROBIOME_REPORT.md` from Scripts 1-4, with the objective sections filled in automatically and the interpretation section left blank, printed alongside this module's database/detection-limit caution as unconditional boilerplate. |

### Verification performed

All shell scripts (`master_microbiome_pipeline.sh`, `master_amplicon_pipeline.sh`, `1_Validate_Classification_Inputs.sh`, both `smoke_test_*.sh`) pass `bash -n`. This session's environment has no `cutadapt`, `kraken2`, `bracken`, `bowtie2`, `Rscript`, or real Python interpreter installed (only a Windows Store stub alias), so the Python scripts (2-5) and both R scripts (`amplicon_dada2_analysis.R`, `04_visualize_microbiome.R`) could not be syntax-checked with `python3 -m py_compile` / `Rscript -e "parse(...)"` or functionally run against real or synthetic fixtures the way Modules 2 and 3 were — they received careful manual line-by-line review only, matching the disclosure Modules 4 and 5 already made in this same environment. Note that `master_amplicon_pipeline.sh` and `amplicon_dada2_analysis.R` are pre-existing delivered material, not authored in this repository this session, so their "verification performed" here is limited to this review pass, not independent construction. Completing the functional testing this module still lacks (a small synthetic Kraken2 report + Bracken table cohort, and a small synthetic FASTQ pair for the amplicon route, run through all scripts end-to-end) on a machine with the actual toolchain — following the Stage 0/Stage 1 static-check and smoke-test pattern documented in this module's source materials' `HPC Validation Test Plan.pdf` (not itself committed to this repo) — is the top item in "still to develop" below.

## How to run

**Shotgun route:**
1. `cd scripts/`
2. `sbatch smoke_test_shotgun_prjna1520859.sh` (set `HOST_REFERENCE` first to exercise host depletion) → `bracken_output/ALL_SAMPLES_MICROBIOME.tsv`
3. `bash 1_Validate_Classification_Inputs.sh $WORKDIR`
4. `python3 2_Assess_Classification_Yield_And_Host_Fraction.py $WORKDIR .` → `classification_yield_summary.tsv`
5. `python3 3_Cohort_Diversity_And_Dominance_Summary.py $WORKDIR/bracken_output/ALL_SAMPLES_MICROBIOME.tsv .` → `cohort_diversity_summary.tsv`
6. `python3 4_Detection_Limit_And_Database_Caveats.py $WORKDIR/bracken_output/ALL_SAMPLES_MICROBIOME.tsv classification_yield_summary.tsv .` → `database_detection_caveats.tsv`
7. `python3 5_Assemble_Microbiome_Report.py --yield-tsv classification_yield_summary.tsv --diversity-tsv cohort_diversity_summary.tsv --caveats-tsv database_detection_caveats.tsv --out MICROBIOME_REPORT.md`
8. Open `MICROBIOME_REPORT.md` and fill in Section 5 by hand, alongside its printed caveat.
9. Optional: `Rscript 04_visualize_microbiome.R` (run from `$WORKDIR`, or with `ALL_SAMPLES_MICROBIOME.tsv` copied alongside it) → `plots/01_microbiome_relative_abundance.png`, `plots/02_microbiome_taxa_heatmap.png`.

**Amplicon route:**
1. `cd scripts/`
2. `sbatch smoke_test_amplicon_prjna917645.sh` → `trimmed_reads/*_R[12]_trimmed.fastq`, `final_amplicon_16S_report.html` (MultiQC)
3. Download a DADA2-formatted SILVA (16S) or UNITE (ITS) reference set (see `amplicon_dada2_analysis.R`'s header for the source), edit `REF_DB_PATH` near the top of the script to point at it (it is not read from the working directory by filename convention — the script's own placeholder is `/path/to/...`), then run `Rscript amplicon_dada2_analysis.R` → `ASV_table.tsv`, `taxonomy.tsv`, `ALL_SAMPLES_AMPLICON.tsv`, `amplicon_composition.pdf`, `quality_profiles.pdf`.

See [`RUNBOOK.md`](RUNBOOK.md) for the full step-by-step execution order, required tools, expected outputs, and troubleshooting.

## Teaching materials

- Lecture deck: not yet provided for this module — still to develop.
- Runbook: [`RUNBOOK.md`](RUNBOOK.md).
- Still to develop per [Section 5 of the curriculum](../../docs/curriculum/modular_bioinformatics_curriculum.md): the outstanding functional (not just syntax/manual-review) testing noted above, a guided student workbook, a dataset card/provenance record for both worked-example BioProjects, and a transfer-task assessment. A same-cohort shotgun-vs-amplicon comparison (see the disclosed dataset/pipeline mismatch above) remains a genuine improvement over this module's current two-BioProject worked example, if `master_microbiome_pipeline.sh` is ever extended to support single-end input.
