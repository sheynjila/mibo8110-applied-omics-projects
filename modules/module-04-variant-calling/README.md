# Module 4 — Variant Calling

**Status: built.**

**Curriculum module:** [Module 4 — Variant Calling](../../docs/curriculum/modular_bioinformatics_curriculum.md).

## Guiding question

Given aligned reads against a reference, which variant calls are real, and which are artifacts of mapping bias, ploidy assumptions, or low coverage that I should filter out?

## Learning outcomes

- Map reads with BWA and understand what reference bias means for variant sensitivity — and measure it directly (mapping rate, breadth of coverage), not just assume it away.
- Call variants with BCFtools (this module's tool; GATK is the curriculum's noted alternative) and design a defensible filter (depth, quality, **and strand bias** — not depth/quality alone).
- Reason explicitly about ploidy assumptions when calling variants in non-diploid or mixed organisms — and detect the specific, real consequence of getting it wrong, not just the concept in the abstract.
- Merge multi-sample calls into a cohort VCF for comparative analysis, and sanity-check the merged call set (Ti/Tv, missingness, heterozygosity) before trusting any downstream clustering built on top of it.

## Portfolio artifact (per curriculum)

A filtered VCF + a QC summary + an explicit statement of the call set's limitations (coverage gaps, reference bias, filter thresholds used). Script 7 below auto-assembles this as `VARIANT_REPORT.md` — sections 1-7 are filled in automatically from the pipeline's own output (Scripts 1-6 plus Tier A's `variant_analysis.R`); section 8 (interpretation) is intentionally left for the analyst to complete by hand, mirroring Module 2's `DE_REPORT.md` and Module 3's `ASSEMBLY_REPORT.md`: no script should decide on its own whether a variant call is biologically real.

## Scripts in this module

Two tiers, in the same "add exactly one idea" style used throughout this repository:

**Tier A — the wet-lab pipeline (pre-existing, unchanged):** turns raw paired-end reads into per-sample filtered VCFs and a merged cohort VCF, plus an SNP-distance/outbreak-clustering analysis.

**Tier B — the call-quality and interpretation pipeline (new this session):** seven scripts that take Tier A's output from "a VCF that ran without erroring" to "a defensible, documented variant-calling report." Each script adds exactly one new idea on top of the previous one.

### Tier A — wet-lab pipeline (existing)

| File | Purpose |
|---|---|
| [`scripts/master_snp_pipeline.sh`](scripts/master_snp_pipeline.sh) | End-to-end Slurm pipeline: SRA download → fastp → BWA-MEM alignment → `samtools markdup` deduplication (with proper read groups) → BCFtools variant calling → fixed-threshold filtering → automated multi-sample cohort VCF merge (`merged_cohort.vcf`). Reference: GCF_000006945.2 (*Salmonella enterica* subsp. *enterica* serovar Typhimurium str. LT2, ASM694v2). |
| [`scripts/variant_analysis.R`](scripts/variant_analysis.R) | SNP/Hamming-distance-based sample comparison (`ape::dist.gene`, not Euclidean) with explicit distance-threshold cluster calling (`outbreak_clusters.tsv`) and a Neighbor-Joining tree — bridges into [Module 5 — Phylogenomics](../module-05-phylogenomics/README.md) territory for outbreak/relatedness interpretation. |
| [`scripts/smoke_test_prjna242847.sh`](scripts/smoke_test_prjna242847.sh) | Three-sample smoke-test wrapper (new this session) — see **Data-provenance note** below for why it exists and does not simply rely on the master script's own built-in default. |

### Data-provenance note (verified via NCBI eutils, same method Module 3 used)

`master_snp_pipeline.sh`'s own fallback (used only when `srr_list.txt` is absent) hardcodes `SRR11092056`/`SRR11092057` as a "default test list." Verifying independently confirms **Module 3's exact finding, in the exact same accessions, reused here**: `SRR11092056` is SARS-CoV-2 metagenomic RNA-seq from bronchoalveolar lavage fluid (Wuhan Institute of Virology, BioProject PRJNA605983 / study SRP249613) — not bacterial WGS, and not alignable in any meaningful way against the Typhimurium LT2 reference this pipeline downloads. Aligning it would not error; it would simply produce near-zero, scattered coverage and a variant call set that is pure alignment noise, which is exactly the failure mode this module's Script 2 exists to catch — but far better to never run the worked example against data that cannot succeed in the first place.

Per this repository's "conserve what's already delivered" policy, `master_snp_pipeline.sh`'s internal fallback is **not edited here**. Instead, `smoke_test_prjna242847.sh` supplies a real, verified, matching worked-example cohort:

| Run | BioSample | Organism | Platform | BioProject | Note |
|---|---|---|---|---|---|
| `SRR40474618` | SAMN62852532 | *S. enterica* subsp. *enterica* serovar Typhimurium | Illumina MiSeq, paired-end | PRJNA242847 (NCBI Pathogen Detection / GenomeTrakr, USDA-FSIS submission) | |
| `SRR40472300` | SAMN62840582 | *S. enterica* subsp. *enterica* serovar Typhimurium | Illumina MiSeq, paired-end | PRJNA242847 | |
| `SRR40472280` | SAMN62840173 | *S. enterica* subsp. *enterica* serovar Typhimurium | Illumina MiSeq, paired-end | PRJNA242847 | |

Three isolates (not two) so the cohort merge and `variant_analysis.R`'s SNP-distance clustering have more than one pair to compare.

### Tier B — call-quality and interpretation pipeline (new)

| # | File | New idea introduced |
|---|---|---|
| 1 | [`scripts/1_Validate_Alignment_And_VCF_Inputs.sh`](scripts/1_Validate_Alignment_And_VCF_Inputs.sh) | Confirm every dedup BAM and filtered VCF is structurally sound (`samtools quickcheck`, indexed, non-empty) **and** that the BAM's `@RG SM:` tag and the VCF's sample column match the sample name **by name**, not by filename assumption — the alignment/calling-layer version of the exact bug class Module 2 Script 1 closed for count matrices. |
| 2 | [`scripts/2_Assess_Mapping_And_Reference_Bias.sh`](scripts/2_Assess_Mapping_And_Reference_Bias.sh) | Measure mapping rate and breadth of coverage (% of the reference at zero depth) per sample directly, and flag samples on both counts as the strongest reference-bias signal — a first look at alignment quality before any variant filter is trusted. |
| 3 | [`scripts/3_Design_Defensible_Variant_Filter.sh`](scripts/3_Design_Defensible_Variant_Filter.sh) | Replace the naive fixed `QUAL>=20, DP>=10, MQ>=30` filter with per-cohort, data-driven thresholds (never looser than the naive ones) **plus a strand-bias exclusion read directly from `INFO/DP4`** — closing the exact gap between this module's own stated learning outcome ("depth, quality, strand bias") and what Tier A's filter actually checks (no strand bias at all). Prints a naive-vs-redesigned PASS-count comparison per sample. |
| 4 | [`scripts/4_Check_Ploidy_Assumptions.sh`](scripts/4_Check_Ploidy_Assumptions.sh) | Detect that `bcftools call -mv` (Tier A, no `--ploidy` flag) defaults to **diploid** against a genome that is actually **haploid** (a single bacterial chromosome), count the resulting biologically-impossible heterozygous calls, and re-call with `--ploidy 1` — the real, quantified consequence behind this module's ploidy-assumption learning outcome, not just the concept stated abstractly. |
| 5 | [`scripts/5_Cohort_VCF_QC_Summary.R`](scripts/5_Cohort_VCF_QC_Summary.R) | Compute cohort-level Ti/Tv ratio (flagging a ratio at or near the 0.5 random-substitution floor as noise-dominated), SNP/indel counts, and per-sample missingness/heterozygosity on the merged cohort VCF — a first read of whether the call set as a whole is trustworthy before `variant_analysis.R`'s clustering is asked to mean anything biologically, and a second, cohort-level check that Script 4's ploidy correction actually propagated forward. |
| 6 | [`scripts/6_Coverage_Gaps_And_Limitations.sh`](scripts/6_Coverage_Gaps_And_Limitations.sh) | Auto-assemble the curriculum's required "coverage gaps, reference bias, filter thresholds used" limitations statement (`LIMITATIONS.md`) as its own first-class, machine-checked step, pulling directly from Scripts 2-4's own output — the objective part of "limitations" a script CAN state as fact, kept separate from Script 7's genuinely subjective interpretation section. |
| 7 | [`scripts/7_Assemble_Variant_Report.R`](scripts/7_Assemble_Variant_Report.R) | Auto-assemble `VARIANT_REPORT.md` from Scripts 1-6 and Tier A's `variant_analysis.R` output (`outbreak_clusters.tsv`, tree/QC plots), with every objective section filled in automatically and a genuine "Interpretation" section left blank for the analyst. |

### Verification performed on Tier B

Every Tier B script was **syntax-checked**: `bash -n` for the four shell scripts (all pass), manual line-by-line review for the two R scripts. Unlike Modules 2 and 3, this session's environment has no `Rscript`, `bcftools`, `samtools`, or `bwa` installed, so these scripts could **not** be functionally run end-to-end against real or synthetic fixtures here — that verification step (synthetic-fixture runs, the way Module 2 tested against noise/seeded-signal count matrices and Module 3 tested against engineered repeat/contamination genomes) is still outstanding and should be completed on a machine or HPC node with the actual toolchain before this module is treated as fully verified, not just written. This is stated plainly rather than claimed away.

## How to run

1. Run the Tier A wet-lab pipeline first: `bash smoke_test_prjna242847.sh` (recommended — verified real accessions) or `sbatch master_snp_pipeline.sh` with your own `srr_list.txt`. Produces `alignment/<SRR>_dedup.bam`, `variants/<SRR>_raw.vcf`, `variants/<SRR>_filtered.vcf`, and `variants/merged_cohort.vcf`.
2. Run Tier B scripts 1→7 in order, each taking the pipeline's `WORKDIR` (default `/scratch/$(whoami)/master_snp_pipeline`) as its argument:
   - `bash 1_Validate_Alignment_And_VCF_Inputs.sh [WORKDIR]`
   - `bash 2_Assess_Mapping_And_Reference_Bias.sh [WORKDIR]` → writes `mapping_bias_summary.tsv`
   - `bash 3_Design_Defensible_Variant_Filter.sh [WORKDIR]` → writes `<SRR>_filtered_v2.vcf` per sample
   - `bash 4_Check_Ploidy_Assumptions.sh [WORKDIR]` → writes `<SRR>_raw_haploid.vcf` per sample
   - **Ploidy-correct re-run:** for a fully defensible result, re-point Script 3 at `<SRR>_raw_haploid.vcf` instead of `<SRR>_raw.vcf`, and re-merge with `bcftools merge` before Script 5 — Script 4 corrects the calling step; it does not rewrite Tier A's merge for you.
   - `Rscript 5_Cohort_VCF_QC_Summary.R variants/merged_cohort.vcf` → writes `cohort_vcf_qc_summary.tsv`, `cohort_vcf_qc_per_sample.tsv`
   - `bash 6_Coverage_Gaps_And_Limitations.sh [WORKDIR]` → writes `LIMITATIONS.md`
   - `Rscript 7_Assemble_Variant_Report.R [WORKDIR]` → writes `VARIANT_REPORT.md`
3. Open `VARIANT_REPORT.md` and fill in Section 8 (interpretation) by hand — the intentional human-judgment step the curriculum requires.

See [`RUNBOOK.md`](RUNBOOK.md) for the full step-by-step execution order, required tools, expected outputs, and troubleshooting.

## Corrections applied (see [changelog](../../docs/corrected_scripts_changelog.md) for full detail)

Automated the cohort VCF merge (previously a manual copy-paste step), added download integrity checks, and a fault-tolerant per-sample loop — all pre-existing corrections to `master_snp_pipeline.sh`/`variant_analysis.R`, unchanged this session. New this session: a verified real worked-example cohort (the original two-accession default is confirmed non-bacterial, see the data-provenance note above), a strand-bias-aware filter redesign (Script 3), an identified-and-corrected ploidy default (Script 4), and a cohort-level QC sanity check (Script 5) — none of which existed in any prior version of this pipeline.

## Teaching materials

- Lecture deck: not yet provided for this module — still to develop.
- Runbook: [`RUNBOOK.md`](RUNBOOK.md).
- Still to develop per [Section 5 of the curriculum](../../docs/curriculum/modular_bioinformatics_curriculum.md): guided student workbook, dataset card/provenance record for PRJNA242847, a transfer-task assessment, and the outstanding functional (not just syntax) testing noted above.
