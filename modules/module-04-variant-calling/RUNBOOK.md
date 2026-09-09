---
output:
  pdf_document: default
  html_document: default
---
# Module 4 — Variant Calling: Runbook

This is an **operational** document — step-by-step execution instructions,
expected outputs, and troubleshooting. For the conceptual "why," see this
module's [`README.md`](README.md). For the full curriculum context, see
[`docs/curriculum/modular_bioinformatics_curriculum.md`](../../docs/curriculum/modular_bioinformatics_curriculum.md).

## Prerequisites

- Slurm environment for Tier A (`master_snp_pipeline.sh`) — environment
  modules (SRA-Toolkit, fastp, BWA, SAMtools, BCFtools, MultiQC) are loaded
  automatically via [`scripts/load_modules.sh`](scripts/load_modules.sh);
  see the root [`RUNBOOK.md`](../../RUNBOOK.md) §4 for how the
  pinned-build / discover-by-name fallback works and how to override a
  module name for your own cluster.
- `samtools` and `bcftools` on `PATH` (outside Slurm too) for Tier B Scripts
  1, 2, 3, 4, and 6 — these are plain shell scripts, not Slurm jobs, and are
  meant to be run interactively against Tier A's output.
- R with `vcfR`, `adegenet`, and `ape` installed — used by Tier A's
  `variant_analysis.R` and Tier B Scripts 5 and 7.
- Module 1 completed or reviewed — this module assumes the same
  "validate/interpret before you trust it" habit taught there, applied one
  level up the pipeline (aligned reads and called variants instead of raw
  reads).

## A real, verified data-provenance issue (read before running the smoke test)

`master_snp_pipeline.sh`'s own built-in fallback `srr_list.txt` (used only
when no `srr_list.txt` is staged) is `SRR11092056`/`SRR11092057`. Verified via
NCBI eutils: `SRR11092056` is SARS-CoV-2 metagenomic RNA-seq (Wuhan Institute
of Virology), not bacterial WGS, and cannot produce a meaningful alignment
against this pipeline's Typhimurium LT2 reference. **Do not rely on the
built-in fallback.** Use `smoke_test_prjna242847.sh`, which stages a real,
verified 3-isolate Typhimurium cohort (PRJNA242847) before invoking
`master_snp_pipeline.sh`. See the README's data-provenance note for the full
verification trail.

## Execution order

| Step | Command | Expected result |
|---|---|---|
| 1 | `bash smoke_test_prjna242847.sh` (or `sbatch master_snp_pipeline.sh` with your own verified `srr_list.txt`) | `variants/<SRR>_raw.vcf`, `variants/<SRR>_filtered.vcf`, `variants/merged_cohort.vcf` |
| 2 | `bash 1_Validate_Alignment_And_VCF_Inputs.sh [WORKDIR]` | Console PASS/FAIL/WARN report per sample; stops hard on a structural or name-mismatch failure |
| 3 | `bash 2_Assess_Mapping_And_Reference_Bias.sh [WORKDIR]` | `mapping_bias_summary.tsv`; console table of % mapped, % properly paired, mean depth, % zero-coverage per sample |
| 4 | `bash 3_Design_Defensible_Variant_Filter.sh [WORKDIR]` | `variants/<SRR>_filtered_v2.vcf`, `variants/<SRR>_strand_biased_sites.tsv`; console naive-vs-v2 PASS-count comparison |
| 5 | `bash 4_Check_Ploidy_Assumptions.sh [WORKDIR]` | `variants/<SRR>_raw_haploid.vcf`; console before/after heterozygous-call count per sample |
| 6 | Re-merge the ploidy-corrected calls: `cd WORKDIR/variants && for f in *_raw_haploid.vcf; do bcftools view -Oz -o ${f%.vcf}.vcf.gz $f && bcftools index -t ${f%.vcf}.vcf.gz; done && bcftools merge *_raw_haploid.vcf.gz -Ov -o merged_cohort_haploid.vcf` | `variants/merged_cohort_haploid.vcf` — the ploidy-correct cohort VCF |
| 7 | `Rscript 5_Cohort_VCF_QC_Summary.R WORKDIR/variants/merged_cohort_haploid.vcf` | `cohort_vcf_qc_summary.tsv`, `cohort_vcf_qc_per_sample.tsv`; console Ti/Tv and heterozygosity report |
| 8 | `bash 6_Coverage_Gaps_And_Limitations.sh [WORKDIR]` | `LIMITATIONS.md` |
| 9 | `cd WORKDIR/variants && Rscript variant_analysis.R` (Tier A, against `merged_cohort_haploid.vcf` — see script header for the required `bcftools merge` setup if not reusing step 6's output) | `variant_quality_qc.pdf`, `phylogenetic_tree.pdf`, `outbreak_clusters.tsv` (if any pair is within `SNP_THRESHOLD`) |
| 10 | `Rscript 7_Assemble_Variant_Report.R WORKDIR` | `VARIANT_REPORT.md`, sections 1-7 auto-filled, section 8 blank |
| 11 | Fill in Section 8 (interpretation) of `VARIANT_REPORT.md` by hand | Completed portfolio artifact |

Step 6 is the one manual hand-off in this module: Script 4 corrects the
*calling* step but does not silently rewrite Tier A's cohort merge for you —
running Scripts 5-7 against the original `merged_cohort.vcf` instead of the
ploidy-corrected re-merge will produce numbers Script 5 itself will flag
(non-zero heterozygosity against a haploid reference), not numbers that are
silently wrong.

## Expected outputs, by step

- **Step 1 (Tier A):** `alignment/<SRR>_dedup.bam[.bai]`, `variants/<SRR>_raw.vcf`, `variants/<SRR>_filtered.vcf`, `variants/merged_cohort.vcf`, `final_multiqc_report.html`.
- **Step 2:** stdout only — validation report.
- **Step 3:** `mapping_bias_summary.tsv`.
- **Step 4:** `variants/<SRR>_filtered_v2.vcf`, `variants/<SRR>_strand_biased_sites.tsv`, `variants/<SRR>_raw_metrics.tsv`.
- **Step 5:** `variants/<SRR>_raw_haploid.vcf`.
- **Step 6:** `variants/merged_cohort_haploid.vcf`.
- **Step 7:** `cohort_vcf_qc_summary.tsv`, `cohort_vcf_qc_per_sample.tsv`.
- **Step 8:** `LIMITATIONS.md`.
- **Step 9 (Tier A):** `variant_quality_qc.pdf`, `phylogenetic_tree.pdf`, `outbreak_clusters.tsv` (conditional).
- **Step 10:** `VARIANT_REPORT.md`.

## Troubleshooting

| Symptom | Likely cause | Fix |
|---|---|---|
| `smoke_test_prjna242847.sh` fails at the `prefetch`/`fasterq-dump` step | SRA-Toolkit module not loaded, or a transient NCBI connectivity issue | Confirm `module load SRA-Toolkit/3.2.0-gompi-2024a` succeeded; retry — `master_snp_pipeline.sh`'s per-sample loop is fault-tolerant and will skip a genuinely failed sample rather than abort the whole run |
| Script 1 fails with `@RG SM: tag does not match filename-derived sample` | A BAM was renamed/copied after alignment without regenerating its read group, or two samples' outputs were cross-copied | Re-run alignment for that sample from `master_snp_pipeline.sh` rather than renaming files by hand — the `-R "@RG\tID:...\tSM:..."` flag is what makes this checkable at all |
| Script 2 flags most/all samples `LOW_MAPPING_RATE` | You are accidentally running this against the pipeline's own broken default accessions (SARS-CoV-2 RNA-seq vs. a bacterial DNA reference) instead of a verified worked-example cohort | Confirm `srr_list.txt` matches `smoke_test_prjna242847.sh`'s accessions or your own verified bacterial WGS list — see the README's data-provenance note |
| Script 3 reports 0 strand-biased sites for every sample | Either genuinely clean data, or too few candidate calls with `INFO/DP4` ALT depth >= 4 to evaluate (very shallow coverage) | Cross-check against Script 2's mean-depth column before treating "0 flagged" as "no strand bias present" |
| Script 4 reports 0 heterozygous calls under the diploid default | Does not confirm ploidy was handled correctly — it means this particular run got lucky, not that the risk doesn't exist. The `*_raw_haploid.vcf` files are still the correct calls to build on | Always continue with the haploid-corrected files regardless of whether this run happened to show 0 het calls |
| Script 5 reports non-zero heterozygosity | You ran it against Tier A's original `merged_cohort.vcf` (diploid-default calls) instead of the ploidy-corrected re-merge from Step 6 | Re-merge `*_raw_haploid.vcf` per Step 6 and re-run Script 5 against that file |
| Script 5 reports Ti/Tv <= 0.6 | Cohort's SNP calls skew toward random noise rather than real biological signal — could be genuinely low-quality data, or evidence you are still running against the wrong (viral RNA-seq vs. bacterial reference) dataset | Confirm the dataset first (see the mapping-rate check above); if the dataset is correct, treat this as a real quality warning, not something to filter past silently |
| `variant_analysis.R` (Tier A) errors on `read.vcfR("merged_cohort.vcf")` | Working directory is not `WORKDIR/variants`, or the file wasn't regenerated after Step 6's re-merge | `cd` into `WORKDIR/variants` first; confirm the intended cohort VCF (original vs. haploid-corrected) actually exists there under the name the script expects |
| Script 7's `VARIANT_REPORT.md` shows "Not available" for a section | The corresponding earlier script was not run, or was run in a different working directory than the one passed to Script 7 | Re-run the missing step with the same `WORKDIR` argument, then re-run Script 7 |
| Script 7's Section 8 is blank | Expected — intentionally left for the analyst to complete by hand | Fill it in based on the actual mapping/coverage flags (Section 2), Ti/Tv (Section 5), and clustering results (Section 6) — not the auto-generated prompt text |

## Completion checklist (portfolio artifact)

The curriculum's Module 4 assessment asks for: *"A filtered VCF + a QC
summary + an explicit statement of the call set's limitations (coverage
gaps, reference bias, filter thresholds used)."* Confirm each is present
before considering the module complete:

- [ ] **Filtered VCF** — `<SRR>_filtered_v2.vcf` (Script 3) or, better,
      the ploidy-corrected `merged_cohort_haploid.vcf` (Step 6) is retained
      alongside the report, not just summarized.
- [ ] **QC summary** — Section 5 of `VARIANT_REPORT.md` reflects a real
      Script 5 run against your actual merged cohort VCF, including the
      Ti/Tv value and per-sample heterozygosity, not a placeholder.
- [ ] **Coverage gaps / reference bias** — Section 2 reflects a real
      Script 2 run (per-sample mapping rate and % zero-coverage), not a
      placeholder.
- [ ] **Filter thresholds used** — Section 3 states the exact per-sample
      thresholds Script 3 actually computed for this cohort, not the naive
      fixed defaults alone.
- [ ] **Ploidy status** — Section 4 states plainly whether the reported
      calls are diploid-default or ploidy-corrected (Script 4), and Section
      5's heterozygosity numbers are consistent with that answer.
- [ ] **Interpretation** — Section 8 is filled in by hand, referencing the
      actual mapping/coverage, filter, ploidy, and Ti/Tv findings above —
      not left as the auto-generated prompt text.
