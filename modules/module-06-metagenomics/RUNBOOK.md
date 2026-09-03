# Module 6 — Metagenomics and Microbiome Profiles: Runbook

This is an **operational** document — step-by-step execution instructions,
expected outputs, and troubleshooting. For the conceptual "why," see this
module's [`README.md`](README.md). For the full curriculum context, see
[`docs/curriculum/modular_bioinformatics_curriculum.md`](../../docs/curriculum/modular_bioinformatics_curriculum.md).

## Prerequisites

**Shotgun route:**
- SRA-Toolkit 3.2.0 (`prefetch`, `fasterq-dump`), fastp 0.23.4, Kraken2
  2.1.2, Bracken 2.8, MultiQC 1.28 — module versions as pinned in
  `master_microbiome_pipeline.sh`'s own `module load` lines.
- A pre-built Kraken2 database (>50GB) — path configured via `KRAKEN_DB`
  inside `master_microbiome_pipeline.sh`.
- Bowtie2 2.5.4 and a host reference genome FASTA — only if `HOST_REFERENCE`
  is set to exercise the optional host-depletion step.
- Bash — used by Script 1.
- Python 3, standard library only — Scripts 2-5 use no third-party
  packages, the same choice Module 5 made.
- R with `ggplot2`, `dplyr`, `tidyr`, `pheatmap`, `RColorBrewer`, `scales`
  — used by `04_visualize_microbiome.R` (optional visualization step, not
  required for the report itself).

**Amplicon route:**
- SRA-Toolkit 3.2.0, FastQC 0.11.9, cutadapt 4.9, MultiQC 1.28 — used by
  `master_amplicon_pipeline.sh`.
- R with the `dada2` and `ggplot2` packages — used by
  `amplicon_dada2_analysis.R`.
- A DADA2-formatted SILVA (16S) or UNITE (ITS) reference FASTA, downloaded
  separately (see `amplicon_dada2_analysis.R`'s header for the source URL).
  This is read from the exact path in that script's own `REF_DB_PATH`
  variable (a placeholder `/path/to/...` by default) — edit that variable
  to point at wherever you actually put the file, rather than relying on a
  working-directory filename convention.

## Execution order — shotgun route

| Step | Command | Expected result |
|---|---|---|
| 0 | `sbatch smoke_test_shotgun_prjna1520859.sh` (optionally `export HOST_REFERENCE=...` first) | `$WORKDIR/bracken_output/ALL_SAMPLES_MICROBIOME.tsv`, `$WORKDIR/kraken_output/<SRR>_report.txt`, `$WORKDIR/final_microbiome_report.html` (MultiQC) |
| 1 | `bash 1_Validate_Classification_Inputs.sh $WORKDIR` | Console PASS/FAIL/WARN report per sample |
| 2 | `python3 2_Assess_Classification_Yield_And_Host_Fraction.py $WORKDIR .` | Console per-sample yield/host report; `classification_yield_summary.tsv` |
| 3 | `python3 3_Cohort_Diversity_And_Dominance_Summary.py $WORKDIR/bracken_output/ALL_SAMPLES_MICROBIOME.tsv .` | Console per-sample diversity/dominance report; `cohort_diversity_summary.tsv` |
| 4 | `python3 4_Detection_Limit_And_Database_Caveats.py $WORKDIR/bracken_output/ALL_SAMPLES_MICROBIOME.tsv classification_yield_summary.tsv .` | Console per-taxon caveat report; `database_detection_caveats.tsv` |
| 5 | `python3 5_Assemble_Microbiome_Report.py --yield-tsv classification_yield_summary.tsv --diversity-tsv cohort_diversity_summary.tsv --caveats-tsv database_detection_caveats.tsv` | `MICROBIOME_REPORT.md`, sections 1-4 auto-filled, section 5 blank |
| 6 | Fill in Section 5 (interpretation) of `MICROBIOME_REPORT.md` by hand, alongside its printed caveat | Completed portfolio artifact |
| 7 (optional) | `Rscript 04_visualize_microbiome.R` (run where `ALL_SAMPLES_MICROBIOME.tsv` or `bracken_output/ALL_SAMPLES_MICROBIOME.tsv` is reachable) | `plots/01_microbiome_relative_abundance.png`, `plots/02_microbiome_taxa_heatmap.png` |

`$WORKDIR` defaults to `/scratch/$(whoami)/master_microbiome_pipeline`
(set inside `master_microbiome_pipeline.sh`).

## Execution order — amplicon route

| Step | Command | Expected result |
|---|---|---|
| 0 | `export AMPLICON_TYPE=16S` (or `ITS`); `sbatch smoke_test_amplicon_prjna917645.sh` | `$WORKDIR/trimmed_reads/<SRR>_R[12]_trimmed.fastq`, `$WORKDIR/qc/<SRR>_cutadapt.log`, `$WORKDIR/qc/<SRR>_1_fastqc.html` (+ trimmed-read FastQC), `$WORKDIR/final_amplicon_<AMPLICON_TYPE>_report.html` (MultiQC) |
| 1 | Download the matching SILVA or UNITE reference set; edit `REF_DB_PATH` near the top of `amplicon_dada2_analysis.R` to point at it; set `AMPLICON_TYPE` in the same script to match Step 0 | Reference FASTA reachable at the exact path in `REF_DB_PATH` |
| 2 | `Rscript amplicon_dada2_analysis.R` (run from `$WORKDIR`, or with `trimmed_reads/` reachable from the working directory) | `quality_profiles.pdf`, `ASV_table.tsv`, `taxonomy.tsv`, `ALL_SAMPLES_AMPLICON.tsv`, `amplicon_composition.pdf` |

`$WORKDIR` defaults to `/scratch/$(whoami)/master_amplicon_pipeline_<AMPLICON_TYPE>`.
This route has no Tier B walkthrough of its own — see README "Relationship
to existing repo content" for why — so its own console output (cutadapt's
per-sample pass rate, DADA2's `filterAndTrim`/`removeBimeraDenovo`
retained-read counts), the MultiQC report, and `quality_profiles.pdf` /
`amplicon_composition.pdf` are this route's QC record. `ALL_SAMPLES_AMPLICON.tsv`
uses the same `Sample`/`Taxon`/`Estimated_Reads`/`Fraction_Total` schema as
the shotgun route's `ALL_SAMPLES_MICROBIOME.tsv` (by the R script's own
design — see its Phase E comment), so Tier B Scripts 3 and 4 above can be
pointed at it directly if you want a diversity/dominance or
database-detection-caveat pass on an amplicon cohort; Scripts 1, 2, and 5
are Kraken2/Bracken-specific and do not apply to amplicon output.

## Expected outputs, by step (shotgun route)

- **Step 0:** `ALL_SAMPLES_MICROBIOME.tsv`, per-sample Kraken2/Bracken files, MultiQC report.
- **Step 1:** console validation report only (no file written).
- **Step 2:** `classification_yield_summary.tsv`.
- **Step 3:** `cohort_diversity_summary.tsv`.
- **Step 4:** `database_detection_caveats.tsv`.
- **Step 5:** `MICROBIOME_REPORT.md`.
- **Step 7 (optional):** `plots/01_microbiome_relative_abundance.png`, `plots/02_microbiome_taxa_heatmap.png`.

## Troubleshooting

| Symptom | Likely cause | Fix |
|---|---|---|
| Step 0 fails at `CRITICAL ERROR: Kraken2 Database not found` | `KRAKEN_DB` inside `master_microbiome_pipeline.sh` does not point at a real, pre-built database on this HPC system | Check your HPC documentation for the shared Kraken2 database path and edit `KRAKEN_DB` before re-running |
| Step 0's shotgun run only produces a single `_1.fastq` per sample, and `fastp`/`kraken2 --paired` then fails | You pointed the shotgun pipeline at single-end SRA data (e.g. PRJNA917645's own shotgun arm — see README's disclosed dataset/pipeline mismatch) | Use `smoke_test_shotgun_prjna1520859.sh`'s verified paired-end cohort, or any other genuinely paired-end shotgun accession, instead |
| Step 1 fails `no *_report.txt files found` | Step 0 has not been run yet, or `$WORKDIR` does not match `master_microbiome_pipeline.sh`'s own `WORKDIR` variable | Re-run Step 0 first; pass the exact same `$WORKDIR` to Script 1 |
| Step 1 reports `[FAIL] ... Bracken output not found` for a sample whose Kraken2 report exists | That sample is in `master_microbiome_pipeline.sh`'s `FAILED_SAMPLES` list from a Bracken-stage error (e.g. Bracken's minimum-read threshold for that sample's read length) | Check that sample's console output/logs from Step 0; re-run just that sample or accept the smaller cohort and note it in `MICROBIOME_REPORT.md` Section 5 |
| Script 2 flags many/most samples `LOW_YIELD` | The Kraken2 database is missing common taxa for this sample type (e.g. a human-only database run against a soil sample), or the sample itself is low-quality/low-complexity | Confirm the database matches the expected community before concluding the sample itself is the problem; this is exactly the "artifact vs. real biology" split this module's guiding question asks about |
| Script 3 flags a sample `DOMINANT` unexpectedly | A single taxon holding most of that sample's reads — could be a real bloom/infection, a genuinely low-complexity sample, or reference-database bias toward one well-represented genome | Cross-check that taxon's identity against the sample's expected biology and against Script 2's classification yield for that sample before treating the dominance as either "real" or "artifact" |
| Script 4 reports many taxa `WEAK_ABSENCE_EVIDENCE` | Several cohort samples have low classification yield (Script 2), so their "missing" taxa rows carry little evidentiary weight | Address the underlying low-yield samples (see the `LOW_YIELD` row above) rather than treating Script 4's output as itself the problem |
| Script 5's `MICROBIOME_REPORT.md` shows "Not available" for a section | The corresponding earlier script's output file was not passed via its `--*-tsv` flag | Re-run the missing step, then re-run Script 5 with the correct file path |
| Script 5's Section 5 is blank | Expected — intentionally left for the analyst, alongside the printed caveat | Fill it in referencing Sections 1-4's actual numbers, and do not remove the caveat text above it |
| Amplicon Step 0's cutadapt log shows a low "Pairs written (passing filters)" percentage | `AMPLICON_TYPE` does not match the actual primers used in the wet lab, or the reads are shotgun (no primer at all) rather than amplicon | Confirm the correct primer pair for this dataset; do not lower `--minimum-length` or drop `--discard-untrimmed` to force a higher pass rate — that would let off-target reads back into DADA2 |
| `amplicon_dada2_analysis.R` stops at `REF_DB_PATH not found at '/path/to/...'` | `REF_DB_PATH` was left at its placeholder value, or the SILVA/UNITE reference FASTA was not actually downloaded to the path it names | Edit `REF_DB_PATH` near the top of the script to the real path of a downloaded, DADA2-formatted reference FASTA (see the script's header for the source URL) and re-run |
| `amplicon_dada2_analysis.R` stops at `No matched forward/reverse trimmed FASTQ pairs found` | `READS_DIR` (`"trimmed_reads"`) isn't reachable from the directory you launched `Rscript` in, or `master_amplicon_pipeline.sh` did not complete for any sample | Run `Rscript amplicon_dada2_analysis.R` from `$WORKDIR` itself, or copy/symlink `trimmed_reads/` next to it; check Step 0's `FAILED_SAMPLES` output if the directory is unexpectedly empty |

## Completion checklist (portfolio artifact)

The curriculum's Module 6 assessment asks for: *"A microbiome profile
report: method chosen and why, QC/contamination handling, abundance
results, and an interpretation that acknowledges database and
detection-limit caveats."* Confirm each is present before considering the
module complete:

- [ ] **Method chosen and why** — Section 1 of `MICROBIOME_REPORT.md`
      states which route (shotgun or amplicon) was run and reflects the
      actual tradeoff for this sample set, not just the boilerplate text.
- [ ] **QC/contamination handling** — Section 2 reflects a real Script 1 +
      Script 2 run against your actual cohort, including any `LOW_YIELD`
      flags and (if used) the actual host-depletion percentages.
- [ ] **Abundance results** — Section 3 reflects a real Script 3 run,
      including any `DOMINANT` flags, not just agreeing/clean samples.
- [ ] **Database and detection-limit caveats** — Section 4 reflects a real
      Script 4 run against your actual cohort table, and the unconditional
      caveat text is not removed or paraphrased away.
- [ ] **Written interpretation** — Section 5 is filled in by hand,
      explicitly weighing which flagged samples/taxa from Sections 2-4
      change how much the abundance results in Section 3 can actually
      support.
