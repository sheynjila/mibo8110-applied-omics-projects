# Module 7 (Elective) — Single-Cell Transcriptomics: Runbook

This is an **operational** document — step-by-step execution instructions,
expected outputs, and troubleshooting. For the conceptual "why," see this
module's [`README.md`](README.md). For the full curriculum context, see
[`docs/curriculum/modular_bioinformatics_curriculum.md`](../../docs/curriculum/modular_bioinformatics_curriculum.md).

## Prerequisites

- SRA-Toolkit, FastQC, STAR, MultiQC — loaded automatically via
  [`scripts/load_modules.sh`](scripts/load_modules.sh) (pinned to
  SRA-Toolkit 3.2.0, FastQC 0.11.9, MultiQC 1.28; STAR's pin is an
  unverified guess, not a tested value, so expect it to fall back to
  `module avail STAR` discovery on most clusters — see the root
  [`RUNBOOK.md`](../../RUNBOOK.md) §4, and set `STAR_MODULE=<exact-name>`
  if discovery finds more than one candidate on yours).
- A 10x cell-barcode whitelist file matching `CHEMISTRY` (`10x_v3` needs a
  16bp-CB/12bp-UMI whitelist — obtainable from Cell Ranger's reference
  bundles or 10x Genomics support). Set `WHITELIST_PATH` before running.
  `master_scrnaseq_pipeline.sh` refuses to run without one (Phase 0).
- R with `Seurat`, `SingleCellExperiment`, `scDblFinder`, and `ggplot2` —
  used by Scripts 1-4.
- Python 3, standard library only — Script 5 uses no third-party packages,
  the same choice Modules 5-6 made.

## Execution order

| Step | Command | Expected result |
|---|---|---|
| 0 | Set `WHITELIST_PATH` to a real 10x v3 whitelist; `sbatch smoke_test_prjna593571.sh` | `$WORKDIR/alignment/<SRR>/Solo.out/Gene/{raw,filtered}/`, `$WORKDIR/final_scrnaseq_report.html` (MultiQC) |
| 1 | `Rscript 1_Load_And_QC_Filter_Cells.R $WORKDIR/alignment/SRR10587809/Solo.out/Gene/filtered/ .` | Console before/after QC report; `qc_filtered.rds`, `qc_summary.tsv`, `qc_metrics_violin.pdf` |
| 2 | `Rscript 2_Detect_And_Remove_Doublets.R qc_filtered.rds .` | Console doublet-rate report; `doublets_removed.rds`, `doublet_summary.tsv` |
| 3 | `Rscript 3_Cluster_Cells.R doublets_removed.rds .` | Console cluster-count/size report; `clustered.rds`, `umap_by_cluster.pdf`, `cluster_sizes.tsv` |
| 4 | `Rscript 4_Identify_Cluster_Markers.R clustered.rds .` | Console top-marker report; `all_markers.tsv`, `top_markers_per_cluster.tsv` |
| 5 | `python3 5_Assemble_SingleCell_Report.py --qc-tsv qc_summary.tsv --doublet-tsv doublet_summary.tsv --cluster-sizes-tsv cluster_sizes.tsv --top-markers-tsv top_markers_per_cluster.tsv` | `SINGLECELL_REPORT.md`, sections 1-4 auto-filled, section 5 blank |
| 6 | Fill in Section 5 (cell-type interpretation) of `SINGLECELL_REPORT.md` by hand, alongside its printed caveat | Completed portfolio artifact |

`$WORKDIR` defaults to `/scratch/$(whoami)/master_scrnaseq_pipeline` (set
inside `master_scrnaseq_pipeline.sh`). Steps 1-5 above process one sample
(`SRR10587809`) end-to-end; repeat Steps 1-5 with `SRR10587810`'s matrix
path to process the second, or merge the two Seurat objects after Script 1
if you want a single combined clustering — this module's scripts process
one sample per invocation and do not merge automatically.

## Expected outputs, by step

- **Step 0:** per-sample STARsolo `raw`/`filtered` matrices, MultiQC report.
- **Step 1:** `qc_filtered.rds`, `qc_summary.tsv`, `qc_metrics_violin.pdf`.
- **Step 2:** `doublets_removed.rds`, `doublet_summary.tsv`.
- **Step 3:** `clustered.rds`, `umap_by_cluster.pdf`, `cluster_sizes.tsv`.
- **Step 4:** `all_markers.tsv`, `top_markers_per_cluster.tsv`.
- **Step 5:** `SINGLECELL_REPORT.md`.

## Troubleshooting

| Symptom | Likely cause | Fix |
|---|---|---|
| Step 0 fails at `CRITICAL ERROR: Barcode whitelist not found` | `WHITELIST_PATH` is unset or points at a nonexistent file (Phase 0's deliberate pre-flight check) | Download the correct whitelist for `CHEMISTRY` (10x v3: 16bp CB + 12bp UMI) and set `WHITELIST_PATH` before re-running |
| Step 0 completes but STARsolo's `filtered/` matrix has near-zero cells/barcodes | A chemistry mismatch — `CHEMISTRY`/`CB_LEN`/`UMI_LEN` don't match how the reads were actually generated (this is a **silent** failure per the script's own methodology notes, not an error) | Confirm the study's actual protocol from its methods/SRA metadata before trusting `CHEMISTRY=10x_v3` on a new dataset — see this module's README for a real example (PRJNA450892/SRR7227261) where the source material's claimed protocol was wrong |
| Step 1 warns `fewer than half the input cells passed QC` | Genuinely poor-quality sample, or `MIN_FEATURES`/`MIN_COUNTS`/`MAX_PCT_MT` don't fit this tissue | Inspect `qc_metrics_violin.pdf` before assuming either explanation; adjust the thresholds at the top of the script only with a stated reason, not to force a target retention percentage |
| Step 2 reports a doublet rate above 20% | High cell-loading concentration for this run, or ambient-RNA-heavy droplets being misclassified as doublets | Check this sample's actual 10x loading target before treating the number as expected; do not silently raise scDblFinder's classification threshold to force a lower rate |
| Step 3 produces 1 cluster (or a cluster count that looks too low/high) | `CLUSTER_RESOLUTION=0.8` does not fit this dataset's actual heterogeneity | Re-run with a different `CLUSTER_RESOLUTION` (Script 3's header explains resolution is a knob for the biological question, not a fact to discover) |
| Script 4 warns some cluster(s) have no significant marker gene | That cluster is not well-separated from its neighbors at the current resolution | Consider a lower `CLUSTER_RESOLUTION` in Script 3 and re-run Scripts 3-4, rather than reporting an empty marker list as-is |
| Script 5's `SINGLECELL_REPORT.md` shows "Not available" for a section | The corresponding earlier script's output file was not passed via its `--*-tsv` flag | Re-run the missing step, then re-run Script 5 with the correct file path |
| Script 5's Section 5 is blank | Expected — intentionally left for the analyst, alongside the printed caveat | Fill it in referencing Sections 1-4's actual numbers and `top_markers_per_cluster.tsv`, and do not remove the caveat text above it |

## Completion checklist (portfolio artifact)

The curriculum's Module 7 assessment asks for: *"A QC'd, clustered
single-cell dataset with marker-gene identification and an interpretation
of what the clusters represent biologically."* Confirm each is present
before considering the module complete:

- [ ] **QC'd dataset** — Section 1 of `SINGLECELL_REPORT.md` reflects a real
      Script 1 run, and `qc_metrics_violin.pdf` was actually inspected, not
      just generated.
- [ ] **Doublet-filtered dataset** — Section 2 reflects a real Script 2 run,
      including the doublet rate, not just a "done" checkbox.
- [ ] **Clustered dataset** — Section 3 reflects a real Script 3 run,
      including any clusters with fewer than 10 cells, not just the total
      count.
- [ ] **Marker-gene identification** — Section 4 reflects a real Script 4
      run, including any cluster with no significant markers.
- [ ] **Written biological interpretation** — Section 5 is filled in by
      hand, explicitly distinguishing what the marker genes support (a
      working cell-type hypothesis) from what clustering alone does not
      establish (a confirmed cell-type identity, or that a cluster boundary
      is not a technical artifact) — the caveat text above it is not
      removed or paraphrased away.
