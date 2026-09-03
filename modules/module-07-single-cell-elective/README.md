# Module 7 (Elective) — Single-Cell Transcriptomics

**Status: built.**

**Curriculum module:** [Module 7 (Elective) — Single-Cell Transcriptomics](../../docs/curriculum/modular_bioinformatics_curriculum.md).

## Guiding question

Given single-cell (10x Genomics-style) reads, how do I get to a trustworthy gene-by-cell count matrix, and which cells/genes should I exclude before any downstream biology?

## Learning outcomes

- Understand barcode/UMI-aware alignment and why it differs fundamentally from bulk RNA-seq alignment (`master_scrnaseq_pipeline.sh`).
- Apply single-cell-specific QC: doublet detection, ambient RNA, and mitochondrial-percentage filtering (Scripts 1-2).
- Cluster cells and identify cluster-defining marker genes (Scripts 3-4), with an interpretation of what the clusters represent biologically left explicitly to the analyst (Script 5).

## Portfolio artifact (per curriculum)

A QC'd, clustered single-cell dataset with marker-gene identification and an interpretation of what the clusters represent biologically. Script 5 below auto-assembles this as `SINGLECELL_REPORT.md` — sections 1-4 are filled in automatically from Scripts 1-4's own output; section 5 (cell-type interpretation) is intentionally left blank, printed alongside a caution that a cluster is a statement about transcriptional similarity, not a confirmed cell type, as unconditional boilerplate rather than optional filler — mirroring Module 2's `DE_REPORT.md`, Module 3's `ASSEMBLY_REPORT.md`, Module 4's `VARIANT_REPORT.md`, Module 5's `PHYLO_REPORT.md`, and Module 6's `MICROBIOME_REPORT.md`.

## Relationship to existing repo content

This module's own draft README (this session's starting point, provided by the user) described itself as **"built (upstream quantification only)"**: `master_scrnaseq_pipeline.sh` — a pre-existing, unchanged, brand-new pipeline (not a correction of an existing script; the domain coverage audit that flagged single-cell RNA-seq as entirely missing from the original 22-script collection led to this being built fresh, the same category as Module 6's amplicon pipeline) — stops at the gene-by-cell count matrix and explicitly does not perform the QC/clustering/marker-identification steps the curriculum's learning outcomes call for. That draft's own "Suggested next step" named exactly what was missing: *"A companion `scrnaseq_seurat_analysis.R` (or Scanpy equivalent) covering QC filtering, clustering, and marker gene identification, to complete this module's portfolio artifact end-to-end."* Scripts 1-4 below are that companion, split into one-idea-per-script Tier B scaffolding (QC filtering → doublet detection → clustering → marker identification) rather than one monolithic script, matching this repository's established Tier B style; Script 5 auto-assembles the portfolio report.

**A disclosed finding about this module's source material, not silently worked around:** two draft scripts handed off for this module (`PRJNA450892_GSE114725_scRNA_processing.sh`, and a near-duplicate misfiled as `PRJNA482620_TCGA_BRCA_alignment.sh` — its filename names a different study than its actual content) both hardcode `SRR7227261` and claim it is BioProject PRJNA450892 / GSE114725, "Single-cell map of immune phenotypes in breast cancer." Verified via NCBI eutils, both claims are wrong, independently: (1) `SRR7227261` actually belongs to BioProject PRJNA438778 ("TCM visualizes trajectories and cell populations from single cell data"), an unrelated ~129K-spot MiSeq dataset far too small to be a real 10x run; (2) GSE114725 itself (confirmed via its own GEO record) is really SRA study SRP148597 / BioProject PRJNA472383, and — more importantly — its actual protocol is **inDrop v2, not 10x Genomics**. Running `master_scrnaseq_pipeline.sh`'s 10x-specific STARsolo logic against inDrop reads would not error; per that script's own methodology notes, a chemistry mismatch is a *silent* failure (near-zero valid cells), exactly the trap a smoke test exists to avoid. This module's worked example (below) uses a different, independently verified real 10x dataset instead. A third file in the same source folder, `PRJNA482620_TCGA_BRCA_downstream.R`, integrates TCGA-BRCA bulk survival analysis (ssGSEA, WGCNA, LASSO Cox, Kaplan-Meier) with single-cell UMAP — genuinely interesting material, but well beyond this module's stated scope (QC → clustering → marker ID); it is not used here and may be worth revisiting for Module 9 (Capstone).

## Worked-example dataset (real, verified via NCBI eutils — same method Modules 3-6 used)

**BioProject PRJNA593571**, "Benchmarking Single-Cell RNA Sequencing Protocols for Cell Atlas Projects" — a deliberate species-mixing ("barnyard") design: human PBMC (60%) + mouse colon (30%) + HEK293T/NIH3T3/MDCK cell lines (10% combined), 10x Chromium V3 (confirmed via read structure: R1=28bp = 16bp cell barcode + 12bp UMI; R2=89bp cDNA), Illumina NovaSeq 6000:

| Run | Sample | Note |
|---|---|---|
| `SRR10587809` | 10XV3_AN4491 | Chromium V3 **with** viability sorting |
| `SRR10587810` | 10XV3_AN4492 | Chromium V3 **without** viability sorting |

A genuine two-condition comparison (viability sorting on/off), not two replicates of one condition. `master_scrnaseq_pipeline.sh`'s own default reference is human-only (GRCh38), so the mouse/cell-line fraction of this barnyard library will simply fail to align — expected, not a bug.

**Why a barnyard design is a bonus for this module specifically:** a droplet containing both human and mouse transcripts is an unambiguous, ground-truth doublet — the classic experimental validation method for any doublet-detection algorithm. Script 2 does **not** use this ground truth (Tier A aligns human-only, so a droplet's mouse content is simply absent from the count matrix Script 2 ever sees) — see Script 2's header. Extending Tier A to a combined human+mouse reference specifically to cross-validate scDblFinder's calls against this dataset's real barnyard ground truth is a genuine, not-yet-built extension — see "still to develop" below.

## Scripts in this module

**Tier A** (unchanged pipeline + a smoke-test wrapper supplying a real, verified `srr_list.txt` cohort in place of `master_scrnaseq_pipeline.sh`'s own hardcoded SRR11092056/SRR11092057 default — already flagged elsewhere in this repository, via NCBI eutils, as SARS-CoV-2 RNA-seq, not single-cell 10x data of any kind):

| File | Purpose |
|---|---|
| [`scripts/master_scrnaseq_pipeline.sh`](scripts/master_scrnaseq_pipeline.sh) | Generic 10x Genomics-style STARsolo pipeline: SRA download → FastQC → barcode/UMI-aware alignment + per-cell quantification (`CHEMISTRY=10x_v3` or `10x_v2`). Stops at the gene-by-cell count matrix. |
| [`scripts/smoke_test_prjna593571.sh`](scripts/smoke_test_prjna593571.sh) | Runs `master_scrnaseq_pipeline.sh` against the verified 2-sample barnyard cohort above. |

**Tier B** (new, one idea per script):

| # | File | New idea introduced |
|---|---|---|
| 1 | [`scripts/1_Load_And_QC_Filter_Cells.R`](scripts/1_Load_And_QC_Filter_Cells.R) | Load STARsolo's filtered matrix into Seurat, compute per-cell QC metrics (genes detected, UMIs, mitochondrial %), and filter out low-quality droplets — a real cell that was dying or ruptured during dissociation, not an empty droplet (STARsolo's own cell-calling already handled that). |
| 2 | [`scripts/2_Detect_And_Remove_Doublets.R`](scripts/2_Detect_And_Remove_Doublets.R) | Run scDblFinder (the field-standard tool, not a hand-rolled high-count heuristic) to flag and remove droplets that captured two cells at once — a failure mode Script 1's own QC metrics cannot catch, since a doublet often looks like a healthy, high-count cell. |
| 3 | [`scripts/3_Cluster_Cells.R`](scripts/3_Cluster_Cells.R) | Normalize → find variable features → PCA → graph-based clustering → UMAP, with an explicit note that clustering resolution is a knob for the biological question being asked, not a fact about the data. |
| 4 | [`scripts/4_Identify_Cluster_Markers.R`](scripts/4_Identify_Cluster_Markers.R) | Compute differentially-expressed marker genes per cluster (`FindAllMarkers`, `p_val_adj`-ranked) — evidence for a cell-type assignment, deliberately stopping short of assigning one. |
| 5 | [`scripts/5_Assemble_SingleCell_Report.py`](scripts/5_Assemble_SingleCell_Report.py) | Auto-assemble `SINGLECELL_REPORT.md` from Scripts 1-4, with the objective sections filled in automatically and the cell-type interpretation section left blank, printed alongside the cluster-vs-cell-type caution as unconditional boilerplate. |

### Verification performed

`master_scrnaseq_pipeline.sh` and `smoke_test_prjna593571.sh` pass `bash -n`. This session's environment has no `STAR`, `Seurat`/`scDblFinder`, `Rscript`, or real Python interpreter installed (only a Windows Store stub alias), so the four R scripts and `5_Assemble_SingleCell_Report.py` could not be syntax-checked with `Rscript -e "parse(...)"` / `python3 -m py_compile`, or functionally run against real or synthetic fixtures — they received careful manual line-by-line review only, matching the disclosure Modules 4-6 already made in this same environment. `master_scrnaseq_pipeline.sh` itself is pre-existing delivered material, not authored in this repository. Completing the functional testing this module still lacks (a small synthetic 10x-format matrix run through Scripts 1-5 end-to-end) on a machine with the actual toolchain is the top item in "still to develop" below.

## How to run

1. `cd scripts/`
2. Obtain a 10x v3 barcode whitelist (see `master_scrnaseq_pipeline.sh`'s header for where) and set `WHITELIST_PATH`.
3. `sbatch smoke_test_prjna593571.sh` → `alignment/<SRR>/Solo.out/Gene/filtered/{matrix.mtx,features.tsv,barcodes.tsv}`
4. `Rscript 1_Load_And_QC_Filter_Cells.R alignment/SRR10587809/Solo.out/Gene/filtered/ .` → `qc_filtered.rds`, `qc_summary.tsv`, `qc_metrics_violin.pdf`
5. `Rscript 2_Detect_And_Remove_Doublets.R qc_filtered.rds .` → `doublets_removed.rds`, `doublet_summary.tsv`
6. `Rscript 3_Cluster_Cells.R doublets_removed.rds .` → `clustered.rds`, `umap_by_cluster.pdf`, `cluster_sizes.tsv`
7. `Rscript 4_Identify_Cluster_Markers.R clustered.rds .` → `all_markers.tsv`, `top_markers_per_cluster.tsv`
8. `python3 5_Assemble_SingleCell_Report.py --qc-tsv qc_summary.tsv --doublet-tsv doublet_summary.tsv --cluster-sizes-tsv cluster_sizes.tsv --top-markers-tsv top_markers_per_cluster.tsv --out SINGLECELL_REPORT.md`
9. Open `SINGLECELL_REPORT.md` and fill in Section 5 by hand, alongside its printed caveat.

See [`RUNBOOK.md`](RUNBOOK.md) for the full step-by-step execution order, required tools, expected outputs, and troubleshooting.

## Teaching materials

- Lecture deck: not yet provided for this module — still to develop.
- Runbook: [`RUNBOOK.md`](RUNBOOK.md).
- Still to develop per [Section 5 of the curriculum](../../docs/curriculum/modular_bioinformatics_curriculum.md): the outstanding functional (not just syntax/manual-review) testing noted above, a guided student workbook, a dataset card/provenance record for PRJNA593571, and a transfer-task assessment. Two genuine extensions beyond this module's current scope: (1) a combined human+mouse Tier A alignment to cross-validate Script 2's scDblFinder calls against this dataset's real barnyard ground truth (see above); (2) `PRJNA482620_TCGA_BRCA_downstream.R`'s bulk-vs-single-cell integrative survival analysis, out of scope here but a plausible fit for Module 9 (Capstone).
