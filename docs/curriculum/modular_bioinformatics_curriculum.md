# Modular Bioinformatics — Applied Omics Curriculum

*From reproducible computing foundations to portfolio-ready biological analysis.*

This is the canonical curriculum reference for MIBO8110. It defines the
learning outcomes, portfolio artifacts, and assessment rubric that every
module's own [README.md](../../README.md) and `RUNBOOK.md` implement in
code. If a module's README says "per curriculum," this is the document it
means.

**Curriculum purpose:** a competency-based, modular program that teaches
students to frame biological questions, select and defend analytical
decisions, execute reproducible workflows, evaluate quality, and
communicate limitations.

| Program feature | Recommended design |
|---|---|
| Audience | Advanced undergraduate, graduate, or laboratory trainee with introductory molecular biology. |
| Delivery | Blended lecture, demonstration, guided laboratory, independent practice, and capstone. |
| Core pathway | Foundations, raw-read QC, bulk RNA-seq, assembly, variant calling, phylogenomics, and metagenomics. |
| Electives | Single-cell transcriptomics and chromatin/regulatory genomics. |
| Evidence of mastery | Reproducible code, QC evidence, interpretation, technical documentation, and oral or written communication. |

## 1. Program Overview

This curriculum develops practical bioinformatics competence through
progressively more independent analysis. It deliberately treats software as
an implementation choice rather than the curriculum's organizing principle.
Students learn to connect biological questions, data structure, statistical
assumptions, quality evidence, and defensible interpretation.

**Design principle:** every module follows the same learning cycle: frame
the question, inspect the data, choose an approach, execute the workflow,
evaluate quality, interpret results, document limitations, and reproduce
the analysis.

### Program goals

- Explain how computation and data mining support hypothesis-driven and hypothesis-generating life-science research.
- Apply statistical reasoning to high-dimensional biological data, including multiple testing and uncertainty.
- Retrieve, organize, validate, and document common biological data types.
- Use command-line, R, or Python tools to investigate biological questions.
- Create reproducible environments, scripts, reports, and provenance records.
- Critically evaluate analytical outputs, database dependence, bias, and study limitations.
- Communicate results to technical and biological audiences.

### Entry prerequisites

| Area | Minimum preparation | Bridging support |
|---|---|---|
| Biology | Central dogma, genomes, genes, variation, and experimental controls. | Pre-course primer and terminology check. |
| Statistics | Distributions, variation, hypothesis tests, confidence, and basic visualization. | Short statistics bootcamp with omics examples. |
| Computing | No advanced programming required; willingness to work with scripts is essential. | Module 0 shell, R/Python, Git, and environments. |
| Research practice | Ability to read a methods section and interpret a figure. | Guided article-reading worksheet. |

### Suggested delivery model

The core sequence can be delivered as a semester course, an intensive
institute, or a series of laboratory rotations. In shorter formats, preserve
Module 0 and the common QC/reproducibility expectations, then let students
select one analytical specialization. Time allocations should be adjusted
to learner preparation and available computing resources.

## 2. Curriculum Architecture

| Stage | Module | Primary portfolio artifact |
|---|---|---|
| Foundation | 0. Computing, data, statistics, and reproducibility | Version-controlled reproducible mini-pipeline |
| Common core | 1. Raw reads and quality control | QC report with justified preprocessing decisions |
| Core analysis | 2. Bulk RNA-seq | Differential-expression report |
| Core analysis | 3. Genome assembly | Evaluated pathogen or microbial assembly |
| Core analysis | 4. Variant calling | Filtered VCF and technical validation summary |
| Core analysis | 5. Phylogenomics | Supported tree and cautious epidemiological interpretation |
| Core analysis | 6. Metagenomics | Community profile and diversity report |
| Elective | 7. Single-cell transcriptomics | Annotated cell atlas with uncertainty statement |
| Elective | 8. Chromatin and regulatory genomics | Peak/accessibility analysis and genome-browser evidence |
| Synthesis | 9. Integrated capstone | Portfolio-ready reproducibility package |

**Recurring assessment pattern:** diagnostic pre-check → instructor
demonstration → guided laboratory → independent transfer task →
interpretation and troubleshooting questions → reproducibility audit →
portfolio artifact and reflection.

---

## MODULE 0 — Computing, Data, Statistics, and Reproducibility

**Guiding question:** How can an analysis be made understandable,
repeatable, portable, and auditable?

**Learning outcomes**
- Navigate and manipulate files safely in Linux.
- Explain and validate FASTA, FASTQ, SAM/BAM, BED, GFF/GTF, and VCF structures.
- Use pipes, redirection, shell variables, loops, and logs in simple workflows.
- Create and document Conda/Mamba or containerized environments.
- Use Git for version history and collaboration.
- Select sensible statistical summaries and visualizations.
- Record data provenance, parameters, software versions, and computing resources.

**Core content and tools**

| Content | Reference implementation | Decision emphasis |
|---|---|---|
| Linux and file operations | Bash, coreutils, grep, awk, sed | Safety, paths, compression, streams, checksums |
| Programming | R or Python notebooks/scripts | Readable code, functions, validation, errors |
| Versioning | Git and a repository host | Commits, README, ignore rules, collaboration |
| Environments | Conda/Mamba; optional containers | Versions, portability, dependency conflicts |
| Compute | WSL, Linux server, cluster, or cloud | CPU, memory, storage, job scheduling, cost awareness |
| Responsible practice | Licensing, privacy, data governance | Permitted use, de-identification, human oversight |

**Assessment:** students submit a repository containing a small pipeline,
environment specification, README, input checks, log, outputs, and a
reproducibility reflection.

> This repo's [Module 0](../../modules/module-00-foundations-reproducibility/README.md)
> adds one more explicit competency beyond the original curriculum text:
> robust HPC environment-module loading (`load_modules.sh`,
> [Script 11](../../modules/module-00-foundations-reproducibility/scripts/11_Robust_HPC_Module_Loading.sh)) —
> because "portable" and "auditable" both break the moment a pinned
> `module load <exact-build>` doesn't exist on a different cluster.

## MODULE 1 — Raw Sequencing Reads and Quality Control

**Guiding question:** Are the reads suitable for the intended biological
analysis, and what preprocessing is justified?

**Learning outcomes**
- Explain FASTQ records and Phred-quality concepts.
- Retrieve a selected public run and verify file integrity.
- Interpret per-base quality, adapter, duplication, length, and composition summaries.
- Distinguish diagnostic QC from automatic trimming.
- Process paired-end reads without breaking read pairing.
- Document preprocessing decisions and their likely consequences.

**Workflow:** define the downstream question and expected sequencing
design → retrieve reads with SRA Toolkit or an approved archive method →
inspect raw reads with FastQC and summarize batches with MultiQC → use
fastp only when evidence supports filtering or adapter removal → repeat QC,
compare before/after metrics, and preserve audit logs.

**Misconceptions to address**
- Higher trimming stringency is not automatically better.
- A passing QC report does not prove biological validity.
- Duplication can reflect technical artifacts or genuine biology.
- Quality thresholds should reflect platform, assay, and downstream purpose.

**Portfolio artifact:** a concise QC report containing the data source,
checksums, diagnostic figures, preprocessing rationale, command log, and
acceptance criteria.

## MODULE 2 — Bulk RNA-seq and Gene Expression

**Guiding question:** Which genes or pathways show evidence of systematic
expression differences after accounting for biological and technical
variation?

**Learning outcomes**
- Validate integer count data and sample metadata.
- Explain biological replication, confounding, covariates, and batch effects.
- Filter low-information features using a documented rule.
- Perform exploratory analysis using transformations, PCA, distances, and sample-level plots.
- Specify designs and contrasts in DESeq2 or edgeR.
- Interpret effect sizes, uncertainty, p-values, and false-discovery-rate adjustment.
- Produce accessible figures and a biologically cautious report.

**Recommended environment:** RStudio or Jupyter with R support. Use
DESeq2 as the primary reference implementation and introduce edgeR as an
alternative. Use ggplot2 and appropriate heatmap/reporting packages for
communication.

| Phase | Required decisions | Evidence |
|---|---|---|
| Design | Replicates, contrasts, covariates, confounding | Design table and rationale |
| Input validation | Counts, identifiers, metadata matching | Validation checks and exclusions |
| Exploration | Transformation, outliers, batch structure | PCA and sample-distance evidence |
| Modeling | Design formula, reference level, contrast | Recorded model and assumptions |
| Inference | Filtering, shrinkage, multiplicity | Effect size and adjusted significance |
| Interpretation | Annotation and pathway context | Claims tied to evidence and limitations |

**Portfolio artifact:** a reproducible report containing metadata
validation, exploratory plots, model specification, differential-expression
results, effect-size visualization, interpretation, and limitations.

## MODULE 3 — Microbial Genome Assembly

**Guiding question:** What genome sequence can be reconstructed from the
reads, and how strong is the evidence for its quality and completeness?

**Learning outcomes**
- Choose among short-read, long-read, and hybrid assembly strategies.
- Run SPAdes for short reads or Flye for long reads with documented parameters.
- Interpret an assembly graph and identify unresolved structures.
- Evaluate contiguity without treating N50 as a complete quality measure.
- Assess completeness, contamination, coverage, and consistency with expected biology.
- Distinguish assembly evaluation from genome annotation.

**Quality framework**

| Dimension | Questions students must answer |
|---|---|
| Contiguity | How fragmented is the assembly, and are large contigs supported? |
| Completeness | Are expected conserved features represented? |
| Correctness | Are there signs of misassembly, inconsistent coverage, or conflicting read support? |
| Contamination | Do taxonomic or compositional signals suggest mixed material? |
| Fitness for purpose | Is the assembly adequate for typing, gene detection, comparative genomics, or closure? |

**Portfolio artifact:** assembly FASTA, environment and commands, assembly
statistics, graph or coverage evidence, quality interpretation, and a clear
fitness-for-purpose conclusion.

## MODULE 4 — Variant Calling

**Guiding question:** Which sequence differences are supported by the
reads, and which filters make the final call set defensible?

**Learning outcomes**
- Select and document an appropriate reference genome.
- Map reads with BWA and evaluate alignment quality.
- Explain coverage, mapping quality, base quality, strand balance, and allele evidence.
- Call variants with BCFtools or GATK using organism-appropriate assumptions.
- Design and justify filters rather than copying thresholds blindly.
- Inspect representative variants and summarize VCF fields.
- Recognize reference bias, ploidy assumptions, and inaccessible regions.

**Required controls and checkpoints**
- Reference identity and version
- Read-group and sample identity checks
- Coverage distribution and low-accessibility regions
- Organism ploidy and expected allele behavior
- Raw versus filtered call counts recorded by the workflow
- Manual inspection of selected high- and low-confidence sites
- Consistent sample and contig naming across files

**Portfolio artifact:** a filtered VCF, callable-region definition, QC
summary, parameter rationale, selected-site inspection evidence, and
limitations statement.

## MODULE 5 — Phylogenomics and Outbreak Context

**Guiding question:** What evolutionary relationships are supported, and
what conclusions are not justified by the tree alone?

**Learning outcomes**
- Construct and validate a comparable sequence or core-variant alignment.
- Identify problematic sites, missingness, and potential recombination.
- Choose a substitution model and infer a tree with IQ-TREE or RAxML-NG.
- Interpret branch support, rooting, scale, and metadata overlays.
- Separate genetic relatedness from direct transmission claims.
- Report uncertainty and sampling limitations.

**Interpretation guardrails**
- A tree depicts relationships under a model; it is not a direct transmission map.
- Sparse or biased sampling can change the apparent context of a cluster.
- Low-quality or recombinant regions may distort inference.
- Epidemiological metadata should complement, not be inferred from, genomic distance alone.
- Branch support and rooting must be shown and explained.

**Portfolio artifact:** alignment, tree file, visualization, model and
support information, metadata statement, and a short outbreak
interpretation with explicit uncertainty.

## MODULE 6 — Metagenomics and Microbiome Profiles

**Guiding question:** Which organisms and functions are represented in a
mixed sample, and how confidently can they be compared across samples?

**Learning outcomes**
- Distinguish amplicon, shotgun profiling, and assembly-based metagenomics.
- Explain negative controls, positive controls, host reads, and contamination risk.
- Generate taxonomic profiles with Kraken2 or MetaPhlAn.
- Describe database dependence and classification uncertainty.
- Import community data into phyloseq or an equivalent framework.
- Calculate and interpret alpha and beta diversity appropriately.
- Explain compositionality and limitations of relative abundance.
- Select a defensible approach to differential abundance.

**Instructional tracks**

| Track | Scope | Suggested tools |
|---|---|---|
| Amplicon | Marker-gene denoising, taxonomy, diversity | DADA2 or mothur; phyloseq |
| Shotgun profiling | Read-level taxonomic and functional profiles | Kraken2/Bracken or MetaPhlAn; optional HUMAnN |
| Assembly-based | Metagenome assembly, bins, evaluation | metaSPAdes/MEGAHIT; binning and quality tools |

**Portfolio artifact:** control-aware community report containing
provenance, classification database/version, QC, abundance and diversity
figures, interpretation, and database/compositional limitations.

## MODULE 7 (elective) — Single-Cell Transcriptomics

**Guiding question:** Which cellular populations and states are supported
by the data without confusing technical structure for biology?

**Learning outcomes**
- Validate cell-by-gene matrices and sample metadata.
- Apply and justify cell- and gene-level QC criteria.
- Recognize empty droplets, doublets, ambient RNA, and mitochondrial-content concerns.
- Normalize, reduce dimensions, cluster, and assess stability.
- Identify markers and annotate cells with evidence and uncertainty.
- Explain batch integration risks and pseudoreplication.
- Distinguish cell-level exploration from sample-level biological inference.

**Reference implementations:** use Seurat in R or Scanpy in Python. Select
one as the teaching implementation so students can focus on concepts, then
provide a translation guide for the alternative ecosystem.

**Portfolio artifact:** annotated cell atlas, QC rationale, marker
evidence, sample-aware interpretation, annotation confidence, and analysis
limitations.

## MODULE 8 (elective) — Chromatin and Regulatory Genomics

**Guiding question:** Which genomic regions show evidence of protein
binding or chromatin accessibility, and how do controls affect
interpretation?

**Learning outcomes**
- Distinguish ChIP-seq, ATAC-seq, and DNA methylation assays.
- Align reads with Bowtie2 and evaluate assay-specific QC.
- Explain controls, duplicate handling, library complexity, and replicate concordance.
- Call enriched regions with MACS2 where appropriate.
- Use BEDTools for interval operations and annotation.
- Visualize normalized tracks in a genome browser.
- Integrate regulatory evidence with expression cautiously.

**Track clarification**

| Assay | Primary signal | Essential teaching emphasis |
|---|---|---|
| ChIP-seq | Protein-associated DNA enrichment | Input/control design, enrichment, replicate agreement |
| ATAC-seq | Accessible chromatin | Fragment pattern, TSS enrichment, accessibility peaks |
| Bisulfite sequencing | DNA methylation state | Conversion efficiency, coverage, methylation context |

**Portfolio artifact:** QC report, peak or methylation output, annotated
regions, browser visualization, replicate/control assessment, and
integrated interpretation.

## MODULE 9 — Integrated Capstone and Professional Portfolio

**Guiding question:** Can the student independently produce a reproducible,
defensible analysis that answers a meaningful biological question?

**Capstone phases:** question and scope proposal → data-management and
computing plan → analysis plan with decision points and acceptance
criteria → pilot analysis using a small subset → full workflow execution →
quality and reproducibility audit → interpretive report and portfolio
packaging → presentation and peer review.

**Required submission package**

| Component | Minimum evidence |
|---|---|
| README | Question, data source, requirements, execution order, outputs |
| Environment | Pinned package/tool versions or container reference |
| Workflow | Scripts or notebooks with parameters and logging |
| Data provenance | Accession, reference versions, checksums, exclusions |
| Quality evidence | Predefined checks and interpretation of failures |
| Results | Machine-readable outputs and publication-ready figures |
| Interpretation | Claims linked to evidence, uncertainty, and limitations |
| Reflection | What changed, what failed, and what would be improved |

---

## 3. Assessment and Rubric Framework

Assessment should reward scientific reasoning and reproducibility, not
merely successful software execution. The following common rubric can be
adapted for each module and the capstone.

| Criterion | Exemplary evidence | Weight |
|---|---|---|
| Question and design | Focused biological question; design addresses replication, controls, and confounding. | 15% |
| Data stewardship | Provenance, identifiers, versions, integrity checks, and exclusions are complete. | 10% |
| Analytical reasoning | Methods and parameters are justified; alternatives and assumptions are considered. | 20% |
| Quality evaluation | Relevant metrics are interpreted and linked to acceptance or remediation decisions. | 15% |
| Reproducibility | Environment, code, order of operations, logs, and outputs can be recreated. | 20% |
| Interpretation | Claims match evidence; uncertainty, bias, and limitations are explicit. | 15% |
| Communication | Figures, captions, organization, and language serve the intended audience. | 5% |

**Recommended formative assessments**
- Predict the effect of changing one parameter before rerunning a workflow.
- Diagnose a deliberately flawed metadata table or file pairing.
- Compare two acceptable tools and defend a selection.
- Interpret a QC failure and propose a bounded remediation plan.
- Review a peer's README for reproducibility gaps.
- Rewrite an overconfident biological claim to match the evidence.

## 4. Computing and Delivery Plan

| Environment | Best use | Constraints to plan for |
|---|---|---|
| Local Windows with WSL | Foundations, scripting, small teaching datasets | Installation support, disk space, memory |
| macOS/Linux workstation | Most small and moderate exercises | Tool architecture and storage differences |
| Galaxy | Guided analysis and accessible workflow learning | Server quotas, queue time, tool versions |
| University cluster | Large raw-read workflows and multi-sample projects | Accounts, scheduler training, file transfer |
| Cloud platform | Scalable or short-term project compute | Cost controls, permissions, data egress |
| Containers | Portable, pinned execution environments | Image storage and platform support |

**Recommended default:** teach concepts with small datasets in WSL, Linux,
or Galaxy. Move full-sized raw-read projects to a managed cluster or cloud
environment. Require the same project structure, provenance, versioning,
and QC evidence in every environment.

> This repo's scripts target a Slurm-managed HPC cluster for every module's
> "full raw-read project" path, and use [`lib/module_loader.sh`](../../lib/module_loader.sh)
> so the pinned tool versions this course was authored against degrade
> gracefully — rather than failing outright — on a cluster with a
> differently-named module tree. See the [root README](../../README.md)
> and [pipeline_orchestration/README.md](../../pipeline_orchestration/README.md)
> for how manual and orchestrated execution both fit into this plan.

**Minimum technical preparation checklist**
- Confirm institutional data-governance requirements.
- Select a supported primary platform and a fallback platform.
- Pre-stage small datasets and reference files.
- Test every environment from a clean learner account.
- Provide an installation check script and troubleshooting guide.
- Define storage quotas, expected runtime, and cleanup rules.
- Provide accessible alternatives for students who cannot install locally.

## 5. Teaching Materials to Develop for Every Module

| Instructor materials | Student materials | Assessment materials |
|---|---|---|
| Lesson plan with timing and misconceptions | Concept primer and glossary | Diagnostic pre-check |
| Annotated demonstration script | Guided workbook | Short knowledge check |
| Expected outputs and troubleshooting notes | Starter code and environment file | Independent transfer task |
| Answer key and facilitation prompts | Dataset card and provenance record | Interpretation rubric |
| Accessibility and adaptation notes | | Submission checklist |
| Reproducibility audit | | |

**Standard module lesson pattern:** Engage (present a biological decision
or flawed result) → Explain (introduce data structure, method,
assumptions, and limitations) → Demonstrate (model a transparent workflow
with commentary) → Practice (complete a guided analysis with checkpoints)
→ Transfer (solve a related problem independently) → Defend (justify
choices and interpret uncertainty) → Package (submit a reproducible
artifact and reflection).

## 6. Selected Reference Frameworks and Resources

These sources support the curriculum's competency-based and
reproducibility-centered design. Tool documentation should be checked and
versioned when individual teaching modules are authored.

- Network for Integrating Bioinformatics into Life Sciences Education (NIBLSE). Core Competencies. <https://qubeshub.org/community/groups/niblse/core_competencies>
- Wilson Sayres, M. A., et al. (2018). Bioinformatics core competencies for undergraduate life sciences education. *PLOS ONE*, 13(6), e0196878. <https://doi.org/10.1371/journal.pone.0196878>
- Bioconductor. Analyzing RNA-seq data with DESeq2. <https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html>
- Bioconductor. edgeR package documentation. <https://bioconductor.org/packages/release/bioc/html/edgeR.html>
- Galaxy Training Network. Training materials and learning pathways. <https://training.galaxyproject.org/>
- Galaxy Training Network. Metagenomics data processing and analysis for microbiome learning pathway. <https://training.galaxyproject.org/training-material/learning-pathways/metagenomics.html>

**Final recommendation:** retain the full program as a modular curriculum
map. For any single course offering, require Modules 0 and 1, select a
coherent set of core analyses, and use the remaining domains as electives
or capstone tracks. This preserves depth while still exposing students to
the broader bioinformatics landscape.
