#!/usr/bin/env bash
# create_module_structure.sh
#
# Master command to scaffold the curriculum folder structure: docs/, and one
# folder per curriculum module (0-9) plus the supplemental case-studies
# folder. Creates directories only — no README.md files are written, so you
# can add each module's README yourself afterward.
#
# Default target is the CURRENT directory — i.e. it assumes you have already
# `cd`ed into your existing repo (e.g. mibo8110-applied-omics-projects) and
# just want the module folders to exist directly inside it, not nested one
# level deeper inside a brand-new subfolder.
#
# Safe to re-run: mkdir -p is idempotent, so running this again after you've
# added files won't overwrite or delete anything.
#
# Usage:
#   cd mibo8110-applied-omics-projects       # your existing repo root (has .git in it)
#   bash create_module_structure.sh          # scaffolds directly into the current directory
#
#   bash create_module_structure.sh /path/to/other-repo   # or target a different existing repo directly

set -euo pipefail

ROOT="${1:-.}"

echo "Scaffolding curriculum structure under: $ROOT (resolved: $(cd "$(dirname "$ROOT")" 2>/dev/null && pwd)/$(basename "$ROOT"))"

# --- cross-cutting docs ---
mkdir -p "$ROOT"
mkdir -p "$ROOT/docs/curriculum"

# --- Module 0 — Computing, Data, Statistics, and Reproducibility ---
mkdir -p "$ROOT/modules/module-00-foundations-reproducibility/scripts"
mkdir -p "$ROOT/modules/module-00-foundations-reproducibility/lecture"

# --- Module 1 — Raw Sequencing Reads and QC ---
mkdir -p "$ROOT/modules/module-01-raw-reads-qc/scripts"

# --- Module 2 — Bulk RNA-seq and Gene Expression ---
mkdir -p "$ROOT/modules/module-02-bulk-rnaseq/scripts"
mkdir -p "$ROOT/modules/module-02-bulk-rnaseq/lecture"

# --- Module 3 — Microbial Genome Assembly ---
mkdir -p "$ROOT/modules/module-03-genome-assembly/scripts"

# --- Module 4 — Variant Calling ---
mkdir -p "$ROOT/modules/module-04-variant-calling/scripts"

# --- Module 5 — Phylogenomics and Outbreak Context ---
mkdir -p "$ROOT/modules/module-05-phylogenomics/scripts"

# --- Module 6 — Metagenomics and Microbiome Profiles ---
mkdir -p "$ROOT/modules/module-06-metagenomics/scripts"

# --- Module 7 (Elective) — Single-Cell Transcriptomics ---
mkdir -p "$ROOT/modules/module-07-single-cell-elective/scripts"

# --- Module 8 (Elective) — Chromatin and Regulatory Genomics ---
mkdir -p "$ROOT/modules/module-08-chromatin-regulatory-elective/scripts"

# --- Module 9 — Integrated Capstone and Professional Portfolio ---
# No scripts/ subfolder here: each capstone project gets its own slug'd
# subfolder created on demand (see modules/module-09-capstone/README.md).
mkdir -p "$ROOT/modules/module-09-capstone"

# --- Supplemental applied case studies (AMR, M. tuberculosis) ---
mkdir -p "$ROOT/modules/supplemental-case-studies/scripts"

echo
echo "Done. Resulting structure:"
find "$ROOT" -type d | sort

echo
echo "Next: add each module's README.md at modules/module-NN-*/README.md"
echo "      (and docs/curriculum/<curriculum-file>.md for the canonical reference)."
