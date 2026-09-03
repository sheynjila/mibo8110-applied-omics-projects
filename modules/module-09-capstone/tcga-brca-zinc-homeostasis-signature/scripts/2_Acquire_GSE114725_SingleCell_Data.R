###############################################################################
# CAPSTONE · SCRIPT 2 of 3 — ACQUIRE GSE114725 SINGLE-CELL SUPPLEMENTARY DATA
# ==============================================================================
# NEW SCRIPT, not part of the delivered material. PRJNA482620_TCGA_BRCA_downstream.R's
# Step 6 calls Read10X("GSE114725_matrix/") without saying how that directory
# gets populated. This script downloads GSE114725's own GEO-deposited
# supplementary files directly - the standard, defensible capstone move of
# reusing a published study's own processed data rather than reprocessing raw
# reads a second time - instead of trying to regenerate a count matrix from
# raw FASTQs (see the disclosed finding below for why that specific path is
# not available to this project anyway).
#
# A DISCLOSED FINDING THIS SCRIPT'S OWN OUTPUT MUST BE CHECKED AGAINST, NOT
# ASSUMED AWAY:
#   GSE114725 (Azizi et al. 2018, Cell - "Single-cell map of diverse immune
#   phenotypes in the breast tumor microenvironment") is the SAME dataset
#   Module 7 (Single-Cell Elective) investigated and found is sequenced with
#   the InDrop protocol, not 10x Genomics (real SRA study SRP148597 /
#   BioProject PRJNA472383 - verified via NCBI eutils; see Module 7's README
#   for the full finding, including two OTHER draft scripts that misattributed
#   this dataset to a different, wrong BioProject/accession entirely).
#   PRJNA482620_TCGA_BRCA_downstream.R's Read10X("GSE114725_matrix/") call
#   assumes the standard 10x trio (barcodes.tsv/features.tsv/matrix.mtx) -
#   this dataset's own GEO supplementary files are NOT confirmed to be in that
#   format (this script's own NCBI eutils check found only a generic "CSV"
#   supplementary-file-type flag on the series, not a per-file listing). This
#   is left as an open, disclosed gap: after running this script, INSPECT
#   whatever actually downloads into geo_supp/ before assuming Read10X() in
#   Script 3 will work unmodified - a genes-by-cells CSV/TXT matrix (a common
#   alternative deposit format for non-10x single-cell studies) needs a
#   different loader (read.csv() + CreateSeuratObject() directly), not
#   Read10X(). Do not "fix" this by silently swapping the loader without
#   first opening the downloaded file and confirming its actual shape.
#
# Install once if needed (Bioconductor, not CRAN):
#   if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
#   BiocManager::install("GEOquery")
###############################################################################

library(GEOquery)

GEO_ACCESSION <- "GSE114725"
OUTDIR <- "geo_supp"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

cat("==================================================================\n")
cat(sprintf(" Downloading %s supplementary files from GEO\n", GEO_ACCESSION))
cat("==================================================================\n")

downloaded <- getGEOSuppFiles(GEO_ACCESSION, baseDir = OUTDIR, makeDirectory = TRUE)
print(downloaded)

cat(sprintf("\nDownloaded %d file(s) to %s/%s/\n", nrow(downloaded), OUTDIR, GEO_ACCESSION))
cat("\n[ACTION REQUIRED] Open the file(s) listed above and confirm their actual format\n")
cat("before running Script 3:\n")
cat("  - If they are the standard 10x trio (barcodes.tsv[.gz], features.tsv[.gz] or genes.tsv[.gz],\n")
cat("    matrix.mtx[.gz]) in one directory: Read10X(\"geo_supp/GSE114725/\") as\n")
cat("    PRJNA482620_TCGA_BRCA_downstream.R's Step 6 already assumes will work unmodified.\n")
cat("  - If they are a plain genes-by-cells CSV/TSV/TXT matrix (the more likely case for an\n")
cat("    InDrop-derived dataset - see this script's header): Script 3's Step 6 needs to be\n")
cat("    changed from Read10X(...) to something like:\n")
cat("        mat <- read.csv(\"geo_supp/GSE114725/<actual filename>\", row.names = 1)\n")
cat("        sc_data <- CreateSeuratObject(counts = as(as.matrix(mat), \"dgCMatrix\"), ...)\n")
cat("    - a genuine, disclosed correction to Tier A that belongs in this project's own README\n")
cat("    changelog once made, not a silent edit.\n")
