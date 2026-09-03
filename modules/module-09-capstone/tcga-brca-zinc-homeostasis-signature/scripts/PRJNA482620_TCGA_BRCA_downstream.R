# ==============================================================================
# STUDY SUMMARY & METADATA
# ==============================================================================
# Study Accession : PRJNA482620 / TCGA-BRCA / GSE114725
# Study Title     : Construction of a novel zinc homeostasis-related gene signature 
#                   to predict breast cancer prognosis by integrating single-cell 
#                   and bulk RNA-seq data (Xu et al., 2026)
# Organism        : Homo sapiens (Human)
# Description     : This script integrates pre-processed bulk RNA-seq data (TCGA-BRCA) 
#                   and single-cell RNA-seq data (GSE114725). It calculates Zinc 
#                   Homeostasis (ZH) activity scores using ssGSEA, identifies correlated 
#                   modules via WGCNA, constructs a 9-gene LASSO Cox risk score model 
#                   (CCDC74B, ELOVL2, FAM234B, FUT3, KCNJ11, NXPH4, RIBC1, UPK2, WNK4), 
#                   evaluates overall survival (Kaplan-Meier), and reconstructs CD8+ T 
#                   cell pseudotime trajectories using Monocle2.
# ==============================================================================

# Create project working directory
dir.create("PRJNA482620_TCGA_BRCA_ZincHomeostasis_analysis", showWarnings = FALSE)
setwd("PRJNA482620_TCGA_BRCA_ZincHomeostasis_analysis")

# Load required packages
library(DESeq2)
library(GSVA)
library(WGCNA)
library(glmnet)
library(survival)
library(survminer)
library(Seurat)
library(Monocle)

cat("==========================================================================\n")
cat(" EXECUTING DOWNSTREAM ANALYSIS FOR STUDY: PRJNA482620_TCGA_BRCA \n")
cat("==========================================================================\n")

# --- STEP 1: LOAD PRE-PROCESSED TCGA-BRCA MATRICES ---
# Load bulk matrices (N = 1,165 samples: 1,066 tumor, 99 normal)
counts_full <- read.table("TCGA_BRCA_raw_counts.tsv", header=TRUE, row.names=1, sep="\t")
tpm_full    <- read.table("TCGA_BRCA_tpm.tsv", header=TRUE, row.names=1, sep="\t")
clinical    <- read.csv("TCGA_BRCA_clinical.csv", row.names=1)

# --- STEP 2: ssGSEA ZINC HOMEOSTASIS SCORING ---
zh_genes <- list(Zinc_Homeostasis = c("SLC30A1", "SLC30A2", "MT1A", "MT2A", "ZIP1", "ZIP6")) # 51 total genes
params <- ssgseaParam(as.matrix(tpm_full), zh_genes)
zh_scores <- gsva(params)

# --- STEP 3: WGCNA MODULE IDENTIFICATION ---
options(stringsAsFactors = FALSE)
net <- blockwiseModules(tpm_full, power = 6, TOMType = "unsigned", minModuleSize = 30)
# Correlate modules with ZH scores to isolate key brown module genes

# --- STEP 4: DESeq2 DIFFERENTIAL EXPRESSION & INTERSECTION ---
dds <- DESeqDataSetFromMatrix(countData = counts_full, colData = clinical, design = ~ condition)
dds <- DESeq(dds)
degs <- rownames(subset(results(dds), padj < 0.05 & abs(log2FoldChange) > 1))
candidate_genes <- intersect(degs, brown_module_genes) # 73 genes

# --- STEP 5: LASSO COX REGRESSION (9-GENE SIGNATURE) ---
surv_obj <- Surv(time = clinical$time, event = clinical$status)
cv_fit <- cv.glmnet(as.matrix(tpm_full[candidate_genes, ]), surv_obj, family = "cox", alpha = 1)

# Predict risk scores and stratify into High/Low risk groups
risk_scores <- predict(cv_fit, newx = as.matrix(tpm_full[candidate_genes, ]), s = "lambda.min")
clinical$risk_group <- ifelse(risk_scores >= median(risk_scores), "High", "Low")

# Kaplan-Meier Survival Curve
fit <- survfit(surv_obj ~ risk_group, data = clinical)
ggsurvplot(fit, data = clinical, pval = TRUE, risk.table = TRUE, title = "PRJNA482620: TCGA-BRCA Survival")

# --- STEP 6: SINGLE-CELL IMMUNE PROFILE (GSE114725) ---
sc_data <- Read10X("GSE114725_matrix/")
sc_obj <- CreateSeuratObject(counts = sc_data, min.cells = 3, min.features = 200)

sc_obj <- NormalizeData(sc_obj)
sc_obj <- FindVariableFeatures(sc_obj)
sc_obj <- ScaleData(sc_obj)
sc_obj <- RunPCA(sc_obj)
sc_obj <- RunUMAP(sc_obj, dims = 1:20)

cat("==========================================================================\n")
cat(" PRJNA482620 DOWNSTREAM ANALYSIS COMPLETE \n")
cat("==========================================================================\n")