rm(list=ls())
gc()

# Load Packages
library(data.table)
library(tidyverse)
library(dplyr)
library(Seurat)
library(ArchR)
library(Signac)
library(DESeq2)

# Parameters
DATA_PATH <- "/ix1/wchen/Shiyue/Projects/2025_01_CE_Reg_P2G/02_P2G/01_DOGMA/Data/"
WORK_PATH <- "/ix1/wchen/Shiyue/Projects/2025_01_CE_Reg_P2G/02_P2G/01_DOGMA/Result/01_CLIPER/01_Input/"

###########################################################################################################################################
# Functions
# Delta method
# Each col is cell
normalize_delta <- function(counts) {

  # counts: gene x cell matrix
  lib_size <- colSums(counts)
  L <- mean(lib_size)

  size_factors <- lib_size / L
  norm_counts <- sweep(counts, 2, size_factors, "/")
  log_norm_counts <- log1p(norm_counts)  # log(x+1)
  
  # filter low expression genes
  n_cells <- ncol(log_norm_counts)
  keep <- rowSums(counts > 0) >= ceiling(0.1 * n_cells) # filter gene express in > 10% cells
  log_norm_counts <- log_norm_counts[keep, , drop = FALSE]
  
#  log_norm_counts <- t(scale(t(log_norm_counts))) # z score
  
  return(log_norm_counts)
}


###########################################################################################################################################
# DESeq2 Normlize
## RNA
normalize_metacell_deseq2 <- function(
    counts                  # gene x metacell
){
  stopifnot(is.matrix(counts) || inherits(counts, "dgCMatrix"))

  keep_cols <- colSums(counts) > 0
  counts <- counts[, keep_cols, drop = FALSE]

  meta_info <- data.frame(dummy = factor(rep("all", ncol(counts))))
  rownames(meta_info) <- colnames(counts)
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = counts, colData = meta_info, design = ~ 1)

  dds <- DESeq2::estimateSizeFactors(dds)
  mat <- SummarizedExperiment::assay(DESeq2::vst(dds, blind = TRUE))
  
  # filter low expression genes
  n_cells <- ncol(counts)
  keep_genes <- rowSums(counts > 0) >= ceiling(0.1 * n_cells) # filter gene express in > 10% cells
  mat <- mat[keep_genes, , drop = FALSE]
  
#  mat <- t(scale(t(mat))) # z score
  
  return(mat)
}

###########################################################################################################################################
# Meta Cell
gene <- fread(paste0(DATA_PATH, "genes_2000.txt"), header = F)
celltype <- fread(paste0(DATA_PATH, "celltype.txt"), header = F)

for(stim in c("Act_IL1B_IL23")){
  
  if(!file.exists(paste0(WORK_PATH, stim))){
    
    system(paste0("mkdir ", paste0(WORK_PATH, stim)))
  }
  
  for(ct in celltype$V1){
    
    dir_path <- file.path(WORK_PATH, stim, ct)
    if (!dir.exists(dir_path)) {
      dir.create(dir_path, recursive = TRUE)
    }
    
    # Load data 
    meta_rna <- readRDS(paste0(DATA_PATH, stim, "/MetaCell/RNA/tcell_dogma_metacell_RNA_", ct, ".rds"))
    rep_rna <- readRDS(paste0(DATA_PATH, stim, "/Pseudo_Rep/RNA/tcell_dogma_metacell_RNA_", ct, ".rds"))
    
    # Normalization
    ## Delta
    meta_rna_delta <- normalize_delta(t(meta_rna))
    rep_rna_delta <- normalize_delta(t(rep_rna))
    
    ## DESeq2
    meta_rna_deseq <- normalize_metacell_deseq2(t(meta_rna))
    rep_rna_deseq <- normalize_metacell_deseq2(t(rep_rna))
    
    # filter genes
    meta_rna_delta <- meta_rna_delta[intersect(gene$V1, rownames(meta_rna_delta)),]
    rep_rna_delta <- rep_rna_delta[intersect(gene$V1, rownames(rep_rna_delta)),]
    
    meta_rna_deseq <- meta_rna_deseq[intersect(gene$V1, rownames(meta_rna_deseq)),]
    rep_rna_deseq <- rep_rna_deseq[intersect(gene$V1, rownames(rep_rna_deseq)),]
    
    ## Transform
    saveRDS(t(meta_rna_delta), paste0(WORK_PATH, stim, "/", ct, "/meta_rna_delta.rds"))
    saveRDS(t(meta_rna_deseq), paste0(WORK_PATH, stim, "/", ct, "/meta_rna_deseq.rds"))
    
    saveRDS(t(rep_rna_delta), paste0(WORK_PATH, stim, "/", ct, "/rep_rna_delta.rds"))
    saveRDS(t(rep_rna_deseq), paste0(WORK_PATH, stim, "/", ct, "/rep_rna_deseq.rds"))
  }
}


###########################################################################################################################################
norm <- rep_rna_delta

# QQ plot
vals <- as.numeric(norm)
qqnorm(vals, main = "Q-Q plot of all normalized values")
qqline(vals, col = "red")

###########################################################################################################################################
# skewness / kurtosis
library(moments)

gene_skew <- apply(norm, 1, skewness)
gene_kurt <- apply(norm, 1, kurtosis)

par(mfrow=c(1,2))
hist(gene_skew, breaks=50, main="Skewness across genes", xlab="Skewness")
hist(gene_kurt, breaks=50, main="Kurtosis across genes", xlab="Kurtosis")
par(mfrow=c(1,1))

###########################################################################################################################################
# Density plot
library(ggplot2)
set.seed(1)

# random select 50 genes
sel_genes <- sample(rownames(norm), 50)
df <- data.frame(
  value = as.vector(norm[sel_genes, ]),
  gene = rep(sel_genes, each = ncol(norm))
)

ggplot(df, aes(x = value)) +
  geom_density(fill="grey", alpha=0.5) +
  facet_wrap(~gene, scales="free_y") +
  theme_bw() +
  labs(title="Density plots of selected genes (z-scored)",
       x="Normalized expression (z)", y="Density")
