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
# Function
# DESeq2 Normlize
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
  
  # filter low expression peaks
  n_cells <- ncol(counts)
  keep_genes <- rowSums(counts > 0) >= ceiling(0.05 * n_cells) # filter peaks in > 5% cells
  mat <- mat[keep_genes, , drop = FALSE]
  
#  mat <- t(scale(t(mat))) # z score
  
  return(mat)
}
  

###########################################################################################################################################
# Delta method
# Each col is cell
normalize_delta <- function(counts) {
  
  # counts: gene x cell matrix
  lib_size <- colSums(counts)
  L <- mean(lib_size)
  
  size_factors <- lib_size / L
  norm_counts <- sweep(counts, 2, size_factors, "/")
  log_norm_counts <- log1p(norm_counts)  # log(x+1)
  
  # filter low expression peaks
  n_cells <- ncol(log_norm_counts)
  keep <- rowSums(counts > 0) >= ceiling(0.05 * n_cells)  # filter peaks in > 5% cells
  log_norm_counts <- log_norm_counts[keep, , drop = FALSE]
  
#  log_norm_counts <- t(scale(t(log_norm_counts))) # z score
  
  return(log_norm_counts)
}

###########################################################################################################################################
# smooth GC-FQ
qsmooth_nogroup <- function(object) {
  
  X <- as.matrix(object)
  rn <- rownames(X); cn <- colnames(X)
  
  X_sorted <- apply(X, 2, sort, partial = NULL, na.last = NA)
  Qref <- rowMeans(X_sorted)
  
  X_norm <- matrix(NA_real_, nrow(X), ncol(X))
  for (i in seq_len(ncol(X))) {
    ref <- Qref
    x   <- X[, i]
    rmin <- rank(x, ties.method = "min")
    dups <- duplicated(rmin)
    if (any(dups)) {
      rrand <- rank(x, ties.method = "random")
      tied.ranks <- unique(rmin[dups])
      for (k in tied.ranks) {
        sel <- rrand[rmin == k]
        ref[sel] <- ave(ref[sel])
      }
    }
    X_norm[, i] <- ref[rmin]
  }
  
  rownames(X_norm) <- rn
  colnames(X_norm) <- cn
  X_norm
}

gc_qsmooth_normalize <- function(counts) {
  stopifnot(is.matrix(counts) || inherits(counts, "dgCMatrix"))
  
  # filter low expression peaks
  n_cells <- ncol(counts)
  keep <- rowSums(counts > 0) >= ceiling(0.05 * n_cells)  # filter peaks in > 5% cells
  counts <- counts[keep, , drop = FALSE]
  counts <- as.matrix(counts)

  keep_autosomes <- function(gr) {
    GenomeInfoDb::seqlevelsStyle(gr) <- "UCSC"
    gr[as.character(GenomicRanges::seqnames(gr)) %in% paste0("chr", 1:22)]
  }
  make_peaks_gr <- function(peaks, bed0 = FALSE) {
    spl1 <- strsplit(peaks, "[:\\-]")
    chr  <- vapply(spl1, `[`, "", 1)
    st   <- as.integer(vapply(spl1, `[`, "", 2))
    en   <- as.integer(vapply(spl1, `[`, "", 3))
    if (isTRUE(bed0)) st <- st + 1L
    GenomicRanges::GRanges(seqnames = chr, ranges = IRanges::IRanges(start = st, end = en))
  }
  make_gc_bins <- function(gcContent, nbins = 20) {
    qs <- stats::quantile(gcContent, probs = seq(0, 1, length.out = nbins + 1), na.rm = TRUE)
    base::cut(gcContent, breaks = qs, include.lowest = TRUE, labels = FALSE)
  }
  gcqn_qsmooth <- function(counts_mat, gcGroups) {
    stopifnot(nrow(counts_mat) == length(gcGroups))
    gcGroups <- as.factor(gcGroups)
    
    out <- matrix(NA_real_, nrow(counts_mat), ncol(counts_mat),
                  dimnames = dimnames(counts_mat))
    
    for (lvl in levels(gcGroups)) {
      idx <- which(gcGroups == lvl)
      X   <- counts_mat[idx, , drop = FALSE]
      
      if (nrow(X) < 2L) {
        out[idx, ] <- as.matrix(X)        
      } else {
        X_qn <- qsmooth_nogroup(as.matrix(X))
        out[idx, ] <- X_qn
      }
    }
    out
  }
  
  peaks_gr_all <- make_peaks_gr(rownames(counts), bed0 = FALSE)
  peaks_gr     <- keep_autosomes(peaks_gr_all)
  
  keep_idx <- which(as.character(GenomicRanges::seqnames(peaks_gr_all)) %in% paste0("chr", 1:22))
  counts   <- counts[keep_idx, , drop = FALSE]
  stopifnot(nrow(counts) == length(peaks_gr))
  
  genome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
  seqs    <- Biostrings::getSeq(genome, peaks_gr)
  gc_frac <- rowSums(Biostrings::letterFrequency(seqs, c("G","C"), as.prob = TRUE))
  stopifnot(length(gc_frac) == nrow(counts))
  
  gc_bins  <- make_gc_bins(gc_frac, nbins = 20)
  counts_smooth <- gcqn_qsmooth(counts, gc_bins)
  
  log_smooth <- log1p(counts_smooth)
#  Z_smooth   <- t(scale(t(log_smooth)))
#  Z_smooth[is.na(Z_smooth)] <- 0 
  
  return(log_smooth)
}

###########################################################################################################################################
# # Peak
# celltype <- fread(paste0(DATA_PATH, "celltype.txt"), header = F)
# for(stim in c("Act_IL1B_IL23")){
#   
#   for(ct in celltype$V1){
#     
#     # Load data 
#     meta_atac <- readRDS(paste0(DATA_PATH, stim, "/MetaCell/Peak/tcell_dogma_metacell_ATAC_", ct, ".rds"))
#     rep_atac <- readRDS(paste0(DATA_PATH, stim, "/Pseudo_Rep/Peak/tcell_dogma_metacell_ATAC_", ct, ".rds"))
#     
#     # Normalization
#     ## Delta
#     meta_atac_delta <- normalize_delta(t(meta_atac))
#     rep_atac_delta <- normalize_delta(t(rep_atac))
#       
#     ## DESeq2
#     meta_atac_deseq <- normalize_metacell_deseq2(t(meta_atac))
#     rep_atac_deseq <- normalize_metacell_deseq2(t(rep_atac))
#     
#     ## smooth GC-FQ
#     meta_atac_gcfq <- gc_qsmooth_normalize(t(meta_atac))
#     rep_atac_gcfq <- gc_qsmooth_normalize(t(rep_atac))
#     
#     # Save
#     saveRDS(meta_atac_delta, paste0(WORK_PATH, stim, "/", ct, "/meta_peak_delta.rds"))
#     saveRDS(rep_atac_delta, paste0(WORK_PATH, stim, "/", ct, "/rep_peak_delta.rds"))
#     saveRDS(meta_atac_deseq, paste0(WORK_PATH, stim, "/", ct, "/meta_peak_deseq.rds"))
#     saveRDS(rep_atac_deseq, paste0(WORK_PATH, stim, "/", ct, "/rep_peak_deseq.rds"))
#     saveRDS(meta_atac_gcfq, paste0(WORK_PATH, stim, "/", ct, "/meta_peak_gcfq.rds"))
#     saveRDS(rep_atac_gcfq, paste0(WORK_PATH, stim, "/", ct, "/rep_peak_gcfq.rds"))
#   }
# }

###########################################################################################################################################
# Tile
celltype <- fread(paste0(DATA_PATH, "celltype.txt"), header = F)
for(stim in c("Act_IL1B_IL23")){
  
  for(ct in celltype$V1){
    
    # Load data 
    meta_atac <- readRDS(paste0(DATA_PATH, stim, "/MetaCell/Tile/tcell_dogma_metacell_tile_", ct, ".rds"))
    rep_atac <- readRDS(paste0(DATA_PATH, stim, "/Pseudo_Rep/Tile/tcell_dogma_metacell_tile_", ct, ".rds"))
    
    # Normalization
    ## Delta
    meta_atac_delta <- normalize_delta(t(meta_atac))
    rep_atac_delta <- normalize_delta(t(rep_atac))
    
    ## DESeq2
    meta_atac_deseq <- normalize_metacell_deseq2(t(meta_atac))
    rep_atac_deseq <- normalize_metacell_deseq2(t(rep_atac))
    
    ## smooth GC-FQ
    meta_atac_gcfq <- gc_qsmooth_normalize(t(meta_atac))
    rep_atac_gcfq <- gc_qsmooth_normalize(t(rep_atac))
    
    # Save
    saveRDS(meta_atac_delta, paste0(WORK_PATH, stim, "/", ct, "/meta_tile_delta.rds"))
    saveRDS(rep_atac_delta, paste0(WORK_PATH, stim, "/", ct, "/rep_tile_delta.rds"))
    saveRDS(meta_atac_deseq, paste0(WORK_PATH, stim, "/", ct, "/meta_tile_deseq.rds"))
    saveRDS(rep_atac_deseq, paste0(WORK_PATH, stim, "/", ct, "/rep_tile_deseq.rds"))
    saveRDS(meta_atac_gcfq, paste0(WORK_PATH, stim, "/", ct, "/meta_tile_gcfq.rds"))
    saveRDS(rep_atac_gcfq, paste0(WORK_PATH, stim, "/", ct, "/rep_tile_gcfq.rds"))
  }
}



