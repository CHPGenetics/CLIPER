rm(list=ls())
gc()

# Load Packages
library(data.table)
library(tidyverse)
library(dplyr)
library(Seurat)
library(ArchR)
library(Signac)
library(SeuratWrappers)
library(Matrix.utils)
library(FNN)
library(Matrix)
library(Rcpp)
library(RcppArmadillo)

# Parameters
DATA_PATH <- "/ix1/rduerr/shared/rduerr_wchen/Shiyue/Chen_Files/2023_07_DOGMA_Revision/SCARlink/Result/"
WORK_PATH <- "/ix1/wchen/Shiyue/Projects/2025_01_CE_Reg_P2G/02_P2G/01_DOGMA/Data/"

setwd(WORK_PATH)

# Load DOGMA-seq data
tcell_dogma <- readRDS("/ix1/wchen/zhongli/RWorkSpace/DOGMA-seq/PBMC/output/tcell_annotated_updated.RDS")
hvg_genes <- VariableFeatures(tcell_dogma, nfeatures = 2000)
write.table(hvg_genes, paste0(WORK_PATH, "genes_2000.txt"), row.names = F, col.names = F, sep = '\t', quote = F)

tcell_dogma$celltype_updated <- as.character(tcell_dogma$celltype_updated)
tcell_dogma$celltype_metacell <- ifelse(tcell_dogma$celltype_updated %in%
                                          c("CD4+ Regulatory (Resting)", "CD4+ Regulatory (Activated)"),
                                        "CD4+ Treg", ifelse(tcell_dogma$celltype_updated %in%
                                                                    c("CD4+ Memory (Resting) - Th1", "CD4+ Memory (Activated) - Th1"),
                                                                  "Th1", ifelse(tcell_dogma$celltype_updated %in%
                                                                                                c("CD4+ Memory (Resting) - Th17", "CD4+ Memory (Activated) - Th17"),
                                                                                              "Th17", ifelse(tcell_dogma$celltype_updated %in%
                                                                                                                             c("CD4+ Memory (Resting) - Other", "CD4+ Memory (Activated) - Other"),
                                                                                                                           "Other CD4+ Memory", ifelse(tcell_dogma$celltype_updated %in%
                                                                                                                                                           c("MAITs (Resting)", "MAITs (Activated)"),
                                                                                                                                                         "MAITs", tcell_dogma$celltype_updated)))))

table(tcell_dogma$celltype_metacell)
Idents(tcell_dogma) <- "celltype_metacell"
tcell_dogma <- subset(x = tcell_dogma, idents = c("CD4+ Memory (Activated) - Tfh",
                                                  "CD4+ Memory (Resting) - Tfh", 
                                                  "CD8+ Regulatory", "Gamma Delta"), invert = TRUE)
tcell_dogma <- subset(x = tcell_dogma, idents = c("MAITs"), invert = TRUE)
table(tcell_dogma$celltype_metacell)

DefaultAssay(tcell_dogma) <- "RNA"
tcell_dogma[["SCT"]] <- NULL

# save cell type info 
tcell_dogma$celltype_metacell <- gsub(" ", "_", tcell_dogma$celltype_metacell)
write.table(unique(tcell_dogma$celltype_metacell), paste0(WORK_PATH, "celltype.txt"),
            row.names = F, col.names = F, sep = '\t', quote = F)

##################################################################################################### 
# Meta cell
sourceCpp("/ix1/wchen/Shiyue/Projects/2025_01_CE_Reg_P2G/02_P2G/01_DOGMA/Code/01_Data_process/KNN_Utils.cpp")
## Function
make_metacell_groups <- function(
    obj,
    reduction = "lsi",
    dims_to_use = 2:30,
    cells_to_use = NULL,
    ct_name = "MC",
    # --- defaults tuned for ~20 cells per group ---
    k = 20,                          # target neighborhood size per seed
    knn_iteration = 500,             # number of seed neighborhoods to propose
    overlap_cutoff = 0.8,            # absolute overlap threshold (= floor(0.8*k))
    seed = 1,
    # post-processing controls
    enforce_unique = TRUE,           # assign each cell to exactly one metacell
    ensure_cover  = TRUE,            # assign any uncovered cells to nearest metacell
    balance_sizes = TRUE,            # split very large groups / merge very small groups
    min_size = NULL,                 # default: ceiling(0.6*k) ~= 12
    max_size = NULL                  # default: 2*k ~= 40
){
  set.seed(seed)
  
  ## Embedding subset
  emb <- Seurat::Embeddings(obj, reduction = reduction)
  emb <- emb[, dims_to_use, drop = FALSE]
  if (!is.null(cells_to_use)) {
    all_cells  <- Seurat::Cells(obj)
    keep_cells <- cells_to_use[cells_to_use %in% all_cells]  # preserves cells_to_use order
    emb <- emb[keep_cells, , drop = FALSE]
  }
  n <- nrow(emb)
  if (n == 0) stop("No cells left after filtering.")
  if (n < 2) stop("Need at least 2 cells to form neighborhoods.")
  
  ## kNN neighborhoods
  k_eff <- min(k, n)
  seeds <- sample(seq_len(n), size = knn_iteration, replace = n < knn_iteration)
  
  nn_index <- FNN::get.knnx(
    data  = emb,
    query = emb[seeds, , drop = FALSE],
    k     = k_eff
  )$nn.index
  # sort each row and keep matrix shape
  nn_index <- t(apply(nn_index, 1, sort))
  storage.mode(nn_index) <- "integer"
  
  ## Overlap pruning (0=keep, -1=drop)
  thr <- max(1L, floor(overlap_cutoff * k_eff))
  keep_flag <- determineOverlapCpp(nn_index, thr)
  nn_index  <- nn_index[keep_flag == 0L, , drop = FALSE]
  
  cell_names <- rownames(emb)
  groups <- lapply(seq_len(nrow(nn_index)), function(i) cell_names[nn_index[i, ]])
  
  ## Enforce unique assignment: greedy take-most-unassigned-first
  if (enforce_unique && length(groups) > 0) {
    # sort once by size (rough proxy for coverage gain)
    ord <- order(vapply(groups, length, integer(1)), decreasing = TRUE)
    groups <- groups[ord]
    
    assigned <- setNames(rep(FALSE, n), cell_names)
    uniq_groups <- vector("list", length(groups))
    kept <- logical(length(groups))
    
    for (g in seq_along(groups)) {
      mem <- groups[[g]]
      new_mem <- mem[!assigned[mem]]
      if (length(new_mem) > 0) {
        uniq_groups[[g]] <- new_mem
        assigned[new_mem] <- TRUE
        kept[g] <- TRUE
      } else {
        uniq_groups[[g]] <- character(0)
      }
    }
    groups <- uniq_groups[kept]
  }
  
  ## Ensure coverage: assign any uncovered cells to nearest group centroid
  if (ensure_cover && length(groups) > 0) {
    covered <- unique(unlist(groups, use.names = FALSE))
    if (length(covered) < n) {
      centers <- t(vapply(groups, function(g) colMeans(emb[g, , drop = FALSE]),
                          numeric(ncol(emb))))
      missing <- setdiff(cell_names, covered)
      if (length(missing) > 0) {
        nn_to_group <- FNN::get.knnx(
          data  = centers,
          query = emb[missing, , drop = FALSE],
          k     = 1
        )$nn.index[,1]
        for (i in seq_along(missing)) {
          groups[[nn_to_group[i]]] <- c(groups[[nn_to_group[i]]], missing[i])
        }
      }
    }
  }
  
  ## Balance sizes: split very large groups; merge very small groups
  if (balance_sizes && length(groups) > 0) {
    tgt_min <- if (is.null(min_size)) ceiling(0.6 * k_eff) else as.integer(min_size)
    tgt_max <- if (is.null(max_size)) 2L * k_eff           else as.integer(max_size)
    
    for (iter in 1:5) {  # a few passes usually suffice
      # split too-large groups
      new_groups <- list()
      for (g in seq_along(groups)) {
        mem <- groups[[g]]; sz <- length(mem)
        if (sz > tgt_max) {
          nsplit <- ceiling(sz / tgt_max)
          km <- stats::kmeans(emb[mem, , drop = FALSE], centers = nsplit, iter.max = 50, nstart = 2)
          for (c in seq_len(nsplit)) {
            sub <- mem[km$cluster == c]
            if (length(sub)) new_groups[[length(new_groups) + 1L]] <- sub
          }
        } else {
          new_groups[[length(new_groups) + 1L]] <- mem
        }
      }
      groups <- new_groups
      
      # recompute centers each pass
      centers <- t(vapply(groups, function(g) colMeans(emb[g, , drop = FALSE]), numeric(ncol(emb))))
      sizes   <- vapply(groups, length, integer(1))
      big_ids   <- which(sizes >= tgt_min)
      small_ids <- which(sizes <  tgt_min)
      
      if (length(small_ids) == 0L) break
      if (length(big_ids) == 0L)   break
      
      # merge each small into nearest big
      for (sid in small_ids) {
        mem <- groups[[sid]]; if (!length(mem)) next
        nn <- FNN::get.knnx(centers[big_ids, , drop = FALSE],
                            matrix(centers[sid, , drop = FALSE], nrow = 1), k = 1)$nn.index[1,1]
        tid <- big_ids[nn]
        groups[[tid]] <- c(groups[[tid]], mem)
        groups[[sid]] <- character(0)
      }
      groups <- Filter(length, groups)
    }
  }
  
  ## Build mapping (no NAs; if any appear, assign to nearest centroid)
  cell2metacell <- setNames(rep(NA_character_, n), cell_names)
  if (length(groups)) {
    for (g in seq_along(groups)) {
      cell2metacell[groups[[g]]] <- paste0(ct_name, "_", g)
    }
  }
  if (anyNA(cell2metacell) && length(groups) > 0) {
    na_cells <- names(cell2metacell)[is.na(cell2metacell)]
    centers <- t(vapply(groups, function(g) colMeans(emb[g, , drop = FALSE]),
                        numeric(ncol(emb))))
    nn <- FNN::get.knnx(
      data  = centers,
      query = emb[na_cells, , drop = FALSE],
      k     = 1
    )$nn.index[,1]
    for (i in seq_along(na_cells)) {
      cell2metacell[na_cells[i]] <- paste0(ct_name, "_", nn[i])
      groups[[nn[i]]] <- c(groups[[nn[i]]], na_cells[i])
    }
  }
  
  ## Diagnostics
  cat("dim(nn_index) =", paste(dim(nn_index), collapse=" x "), "\n")
  lens <- vapply(groups, length, integer(1))
  if (length(lens)) print(summary(lens))
  cat("#groups:", length(groups),
      " n_assigned:", sum(!is.na(cell2metacell)),
      " n_cells:", n, "\n")
  
  list(
    metacell_members = groups,
    cell2metacell    = cell2metacell,
    seeds_used       = rownames(emb)[seeds],
    k                = k_eff,
    knn_iteration    = knn_iteration,
    overlap_threshold= thr)
}

Idents(tcell_dogma) <- "condition"
for(stim in c("Act_IL1B_IL23", "Act_IL1B_IL23_PGE2", "Act_IL1B_IL23_TGFB", "Act_IL1B_IL23_PGE2_TGFB")){
  
  barcode_map_list <- list()
  # subset data
  tcell_dogma_subset <- subset(x = tcell_dogma, idents = stim)
  
  for (ct in unique(tcell_dogma_subset$celltype_metacell)){
    
    subset_obj <- subset(tcell_dogma_subset, subset = celltype_metacell == ct)
    
    subset_obj[["lsi"]] <- NULL
    subset_obj <- RunTFIDF(subset_obj)
    subset_obj <- FindTopFeatures(subset_obj, min.cutoff = "q0")
    subset_obj <- RunSVD(subset_obj)
    
    mc <- make_metacell_groups(obj = subset_obj, ct_name = ct)
    subset_obj$metacell_id <- mc$cell2metacell[Seurat::Cells(subset_obj)]
    metacell_id <- subset_obj$metacell_id
    
    # Save metacell info
    meta <- subset_obj@meta.data[, c("gex_barcode", "atac_barcode", "metacell_id")]
    saveRDS(meta, paste0(WORK_PATH, stim, "/MetaCell/metacell_", ct, ".rds"))
    
    # RNA aggregation
    rna_counts <- GetAssayData(subset_obj, assay = "RNA", slot = "counts")
    meta_rna <- aggregate.Matrix(t(rna_counts), groupings = metacell_id, fun = "sum")
    
    if(!file.exists(paste0(WORK_PATH, stim, "/MetaCell/RNA"))){
      
      system(paste0("mkdir ", WORK_PATH, stim, "/MetaCell/RNA"))
    }
    saveRDS(meta_rna, paste0(WORK_PATH, stim, "/MetaCell/RNA/tcell_dogma_metacell_RNA_", ct, ".rds"))
    
    barcode_map_list[[ct]] <- data.frame(barcode = colnames(subset_obj),
                                         metacell = metacell_id, celltype = ct)
    
    print(ct)
  }
  
  barcode_to_metacell <- do.call(rbind, barcode_map_list)
  saveRDS(barcode_to_metacell, paste0(WORK_PATH, stim, "/MetaCell/tcell_dogma_metacell_barcodes.rds"))
}




