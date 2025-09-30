rm(list=ls())
gc()

# Load Packages
library(data.table)
library(tidyverse)
library(dplyr)
library(Seurat)
library(ArchR)
library(Signac)
library(BSgenome.Hsapiens.UCSC.hg38)
library(SeuratWrappers)
library(Matrix.utils)
library(Matrix)

# Parameters
DATA_PATH <- "/ix1/rduerr/shared/rduerr_wchen/Shiyue/Chen_Files/2023_07_DOGMA_Revision/SCARlink/Result/"
WORK_PATH <- "/ix1/wchen/Shiyue/Projects/2025_01_CE_Reg_P2G/02_P2G/01_DOGMA/Data/"

setwd(WORK_PATH)

celltype <- fread(paste0(WORK_PATH, "celltype.txt"), header = F)
# Load ArchR project
scatac.object <- loadArchRProject(paste0(DATA_PATH, 'DOGMA_ATAC/'))

# Combine cell types
cd <- ArchR::getCellColData(scatac.object, select = "celltype_updated")
cells <- rownames(cd)
ctu <- as.character(cd$celltype_updated)

ct_meta <- dplyr::case_when(
  ctu %in% c("CD4+ Regulatory (Resting)", "CD4+ Regulatory (Activated)") ~ "CD4+ Treg",
  ctu %in% c("CD4+ Memory (Resting) - Th1",  "CD4+ Memory (Activated) - Th1")  ~ "Th1",
  ctu %in% c("CD4+ Memory (Resting) - Th17", "CD4+ Memory (Activated) - Th17") ~ "Th17",
  ctu %in% c("CD4+ Memory (Resting) - Other","CD4+ Memory (Activated) - Other")~ "Other CD4+ Memory",
  ctu %in% c("MAITs (Resting)", "MAITs (Activated)") ~ "MAITs",
  TRUE ~ ctu
)

ct_meta <- gsub(" ", "_", ct_meta)
names(ct_meta) <- cells
scatac.object <- ArchR::addCellColData(scatac.object, data = ct_meta, 
                                       cells = names(ct_meta), name = "celltype_metacell")

#####################################################################################################
# ArchR pseudo-replicates
meta <- as.data.frame(scatac.object@cellColData)
meta$cell <- rownames(meta)

#for(stim in c("Act_IL1B_IL23", "Act_IL1B_IL23_PGE2", "Act_IL1B_IL23_TGFB", "Act_IL1B_IL23_PGE2_TGFB")){
for(stim in c("Act_IL1B_IL23_TGFB", "Act_IL1B_IL23_PGE2_TGFB")){
  
  # Filter cells
  cells <- meta$cell[ meta[["condition"]] == stim & meta[["celltype_metacell"]] %in% celltype$V1]
  if (length(cells) == 0L) {
    message(sprintf("[SKIP] %s / %s: 0 cells", stim, ct))
    next
  }
  
  outdir <- paste0("ArchR_subset_", stim)
  system(paste0("mkdir ", outdir))
  
  proj_sub <- subsetArchRProject(
    ArchRProj = scatac.object,
    cells = cells,
    outputDirectory = outdir,
    force = TRUE)
  
  # Tiles
  tile_mt <- readRDS(paste0(WORK_PATH, "tcell_dogma_filtered_tile_", stim, ".rds"))
  
  barcode_map_list <- list()
  for(ct in celltype$V1){
    
    # Subset
    outdir_ct <- paste0("temp/ArchR_subset_", stim, "_", ct)
    system(paste0("mkdir ", outdir_ct))
    
    cells <- meta$cell[ meta[["condition"]] == stim & meta[["celltype_metacell"]] == ct]
    
    proj_sub_ct <- subsetArchRProject(
      ArchRProj = proj_sub,
      cells = cells,
      outputDirectory = outdir_ct,
      force = TRUE)
    
    # Pseudo Rep
    n_cell <- length(cells)
    min_rep <- floor(length(cells)/20)
    pseu_rep <- addGroupCoverages(ArchRProj = proj_sub_ct, groupBy = "celltype_metacell", minReplicates = min_rep, useLabels = FALSE, 
                                  minCells = 20, maxCells = 100, maxReplicates = 1000, returnGroups = TRUE, force = TRUE)
    barcode_to_rep <- lapply(names(pseu_rep), function(ct) {
      lapply(names(pseu_rep[[ct]]), function(repname) {
        data.frame(
          barcode   = pseu_rep[[ct]][[repname]],
          replicate = repname,
          cell_type = ct,
          stringsAsFactors = FALSE
        )
      }) %>% bind_rows()
    }) %>% bind_rows()
    
    barcode_to_rep$replicate <- paste0(barcode_to_rep$cell_type, "_", barcode_to_rep$replicate)
    barcode_map_list[[ct]] <- barcode_to_rep
    saveRDS(barcode_to_rep, paste0(WORK_PATH, stim, "/Pseudo_Rep/metacell_", ct, ".rds"))
    
    # Get cells
    archr_barcode <- barcode_to_rep$barcode
    
    # Subset Tile data
    mat_ct_sub <- tile_mt[,archr_barcode]
    meta_tile <- aggregate.Matrix(t(mat_ct_sub), groupings = barcode_to_rep$replicate, fun = "sum")
    saveRDS(meta_tile, paste0(WORK_PATH, stim, "/Pseudo_Rep/Tile/tcell_dogma_metacell_tile_", ct, ".rds"))
    
    rm(proj_sub_ct); gc()
    try(unlink(outdir_ct, recursive = TRUE, force = TRUE), silent = TRUE)
    system(paste0("rm -r ", outdir_ct))
  }
  
  barcode_to_metacell <- do.call(rbind, barcode_map_list)
  saveRDS(barcode_to_metacell, paste0(WORK_PATH, stim, "/Pseudo_Rep/tcell_dogma_metacell_barcodes.rds"))
  
  rm(proj_sub); gc()
  try(unlink(outdir, recursive = TRUE, force = TRUE), silent = TRUE)
}

#####################################################################################################
# ArchR 500bp Peaks
peaks <- readRDS('/ix1/wchen/zhongli/RWorkSpace/DOGMA-seq/PBMC/output/pseudo_bulk/peak.RDS')
for(stim in c("Act_IL1B_IL23", "Act_IL1B_IL23_PGE2", "Act_IL1B_IL23_TGFB", "Act_IL1B_IL23_PGE2_TGFB")){

  for(ct in celltype$V1){

    # Load barcode-to-metacell mapping table
    barcode_to_metacell <- readRDS(paste0(WORK_PATH, stim, "/Pseudo_Rep/metacell_", ct, ".rds"))
    barcode_to_metacell$barcode <- gsub("#", "_", barcode_to_metacell$barcode)

    # Subset Peak data
    peaks_subset <- peaks[,barcode_to_metacell$barcode]
    meta_atac <- aggregate.Matrix(t(peaks_subset), groupings = barcode_to_metacell$replicate, fun = "sum")
    saveRDS(meta_atac, paste0(WORK_PATH, stim, "/Pseudo_Rep/Peak/tcell_dogma_metacell_ATAC_", ct, ".rds"))
    
    print(paste0(stim, ct))
  }
}


