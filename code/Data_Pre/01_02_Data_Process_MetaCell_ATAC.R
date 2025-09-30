rm(list=ls())
gc()

# Load Packages
library(data.table)
library(tidyverse)
library(dplyr)
library(Seurat)
library(ArchR)
library(Signac)
library(Matrix.utils)

# # Parameters
# DATA_PATH <- "/ix1/wchen/zhongli/RWorkSpace/DOGMA-seq/PBMC/output/"
# WORK_PATH <- "/ix1/wchen/Shiyue/Projects/2025_01_CE_Reg_P2G/02_P2G/01_DOGMA/Data/"
# 
# setwd(WORK_PATH)
# 
# celltype <- fread(paste0(WORK_PATH, "celltype.txt"), header = F)
# peaks <- readRDS(paste0(DATA_PATH, '/pseudo_bulk/peak.RDS'))
# 
# # ArchR 500bp peaks
# for(stim in c("Act_IL1B_IL23", "Act_IL1B_IL23_PGE2", "Act_IL1B_IL23_TGFB", "Act_IL1B_IL23_PGE2_TGFB")){
#   
#   for(ct in celltype$V1){
#     
#     # Load barcode-to-metacell mapping table
#     barcode_to_metacell <- readRDS(paste0(WORK_PATH, stim, "/MetaCell/metacell_", ct, ".rds"))
#     
#     # Subset Peak data
#     peaks_subset <- peaks[,rownames(barcode_to_metacell)]
#     meta_atac <- aggregate.Matrix(t(peaks_subset), groupings = barcode_to_metacell$metacell_id, fun = "sum")
#     saveRDS(meta_atac, paste0(WORK_PATH, stim, "/MetaCell/Peak/tcell_dogma_metacell_ATAC_", ct, ".rds"))
#   }
# }

#####################################################################################################
# ArchR tile
addArchRThreads(threads = 8)
# Parameters
DATA_PATH <- "/ix1/rduerr/shared/rduerr_wchen/Shiyue/Chen_Files/2023_07_DOGMA_Revision/SCARlink/Result/"
WORK_PATH <- "/ix1/wchen/Shiyue/Projects/2025_01_CE_Reg_P2G/02_P2G/01_DOGMA/Data/"

setwd(WORK_PATH)

celltype <- fread(paste0(WORK_PATH, "celltype.txt"), header = F)
# Load ArchR project
scatac.object <- loadArchRProject(paste0(DATA_PATH, 'DOGMA_ATAC/'))

# Extract tiles
for(stim in c("Act_IL1B_IL23", "Act_IL1B_IL23_PGE2", "Act_IL1B_IL23_TGFB", "Act_IL1B_IL23_PGE2_TGFB")[-1]){
  
  barcodes_stim <- readRDS(paste0(WORK_PATH, stim, "/MetaCell/tcell_dogma_metacell_barcodes.rds"))
  
  # Format ArchR barcodes from "04191#CTCACAACAGCAAGTG-1" to "04191_CTCACAACAGCAAGTG-1"
  archr_barcodes <- rownames(scatac.object@cellColData)
  formatted_barcodes <- gsub("#", "_", archr_barcodes)  # Replace '#' with '_'
  names(formatted_barcodes) <- archr_barcodes  # Keep original as names for mapping
  
  # Filter barcodes present in both ArchR and RNA-derived metacells
  barcodes_stim <- barcodes_stim[barcodes_stim$barcode %in% formatted_barcodes, ]
  
  # Replace formatted barcode with original ArchR barcode (for pulling matrices)
  barcodes_stim$archr_barcode <- names(formatted_barcodes)[match(barcodes_stim$barcode, formatted_barcodes)]
  
  archr_cells <- barcodes_stim$archr_barcode
  
  setwd(paste0(WORK_PATH, stim, "/MetaCell/Tile/"))
  # Subset ArchR project to current cell type only
  proj_sub <- subsetArchRProject(
    ArchRProj = scatac.object,
    cells = archr_cells,
    dropCells = TRUE,
    logFile = createLogFile(paste0("subset_", stim)),
    force = TRUE)
  
  proj_sub <- addTileMatrix(proj_sub, binarize = FALSE, force = TRUE)
  
  # Load the TileMatrix from subset project
  tile_mat_ct <- getMatrixFromProject(
    proj_sub,
    useMatrix = "TileMatrix",
    binarize = FALSE)
  
  # Extract sparse matrix from SummarizedExperiment object
  mat_ct <- assays(tile_mat_ct)$TileMatrix
  rownames(mat_ct) <- paste0(tile_mat_ct@elementMetadata$seqnames, ":",
                             tile_mat_ct@elementMetadata$start, "-", tile_mat_ct@elementMetadata$start+500)
  
  saveRDS(mat_ct, file = file.path(WORK_PATH, paste0("tcell_dogma_filtered_tile_", stim, ".rds")))
  
  # Tile process
  for(ct in celltype$V1){
    
    # Load barcode-to-metacell mapping table
    barcode_to_metacell <- readRDS(paste0(WORK_PATH, stim, "/MetaCell/metacell_", ct, ".rds"))
    
    # Format ArchR barcodes from "04191#CTCACAACAGCAAGTG-1" to "04191_CTCACAACAGCAAGTG-1"
    archr_barcodes <- rownames(scatac.object@cellColData)
    formatted_barcodes <- gsub("#", "_", archr_barcodes)  # Replace '#' with '_'
    names(formatted_barcodes) <- archr_barcodes  # Keep original as names for mapping
    
    # Filter barcodes present in both ArchR and RNA-derived metacells
    barcode_to_metacell$barcode <- rownames(barcode_to_metacell)
    barcode_to_metacell <- barcode_to_metacell[barcode_to_metacell$barcode %in% formatted_barcodes, ]
    
    # Replace formatted barcode with original ArchR barcode (for pulling matrices)
    barcode_to_metacell$archr_barcode <- names(formatted_barcodes)[match(barcode_to_metacell$barcode, formatted_barcodes)]
    
    # Subset Tile data
    mat_ct_sub <- mat_ct[,barcode_to_metacell$archr_barcode]
    meta_tile <- aggregate.Matrix(t(mat_ct_sub), groupings = barcode_to_metacell$metacell_id, fun = "sum")
    saveRDS(meta_tile, paste0(WORK_PATH, stim, "/MetaCell/Tile/tcell_dogma_metacell_tile_", ct, ".rds"))
  }
  
  rm(proj_sub, tile_mat_ct, mat_ct)
  gc()
}








