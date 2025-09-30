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
tcell_dogma$celltype_metacell <- gsub(" ", "_", tcell_dogma$celltype_metacell)

celltype <- fread(paste0(WORK_PATH, "celltype.txt"), header = F)

Idents(tcell_dogma) <- "condition"
for(stim in c("Act_IL1B_IL23", "Act_IL1B_IL23_PGE2", "Act_IL1B_IL23_TGFB", "Act_IL1B_IL23_PGE2_TGFB")){
  
  # subset data
  tcell_dogma_subset <- subset(x = tcell_dogma, idents = stim)
  
  for (ct in celltype$V1){
    
    subset_obj <- subset(tcell_dogma_subset, subset = celltype_metacell == ct)
    
    mc <- readRDS(paste0(WORK_PATH, stim, "/Pseudo_Rep/metacell_", ct, ".rds"))
    mc$barcode <- gsub("#", "_", mc$barcode)
    
    subset_obj <- subset_obj[, Seurat::Cells(subset_obj) %in% mc$barcode]
    
    # RNA aggregation
    rna_counts <- GetAssayData(subset_obj, assay = "RNA", slot = "counts")
    rna_counts <- rna_counts[,mc$barcode]
    meta_rna <- aggregate.Matrix(t(rna_counts), groupings = mc$replicate, fun = "sum")
    
    if(!file.exists(paste0(WORK_PATH, stim, "/Pseudo_Rep/RNA"))){
      
      system(paste0("mkdir ", WORK_PATH, stim, "/Pseudo_Rep/RNA"))
    }
    saveRDS(meta_rna, paste0(WORK_PATH, stim, "/Pseudo_Rep/RNA/tcell_dogma_metacell_RNA_", ct, ".rds"))
    
    print(paste0(stim, ct))
  }
}

