rm(list=ls())
gc()

# Load Packages
library(data.table)
library(tidyverse)
library(dplyr)
library(Seurat)
library(ArchR)
library(Signac)
library(GenomicRanges)
library(stringr)

# Parameters
DATA_PATH <- "/ix1/wchen/Shiyue/Projects/2025_01_CE_Reg_P2G/02_P2G/01_DOGMA/Data/"
WORK_PATH <- "/ix1/wchen/Shiyue/Projects/2025_01_CE_Reg_P2G/02_P2G/01_DOGMA/Result/01_CLIPER/01_Input/"

###########################################################################################################################################
# Functions
# Step 1: Parse tile names into GRanges
parse_tiles <- function(tile_names) {
  
  tile_df <- str_match(tile_names, "^(chr[^:-]+)[:\\-](\\d+)-(\\d+)")
  
  # Remove rows with NA
  valid_idx <- complete.cases(tile_df)
  tile_df <- tile_df[valid_idx, ]
  tile_names_valid <- tile_names[valid_idx]
  
  start_pos <- as.numeric(tile_df[, 3])
  end_pos   <- as.numeric(tile_df[, 4])
  valid_range <- which(!is.na(start_pos) & !is.na(end_pos) & end_pos >= start_pos)
  tile_df <- tile_df[valid_range, ]
  tile_names_valid <- tile_names_valid[valid_range]
  start_pos <- start_pos[valid_range]
  end_pos   <- end_pos[valid_range]
  
  tile_gr <- GRanges(
    seqnames = tile_df[, 2],
    ranges = IRanges(start = start_pos, end = end_pos),
    tile_name = tile_names_valid)
  
  return(tile_gr)
}

# Step 2: Convert gene_ref into ±250kb GRanges
expand_gene_window <- function(gene_df, flank = 250000) {
  
  gene_gr <- GRanges(seqnames = paste0("chr", gene_df$chromosome_name), 
                     IRanges(start = pmax(1, gene_df$start_position - flank), 
                             end = gene_df$end_position + flank),
                     gene = gene_ref$hgnc_symbol)
  
  return(gene_gr)
}

# Step 3: Match tiles to expanded gene windows
get_tiles_near_genes <- function(tile_names, gene_ref, flank = 250000) {
  
  tile_gr <- parse_tiles(tile_names)
  gene_gr <- expand_gene_window(gene_ref, flank)
  hits <- findOverlaps(tile_gr, gene_gr)
  matched_tiles <- tile_gr[queryHits(hits)]$tile_name
  
  return(unique(matched_tiles))
}

###########################################################################################################################################
# Load data
celltype <- fread(paste0(DATA_PATH, "celltype.txt"), header = F)
gene <- fread(paste0(DATA_PATH, "genes_2000.txt"), header = F)
genes_hg38 <- readRDS("/ix1/wchen/Shiyue/Projects/2025_01_CE_Reg/Data/genes_hg38.rds")

for(stim in c("Act_IL1B_IL23")){
    
    for(ct in celltype$V1[9:10]){
      
      for(cell in c("meta", "rep")){
        
        #for(atac in c("peak", "tile")){
        for(atac in c("tile")){
          
          for(rna_norm in c("delta", "deseq")){
            
            for(atac_norm in c("delta", "deseq", "gcfq")){
              
              meta_rna <- readRDS(paste0(WORK_PATH, stim, "/", ct, "/", cell, "_rna_", rna_norm, ".rds"))
              meta_atac <- readRDS(paste0(WORK_PATH, stim, "/", ct, "/", cell, "_", atac, "_", atac_norm, ".rds"))
              
              cell_name <- intersect(rownames(meta_rna), colnames(meta_atac))
              meta_rna <- meta_rna[cell_name,]
              meta_atac <- meta_atac[,cell_name]
              meta_atac <- t(meta_atac)
              
              cliper_input <- list()
              
              for(gene_select in gene$V1){
                
                gene_ref <- genes_hg38 %>% filter(hgnc_symbol %in% c(gene_select))
                gene_ref$chromosome_name <- as.character(gene_ref$chromosome_name)
                
                if(gene_select %in% colnames(meta_rna) & nrow(gene_ref) > 0){
                  
                  peaks_near_genes <- get_tiles_near_genes(tile_names = colnames(meta_atac),
                                                           gene_ref = gene_ref,
                                                           flank = 250000)
                  
                  if (is.null(peaks_near_genes) || length(peaks_near_genes) == 0) next
                  meta_atac_select <- meta_atac[,peaks_near_genes]
                  
                  inputs <- cbind(y = meta_rna[, gene_select, drop = TRUE],
                                  as.data.frame(meta_atac_select, stringsAsFactors = FALSE))
                  cliper_input[[gene_select]] <- inputs
                }
              }
              
              saveRDS(cliper_input, paste0(WORK_PATH, stim, "/", ct, "/CLIPER_input_", cell, "_rna_", rna_norm, "_", atac, "_", atac_norm, ".rds"))
              print(paste0(ct, "_", cell, "_rna_", rna_norm, "_", atac, "_", atac_norm))
          }
        }
      }
    }
  }
}


###########################################################################################################################################
# Check Residual
# Load data
cliper_input <- readRDS("/ix1/wchen/Shiyue/Projects/2025_01_CE_Reg_P2G/02_P2G/01_DOGMA/Result/01_CLIPER/01_Input/Act_IL1B_IL23/Th17/CLIPER_input_meta_rna_deseq_peak_gcfq.rds")

residuals_list <- list()

for (g in names(cliper_input)) {
  
  df <- cliper_input[[g]]
  
  df <- scale(df, center = T, scale = F)
  
  fit <- lm(y ~ ., data = data.frame(df))
  res <- residuals(fit)
  residuals_list[[g]] <- res
}

res_mat <- do.call(cbind, residuals_list)
colnames(res_mat) <- names(residuals_list)
res_mat <- t(res_mat)

# QQ plot
vals <- as.numeric(res_mat)
qqnorm(vals, main = "Q-Q plot of all residuals")
qqline(vals, col = "red")

###########################################################################################################################################
# skewness / kurtosis
library(moments)

gene_skew <- apply(res_mat, 1, skewness)
gene_kurt <- apply(res_mat, 1, kurtosis)

par(mfrow=c(1,2))
hist(gene_skew, breaks=50, main="Skewness across genes", xlab="Skewness")
hist(gene_kurt, breaks=50, main="Kurtosis across genes", xlab="Kurtosis")
par(mfrow=c(1,1))

###########################################################################################################################################
# Density plot
library(ggplot2)
set.seed(1)

# random select 50 genes
sel_genes <- sample(rownames(res_mat), 36)
df <- data.frame(
  value = as.vector(res_mat[sel_genes, ]),
  gene = rep(sel_genes, each = ncol(res_mat))
)

ggplot(df, aes(x = value)) +
  geom_density(fill="grey", alpha=0.5) +
  facet_wrap(~gene, scales="free_y") +
  theme_bw() +
  labs(title="Density plots of selected genes (z-scored)",
       x="Normalized expression (z)", y="Density")

###########################################################################################################################################
normality_results <- list()

for (g in names(cliper_input)) {
  
  res <- res_mat[g,]
  
  if(sum(res) > 0){
    
    # Shapiro–Wilk
    shapiro_p <- NA
    if (length(res) >= 3 && length(res) <= 5000) {
      shapiro_p <- shapiro.test(res)$p.value
    }
    
    # KS test
    ks_p <- NA
    if (length(res) > 10) {
      ks_p <- ks.test(
        scale(res), "pnorm" 
      )$p.value
    }
    
    normality_results[[g]] <- data.frame(
      gene = g,
      shapiro_p = shapiro_p,
      ks_p = ks_p
    )
  }
}

normality_df <- do.call(rbind, normality_results)

nrow(normality_df)
sum(normality_df$ks_p > 0.1)
sum(normality_df$shapiro_p > 0.1)

# Summay Meta Cell Number
metacell <- data.frame()

for (ct in celltype$V1) {
  cliper_input <- readRDS(paste0(WORK_PATH, "Act_IL1B_IL23/", ct,
                                 "/CLIPER_input_meta_rna_deseq_peak_gcfq.rds"))
  
  ncell <- nrow(cliper_input[[1]])
  
  npeak_vec <- sapply(cliper_input, ncol)
  
  metacell <- rbind(
    metacell,
    data.frame(
      celltype = ct,
      ncell    = ncell,
      npeak    = npeak_vec
    )
  )
  
  print(ct)
}

metacell$p_n <- metacell$npeak / metacell$ncell

png(paste0(PLOT_PATH, "05_MetaCell_N_P.png"), width = 1800, height = 1800, res = 300)
ggplot(metacell, aes(x = celltype, y = npeak/ncell)) +
  geom_boxplot(fill = "lightblue") +
  labs(x = "Cell type", y = "npeak / ncell",
       title = "Distribution of npeak/ncell across cell types") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 55, hjust = 1))
dev.off()

ord <- metacell %>%
  dplyr::as_tibble() %>%
  dplyr::group_by(celltype) %>%
  dplyr::summarise(med = median(ratio, na.rm = TRUE), .groups = "drop") %>%
  dplyr::arrange(dplyr::desc(med)) %>%
  dplyr::pull(.data$celltype)

desired_levels <- c(
  "CD8+ Naive (Activated)",
  "CD8+ Naive (Resting)",
  "CD8+ Memory (Activated)",
  "CD8+ Memory (Resting)",
  "CD4+ Naive (Activated)",
  "CD4+ Naive (Resting)",
  "CD4+ Treg",
  "Th1",
  "Th17",
  "Other CD4+ Memory"
)
desired_levels <- gsub(" ", "_", desired_levels)

metacell <- metacell %>%
  mutate(celltype = gsub(" ", "_", celltype),
         celltype = factor(celltype, levels = desired_levels),
         ratio = npeak / ncell)

# UMAP
cell_order <- desired_levels
cell_colors <- c("#2874A6", "#FFD166", "#FFC8A2", "#D46A6A", 
                 "#40B5AD", "#B0DFE6", "#CF9FFF", "#FF9D5C", "#9FE2BF", "#F3B0C3")
names(cell_colors) <- cell_order

PLOT_PATH <- "/ix1/wchen/Shiyue/Projects/2025_01_CE_Reg_P2G/02_P2G/01_DOGMA/Result/01_CLIPER/Plots/Simulation/"

png(paste0(PLOT_PATH, "05_MetaCell_N_P.png"), width = 2400, height = 1400, res = 300)
ggplot(metacell, aes(x = celltype, y = ratio, fill = celltype)) +
  geom_violin(color = "grey30", alpha = 0.9, trim = TRUE) +
  geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white", color = "grey30") +
  stat_summary(fun = median, geom = "point", shape = 23, size = 2.5,
               fill = "white", color = "black") +
  scale_fill_manual(values = cell_colors, drop = FALSE) +
  labs(x = NULL, y = "npeak / ncell",
       title = "Ratio of npeak/ncell across cell types") +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 55, hjust = 1),
    panel.grid.major.x = element_blank(),
    plot.title = element_text(face = "bold"),
    legend.position = "none"
  )
dev.off()
