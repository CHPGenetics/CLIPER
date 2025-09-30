rm(list=ls())
gc()

# Load Packages
library(tidyverse)
library(ggplot2)
library(mclust)
library(clue)
library(coda)
library(cluster)
library(glmnet)
library(RColorBrewer)
library(MASS)
library(MCMCpack)
library(igraph)
library(data.table)

DATA_PATH <- "/ix1/wchen/Shiyue/Projects/2025_01_CE_Reg_P2G/02_P2G/01_DOGMA/Result/01_CLIPER/01_Input/"
WORK_PATH <- "/ix1/wchen/Shiyue/Projects/2025_01_CE_Reg_P2G/02_P2G/01_DOGMA/Result/01_CLIPER/02_Output/wo_b_mu0/"

## Based on pair matrix
find_clusters <- function(matrix) {
  p <- nrow(matrix)
  visited <- rep(FALSE, p)
  clusters <- rep(NA, p)
  cluster_id <- 0
  
  for (i in 1:p) {
    if (!visited[i]) {
      cluster_id <- cluster_id + 1
      queue <- c(i)
      while (length(queue) > 0) {
        node <- queue[1]
        queue <- queue[-1]
        if (!visited[node]) {
          visited[node] <- TRUE
          clusters[node] <- cluster_id
          neighbors <- which(matrix[node, ] == 1 & !visited)
          queue <- c(queue, neighbors)
        }
      }
    }
  }
  return(clusters)
}

# Clustering Accuracy
calculate_accuracy <- function(true_matrix, pred_matrix) {
  
  if (!all(dim(true_matrix) == dim(pred_matrix))) {
    stop("The matrices must have the same dimensions.")
  }
  
  # if (!all(true_matrix == t(true_matrix)) || !all(pred_matrix == t(pred_matrix))) {
  #   stop("Both matrices should be symmetric for pairwise comparison.")
  # }
  
  upper_tri_indices <- upper.tri(true_matrix)
  true_values <- true_matrix[upper_tri_indices]
  pred_values <- pred_matrix[upper_tri_indices]
  
  TP <- sum(true_values == 1 & pred_values == 1)
  TN <- sum(true_values == 0 & pred_values == 0)
  FP <- sum(true_values == 0 & pred_values == 1)
  FN <- sum(true_values == 1 & pred_values == 0)
  
  accuracy <- (TP + TN) / (TP + TN + FP + FN)
  
  cat("True Positives (TP):", TP, "\n")
  cat("True Negatives (TN):", TN, "\n")
  cat("False Positives (FP):", FP, "\n")
  cat("False Negatives (FN):", FN, "\n")
  cat("Accuracy:", accuracy*100, "%\n")
  
  return(accuracy)
}

celltype <- fread("/ix1/wchen/Shiyue/Projects/2025_01_CE_Reg_P2G/02_P2G/01_DOGMA/Data/celltype.txt", header = F)
genes <- fread("/ix1/wchen/Shiyue/Projects/2025_01_CE_Reg_P2G/02_P2G/01_DOGMA/Data/genes_2000.txt", header = F)
p1 <- 0.5
K = 5

batch_list <- c()
# for(cell in c("meta", "rep")){
for(cell in c("meta")){
  
  for(atac in c("peak")){

    for(rna_norm in c("delta", "deseq")){
      
      for(atac_norm in c("delta", "deseq", "gcfq")){
        
        batch_list <- c(batch_list, paste0(cell, "_rna_", rna_norm, "_", atac, "_", atac_norm))
      }
    }
  }
}

for(stim in c("Act_IL1B_IL23")){
  
  for(batch in batch_list){
    
    for(center in c("Center", "No_Center")){
      
      Results <- data.frame()
      non_cluster1_all <- data.frame()
      
      for(ct in celltype$V1){
        
        print(ct)
        for(gene in genes$V1){
          
          PATH <- Sys.glob(paste0(WORK_PATH, center, "/", ct, "/mcmc_sim_", batch, "_", gene, "*.rds"))
          
          if(length(PATH) == 1){
            
            result <- readRDS(PATH)

            y <- result$input_y
            
            posterior_m <- result$posterior_m
            posterior_y <- result$posterior_y
            posterior_b <- result$posterior_b
            predicted_cluster <- apply(posterior_m, 1, which.max)
            
            # R square
            R2 <- 1 - sum((y - posterior_y)^2) / sum((y - mean(y))^2)
            
            # Cluster table summary
            cluster_counts <- paste0(table(predicted_cluster), collapse = ",")
            
            # Predictors not in cluster 1
            idx_non_cluster1 <- which(predicted_cluster != 1)
            if (length(idx_non_cluster1) > 0) {
              tmp_df <- data.frame(
                Peak = result$peaks[idx_non_cluster1],
                Gene = gene,
                Cell_Type = ct,
                Cluster = predicted_cluster[idx_non_cluster1],
                Posterior_b = posterior_b[predicted_cluster[idx_non_cluster1]],
                stringsAsFactors = FALSE)
              non_cluster1_all <- rbind(non_cluster1_all, tmp_df)
            }
            
            non_cluster1_names <- if (length(idx_non_cluster1) > 0) {
              paste(result$peaks[idx_non_cluster1], collapse = ";")
            } else {
              NA
            }
            
            # Summary info
            info <- c(gene = gene, ct = ct, R2 = R2, 
                      K_predicted = length(unique(predicted_cluster)), 
                      cluster = paste0(table(predicted_cluster), collapse = ","))
            Results <- rbind(Results, info)
          }
        }
        
      }
      
      colnames(Results) <- c("Gene", "Cell_Type", "R2", "K_predicted", "Cluster_n")
      saveRDS(Results, paste0(WORK_PATH, center, "/Results_summary_", batch, ".rds"))
      saveRDS(non_cluster1_all, paste0(WORK_PATH, center, "/non_cluster1_all", batch, ".rds"))
    }
  }
}

Results$Other_Clusters_Total <- sapply(Results$Cluster_n, function(x) {
  nums <- as.numeric(unlist(strsplit(x, ",")))
  if (length(nums) <= 1) {
    return(0)
  } else {
    return(sum(nums[-1]))
  }
})


for(ct in celltype$V1){
  
  Result_sub <- Results %>% filter(Cell_Type == ct)
  non_cluster1_all_sub <- non_cluster1_all %>% filter(Cell_Type == ct)
  
  Result_sub <- Result_sub %>%
    rowwise() %>%
    mutate(
      assoc_values = list(as.numeric(unlist(strsplit(Cluster_n, ",")))[-1]),
      assoc_min = ifelse(length(assoc_values) > 0, min(assoc_values), 0),
      assoc_max = ifelse(length(assoc_values) > 0, max(assoc_values), 0),
      assoc_sum = ifelse(length(assoc_values) > 0, sum(assoc_values), 0)
    ) %>%
    ungroup() %>%
    dplyr::select(-assoc_values)
  
  cat("Cell Type:", ct, "\n")
  cat("  Gene Number:", length(unique(Result_sub$Gene)), "\n")
  cat("  Genes with Significant Associations:", 
      sum(Result_sub$assoc_sum != 0), 
      "out of", length(unique(Result_sub$Gene)), "\n")
  cat("  Association Number (Mean (Range)):", 
      round(mean(Result_sub$assoc_sum)), 
      "(", min(Result_sub$assoc_sum), "-", max(Result_sub$assoc_sum), ")\n")
  cat("  Positive Associations:", sum(non_cluster1_all_sub$Posterior_b > 0), "\n")
  cat("  Negative Associations:", sum(non_cluster1_all_sub$Posterior_b < 0), "\n")
  cat("\n")
}


non_cluster1_all <- readRDS(paste0("/ix1/wchen/Shiyue/Projects/2025_01_CE_Reg_P2G/02_P2G/01_DOGMA/Result/01_CLIPER/02_Output/wo_b_mu0/Center", "/non_cluster1_allmeta_rna_deseq_peak_gcfq.rds"))

mean_beta <- mean(non_cluster1_all$Posterior_b, na.rm = TRUE)
ggplot(non_cluster1_all, aes(x = Posterior_b)) +
  geom_histogram(
    bins = 30,  
    color = "black",
    fill = "#2874A6",
    alpha = 0.7) +
  geom_vline(aes(xintercept = mean_beta), 
             color = "red", linetype = "dashed", size = 1) +
  annotate("text",
           x = mean_beta, 
           y = 40000, 
           label = paste0("Mean = ", round(mean_beta, 2)),
           hjust = -0.1, color = "red") +
  labs(
    title = "Distribution of Beta",
    x = "Beta",
    y = "Count"
  ) +
  theme_classic()

means_df <- non_cluster1_all %>%
  group_by(Cell_Type) %>%
  summarise(mean_beta = mean(Posterior_b, na.rm = TRUE))

pastel_colors <- c("#2874A6", "#FFD166", "#FFC8A2", "#D46A6A", "#40B5AD",
                   "#B0DFE6", "#CF9FFF", "#FF9D5C", "#77DD77", "#F3B0C3")
ggplot(non_cluster1_all, aes(x = Posterior_b)) +
  geom_histogram(aes(fill = Cell_Type),
                 bins = 30,
                 color = "white",
                 alpha = 0.9) +
  geom_vline(data = means_df,
             aes(xintercept = mean_beta, color = "black"),
             linetype = "dashed", size = 0.8, show.legend = FALSE) +
  geom_text(data = means_df,
            aes(x = mean_beta, y = Inf,
                label = paste0("Mean=", round(mean_beta, 2)),
                color = "black"),
            vjust = 1.2, size = 3, show.legend = FALSE) +
  facet_wrap(~ Cell_Type, scales = "free_y") +
  scale_fill_manual(values = pastel_colors) +
  scale_color_manual(values = pastel_colors) +
  labs(title = "Beta distribution by cell type",
       x = "Beta", y = "Count") +
  theme_classic() +
  theme(legend.position = "none")

df_pie <- non_cluster1_all %>%
  dplyr::mutate(Sign = ifelse(Posterior_b > 0, "Positive", "Negative")) %>%
  dplyr::count(Cell_Type, Sign, name = "n") %>%         
  dplyr::group_by(Cell_Type) %>%
  dplyr::mutate(prop  = n / sum(n),
                label = paste0(round(prop * 100, 1), "%")) %>%
  dplyr::ungroup()

ggplot(df_pie, aes(x = "", y = prop, fill = Sign)) +
  geom_col(width = 1, color = "white") +
  coord_polar(theta = "y") +
  facet_wrap(~ Cell_Type) +
  geom_text(aes(label = label), 
            position = position_stack(vjust = 0.5), 
            color = "black", size = 3) +
  scale_fill_manual(values = c("Positive" = "#77DD77", "Negative" = "#FF9D5C")) +
  labs(title = "Proportion of positive/negative Posterior_b by cell type",
       fill = "Sign") +
  theme_void() +
  theme(strip.text = element_text(size = 10, face = "bold"))




