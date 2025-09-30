rm(list=ls())
gc()

# Load Packages
library(Rcpp)
library(RcppArmadillo)
library(MCMCpack)
library(MASS)
library(tidyverse)
library(data.table)
library(ggplot2)
library(pROC)
library(mclust)
library(clue)
library(coda)
library(cluster)
library(glmnet)

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)
## How to run?
## Rscript script.R 100 20 3 5000 1000 2 2 ... "/path/to/output"

if (length(args) < 9) {
  stop("Usage: Rscript script.R <n> <p> <K> <n_iter> <burn_in> <tau2> <rho0> <outpath>")
}

gene <- args[1]
ct <- args[2]
K <- as.numeric(args[3])
n_iter <- as.numeric(args[4])
burn_in <- as.numeric(args[5])
tau2 <- as.numeric(args[6])
rho0 <- as.numeric(args[7])
p1 <- as.numeric(args[8])
add_b_mu <- as.logical(args[9])
outpath <- args[10]
batch <- args[11]
center <- as.logical(args[12])

sourceCpp("/ix1/wchen/Shiyue/Projects/2025_01_CE_Reg_P2G/02_P2G/01_DOGMA/Code/02_CLIPER/01_CER_PCGS.cpp")

CER_PCGS <- function(X, y, K, n_iter = 10000, burn_in = 5000, 
                     tau2 = 10000, rho0 = 1, p1 = 0.5, add_b_mu = FALSE, predict = FALSE, center = TRUE) {
  
  p <- ncol(X)
  n <- nrow(X)
  
  if(p <= K){
    
    K <- max(p - 1, 1)
  }
  
  if(center){
    
    X <- scale(X, center = T, scale = F)
    y <- scale(y, center = T, scale = F)
  }
  
  selected_idx <- seq_len(p)
  p_cap <- ceiling(1.5 * n)
  if (p > p_cap) {
    # Drop zero-variance columns first (cor() would yield NA)
    sds <- apply(X, 2, sd)
    nonzero <- which(sds > 0 & !is.na(sds))
    
    # Compute marginal correlations
    cors <- suppressWarnings(stats::cor(X[, nonzero, drop = FALSE], y))
    # If y has zero variance, cor gives NA everywhere
    if (all(is.na(cors))) stop("y has zero variance or correlations could not be computed.")
    
    # Rank by absolute correlation
    ord <- order(abs(cors), decreasing = TRUE, na.last = NA)
    keep_m <- min(p_cap, length(ord))
    keep_in_nonzero <- ord[seq_len(keep_m)]
    selected_idx <- nonzero[keep_in_nonzero]
    
    # Subset X and update p
    X <- X[, selected_idx, drop = FALSE]
    p <- ncol(X)
  }
  
  # K-Means Clustering on Ridge Regression Coefficients
  ridge_cv <- cv.glmnet(X, y, alpha = 0, lambda = exp(seq(-5, 5, length = 100)))
  ridge_lambda <- ridge_cv$lambda.min
  ridge_fit <- glmnet(X, y, alpha = 0, lambda = ridge_lambda)
  beta_hat <- as.vector(coef(ridge_fit))[-1]  # Remove intercept
  
  # Perform K-Means Clustering
  kmeans_result <- kmeans(beta_hat, centers = K, nstart = 10)
  
  # Identify the cluster whose center is closest to zero
  closest_cluster <- which.min(abs(kmeans_result$centers))
  
  # Create assignment matrix
  m <- matrix(0, nrow = p, ncol = K)
  for (j in 1:p) {
    m[j, kmeans_result$cluster[j]] <- 1
  }
  
  # Swap the closest-to-zero cluster with the first column
  if (closest_cluster != 1) {
    m[, c(1, closest_cluster)] <- m[, c(closest_cluster, 1)]
  }
  
  # Ensure no empty clusters
  for (k in 1:K) {
    if (sum(m[, k]) == 0) {
      # Find the closest non-empty cluster
      non_empty_clusters <- which(colSums(m) > 1)
      closest_non_empty <- non_empty_clusters[which.min(abs(kmeans_result$centers[non_empty_clusters] - kmeans_result$centers[k]))]
      
      # Transfer one point from the closest non-empty cluster
      j_to_transfer <- which(m[, closest_non_empty] == 1)[1]
      m[j_to_transfer, ] <- 0
      m[j_to_transfer, k] <- 1
    }
  }
  
  res <- CER_PCGS_rcpp(X = X, 
                       y = y, 
                       K = K, 
                       n_iter = n_iter, 
                       burn_in = burn_in, 
                       tau2 = tau2, 
                       rho0 = rho0, 
                       p1 = p1, 
                       add_b_mu = add_b_mu,
                       predict = predict, 
                       m_init = m)
  
  n_post <- dim(res$m_samples)[3L]
  clusters <- apply(res$m_samples, c(1, 3), which.max)  # p x n_post
  
  posterior_pair <- Reduce(`+`, lapply(seq_len(n_post), function(s) {
    outer(clusters[, s], clusters[, s], `==`)
  })) / n_post
  
  T_keep <- dim(simplify2array(res$m_samples))[3]
  b_keep <- tail(res$b_samples, T_keep)
  clusters <- apply(simplify2array(res$m_samples), c(1, 3), which.max)  # p x T_keep
  
  effects <- matrix(
    b_keep[cbind(rep(1:T_keep, each = p), as.vector(clusters))],
    nrow = p, ncol = T_keep)
  posterior_predictor_b <- rowMeans(effects)
  
  res$m_samples <- NULL
  gc()
  
  return(list(posterior_sigma2 = res$posterior_sigma2, posterior_m = res$posterior_m, 
              posterior_mu = res$posterior_mu, posterior_b = res$posterior_b, posterior_pair = posterior_pair,
              b_samples = res$b_samples, y_samples = res$y_samples, posterior_y = res$posterior_y, 
              posterior_predictor_b = posterior_predictor_b, input_y = y, peaks = colnames(X)))
}

DATA_PATH <- "/ix1/wchen/Shiyue/Projects/2025_01_CE_Reg_P2G/02_P2G/01_DOGMA/Result/01_CLIPER/01_Input/Act_IL1B_IL23/"

cliper_input <- readRDS(paste0(DATA_PATH, ct, "/CLIPER_input_", batch, ".rds"))
Input_data <- cliper_input[[gene]]

if(!is.null(Input_data) && ncol(Input_data) > 3){
  
  X <- as.matrix(subset(Input_data, select = -y))
  y <- Input_data$y
  
  p <- ncol(X)
  n <- nrow(X)
  
  print(paste0("gene: ", gene))
  print(paste("celltype:", ct))
  print(paste("n:", n))
  print(paste("p:", p))
  print(paste("K:", K))
  print(paste("n_iter:", n_iter))
  print(paste("burn_in:", burn_in))
  print(paste("tau2:", tau2))
  print(paste("rho0:", rho0))
  print(paste("p1:", p1))
  print(paste("add_b_mu:", add_b_mu))
  print(paste("outpath:", outpath))
  print(paste("batch:", batch))
  print(paste("center:", center))
  
  outfile <- file.path(outpath, paste0("mcmc_sim_", batch, "_", gene, "_", 
                                       K, "_", tau2, "_", rho0, "_", p1, ".rds"))
  if (!file.exists(outfile)) {
    
    set.seed(20250302)
    # Run MCMC
    start <- Sys.time()
    mcmc_results <- CER_PCGS(X, y, K, n_iter, burn_in, tau2, rho0, p1, add_b_mu = add_b_mu, predict = TRUE, center = center)
    end <- Sys.time()
    print(end-start)
    
    saveRDS(mcmc_results, outfile)
  }

} else{
  
  message(sprintf("Too few peaks for %s in %s", gene, ct))
}


