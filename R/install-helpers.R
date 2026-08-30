cliper_install_deps <- function(include_seurat_wrappers = TRUE, ask = TRUE) {
  cran <- c(
    "clue", "cluster", "coda", "cols4all", "data.table", "dplyr", "FNN", "ggplot2", "glmnet",
    "MASS", "Matrix", "MCMCpack", "mclust", "Rcpp",
    "RcppArmadillo", "patchwork", "remotes", "scales", "Seurat", "Signac", "stringr",
    "tibble", "tidyverse"
  )
  bioc <- c(
    "BiocGenerics", "Biostrings", "BSgenome.Hsapiens.UCSC.hg38", "DESeq2", "edgeR",
    "GenomeInfoDb", "GenomicRanges", "IRanges", "S4Vectors",
    "rtracklayer", "SummarizedExperiment"
  )

  missing_cran <- cran[!vapply(cran, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_cran) > 0) {
    install.packages(missing_cran)
  }

  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }

  missing_bioc <- bioc[!vapply(bioc, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_bioc) > 0) {
    BiocManager::install(missing_bioc, ask = ask)
  }

  if (isTRUE(include_seurat_wrappers) &&
      !requireNamespace("SeuratWrappers", quietly = TRUE)) {
    remotes::install_github("satijalab/seurat-wrappers")
  }

  invisible(TRUE)
}
