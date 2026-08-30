library(data.table)
library(tidyverse)
library(dplyr)
library(Seurat)
library(ArchR)
library(Signac)
library(SeuratWrappers)
library(FNN)
library(Matrix)
library(Rcpp)
library(RcppArmadillo)
library(GenomicRanges)
library(IRanges)
library(stringr)
library(MCMCpack)
library(MASS)
library(mclust)
library(clue)
library(coda)
library(cluster)
library(glmnet)
library(cols4all)
library(ggplot2)
library(patchwork)
library(rtracklayer)

.cliper_source_dir <- local({
  frames <- sys.frames()
  ofiles <- vapply(
    frames,
    function(x) {
      if (!is.null(x$ofile)) x$ofile else NA_character_
    },
    character(1)
  )
  ofiles <- ofiles[!is.na(ofiles)]
  if (length(ofiles) > 0) {
    dirname(normalizePath(ofiles[[length(ofiles)]], mustWork = FALSE))
  } else {
    getwd()
  }
})

sourceCpp(file.path(.cliper_source_dir, "KNN.cpp"))

make_metacell_groups <- function(
    obj,
    reduction = "lsi",
    dims_to_use = 2:30,
    cells_to_use = NULL,
    celltype_col = NULL,
    celltype_value = NULL,
    prefix = NULL,
    k = 10,
    knn_iteration = 500,
    overlap_cutoff = 0.3,
    seed = 1,
    balance_sizes = TRUE,
    min_size = NULL,
    max_size = NULL,
    target_metacells = NULL,
    max_metacells = 500,
    min_metacells = 200,
    target_cells_per_metacell = 15,
    min_cells_per_metacell = 10,
    candidate_multiplier = 5,
    adaptive_k = TRUE,
    balance_final_metacells = FALSE,
    force_target_metacells = FALSE,
    target_memberships_per_cell = 1.2,
    selection_method = c("coverage", "diverse"),
    allow_overlap = TRUE,
    min_coverage_warning = 0.75,
    max_mean_reuse_warning = 3,
    verbose = TRUE
) {
  set.seed(seed)
  
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Package 'Seurat' is required.")
  }
  if (!requireNamespace("FNN", quietly = TRUE)) {
    stop("Package 'FNN' is required.")
  }
  
  if (!exists("determineOverlapCpp", mode = "function")) {
    stop("Function 'determineOverlapCpp' is not found. Please source/compile KNN.cpp first.")
  }
  
  selection_method <- match.arg(selection_method)
  if (!is.logical(balance_final_metacells) || length(balance_final_metacells) != 1L) {
    stop("balance_final_metacells must be TRUE or FALSE.")
  }
  if (!is.logical(force_target_metacells) || length(force_target_metacells) != 1L) {
    stop("force_target_metacells must be TRUE or FALSE.")
  }
  if (!is.finite(target_memberships_per_cell) || target_memberships_per_cell <= 0) {
    stop("target_memberships_per_cell must be positive.")
  }
  
  meta <- obj@meta.data
  all_cells <- Seurat::Cells(obj)
  
  if (is.null(cells_to_use)) {
    cells_to_use <- all_cells
  } else {
    cells_to_use <- intersect(cells_to_use, all_cells)
  }
  
  if (length(cells_to_use) == 0L) {
    stop("No cells left after filtering cells_to_use.")
  }
  
  if (!reduction %in% names(obj@reductions)) {
    stop("Reduction '", reduction, "' is not present in the Seurat object.")
  }
  
  emb_all <- Seurat::Embeddings(obj, reduction = reduction)
  cells_to_use <- intersect(cells_to_use, rownames(emb_all))
  
  if (length(cells_to_use) == 0L) {
    stop("No selected cells are present in the requested embedding.")
  }
  
  dims_to_use <- dims_to_use[dims_to_use >= 1L & dims_to_use <= ncol(emb_all)]
  if (length(dims_to_use) == 0L) {
    stop("No valid dimensions in dims_to_use.")
  }
  
  emb_all <- emb_all[, dims_to_use, drop = FALSE]
  
  if (!is.null(max_metacells) && !is.finite(max_metacells)) {
    stop("max_metacells must be finite or NULL.")
  }
  if (!is.null(target_metacells) && !is.finite(target_metacells)) {
    stop("target_metacells must be finite or NULL.")
  }
  if (!is.null(min_metacells) && !is.finite(min_metacells)) {
    stop("min_metacells must be finite or NULL.")
  }
  if (!is.null(target_metacells)) target_metacells <- as.integer(target_metacells)
  if (!is.null(max_metacells)) max_metacells <- as.integer(max_metacells)
  if (!is.null(min_metacells)) min_metacells <- as.integer(min_metacells)
  if (!is.finite(target_cells_per_metacell) || target_cells_per_metacell <= 0) {
    stop("target_cells_per_metacell must be positive.")
  }
  min_cells_per_metacell <- as.integer(min_cells_per_metacell)
  if (!is.finite(min_cells_per_metacell) || min_cells_per_metacell < 2L) {
    stop("min_cells_per_metacell must be an integer >= 2.")
  }
  if (!is.finite(overlap_cutoff) || overlap_cutoff < 0 || overlap_cutoff > 1) {
    stop("overlap_cutoff must be between 0 and 1.")
  }
  if (!is.finite(candidate_multiplier) || candidate_multiplier < 1) {
    stop("candidate_multiplier must be >= 1.")
  }
  
  if (is.null(celltype_col)) {
    if (is.null(celltype_value)) {
      if (!is.null(prefix)) {
        ct_vec <- rep(prefix, length(cells_to_use))
      } else {
        ct_vec <- rep("metacell", length(cells_to_use))
      }
    } else {
      ct_vec <- rep(as.character(celltype_value), length(cells_to_use))
    }
    names(ct_vec) <- cells_to_use
  } else {
    if (!celltype_col %in% colnames(meta)) {
      stop("celltype_col must be a column name in obj@meta.data.")
    }
    
    ct_vec <- meta[cells_to_use, celltype_col, drop = TRUE]
    ct_vec <- as.character(ct_vec)
    names(ct_vec) <- cells_to_use
    
    keep_non_na <- !is.na(ct_vec)
    if (any(!keep_non_na)) {
      warning("Some cells have NA cell type labels and will be removed.")
      cells_to_use <- cells_to_use[keep_non_na]
      ct_vec <- ct_vec[keep_non_na]
    }
    
    if (!is.null(celltype_value)) {
      keep_value <- ct_vec == as.character(celltype_value)
      cells_to_use <- cells_to_use[keep_value]
      ct_vec <- ct_vec[keep_value]
    }
    
    if (length(cells_to_use) == 0L) {
      stop("No cells left after applying celltype_col/celltype_value filtering.")
    }
  }
  
  make_safe_name <- function(x) {
    x <- as.character(x)
    x <- gsub("[^A-Za-z0-9_]+", "_", x)
    x <- gsub("_+", "_", x)
    x <- gsub("^_|_$", "", x)
    ifelse(nchar(x) == 0L, "metacell", x)
  }
  
  sqdist_to_point <- function(mat, point) {
    rowSums((mat - matrix(point, nrow(mat), length(point), byrow = TRUE))^2)
  }
  
  select_diverse_groups <- function(groups, emb_work, max_keep) {
    if (is.null(max_keep) || length(groups) <= max_keep) {
      return(groups)
    }
    
    centers <- t(vapply(
      groups,
      function(g) colMeans(emb_work[g, , drop = FALSE]),
      numeric(ncol(emb_work))
    ))
    
    selected <- integer(max_keep)
    selected[1L] <- sample(seq_along(groups), 1L)
    dmin <- sqdist_to_point(centers, centers[selected[1L], ])
    
    if (max_keep >= 2L) {
      for (ii in 2:max_keep) {
        remaining <- setdiff(seq_along(groups), selected[seq_len(ii - 1L)])
        if (length(remaining) == 0L) {
          selected <- selected[seq_len(ii - 1L)]
          break
        }
        selected[ii] <- remaining[which.max(dmin[remaining])]
        dnew <- sqdist_to_point(centers, centers[selected[ii], ])
        dmin <- pmin(dmin, dnew)
      }
    }
    
    groups[selected]
  }
  
  select_coverage_groups <- function(groups, emb_work, max_keep) {
    if (is.null(max_keep) || length(groups) <= max_keep) {
      return(groups)
    }
    
    cell_names_local <- rownames(emb_work)
    groups_idx <- lapply(groups, function(g) match(g, cell_names_local))
    groups_idx <- lapply(groups_idx, function(g) g[!is.na(g)])
    keep_nonempty <- vapply(groups_idx, length, integer(1)) > 0L
    groups <- groups[keep_nonempty]
    groups_idx <- groups_idx[keep_nonempty]
    
    if (length(groups) <= max_keep) {
      return(groups)
    }
    
    centers <- t(vapply(
      groups_idx,
      function(g) colMeans(emb_work[g, , drop = FALSE]),
      numeric(ncol(emb_work))
    ))
    
    covered <- rep(FALSE, nrow(emb_work))
    selected <- integer(0)
    remaining <- seq_along(groups_idx)
    dmin <- rep(Inf, length(groups_idx))
    
    for (ii in seq_len(max_keep)) {
      if (length(remaining) == 0L) break
      new_counts <- vapply(
        remaining,
        function(r) sum(!covered[groups_idx[[r]]]),
        integer(1)
      )
      best_pool <- remaining[new_counts == max(new_counts)]
      
      if (length(selected) == 0L) {
        best <- best_pool[1L]
      } else if (length(best_pool) > 1L) {
        best <- best_pool[which.max(dmin[best_pool])]
      } else {
        best <- best_pool[1L]
      }
      
      selected <- c(selected, best)
      covered[groups_idx[[best]]] <- TRUE
      remaining <- setdiff(remaining, best)
      dnew <- sqdist_to_point(centers, centers[best, ])
      dmin <- pmin(dmin, dnew)
    }
    
    groups[selected]
  }
  
  pairwise_overlap_summary <- function(groups) {
    if (length(groups) <= 1L) {
      return(list(max_overlap_cells = 0L, max_overlap_fraction = 0))
    }
    max_count <- 0L
    max_frac <- 0
    lens <- vapply(groups, length, integer(1))
    for (i in seq_len(length(groups) - 1L)) {
      gi <- groups[[i]]
      for (j in seq.int(i + 1L, length(groups))) {
        ov <- length(intersect(gi, groups[[j]]))
        if (ov > max_count) max_count <- ov
        denom <- min(lens[i], lens[j])
        frac <- if (denom > 0L) ov / denom else 0
        if (frac > max_frac) max_frac <- frac
      }
    }
    list(max_overlap_cells = max_count, max_overlap_fraction = max_frac)
  }
  
  reuse_summary <- function(long_map, n_total_cells) {
    if (nrow(long_map) == 0L) {
      return(list(
        mean_reuse_covered = 0,
        median_reuse_covered = 0,
        q95_reuse_covered = 0,
        max_reuse = 0L,
        mean_reuse_all = 0
      ))
    }
    reuse <- as.integer(table(long_map$barcode))
    list(
      mean_reuse_covered = mean(reuse),
      median_reuse_covered = stats::median(reuse),
      q95_reuse_covered = as.numeric(stats::quantile(reuse, 0.95, names = FALSE)),
      max_reuse = max(reuse),
      mean_reuse_all = nrow(long_map) / n_total_cells
    )
  }
  
  resolve_target_ct <- function(n) {
    if (isTRUE(balance_final_metacells)) {
      target_ct <- if (!is.null(target_metacells)) {
        as.integer(target_metacells)
      } else {
        500L
      }
    } else {
      if (!is.null(target_metacells)) {
        target_ct <- as.integer(target_metacells)
      } else {
        target_ct <- as.integer(ceiling(n / target_cells_per_metacell))
      }
      if (!is.null(min_metacells)) {
        target_ct <- max(target_ct, as.integer(min_metacells))
      }
      if (!is.null(max_metacells)) {
        target_ct <- min(target_ct, as.integer(max_metacells))
      }
    }
    target_ct <- min(max(1L, target_ct), n)
    target_ct
  }
  
  make_long_map <- function(groups, ct) {
    lens <- vapply(groups, length, integer(1))
    if (length(groups) == 0L || sum(lens) == 0L) {
      return(data.frame(
        barcode = character(0),
        metacell = character(0),
        celltype = character(0),
        stringsAsFactors = FALSE
      ))
    }
    data.frame(
      barcode = unlist(groups, use.names = FALSE),
      metacell = rep(names(groups), lens),
      celltype = ct,
      stringsAsFactors = FALSE
    )
  }
  
  make_primary_map <- function(groups, cell_names) {
    long <- make_long_map(groups, ct = "tmp")
    primary <- setNames(rep(NA_character_, length(cell_names)), cell_names)
    if (nrow(long) > 0L) {
      first_hit <- !duplicated(long$barcode)
      primary[long$barcode[first_hit]] <- long$metacell[first_hit]
    }
    primary
  }
  
  make_one_ct <- function(ct, ct_cells, seed_ct) {
    set.seed(seed_ct)
    
    emb <- emb_all[ct_cells, , drop = FALSE]
    n <- nrow(emb)
    
    if (n < min_cells_per_metacell) {
      warning(
        "Skipping group '", ct, "' because it has fewer than min_cells_per_metacell cells."
      )
      return(NULL)
    }
    
    cell_names <- rownames(emb)
    emb <- as.matrix(emb)
    storage.mode(emb) <- "double"
    
    emb_scale <- apply(emb, 2, stats::sd)
    emb_scale[emb_scale == 0 | !is.finite(emb_scale)] <- 1
    emb_work <- scale(emb, center = TRUE, scale = emb_scale)
    emb_work <- as.matrix(emb_work)
    emb_work[!is.finite(emb_work)] <- 0
    rownames(emb_work) <- cell_names
    
    target_ct <- resolve_target_ct(n)
    max_ct <- target_ct
    
    k_base <- as.integer(k)
    if (!is.finite(k_base) || k_base < 2L) {
      stop("k must be an integer >= 2.")
    }
    if (isTRUE(adaptive_k) && max_ct > 0L) {
      if (isTRUE(balance_final_metacells)) {
        k_base <- max(
          k_base,
          ceiling(target_memberships_per_cell * n / max_ct)
        )
      } else {
        k_base <- max(k_base, ceiling(n / max_ct))
      }
    }
    k_eff <- min(max(k_base, min_cells_per_metacell), n)
    
    if (!is.null(max_size)) {
      k_eff <- min(k_eff, as.integer(max_size), n)
      k_eff <- max(k_eff, min_cells_per_metacell)
    }
    if (!is.null(min_size)) {
      k_eff <- max(k_eff, as.integer(min_size))
      k_eff <- min(k_eff, n)
    }
    
    if (k_eff < min_cells_per_metacell) {
      warning("Skipping group '", ct, "' because k_eff < min_cells_per_metacell.")
      return(NULL)
    }
    
    requested_candidates <- max(
      as.integer(knn_iteration),
      as.integer(ceiling(target_ct * candidate_multiplier))
    )
    n_seed <- min(requested_candidates, n)
    n_seed <- max(1L, n_seed)
    
    seeds <- integer(n_seed)
    seeds[1L] <- sample(seq_len(n), 1L)
    dmin <- sqdist_to_point(emb_work, emb_work[seeds[1L], ])
    
    if (n_seed >= 2L) {
      for (i in 2:n_seed) {
        remaining <- setdiff(seq_len(n), seeds[seq_len(i - 1L)])
        if (length(remaining) == 0L) {
          seeds <- seeds[seq_len(i - 1L)]
          break
        }
        next_seed <- remaining[which.max(dmin[remaining])]
        seeds[i] <- next_seed
        d_new <- sqdist_to_point(emb_work, emb_work[next_seed, ])
        dmin <- pmin(dmin, d_new)
      }
    }
    
    seeds <- unique(seeds[seeds > 0L])
    
    nn_index <- FNN::get.knnx(
      data = emb_work,
      query = emb_work[seeds, , drop = FALSE],
      k = k_eff
    )$nn.index
    
    if (is.null(dim(nn_index))) {
      nn_index <- matrix(nn_index, nrow = 1L)
    }
    
    nn_index <- t(apply(nn_index, 1, sort))
    storage.mode(nn_index) <- "integer"
    
    candidate_groups <- lapply(seq_len(nrow(nn_index)), function(i) {
      unique(cell_names[nn_index[i, ]])
    })
    
    keep_len <- vapply(candidate_groups, length, integer(1)) >= min_cells_per_metacell
    candidate_groups <- candidate_groups[keep_len]
    nn_index <- nn_index[keep_len, , drop = FALSE]
    seeds <- seeds[keep_len]
    
    group_key <- vapply(candidate_groups, function(g) paste(g, collapse = "|"), character(1))
    keep_unique <- !duplicated(group_key)
    candidate_groups <- candidate_groups[keep_unique]
    nn_index <- nn_index[keep_unique, , drop = FALSE]
    seeds <- seeds[keep_unique]
    
    if (length(candidate_groups) == 0L) {
      warning("No candidate metacells for group '", ct, "'.")
      return(NULL)
    }
    
    overlap_threshold <- max(1L, floor(overlap_cutoff * k_eff))
    keep_flag <- determineOverlapCpp(nn_index, overlap_threshold)
    keep_idx <- which(keep_flag == 0L)
    
    if (length(keep_idx) == 0L) {
      keep_idx <- 1L
    }
    
    groups <- candidate_groups[keep_idx]
    kept_seed_ids <- seeds[keep_idx]
    
    if (isTRUE(force_target_metacells) && length(groups) < max_ct && length(candidate_groups) >= max_ct) {
      groups <- candidate_groups
    }
    
    if (length(groups) > max_ct) {
      if (identical(selection_method, "coverage")) {
        groups <- select_coverage_groups(groups, emb_work, max_ct)
      } else {
        groups <- select_diverse_groups(groups, emb_work, max_ct)
      }
    }
    
    ct_safe <- make_safe_name(ct)
    names(groups) <- paste0(ct_safe, "_", seq_along(groups))
    
    if (!isTRUE(allow_overlap)) {
      assigned <- setNames(rep(FALSE, n), cell_names)
      uniq_groups <- list()
      for (g in seq_along(groups)) {
        mem <- groups[[g]]
        new_mem <- mem[!assigned[mem]]
        if (length(new_mem) > 0L) {
          uniq_groups[[length(uniq_groups) + 1L]] <- new_mem
          assigned[new_mem] <- TRUE
        }
      }
      groups <- uniq_groups
      names(groups) <- paste0(ct_safe, "_", seq_along(groups))
      
      missing <- names(assigned)[!assigned]
      if (length(missing) > 0L && length(groups) > 0L) {
        centers_final <- t(vapply(
          groups,
          function(g) colMeans(emb_work[g, , drop = FALSE]),
          numeric(ncol(emb_work))
        ))
        nn <- FNN::get.knnx(
          data = centers_final,
          query = emb_work[missing, , drop = FALSE],
          k = 1L
        )$nn.index[, 1L]
        for (ii in seq_along(missing)) {
          groups[[nn[ii]]] <- c(groups[[nn[ii]]], missing[ii])
        }
      }
    }
    
    lens <- vapply(groups, length, integer(1))
    ov_sum <- pairwise_overlap_summary(groups)
    long_map <- make_long_map(groups, ct)
    primary_map <- make_primary_map(groups, cell_names)
    n_covered <- length(unique(long_map$barcode))
    reuse_sum <- reuse_summary(long_map, n_total_cells = n)
    
    if (verbose && n_covered / n < min_coverage_warning) {
      warning(
        "Low metacell coverage for '", ct, "': ",
        round(n_covered / n, 3),
        ". Consider increasing target_memberships_per_cell, increasing target_metacells, or using selection_method = 'coverage'."
      )
    }
    if (verbose && reuse_sum$mean_reuse_covered > max_mean_reuse_warning) {
      warning(
        "High mean cell reuse among covered cells for '", ct, "': ",
        round(reuse_sum$mean_reuse_covered, 2),
        ". Consider decreasing target_memberships_per_cell, decreasing target_metacells, or tightening overlap_cutoff."
      )
    }
    
    list(
      metacell_members = groups,
      cell2metacell = primary_map,
      cell2metacell_long = long_map,
      seeds_used = cell_names[intersect(kept_seed_ids, seq_len(n))],
      k = k_eff,
      knn_iteration = n_seed,
      overlap_cutoff = overlap_cutoff,
      overlap_threshold = overlap_threshold,
      target_metacells = target_ct,
      max_metacells = max_ct,
      min_cells_per_metacell = min_cells_per_metacell,
      force_target_metacells = force_target_metacells,
      n_cells = n,
      n_cells_covered = n_covered,
      cell_coverage_fraction = n_covered / n,
      mean_reuse_covered = reuse_sum$mean_reuse_covered,
      median_reuse_covered = reuse_sum$median_reuse_covered,
      q95_reuse_covered = reuse_sum$q95_reuse_covered,
      max_reuse = reuse_sum$max_reuse,
      mean_reuse_all = reuse_sum$mean_reuse_all,
      n_metacells = length(groups),
      size_summary = summary(lens),
      max_pairwise_overlap_cells = ov_sum$max_overlap_cells,
      max_pairwise_overlap_fraction = ov_sum$max_overlap_fraction
    )
  }
  
  celltypes <- unique(ct_vec)
  results_by_ct <- list()
  
  for (i in seq_along(celltypes)) {
    ct <- celltypes[i]
    ct_cells <- names(ct_vec)[ct_vec == ct]
    
    if (verbose) {
      cat("\nConstructing CLIPER metacells for:", ct, "\n")
      cat("n_cells:", length(ct_cells), "\n")
    }
    
    res <- make_one_ct(ct = ct, ct_cells = ct_cells, seed_ct = seed + i - 1L)
    
    if (!is.null(res)) {
      results_by_ct[[ct]] <- res
      
      if (verbose) {
        lens <- vapply(res$metacell_members, length, integer(1))
        cat("#metacells:", length(res$metacell_members), "\n")
        cat("k_eff:", res$k, "\n")
        cat("covered cells:", res$n_cells_covered, "/", res$n_cells,
            "(", round(res$cell_coverage_fraction, 3), ")\n")
        cat("max pairwise overlap fraction:", round(res$max_pairwise_overlap_fraction, 3), "\n")
        print(summary(lens))
      }
    }
  }
  
  if (length(results_by_ct) == 0L) {
    stop("No metacells were constructed.")
  }
  
  metacell_members <- unlist(
    unname(lapply(results_by_ct, function(x) x$metacell_members)),
    recursive = FALSE,
    use.names = TRUE
  )
  
  if (is.null(names(metacell_members)) || any(is.na(names(metacell_members))) || any(names(metacell_members) == "")) {
    names(metacell_members) <- paste0("metacell_", seq_along(metacell_members))
  }
  
  cell2metacell <- unlist(
    unname(lapply(results_by_ct, function(x) x$cell2metacell)),
    use.names = TRUE
  )
  
  cell2metacell_long <- do.call(
    rbind,
    unname(lapply(results_by_ct, function(x) x$cell2metacell_long))
  )
  
  seeds_used <- unlist(
    lapply(results_by_ct, function(x) x$seeds_used),
    use.names = FALSE
  )
  
  final_sizes <- vapply(metacell_members, length, integer(1))
  
  qc <- data.frame(
    metacell = names(metacell_members),
    size = as.integer(final_sizes),
    stringsAsFactors = FALSE
  )
  qc$group <- sub("_[0-9]+$", "", qc$metacell)
  
  qc_by_ct <- data.frame(
    group = names(results_by_ct),
    n_cells = vapply(results_by_ct, function(x) x$n_cells, integer(1)),
    n_cells_covered = vapply(results_by_ct, function(x) x$n_cells_covered, integer(1)),
    cell_coverage_fraction = vapply(results_by_ct, function(x) x$cell_coverage_fraction, numeric(1)),
    mean_reuse_covered = vapply(results_by_ct, function(x) x$mean_reuse_covered, numeric(1)),
    median_reuse_covered = vapply(results_by_ct, function(x) x$median_reuse_covered, numeric(1)),
    q95_reuse_covered = vapply(results_by_ct, function(x) x$q95_reuse_covered, numeric(1)),
    max_reuse = vapply(results_by_ct, function(x) x$max_reuse, integer(1)),
    mean_reuse_all = vapply(results_by_ct, function(x) x$mean_reuse_all, numeric(1)),
    n_metacells = vapply(results_by_ct, function(x) x$n_metacells, integer(1)),
    k_eff = vapply(results_by_ct, function(x) as.integer(x$k), integer(1)),
    overlap_cutoff = vapply(results_by_ct, function(x) x$overlap_cutoff, numeric(1)),
    max_pairwise_overlap_fraction = vapply(results_by_ct, function(x) x$max_pairwise_overlap_fraction, numeric(1)),
    stringsAsFactors = FALSE
  )
  
  if (verbose) {
    cat("\nCLIPER metacell construction summary\n")
    cat("#metacells:", length(metacell_members), "\n")
    cat("long barcode-metacell memberships:", nrow(cell2metacell_long), "\n")
    print(summary(final_sizes))
    print(qc_by_ct)
  }
  
  list(
    metacell_members = metacell_members,
    cell2metacell = cell2metacell,
    cell2metacell_long = cell2metacell_long,
    seeds_used = seeds_used,
    k = k,
    knn_iteration = knn_iteration,
    overlap_cutoff = overlap_cutoff,
    target_metacells = target_metacells,
    max_metacells = max_metacells,
    min_metacells = min_metacells,
    target_cells_per_metacell = target_cells_per_metacell,
    min_cells_per_metacell = min_cells_per_metacell,
    allow_overlap = allow_overlap,
    balance_final_metacells = balance_final_metacells,
    force_target_metacells = force_target_metacells,
    target_memberships_per_cell = target_memberships_per_cell,
    selection_method = selection_method,
    overlap_threshold = vapply(results_by_ct, function(x) as.integer(x$overlap_threshold), integer(1)),
    results_by_ct = results_by_ct,
    qc = qc,
    qc_by_ct = qc_by_ct
  )
}

aggregate_metacell_counts <- function(
    counts,
    cell2metacell = NULL,
    metacell_members = NULL,
    cells_are_columns = TRUE
) {
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("Package 'Matrix' is required.")
  }
  
  if (!is.null(metacell_members)) {
    if (length(metacell_members) == 0L) {
      stop("metacell_members is empty.")
    }
    if (is.null(names(metacell_members)) || any(names(metacell_members) == "")) {
      names(metacell_members) <- paste0("metacell_", seq_along(metacell_members))
    }
    
    cell_names <- if (cells_are_columns) colnames(counts) else rownames(counts)
    if (is.null(cell_names)) {
      stop("counts must have cell names in colnames(counts) or rownames(counts).")
    }
    
    metacell_members <- lapply(metacell_members, function(g) {
      unique(intersect(as.character(g), cell_names))
    })
    keep_groups <- vapply(metacell_members, length, integer(1)) > 0L
    metacell_members <- metacell_members[keep_groups]
    
    if (length(metacell_members) == 0L) {
      stop("No metacell members overlap the count matrix cell names.")
    }
    
    lens <- vapply(metacell_members, length, integer(1))
    i <- match(unlist(metacell_members, use.names = FALSE), cell_names)
    j <- rep(seq_along(metacell_members), lens)
    
    group_mat <- Matrix::sparseMatrix(
      i = i,
      j = j,
      x = 1,
      dims = c(length(cell_names), length(metacell_members)),
      dimnames = list(cell_names, names(metacell_members))
    )
    
    if (cells_are_columns) {
      counts_use <- counts[, cell_names, drop = FALSE]
      out <- counts_use %*% group_mat
    } else {
      counts_use <- counts[cell_names, , drop = FALSE]
      out <- Matrix::t(group_mat) %*% counts_use
    }
    
    return(out)
  }
  
  if (is.null(cell2metacell)) {
    stop("Provide either cell2metacell or metacell_members.")
  }
  
  if (cells_are_columns) {
    keep_cells <- intersect(colnames(counts), names(cell2metacell))
    counts <- counts[, keep_cells, drop = FALSE]
    group <- cell2metacell[keep_cells]
    
    mm <- Matrix::sparse.model.matrix(~ 0 + factor(group))
    colnames(mm) <- levels(factor(group))
    
    out <- counts %*% mm
  } else {
    keep_cells <- intersect(rownames(counts), names(cell2metacell))
    counts <- counts[keep_cells, , drop = FALSE]
    group <- cell2metacell[keep_cells]
    
    mm <- Matrix::sparse.model.matrix(~ 0 + factor(group))
    colnames(mm) <- levels(factor(group))
    
    out <- Matrix::t(mm) %*% counts
  }
  
  out
}

normalize_delta <- function(
    counts,
    qc_cutoff = 0.05
) {
  
  lib_size <- colSums(counts)
  L <- mean(lib_size)
  
  size_factors <- lib_size / L
  norm_counts <- sweep(counts, 2, size_factors, "/")
  log_norm_counts <- log1p(norm_counts)
  
  n_cells <- ncol(log_norm_counts)
  keep <- rowSums(counts > 2) >= ceiling(qc_cutoff * n_cells)
  log_norm_counts <- log_norm_counts[keep, , drop = FALSE]
  
  
  return(log_norm_counts)
}

normalize_metacell_deseq2 <- function(
    counts,
    qc_cutoff = 0.05,
    min_total_count_per_cell = 1
){
  stopifnot(is.matrix(counts) || inherits(counts, "dgCMatrix"))
  
  keep_cols <- colSums(counts) >= min_total_count_per_cell
  counts <- counts[, keep_cols, drop = FALSE]
  if (ncol(counts) == 0) stop("All columns dropped (all metacells have 0 counts).")
  
  n_cells <- ncol(counts)
  keep_rows <- rowSums(counts > 2) >= ceiling(qc_cutoff * n_cells)
  counts <- counts[keep_rows, , drop = FALSE]
  if (nrow(counts) == 0) stop("All rows dropped by qc_cutoff filter.")
  
  if (inherits(counts, "dgCMatrix")) {
    nnz <- counts@x
    if (any(abs(nnz - round(nnz)) > 0)) {
      stop("Input 'counts' contains non-integers (sparse matrix). ",
           "DESeq2 requires raw integer counts. ",
           "This usually means you passed normalized/averaged data instead of summed counts.")
    }
    counts@x <- as.numeric(round(counts@x))
  } else {
    if (any(abs(counts - round(counts)) > 0, na.rm = TRUE)) {
      stop("Input 'counts' contains non-integers. ",
           "DESeq2 requires raw integer counts. ",
           "This usually means you passed normalized/averaged data instead of summed counts.")
    }
    storage.mode(counts) <- "integer"
  }
  
  meta_info <- data.frame(dummy = factor(rep("all", ncol(counts))))
  rownames(meta_info) <- colnames(counts)
  
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = counts,
    colData = meta_info,
    design = ~ 1
  )
  
  dds <- DESeq2::estimateSizeFactors(dds)
  
  base_mean <- Matrix::rowMeans(DESeq2::counts(dds, normalized = TRUE))
  
  if (sum(base_mean > 5, na.rm = TRUE) < 1000) {
    message("[CLIPER] Using varianceStabilizingTransformation() because too few genes have mean normalized count > 5.")
    mat <- SummarizedExperiment::assay(
      DESeq2::varianceStabilizingTransformation(dds, blind = TRUE)
    )
  } else {
    mat <- SummarizedExperiment::assay(
      DESeq2::vst(dds, blind = TRUE)
    )
  }
  
  mat
}

normalize_metacell_logcpm <- function(
    counts,
    min_total_count_per_cell = 1,
    qc_cutoff = 0.05,
    norm_method = "TMM",
    prior_count = 1
){
  keep_cols <- Matrix::colSums(counts) >= min_total_count_per_cell
  counts <- counts[, keep_cols, drop = FALSE]
  if (ncol(counts) == 0) stop("All columns dropped (all metacells have 0 counts).")
  
  n_cells <- ncol(counts)
  keep_rows <- Matrix::rowSums(counts > 2) >= ceiling(qc_cutoff * n_cells)
  counts <- counts[keep_rows, , drop = FALSE]
  if (nrow(counts) == 0) stop("All rows dropped by qc_cutoff filter.")
  
  y <- edgeR::DGEList(counts = counts)
  if (!identical(norm_method, "none")) {
    y <- edgeR::calcNormFactors(y, method = norm_method)
  }
  mat <- edgeR::cpm(y, log = TRUE, prior.count = prior_count)
  
  mat
}

gc_fq_normalize <- function(
    counts,
    genome = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38,
    nbins = 20,
    qc_cutoff = 0.05,
    keep_chr = c(paste0("chr", 1:22), "chrX"),
    bed0 = FALSE,
    return_log1p = TRUE,
    verbose = TRUE
){
  stopifnot(is.matrix(counts) || inherits(counts, "dgCMatrix"))
  if (is.null(rownames(counts))) stop("counts must have rownames as peak strings like chr-start-end or chr:start-end")
  
  make_peaks_gr <- function(peaks, bed0 = FALSE) {
    spl1 <- strsplit(peaks, "[:\\-]")
    chr  <- vapply(spl1, `[`, "", 1)
    st   <- as.integer(vapply(spl1, `[`, "", 2))
    en   <- as.integer(vapply(spl1, `[`, "", 3))
    if (isTRUE(bed0)) st <- st + 1L
    GenomicRanges::GRanges(seqnames = chr, ranges = IRanges::IRanges(start = st, end = en))
  }
  
  qn_nogroup <- function(X) {
    X <- as.matrix(X)
    rn <- rownames(X); cn <- colnames(X)
    
    Qref <- rowMeans(apply(X, 2, sort, na.last = NA))
    
    out <- matrix(NA_real_, nrow(X), ncol(X), dimnames = list(rn, cn))
    for (j in seq_len(ncol(X))) {
      x <- X[, j]
      
      o <- order(x, na.last = TRUE)
      y <- rep(NA_real_, length(x))
      y[o] <- Qref
      
      ok <- !is.na(x)
      if (any(ok)) {
        y_ok <- y[ok]
        x_ok <- x[ok]
        y[ok] <- ave(y_ok, x_ok, FUN = mean)
      }
      out[, j] <- y
    }
    out
  }
  
  make_gc_bins <- function(gc, nbins = 20) {
    qs <- stats::quantile(gc, probs = seq(0, 1, length.out = nbins + 1), na.rm = TRUE, names = FALSE)
    qs <- unique(as.numeric(qs))
    if (length(qs) < 3L) stop("GC-content has too few unique values to make bins. Try smaller nbins or check peaks.")
    base::cut(gc, breaks = qs, include.lowest = TRUE, labels = FALSE)
  }
  
  n_cells <- ncol(counts)
  keep_detect <- Matrix::rowSums(counts > 0) >= ceiling(qc_cutoff * n_cells)
  if (verbose) message("Keep by detection: ", sum(keep_detect), " / ", length(keep_detect))
  counts <- counts[keep_detect, , drop = FALSE]
  
  peaks_gr_all <- make_peaks_gr(rownames(counts), bed0 = bed0)
  GenomeInfoDb::seqlevelsStyle(peaks_gr_all) <- "UCSC"
  
  keep_chr_idx <- as.character(GenomicRanges::seqnames(peaks_gr_all)) %in% keep_chr
  if (verbose) message("Keep by chr: ", sum(keep_chr_idx), " / ", length(keep_chr_idx))
  
  peaks_gr <- peaks_gr_all[keep_chr_idx]
  counts  <- counts[keep_chr_idx, , drop = FALSE]
  
  seqs <- Biostrings::getSeq(genome, peaks_gr)
  gc_mat  <- Biostrings::letterFrequency(seqs, letters = c("G", "C"), as.prob = TRUE)
  gc_frac <- rowSums(gc_mat)
  stopifnot(length(gc_frac) == nrow(counts))
  
  gc_bins <- make_gc_bins(gc_frac, nbins = nbins)
  gc_bins <- as.factor(gc_bins)
  
  if (verbose) {
    tab <- table(gc_bins)
    message("GC bins: ", length(tab), " (min bin size = ", min(tab), ")")
  }
  
  X <- as.matrix(counts)
  
  out <- matrix(NA_real_, nrow(X), ncol(X), dimnames = dimnames(X))
  for (lvl in levels(gc_bins)) {
    idx <- which(gc_bins == lvl)
    if (length(idx) <= 1L) {
      out[idx, ] <- X[idx, , drop = FALSE]
    } else {
      out[idx, ] <- qn_nogroup(X[idx, , drop = FALSE])
    }
  }
  
  if (return_log1p) {
    out <- log1p(out)
  }
  
  
  return(out)
}

data_norm <- function(
    RNA_meta,
    ATAC_meta,
    norm_RNA = c("DESeq2", "Delta", "CPM"),
    norm_ATAC = c("Delta", "DESeq2", "CPM", "GC-FQ"),
    qc_cutoff_RNA = 0.1,
    qc_cutoff_ATAC = 0.05
){
  norm_RNA <- match.arg(norm_RNA)
  norm_ATAC <- match.arg(norm_ATAC)
  if (!is.finite(qc_cutoff_RNA) || qc_cutoff_RNA < 0 || qc_cutoff_RNA > 1) {
    stop("qc_cutoff_RNA must be between 0 and 1.")
  }
  if (!is.finite(qc_cutoff_ATAC) || qc_cutoff_ATAC < 0 || qc_cutoff_ATAC > 1) {
    stop("qc_cutoff_ATAC must be between 0 and 1.")
  }
  
  if (norm_RNA == "DESeq2") {
    meta_rna_norm <- normalize_metacell_deseq2(t(RNA_meta), qc_cutoff = qc_cutoff_RNA)
  } else if (norm_RNA == "Delta") {
    meta_rna_norm <- normalize_delta(t(RNA_meta), qc_cutoff = qc_cutoff_RNA)
  } else if (norm_RNA == "CPM") {
    meta_rna_norm <- normalize_metacell_logcpm(t(RNA_meta), qc_cutoff = qc_cutoff_RNA)
  }
  
  if (norm_ATAC == "DESeq2") {
    meta_atac_norm <- normalize_metacell_deseq2(t(ATAC_meta), qc_cutoff = qc_cutoff_ATAC)
  } else if (norm_ATAC == "Delta") {
    meta_atac_norm <- normalize_delta(t(ATAC_meta), qc_cutoff = qc_cutoff_ATAC)
  } else if (norm_ATAC == "CPM") {
    meta_atac_norm <- normalize_metacell_logcpm(t(ATAC_meta), qc_cutoff = qc_cutoff_ATAC)
  } else if (norm_ATAC == "GC-FQ") {
    meta_atac_norm <- gc_fq_normalize(t(ATAC_meta), qc_cutoff = qc_cutoff_ATAC)
  }
  
  list(meta_rna = meta_rna_norm, meta_atac = meta_atac_norm)
}


score_cells_for_downsampling <- function(
    meta,
    cells,
    quality_cols = NULL,
    higher_is_better_cols = NULL,
    lower_is_better_cols = NULL
) {
  cells <- intersect(cells, rownames(meta))
  if (length(cells) == 0L) {
    return(setNames(numeric(0), character(0)))
  }
  
  if (is.null(quality_cols) || length(quality_cols) == 0L) {
    return(setNames(rep(NA_real_, length(cells)), cells))
  }
  
  quality_cols <- intersect(quality_cols, colnames(meta))
  if (length(quality_cols) == 0L) {
    warning("None of quality_cols are present in obj@meta.data; skipping quality-score filtering.")
    return(setNames(rep(NA_real_, length(cells)), cells))
  }
  
  meta_qc <- meta[cells, quality_cols, drop = FALSE]
  
  infer_lower_better <- function(x) {
    x_low <- tolower(x)
    grepl("nucleosome|percent\\.?mt|pct\\.?mt|mito|blacklist|doublet", x_low)
  }
  
  if (is.null(lower_is_better_cols)) {
    lower_is_better_cols <- quality_cols[vapply(quality_cols, infer_lower_better, logical(1))]
  } else {
    lower_is_better_cols <- intersect(lower_is_better_cols, quality_cols)
  }
  
  if (is.null(higher_is_better_cols)) {
    higher_is_better_cols <- setdiff(quality_cols, lower_is_better_cols)
  } else {
    higher_is_better_cols <- intersect(higher_is_better_cols, quality_cols)
  }
  
  score_mat <- lapply(quality_cols, function(cc) {
    x <- suppressWarnings(as.numeric(meta_qc[[cc]]))
    if (all(is.na(x))) {
      return(rep(NA_real_, length(x)))
    }
    
    if (all(x >= 0, na.rm = TRUE)) {
      x_use <- log1p(x)
    } else {
      x_use <- x
    }
    
    sx <- stats::sd(x_use, na.rm = TRUE)
    if (!is.finite(sx) || sx == 0) {
      z <- rep(0, length(x_use))
    } else {
      z <- as.numeric(scale(x_use))
    }
    
    if (cc %in% lower_is_better_cols) {
      z <- -z
    }
    z
  })
  
  score_mat <- do.call(cbind, score_mat)
  score <- rowMeans(score_mat, na.rm = TRUE)
  score[!is.finite(score)] <- NA_real_
  names(score) <- cells
  score
}

sample_cells_for_metacell <- function(
    meta,
    cells,
    max_cells = 10000,
    strata_cols = NULL,
    quality_cols = NULL,
    quality_filter_quantile = NULL,
    higher_is_better_cols = NULL,
    lower_is_better_cols = NULL,
    min_cells_per_stratum = 20,
    seed = 1,
    verbose = TRUE
) {
  set.seed(seed)
  
  cells <- unique(intersect(cells, rownames(meta)))
  n_original <- length(cells)
  
  if (n_original == 0L) {
    return(list(
      cells = character(0),
      summary = data.frame(
        n_cells_original = 0L,
        n_cells_after_quality_filter = 0L,
        n_cells_sampled = 0L,
        max_cells = ifelse(is.null(max_cells), NA_integer_, as.integer(max_cells)),
        quality_filter_used = FALSE,
        downsample_used = FALSE,
        stringsAsFactors = FALSE
      )
    ))
  }
  
  if (is.null(max_cells) || !is.finite(max_cells)) {
    max_cells <- n_original
  }
  max_cells <- as.integer(max_cells)
  if (max_cells < 1L) {
    stop("max_cells must be positive, finite, or NULL.")
  }
  
  strata_cols <- strata_cols[!is.na(strata_cols) & nzchar(strata_cols)]
  strata_cols <- intersect(strata_cols, colnames(meta))
  
  df <- data.frame(
    cell = cells,
    stringsAsFactors = FALSE
  )
  
  if (length(strata_cols) == 0L) {
    df$.stratum <- "all"
  } else {
    strata_df <- meta[cells, strata_cols, drop = FALSE]
    strata_df[] <- lapply(strata_df, function(x) {
      x <- as.character(x)
      x[is.na(x) | x == ""] <- "NA"
      x
    })
    df$.stratum <- apply(strata_df, 1, paste, collapse = "__")
  }
  
  min_cells_per_stratum <- as.integer(min_cells_per_stratum)
  if (!is.finite(min_cells_per_stratum) || min_cells_per_stratum < 1L) {
    min_cells_per_stratum <- 1L
  }
  if (min_cells_per_stratum > 1L) {
    tab0 <- table(df$.stratum)
    small_strata <- names(tab0)[tab0 < min_cells_per_stratum]
    if (length(small_strata) > 0L && length(small_strata) < length(tab0)) {
      df$.stratum[df$.stratum %in% small_strata] <- "small_strata_pooled"
    }
  }
  
  quality_filter_used <- FALSE
  if (!is.null(quality_filter_quantile)) {
    if (!is.finite(quality_filter_quantile) || quality_filter_quantile < 0 || quality_filter_quantile >= 1) {
      stop("quality_filter_quantile must be in [0, 1), or NULL.")
    }
    
    qc_score <- score_cells_for_downsampling(
      meta = meta,
      cells = df$cell,
      quality_cols = quality_cols,
      higher_is_better_cols = higher_is_better_cols,
      lower_is_better_cols = lower_is_better_cols
    )
    
    if (any(is.finite(qc_score))) {
      cutoff <- stats::quantile(qc_score, probs = quality_filter_quantile, na.rm = TRUE, names = FALSE)
      keep_qc <- is.na(qc_score[df$cell]) | qc_score[df$cell] >= cutoff
      
      if (sum(keep_qc) >= min(max_cells, n_original)) {
        df <- df[keep_qc, , drop = FALSE]
        quality_filter_used <- TRUE
      } else if (verbose) {
        warning(
          "Skipping quality_filter_quantile because it would leave fewer cells than needed for max_cells."
        )
      }
    }
  }
  
  n_after_qc <- nrow(df)
  
  if (n_after_qc <= max_cells) {
    return(list(
      cells = df$cell,
      summary = data.frame(
        n_cells_original = n_original,
        n_cells_after_quality_filter = n_after_qc,
        n_cells_sampled = n_after_qc,
        max_cells = max_cells,
        n_strata = length(unique(df$.stratum)),
        strata_cols = ifelse(length(strata_cols) == 0L, "none", paste(strata_cols, collapse = ",")),
        quality_filter_used = quality_filter_used,
        downsample_used = FALSE,
        stringsAsFactors = FALSE
      )
    ))
  }
  
  tab <- table(df$.stratum)
  strata <- names(tab)
  n_strata <- length(strata)
  
  if (n_strata > max_cells) {
    sampled <- sample(df$cell, size = max_cells, replace = FALSE)
  } else {
    target_raw <- as.numeric(tab) / sum(tab) * max_cells
    target_n <- floor(target_raw)
    
    target_n[target_n == 0L & as.integer(tab) > 0L] <- 1L
    target_n <- pmin(target_n, as.integer(tab))
    
    while (sum(target_n) < max_cells) {
      capacity <- as.integer(tab) - target_n
      if (all(capacity <= 0L)) break
      residual <- target_raw - floor(target_raw)
      residual[capacity <= 0L] <- -Inf
      jj <- which.max(residual)
      target_n[jj] <- target_n[jj] + 1L
    }
    
    while (sum(target_n) > max_cells) {
      removable <- target_n > 1L
      if (!any(removable)) break
      jj <- which.max(ifelse(removable, target_n, -Inf))
      target_n[jj] <- target_n[jj] - 1L
    }
    
    sampled <- unlist(lapply(seq_along(strata), function(ii) {
      cc <- df$cell[df$.stratum == strata[ii]]
      nn <- min(length(cc), target_n[ii])
      if (nn <= 0L) return(character(0))
      sample(cc, size = nn, replace = FALSE)
    }), use.names = FALSE)
    
    if (length(sampled) < max_cells) {
      add_pool <- setdiff(df$cell, sampled)
      if (length(add_pool) > 0L) {
        sampled <- c(sampled, sample(add_pool, size = min(max_cells - length(sampled), length(add_pool))))
      }
    }
    if (length(sampled) > max_cells) {
      sampled <- sample(sampled, size = max_cells, replace = FALSE)
    }
  }
  
  if (verbose) {
    message(
      "Representative downsampling for metacells: ",
      n_original, " -> ", length(sampled), " cells",
      if (length(strata_cols) > 0L) paste0("; strata = ", paste(strata_cols, collapse = ",")) else ""
    )
  }
  
  list(
    cells = sampled,
    summary = data.frame(
      n_cells_original = n_original,
      n_cells_after_quality_filter = n_after_qc,
      n_cells_sampled = length(sampled),
      max_cells = max_cells,
      n_strata = n_strata,
      strata_cols = ifelse(length(strata_cols) == 0L, "none", paste(strata_cols, collapse = ",")),
      quality_filter_used = quality_filter_used,
      downsample_used = TRUE,
      stringsAsFactors = FALSE
    )
  )
}

Create_Signac_CLIPER_obj <- function(
    signac_obj = NULL,
    RNA = "RNA",
    ATAC = "Peaks",
    Barcodes_col = "orig.ident",
    norm_RNA = "DESeq2",
    norm_ATAC = "GC-FQ",
    qc_cutoff_RNA = 0.1,
    qc_cutoff_ATAC = 0.05,
    celltype = "celltype",
    target_metacells = NULL,
    max_metacells = 500,
    min_metacells = 200,
    target_cells_per_metacell = 15,
    min_cells_per_metacell = 10,
    overlap_cutoff = 0.3,
    k = 10,
    knn_iteration = 500,
    candidate_multiplier = 5,
    adaptive_k = TRUE,
    balance_final_metacells = FALSE,
    force_target_metacells = FALSE,
    target_memberships_per_cell = 1.2,
    selection_method = c("coverage", "diverse"),
    allow_overlap = TRUE,
    min_coverage_warning = 0.75,
    max_mean_reuse_warning = 3,
    min_cells_per_celltype = 200,
    downsample_large_celltypes = TRUE,
    max_cells_per_celltype = 10000,
    downsample_strata_cols = NULL,
    downsample_quality_cols = NULL,
    downsample_quality_filter_quantile = NULL,
    downsample_higher_is_better_cols = NULL,
    downsample_lower_is_better_cols = NULL,
    downsample_min_cells_per_stratum = 20,
    seed = 1
) {
  
  selection_method <- match.arg(selection_method)
  obj_meta <- signac_obj@meta.data
  cts <- unique(stats::na.omit(as.character(obj_meta[[celltype]])))
  
  if (isTRUE(downsample_large_celltypes)) {
    if (is.null(max_cells_per_celltype) || !is.finite(max_cells_per_celltype) || max_cells_per_celltype < 1L) {
      stop("max_cells_per_celltype must be a positive finite number when downsample_large_celltypes = TRUE.")
    }
    max_cells_per_celltype <- as.integer(max_cells_per_celltype)
  }
  
  if (is.null(downsample_strata_cols)) {
    downsample_strata_cols_use <- Barcodes_col
  } else {
    downsample_strata_cols_use <- downsample_strata_cols
  }
  downsample_strata_cols_use <- intersect(downsample_strata_cols_use, colnames(obj_meta))
  
  metacell <- setNames(
    lapply(seq_along(cts), function(ii) {
      ct <- cts[ii]
      
      cells_ct_all <- rownames(obj_meta)[as.character(obj_meta[[celltype]]) == ct]
      cells_ct_all <- cells_ct_all[!is.na(cells_ct_all)]
      
      if (length(cells_ct_all) < min_cells_per_celltype) {
        warning("Skipping ", ct, ": fewer than ", min_cells_per_celltype, " cells.")
        return(NULL)
      }
      
      if (isTRUE(downsample_large_celltypes)) {
        sampled_ct <- sample_cells_for_metacell(
          meta = obj_meta,
          cells = cells_ct_all,
          max_cells = max_cells_per_celltype,
          strata_cols = downsample_strata_cols_use,
          quality_cols = downsample_quality_cols,
          quality_filter_quantile = downsample_quality_filter_quantile,
          higher_is_better_cols = downsample_higher_is_better_cols,
          lower_is_better_cols = downsample_lower_is_better_cols,
          min_cells_per_stratum = downsample_min_cells_per_stratum,
          seed = seed + ii - 1L,
          verbose = FALSE
        )
        cells_ct <- sampled_ct$cells
        sampling_summary <- sampled_ct$summary
      } else {
        cells_ct <- cells_ct_all
        sampling_summary <- data.frame(
          n_cells_original = length(cells_ct_all),
          n_cells_after_quality_filter = length(cells_ct_all),
          n_cells_sampled = length(cells_ct_all),
          max_cells = NA_integer_,
          n_strata = NA_integer_,
          strata_cols = "none",
          quality_filter_used = FALSE,
          downsample_used = FALSE,
          stringsAsFactors = FALSE
        )
      }
      
      if (length(cells_ct) < min_cells_per_celltype) {
        warning("Skipping ", ct, ": fewer than ", min_cells_per_celltype, " cells after downsampling/quality filtering.")
        return(NULL)
      }
      
      subset_obj <- subset(signac_obj, cells = cells_ct)
      
      DefaultAssay(subset_obj) <- ATAC
      
      if ("lsi" %in% names(subset_obj@reductions)) {
        subset_obj[["lsi"]] <- NULL
      }
      subset_obj <- RunTFIDF(subset_obj)
      subset_obj <- FindTopFeatures(subset_obj, min.cutoff = "q0")
      subset_obj <- RunSVD(subset_obj)
      
      mc <- make_metacell_groups(
        obj = subset_obj,
        reduction = "lsi",
        dims_to_use = 2:30,
        celltype_col = NULL,
        celltype_value = ct,
        k = k,
        knn_iteration = knn_iteration,
        overlap_cutoff = overlap_cutoff,
        target_metacells = target_metacells,
        max_metacells = max_metacells,
        min_metacells = min_metacells,
        target_cells_per_metacell = target_cells_per_metacell,
        min_cells_per_metacell = min_cells_per_metacell,
        candidate_multiplier = candidate_multiplier,
        adaptive_k = adaptive_k,
        balance_final_metacells = balance_final_metacells,
        force_target_metacells = force_target_metacells,
        target_memberships_per_cell = target_memberships_per_cell,
        selection_method = selection_method,
        allow_overlap = allow_overlap,
        min_coverage_warning = min_coverage_warning,
        max_mean_reuse_warning = max_mean_reuse_warning,
        seed = seed + ii - 1L,
        verbose = FALSE
      )
      
      subset_obj$metacell_id_primary <- unname(mc$cell2metacell[colnames(subset_obj)])
      
      RNA_counts  <- GetAssayData(subset_obj, assay = RNA,  layer = "counts")
      ATAC_counts <- GetAssayData(subset_obj, assay = ATAC, layer = "counts")
      
      sc_gene_detect <- Matrix::rowMeans(RNA_counts > 0)
      keep_genes_sc <- names(sc_gene_detect)[sc_gene_detect >= qc_cutoff_RNA]
      RNA_counts <- RNA_counts[keep_genes_sc, , drop = FALSE]
      
      meta_rna <- t(aggregate_metacell_counts(
        counts = RNA_counts,
        metacell_members = mc$metacell_members,
        cells_are_columns = TRUE
      ))
      
      meta_atac <- t(aggregate_metacell_counts(
        counts = ATAC_counts,
        metacell_members = mc$metacell_members,
        cells_are_columns = TRUE
      ))
      
      data_normalized <- data_norm(
        RNA_meta  = meta_rna,
        ATAC_meta = meta_atac,
        norm_RNA  = norm_RNA,
        norm_ATAC = norm_ATAC,
        qc_cutoff_RNA = qc_cutoff_RNA,
        qc_cutoff_ATAC = qc_cutoff_ATAC
      )
      
      message(ct)
      
      barcode_map <- mc$cell2metacell_long
      if (!Barcodes_col %in% colnames(subset_obj@meta.data)) {
        barcode_map[[Barcodes_col]] <- NA_character_
      } else {
        barcode_map[[Barcodes_col]] <- subset_obj@meta.data[barcode_map$barcode, Barcodes_col, drop = TRUE]
      }
      
      list(
        meta_rna          = data_normalized[["meta_rna"]],
        meta_atac         = data_normalized[["meta_atac"]],
        barcode_map       = barcode_map,
        metacell_sampling = sampling_summary,
        metacell_members  = mc$metacell_members,
        metacell_qc       = mc$qc,
        metacell_qc_by_ct = mc$qc_by_ct,
        metacell_settings = list(
          target_metacells = target_metacells,
          max_metacells = max_metacells,
          min_metacells = min_metacells,
          qc_cutoff_RNA = qc_cutoff_RNA,
          qc_cutoff_ATAC = qc_cutoff_ATAC,
          target_cells_per_metacell = target_cells_per_metacell,
          resolved_target_metacells = mc$results_by_ct[[ct]]$target_metacells,
          min_cells_per_metacell = min_cells_per_metacell,
          overlap_cutoff = overlap_cutoff,
          k = k,
          k_eff = mc$results_by_ct[[ct]]$k,
          knn_iteration = knn_iteration,
          candidate_multiplier = candidate_multiplier,
          adaptive_k = adaptive_k,
          balance_final_metacells = balance_final_metacells,
          force_target_metacells = force_target_metacells,
          target_memberships_per_cell = target_memberships_per_cell,
          selection_method = selection_method,
          allow_overlap = allow_overlap,
          downsample_large_celltypes = downsample_large_celltypes,
          max_cells_per_celltype = max_cells_per_celltype,
          downsample_strata_cols = downsample_strata_cols_use,
          downsample_quality_cols = downsample_quality_cols,
          downsample_quality_filter_quantile = downsample_quality_filter_quantile
        )
      )
    }),
    cts
  )
  
  metacell <- Filter(Negate(is.null), metacell)
  if (length(metacell) == 0L) {
    stop("No cell types produced a CLIPER object. Check cell counts and celltype labels.")
  }
  
  metacell
}

make_metacell_groups_from_emb <- function(
    emb,
    ct_name = "celltype",
    k = 10,
    knn_iteration = 500,
    overlap_cutoff = 0.3,
    target_metacells = NULL,
    max_metacells = 500,
    min_metacells = 200,
    target_cells_per_metacell = 15,
    min_cells_per_metacell = 10,
    candidate_multiplier = 5,
    adaptive_k = TRUE,
    balance_final_metacells = FALSE,
    force_target_metacells = FALSE,
    target_memberships_per_cell = 1.2,
    selection_method = c("coverage", "diverse"),
    allow_overlap = TRUE,
    seed = 1,
    verbose = FALSE
){
  if (!exists("determineOverlapCpp")) {
    stop("determineOverlapCpp not found in environment. ",
         "Please load the package / source the C++ function you used in Signac.")
  }
  if (!requireNamespace("FNN", quietly = TRUE)) {
    stop("Package 'FNN' is required.")
  }
  selection_method <- match.arg(selection_method)
  if (!is.logical(force_target_metacells) || length(force_target_metacells) != 1L) {
    stop("force_target_metacells must be TRUE or FALSE.")
  }
  if (!is.finite(target_memberships_per_cell) || target_memberships_per_cell <= 0) {
    stop("target_memberships_per_cell must be positive.")
  }
  set.seed(seed)
  
  if (is.null(rownames(emb)) || any(rownames(emb) == "")) {
    stop("emb must have non-empty rownames as cell names.")
  }
  
  emb <- as.matrix(emb)
  storage.mode(emb) <- "double"
  emb[!is.finite(emb)] <- 0
  
  n <- nrow(emb)
  if (n < min_cells_per_metacell) {
    stop("Need at least min_cells_per_metacell cells to form metacells.")
  }
  
  if (isTRUE(balance_final_metacells)) {
    target_ct <- if (!is.null(target_metacells)) as.integer(target_metacells) else 500L
  } else {
    if (!is.null(target_metacells)) {
      target_ct <- as.integer(target_metacells)
    } else {
      target_ct <- as.integer(ceiling(n / target_cells_per_metacell))
    }
    if (!is.null(min_metacells)) target_ct <- max(target_ct, as.integer(min_metacells))
    if (!is.null(max_metacells)) target_ct <- min(target_ct, as.integer(max_metacells))
  }
  target_ct <- min(max(1L, target_ct), n)
  max_ct <- target_ct
  
  k_base <- as.integer(k)
  if (isTRUE(adaptive_k) && max_ct > 0L) {
    if (isTRUE(balance_final_metacells)) {
      k_base <- max(k_base, ceiling(target_memberships_per_cell * n / max_ct))
    } else {
      k_base <- max(k_base, ceiling(n / max_ct))
    }
  }
  k_eff <- min(max(k_base, as.integer(min_cells_per_metacell)), n)
  
  n_seed <- min(
    max(as.integer(knn_iteration), as.integer(ceiling(target_ct * candidate_multiplier))),
    n
  )
  
  sqdist_to_point <- function(mat, point) {
    rowSums((mat - matrix(point, nrow(mat), length(point), byrow = TRUE))^2)
  }
  
  seeds <- integer(n_seed)
  seeds[1L] <- sample(seq_len(n), 1L)
  dmin <- sqdist_to_point(emb, emb[seeds[1L], ])
  if (n_seed >= 2L) {
    for (i in 2:n_seed) {
      remaining <- setdiff(seq_len(n), seeds[seq_len(i - 1L)])
      if (length(remaining) == 0L) {
        seeds <- seeds[seq_len(i - 1L)]
        break
      }
      seeds[i] <- remaining[which.max(dmin[remaining])]
      dmin <- pmin(dmin, sqdist_to_point(emb, emb[seeds[i], ]))
    }
  }
  
  nn_index <- FNN::get.knnx(
    data  = emb,
    query = emb[seeds, , drop = FALSE],
    k     = k_eff
  )$nn.index
  nn_index <- t(apply(nn_index, 1, sort))
  storage.mode(nn_index) <- "integer"
  
  select_groups_from_emb <- function(groups, max_keep, method = c("coverage", "diverse")) {
    method <- match.arg(method)
    if (length(groups) <= max_keep) return(groups)
    centers <- t(vapply(groups, function(g) colMeans(emb[g, , drop = FALSE]), numeric(ncol(emb))))
    selected <- integer(0)
    remaining <- seq_along(groups)
    dmin <- rep(Inf, length(groups))
    covered <- setNames(rep(FALSE, n), rownames(emb))
    
    for (ii in seq_len(max_keep)) {
      if (length(remaining) == 0L) break
      if (identical(method, "coverage")) {
        new_counts <- vapply(remaining, function(r) sum(!covered[groups[[r]]]), integer(1))
        best_pool <- remaining[new_counts == max(new_counts)]
      } else {
        best_pool <- remaining
      }
      if (length(selected) == 0L) {
        best <- best_pool[1L]
      } else if (length(best_pool) > 1L) {
        best <- best_pool[which.max(dmin[best_pool])]
      } else {
        best <- best_pool[1L]
      }
      selected <- c(selected, best)
      covered[groups[[best]]] <- TRUE
      remaining <- setdiff(remaining, best)
      dnew <- sqdist_to_point(centers, centers[best, ])
      dmin <- pmin(dmin, dnew)
    }
    groups[selected]
  }
  
  thr <- max(1L, floor(overlap_cutoff * k_eff))
  keep_flag <- determineOverlapCpp(nn_index, thr)
  keep_idx <- which(keep_flag == 0L)
  if (length(keep_idx) == 0L) keep_idx <- 1L
  nn_index_pruned <- nn_index[keep_idx, , drop = FALSE]
  
  cell_names <- rownames(emb)
  groups <- lapply(seq_len(nrow(nn_index_pruned)), function(i) unique(cell_names[nn_index_pruned[i, ]]))
  groups <- groups[vapply(groups, length, integer(1)) >= min_cells_per_metacell]
  
  if (isTRUE(force_target_metacells) && length(groups) < max_ct && nrow(nn_index) >= max_ct) {
    groups <- lapply(seq_len(nrow(nn_index)), function(i) unique(cell_names[nn_index[i, ]]))
    groups <- groups[vapply(groups, length, integer(1)) >= min_cells_per_metacell]
  }
  
  if (length(groups) > max_ct) {
    groups <- select_groups_from_emb(groups, max_ct, method = selection_method)
  }
  
  if (!isTRUE(allow_overlap)) {
    assigned <- setNames(rep(FALSE, n), cell_names)
    uniq_groups <- list()
    for (g in seq_along(groups)) {
      mem <- groups[[g]]
      new_mem <- mem[!assigned[mem]]
      if (length(new_mem) > 0L) {
        uniq_groups[[length(uniq_groups) + 1L]] <- new_mem
        assigned[new_mem] <- TRUE
      }
    }
    groups <- uniq_groups
  }
  
  names(groups) <- paste0(ct_name, "_", seq_along(groups))
  
  lens <- vapply(groups, length, integer(1))
  long_map <- data.frame(
    barcode = unlist(groups, use.names = FALSE),
    metacell = rep(names(groups), lens),
    celltype = ct_name,
    stringsAsFactors = FALSE
  )
  cell2metacell <- setNames(rep(NA_character_, n), cell_names)
  first_hit <- !duplicated(long_map$barcode)
  cell2metacell[long_map$barcode[first_hit]] <- long_map$metacell[first_hit]
  n_covered <- length(unique(long_map$barcode))
  reuse_counts <- as.integer(table(long_map$barcode))
  mean_reuse_covered <- if (length(reuse_counts) > 0L) mean(reuse_counts) else 0
  median_reuse_covered <- if (length(reuse_counts) > 0L) stats::median(reuse_counts) else 0
  q95_reuse_covered <- if (length(reuse_counts) > 0L) as.numeric(stats::quantile(reuse_counts, 0.95, names = FALSE)) else 0
  max_reuse <- if (length(reuse_counts) > 0L) max(reuse_counts) else 0L
  
  if (verbose) {
    message(ct_name, ": #groups=", length(groups),
            " n_cells=", n,
            " k_eff=", k_eff,
            " covered_cells=", n_covered,
            " coverage=", round(n_covered / n, 3),
            " mean_reuse_covered=", round(mean_reuse_covered, 2),
            " size_summary=", paste(capture.output(summary(lens)), collapse="; "))
  }
  
  list(
    metacell_members = groups,
    cell2metacell = cell2metacell,
    cell2metacell_long = long_map,
    seeds_used = rownames(emb)[seeds],
    k = k_eff,
    knn_iteration = n_seed,
    overlap_cutoff = overlap_cutoff,
    overlap_threshold = thr,
    target_metacells = target_ct,
    max_metacells = max_ct,
    min_metacells = min_metacells,
    target_cells_per_metacell = target_cells_per_metacell,
    min_cells_per_metacell = min_cells_per_metacell,
    balance_final_metacells = balance_final_metacells,
    force_target_metacells = force_target_metacells,
    target_memberships_per_cell = target_memberships_per_cell,
    selection_method = selection_method,
    allow_overlap = allow_overlap,
    n_cells = n,
    n_cells_covered = n_covered,
    cell_coverage_fraction = n_covered / n,
    mean_reuse_covered = mean_reuse_covered,
    median_reuse_covered = median_reuse_covered,
    q95_reuse_covered = q95_reuse_covered,
    max_reuse = max_reuse,
    mean_reuse_all = nrow(long_map) / n
  )
}

parse_peaks <- function(peak_names) {
  
  peak_names_std <- gsub("\\s+", "", as.character(peak_names))
  peak_names_std <- sub(
    "^(chr[^_:-]+)_([0-9]+)_([0-9]+)$",
    "\\1-\\2-\\3",
    peak_names_std
  )
  peak_df <- stringr::str_match(peak_names_std, "^(chr[^:-]+)[:\\-](\\d+)-(\\d+)$")
  
  valid_idx <- stats::complete.cases(peak_df)
  if (!any(valid_idx)) {
    stop("No valid peaks parsed. Expected format like 'chr1:100-200'.")
  }
  
  peak_df <- peak_df[valid_idx, , drop = FALSE]
  peak_names_valid <- peak_names[valid_idx]
  invalid_chr <- is.na(peak_df[, 2]) |
    peak_df[, 2] %in% c("chr", "chrNA", "chrNaN", "chrNULL") |
    grepl("^chr(NA|NaN|null)$", peak_df[, 2], ignore.case = TRUE)
  if (any(invalid_chr)) {
    stop(
      "Parsed peak IDs contain invalid chromosome labels such as chrNA. ",
      "Check upstream region_id construction."
    )
  }
  
  start_pos <- as.numeric(peak_df[, 3])
  end_pos   <- as.numeric(peak_df[, 4])
  
  keep <- !is.na(start_pos) & !is.na(end_pos) & end_pos >= start_pos
  if (!any(keep)) stop("All parsed peaks had invalid start/end.")
  
  peak_df <- peak_df[keep, , drop = FALSE]
  peak_names_valid <- peak_names_valid[keep]
  start_pos <- start_pos[keep]
  end_pos   <- end_pos[keep]
  
  peak_gr <- GRanges(
    seqnames = peak_df[, 2],
    ranges   = IRanges(start = start_pos, end = end_pos)
  )
  
  mcols(peak_gr)$peak_name <- peak_names_valid
  
  names(peak_gr) <- peak_names_valid
  
  peak_gr
}

make_gene_window_for_selected <- function(gr_anno,
                                          genes,
                                          flank = 250000,
                                          gene_col = "gene_name") {
  
  if (!inherits(gr_anno, "GRanges")) stop("gr_anno must be a GRanges.")
  if (!gene_col %in% colnames(mcols(gr_anno))) {
    stop(paste0("gene_col '", gene_col, "' not found in gr_anno mcols. ",
                "Available: ", paste(colnames(mcols(gr_anno)), collapse = ", ")))
  }
  
  g <- as.character(mcols(gr_anno)[[gene_col]])
  keep <- !is.na(g) & g %in% genes
  if (!any(keep)) stop("No annotation rows matched the selected genes.")
  
  gr_sub <- gr_anno[keep]
  g_sub  <- as.character(mcols(gr_sub)[[gene_col]])
  
  grl <- split(gr_sub, g_sub)
  red <- suppressWarnings(reduce(grl))
  
  spans1 <- endoapply(red, function(gr) {
    GRanges(
      seqnames = seqnames(gr)[1],
      ranges   = IRanges(start = min(start(gr)), end = max(end(gr))),
      strand   = "*"
    )
  })
  
  gene_span <- unlist(spans1, use.names = TRUE)
  mcols(gene_span)$gene <- names(gene_span)
  
  gene_window <- GRanges(
    seqnames = seqnames(gene_span),
    ranges   = IRanges(
      start = pmax(1, start(gene_span) - flank),
      end   = end(gene_span) + flank
    )
  )
  mcols(gene_window)$gene <- mcols(gene_span)$gene
  names(gene_window) <- mcols(gene_window)$gene
  
  gene_window
}

get_peaks_by_gene_gr <- function(peak_names, gene_window_gr) {
  
  peak_gr <- parse_peaks(peak_names)
  
  hits <- findOverlaps(peak_gr, gene_window_gr, ignore.strand = TRUE)
  if (length(hits) == 0) {
    return(setNames(list(), character(0)))
  }
  
  peaks_hit <- mcols(peak_gr)$peak_name[queryHits(hits)]
  genes_hit <- mcols(gene_window_gr)$gene[subjectHits(hits)]
  
  out <- split(peaks_hit, genes_hit)
  out <- lapply(out, unique)
  
  out
}

sourceCpp(file.path(.cliper_source_dir, "CER_PCGS.cpp"))
CER_PCGS <- function(X, y, K, n_iter = 5000, burn_in = 1000, alpha_conc = 5,
                     tau2 = 10000, rho0 = 1, p1 = 0.7, add_b_mu = FALSE, scale = TRUE,
                     scale_y = FALSE) {
  
  p <- ncol(X)
  n <- nrow(X)
  
  if (p < 3) stop("Need at least 3 predictors for CLIPER with null/effect clusters.")
  K <- min(K, p - 1)
  if (K < 2) stop("K must be at least 2.")
  
  if(p <= K){
    
    K <- p - 1
  }
  
  sds0 <- apply(X, 2, sd, na.rm = TRUE)
  keep0 <- is.finite(sds0) & sds0 > 0
  X <- X[, keep0, drop = FALSE]
  if (ncol(X) < 3) stop("Fewer than 3 nonzero-variance predictors after filtering.")
  p <- ncol(X)
  K <- min(K, p - 1)
  if (K < 2) stop("K must be at least 2 after filtering.")
  
  X <- scale(X, center = T, scale = scale)
  y <- scale(y, center = TRUE, scale = scale_y)
  
  selected_idx <- seq_len(p)
  p_cap <- ceiling(1.2 * n)
  if (p > p_cap) {
    sds <- apply(X, 2, sd)
    nonzero <- which(sds > 0 & !is.na(sds))
    
    cors <- suppressWarnings(stats::cor(X[, nonzero, drop = FALSE], y))
    if (all(is.na(cors))) stop("y has zero variance or correlations could not be computed.")
    
    ord <- order(abs(cors), decreasing = TRUE, na.last = NA)
    keep_m <- min(p_cap, length(ord))
    keep_in_nonzero <- ord[seq_len(keep_m)]
    selected_idx <- nonzero[keep_in_nonzero]
    
    X <- X[, selected_idx, drop = FALSE]
    p <- ncol(X)
    K <- min(K, p - 1)
  }
  
  ridge_cv <- cv.glmnet(X, y, alpha = 0, lambda = exp(seq(-5, 5, length = 100)))
  ridge_lambda <- ridge_cv$lambda.min
  ridge_fit <- glmnet(X, y, alpha = 0, lambda = ridge_lambda)
  beta_hat <- as.vector(coef(ridge_fit))[-1]
  
  n_beta_unique <- length(unique(beta_hat[is.finite(beta_hat)]))
  K <- min(K, p - 1, n_beta_unique)
  if (K < 2) stop("Fewer than 2 distinct ridge coefficients for kmeans initialization.")
  kmeans_result <- kmeans(beta_hat, centers = K, nstart = 10)
  
  centers <- as.numeric(kmeans_result$centers)
  cluster0 <- kmeans_result$cluster
  
  closest_cluster <- which.min(abs(centers))
  
  m <- matrix(0, nrow = p, ncol = K)
  m[cbind(seq_len(p), cluster0)] <- 1
  
  if (closest_cluster != 1) {
    swap <- c(1, closest_cluster)
    m[, swap] <- m[, rev(swap)]
    centers[swap] <- centers[rev(swap)]
  }
  
  for (k in 1:K) {
    if (sum(m[, k]) == 0) {
      donors <- which(colSums(m) > 1)
      if (!length(donors)) donors <- which(colSums(m) > 0)
      closest_non_empty <- donors[
        which.min(abs(centers[donors] - centers[k]))
      ]
      
      j_to_transfer <- which(m[, closest_non_empty] == 1)[1]
      m[j_to_transfer, ] <- 0
      m[j_to_transfer, k] <- 1
    }
  }
  
  res <- CER_PCGS_rcpp(
    X = X, y = y, K = K, n_iter = n_iter, burn_in = burn_in, alpha_conc = alpha_conc,
    tau2 = tau2, rho0 = rho0, p1 = p1, add_b_mu = add_b_mu,
    predict = FALSE, m_init = m
  )
  
  list(
    posterior_m      = res$posterior_m,
    posterior_b      = res$posterior_b,
    posterior_b_pred = res$posterior_b_pred,
    
    posterior_b_pred_sd   = res$posterior_b_pred_sd,
    posterior_b_pred_q025 = res$posterior_b_pred_q025,
    posterior_b_pred_q975 = res$posterior_b_pred_q975,
    
    posterior_y      = res$posterior_y,
    b_samples        = res$b_samples,
    input_peaks      = colnames(X),
    y_cliper         = y
  )
}

Run_CLIPER <- function(
    cliper_obj = NULL,
    flank = 500000,
    gene_list = NULL,
    gr_anno = NULL,
    p1 = 0.7, K = 5, n_iter = 10000, burn_in = 5000, alpha_conc = 5,
    tau2 = 5000, rho0 = 0.5,
    add_b_mu = FALSE, scale = TRUE, scale_y = FALSE, seed = 2001
){
  
  if (is.null(cliper_obj) || length(cliper_obj) == 0) stop("cliper_obj is NULL/empty.")
  if (is.null(gr_anno)) stop("Annotation is NULL. Add gene annotation first.")
  if (is.null(gene_list) || length(gene_list) == 0) stop("gene_list is NULL/empty.")
  if (!any(gene_list %in% gr_anno$gene_name)) stop("None of the genes in gene_list are in gr_anno$gene_name.")
  
  gene_window <- make_gene_window_for_selected(
    gr_anno  = gr_anno,
    genes    = gene_list,
    flank    = flank,
    gene_col = "gene_name"
  )
  
  cts <- names(cliper_obj)
  
  cliper_result <- setNames(
    lapply(seq_along(cts), function(ct_idx) {
      
      ct <- cts[ct_idx]
      message(sprintf(
        "[CLIPER] Cell type %d / %d : %s",
        ct_idx, length(cts), ct
      ))
      
      meta_atac <- t(cliper_obj[[ct]]$meta_atac)
      meta_rna  <- t(cliper_obj[[ct]]$meta_rna)
      
      peaks_by_gene <- get_peaks_by_gene_gr(
        peak_names     = colnames(meta_atac),
        gene_window_gr = gene_window
      )
      
      genes_to_test <- names(peaks_by_gene)
      n_gene <- length(genes_to_test)
      
      out_list <- list()
      all_list <- list()
      info_list <- list()
      out_i <- 0L
      all_i <- 0L
      info_i <- 0L
      
      for (g_idx in seq_along(genes_to_test)) {
        
        gene_select <- genes_to_test[g_idx]
        
        if (g_idx %% 10 == 0 || g_idx == 1 || g_idx == n_gene) {
          message(sprintf(
            "  └─ gene %d / %d (%s)",
            g_idx, n_gene, gene_select
          ))
        }
        
        if (!gene_select %in% colnames(meta_rna)) next
        
        peaks_near_gene <- peaks_by_gene[[gene_select]]
        if (is.null(peaks_near_gene) || length(peaks_near_gene) == 0) next
        
        meta_atac_select <- meta_atac[, peaks_near_gene, drop = FALSE]
        if (ncol(meta_atac_select) <= 2) next
        
        X <- meta_atac_select
        
        common_mc <- intersect(rownames(X), rownames(meta_rna))
        if (length(common_mc) < 3) next
        X_use <- X[common_mc, , drop = FALSE]
        y <- meta_rna[common_mc, gene_select, drop = TRUE]
        
        set.seed(seed + g_idx)
        result <- CER_PCGS(
          X = X_use, y = y,
          K = K,n_iter = n_iter, burn_in = burn_in,
          alpha_conc = alpha_conc,
          tau2 = tau2,rho0 = rho0, p1 = p1,
          add_b_mu = add_b_mu,
          scale = scale,
          scale_y = scale_y
        )
        
        posterior_m      <- result$posterior_m
        posterior_b      <- result$posterior_b
        posterior_b_pred <- result$posterior_b_pred
        input_peaks      <- result$input_peaks
        y_cliper         <- result$y_cliper
        posterior_b_pred_sd   <- result$posterior_b_pred_sd
        posterior_b_pred_q025 <- result$posterior_b_pred_q025
        posterior_b_pred_q975 <- result$posterior_b_pred_q975
        
        predicted_cluster <- apply(posterior_m, 1, which.max)
        
        all_i <- all_i + 1L
        all_list[[all_i]] <- data.frame(
          Peak = input_peaks,
          Gene = gene_select,
          Cell_Type = ct,
          Cluster = predicted_cluster,
          Posterior_b = posterior_b_pred,
          Posterior_b_sd = posterior_b_pred_sd,
          Beta_q025 = posterior_b_pred_q025,
          Beta_q975 = posterior_b_pred_q975,
          PPIP = 1 - posterior_m[, 1],
          stringsAsFactors = FALSE
        )
        
        idx_non_cluster1 <- which(predicted_cluster != 1)
        if (length(idx_non_cluster1) > 0) {
          out_i <- out_i + 1L
          out_list[[out_i]] <- data.frame(
            Peak = input_peaks[idx_non_cluster1],
            Gene = gene_select,
            Cell_Type = ct,
            Cluster = predicted_cluster[idx_non_cluster1],
            Posterior_b = posterior_b_pred[idx_non_cluster1],
            Posterior_b_sd = posterior_b_pred_sd[idx_non_cluster1],
            Beta_q025 = posterior_b_pred_q025[idx_non_cluster1],
            Beta_q975 = posterior_b_pred_q975[idx_non_cluster1],
            PPIP = 1 - posterior_m[idx_non_cluster1, 1],
            stringsAsFactors = FALSE
          )
        }
        
        info_i <- info_i + 1L
        info_list[[info_i]] <- data.frame(
          gene = gene_select,
          ct = ct,
          K_predicted = length(unique(predicted_cluster)),
          cluster = paste0(table(predicted_cluster), collapse = ","),
          stringsAsFactors = FALSE
        )
      }
      
      summary_all <- if (length(all_list)) dplyr::bind_rows(all_list) else data.frame()
      cliper_summary <- if (length(out_list)) dplyr::bind_rows(out_list) else data.frame()
      
      message(ct)
      
      list(
        summary_all = summary_all,
        cliper_summary = cliper_summary,
        summary_info   = if (length(info_list)) dplyr::bind_rows(info_list) else data.frame()
      )
    }),
    cts
  )
  
  cliper_result
}


compare_p2g <- function(
    cliper_1,
    cliper_2,
    name_1 = "condition1",
    name_2 = "condition2",
    by = c("Peak", "Gene", "Cell_Type"),
    pip_cutoff = 0.5,
    pip_rule = c("min", "max"),
    decision_rule = c("delta_ci", "beta_ci_overlap"),
    delta_cutoff = 0
) {
  pip_rule <- match.arg(pip_rule)
  decision_rule <- match.arg(decision_rule)
  
  required_cols <- c(by, "Posterior_b", "Beta_q025", "Beta_q975", "PPIP")
  
  missing_1 <- setdiff(required_cols, colnames(cliper_1))
  missing_2 <- setdiff(required_cols, colnames(cliper_2))
  
  if (length(missing_1) > 0) {
    stop("cliper_1 is missing columns: ", paste(missing_1, collapse = ", "))
  }
  if (length(missing_2) > 0) {
    stop("cliper_2 is missing columns: ", paste(missing_2, collapse = ", "))
  }
  
  x1 <- cliper_1[, required_cols, drop = FALSE]
  x2 <- cliper_2[, required_cols, drop = FALSE]
  
  colnames(x1)[match(c("Posterior_b", "Beta_q025", "Beta_q975", "PPIP"), colnames(x1))] <-
    paste0(c("Beta_mean", "Beta_q025", "Beta_q975", "PPIP"), "_", name_1)
  
  colnames(x2)[match(c("Posterior_b", "Beta_q025", "Beta_q975", "PPIP"), colnames(x2))] <-
    paste0(c("Beta_mean", "Beta_q025", "Beta_q975", "PPIP"), "_", name_2)
  
  out <- dplyr::inner_join(x1, x2, by = by)
  
  b1 <- out[[paste0("Beta_mean_", name_1)]]
  l1 <- out[[paste0("Beta_q025_", name_1)]]
  u1 <- out[[paste0("Beta_q975_", name_1)]]
  p1 <- out[[paste0("PPIP_", name_1)]]
  
  b2 <- out[[paste0("Beta_mean_", name_2)]]
  l2 <- out[[paste0("Beta_q025_", name_2)]]
  u2 <- out[[paste0("Beta_q975_", name_2)]]
  p2 <- out[[paste0("PPIP_", name_2)]]
  
  out$Delta_beta <- b1 - b2
  out$Abs_delta_beta <- abs(out$Delta_beta)
  
  out$Delta_PPIP <- p1 - p2
  out$Abs_delta_PPIP <- abs(out$Delta_PPIP)
  
  se1 <- (u1 - l1) / (2 * 1.96)
  se2 <- (u2 - l2) / (2 * 1.96)
  se_delta <- sqrt(se1^2 + se2^2)
  
  out$Beta_se_approx_1 <- se1
  out$Beta_se_approx_2 <- se2
  out$Delta_beta_se <- se_delta
  
  out$Delta_beta_q025 <- out$Delta_beta - 1.96 * se_delta
  out$Delta_beta_q975 <- out$Delta_beta + 1.96 * se_delta
  
  out$Delta_CI_cross0 <- out$Delta_beta_q025 <= 0 & out$Delta_beta_q975 >= 0
  out$Delta_CI_exclude0 <- !out$Delta_CI_cross0
  
  out$Beta_CI_overlap <- pmax(l1, l2) <= pmin(u1, u2)
  
  out$Direction <- ifelse(
    out$Delta_beta_q025 > 0,
    paste0(name_1, "_stronger"),
    ifelse(
      out$Delta_beta_q975 < 0,
      paste0(name_2, "_stronger"),
      "uncertain"
    )
  )
  
  out$Max_PPIP <- pmax(p1, p2)
  out$Min_PPIP <- pmin(p1, p2)
  
  out$Pass_PPIP <- if (pip_rule == "min") {
    out$Min_PPIP >= pip_cutoff
  } else {
    out$Max_PPIP >= pip_cutoff
  }
  
  out$Pass_delta <- out$Abs_delta_beta >= delta_cutoff
  
  out$Differential <- if (decision_rule == "delta_ci") {
    out$Pass_PPIP & out$Pass_delta & out$Delta_CI_exclude0
  } else {
    out$Pass_PPIP & out$Pass_delta & !out$Beta_CI_overlap
  }
  
  out <- out[order(
    out$Differential,
    out$Abs_delta_beta,
    out$Min_PPIP,
    decreasing = TRUE
  ), ]
  
  rownames(out) <- NULL
  
  out
}


# Robust region ID -> GRanges
# Accepts strings like "chr1-100-200", "chr1:100-200", "chr1_100_200", with spaces removed
ids_to_granges_safe <- function(region_ids) {
  region_ids <- as.character(region_ids)
  
  region_ids2 <- region_ids %>%
    stringr::str_replace_all(":", "-") %>%
    stringr::str_replace_all("_", "-") %>%
    stringr::str_replace_all("\\s+", "")
  
  gr <- Signac::StringToGRanges(region_ids2, sep = c("-", "-"))
  
  S4Vectors::mcols(gr)$region_id <- region_ids
  
  ok <- !is.na(as.character(GenomicRanges::seqnames(gr))) &
    !is.na(GenomicRanges::start(gr)) &
    !is.na(GenomicRanges::end(gr))
  
  if (any(!ok)) {
    message("Dropping ", sum(!ok),
            " region IDs that cannot be parsed to GRanges (check ID format).")
    gr <- gr[ok]
  }
  gr
}

# Generic feature intervals -> assay region IDs (label-aware)
# Example usage: map CRISPR intervals (label=gene) to ATAC peaks; but can be used for any feature intervals
map_feature_intervals_to_assay_regions <- function(
    obj,
    feature_df,
    label_col = "target_gene_short",
    chr_col  = "chr.candidate_enhancer",
    start_col = "start.candidate_enhancer",
    end_col   = "stop.candidate_enhancer",
    assay_name = "ATAC",
    min_overlap_bp = 1
) {
  assay_region_ids <- rownames(obj[[assay_name]])
  assay_gr <- ids_to_granges_safe(assay_region_ids)
  if (length(assay_gr) == 0) {
    stop("Assay region GRanges is empty after parsing. Check assay feature names.")
  }
  
  names(assay_gr) <- S4Vectors::mcols(assay_gr)$region_id
  
  feat <- feature_df %>%
    dplyr::transmute(
      label = .data[[label_col]],
      chr0  = as.character(.data[[chr_col]]),
      start = as.integer(.data[[start_col]]),
      end   = as.integer(.data[[end_col]])
    ) %>%
    dplyr::filter(!is.na(label), !is.na(chr0), !is.na(start), !is.na(end)) %>%
    dplyr::distinct()
  
  feat <- feat %>%
    dplyr::mutate(chr = ifelse(stringr::str_detect(chr0, "^chr"), chr0, paste0("chr", chr0))) %>%
    dplyr::select(-chr0)
  
  feature_gr <- GenomicRanges::GRanges(
    seqnames = feat$chr,
    ranges = IRanges::IRanges(start = feat$start, end = feat$end),
    label = feat$label
  )
  
  ok2 <- !is.na(as.character(GenomicRanges::seqnames(feature_gr)))
  if (any(!ok2)) {
    message("Dropping ", sum(!ok2), " feature intervals with NA seqnames.")
    feature_gr <- feature_gr[ok2]
  }
  
  hits <- GenomicRanges::findOverlaps(
    feature_gr, assay_gr,
    ignore.strand = TRUE,
    minoverlap = min_overlap_bp
  )
  
  if (length(hits) == 0) {
    return(tibble::tibble(label = character(), region_id = character()))
  }
  
  tibble::tibble(
    label = S4Vectors::mcols(feature_gr)$label[S4Vectors::queryHits(hits)],
    region_id = S4Vectors::mcols(assay_gr)$region_id[S4Vectors::subjectHits(hits)]
  ) %>%
    dplyr::distinct()
}

# Enrichment
# The cis region is constructed from gr_anno for the requested genes using
# gene body +/- cis_flank. For a fair comparison across methods, use the same
# genes, cis_flank, external_background, and external_targets in every call.

.enrich_gr_key <- function(gr) {
  paste0(
    as.character(GenomicRanges::seqnames(gr)), ":",
    BiocGenerics::start(gr), "-", BiocGenerics::end(gr)
  )
}

.enrich_as_gr <- function(x, label, deduplicate = TRUE) {
  if (inherits(x, "GRanges")) {
    gr <- x
  } else {
    if (!exists("parse_peaks", mode = "function")) {
      stop("parse_peaks() must be defined before using ", label, ".")
    }
    x <- unique(as.character(x))
    x <- x[!is.na(x) & nzchar(x)]
    if (length(x) == 0L) stop(label, " is empty.")
    gr <- parse_peaks(x)
  }
  
  if (length(gr) == 0L) stop(label, " cannot be parsed to GRanges.")
  
  if (isTRUE(deduplicate)) {
    gr <- gr[!duplicated(.enrich_gr_key(gr))]
  }
  gr
}

.enrich_as_named_gr <- function(x, label) {
  if (inherits(x, "GRanges")) {
    gr <- x
    
    if ("gene" %in% colnames(S4Vectors::mcols(gr))) {
      gene <- as.character(S4Vectors::mcols(gr)$gene)
    } else if (!is.null(names(gr)) && all(nzchar(names(gr)))) {
      gene <- as.character(names(gr))
    } else {
      stop(
        label,
        " must have a 'gene' metadata column or non-empty names."
      )
    }
    
    region_id <- .enrich_gr_key(gr)
  } else {
    if (
      is.null(names(x)) ||
      any(is.na(names(x))) ||
      any(names(x) == "")
    ) {
      stop(label, " must be a named character vector with names = gene.")
    }
    
    gene <- as.character(names(x))
    region_id <- as.character(x)
    keep <- !is.na(gene) & nzchar(gene) & !is.na(region_id) & nzchar(region_id)
    gene <- gene[keep]
    region_id <- region_id[keep]
    if (length(region_id) == 0L) stop(label, " is empty after filtering.")
    
    region_std <- gsub("\\s+", "", region_id)
    region_std <- gsub("[:_]", "-", region_std)
    mat <- stringr::str_match(
      region_std,
      "^(chr[^-]+)-([0-9]+)-([0-9]+)$"
    )
    ok <- !is.na(mat[, 1])
    
    if (any(!ok)) {
      message(
        "Dropping ", sum(!ok), " ", label,
        " region IDs that cannot be parsed."
      )
    }
    
    gene <- gene[ok]
    region_id <- region_id[ok]
    mat <- mat[ok, , drop = FALSE]
    if (length(region_id) == 0L) {
      stop(label, " contains no valid genomic regions.")
    }
    
    start_pos <- as.integer(mat[, 3])
    end_pos <- as.integer(mat[, 4])
    valid_range <- !is.na(start_pos) & !is.na(end_pos) & end_pos >= start_pos
    
    gene <- gene[valid_range]
    region_id <- region_id[valid_range]
    mat <- mat[valid_range, , drop = FALSE]
    start_pos <- start_pos[valid_range]
    end_pos <- end_pos[valid_range]
    
    if (length(region_id) == 0L) {
      stop(label, " contains no valid genomic ranges.")
    }
    
    gr <- GenomicRanges::GRanges(
      seqnames = mat[, 2],
      ranges = IRanges::IRanges(start = start_pos, end = end_pos)
    )
  }
  
  region_key <- .enrich_gr_key(gr)
  pair_id <- paste(gene, region_key, sep = "||")
  
  S4Vectors::mcols(gr)$region_id <- region_id
  S4Vectors::mcols(gr)$region_key <- region_key
  S4Vectors::mcols(gr)$gene <- gene
  S4Vectors::mcols(gr)$pair_id <- pair_id
  
  gr[!duplicated(pair_id)]
}

.enrich_make_cis_windows <- function(
    gr_anno,
    genes,
    cis_flank,
    gene_col = "gene_name"
) {
  if (!exists("make_gene_window_for_selected", mode = "function")) {
    stop(
      "make_gene_window_for_selected() must be defined before running ",
      "the enrichment functions."
    )
  }
  if (!inherits(gr_anno, "GRanges")) stop("gr_anno must be a GRanges.")
  if (is.null(genes) || length(genes) == 0L) stop("genes cannot be NULL/empty.")
  if (!is.finite(cis_flank) || length(cis_flank) != 1L || cis_flank < 0) {
    stop("cis_flank must be one non-negative finite number.")
  }
  
  genes <- unique(as.character(genes))
  genes <- genes[!is.na(genes) & nzchar(genes)]
  if (length(genes) == 0L) stop("No valid genes remain.")
  
  cis_by_gene <- make_gene_window_for_selected(
    gr_anno = gr_anno,
    genes = genes,
    flank = as.integer(cis_flank),
    gene_col = gene_col
  )
  
  if (length(cis_by_gene) == 0L) stop("No cis regions were constructed.")
  
  if (!"gene" %in% colnames(S4Vectors::mcols(cis_by_gene))) {
    if (!is.null(names(cis_by_gene)) && all(nzchar(names(cis_by_gene)))) {
      S4Vectors::mcols(cis_by_gene)$gene <- names(cis_by_gene)
    } else {
      stop("Constructed cis regions do not contain gene labels.")
    }
  }
  
  cis_by_gene
}

.enrich_overlap_flags <- function(query_gr, subject_gr, min_overlap_bp = 1L) {
  out <- rep(FALSE, length(query_gr))
  if (length(query_gr) == 0L || length(subject_gr) == 0L) return(out)
  
  hits <- GenomicRanges::findOverlaps(
    query_gr,
    subject_gr,
    minoverlap = as.integer(min_overlap_bp),
    ignore.strand = TRUE
  )
  if (length(hits) > 0L) {
    out[unique(S4Vectors::queryHits(hits))] <- TRUE
  }
  out
}

.enrich_pair_overlap_flags <- function(query_gr, subject_gr, min_overlap_bp = 1L) {
  out <- rep(FALSE, length(query_gr))
  if (length(query_gr) == 0L || length(subject_gr) == 0L) return(out)
  
  hits <- GenomicRanges::findOverlaps(
    query_gr,
    subject_gr,
    minoverlap = as.integer(min_overlap_bp),
    ignore.strand = TRUE
  )
  if (length(hits) == 0L) return(out)
  
  q_gene <- as.character(
    S4Vectors::mcols(query_gr)$gene[S4Vectors::queryHits(hits)]
  )
  s_gene <- as.character(
    S4Vectors::mcols(subject_gr)$gene[S4Vectors::subjectHits(hits)]
  )
  hits <- hits[q_gene == s_gene]
  
  if (length(hits) > 0L) {
    out[unique(S4Vectors::queryHits(hits))] <- TRUE
  }
  out
}

.enrich_target_flags <- function(
    query_gr,
    target_gr,
    match = c("exact", "overlap"),
    min_overlap_bp = 1L,
    gene_matched = FALSE
) {
  match <- match.arg(match)
  
  if (match == "exact") {
    if (isTRUE(gene_matched)) {
      return(
        S4Vectors::mcols(query_gr)$pair_id %in%
          S4Vectors::mcols(target_gr)$pair_id
      )
    }
    return(.enrich_gr_key(query_gr) %in% .enrich_gr_key(target_gr))
  }
  
  if (isTRUE(gene_matched)) {
    .enrich_pair_overlap_flags(query_gr, target_gr, min_overlap_bp)
  } else {
    .enrich_overlap_flags(query_gr, target_gr, min_overlap_bp)
  }
}

.enrich_safe_ratio <- function(num, den) {
  if (
    length(num) != 1L || length(den) != 1L ||
    !is.finite(num) || !is.finite(den) || den <= 0
  ) {
    return(NA_real_)
  }
  num / den
}

.enrich_safe_enrichment <- function(
    G_selected,
    K_selected,
    G_background,
    K_background
) {
  p_selected <- .enrich_safe_ratio(G_selected, K_selected)
  p_background <- .enrich_safe_ratio(G_background, K_background)
  
  if (is.finite(p_selected) && is.finite(p_background) && p_background > 0) {
    p_selected / p_background
  } else {
    NA_real_
  }
}

summarize_trait_bootstrap_enrichment <- function(
    dt,
    group_cols,
    trait_col = "Phenotype",
    value_col = "enrichment",
    n_boot = 1000L,
    conf_level = 0.95,
    seed = 1L,
    min_traits = 2L,
    status_col = "status",
    ok_status = "ok"
) {
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Package data.table is required.")
  }
  if (is.null(dt) || nrow(dt) == 0L) {
    return(data.table::data.table())
  }
  x <- data.table::as.data.table(dt)
  required_cols <- unique(c(group_cols, trait_col, value_col))
  missing_cols <- setdiff(required_cols, colnames(x))
  if (length(missing_cols) > 0L) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  if (!is.null(status_col) && status_col %in% colnames(x)) {
    x <- x[get(status_col) == ok_status]
  }
  x[, (value_col) := suppressWarnings(as.numeric(get(value_col)))]
  x <- x[!is.na(get(trait_col)) & nzchar(as.character(get(trait_col))) &
           is.finite(get(value_col))]
  if (nrow(x) == 0L) return(data.table::data.table())

  by_cols <- unique(c(group_cols, trait_col))
  trait_dt <- x[
    ,
    .(
      trait_enrichment = mean(get(value_col), na.rm = TRUE),
      n_rows_per_trait = .N
    ),
    by = by_cols
  ]

  alpha <- (1 - conf_level) / 2
  old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  } else {
    NULL
  }
  on.exit({
    if (is.null(old_seed)) {
      if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    } else {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(seed)

  trait_dt[
    ,
    {
      values <- trait_enrichment[is.finite(trait_enrichment)]
      n_traits <- length(values)
      if (n_traits >= min_traits && n_boot > 0L) {
        boot <- replicate(
          as.integer(n_boot),
          mean(sample(values, size = n_traits, replace = TRUE)),
          simplify = TRUE
        )
        ci_low <- unname(stats::quantile(boot, probs = alpha, na.rm = TRUE, names = FALSE))
        ci_high <- unname(stats::quantile(boot, probs = 1 - alpha, na.rm = TRUE, names = FALSE))
      } else {
        ci_low <- NA_real_
        ci_high <- NA_real_
      }
      .(
        n_traits = n_traits,
        n_input_rows = sum(n_rows_per_trait),
        mean_enrichment = mean(values, na.rm = TRUE),
        median_enrichment = stats::median(values, na.rm = TRUE),
        sd_enrichment = stats::sd(values, na.rm = TRUE),
        se_enrichment = stats::sd(values, na.rm = TRUE) / sqrt(n_traits),
        bootstrap_ci_low = ci_low,
        bootstrap_ci_high = ci_high,
        bootstrap_n = as.integer(n_boot),
        bootstrap_conf_level = conf_level,
        bootstrap_unit = trait_col,
        bootstrap_weighting = "equal_trait"
      )
    },
    by = group_cols
  ]
}

.enrich_fisher_selected_vs_rest <- function(
    G_selected,
    K_selected,
    G_rest,
    K_rest,
    alternative = "greater",
    conf_level = 0.95,
    target_label = "external_target"
) {
  mat <- matrix(
    c(
      G_selected, K_selected - G_selected,
      G_rest, K_rest - G_rest
    ),
    nrow = 2L,
    byrow = TRUE,
    dimnames = list(
      c("selected", "cis_rest"),
      c(target_label, paste0("non_", target_label))
    )
  )
  
  out <- list(
    odds_ratio = NA_real_,
    or_ci_lower = NA_real_,
    or_ci_upper = NA_real_,
    p_value = NA_real_,
    alternative = alternative,
    conf_level = conf_level,
    table = mat
  )
  
  if (
    K_selected > 0L && K_rest > 0L &&
    G_selected >= 0L && G_selected <= K_selected &&
    G_rest >= 0L && G_rest <= K_rest
  ) {
    ft <- suppressWarnings(
      stats::fisher.test(
        mat,
        alternative = alternative,
        conf.level = conf_level
      )
    )
    out$odds_ratio <- unname(ft$estimate)
    out$or_ci_lower <- unname(ft$conf.int[1])
    out$or_ci_upper <- unname(ft$conf.int[2])
    out$p_value <- ft$p.value
  }
  
  out
}


## Variant-level enrichment without target-gene matching
##    Examples: GWAS variants
# -----------------------------------------------------------------------------
compute_region_enrichment_cisRest_externalBG <- function(
    selected_regions,
    external_background,
    external_targets,
    gr_anno,
    genes = NULL,
    cis_flank = 500000,
    gene_col = "gene_name",
    min_overlap_bp = 1,
    restrict_selected_to_cis = TRUE,
    target_match = c("exact", "overlap"),
    fisher_alternative = c("two.sided", "greater", "less"),
    conf_level = 0.95,
    verbose = TRUE,
    return_variant_tables = FALSE,
    return_cis_regions = FALSE
) {
  target_match <- match.arg(target_match)
  fisher_alternative <- match.arg(fisher_alternative)
  
  if (length(selected_regions) == 0L) stop("selected_regions is empty.")
  if (length(external_background) == 0L) stop("external_background is empty.")
  if (length(external_targets) == 0L) stop("external_targets is empty.")
  if (!is.finite(min_overlap_bp) || min_overlap_bp < 1L) {
    stop("min_overlap_bp must be at least 1.")
  }
  if (!is.finite(conf_level) || conf_level <= 0 || conf_level >= 1) {
    stop("conf_level must be between 0 and 1.")
  }
  
  if (is.null(genes)) {
    selected_names <- names(selected_regions)
    selected_names <- selected_names[!is.na(selected_names) & nzchar(selected_names)]
    if (length(selected_names) == 0L) {
      stop(
        "Provide genes explicitly. They are required to construct the cis ",
        "background when selected_regions is not named by gene."
      )
    }
    genes <- unique(selected_names)
    warning(
      "genes was inferred from names(selected_regions). For comparisons across ",
      "methods, explicitly provide the same genes in every call."
    )
  }
  
  cis_by_gene <- .enrich_make_cis_windows(
    gr_anno = gr_anno,
    genes = genes,
    cis_flank = cis_flank,
    gene_col = gene_col
  )
  cis_union <- GenomicRanges::reduce(cis_by_gene, ignore.strand = TRUE)
  
  SEL0 <- .enrich_as_gr(selected_regions, "selected_regions")
  EXT_BG <- .enrich_as_gr(external_background, "external_background")
  EXT_TG <- .enrich_as_gr(external_targets, "external_targets")
  
  selected_in_cis <- .enrich_overlap_flags(SEL0, cis_union, min_overlap_bp)
  SEL <- if (isTRUE(restrict_selected_to_cis)) {
    SEL0[selected_in_cis]
  } else {
    SEL0
  }
  if (length(SEL) == 0L) {
    stop("No selected regions overlap the constructed cis regions.")
  }
  
  # Eligible variants are restricted to the union of the requested cis regions.
  # A variant is counted once even if several gene windows overlap it.
  in_cis <- .enrich_overlap_flags(EXT_BG, cis_union, min_overlap_bp)
  CIS_VARIANTS <- EXT_BG[in_cis]
  if (length(CIS_VARIANTS) == 0L) {
    stop("No external_background variants overlap the constructed cis regions.")
  }
  
  is_selected <- .enrich_overlap_flags(CIS_VARIANTS, SEL, min_overlap_bp)
  is_target <- .enrich_target_flags(
    query_gr = CIS_VARIANTS,
    target_gr = EXT_TG,
    match = target_match,
    min_overlap_bp = min_overlap_bp,
    gene_matched = FALSE
  )
  is_rest <- !is_selected
  
  K_selected <- sum(is_selected)
  G_selected <- sum(is_selected & is_target)
  K_background <- length(CIS_VARIANTS)
  G_background <- sum(is_target)
  K_rest <- sum(is_rest)
  G_rest <- sum(is_rest & is_target)
  
  density_selected <- .enrich_safe_ratio(G_selected, K_selected)
  density_background <- .enrich_safe_ratio(G_background, K_background)
  density_rest <- .enrich_safe_ratio(G_rest, K_rest)
  enrichment <- .enrich_safe_enrichment(
    G_selected,
    K_selected,
    G_background,
    K_background
  )
  
  fisher <- .enrich_fisher_selected_vs_rest(
    G_selected = G_selected,
    K_selected = K_selected,
    G_rest = G_rest,
    K_rest = K_rest,
    alternative = fisher_alternative,
    conf_level = conf_level,
    target_label = "external_target"
  )
  fisher$inference_note <- paste(
    "Fisher's exact test compares selected variants with the non-overlapping",
    "cis-rest variants. Its odds ratio is distinct from the main enrichment",
    "estimate, whose denominator is the full cis region."
  )
  
  if (K_selected == 0L) {
    warning("No eligible external-background variants overlap selected_regions.")
  }
  if (G_background == 0L) {
    warning("No external targets were found in the full cis background.")
  }
  
  if (verbose) {
    message("Variant enrichment using the full cis-region background")
    message(
      "Genes: ", length(unique(as.character(genes))),
      " | cis flank: ", format(cis_flank, big.mark = ","), " bp"
    )
    message(
      "Selected regions in cis: ", length(SEL),
      " | eligible cis variants: ", K_background
    )
    message(
      "Selected density: ", G_selected, "/", K_selected,
      " = ", signif(density_selected, 4)
    )
    message(
      "Full cis density: ", G_background, "/", K_background,
      " = ", signif(density_background, 4)
    )
    message("Enrichment: ", signif(enrichment, 4))
    message(
      "Fisher selected-vs-cis-rest OR: ", signif(fisher$odds_ratio, 4),
      " | p = ", signif(fisher$p_value, 4)
    )
  }
  
  out <- list(
    counts = list(
      genes_requested = length(unique(as.character(genes))),
      genes_with_cis_region = length(unique(S4Vectors::mcols(cis_by_gene)$gene)),
      selected_region_n = length(SEL),
      cis_region_n_before_union = length(cis_by_gene),
      cis_region_n_after_union = length(cis_union),
      K_selected = K_selected,
      G_selected = G_selected,
      K_background = K_background,
      G_background = G_background,
      K_cis_rest = K_rest,
      G_cis_rest = G_rest
    ),
    density = list(
      selected = density_selected,
      cis_background = density_background,
      cis_rest = density_rest
    ),
    enrichment = enrichment,
    risk_ratio = list(
      estimate = enrichment,
      definition = paste(
        "target-variant density in selected regions /",
        "target-variant density in the full cis-region background"
      )
    ),
    fisher = fisher,
    settings = list(
      cis_flank = cis_flank,
      gene_col = gene_col,
      selected_included_in_background_denominator = TRUE,
      unit = "unique external-background variant",
      target_match = target_match,
      min_overlap_bp = min_overlap_bp
    )
  )
  
  if (isTRUE(return_variant_tables)) {
    out$variant_table <- data.frame(
      variant_id = .enrich_gr_key(CIS_VARIANTS),
      in_selected = is_selected,
      in_cis_rest = is_rest,
      external_target = is_target,
      stringsAsFactors = FALSE
    )
  }
  if (isTRUE(return_cis_regions)) {
    out$cis_regions_by_gene <- cis_by_gene
    out$cis_regions_union <- cis_union
  }
  
  out
}


## Region-gene-pair enrichment with target-gene matching
##    Examples: eQTL, CRISPRi, HiChIP, promoter-capture Hi-C
# -----------------------------------------------------------------------------
compute_region_enrichment_cisRest_externalBG_geneMatched <- function(
    selected_regions,
    external_background,
    external_targets,
    gr_anno,
    genes = NULL,
    cis_flank = 500000,
    gene_col = "gene_name",
    min_overlap_bp = 1,
    restrict_selected_to_cis = TRUE,
    external_target_match = c("exact", "overlap"),
    fisher_alternative = c("two.sided", "greater", "less"),
    conf_level = 0.95,
    n_boot = 1000,
    bootstrap_seed = NULL,
    verbose = TRUE,
    return_pair_tables = FALSE,
    return_cis_regions = FALSE
) {
  external_target_match <- match.arg(external_target_match)
  fisher_alternative <- match.arg(fisher_alternative)
  
  if (length(selected_regions) == 0L) stop("selected_regions is empty.")
  if (length(external_background) == 0L) stop("external_background is empty.")
  if (length(external_targets) == 0L) stop("external_targets is empty.")
  if (!is.finite(min_overlap_bp) || min_overlap_bp < 1L) {
    stop("min_overlap_bp must be at least 1.")
  }
  if (!is.finite(conf_level) || conf_level <= 0 || conf_level >= 1) {
    stop("conf_level must be between 0 and 1.")
  }
  if (!is.finite(n_boot) || n_boot < 0) {
    stop("n_boot must be non-negative.")
  }
  n_boot <- as.integer(n_boot)
  
  SEL0   <- .enrich_as_named_gr(selected_regions, "selected_regions")
  EXT_BG <- .enrich_as_named_gr(external_background, "external_background")
  EXT_TG <- .enrich_as_named_gr(external_targets, "external_targets")
  
  if (is.null(genes)) {
    genes <- unique(as.character(S4Vectors::mcols(EXT_BG)$gene))
    warning(
      "genes was inferred from external_background. For comparisons across ",
      "methods, explicitly provide the same genes in every call."
    )
  }
  genes <- unique(as.character(genes))
  
  cis_by_gene <- .enrich_make_cis_windows(
    gr_anno = gr_anno,
    genes = genes,
    cis_flank = cis_flank,
    gene_col = gene_col
  )
  
  # Add the metadata required by gene-matched overlap functions.
  cis_gene <- as.character(S4Vectors::mcols(cis_by_gene)$gene)
  cis_region_key <- .enrich_gr_key(cis_by_gene)
  S4Vectors::mcols(cis_by_gene)$region_id <- cis_region_key
  S4Vectors::mcols(cis_by_gene)$region_key <- cis_region_key
  S4Vectors::mcols(cis_by_gene)$pair_id <- paste(
    cis_gene,
    cis_region_key,
    sep = "||"
  )
  
  selected_in_cis <- .enrich_pair_overlap_flags(
    SEL0,
    cis_by_gene,
    min_overlap_bp
  )
  SEL <- if (isTRUE(restrict_selected_to_cis)) {
    SEL0[selected_in_cis]
  } else {
    SEL0
  }
  if (length(SEL) == 0L) {
    stop("No selected region-gene pairs overlap their corresponding cis regions.")
  }
  
  # Analysis units: eligible external region-gene pairs in the gene-specific cis regions.
  in_cis <- .enrich_pair_overlap_flags(
    EXT_BG,
    cis_by_gene,
    min_overlap_bp
  )
  UNITS <- EXT_BG[in_cis]
  if (length(UNITS) == 0L) {
    stop(
      "No external_background region-gene pairs overlap their corresponding ",
      "cis regions."
    )
  }
  
  # A unit is selected only if it overlaps a selected region linked to the same gene.
  selected_flag <- .enrich_pair_overlap_flags(
    UNITS,
    SEL,
    min_overlap_bp
  )
  
  # A unit is a target only if it overlaps a target region linked to the same gene.
  target_flag <- .enrich_target_flags(
    query_gr = UNITS,
    target_gr = EXT_TG,
    match = external_target_match,
    min_overlap_bp = min_overlap_bp,
    gene_matched = TRUE
  )
  
  background_flag <- rep(TRUE, length(UNITS))
  rest_flag <- !selected_flag
  
  genes_in_units <- sort(unique(as.character(S4Vectors::mcols(UNITS)$gene)))
  gene_rows <- lapply(genes_in_units, function(g) {
    idx <- which(as.character(S4Vectors::mcols(UNITS)$gene) == g)
    sel_idx <- idx[selected_flag[idx]]
    rest_idx <- idx[rest_flag[idx]]
    
    K_selected <- length(sel_idx)
    G_selected <- sum(target_flag[sel_idx])
    K_background <- length(idx)
    G_background <- sum(target_flag[idx])
    K_rest <- length(rest_idx)
    G_rest <- sum(target_flag[rest_idx])
    
    density_selected <- .enrich_safe_ratio(G_selected, K_selected)
    density_background <- .enrich_safe_ratio(G_background, K_background)
    density_rest <- .enrich_safe_ratio(G_rest, K_rest)
    enrichment_gene <- .enrich_safe_enrichment(
      G_selected,
      K_selected,
      G_background,
      K_background
    )
    
    data.frame(
      gene = g,
      K_selected = K_selected,
      G_selected = G_selected,
      K_background = K_background,
      G_background = G_background,
      K_cis_rest = K_rest,
      G_cis_rest = G_rest,
      density_selected = density_selected,
      density_background = density_background,
      density_cis_rest = density_rest,
      enrichment_gene = enrichment_gene,
      has_selected_background_pair = K_selected > 0,
      has_background_target = G_background > 0,
      valid_gene_ratio = is.finite(enrichment_gene),
      stringsAsFactors = FALSE
    )
  })
  gene_table <- do.call(rbind, gene_rows)
  
  K_selected <- sum(selected_flag)
  G_selected <- sum(selected_flag & target_flag)
  K_background <- sum(background_flag)
  G_background <- sum(background_flag & target_flag)
  K_rest <- sum(rest_flag)
  G_rest <- sum(rest_flag & target_flag)
  
  density_selected <- .enrich_safe_ratio(G_selected, K_selected)
  density_background <- .enrich_safe_ratio(G_background, K_background)
  density_rest <- .enrich_safe_ratio(G_rest, K_rest)
  
  enrichment <- .enrich_safe_enrichment(
    G_selected,
    K_selected,
    G_background,
    K_background
  )
  
  enrichment_definition <- paste(
    "pooled selected density / pooled full cis density,",
    "using gene-matched external-background region-gene pairs"
  )
  
  fisher <- .enrich_fisher_selected_vs_rest(
    G_selected = G_selected,
    K_selected = K_selected,
    G_rest = G_rest,
    K_rest = K_rest,
    alternative = fisher_alternative,
    conf_level = conf_level,
    target_label = "external_target_gene_matched"
  )
  fisher$inference_note <- paste(
    "Fisher's exact test compares selected pairs with the non-overlapping",
    "cis-rest pairs. Its odds ratio is distinct from the main enrichment",
    "estimate, whose denominator contains all eligible cis pairs."
  )
  
  bootstrap_values <- numeric(0)
  bootstrap_ci_lower <- NA_real_
  bootstrap_ci_upper <- NA_real_
  bootstrap_n_valid <- 0L
  
  if (n_boot > 0L && nrow(gene_table) >= 2L) {
    if (!is.null(bootstrap_seed)) set.seed(bootstrap_seed)
    
    bootstrap_values <- replicate(n_boot, {
      sampled <- sample(
        seq_len(nrow(gene_table)),
        size = nrow(gene_table),
        replace = TRUE
      )
      .enrich_safe_enrichment(
        sum(gene_table$G_selected[sampled]),
        sum(gene_table$K_selected[sampled]),
        sum(gene_table$G_background[sampled]),
        sum(gene_table$K_background[sampled])
      )
    })
    
    bootstrap_values <- bootstrap_values[is.finite(bootstrap_values)]
    bootstrap_n_valid <- length(bootstrap_values)
    
    if (bootstrap_n_valid > 0L) {
      alpha <- (1 - conf_level) / 2
      bootstrap_ci_lower <- unname(
        stats::quantile(bootstrap_values, probs = alpha, na.rm = TRUE)
      )
      bootstrap_ci_upper <- unname(
        stats::quantile(bootstrap_values, probs = 1 - alpha, na.rm = TRUE)
      )
    }
  }
  
  if (K_selected == 0L) {
    warning("No eligible external-background pairs overlap selected_regions.")
  }
  if (G_background == 0L) {
    warning("No external targets were found in the full cis background.")
  }
  
  if (verbose) {
    message("Gene-matched enrichment")
    message(
      "Genes: ", length(unique(as.character(genes))),
      " | cis flank: ", format(cis_flank, big.mark = ","), " bp"
    )
    message(
      "Selected density: ", G_selected, "/", K_selected,
      " = ", signif(density_selected, 4)
    )
    message(
      "Full cis density: ", G_background, "/", K_background,
      " = ", signif(density_background, 4)
    )
    message(
      "Enrichment = Selected density / Full cis density = ",
      signif(enrichment, 4)
    )
    message(
      "Fisher OR (selected vs cis-rest): ",
      signif(fisher$odds_ratio, 4),
      " | ", conf_level * 100, "% CI [",
      signif(fisher$or_ci_lower, 4), ", ",
      signif(fisher$or_ci_upper, 4), "]",
      " | p = ", signif(fisher$p_value, 4),
      " | alternative = ", fisher_alternative
    )
    if (bootstrap_n_valid > 0L) {
      message(
        conf_level * 100, "% gene-bootstrap CI for pooled enrichment [",
        signif(bootstrap_ci_lower, 4), ", ",
        signif(bootstrap_ci_upper, 4), "]"
      )
    }
  }
  
  out <- list(
    counts = list(
      genes_requested = length(unique(as.character(genes))),
      genes_with_cis_region = length(unique(S4Vectors::mcols(cis_by_gene)$gene)),
      selected_region_gene_n = length(SEL),
      external_background_pair_n_in_cis = length(UNITS),
      K_selected = K_selected,
      G_selected = G_selected,
      K_background = K_background,
      G_background = G_background,
      K_cis_rest = K_rest,
      G_cis_rest = G_rest,
      genes_total = nrow(gene_table),
      genes_with_selected_external_background_pair = sum(gene_table$has_selected_background_pair),
      genes_with_background_target = sum(gene_table$has_background_target),
      genes_with_finite_gene_ratio = sum(gene_table$valid_gene_ratio)
    ),
    density = list(
      selected = density_selected,
      full_cis = density_background,
      cis_rest = density_rest,
      selected_pooled = density_selected,
      cis_background_pooled = density_background,
      cis_rest_pooled = density_rest
    ),
    enrichment = enrichment,
    risk_ratio = list(
      estimate = enrichment,
      definition = enrichment_definition,
      bootstrap_ci_lower = bootstrap_ci_lower,
      bootstrap_ci_upper = bootstrap_ci_upper,
      bootstrap_n_valid = bootstrap_n_valid,
      conf_level = conf_level
    ),
    fisher = fisher,
    gene_table = gene_table,
    settings = list(
      cis_flank = cis_flank,
      gene_col = gene_col,
      selected_included_in_background_denominator = TRUE,
      unit = "external-background region-gene pair",
      aggregation = "pooled",
      target_match = external_target_match,
      min_overlap_bp = min_overlap_bp
    )
  )
  
  if (isTRUE(return_pair_tables)) {
    out$analysis_unit_table <- data.frame(
      gene = as.character(S4Vectors::mcols(UNITS)$gene),
      region_id = as.character(S4Vectors::mcols(UNITS)$region_id),
      pair_id = as.character(S4Vectors::mcols(UNITS)$pair_id),
      in_selected = selected_flag,
      in_cis_rest = rest_flag,
      external_target = target_flag,
      stringsAsFactors = FALSE
    )
  }
  if (isTRUE(return_cis_regions)) {
    out$cis_regions_by_gene <- cis_by_gene
  }
  
  out
}

liftover_snp_hg19_to_hg38 <- function(
    df,
    chr_col,
    pos_col,
    chain_path = CHAIN_FILE,
    seqlevels_style = "UCSC",
    keep_standard = TRUE,
    verbose = TRUE
) {
  stopifnot(is.data.frame(df))
  for (nm in c(chr_col, pos_col)) {
    if (!nm %in% colnames(df)) stop("Column not found: ", nm)
  }
  if (!file.exists(chain_path)) stop("chain_path not found: ", chain_path)
  
  dat <- df %>%
    mutate(
      .row_id = dplyr::row_number(),
      .chr = as.character(.data[[chr_col]]),
      .pos = suppressWarnings(as.integer(.data[[pos_col]]))
    ) %>%
    filter(!is.na(.chr), !is.na(.pos), .pos > 0)
  
  if (nrow(dat) == 0L) stop("No valid SNP positions after cleaning GWAS coordinates.")
  
  gr37 <- makeGRangesFromDataFrame(
    dat,
    seqnames.field = ".chr",
    start.field = ".pos",
    end.field = ".pos",
    keep.extra.columns = TRUE
  )
  seqlevelsStyle(gr37) <- seqlevels_style
  
  if (keep_standard) {
    gr37 <- keepStandardChromosomes(gr37, pruning.mode = "coarse")
  }
  
  chain <- import.chain(chain_path)
  gr38_list <- liftOver(gr37, chain)
  n_map <- elementNROWS(gr38_list)
  idx_unique <- which(n_map == 1L)
  
  gr38_unique <- unlist(gr38_list[idx_unique], use.names = FALSE)
  genome(gr38_unique) <- "hg38"
  
  row_id <- mcols(gr37)$.row_id
  df_hg38 <- as.data.table(df[row_id[idx_unique], , drop = FALSE])
  
  original_chr_col <- paste0(chr_col, "_GRCh37")
  original_pos_col <- paste0(pos_col, "_GRCh37")
  df_hg38[, (original_chr_col) := get(chr_col)]
  df_hg38[, (original_pos_col) := get(pos_col)]
  df_hg38[, (chr_col) := as.character(seqnames(gr38_unique))]
  df_hg38[, (pos_col) := start(gr38_unique)]
  
  diagnostics <- data.table(
    row_id = row_id,
    chromosome_GRCh37 = as.character(seqnames(gr37)),
    base_pair_location_GRCh37 = start(gr37),
    n_hg38_maps = n_map,
    liftover_status = fifelse(
      n_map == 1L,
      "unique",
      fifelse(n_map == 0L, "unmapped", "multi_mapped")
    )
  )
  
  if (verbose) {
    message(
      paste0(
        "Input rows: ", nrow(df), "\n",
        "Valid SNPs used: ", length(gr37), "\n",
        "Unique hg38 maps: ", length(idx_unique), "\n",
        "Multi-mapped SNPs: ", sum(n_map > 1L), "\n",
        "Unmapped SNPs: ", sum(n_map == 0L)
      )
    )
  }
  
  list(df_hg38 = df_hg38, diagnostics = diagnostics)
}

bind_cliper_summary <- function(
    cliper_output,
    use = c("summary_all", "cliper_summary")
) {
  use <- match.arg(use)
  
  if (is.data.frame(cliper_output)) {
    out <- as.data.table(cliper_output)
    return(out)
  }
  
  rbindlist(
    lapply(names(cliper_output), function(ct) {
      x <- cliper_output[[ct]]
      
      if (is.list(x) && use %in% names(x)) {
        x <- x[[use]]
      }
      
      if (!is.data.frame(x) || nrow(x) == 0) {
        return(NULL)
      }
      
      x <- as.data.table(x)
      
      if (!"Cell_Type" %in% colnames(x)) {
        x[, Cell_Type := ct]
      }
      
      x
    }),
    fill = TRUE
  )
}


fix_chr_style <- function(chr, target_chr) {
  chr <- as.character(chr)
  target_has_chr <- grepl("^chr", as.character(target_chr))
  
  if (target_has_chr) {
    paste0("chr", sub("^chr", "", chr))
  } else {
    sub("^chr", "", chr)
  }
}


make_plot_region <- function(
    region = NULL,
    gene,
    gr_anno,
    flank = 250000,
    gene_upstream = NULL,
    gene_downstream = NULL,
    gene_col = "gene_name"
) {
  if (is.null(region)) {
    if (!is.null(gene_upstream) || !is.null(gene_downstream)) {
      if (!inherits(gr_anno, "GRanges")) {
        stop("gr_anno must be a GRanges when using gene_upstream/gene_downstream.")
      }
      if (!gene_col %in% colnames(S4Vectors::mcols(gr_anno))) {
        stop("gene_col not found in gr_anno: ", gene_col)
      }
      
      gene_upstream <- if (is.null(gene_upstream)) flank else as.numeric(gene_upstream)[1]
      gene_downstream <- if (is.null(gene_downstream)) flank else as.numeric(gene_downstream)[1]
      if (!is.finite(gene_upstream) || gene_upstream < 0) {
        stop("gene_upstream should be a non-negative number.")
      }
      if (!is.finite(gene_downstream) || gene_downstream < 0) {
        stop("gene_downstream should be a non-negative number.")
      }
      
      gene_hits <- gr_anno[as.character(S4Vectors::mcols(gr_anno)[[gene_col]]) == gene]
      if (length(gene_hits) == 0L) {
        stop("Cannot find gene in gr_anno: ", gene)
      }
      
      gene_chr <- unique(as.character(GenomicRanges::seqnames(gene_hits)))
      if (length(gene_chr) != 1L) {
        stop("Gene maps to multiple chromosomes in gr_anno: ", gene)
      }
      
      region_start <- max(1, min(GenomicRanges::start(gene_hits), na.rm = TRUE) - gene_upstream)
      region_end <- max(GenomicRanges::end(gene_hits), na.rm = TRUE) + gene_downstream
      region_gr <- GenomicRanges::GRanges(
        seqnames = gene_chr,
        ranges = IRanges::IRanges(start = as.integer(region_start), end = as.integer(region_end))
      )
      names(region_gr) <- gene
      return(region_gr)
    }
    
    region_gr <- make_gene_window_for_selected(
      gr_anno = gr_anno,
      genes = gene,
      flank = flank,
      gene_col = gene_col
    )
    
    region_gr <- region_gr[gene]
    
    if (length(region_gr) != 1L) {
      stop("Cannot uniquely define region for gene: ", gene)
    }
    
    return(region_gr)
  }
  
  if (inherits(region, "GRanges")) {
    if (length(region) != 1L) {
      stop("region should be a single GRanges interval.")
    }
    return(region)
  }
  
  region_gr <- ids_to_granges_safe(region)
  
  if (length(region_gr) != 1L) {
    stop("region should resolve to exactly one interval.")
  }
  
  region_gr
}


prepare_cliper_track <- function(
    cliper_output,
    gene,
    region_gr,
    use = "summary_all",
    gene_col = "Gene",
    peak_col = "Peak",
    celltype_col = "Cell_Type",
    ppip_col = "PPIP",
    beta_col = "Posterior_b"
) {
  cliper_dt <- bind_cliper_summary(cliper_output, use = use)
  
  need <- c(gene_col, peak_col, celltype_col, ppip_col, beta_col)
  miss <- setdiff(need, colnames(cliper_dt))
  if (length(miss) > 0) {
    stop("Missing columns in CLIPER output: ", paste(miss, collapse = ", "))
  }
  
  cliper_dt <- as.data.table(cliper_dt)
  cliper_dt <- cliper_dt[get(gene_col) == gene]
  
  if (nrow(cliper_dt) == 0) {
    stop("No CLIPER rows found for gene: ", gene)
  }
  
  # Important: parse unique peaks only, otherwise merge becomes many-to-many
  unique_peaks <- unique(cliper_dt[[peak_col]])
  peak_gr <- ids_to_granges_safe(unique_peaks)
  
  peak_dt <- data.table(
    Peak_tmp = S4Vectors::mcols(peak_gr)$region_id,
    peak_chr = as.character(GenomicRanges::seqnames(peak_gr)),
    peak_start = GenomicRanges::start(peak_gr),
    peak_end = GenomicRanges::end(peak_gr)
  )
  
  setnames(peak_dt, "Peak_tmp", peak_col)
  peak_dt <- unique(peak_dt, by = peak_col)
  
  # This should now be many-to-one, not many-to-many
  cliper_dt <- merge(
    cliper_dt,
    peak_dt,
    by = peak_col,
    all.x = TRUE,
    all.y = FALSE,
    sort = FALSE
  )
  
  target_chr <- as.character(GenomicRanges::seqnames(region_gr))
  cliper_dt[, peak_chr := fix_chr_style(peak_chr, target_chr)]
  
  cliper_dt <- cliper_dt[
    peak_chr == target_chr &
      peak_end >= GenomicRanges::start(region_gr) &
      peak_start <= GenomicRanges::end(region_gr)
  ]
  
  if (nrow(cliper_dt) == 0) {
    stop("No CLIPER peaks overlap the requested plotting region.")
  }
  
  cliper_dt[, peak_mid := (peak_start + peak_end) / 2]
  cliper_dt[, celltype := as.character(get(celltype_col))]
  cliper_dt[, PPIP_plot := as.numeric(get(ppip_col))]
  cliper_dt[, Posterior_b_plot := as.numeric(get(beta_col))]
  
  cliper_dt
}

prepare_gwas_track <- function(
    gwas_finemap,
    region_gr,
    chr_col = "Chrom",
    pos_col = "Position",
    pip_col = "PIP"
) {
  target_chr <- as.character(seqnames(region_gr))
  
  if (inherits(gwas_finemap, "GRanges")) {
    gwas_dt <- as.data.table(as.data.frame(gwas_finemap))
    gwas_dt[, chr := as.character(seqnames)]
    gwas_dt[, pos := start]
  } else {
    gwas_dt <- as.data.table(gwas_finemap)
    
    if (!chr_col %in% colnames(gwas_dt)) {
      stop("chr_col not found in gwas_finemap: ", chr_col)
    }
    if (!pos_col %in% colnames(gwas_dt)) {
      stop("pos_col not found in gwas_finemap: ", pos_col)
    }
    
    gwas_dt[, chr := as.character(get(chr_col))]
    gwas_dt[, pos := as.integer(get(pos_col))]
  }
  
  if (!pip_col %in% colnames(gwas_dt)) {
    stop("pip_col not found in gwas_finemap: ", pip_col)
  }
  
  gwas_dt[, chr := fix_chr_style(chr, target_chr)]
  gwas_dt[, PIP_plot := as.numeric(get(pip_col))]
  
  gwas_dt[
    chr == target_chr &
      pos >= start(region_gr) &
      pos <= end(region_gr) &
      !is.na(PIP_plot)
  ]
}


prepare_gene_track <- function(
    gr_anno,
    region_gr,
    target_gene,
    gene_col = "gene_name",
    max_genes = 25
) {
  if (!inherits(gr_anno, "GRanges")) {
    stop("gr_anno must be a GRanges.")
  }
  
  if (!gene_col %in% colnames(S4Vectors::mcols(gr_anno))) {
    stop("gene_col not found in gr_anno: ", gene_col)
  }
  
  target_chr <- as.character(GenomicRanges::seqnames(region_gr))
  region_start <- GenomicRanges::start(region_gr)
  region_end <- GenomicRanges::end(region_gr)
  
  # Do NOT use as.data.frame(gr_anno), because duplicated names(gr_anno)
  # can trigger duplicate row.names errors.
  gene_dt <- data.table(
    chr = fix_chr_style(as.character(GenomicRanges::seqnames(gr_anno)), target_chr),
    start = GenomicRanges::start(gr_anno),
    end = GenomicRanges::end(gr_anno),
    strand = as.character(GenomicRanges::strand(gr_anno)),
    gene_name_plot = as.character(S4Vectors::mcols(gr_anno)[[gene_col]])
  )
  
  gene_dt <- gene_dt[
    !is.na(gene_name_plot) &
      gene_name_plot != "" &
      chr == target_chr &
      end >= region_start &
      start <= region_end
  ]
  
  if (nrow(gene_dt) == 0L) {
    warning("No genes found in the plotting region.")
    return(data.table(
      gene_name_plot = character(),
      start = integer(),
      end = integer(),
      strand = character(),
      gene_mid = numeric(),
      is_target = logical(),
      y = numeric()
    ))
  }
  
  gene_track_dt <- gene_dt[
    ,
    .(
      start = min(start, na.rm = TRUE),
      end = max(end, na.rm = TRUE),
      strand = {
        ss <- strand[!is.na(strand) & strand != "*"]
        if (length(ss) > 0) ss[1] else "*"
      }
    ),
    by = gene_name_plot
  ]
  
  gene_track_dt[, gene_mid := (start + end) / 2]
  gene_track_dt[, is_target := gene_name_plot == target_gene]
  
  if (nrow(gene_track_dt) > max_genes) {
    region_mid <- (region_start + region_end) / 2
    
    gene_track_dt[, dist_to_region_mid := abs(gene_mid - region_mid)]
    
    gene_track_dt <- gene_track_dt[
      order(!is_target, dist_to_region_mid)
    ][seq_len(min(max_genes, .N))]
    
    gene_track_dt[, dist_to_region_mid := NULL]
  }
  
  gene_track_dt[, y := ((seq_len(.N) - 1L) %% 4L) + 1L]
  
  gene_track_dt
}

make_binned_coverage <- function(
    object,
    region_gr,
    celltype_col,
    assay_atac = "peaks",
    celltypes = NULL,
    bin_size = 100
) {
  meta <- as.data.table(object@meta.data, keep.rownames = "cell")
  
  if (!celltype_col %in% colnames(meta)) {
    stop("celltype_col not found in object metadata: ", celltype_col)
  }
  
  if (!assay_atac %in% names(object@assays)) {
    stop(
      "assay_atac not found: ", assay_atac,
      ". Available assays: ", paste(names(object@assays), collapse = ", ")
    )
  }
  
  meta[, celltype := as.character(get(celltype_col))]
  meta <- meta[!is.na(celltype)]
  
  if (is.null(celltypes)) {
    celltypes <- unique(meta$celltype)
  }
  
  celltypes <- as.character(celltypes)
  meta <- meta[celltype %in% celltypes]
  celltypes <- celltypes[celltypes %in% unique(meta$celltype)]
  
  if (length(celltypes) == 0L) {
    stop("No selected celltypes are present in object metadata.")
  }
  
  target_chr <- as.character(GenomicRanges::seqnames(region_gr))
  region_start <- GenomicRanges::start(region_gr)
  region_end <- GenomicRanges::end(region_gr)
  
  starts <- seq(region_start, region_end, by = bin_size)
  ends <- pmin(starts + bin_size - 1L, region_end)
  
  bins <- GenomicRanges::GRanges(
    seqnames = target_chr,
    ranges = IRanges::IRanges(start = starts, end = ends)
  )
  
  names(bins) <- paste0("bin_", seq_along(bins))
  
  fragments <- Signac::Fragments(object[[assay_atac]])
  
  if (length(fragments) == 0L) {
    stop("No fragments found in object[[", assay_atac, "]].")
  }
  
  cells_use <- intersect(meta$cell, colnames(object))
  
  if (length(cells_use) == 0L) {
    stop("No cells overlap between metadata and object colnames.")
  }
  
  mat <- Signac::FeatureMatrix(
    fragments = fragments,
    features = bins,
    cells = cells_use
  )
  
  # Critical fix: FeatureMatrix may not preserve names(bins)
  if (nrow(mat) != length(bins)) {
    stop("FeatureMatrix returned a different number of rows than bins.")
  }
  
  rownames(mat) <- names(bins)
  
  meta_use <- meta[match(colnames(mat), cell), .(cell, celltype)]
  keep <- !is.na(meta_use$celltype)
  
  mat <- mat[, keep, drop = FALSE]
  meta_use <- meta_use[keep]
  
  group <- factor(meta_use$celltype, levels = celltypes)
  
  keep2 <- !is.na(group)
  mat <- mat[, keep2, drop = FALSE]
  group <- group[keep2]
  
  if (ncol(mat) == 0L) {
    stop("No cells remain after matching cell type labels.")
  }
  
  group_mat <- Matrix::sparseMatrix(
    i = seq_along(group),
    j = as.integer(group),
    x = 1,
    dims = c(length(group), length(celltypes)),
    dimnames = list(colnames(mat), celltypes)
  )
  
  cov_by_ct <- mat %*% group_mat
  
  n_by_ct <- as.numeric(tabulate(as.integer(group), nbins = length(celltypes)))
  names(n_by_ct) <- celltypes
  
  cov_by_ct <- sweep(
    as.matrix(cov_by_ct),
    2,
    pmax(n_by_ct, 1),
    "/"
  )
  
  rownames(cov_by_ct) <- names(bins)
  
  cov_dt <- as.data.table(as.table(cov_by_ct))
  setnames(cov_dt, c("bin", "celltype", "coverage"))
  
  cov_dt[, bin := as.character(bin)]
  cov_dt[, celltype := as.character(celltype)]
  cov_dt[, coverage := as.numeric(coverage)]
  
  bin_dt <- data.table(
    bin = names(bins),
    x = (GenomicRanges::start(bins) + GenomicRanges::end(bins)) / 2,
    bin_start = GenomicRanges::start(bins),
    bin_end = GenomicRanges::end(bins)
  )
  
  cov_dt <- merge(cov_dt, bin_dt, by = "bin", all.x = TRUE, sort = FALSE)
  
  if (any(is.na(cov_dt$x))) {
    stop("Coverage bins failed to merge with bin coordinates. This should not happen after rownames(mat) <- names(bins).")
  }
  
  cov_dt[, coverage_scaled := {
    mx <- max(coverage, na.rm = TRUE)
    if (is.finite(mx) && mx > 0) coverage / mx else rep(0, .N)
  }, by = celltype]
  
  cov_dt[, coverage_scaled := as.numeric(coverage_scaled)]
  cov_dt[!is.finite(coverage_scaled), coverage_scaled := 0]
  
  cov_dt
}

get_c4a_colors <- function(
    palette,
    n,
    type = "cat",
    reverse = FALSE
) {
  if (!requireNamespace("cols4all", quietly = TRUE)) {
    stop("Package 'cols4all' is required. Please install it with install.packages('cols4all').")
  }
  
  cols <- tryCatch(
    cols4all::c4a(
      palette = palette,
      n = n,
      type = type,
      reverse = reverse
    ),
    error = function(e1) {
      tryCatch(
        cols4all::c4a(
          palette,
          n = n,
          type = type,
          reverse = reverse
        ),
        error = function(e2) {
          stop(
            "Failed to get colors from cols4all palette: ",
            palette,
            "\nOriginal error: ",
            conditionMessage(e2)
          )
        }
      )
    }
  )
  
  as.character(cols)
}


get_c4a_celltype_colors <- function(
    celltypes,
    palette = "carto.safe",
    reverse = FALSE
) {
  celltypes <- as.character(celltypes)
  
  cols <- get_c4a_colors(
    palette = palette,
    n = length(celltypes),
    type = "cat",
    reverse = reverse
  )
  
  setNames(cols, celltypes)
}

make_one_celltype_block <- function(
    ct,
    cov_dt,
    cliper_plot_dt,
    region_start,
    region_end,
    ct_color,
    beta_cols,
    beta_limits,
    show_x = FALSE,
    show_legend = FALSE,
    show_peak_labels = FALSE,
    highlighted_peak_regions = NULL,
    coverage_alpha = 0.45,
    ppip_point_size = 2.0,
    highlight_peak_band_width = NULL,
    highlight_peak_fill = "#9ECAE1",
    highlight_alpha = 0.35,
    plot_base_size = 12,
    celltype_label_size = 12,
    axis_linewidth = 0.45
) {
  ct <- as.character(ct)
  
  cov_sub <- cov_dt[as.character(celltype) == ct]
  ppip_sub <- cliper_plot_dt[as.character(celltype) == ct]
  ppip_high <- ppip_sub[highlight_peak %in% TRUE]
  peak_high_regions <- unique(ppip_high[, .(
    peak_mid,
    Peak
  )])
  if (nrow(peak_high_regions) > 0L) {
    peak_high_regions[, xmin := pmax(region_start, peak_mid - highlight_peak_band_width / 2)]
    peak_high_regions[, xmax := pmin(region_end, peak_mid + highlight_peak_band_width / 2)]
  } else {
    peak_high_regions[, xmin := numeric()]
    peak_high_regions[, xmax := numeric()]
  }
  coverage_high_regions <- if (show_peak_labels && !is.null(highlighted_peak_regions)) {
    highlighted_peak_regions
  } else {
    peak_high_regions
  }
  if (nrow(coverage_high_regions) > 0L) {
    coverage_high_regions[, label_x := pmin(pmax(peak_mid, region_start), region_end)]
  } else {
    coverage_high_regions[, label_x := numeric()]
  }
  
  common_theme <- theme_classic(base_size = plot_base_size) +
    theme(
      panel.grid = element_blank(),
      plot.margin = margin(1.5, 2, 1.5, 2),
      axis.title.x = element_blank(),
      axis.line.x = element_blank(),
      axis.ticks.x = if (show_x) element_line(linewidth = axis_linewidth) else element_blank(),
      axis.text.x = if (show_x) element_text(size = plot_base_size - 1) else element_blank(),
      axis.text.y = element_text(size = plot_base_size - 2),
      axis.ticks.y = element_line(linewidth = axis_linewidth),
      axis.line.y = element_line(linewidth = axis_linewidth)
    )
  
  p_cov_ct <- ggplot(cov_sub, aes(x = x, y = coverage_scaled)) +
    geom_rect(
      data = coverage_high_regions,
      aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
      inherit.aes = FALSE,
      fill = highlight_peak_fill,
      alpha = highlight_alpha,
      color = NA
    ) +
    geom_area(
      fill = ct_color,
      alpha = coverage_alpha,
      linewidth = 0
    ) +
    geom_line(
      color = ct_color,
      linewidth = 0.30,
      alpha = 0.95
    ) +
    geom_text(
      data = if (show_peak_labels) coverage_high_regions else coverage_high_regions[0],
      aes(x = label_x, y = 0.98, label = Peak),
      inherit.aes = FALSE,
      color = "black",
      size = max(3.0, plot_base_size * 0.28),
      fontface = "bold",
      angle = 0,
      hjust = 0.5,
      vjust = 1
    ) +
    scale_y_continuous(
      breaks = c(0, 1),
      limits = c(0, 1),
      expand = expansion(mult = c(0, 0.03))
    ) +
    coord_cartesian(
      xlim = c(region_start, region_end),
      expand = FALSE,
      clip = "off"
    ) +
    labs(
      x = NULL,
      y = ct
    ) +
    common_theme +
    theme(
      plot.margin = margin(1.5, 2, 4, 2),
      axis.title.y = element_text(
        angle = 0,
        hjust = 1,
        vjust = 0.62,
        size = celltype_label_size,
        face = "bold",
        color = ct_color,
        margin = margin(r = 11)
      )
    )
  
  if (nrow(ppip_sub) == 0L) {
    return(p_cov_ct)
  }
  
  p_ppip_ct <- ggplot(ppip_sub) +
    geom_rect(
      data = peak_high_regions,
      aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
      inherit.aes = FALSE,
      fill = highlight_peak_fill,
      alpha = highlight_alpha,
      color = NA
    ) +
    geom_segment(
      aes(
        x = peak_start,
        xend = peak_end,
        y = 0.02,
        yend = 0.02
      ),
      linewidth = 0.30,
      alpha = 0.40
    ) +
    geom_point(
      aes(
        x = peak_mid,
        y = PPIP_plot,
        fill = Posterior_b_plot
      ),
      shape = 21,
      color = "black",
      stroke = 0.30,
      size = ppip_point_size,
      alpha = 0.95
    ) +
    scale_fill_gradientn(
      colours = beta_cols,
      limits = beta_limits,
      oob = scales::squish,
      name = "Posterior b",
      guide = if (show_legend) {
        guide_colorbar(
          title.position = "top",
          order = 3,
          barheight = unit(44, "pt"),
          barwidth = unit(6, "pt")
        )
      } else {
        "none"
      }
    ) +
    scale_y_continuous(
      breaks = c(0, 1),
      limits = c(-0.05, 1),
      expand = expansion(mult = c(0, 0.03))
    ) +
    coord_cartesian(
      xlim = c(region_start, region_end),
      expand = FALSE,
      clip = "off"
    ) +
    labs(
      x = NULL,
      y = "CLIPER\nPPIP"
    ) +
    common_theme +
    theme(
      axis.line.x = element_line(linewidth = axis_linewidth),
      axis.ticks.x = element_line(linewidth = axis_linewidth),
      axis.text.x = element_blank(),
      plot.margin = margin(5, 2, 4, 2),
      axis.title.y = element_text(
        angle = 0,
        hjust = 1,
        vjust = 0.5,
        size = plot_base_size - 1,
        margin = margin(r = 11)
      )
    )
  
  p_cov_ct / p_ppip_ct +
    patchwork::plot_layout(heights = c(1.00, 0.67))
}


plot_cliper_p2g <- function(
    object,
    cliper_output,
    gene,
    gr_anno,
    gwas_finemap = NULL,
    region = NULL,
    flank = 250000,
    gene_upstream = NULL,
    gene_downstream = NULL,
    celltype_col = "celltype",
    celltypes = NULL,
    assay_atac = "ATAC",
    assay_rna = "RNA",
    gene_anno_col = "gene_name",
    gwas_chr_col = "Chrom",
    gwas_pos_col = "Position",
    gwas_pip_col = "PIP",
    bin_size = 100,
    use_cliper = "summary_all",
    max_genes = 25,
    c4a_celltype_palette = "carto.safe",
    c4a_beta_palette = "hcl.blue_red",
    c4a_celltype_reverse = FALSE,
    c4a_beta_reverse = FALSE,
    expression_width = NULL,
    browser_width = 9.50,
    coverage_alpha = 0.45,
    ppip_point_size = 2.0,
    highlight_peaks = NULL,
    highlight_snps = NULL,
    gwas_snp_col = NULL,
    highlight_peak_band_width = NULL,
    highlight_peak_fill = "#9ECAE1",
    highlight_snp_fill = "#E64B35",
    highlight_alpha = 0.35,
    dotplot_cols = c("grey90", "#B2182B"),
    dotplot_dot_scale = 6,
    dotplot_scale_min = 0,
    dotplot_scale_max = 100,
    plot_base_size = 12,
    celltype_label_size = 12,
    axis_linewidth = 0.45,
    return_data = FALSE
) {
  if (!requireNamespace("cols4all", quietly = TRUE)) {
    stop("Package 'cols4all' is required. Please install it with install.packages('cols4all').")
  }
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("Package 'patchwork' is required.")
  }
  if (!requireNamespace("Signac", quietly = TRUE)) {
    stop("Package 'Signac' is required.")
  }
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Package 'Seurat' is required.")
  }
  if (is.null(expression_width)) {
    expression_width <- browser_width / 5
  }
  expression_width <- as.numeric(expression_width)[1]
  if (!is.finite(expression_width) || expression_width <= 0) {
    stop("expression_width should be a positive number.")
  }
  plot_gwas <- !is.null(gwas_finemap)
  region_gr <- make_plot_region(
    region = region,
    gene = gene,
    gr_anno = gr_anno,
    flank = flank,
    gene_upstream = gene_upstream,
    gene_downstream = gene_downstream,
    gene_col = gene_anno_col
  )
  
  target_chr <- as.character(GenomicRanges::seqnames(region_gr))
  region_start <- GenomicRanges::start(region_gr)
  region_end <- GenomicRanges::end(region_gr)
  region_width <- region_end - region_start
  
  if (is.null(highlight_peak_band_width)) {
    highlight_peak_band_width <- max(
      1,
      region_width * ppip_point_size / (browser_width * 25.4)
    )
  }
  highlight_peak_band_width <- as.numeric(highlight_peak_band_width)[1]
  if (!is.finite(highlight_peak_band_width) || highlight_peak_band_width <= 0) {
    stop("highlight_peak_band_width should be a positive number.")
  }
  
  region_string <- paste0(target_chr, ":", region_start, "-", region_end)
  
  cliper_dt <- prepare_cliper_track(
    cliper_output = cliper_output,
    gene = gene,
    region_gr = region_gr,
    use = use_cliper
  )
  
  meta <- data.table::as.data.table(object@meta.data, keep.rownames = "cell")
  
  if (!celltype_col %in% colnames(meta)) {
    stop("celltype_col not found in object metadata: ", celltype_col)
  }
  
  meta[, celltype := as.character(get(celltype_col))]
  
  if (is.null(celltypes)) {
    celltype_raw <- object@meta.data[[celltype_col]]
    if (is.factor(celltype_raw)) {
      celltypes <- levels(celltype_raw)
      celltypes <- celltypes[celltypes %in% meta[!is.na(celltype), unique(celltype)]]
    } else {
      celltypes <- sort(unique(meta[!is.na(celltype), celltype]))
    }
  }
  
  celltypes <- as.character(celltypes)
  
  if (length(celltypes) == 0L) {
    stop("No cell types found in object metadata.")
  }
  
  cliper_dt <- cliper_dt[celltype %in% celltypes]
  
  if (nrow(cliper_dt) == 0L) {
    stop("No CLIPER rows remain after filtering to selected celltypes.")
  }
  
  cliper_dt[, celltype := factor(as.character(celltype), levels = celltypes)]
  
  cov_dt <- make_binned_coverage(
    object = object,
    region_gr = region_gr,
    celltype_col = celltype_col,
    assay_atac = assay_atac,
    celltypes = celltypes,
    bin_size = bin_size
  )
  
  cov_dt[, celltype := factor(as.character(celltype), levels = celltypes)]
  
  gwas_dt <- NULL
  if (plot_gwas) {
    gwas_dt <- prepare_gwas_track(
      gwas_finemap = gwas_finemap,
      region_gr = region_gr,
      chr_col = gwas_chr_col,
      pos_col = gwas_pos_col,
      pip_col = gwas_pip_col
    )
  }
  
  gene_track_dt <- prepare_gene_track(
    gr_anno = gr_anno,
    region_gr = region_gr,
    target_gene = gene,
    gene_col = gene_anno_col,
    max_genes = max_genes
  )
  
  cliper_plot_dt <- cliper_dt[
    !is.na(PPIP_plot) &
      !is.na(Posterior_b_plot)
  ]
  
  if (nrow(cliper_plot_dt) == 0L) {
    stop("No non-missing CLIPER PPIP / Posterior_b values to plot.")
  }
  
  highlight_peaks <- unique(as.character(highlight_peaks))
  highlight_peaks <- highlight_peaks[!is.na(highlight_peaks) & highlight_peaks != ""]
  cliper_plot_dt[, peak_coord := paste0(peak_chr, ":", peak_start, "-", peak_end)]
  cliper_plot_dt[, peak_coord_nochr := paste0(sub("^chr", "", peak_chr), ":", peak_start, "-", peak_end)]
  cliper_plot_dt[, highlight_peak := FALSE]
  if (length(highlight_peaks) > 0L) {
    cliper_plot_dt[
      as.character(Peak) %in% highlight_peaks |
        peak_coord %in% highlight_peaks |
        peak_coord_nochr %in% highlight_peaks,
      highlight_peak := TRUE
    ]
  }
  
  if (plot_gwas) {
    highlight_snps <- unique(as.character(highlight_snps))
    highlight_snps <- highlight_snps[!is.na(highlight_snps) & highlight_snps != ""]
    gwas_dt[, snp_coord := paste0(chr, ":", pos)]
    gwas_dt[, snp_coord_nochr := paste0(sub("^chr", "", chr), ":", pos)]
    gwas_dt[, highlight_snp := FALSE]
    if (length(highlight_snps) > 0L) {
      snp_cols <- if (!is.null(gwas_snp_col)) {
        if (!gwas_snp_col %in% colnames(gwas_dt)) {
          stop("gwas_snp_col not found in gwas_finemap: ", gwas_snp_col)
        }
        gwas_snp_col
      } else {
        intersect(
          c("SNP", "snp", "rsid", "RSID", "rsID", "Name", "name", "Index", "index",
            "variant_id", "Variant", "MarkerName", "ID", "id"),
          colnames(gwas_dt)
        )
      }
      
      gwas_dt[
        snp_coord %in% highlight_snps |
          snp_coord_nochr %in% highlight_snps |
          as.character(pos) %in% highlight_snps,
        highlight_snp := TRUE
      ]
      
      if (length(snp_cols) > 0L) {
        gwas_dt[
          Reduce(
            `|`,
            lapply(snp_cols, function(cc) as.character(get(cc)) %in% highlight_snps)
          ),
          highlight_snp := TRUE
        ]
      }
    }
    
    snp_label_cols <- if (!is.null(gwas_snp_col) && gwas_snp_col %in% colnames(gwas_dt)) {
      gwas_snp_col
    } else {
      intersect(
        c("SNP", "snp", "rsid", "RSID", "rsID", "Name", "name", "Index", "index",
          "variant_id", "Variant", "MarkerName", "ID", "id"),
        colnames(gwas_dt)
      )
    }
    gwas_dt[, highlight_snp_label := fifelse(
      highlight_snp %in% TRUE,
      if (length(snp_label_cols) > 0L) as.character(get(snp_label_cols[1])) else snp_coord,
      NA_character_
    )]
  }
  
  ct_cols <- get_c4a_celltype_colors(
    celltypes = celltypes,
    palette = c4a_celltype_palette,
    reverse = c4a_celltype_reverse
  )
  
  beta_cols <- get_c4a_colors(
    palette = c4a_beta_palette,
    n = 11,
    type = "div",
    reverse = c4a_beta_reverse
  )
  
  beta_abs <- max(abs(cliper_plot_dt$Posterior_b_plot), na.rm = TRUE)
  if (!is.finite(beta_abs) || beta_abs == 0) {
    beta_abs <- 1
  }
  beta_limits <- c(-beta_abs, beta_abs)
  ppip_legend_i <- match(
    TRUE,
    celltypes %in% unique(as.character(cliper_plot_dt$celltype))
  )
  celltypes_with_ppip <- unique(as.character(cliper_plot_dt$celltype))
  n_celltypes_with_ppip <- sum(celltypes %in% celltypes_with_ppip)
  browser_height <- max(4.5, (length(celltypes) + 0.67 * n_celltypes_with_ppip) * 0.70)
  highlighted_peak_regions <- unique(cliper_plot_dt[
    highlight_peak %in% TRUE,
    .(
      peak_mid,
      Peak
    )
  ])
  if (nrow(highlighted_peak_regions) > 0L) {
    highlighted_peak_regions[, xmin := pmax(region_start, peak_mid - highlight_peak_band_width / 2)]
    highlighted_peak_regions[, xmax := pmin(region_end, peak_mid + highlight_peak_band_width / 2)]
  } else {
    highlighted_peak_regions[, xmin := numeric()]
    highlighted_peak_regions[, xmax := numeric()]
  }
  
  object_expr <- object
  
  cells_expr <- rownames(object_expr@meta.data)[
    as.character(object_expr@meta.data[[celltype_col]]) %in% celltypes
  ]
  
  if (length(cells_expr) == 0L) {
    stop("No cells remain for expression plot after filtering to selected celltypes.")
  }
  
  object_expr <- subset(object_expr, cells = cells_expr)
  
  object_expr@meta.data[[celltype_col]] <- factor(
    as.character(object_expr@meta.data[[celltype_col]]),
    levels = celltypes
  )
  
  object_expr@meta.data[[celltype_col]] <- droplevels(
    object_expr@meta.data[[celltype_col]]
  )
  
  Seurat::DefaultAssay(object_expr) <- assay_rna
  
  p_expr <- Seurat::DotPlot(
    object = object_expr,
    features = gene,
    group.by = celltype_col,
    assay = assay_rna,
    cols = dotplot_cols,
    dot.scale = dotplot_dot_scale,
    scale.min = dotplot_scale_min,
    scale.max = dotplot_scale_max
  ) +
    ggplot2::scale_y_discrete(limits = rev(celltypes)) +
    ggplot2::labs(
      title = paste0(gene, "\nGene expression"),
      x = NULL,
      y = NULL
    ) +
    ggplot2::guides(
      color = ggplot2::guide_colorbar(order = 1),
      colour = ggplot2::guide_colorbar(order = 1),
      size = ggplot2::guide_legend(order = 2)
    ) +
    ggplot2::theme_classic(base_size = plot_base_size) +
    ggplot2::theme(
      legend.position = "right",
      plot.title = ggplot2::element_text(size = plot_base_size, face = "bold", hjust = 0.5),
      axis.title = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      axis.line.x = ggplot2::element_blank(),
      axis.line.y = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      strip.background = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(2, 2, 2, 2)
    )
  
  ct_blocks <- lapply(seq_along(celltypes), function(i) {
    ct <- as.character(celltypes[i])
    
    make_one_celltype_block(
      ct = ct,
      cov_dt = cov_dt,
      cliper_plot_dt = cliper_plot_dt,
      region_start = region_start,
      region_end = region_end,
      ct_color = unname(ct_cols[ct]),
      beta_cols = beta_cols,
      beta_limits = beta_limits,
      show_x = FALSE,
      show_legend = i == ppip_legend_i,
      show_peak_labels = i == 1,
      highlighted_peak_regions = highlighted_peak_regions,
      coverage_alpha = coverage_alpha,
      ppip_point_size = ppip_point_size,
      highlight_peak_band_width = highlight_peak_band_width,
      highlight_peak_fill = highlight_peak_fill,
      highlight_alpha = highlight_alpha,
      plot_base_size = plot_base_size,
      celltype_label_size = celltype_label_size,
      axis_linewidth = axis_linewidth
    )
  })
  
  p_browser <- patchwork::wrap_plots(
    ct_blocks,
    ncol = 1
  ) +
    patchwork::plot_layout(
      heights = ifelse(celltypes %in% celltypes_with_ppip, 1.67, 1.00)
    )
  
  peak_region_dt <- unique(cliper_plot_dt[, .(
    peak_start,
    peak_end,
    peak_mid,
    highlight_peak
  )])
  peak_region_dt[, y := 1]
  
  p_peaks <- ggplot2::ggplot(peak_region_dt) +
    ggplot2::geom_rect(
      data = unique(peak_region_dt[highlight_peak %in% TRUE, .(xmin = peak_start, xmax = peak_end)]),
      ggplot2::aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
      inherit.aes = FALSE,
      fill = highlight_peak_fill,
      alpha = highlight_alpha,
      color = NA
    ) +
    ggplot2::geom_segment(
      ggplot2::aes(
        x = peak_start,
        xend = peak_end,
        y = y,
        yend = y,
        linewidth = highlight_peak
      ),
      color = "#5A5A5A",
      lineend = "round"
    ) +
    ggplot2::scale_linewidth_manual(values = c(`FALSE` = 0.55, `TRUE` = 1.1), guide = "none") +
    ggplot2::coord_cartesian(
      xlim = c(region_start, region_end),
      ylim = c(0.75, 1.25),
      expand = FALSE,
      clip = "off"
    ) +
    ggplot2::labs(x = NULL, y = "Peaks\nRegion") +
    ggplot2::theme_classic(base_size = plot_base_size) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      axis.line = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_text(
        angle = 0,
        hjust = 1,
        vjust = 0.5,
        size = plot_base_size - 1,
        margin = ggplot2::margin(r = 7)
      ),
      plot.margin = ggplot2::margin(1, 2, 1, 2)
    )
  
  if (plot_gwas) {
    p_gwas <- ggplot2::ggplot(gwas_dt, ggplot2::aes(x = pos, y = PIP_plot)) +
      ggplot2::geom_segment(
        ggplot2::aes(xend = pos, y = 0, yend = PIP_plot),
        alpha = 0.45,
        linewidth = 0.35
      ) +
      ggplot2::geom_point(
        shape = 21,
        fill = "white",
        color = "black",
        stroke = 0.25,
        size = 1.85,
        alpha = 0.90
      ) +
      ggplot2::geom_point(
        data = gwas_dt[highlight_snp %in% TRUE],
        ggplot2::aes(x = pos, y = PIP_plot),
        inherit.aes = FALSE,
        shape = 21,
        fill = highlight_snp_fill,
        color = "#B2182B",
        stroke = 0.55,
        size = 2.8,
        alpha = 1
      ) +
      ggplot2::geom_text(
        data = gwas_dt[highlight_snp %in% TRUE],
        ggplot2::aes(
          x = pos,
          y = PIP_plot,
          label = highlight_snp_label
        ),
        inherit.aes = FALSE,
        color = "#B2182B",
        size = max(3.0, plot_base_size * 0.28),
        fontface = "bold",
        nudge_y = max(1, gwas_dt$PIP_plot, na.rm = TRUE) * 0.06,
        angle = 0,
        hjust = 0.5,
        vjust = 0
      ) +
      ggplot2::coord_cartesian(
        xlim = c(region_start, region_end),
        ylim = c(0, max(1, gwas_dt$PIP_plot, na.rm = TRUE) * 1.15),
        expand = FALSE,
        clip = "off"
      ) +
      ggplot2::labs(
        x = NULL,
        y = "GWAS\nPIP"
      ) +
      ggplot2::theme_classic(base_size = plot_base_size) +
      ggplot2::theme(
        panel.grid = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank(),
        axis.title.y = ggplot2::element_text(
          angle = 0,
          hjust = 1,
          vjust = 0.5,
          size = plot_base_size - 1,
          margin = ggplot2::margin(r = 7)
        ),
        axis.text.y = ggplot2::element_text(size = plot_base_size - 2),
        axis.ticks.y = ggplot2::element_line(linewidth = axis_linewidth),
        axis.line.y = ggplot2::element_line(linewidth = axis_linewidth),
        axis.line.x = ggplot2::element_blank(),
        plot.margin = ggplot2::margin(2, 2, 2, 2)
      )
  }
  
  p_gene <- ggplot2::ggplot(gene_track_dt) +
    ggplot2::geom_segment(
      data = gene_track_dt[strand != "-"],
      ggplot2::aes(
        x = start,
        xend = end,
        y = y,
        yend = y,
        linewidth = is_target,
        color = is_target
      ),
      arrow = ggplot2::arrow(
        length = grid::unit(0.055, "inches"),
        ends = "last",
        type = "closed"
      ),
      lineend = "round"
    ) +
    ggplot2::geom_segment(
      data = gene_track_dt[strand == "-"],
      ggplot2::aes(
        x = start,
        xend = end,
        y = y,
        yend = y,
        linewidth = is_target,
        color = is_target
      ),
      arrow = ggplot2::arrow(
        length = grid::unit(0.055, "inches"),
        ends = "first",
        type = "closed"
      ),
      lineend = "round"
    ) +
    ggplot2::geom_text(
      ggplot2::aes(
        x = gene_mid,
        y = y + 0.22,
        label = gene_name_plot,
        fontface = ifelse(is_target, "bold", "plain"),
        color = is_target
      ),
      size = max(3.1, plot_base_size * 0.32)
    ) +
    ggplot2::scale_linewidth_manual(
      values = c(`FALSE` = 0.45, `TRUE` = 1.05),
      guide = "none"
    ) +
    ggplot2::scale_color_manual(
      values = c(`FALSE` = "#5F6B7A", `TRUE` = "#B2182B"),
      guide = "none"
    ) +
    ggplot2::coord_cartesian(
      xlim = c(region_start, region_end),
      expand = FALSE,
      clip = "off"
    ) +
    ggplot2::labs(
      x = paste0(target_chr, " position"),
      y = NULL
    ) +
    ggplot2::theme_classic(base_size = plot_base_size) +
    ggplot2::theme(
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      axis.line.y = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(size = plot_base_size - 1),
      axis.ticks.x = ggplot2::element_line(linewidth = axis_linewidth),
      axis.line.x = ggplot2::element_line(linewidth = axis_linewidth),
      panel.grid = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(2, 2, 2, 2)
    )
  
  p_top <- (p_browser | p_expr) +
    patchwork::plot_layout(widths = c(browser_width, expression_width))
  p_peaks_row <- (p_peaks | patchwork::plot_spacer()) +
    patchwork::plot_layout(widths = c(browser_width, expression_width))
  p_gene_row <- (p_gene | patchwork::plot_spacer()) +
    patchwork::plot_layout(widths = c(browser_width, expression_width))
  
  if (plot_gwas) {
    p_gwas_row <- (p_gwas | patchwork::plot_spacer()) +
      patchwork::plot_layout(widths = c(browser_width, expression_width))
    final_plot <- p_top /
      p_peaks_row /
      p_gwas_row /
      p_gene_row +
      patchwork::plot_layout(
        heights = c(
          browser_height,
          0.28,
          0.95,
          0.80
        ),
        guides = "collect"
      )
  } else {
    final_plot <- p_top /
      p_peaks_row /
      p_gene_row +
      patchwork::plot_layout(
        heights = c(
          browser_height,
          0.28,
          0.80
        ),
        guides = "collect"
      )
  }
  
  final_plot <- final_plot +
    patchwork::plot_annotation(
      title = paste0(gene, " CLIPER peak-to-gene fine-mapping"),
      subtitle = region_string
    ) &
    ggplot2::theme(
      legend.position = "right",
      plot.title = ggplot2::element_text(face = "bold", size = 12),
      plot.subtitle = ggplot2::element_text(size = 9)
    )
  
  if (return_data) {
    return(list(
      plot = final_plot,
      region = region_gr,
      cliper = cliper_dt,
      cliper_plot_dt = cliper_plot_dt,
      coverage = cov_dt,
      gwas = gwas_dt,
      peaks = peak_region_dt,
      genes = gene_track_dt,
      celltype_colors = ct_cols,
      beta_colors = beta_cols
    ))
  }
  
  final_plot
}
