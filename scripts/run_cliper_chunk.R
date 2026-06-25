#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

usage <- function() {
  cat(
    "Usage:\n",
    "  Rscript scripts/run_cliper_chunk.R run \\\n",
    "    --cliper_obj cliper_obj.rds --genes genes.txt --gr_anno gr_anno.rds \\\n",
    "    --outdir cliper_chunks --chunk_id 1 --n_chunks 10 [Run_CLIPER options]\n",
    "\n",
    "  Rscript scripts/run_cliper_chunk.R merge \\\n",
    "    --outdir cliper_chunks --merged_rds cliper_merged.rds \\\n",
    "    --summary_all_csv summary_all.csv --summary_csv cliper_summary.csv \\\n",
    "    --select_csv cliper_select.csv --info_csv summary_info.csv\n",
    "\n",
    "Common Run_CLIPER options:\n",
    "  --flank 500000 --p1 0.8 --K 5 --n_iter 10000 --burn_in 5000\n",
    "  --alpha_conc 5 --tau2 5000 --rho0 0.5 --posterior_b_cutoff 0.1 --pip_cutoff 0.8\n",
    "  --no_scale to disable X scaling; --scale_y to scale y after centering\n",
    sep = ""
  )
  quit(status = 1)
}

if (length(args) < 1) usage()
mode <- args[[1]]
args <- args[-1]

parse_args <- function(x) {
  out <- list()
  i <- 1L
  while (i <= length(x)) {
    key <- x[[i]]
    if (!startsWith(key, "--")) {
      stop("Unexpected argument: ", key, call. = FALSE)
    }
    key <- sub("^--", "", key)
    if (i == length(x) || startsWith(x[[i + 1L]], "--")) {
      out[[key]] <- TRUE
      i <- i + 1L
    } else {
      out[[key]] <- x[[i + 1L]]
      i <- i + 2L
    }
  }
  out
}

opt <- parse_args(args)

need <- function(name) {
  if (is.null(opt[[name]]) || !nzchar(opt[[name]])) {
    stop("Missing required option --", name, call. = FALSE)
  }
  opt[[name]]
}

as_num <- function(name, default) {
  if (is.null(opt[[name]])) return(default)
  as.numeric(opt[[name]])
}

as_int <- function(name, default) {
  as.integer(as_num(name, default))
}

read_gene_list <- function(path) {
  genes <- readLines(path, warn = FALSE)
  genes <- trimws(genes)
  genes <- genes[nzchar(genes)]
  unique(genes)
}

split_genes <- function(genes, chunk_id, n_chunks) {
  if (length(genes) == 0L) {
    return(character(0))
  }
  idx <- seq_along(genes)
  keep <- ((idx - 1L) %% n_chunks) + 1L == chunk_id
  genes[keep]
}

bind_or_empty <- function(x) {
  x <- Filter(function(z) !is.null(z) && nrow(z) > 0L, x)
  if (length(x) == 0L) {
    return(data.frame())
  }
  dplyr::bind_rows(x)
}

merge_results <- function(files) {
  chunk_results <- lapply(files, readRDS)
  cts <- sort(unique(unlist(lapply(chunk_results, names), use.names = FALSE)))

  merged <- setNames(vector("list", length(cts)), cts)
  for (ct in cts) {
    per_ct <- lapply(chunk_results, function(x) x[[ct]])
    per_ct <- Filter(Negate(is.null), per_ct)

    merged[[ct]] <- list(
      summary_all = bind_or_empty(lapply(per_ct, `[[`, "summary_all")),
      cliper_summary = bind_or_empty(lapply(per_ct, `[[`, "cliper_summary")),
      cliper_select = bind_or_empty(lapply(per_ct, `[[`, "cliper_select")),
      summary_info = bind_or_empty(lapply(per_ct, `[[`, "summary_info"))
    )
  }

  merged
}

if (identical(mode, "run")) {
  suppressPackageStartupMessages(library(CLIPER))

  cliper_obj_path <- need("cliper_obj")
  genes_path <- need("genes")
  gr_anno_path <- need("gr_anno")
  outdir <- need("outdir")
  chunk_id <- as_int("chunk_id", NA_integer_)
  n_chunks <- as_int("n_chunks", NA_integer_)

  if (!is.finite(chunk_id) || !is.finite(n_chunks) || chunk_id < 1L || n_chunks < 1L || chunk_id > n_chunks) {
    stop("--chunk_id must be between 1 and --n_chunks.", call. = FALSE)
  }

  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

  genes <- read_gene_list(genes_path)
  genes_chunk <- split_genes(genes, chunk_id, n_chunks)

  message("[CLIPER chunk ", chunk_id, "/", n_chunks, "] genes: ", length(genes_chunk))

  out_file <- file.path(outdir, sprintf("cliper_chunk_%03d_of_%03d.rds", chunk_id, n_chunks))
  gene_file <- file.path(outdir, sprintf("genes_chunk_%03d_of_%03d.txt", chunk_id, n_chunks))
  writeLines(genes_chunk, gene_file)

  if (length(genes_chunk) == 0L) {
    saveRDS(list(), out_file)
    quit(status = 0)
  }

  cliper_obj <- readRDS(cliper_obj_path)
  gr_anno <- readRDS(gr_anno_path)

  set.seed(as_int("seed", 2001L) + chunk_id)

  result <- Run_CLIPER(
    cliper_obj = cliper_obj,
    flank = as_num("flank", 500000),
    gene_list = genes_chunk,
    gr_anno = gr_anno,
    p1 = as_num("p1", 0.8),
    K = as_int("K", 5L),
    n_iter = as_int("n_iter", 10000L),
    burn_in = as_int("burn_in", 5000L),
    alpha_conc = as_num("alpha_conc", 5),
    tau2 = as_num("tau2", 5000),
    rho0 = as_num("rho0", 0.5),
    add_b_mu = isTRUE(opt[["add_b_mu"]]),
    scale = !isTRUE(opt[["no_scale"]]),
    scale_y = isTRUE(opt[["scale_y"]]),
    seed = as_int("seed", 2001L) + chunk_id,
    posterior_b_cutoff = as_num("posterior_b_cutoff", 0.1),
    pip_cutoff = as_num("pip_cutoff", 0.8)
  )

  saveRDS(result, out_file)
  message("[CLIPER chunk ", chunk_id, "/", n_chunks, "] saved: ", out_file)
  quit(status = 0)
}

if (identical(mode, "merge")) {
  suppressPackageStartupMessages(library(dplyr))

  outdir <- need("outdir")
  merged_rds <- need("merged_rds")
  summary_all_csv <- opt[["summary_all_csv"]]
  summary_csv <- opt[["summary_csv"]]
  select_csv <- opt[["select_csv"]]
  info_csv <- opt[["info_csv"]]

  files <- list.files(outdir, pattern = "^cliper_chunk_.*\\.rds$", full.names = TRUE)
  if (length(files) == 0L) {
    stop("No chunk RDS files found in ", outdir, call. = FALSE)
  }
  files <- sort(files)

  merged <- merge_results(files)
  saveRDS(merged, merged_rds)
  message("[CLIPER merge] saved: ", merged_rds)

  all_summary_all <- bind_or_empty(lapply(merged, `[[`, "summary_all"))
  all_summary <- bind_or_empty(lapply(merged, `[[`, "cliper_summary"))
  all_select <- bind_or_empty(lapply(merged, `[[`, "cliper_select"))
  all_info <- bind_or_empty(lapply(merged, `[[`, "summary_info"))

  if (!is.null(summary_all_csv)) {
    utils::write.csv(all_summary_all, summary_all_csv, row.names = FALSE)
    message("[CLIPER merge] saved: ", summary_all_csv)
  }
  if (!is.null(summary_csv)) {
    utils::write.csv(all_summary, summary_csv, row.names = FALSE)
    message("[CLIPER merge] saved: ", summary_csv)
  }
  if (!is.null(select_csv)) {
    utils::write.csv(all_select, select_csv, row.names = FALSE)
    message("[CLIPER merge] saved: ", select_csv)
  }
  if (!is.null(info_csv)) {
    utils::write.csv(all_info, info_csv, row.names = FALSE)
    message("[CLIPER merge] saved: ", info_csv)
  }

  quit(status = 0)
}

usage()
