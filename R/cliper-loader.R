.cliper_legacy_env <- new.env(parent = globalenv())
.cliper_state <- new.env(parent = emptyenv())
.cliper_state$legacy_loaded <- FALSE

.cliper_legacy_file <- function() {
  legacy_file <- system.file("legacy", "01_Functions.R", package = "CLIPER")
  if (!nzchar(legacy_file)) {
    legacy_file <- file.path(getwd(), "inst", "legacy", "01_Functions.R")
  }
  legacy_file
}

cliper_load <- function(verbose = TRUE) {
  if (isTRUE(.cliper_state$legacy_loaded)) {
    return(invisible(.cliper_legacy_env))
  }

  legacy_file <- .cliper_legacy_file()
  if (!file.exists(legacy_file)) {
    stop("Cannot find bundled CLIPER legacy source: ", legacy_file, call. = FALSE)
  }

  if (isTRUE(verbose)) {
    message("[CLIPER] Loading bundled algorithm code")
  }

  assign(
    "library",
    function(package, ..., character.only = FALSE) {
      pkg <- if (isTRUE(character.only)) package else as.character(substitute(package))
      if (identical(pkg, "ArchR")) {
        return(invisible(TRUE))
      }
      base::library(pkg, ..., character.only = TRUE)
    },
    envir = .cliper_legacy_env
  )

  oldwd <- getwd()
  on.exit(setwd(oldwd), add = TRUE)
  setwd(dirname(legacy_file))
  sys.source(legacy_file, envir = .cliper_legacy_env, keep.source = FALSE)

  .cliper_state$legacy_loaded <- TRUE

  if (isTRUE(verbose)) {
    message("[CLIPER] Ready")
  }

  invisible(.cliper_legacy_env)
}

cliper_functions <- function() {
  cliper_load(verbose = FALSE)
  sort(ls(envir = .cliper_legacy_env, pattern = "^[A-Za-z0-9_.]+$"))
}

.cliper_call <- function(name, ..., .cliper_verbose = FALSE) {
  cliper_load(verbose = .cliper_verbose)
  fn <- get(name, envir = .cliper_legacy_env, inherits = FALSE)
  fn(...)
}

make_metacell_groups <- function(..., .cliper_verbose = FALSE) {
  .cliper_call("make_metacell_groups", ..., .cliper_verbose = .cliper_verbose)
}

normalize_delta <- function(..., .cliper_verbose = FALSE) {
  .cliper_call("normalize_delta", ..., .cliper_verbose = .cliper_verbose)
}

normalize_metacell_deseq2 <- function(..., .cliper_verbose = FALSE) {
  .cliper_call("normalize_metacell_deseq2", ..., .cliper_verbose = .cliper_verbose)
}

normalize_metacell_logcpm <- function(..., .cliper_verbose = FALSE) {
  .cliper_call("normalize_metacell_logcpm", ..., .cliper_verbose = .cliper_verbose)
}

gc_fq_normalize <- function(..., .cliper_verbose = FALSE) {
  .cliper_call("gc_fq_normalize", ..., .cliper_verbose = .cliper_verbose)
}

data_norm <- function(..., .cliper_verbose = FALSE) {
  .cliper_call("data_norm", ..., .cliper_verbose = .cliper_verbose)
}

aggregate_metacell_counts <- function(..., .cliper_verbose = FALSE) {
  .cliper_call("aggregate_metacell_counts", ..., .cliper_verbose = .cliper_verbose)
}

Create_Signac_CLIPER_obj <- function(..., .cliper_verbose = TRUE) {
  if (isTRUE(.cliper_verbose)) message("[CLIPER] Creating Signac CLIPER object")
  out <- .cliper_call("Create_Signac_CLIPER_obj", ..., .cliper_verbose = .cliper_verbose)
  if (isTRUE(.cliper_verbose)) message("[CLIPER] Signac CLIPER object created")
  out
}

make_metacell_groups_from_emb <- function(..., .cliper_verbose = FALSE) {
  .cliper_call("make_metacell_groups_from_emb", ..., .cliper_verbose = .cliper_verbose)
}

parse_peaks <- function(..., .cliper_verbose = FALSE) {
  .cliper_call("parse_peaks", ..., .cliper_verbose = .cliper_verbose)
}

make_gene_window_for_selected <- function(..., .cliper_verbose = FALSE) {
  .cliper_call("make_gene_window_for_selected", ..., .cliper_verbose = .cliper_verbose)
}

get_peaks_by_gene_gr <- function(..., .cliper_verbose = FALSE) {
  .cliper_call("get_peaks_by_gene_gr", ..., .cliper_verbose = .cliper_verbose)
}

CER_PCGS <- function(..., .cliper_verbose = FALSE) {
  .cliper_call("CER_PCGS", ..., .cliper_verbose = .cliper_verbose)
}

Run_CLIPER <- function(..., .cliper_verbose = TRUE) {
  if (isTRUE(.cliper_verbose)) message("[CLIPER] Running peak-to-gene model")
  out <- .cliper_call("Run_CLIPER", ..., .cliper_verbose = .cliper_verbose)
  if (isTRUE(.cliper_verbose)) message("[CLIPER] Run complete")
  out
}

compare_p2g <- function(..., .cliper_verbose = FALSE) {
  .cliper_call("compare_p2g", ..., .cliper_verbose = .cliper_verbose)
}

ids_to_granges_safe <- function(..., .cliper_verbose = FALSE) {
  .cliper_call("ids_to_granges_safe", ..., .cliper_verbose = .cliper_verbose)
}

map_feature_intervals_to_assay_regions <- function(..., .cliper_verbose = FALSE) {
  .cliper_call("map_feature_intervals_to_assay_regions", ..., .cliper_verbose = .cliper_verbose)
}

compute_region_enrichment_cisRest_externalBG <- function(..., .cliper_verbose = FALSE) {
  .cliper_call("compute_region_enrichment_cisRest_externalBG", ..., .cliper_verbose = .cliper_verbose)
}
