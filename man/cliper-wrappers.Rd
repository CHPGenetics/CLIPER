\name{cliper-wrappers}
\alias{cliper_load}
\alias{cliper_functions}
\alias{make_metacell_groups}
\alias{normalize_delta}
\alias{normalize_metacell_deseq2}
\alias{normalize_metacell_logcpm}
\alias{gc_fq_normalize}
\alias{data_norm}
\alias{aggregate_metacell_counts}
\alias{Create_Signac_CLIPER_obj}
\alias{make_metacell_groups_from_emb}
\alias{parse_peaks}
\alias{make_gene_window_for_selected}
\alias{get_peaks_by_gene_gr}
\alias{CER_PCGS}
\alias{Run_CLIPER}
\alias{compare_p2g}
\title{CLIPER wrapper functions}
\usage{
cliper_load(verbose = TRUE)
cliper_functions()
make_metacell_groups(..., .cliper_verbose = FALSE)
normalize_delta(..., .cliper_verbose = FALSE)
normalize_metacell_deseq2(..., .cliper_verbose = FALSE)
normalize_metacell_logcpm(..., .cliper_verbose = FALSE)
gc_fq_normalize(..., .cliper_verbose = FALSE)
data_norm(..., .cliper_verbose = FALSE)
aggregate_metacell_counts(..., .cliper_verbose = FALSE)
Create_Signac_CLIPER_obj(..., .cliper_verbose = TRUE)
make_metacell_groups_from_emb(..., .cliper_verbose = FALSE)
parse_peaks(..., .cliper_verbose = FALSE)
make_gene_window_for_selected(..., .cliper_verbose = FALSE)
get_peaks_by_gene_gr(..., .cliper_verbose = FALSE)
CER_PCGS(..., .cliper_verbose = FALSE)
Run_CLIPER(..., .cliper_verbose = TRUE)
compare_p2g(..., .cliper_verbose = FALSE)
}
\arguments{
\item{...}{Arguments passed through to the bundled CLIPER legacy function.}
\item{.cliper_verbose}{Logical; whether to print package-level progress messages.}
}
\description{
These exported functions forward calls to the bundled CLIPER legacy code.
}
\value{
The return value of the corresponding CLIPER legacy function.
}
