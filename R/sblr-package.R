#' sblr: Scalable Bayesian Linear Regression
#'
#' Experimental STBLR and MTBLR models for individual-level packed-BED data
#' and summary statistics represented by CSR sparse-LD or block-eigen
#' operators. Canonical prior families are BayesC, BayesR, and
#' annotation-informed BayesRC. See [stblr_csr()], [stblr_csr_annot()],
#' [stblr_block_eigen()], [stblr_bed()], [mtblr_csr()],
#' [mtblr_block_eigen()], and [mtblr_bed()].
#'
#' The public `s` model prefix denotes summary-statistics data and is
#' independent of the optional MAF-dependent `maf_effect_s` effect scale.
#'
#' @useDynLib sblr, .registration = TRUE
#' @importFrom Rcpp evalCpp
"_PACKAGE"
