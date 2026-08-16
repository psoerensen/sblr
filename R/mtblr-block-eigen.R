#' Fit Cheng Multi-Trait BayesC-Pi with Block-Eigen Operators
#'
#' Fits the corrected Cheng MT-BayesC-Pi model from independent summary-
#' statistic providers backed by block-eigen cross-product operators. Each
#' provider contributes to exactly one declared trait and may have its own
#' marker subset, order, sample size, residual scale, blocks, and retained
#' eigen rank. Sample overlap is not modelled.
#'
#' Full-rank blocks are exact relative to the declared block-diagonal operator.
#' Retained-rank blocks target the retained reconstruction; omitted cross-block
#' terms are never reconstructed. Individual-level predictions and full
#' residual-covariance draws are unavailable for this likelihood.
#'
#' @inheritParams mtblr_csr
#' @param operator_resources A nonempty list of block-eigen resource
#'   descriptors. Each descriptor must contain `resource_id`, `marker_ids`,
#'   `alleles`, genotype coding and scaling metadata, an approximation
#'   declaration, provenance, and `blocks`. Each block contains `marker_ids`,
#'   `eigenvectors`, and `eigenvalues`.
#' @return An object of class `mtblr_fit`. The validated `blr_raw` version 2
#'   object is available according to the standard output contract.
#' @export
mtblr_block_eigen <- function(
    providers, operator_resources, global_marker_ids, global_alleles,
    trait_ids, method = "bayesc", initial_marker_covariance,
    marker_covariance_prior_degrees_of_freedom,
    marker_covariance_prior_scale,
    initial_activity_pattern_probability = NULL,
    activity_pattern_dirichlet_prior = NULL,
    nit = 1000L, nburn = 500L, nthin = 1L, seed = 1,
    nchains = 1L, ncores = 1L, chain_seeds = NULL,
    keep_chains = FALSE, convergence = c("core", "none"),
    convergence_control = NULL, keep_traces = TRUE,
    memory_limit_bytes = 256 * 1024^2) {
  .blr_validate_exact_public_call(sys.call(), sys.function(),
                                  "mtblr_block_eigen()")
  .mtblr_summary_public_fit(
    providers, operator_resources, global_marker_ids, global_alleles,
    trait_ids, "block_eigen", "mtblr_block_eigen", method,
    initial_marker_covariance,
    marker_covariance_prior_degrees_of_freedom,
    marker_covariance_prior_scale,
    initial_activity_pattern_probability,
    activity_pattern_dirichlet_prior,
    nit, nburn, nthin, seed, nchains, ncores, chain_seeds,
    keep_chains, convergence, convergence_control, keep_traces,
    memory_limit_bytes)
}
