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
    component_scales = NULL, initial_scale_probability = NULL,
    scale_dirichlet_prior = NULL, marker_multipliers = NULL,
    nit = 1000L, nburn = 500L, nthin = 1L, seed = 1,
    nchains = 1L, ncores = 1L, chain_seeds = NULL,
    keep_chains = FALSE, convergence = c("core", "none"),
    convergence_control = NULL, keep_traces = TRUE,
    memory_limit_bytes = 256 * 1024^2) {
  .blr_validate_exact_public_call(sys.call(), sys.function(),
                                  "mtblr_block_eigen()")
  .mtblr_summary_public_fit(
    providers = providers, operator_resources = operator_resources,
    global_marker_ids = global_marker_ids, global_alleles = global_alleles,
    trait_ids = trait_ids, representation = "block_eigen",
    interface = "mtblr_block_eigen", method = method,
    initial_marker_covariance = initial_marker_covariance,
    marker_covariance_prior_degrees_of_freedom =
      marker_covariance_prior_degrees_of_freedom,
    marker_covariance_prior_scale = marker_covariance_prior_scale,
    initial_activity_pattern_probability =
      initial_activity_pattern_probability,
    activity_pattern_dirichlet_prior = activity_pattern_dirichlet_prior,
    component_scales = component_scales,
    initial_scale_probability = initial_scale_probability,
    scale_dirichlet_prior = scale_dirichlet_prior,
    marker_multipliers = marker_multipliers,
    nit = nit, nburn = nburn, nthin = nthin, seed = seed,
    nchains = nchains, ncores = ncores, chain_seeds = chain_seeds,
    keep_chains = keep_chains, convergence = convergence,
    convergence_control = convergence_control, keep_traces = keep_traces,
    memory_limit_bytes = memory_limit_bytes)
}
