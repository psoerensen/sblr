#' Fit Cheng multi-trait BayesC or BayesR with sparse LD
#'
#' Fits the supported independent-provider summary-statistic Cheng
#' MT-BayesC-Pi or factorized pattern-by-scale MT-BayesR model using CSR
#' cross-product operators. Providers may have distinct marker subsets,
#' orders, sample sizes, residual scales, and LD populations. Missing provider
#' markers contribute no likelihood information.
#'
#' @param providers Named provider descriptors containing provider and trait
#'   IDs, operator resource ID, named score, sample size, positive fixed
#'   residual scale, independent-summary regime, and effect scale.
#' @param operator_resources CSR resource descriptors with stable IDs, marker
#'   and allele metadata, coding/scaling metadata, approximation declaration,
#'   provenance, CSR storage, and cross-product diagonal.
#' @param global_marker_ids Ordered unique global marker identifiers.
#' @param global_alleles Global allele metadata aligned to
#'   `global_marker_ids`.
#' @param trait_ids Ordered unique trait identifiers.
#' @param method Either `"bayesc"` or `"bayesr"`.
#' @param initial_marker_covariance Initial finite SPD marker covariance.
#' @param marker_covariance_prior_degrees_of_freedom Inverse-Wishart degrees
#'   of freedom for the marker covariance.
#' @param marker_covariance_prior_scale Finite SPD inverse-Wishart scale matrix.
#' @param initial_activity_pattern_probability Optional initial activity-
#'   pattern simplex.
#' @param activity_pattern_dirichlet_prior Optional positive activity-pattern
#'   Dirichlet prior.
#' @param component_scales Positive, unique, increasing BayesR variance
#'   multipliers. The null state is not included.
#' @param initial_scale_probability Optional conditional simplex over positive
#'   component scales.
#' @param scale_dirichlet_prior Optional positive scale Dirichlet prior.
#' @param marker_multipliers Optional positive marker-specific multipliers.
#' @param nit,nburn,nthin Post-burn iterations, burn-in iterations, and
#'   thinning interval.
#' @param seed Master seed used by the Phase 3 seed contract.
#' @param nchains Number of complete joint chains.
#' @param ncores Requested chain workers.
#' @param chain_seeds Optional exact base seed per chain.
#' @param keep_chains Retain compact chain records in the formatted fit.
#' @param convergence Convergence capture policy, `"core"` or `"none"`.
#' @param convergence_control Reserved; must be `NULL`.
#' @param keep_traces Capture the established unthinned convergence state.
#' @param memory_limit_bytes Incremental fit-memory limit checked before native
#'   sampling.
#' @return An `mtblr_fit`; validated raw-v2 is available as
#'   `attr(fit, "blr_raw")`.
#' @export
mtblr_csr <- function(
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
  .blr_validate_exact_public_call(sys.call(), sys.function(), "mtblr_csr()")
  .mtblr_summary_public_fit(
    providers = providers, operator_resources = operator_resources,
    global_marker_ids = global_marker_ids, global_alleles = global_alleles,
    trait_ids = trait_ids, representation = "csr", interface = "mtblr_csr",
    method = method, initial_marker_covariance = initial_marker_covariance,
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
