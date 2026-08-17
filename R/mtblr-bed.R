#' Fit Cheng multi-trait BayesC or BayesR from packed BED genotypes
#'
#' Fits the supported common-sample Cheng MT-BayesC-Pi or factorized
#' pattern-by-scale MT-BayesR model. Phenotypes must be complete and centred
#' or residualized as appropriate for the analysis. Independent-provider
#' summary-statistic analyses use [mtblr_csr()] or [mtblr_block_eigen()].
#'
#' @param y Complete numeric phenotype matrix with one column per trait.
#' @param Glist One packed-BED genotype list shared by all traits.
#' @param trait_ids Unique stable trait IDs in phenotype-column order.
#' @param sample_ids Unique stable sample IDs in phenotype-row order.
#' @param method Either `"bayesc"` or `"bayesr"`.
#' @param residual_covariance_policy Either `"fixed_full"` or
#'   `"sampled_full"`.
#' @param fixed_residual_covariance Fixed full residual covariance for
#'   `fixed_full`; otherwise `NULL`.
#' @param initial_marker_covariance Initial full marker-effect covariance.
#' @param marker_covariance_prior_degrees_of_freedom Inverse-Wishart degrees
#'   of freedom for the marker covariance.
#' @param marker_covariance_prior_scale Inverse-Wishart scale matrix for the
#'   marker covariance.
#' @param initial_activity_pattern_probability Optional initial simplex over
#'   the canonical complete binary activity-pattern order.
#' @param activity_pattern_dirichlet_prior Optional positive Dirichlet prior
#'   in the same activity-pattern order.
#' @param component_scales,initial_scale_probability,scale_dirichlet_prior,marker_multipliers
#'   BayesR positive variance multipliers, their optional initial and
#'   Dirichlet probability vectors, and positive marker-specific multipliers.
#'   These arguments must be `NULL` for `method = "bayesc"`.
#' @param initial_residual_covariance Initial full residual covariance for
#'   `sampled_full`; otherwise `NULL`.
#' @param residual_covariance_prior_degrees_of_freedom Inverse-Wishart degrees
#'   of freedom for sampled residual covariance.
#' @param residual_covariance_prior_scale Inverse-Wishart scale matrix for
#'   sampled residual covariance.
#' @param chr,cls,rows,block_size Packed-BED marker and sample selection
#'   controls, following [stblr_bed()].
#' @param nit,nburn,nthin Post-burn iterations, burn-in iterations, and
#'   thinning interval.
#' @param seed Master seed resolved through seed-contract version 1.
#' @param nchains Number of complete joint chains.
#' @param ncores Requested chain workers.
#' @param chain_seeds Optional signed-int32 or uint32-compatible base seed per
#'   chain.
#' @param keep_chains Retain compact chain-final records in the formatted fit.
#' @param convergence Either `"core"` or `"none"`.
#' @param convergence_control Reserved; must be `NULL`.
#' @param keep_traces Retain unthinned observational convergence traces.
#' @param memory_limit_bytes Incremental fit-allocation limit checked before
#'   provider construction and native sampling. `NULL` or positive `Inf`
#'   disables the finite limit.
#' @return An `mtblr_fit` and `blr_fit`. The validated `blr_raw` version 2
#'   object is available as `attr(fit, "blr_raw")`.
#' @export
mtblr_bed <- function(
    y, Glist, trait_ids = colnames(as.matrix(y)),
    sample_ids = rownames(as.matrix(y)), method = "bayesc",
    residual_covariance_policy = "fixed_full",
    fixed_residual_covariance = NULL,
    initial_marker_covariance,
    marker_covariance_prior_degrees_of_freedom,
    marker_covariance_prior_scale,
    initial_activity_pattern_probability = NULL,
    activity_pattern_dirichlet_prior = NULL,
    component_scales = NULL, initial_scale_probability = NULL,
    scale_dirichlet_prior = NULL, marker_multipliers = NULL,
    initial_residual_covariance = NULL,
    residual_covariance_prior_degrees_of_freedom = NULL,
    residual_covariance_prior_scale = NULL,
    chr = NULL, cls = NULL, rows = NULL, block_size = 1000L,
    nit = 1000L, nburn = 500L, nthin = 1L, seed = 1,
    nchains = 1L, ncores = 1L, chain_seeds = NULL,
    keep_chains = FALSE, convergence = c("core", "none"),
    convergence_control = NULL, keep_traces = TRUE,
    memory_limit_bytes = 256 * 1024^2) {
  .blr_validate_exact_public_call(sys.call(), sys.function(), "mtblr_bed()")
  method <- .blr_character_scalar(method, "method", c("bayesc", "bayesr"))
  residual_covariance_policy <- .blr_character_scalar(
    residual_covariance_policy, "residual_covariance_policy",
    c("fixed_full", "sampled_full"))
  controls <- .blr_joint_mt_controls(
    nit, nburn, nthin, seed, nchains, ncores, chain_seeds,
    keep_chains, convergence, convergence_control, keep_traces,
    "mtblr_bed()")
  phenotype <- as.matrix(y)
  if (!is.numeric(phenotype) || nrow(phenotype) <= 1L ||
      ncol(phenotype) < 2L || any(!is.finite(phenotype))) {
    stop("y must be a complete finite N x T numeric matrix with T >= 2.",
         call. = FALSE)
  }
  trait_ids <- .blr_ids(trait_ids, "trait_ids")
  if (length(trait_ids) != ncol(phenotype) ||
      (!is.null(colnames(phenotype)) &&
       !identical(colnames(phenotype), trait_ids))) {
    stop("trait_ids must exactly match phenotype column order.", call. = FALSE)
  }
  sample_ids <- .blr_ids(sample_ids, "sample_ids")
  if (length(sample_ids) != nrow(phenotype) ||
      (!is.null(rownames(phenotype)) &&
       !identical(rownames(phenotype), sample_ids))) {
    stop("sample_ids must exactly match phenotype row order.", call. = FALSE)
  }
  dimnames(phenotype) <- list(sample_ids, trait_ids)
  raw <- .blr_cheng_mt_bayesc_bed_qualification(
    y = phenotype, Glist = Glist,
    fixed_residual_covariance = fixed_residual_covariance,
    initial_marker_covariance = initial_marker_covariance,
    marker_covariance_prior_df =
      marker_covariance_prior_degrees_of_freedom,
    marker_covariance_prior_scale = marker_covariance_prior_scale,
    initial_activity_pattern_probability =
      initial_activity_pattern_probability,
    activity_pattern_dirichlet_prior = activity_pattern_dirichlet_prior,
    method = method, component_scales = component_scales,
    initial_scale_probability = initial_scale_probability,
    scale_dirichlet_prior = scale_dirichlet_prior,
    marker_multipliers = marker_multipliers,
    update_marker_covariance = TRUE,
    update_activity_pattern_probability = TRUE,
    burn_in_iterations = controls$chain$nburn,
    sampling_iterations = controls$chain$nit,
    thin_interval = controls$chain$nthin,
    chains = controls$chain$nchains, cores = controls$chain$ncores,
    seed = controls$chain$seed,
    chain_seeds = controls$chain$chain_seeds_requested,
    keep_traces = controls$keep_traces &&
      identical(controls$convergence, "core"),
    chr = chr, cls = cls, rows = rows, block_size = block_size,
    residual_covariance_policy = residual_covariance_policy,
    initial_residual_covariance = initial_residual_covariance,
    residual_covariance_prior_df =
      residual_covariance_prior_degrees_of_freedom,
    residual_covariance_prior_scale = residual_covariance_prior_scale,
    memory_limit_bytes = memory_limit_bytes)
  raw <- .blr_phase5b_promote_raw(raw)
  .blr_format_cheng_mt_raw_v2(
    raw, keep_chains = controls$chain$keep_chains)
}
