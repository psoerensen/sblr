#' Fit ST-BLR with Fixed Marker-Specific Priors
#'
#' Fits the CSR ST-BLR sampler using fixed marker-specific inclusion
#' probabilities and/or marker-effect variance multipliers. Marker-specific
#' priors can be supplied directly or derived from a marker annotation matrix
#' and fixed annotation coefficients.
#'
#' @param stats Sufficient statistics returned by [bed_xtx_xty()].
#' @param ld_prefix Prefix of the disk-backed CSR LD files.
#' @param A Optional `m` by `K` marker annotation matrix. Rows must correspond
#'   to the markers in `stats` and the sparse LD files.
#' @param n Sample size. Defaults to `stats$n` when available.
#' @param m Number of markers. Inferred from `stats` when omitted.
#' @param pi_marker Backward-compatible initial inclusion probability.
#' @param pi_init,pi_vb_init,pi_prior_mean,pi_prior_strength,pi_prior_a,pi_prior_b
#'   Inclusion-probability and marker-variance prior controls.
#' @param h2 Initial heritability.
#' @param nub,nue Prior degrees of freedom.
#' @param B,E Optional initial marker-effect and residual covariance matrices.
#' @param ssb_prior,sse_prior Optional prior scale matrices.
#' @param updateB,updateE,updatePi Logical sampler update controls.
#' @param adjE Residual adjustment factor.
#' @param nit,nburn,nthin MCMC iteration controls.
#' @param ncores Number of OpenMP threads.
#' @param seed Sampler seed.
#' @param b_init,d_init Optional initial marker effects and inclusion states.
#' @param use_d_init Use the supplied initial inclusion states.
#' @param r_init Optional initial residual state.
#' @param use_r_init Use the supplied initial residual state.
#' @param rebuild_r_before_updateE Rebuild the residual before updating `E`.
#' @param add_intercept Add an intercept column to `A`.
#' @param standardize_annotations Standardize non-binary annotation columns.
#' @param center_binary_annotations Center and standardize binary annotations.
#' @param use_pi_marker Use marker-specific inclusion probabilities.
#' @param use_vb_multiplier Use marker-specific variance multipliers.
#' @param fixed_pi_marker Optional list of fixed marker-specific inclusion
#'   probabilities, one length-`m` vector per trait.
#' @param fixed_vb_multiplier Optional list of fixed marker-specific variance
#'   multipliers, one length-`m` vector per trait.
#' @param beta_pi,beta_vb Optional fixed annotation coefficients used to derive
#'   marker-specific inclusion probabilities and variance multipliers from `A`.
#' @param pi_min,pi_max Bounds for annotation-derived inclusion probabilities.
#' @param vb_multiplier_min,vb_multiplier_max Bounds for annotation-derived
#'   variance multipliers.
#' @return A formatted ST-BLR fit with the resolved fixed priors in `input`.
#' @export
stblr_csr_prior_annot <- function(
  stats,
  ld_prefix,
  A = NULL,
  n = NULL,
  m = NULL,
  pi_marker = 0.001,
  pi_init = NULL,
  pi_vb_init = NULL,
  pi_prior_mean = NULL,
  pi_prior_strength = NULL,
  pi_prior_a = NULL,
  pi_prior_b = NULL,
  h2 = 0.5,
  nub = 4,
  nue = 4,
  B = NULL,
  E = NULL,
  ssb_prior = NULL,
  sse_prior = NULL,
  updateB = TRUE,
  updateE = TRUE,
  updatePi = TRUE,
  adjE = 0.9,
  nit = 1000,
  nburn = 100,
  nthin = 1,
  ncores = 3,
  seed = 10,
  b_init = NULL,
  d_init = NULL,
  use_d_init = FALSE,
  r_init = NULL,
  use_r_init = FALSE,
  rebuild_r_before_updateE = FALSE,
  add_intercept = FALSE,
  standardize_annotations = TRUE,
  center_binary_annotations = FALSE,
  use_pi_marker = FALSE,
  use_vb_multiplier = FALSE,
  fixed_pi_marker = NULL,
  fixed_vb_multiplier = NULL,
  beta_pi = NULL,
  beta_vb = NULL,
  pi_min = 1e-8,
  pi_max = 0.5,
  vb_multiplier_min = 1e-3,
  vb_multiplier_max = 1e3
) {
 dims <- .stblr_get_nt_m_names(stats, n = n, m = m)
 nt <- dims$nt
 n <- dims$n
 m <- dims$m
 trait_names <- dims$trait_names
 variable_names <- dims$variable_names

 .stblr_validate_stats(stats, nt = nt, m = m)

 arch <- .stblr_resolve_architecture(
  pi_marker = pi_marker,
  pi_init = pi_init,
  pi_vb_init = pi_vb_init,
  pi_prior_mean = pi_prior_mean,
  pi_prior_strength = pi_prior_strength,
  pi_prior_a = pi_prior_a,
  pi_prior_b = pi_prior_b
 )

 pri <- .stblr_make_csr_variance_priors(
  stats = stats,
  n = n,
  m = m,
  nt = nt,
  h2 = h2,
  nub = nub,
  nue = nue,
  pi_vb_init = arch$pi_vb_init,
  pi_prior_mean = arch$pi_prior_mean,
  trait_names = trait_names,
  B = B,
  E = E,
  ssb_prior = ssb_prior,
  sse_prior = sse_prior
 )

 state <- .stblr_init_marker_state(
  nt = nt,
  m = m,
  b_init = b_init,
  d_init = d_init
 )
 r_init <- .stblr_init_r_state(
  stats,
  nt = nt,
  m = m,
  use_r_init = use_r_init,
  r_init = r_init
 )

 A <- .stblr_prepare_annotation_matrix(
  A = A,
  m = m,
  variable_names = variable_names,
  add_intercept = add_intercept,
  standardize = standardize_annotations,
  center_binary = center_binary_annotations
 )

 if (!is.null(fixed_pi_marker)) {
  use_pi_marker <- TRUE
  pi_marker_list <- fixed_pi_marker
 } else if (use_pi_marker || !is.null(beta_pi)) {
  ann_prior <- .stblr_make_prior_from_annotations(
   A = A,
   nt = nt,
   pi_base = arch$pi_init,
   beta_pi = beta_pi,
   beta_vb = beta_vb,
   pi_min = pi_min,
   pi_max = pi_max,
   vb_multiplier_min = vb_multiplier_min,
   vb_multiplier_max = vb_multiplier_max
  )
  pi_marker_list <- ann_prior$pi_marker
 } else {
  pi_marker_list <- lapply(seq_len(nt), function(i) rep(arch$pi_init, m))
 }

 if (!is.null(fixed_vb_multiplier)) {
  use_vb_multiplier <- TRUE
  vb_multiplier_list <- fixed_vb_multiplier
 } else if (use_vb_multiplier || !is.null(beta_vb)) {
  if (!exists("ann_prior", inherits = FALSE)) {
   ann_prior <- .stblr_make_prior_from_annotations(
    A = A,
    nt = nt,
    pi_base = arch$pi_init,
    beta_pi = beta_pi,
    beta_vb = beta_vb,
    pi_min = pi_min,
    pi_max = pi_max,
    vb_multiplier_min = vb_multiplier_min,
    vb_multiplier_max = vb_multiplier_max
   )
  }
  vb_multiplier_list <- ann_prior$vb_multiplier
 } else {
  vb_multiplier_list <- lapply(seq_len(nt), function(i) rep(1, m))
 }

 raw_fit <- stblr_cpg_omp_csr_prior(
  wy = stats$wy,
  ww = stats$ww,
  yy = stats$yy,
  b_init = state$b_init,
  d_init = state$d_init,
  use_d_init = use_d_init,
  r_init = r_init,
  use_r_init = use_r_init,
  rebuild_r_before_updateE = rebuild_r_before_updateE,
  ld_prefix = ld_prefix,
  B = pri$B,
  E = pri$E,
  ssb_prior = pri$ssb_prior_list,
  sse_prior = pri$sse_prior_list,
  pi = arch$pi,
  use_pi_marker = use_pi_marker,
  pi_marker = pi_marker_list,
  use_vb_multiplier = use_vb_multiplier,
  vb_multiplier = vb_multiplier_list,
  nub = nub,
  nue = nue,
  updateB = updateB,
  updateE = updateE,
  updatePi = updatePi,
  adjE = adjE,
  n = rep(as.integer(n), nt),
  nit = as.integer(nit),
  nburn = as.integer(nburn),
  nthin = as.integer(nthin),
  pi_prior_a = arch$pi_prior_a,
  pi_prior_b = arch$pi_prior_b,
  ncores = as.integer(ncores),
  seed = as.integer(seed)
 )

 fit <- .format_stblr_fit(
  raw_fit,
  nt = nt,
  m = m,
  trait_names = trait_names,
  variable_names = variable_names
 )

 fit$input <- c(
  list(
   model = "prior",
   n = n,
   m = m,
   nt = nt,
   A = A,
   annotation_names = colnames(A),
   use_pi_marker = use_pi_marker,
   use_vb_multiplier = use_vb_multiplier,
   pi_marker = if (use_pi_marker) pi_marker_list else NULL,
   vb_multiplier = if (use_vb_multiplier) vb_multiplier_list else NULL,
   h2 = h2,
   nub = nub,
   nue = nue,
   vy = pri$vy,
   B = pri$B,
   E = pri$E,
   ssb_prior = pri$ssb_prior,
   sse_prior = pri$sse_prior,
   updateB = updateB,
   updateE = updateE,
   updatePi = updatePi,
   adjE = adjE,
   nit = nit,
   nburn = nburn,
   nthin = nthin,
   ncores = ncores,
   seed = seed,
   use_d_init = use_d_init,
   use_r_init = use_r_init,
   rebuild_r_before_updateE = rebuild_r_before_updateE,
   ld_prefix = ld_prefix
  ),
  arch
 )
 .standardize_stblr_annotation_fit(fit, "prior")
}
