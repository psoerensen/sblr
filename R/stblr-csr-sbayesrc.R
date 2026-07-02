#' Fit an SBayesRC-Style Annotation-Dependent BayesR Prior
#'
#' Fits a CSR ST-BLR model with an SBayesRC-style annotation-dependent BayesR
#' prior. Marker mixture probabilities are modeled through probit
#' stick-breaking regressions on an annotation matrix. The core prior family is
#' closely aligned with the published SBayesRC prior structure, but the exposed
#' default gamma grid and CSR sparse-LD likelihood differ from the published
#' GCTB SBayesRC implementation.
#'
#' @param stats Sufficient statistics returned by [bed_xtx_xty()].
#' @param ld_prefix Prefix of the disk-backed CSR LD files.
#' @param A An `m` by `K` numeric marker annotation matrix. Rows must correspond
#'   to the markers in `stats` and the sparse LD files.
#' @param Glist Optional genotype/LD metadata. Required when `selection_s` is
#'   non-`NULL` so MAF can be aligned to the CSR LD marker order using
#'   `Glist$rsidsLD`, `Glist$rsids`, and `Glist$maf`.
#' @param n Sample size. Defaults to `stats$n` when available.
#' @param m Number of markers. Inferred from `stats` when omitted.
#' @param gamma Numeric mixture variance multipliers. The first value must be
#'   zero and all remaining values must be positive.
#' @param pi_marker Backward-compatible initial active-marker probability.
#' @param pi_init,pi_vb_init,pi_prior_mean,pi_prior_strength,pi_prior_a,pi_prior_b
#'   Active-marker probability and marker-variance prior controls.
#' @param h2 Initial heritability.
#' @param selection_s Optional fixed global BayesS-style MAF-scaling parameter.
#'   If supplied, non-null component prior variances are scaled by
#'   `h^(selection_s + 1)`, where `h = 2p(1-p)`. `selection_s = NULL` preserves
#'   the ordinary SBayesRC prior. Sampler-level estimation of `selection_s` is
#'   not implemented.
#' @param nub,nue Prior degrees of freedom.
#' @param B,E Optional initial marker-effect and residual covariance matrices.
#' @param ssb_prior,sse_prior Optional prior scale matrices.
#' @param updateAlpha,updateB,updateE Logical sampler update controls.
#' @param adjE Residual adjustment factor.
#' @param nit,nburn,nthin MCMC iteration controls.
#' @param ncores Number of OpenMP threads.
#' @param seed Sampler seed.
#' @param nchains Number of independent MCMC chains.
#' @param chain_seeds Optional integer seeds, one per chain.
#' @param keep_chains Logical; return compact per-chain summaries.
#' @param updateLDswap Logical; attempt optional active/null LD-swap
#'   Metropolis-Hastings moves.
#' @param ld_swap_prob Probability per MCMC iteration of attempting LD-swap
#'   moves when `updateLDswap = TRUE`.
#' @param ld_swap_r2 Minimum LD r-squared for candidate swap partners.
#' @param ld_swap_max_friends Maximum number of high-LD friends stored per
#'   marker for swap proposals.
#' @param ld_swap_moves Number of swap attempts when LD-swap is triggered.
#' @param b_init Optional initial marker effects.
#' @param comp_init Optional initial zero-based mixture-component indices, one
#'   length-`m` vector per trait.
#' @param use_comp_init Use the supplied initial component indices.
#' @param r_init Optional initial residual state.
#' @param use_r_init Use the supplied initial residual state.
#' @param rebuild_r_before_updateE Rebuild the residual before updating `E`.
#' @param add_intercept Add an intercept column to `A` when one is not present.
#' @param standardize_annotations Standardize non-binary annotation columns.
#' @param center_binary_annotations Center and standardize binary annotations.
#' @param active_comp_weights Optional weights distributing `pi_init` across
#'   the active mixture components.
#' @param alpha_init Optional initial annotation coefficient matrix.
#' @param sigmaSqAlpha_init Optional initial annotation-effect variances.
#' @param intercept_flat Use a flat prior for the intercept coefficients.
#' @param sigmaSqAlpha_a,sigmaSqAlpha_b Annotation-effect variance prior
#'   controls.
#' @param pi_floor Lower probability bound used by the sampler.
#' @param alpha_update_every Iterations between annotation-effect updates.
#'
#' @return A formatted SBayesRC-style ST-BLR fit with annotation and
#'   mixture-component posterior summaries, including `comp_prob`,
#'   `dm_component_mean`, `alpha`, `sigmaSqAlpha`, and `ncomp` where available.
#'   The fit includes `vle` and `vld` traces using the same definitions and
#'   formatting conventions as annotation-unaware CSR fits.
#' @export
stblr_csr_sbayesrc_generic <- function(
  stats,
  ld_prefix,
  A,
  n = NULL,
  m = NULL,
  gamma = c(0, 0.01, 0.1, 1),
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
  updateAlpha = TRUE,
  updateB = TRUE,
  updateE = TRUE,
  adjE = 0.9,
  nit = 1000,
  nburn = 100,
  nthin = 1,
  ncores = 3,
  seed = 10,
  nchains = 1L,
  chain_seeds = NULL,
  keep_chains = FALSE,
  updateLDswap = FALSE,
  ld_swap_prob = 0.05,
  ld_swap_r2 = 0.8,
  ld_swap_max_friends = 50L,
  ld_swap_moves = 1L,
  b_init = NULL,
  comp_init = NULL,
  use_comp_init = FALSE,
  r_init = NULL,
  use_r_init = FALSE,
  rebuild_r_before_updateE = FALSE,
  add_intercept = TRUE,
  standardize_annotations = TRUE,
  center_binary_annotations = FALSE,
  active_comp_weights = NULL,
  alpha_init = NULL,
  sigmaSqAlpha_init = NULL,
  intercept_flat = TRUE,
  sigmaSqAlpha_a = 2,
  sigmaSqAlpha_b = 2,
  pi_floor = 1e-12,
  alpha_update_every = 10,
  selection_s = NULL,
  Glist = NULL
) {
 dims <- .stblr_get_nt_m_names(stats, n = n, m = m)
 nt <- dims$nt
 n <- dims$n
 m <- dims$m
 trait_names <- dims$trait_names
 variable_names <- dims$variable_names

 .stblr_validate_stats(stats, nt = nt, m = m)

 if (!is.numeric(nchains) || length(nchains) != 1L ||
     !is.finite(nchains) || nchains < 1 || nchains != floor(nchains)) {
  stop("nchains must be a positive integer scalar.")
 }
 nchains <- as.integer(nchains)
 if (!is.logical(keep_chains) || length(keep_chains) != 1L ||
     is.na(keep_chains)) {
  stop("keep_chains must be TRUE or FALSE.")
 }
 if (!is.null(chain_seeds)) {
  if (!is.numeric(chain_seeds) || length(chain_seeds) != nchains ||
      anyNA(chain_seeds) || any(!is.finite(chain_seeds)) ||
      any(chain_seeds != floor(chain_seeds))) {
   stop("chain_seeds must be NULL or an integer/numeric vector of length nchains.")
  }
  chain_seeds <- as.integer(chain_seeds)
 } else {
  chain_seeds <- integer()
 }
 .validate_ld_swap_args(
  updateLDswap, ld_swap_prob, ld_swap_r2, ld_swap_max_friends, ld_swap_moves
 )
 selection_s_info <- .stblr_prepare_csr_selection_s(
  selection_s = selection_s,
  Glist = Glist,
  m = m,
  scheduled = FALSE,
  backend = "sbayesrc"
 )

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

 state <- .stblr_init_marker_state(nt = nt, m = m, b_init = b_init)
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
 if (ncol(A) == 0) {
  stop("stblr_csr_sbayesrc_generic() requires at least one annotation column in A.")
 }

 alpha <- make_sbayesrc_alpha_init(
  A = A,
  gamma = gamma,
  pi_init = arch$pi_init,
  active_comp_weights = active_comp_weights,
  alpha_init = alpha_init,
  sigmaSqAlpha_init = sigmaSqAlpha_init
 )
 gamma <- alpha$gamma
 Kgamma <- length(gamma)

 if (is.null(comp_init)) {
  comp_init <- lapply(seq_len(nt), function(i) rep(0L, m))
 }
 if (!is.list(comp_init) || length(comp_init) != nt ||
     any(lengths(comp_init) != m)) {
  stop("comp_init must be a list of length nt, each element length m.")
 }
 valid_comp <- vapply(
 comp_init,
  function(x) {
   is.numeric(x) && all(is.finite(x)) && all(x == floor(x)) &&
    all(x >= 0 & x < Kgamma)
  },
  logical(1)
 )
 if (!all(valid_comp)) {
  stop("comp_init values must be finite zero-based component indices in gamma.")
 }
 comp_init <- lapply(comp_init, as.integer)

 for (nm in c("sigmaSqAlpha_a", "sigmaSqAlpha_b")) {
  value <- get(nm)
  if (!is.numeric(value) || length(value) != 1 || !is.finite(value) ||
      value <= 0) {
   stop(nm, " must be a positive finite scalar.")
  }
 }
 if (!is.numeric(pi_floor) || length(pi_floor) != 1 ||
     !is.finite(pi_floor) || pi_floor <= 0 || pi_floor >= 1) {
  stop("pi_floor must be a finite scalar in (0, 1).")
 }
 if (!is.numeric(alpha_update_every) || length(alpha_update_every) != 1 ||
     !is.finite(alpha_update_every) || alpha_update_every <= 0 ||
     alpha_update_every != as.integer(alpha_update_every)) {
  stop("alpha_update_every must be a positive integer.")
 }

 raw_fit <- stblr_cpg_omp_csr_sbayesrc(
  wy = stats$wy,
  ww = stats$ww,
  yy = stats$yy,
  b_init = state$b_init,
  comp_init = comp_init,
  use_comp_init = use_comp_init,
  r_init = r_init,
  use_r_init = use_r_init,
  rebuild_r_before_updateE = rebuild_r_before_updateE,
  ld_prefix = ld_prefix,
  B = pri$B,
  E = pri$E,
  ssb_prior = pri$ssb_prior_list,
  sse_prior = pri$sse_prior_list,
  A = A,
  gamma = gamma,
  alpha_init = alpha$alpha_init,
  sigmaSqAlpha_init = alpha$sigmaSqAlpha_init,
  intercept_flat = intercept_flat,
  sigmaSqAlpha_a = sigmaSqAlpha_a,
  sigmaSqAlpha_b = sigmaSqAlpha_b,
  pi_floor = pi_floor,
  nub = nub,
  nue = nue,
  updateAlpha = updateAlpha,
  updateB = updateB,
  updateE = updateE,
  alpha_update_every = as.integer(alpha_update_every),
  adjE = adjE,
  n = rep(as.integer(n), nt),
  nit = as.integer(nit),
  nburn = as.integer(nburn),
  nthin = as.integer(nthin),
  ncores = as.integer(ncores),
  seed = as.integer(seed),
  nchains = nchains,
  keep_chains = keep_chains,
  chain_seeds = chain_seeds,
  updateLDswap = updateLDswap,
  ld_swap_prob = ld_swap_prob,
  ld_swap_r2 = ld_swap_r2,
  ld_swap_max_friends = as.integer(ld_swap_max_friends),
  ld_swap_moves = as.integer(ld_swap_moves),
  selection_s_prior_scale = selection_s_info$prior_scale
 )

 fit <- format_sbayesrc_csr_fit(
  fit = raw_fit,
  nt = nt,
  m = m,
  gamma = gamma,
  n_anno = ncol(A),
  trait_names = trait_names,
  variable_names = variable_names,
  annotation_names = colnames(A),
  nchains = nchains,
  keep_chains = keep_chains
 )

 fit$input <- c(
  list(
   model = "sbayesrc",
   n = n,
   m = m,
   nt = nt,
   A = A,
   annotation_names = colnames(A),
   gamma = gamma,
   component_names = colnames(fit$ncomp),
   active_comp_weights = alpha$active_comp_weights,
   component_prob_init = alpha$component_prob_init,
   step_prob_init = alpha$step_prob_init,
   alpha_init = alpha$alpha_init,
   sigmaSqAlpha_init = alpha$sigmaSqAlpha_init,
   comp_init = comp_init,
   use_comp_init = use_comp_init,
   intercept_flat = intercept_flat,
   sigmaSqAlpha_a = sigmaSqAlpha_a,
   sigmaSqAlpha_b = sigmaSqAlpha_b,
   pi_floor = pi_floor,
   updateAlpha = updateAlpha,
   alpha_update_every = alpha_update_every,
   add_intercept = add_intercept,
   standardize_annotations = standardize_annotations,
   center_binary_annotations = center_binary_annotations,
   h2 = h2,
   selection_s = selection_s_info$selection_s,
   selection_s_fixed = selection_s_info$fixed,
   selection_s_exponent = selection_s_info$exponent,
   selection_s_scale = "standardized_genotype_effect",
   nub = nub,
   nue = nue,
   vy = pri$vy,
   B = pri$B,
   E = pri$E,
   ssb_prior = pri$ssb_prior,
   sse_prior = pri$sse_prior,
   updateB = updateB,
   updateE = updateE,
   adjE = adjE,
   nit = nit,
   nburn = nburn,
   nthin = nthin,
   ncores = ncores,
   seed = seed,
   nchains = nchains,
   keep_chains = keep_chains,
   chain_seeds = if (length(chain_seeds)) chain_seeds else NULL,
   updateLDswap = updateLDswap,
   ld_swap_prob = ld_swap_prob,
   ld_swap_r2 = ld_swap_r2,
   ld_swap_max_friends = as.integer(ld_swap_max_friends),
   ld_swap_moves = as.integer(ld_swap_moves),
   chain_seed_rule = if (length(chain_seeds)) {
    "chain_seeds[chain] + 1000003 * (trait + 1)"
   } else if (nchains == 1L) {
    "seed + 1000003 * (trait + 1)"
   } else {
    "seed + 1000003 * (trait + 1) + 9176 * (chain + 1)"
   },
   use_r_init = use_r_init,
   rebuild_r_before_updateE = rebuild_r_before_updateE,
   ld_prefix = ld_prefix
  ),
  arch
 )

 .standardize_stblr_annotation_fit(fit, "sbayesrc")
}
