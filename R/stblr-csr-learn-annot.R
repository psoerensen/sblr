.format_csr_annot_fit <- function(
  fit,
  nt,
  m,
  annotation_names,
  trait_names = NULL,
  variable_names = NULL
) {
 if (length(fit) < 22) {
  stop(".format_csr_annot_fit() expects the 22-slot annotation CSR return object.")
 }

 names(fit)[1:22] <- c(
  "bm", "dm", "wy", "r", "b", "d", "o",
  "vbs", "vgs", "ves",
  "covb", "covg", "cove",
  "vb", "vg", "ve",
  "pi", "pim", "eta_pi", "eta_vb",
  "vle", "vld"
 )

 if (is.null(trait_names)) trait_names <- paste0("T", seq_len(nt))
 if (is.null(variable_names)) variable_names <- paste0("V", seq_len(m))

 for (i in 1:7) {
  fit[[i]] <- as.matrix(as.data.frame(fit[[i]]))
  rownames(fit[[i]]) <- variable_names
  colnames(fit[[i]]) <- trait_names
 }

 for (i in 8:10) {
  fit[[i]] <- as.matrix(as.data.frame(fit[[i]]))
  rownames(fit[[i]]) <- paste0("Iter", seq_len(nrow(fit[[i]])))
  colnames(fit[[i]]) <- trait_names
 }

 for (i in 11:16) {
  fit[[i]] <- matrix(unlist(fit[[i]]), ncol = nt, byrow = TRUE)
  rownames(fit[[i]]) <- trait_names
  colnames(fit[[i]]) <- trait_names
 }

 fit$pi <- matrix(unlist(fit$pi), ncol = 2, byrow = TRUE)
 rownames(fit$pi) <- trait_names
 colnames(fit$pi) <- c("pi0", "pi1")

 fit$pim <- matrix(unlist(fit$pim), ncol = 2, byrow = TRUE)
 rownames(fit$pim) <- trait_names
 colnames(fit$pim) <- c("pi0", "pi1")

 eta_pi <- matrix(unlist(fit$eta_pi), nrow = nt, byrow = TRUE)
 eta_vb <- matrix(unlist(fit$eta_vb), nrow = nt, byrow = TRUE)
 rownames(eta_pi) <- rownames(eta_vb) <- trait_names
 colnames(eta_pi) <- colnames(eta_vb) <- annotation_names

 vle <- as.matrix(as.data.frame(fit$vle))
 vld <- as.matrix(as.data.frame(fit$vld))
 rownames(vle) <- paste0("Iter", seq_len(nrow(vle)))
 rownames(vld) <- paste0("Iter", seq_len(nrow(vld)))
 colnames(vle) <- trait_names
 colnames(vld) <- trait_names

 out <- fit[1:18]
 out$eta_pi <- eta_pi
 out$eta_vb <- eta_vb
 out$vle <- vle
 out$vld <- vld

 if (sum(diag(out$covb)) > 0) out$rb <- cov2cor(out$covb)
 if (sum(diag(out$covg)) > 0) out$rg <- cov2cor(out$covg)
 if (sum(diag(out$cove)) > 0) out$re <- cov2cor(out$cove)

 out
}

#' Fit ST-BLR with Learned Annotation Effects
#'
#' Fits the CSR ST-BLR sampler while learning annotation effects on
#' marker-specific inclusion probabilities and, optionally, marker-effect
#' variance multipliers.
#'
#' @param stats Sufficient statistics returned by [bed_xtx_xty()].
#' @param ld_prefix Prefix of the disk-backed CSR LD files.
#' @param A An `m` by `K` marker annotation matrix. Rows must correspond to the
#'   markers in `stats` and the sparse LD files.
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
#' @param learn_pi_annot Learn annotation effects on inclusion probabilities.
#' @param learn_vb_annot Learn annotation effects on variance multipliers.
#' @param eta_pi_init,eta_vb_init Initial `K` by `nt` annotation-effect
#'   matrices.
#' @param sigma_eta_pi,sigma_eta_vb Annotation-effect prior scales.
#' @param rw_sd_eta_pi,rw_sd_eta_vb Annotation-effect proposal scales.
#' @param annot_update_every Iterations between annotation-effect updates.
#' @param pi_min,pi_max Bounds for marker-specific inclusion probabilities.
#' @param vb_multiplier_min,vb_multiplier_max Bounds for marker-specific
#'   variance multipliers.
#' @return A formatted ST-BLR fit with learned annotation effects.
#' @export
stblr_csr_learn_annot <- function(
  stats,
  ld_prefix,
  A,
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
  learn_pi_annot = TRUE,
  learn_vb_annot = FALSE,
  eta_pi_init = NULL,
  eta_vb_init = NULL,
  sigma_eta_pi = 1,
  sigma_eta_vb = 1,
  rw_sd_eta_pi = 0.05,
  rw_sd_eta_vb = 0.05,
  annot_update_every = 10,
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
 K <- ncol(A)
 if (K == 0) {
  stop("stblr_csr_learn_annot() requires at least one annotation column in A.")
 }

 if (is.null(eta_pi_init)) eta_pi_init <- matrix(0, nrow = K, ncol = nt)
 if (is.null(eta_vb_init)) eta_vb_init <- matrix(0, nrow = K, ncol = nt)

 eta_pi_init <- as.matrix(eta_pi_init)
 eta_vb_init <- as.matrix(eta_vb_init)
 storage.mode(eta_pi_init) <- "double"
 storage.mode(eta_vb_init) <- "double"

 if (!all(dim(eta_pi_init) == c(K, nt))) {
  stop("eta_pi_init must be ncol(A) x nt.")
 }
 if (!all(dim(eta_vb_init) == c(K, nt))) {
  stop("eta_vb_init must be ncol(A) x nt.")
 }

 raw_fit <- stblr_cpg_omp_csr_annot(
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
  A = A,
  learn_pi_annot = learn_pi_annot,
  learn_vb_annot = learn_vb_annot,
  eta_pi_init = eta_pi_init,
  eta_vb_init = eta_vb_init,
  sigma_eta_pi = sigma_eta_pi,
  sigma_eta_vb = sigma_eta_vb,
  rw_sd_eta_pi = rw_sd_eta_pi,
  rw_sd_eta_vb = rw_sd_eta_vb,
  annot_update_every = as.integer(annot_update_every),
  pi_min = pi_min,
  pi_max = pi_max,
  vb_multiplier_min = vb_multiplier_min,
  vb_multiplier_max = vb_multiplier_max,
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

 fit <- .format_csr_annot_fit(
  fit = raw_fit,
  nt = nt,
  m = m,
  annotation_names = colnames(A),
  trait_names = trait_names,
  variable_names = variable_names
 )

 fit$input <- c(
  list(
   model = "annot",
   n = n,
   m = m,
   nt = nt,
   A = A,
   annotation_names = colnames(A),
   learn_pi_annot = learn_pi_annot,
   learn_vb_annot = learn_vb_annot,
   eta_pi_init = eta_pi_init,
   eta_vb_init = eta_vb_init,
   sigma_eta_pi = sigma_eta_pi,
   sigma_eta_vb = sigma_eta_vb,
   rw_sd_eta_pi = rw_sd_eta_pi,
   rw_sd_eta_vb = rw_sd_eta_vb,
   annot_update_every = annot_update_every,
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
 fit
}
