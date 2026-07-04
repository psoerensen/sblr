.format_csr_group_annot_fit <- function(
  fit,
  nt,
  m,
  ngroup,
  group_names = NULL,
  trait_names = NULL,
  variable_names = NULL
) {
 if (length(fit) < 26) {
  stop(".format_csr_group_annot_fit() expects the 26-slot group CSR return object.")
 }

 names(fit)[1:26] <- c(
  "bm", "dm", "wy", "r", "b", "d", "o",
  "vbs", "vgs", "ves",
  "covb", "covg", "cove",
  "vb", "vg", "ve",
  "pi", "pim", "diagnostics", "pimarker",
  "vle", "vld",
  "group_pi", "group_vb_multiplier", "group_nincluded", "group_size"
 )

 if (is.null(trait_names)) trait_names <- paste0("T", seq_len(nt))
 if (is.null(variable_names)) variable_names <- paste0("V", seq_len(m))
 if (is.null(group_names)) group_names <- paste0("G", seq_len(ngroup))

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

 fit$diagnostics <- matrix(unlist(fit$diagnostics), ncol = 4, byrow = TRUE)
 rownames(fit$diagnostics) <- trait_names
 colnames(fit$diagnostics) <- c(
  "log_cpo", "mean_log_cpo", "seconds_mean", "seconds_max"
 )

 fit$pimarker <- matrix(unlist(fit$pimarker), ncol = 2, byrow = TRUE)
 rownames(fit$pimarker) <- trait_names
 colnames(fit$pimarker) <- c("nsamples", "n")

 vle <- as.matrix(as.data.frame(fit$vle))
 vld <- as.matrix(as.data.frame(fit$vld))
 rownames(vle) <- paste0("Iter", seq_len(nrow(vle)))
 rownames(vld) <- paste0("Iter", seq_len(nrow(vld)))
 colnames(vle) <- trait_names
 colnames(vld) <- trait_names

 group_pi <- matrix(
  unlist(fit$group_pi), nrow = nt, ncol = ngroup, byrow = TRUE
 )
 group_vb_multiplier <- matrix(
  unlist(fit$group_vb_multiplier), nrow = nt, ncol = ngroup, byrow = TRUE
 )
 group_nincluded <- matrix(
  unlist(fit$group_nincluded), nrow = nt, ncol = ngroup, byrow = TRUE
 )
 group_size <- matrix(
  unlist(fit$group_size), nrow = nt, ncol = ngroup, byrow = TRUE
 )

 rownames(group_pi) <- rownames(group_vb_multiplier) <-
  rownames(group_nincluded) <- rownames(group_size) <- trait_names
 colnames(group_pi) <- colnames(group_vb_multiplier) <-
  colnames(group_nincluded) <- colnames(group_size) <- group_names

 out <- fit[1:20]
 out$vle <- vle
 out$vld <- vld
 out$group_pi <- group_pi
 out$group_vb_multiplier <- group_vb_multiplier
 out$group_nincluded <- group_nincluded
 out$group_size <- group_size
 out$log_cpo <- out$diagnostics[, "log_cpo"]
 out$mean_log_cpo <- out$diagnostics[, "mean_log_cpo"]

 has_ld_swap <- length(fit) >= 27L &&
  length(fit[[27L]]) == nt &&
  all(vapply(fit[[27L]], function(x) {
   length(x) >= 4L && isTRUE(all.equal(as.numeric(x[[4L]]), 1))
  }, logical(1)))
 if (has_ld_swap) {
  names(fit)[27L] <- "ld_swap_raw"
  ld_swap <- matrix(unlist(fit$ld_swap_raw, use.names = FALSE), ncol = 4, byrow = TRUE)
  ld_swap <- ld_swap[, 1:3, drop = FALSE]
  rownames(ld_swap) <- trait_names
  colnames(ld_swap) <- c("attempted", "accepted", "acceptance_rate")
  out$ld_swap <- as.data.frame(ld_swap)
 }

 summary_start <- if (has_ld_swap) 28L else 27L
 summary_end <- summary_start + 5L
 if (length(fit) >= summary_end) {
  names(fit)[summary_start:summary_end] <- c(
   "bm_sd", "bm_min", "bm_max", "dm_sd", "dm_min", "dm_max"
  )
  for (nm in c("bm_sd", "bm_min", "bm_max", "dm_sd", "dm_min", "dm_max")) {
   out[[nm]] <- as.matrix(as.data.frame(fit[[nm]]))
   rownames(out[[nm]]) <- variable_names
   colnames(out[[nm]]) <- trait_names
  }
 }

 chain_start <- summary_end + 1L
 if (length(fit) >= (chain_start + 5L)) {
  names(fit)[chain_start:(chain_start + 5L)] <- c(
   "chain_dm_raw", "chain_bm_raw", "chain_ld_swap_raw", "chain_group_pi_raw",
   "chain_group_vb_raw", "chain_group_nincluded_raw"
  )
  chains <- setNames(vector("list", nt), trait_names)
  ld_swap_chains <- list()
  for (tt in seq_len(nt)) {
   dm_raw <- as.numeric(fit$chain_dm_raw[[tt]])
   bm_raw <- as.numeric(fit$chain_bm_raw[[tt]])
   diag_raw <- as.numeric(fit$chain_ld_swap_raw[[tt]])
   nchains <- if (m > 0L) length(dm_raw) %/% m else 0L
   trait_chains <- setNames(vector("list", nchains), paste0("chain", seq_len(nchains)))
   trait_ld <- matrix(NA_real_, nrow = nchains, ncol = 3)
   colnames(trait_ld) <- c("attempted", "accepted", "acceptance_rate")
   rownames(trait_ld) <- names(trait_chains)
   for (cc in seq_len(nchains)) {
    idx <- ((cc - 1L) * m + 1L):(cc * m)
    gidx <- ((cc - 1L) * ngroup + 1L):(cc * ngroup)
    trait_chains[[cc]] <- list(
     dm = stats::setNames(dm_raw[idx], variable_names),
     bm = stats::setNames(bm_raw[idx], variable_names),
     group_pi = stats::setNames(as.numeric(fit$chain_group_pi_raw[[tt]])[gidx], group_names),
     group_vb_multiplier = stats::setNames(as.numeric(fit$chain_group_vb_raw[[tt]])[gidx], group_names),
     group_nincluded = stats::setNames(as.numeric(fit$chain_group_nincluded_raw[[tt]])[gidx], group_names)
    )
    didx <- ((cc - 1L) * 4L + 1L):((cc - 1L) * 4L + 4L)
    if (length(diag_raw) >= max(didx)) {
     trait_chains[[cc]]$ld_swap <- stats::setNames(
      diag_raw[didx][1:3], c("attempted", "accepted", "acceptance_rate")
     )
     trait_ld[cc, ] <- diag_raw[didx][1:3]
    }
   }
   chains[[tt]] <- trait_chains
   ld_swap_chains[[trait_names[tt]]] <- as.data.frame(trait_ld)
  }
  out$chains <- chains
  out$ld_swap_chains <- ld_swap_chains
 }

 if (sum(diag(out$covb)) > 0) out$rb <- cov2cor(out$covb)
 if (sum(diag(out$covg)) > 0) out$rg <- cov2cor(out$covg)
 if (sum(diag(out$cove)) > 0) out$re <- cov2cor(out$cove)

 .stblr_ensure_ld_swap_fields(out)
}

#' Fit ST-BLR with Group-Level Annotation Priors
#'
#' Fits the CSR ST-BLR sampler using group-level inclusion probabilities and
#' marker-effect variance multipliers.
#'
#' @param stats Sufficient statistics returned by [bed_xtx_xty()].
#' @param ld_prefix Prefix of the disk-backed CSR LD files.
#' @param group A length-`m` marker group vector. Named vectors are reordered
#'   to match marker names when possible. Defaults to one group.
#' @param n Sample size. Defaults to `stats$n` when available.
#' @param m Number of markers. Inferred from `stats` when omitted.
#' @param group_names Optional group ordering.
#' @param pi_marker Backward-compatible initial inclusion probability.
#' @param pi_init,pi_vb_init,pi_prior_mean,pi_prior_strength,pi_prior_a,pi_prior_b
#'   Global inclusion-probability and marker-variance prior controls.
#' @param h2 Initial heritability.
#' @param nub,nue Prior degrees of freedom.
#' @param B,E Optional initial marker-effect and residual covariance matrices.
#' @param ssb_prior,sse_prior Optional prior scale matrices.
#' @param updateB,updateE,updatePi Logical sampler update controls.
#' @param adjE Residual adjustment factor.
#' @param nit,nburn,nthin MCMC iteration controls.
#' @param ncores Number of OpenMP threads.
#' @param seed Sampler seed.
#' @param nchains Number of independent MCMC chains.
#' @param keep_chains Logical; return compact per-chain summaries.
#' @param chain_seeds Optional numeric or integer vector of length `nchains`.
#' @param updateLDswap Logical; attempt active/null LD-swap
#'   Metropolis-Hastings moves. Group-specific inclusion probabilities and
#'   variance multipliers are included in the MH prior ratio.
#' @param ld_swap_prob Probability per MCMC iteration of attempting LD-swap
#'   moves when `updateLDswap = TRUE`.
#' @param ld_swap_r2 Minimum LD r-squared for candidate swap partners.
#' @param ld_swap_max_friends Maximum number of high-LD friends stored per
#'   marker for swap proposals.
#' @param ld_swap_moves Number of swap attempts when LD-swap is triggered.
#' @param b_init,d_init Optional initial marker effects and inclusion states.
#' @param use_d_init Use the supplied initial inclusion states.
#' @param r_init Optional initial residual state.
#' @param use_r_init Use the supplied initial residual state.
#' @param rebuild_r_before_updateE Rebuild the residual before updating `E`.
#' @param group_pi_init Initial group inclusion probabilities.
#' @param group_vb_multiplier_init Initial group variance multipliers.
#' @param pi_group_prior_mean,pi_group_prior_strength Group beta-prior controls.
#' @param pi_group_prior_a,pi_group_prior_b Optional group beta-prior shapes.
#' @param updateGroupVb Update group variance multipliers.
#' @param nub_group Group variance-prior degrees of freedom.
#' @param ssb_group_prior Group variance-prior scale.
#' @param normalize_group_vb Normalize group variance multipliers.
#' @return A formatted ST-BLR fit with group-level posterior summaries. The
#'   fit includes `vle` and `vld` traces using the same definitions and
#'   formatting conventions as annotation-unaware CSR fits.
#' @export
stblr_csr_group_annot <- function(
  stats,
  ld_prefix,
  group = NULL,
  n = NULL,
  m = NULL,
  group_names = NULL,
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
  nchains = 1L,
  keep_chains = FALSE,
  chain_seeds = NULL,
  updateLDswap = FALSE,
  ld_swap_prob = 0.05,
  ld_swap_r2 = 0.8,
  ld_swap_max_friends = 50L,
  ld_swap_moves = 1L,
  b_init = NULL,
  d_init = NULL,
  use_d_init = FALSE,
  r_init = NULL,
  use_r_init = FALSE,
  rebuild_r_before_updateE = FALSE,
  group_pi_init = NULL,
  group_vb_multiplier_init = NULL,
  pi_group_prior_mean = NULL,
  pi_group_prior_strength = NULL,
  pi_group_prior_a = NULL,
  pi_group_prior_b = NULL,
  updateGroupVb = FALSE,
  nub_group = 4,
  ssb_group_prior = 1,
  normalize_group_vb = TRUE
) {
 .validate_ld_swap_args(
  updateLDswap, ld_swap_prob, ld_swap_r2, ld_swap_max_friends, ld_swap_moves
 )
 chain_args <- .stblr_validate_annotation_chain_args(
  nchains = nchains,
  keep_chains = keep_chains,
  chain_seeds = chain_seeds
 )
 nchains <- chain_args$nchains
 keep_chains <- chain_args$keep_chains
 chain_seeds <- chain_args$chain_seeds

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

 group_info <- .stblr_prepare_group_index(
  group = group,
  m = m,
  variable_names = variable_names,
  group_names = group_names
 )

 group_priors <- .stblr_make_group_priors(
  ngroup = group_info$ngroup,
  nt = nt,
  pi_init = arch$pi_init,
  group_pi_init = group_pi_init,
  group_vb_multiplier_init = group_vb_multiplier_init,
  pi_group_prior_mean = pi_group_prior_mean,
  pi_group_prior_strength = pi_group_prior_strength,
  pi_group_prior_a = pi_group_prior_a,
  pi_group_prior_b = pi_group_prior_b,
  group_names = group_info$group_names,
  trait_names = trait_names
 )

 raw_fit <- stblr_cpg_omp_csr_group_annot(
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
  group_index = as.integer(group_info$group_index0),
  ngroup = as.integer(group_info$ngroup),
  group_pi_init = group_priors$group_pi_init,
  pi_group_prior_a = group_priors$pi_group_prior_a,
  pi_group_prior_b = group_priors$pi_group_prior_b,
  group_vb_multiplier_init = group_priors$group_vb_multiplier_init,
  updateGroupVb = updateGroupVb,
  nub_group = nub_group,
  ssb_group_prior = ssb_group_prior,
  normalize_group_vb = normalize_group_vb,
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
  ncores = as.integer(ncores),
  seed = as.integer(seed),
  nchains = nchains,
  keep_chains = keep_chains,
  chain_seeds = chain_seeds,
  updateLDswap = updateLDswap,
  ld_swap_prob = ld_swap_prob,
  ld_swap_r2 = ld_swap_r2,
  ld_swap_max_friends = as.integer(ld_swap_max_friends),
  ld_swap_moves = as.integer(ld_swap_moves)
 )

 if (.is_stblr_raw_v1(raw_fit)) {
  raw_fit$group$group_names <- group_info$group_names
  fit <- .format_stblr_raw_v1(
   raw_fit,
   trait_names = trait_names,
   variable_names = variable_names
  )
 } else {
  fit <- .format_csr_group_annot_fit(
   fit = raw_fit,
   nt = nt,
   m = m,
   ngroup = group_info$ngroup,
   group_names = group_info$group_names,
   trait_names = trait_names,
   variable_names = variable_names
  )
 }

 fit$input <- c(
  list(
   model = "group",
   n = n,
   m = m,
   nt = nt,
   group = group_info$group,
   group_index0 = group_info$group_index0,
   ngroup = group_info$ngroup,
   group_names = group_info$group_names,
   group_size = group_info$group_size,
   group_pi_init = group_priors$group_pi_init_matrix,
   group_vb_multiplier_init = group_priors$group_vb_multiplier_init_matrix,
   pi_group_prior_mean = group_priors$pi_group_prior_mean,
   pi_group_prior_strength = group_priors$pi_group_prior_strength,
   pi_group_prior_a = group_priors$pi_group_prior_a,
   pi_group_prior_b = group_priors$pi_group_prior_b,
   updateGroupVb = updateGroupVb,
   nub_group = nub_group,
   ssb_group_prior = ssb_group_prior,
   normalize_group_vb = normalize_group_vb,
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
   nchains = nchains,
   keep_chains = keep_chains,
   chain_seeds = if (length(chain_seeds)) chain_seeds else NULL,
   chain_seed_rule = if (length(chain_seeds)) {
    "chain_seeds[chain] + 1000003 * (trait + 1)"
   } else if (nchains == 1L) {
    "seed + 1000003 * (trait + 1)"
   } else {
    "seed + 1000003 * (trait + 1) + 9176 * (chain + 1)"
   },
   updateLDswap = updateLDswap,
   ld_swap_prob = ld_swap_prob,
   ld_swap_r2 = ld_swap_r2,
   ld_swap_max_friends = as.integer(ld_swap_max_friends),
   ld_swap_moves = as.integer(ld_swap_moves),
   use_d_init = use_d_init,
   use_r_init = use_r_init,
   rebuild_r_before_updateE = rebuild_r_before_updateE,
   ld_prefix = ld_prefix
  ),
  arch
 )
 .standardize_stblr_annotation_fit(fit, "group")
}
