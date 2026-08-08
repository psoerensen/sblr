# Internal Phase-5C SBayesRC-S-EM helpers. The genomic E-step is inherited
# from SBayesRC-EM; the shared-delta outer layer targets a distinct MAP model.

.sbayesrc_s_em_log_sum_exp <- function(x) {
 x <- as.numeric(x)
 maximum <- max(x)
 maximum + log(sum(exp(x - maximum)))
}

.sbayesrc_s_em_model_key <- function(delta) {
 paste0(as.integer(delta), collapse = "")
}

.sbayesrc_s_em_validate_hyperparameters <- function(pi_a, tau2, stick_count) {
 if (!is.numeric(pi_a) || length(pi_a) != 1L || !is.finite(pi_a) ||
     pi_a <= 0 || pi_a >= 1) {
  stop("pi_a must be a finite scalar in (0, 1).")
 }
 tau2 <- as.numeric(tau2)
 if (length(tau2) != stick_count || any(!is.finite(tau2)) ||
     any(tau2 <= 0)) {
  stop("tau2 must contain one positive finite slab variance per stick.")
 }
 list(pi_a = as.numeric(pi_a), tau2 = tau2)
}

.sbayesrc_s_em_stick_laplace <- function(
  A, eligible, success, selected, tau2, intercept_mean,
  intercept_precision, start, probability_floor = 1e-12, maxit = 500L
) {
 columns <- c(1L, which(selected) + 1L)
 design <- A[, columns, drop = FALSE]
 prior_mean <- c(intercept_mean, rep(0, sum(selected)))
 prior_precision <- c(intercept_precision, rep(1 / tau2, sum(selected)))
 start <- as.numeric(start[columns])
 if (length(start) != ncol(design) || any(!is.finite(start))) {
  start <- prior_mean
 }
 log_joint <- function(parameter) {
  eta <- drop(design %*% parameter)
  probability <- .sbayesrc_mcem_clamp_probability(
   stats::pnorm(eta), floor = probability_floor
  )
  log_likelihood <- sum(
   success * log(probability) +
    (eligible - success) * log1p(-probability)
  )
  log_prior <- 0.5 * sum(log(prior_precision)) -
   0.5 * length(parameter) * log(2 * pi) -
   0.5 * sum(prior_precision * (parameter - prior_mean)^2)
  log_likelihood + log_prior
 }
 fit <- stats::optim(
  par = start, fn = function(parameter) -log_joint(parameter),
  method = "BFGS", hessian = TRUE,
  control = list(maxit = as.integer(maxit), reltol = 1e-12)
 )
 hessian <- (fit$hessian + t(fit$hessian)) / 2
 eigenvalue <- eigen(hessian, symmetric = TRUE, only.values = TRUE)$values
 if (any(!is.finite(eigenvalue)) || min(eigenvalue) <= 1e-10) {
  hessian <- hessian + diag(max(1e-8, 1e-8 - min(eigenvalue)), nrow(hessian))
 }
 chol_hessian <- chol(hessian)
 log_det_hessian <- 2 * sum(log(diag(chol_hessian)))
 covariance <- chol2inv(chol_hessian)
 log_marginal <- log_joint(fit$par) +
  0.5 * length(fit$par) * log(2 * pi) - 0.5 * log_det_hessian
 full_mode <- rep(0, ncol(A))
 full_mode[columns] <- fit$par
 full_covariance <- matrix(0, ncol(A), ncol(A))
 full_covariance[columns, columns] <- covariance
 list(
  log_marginal = log_marginal,
  mode = full_mode,
  covariance = full_covariance,
  convergence = fit$convergence,
  objective = fit$value
 )
}

.sbayesrc_s_em_model_laplace <- function(
  A, responsibility, delta, pi_a, tau2, intercept_prior_resolved,
  alpha_start, probability_floor = 1e-12, maxit = 500L
) {
 A <- as.matrix(A)
 responsibility <- as.matrix(responsibility)
 alpha_start <- as.matrix(alpha_start)
 delta <- as.integer(delta)
 stick_count <- ncol(responsibility) - 1L
 selectable_count <- ncol(A) - 1L
 if (length(delta) != selectable_count || any(!delta %in% 0:1) ||
     !identical(dim(alpha_start), c(ncol(A), stick_count))) {
  stop("delta or alpha_start has incompatible dimensions.")
 }
 hyper <- .sbayesrc_s_em_validate_hyperparameters(pi_a, tau2, stick_count)
 intercept <- .sbayesrc_mcem_intercept_prior(
  intercept_prior_resolved, stick_count
 )
 soft <- .sbayesrc_mcem_soft_stick_information(responsibility)
 sticks <- vector("list", stick_count)
 alpha_mode <- matrix(0, ncol(A), stick_count,
                      dimnames = dimnames(alpha_start))
 log_marginal <- 0
 for (stick in seq_len(stick_count)) {
  sticks[[stick]] <- .sbayesrc_s_em_stick_laplace(
   A, soft$eligible[, stick], soft$success[, stick], delta == 1L,
   hyper$tau2[stick], intercept$mean[stick], intercept$precision[stick],
   alpha_start[, stick], probability_floor, maxit
  )
  alpha_mode[, stick] <- sticks[[stick]]$mode
  log_marginal <- log_marginal + sticks[[stick]]$log_marginal
 }
 included <- sum(delta)
 log_model_prior <- included * log(hyper$pi_a) +
  (selectable_count - included) * log1p(-hyper$pi_a)
 list(
  delta = delta,
  log_weight = log_marginal + log_model_prior,
  log_marginal = log_marginal,
  log_model_prior = log_model_prior,
  alpha_mode = alpha_mode,
  sticks = sticks,
  convergence = max(vapply(sticks, `[[`, integer(1L), "convergence"))
 )
}

.sbayesrc_s_em_selection_update <- function(
  A, responsibility, delta_start, alpha_start, pi_a, tau2,
  intercept_prior_resolved, sweeps = 2000L, burn = 500L,
  seed = 1L, probability_floor = 1e-12, maxit = 500L
) {
 selectable_count <- ncol(A) - 1L
 delta <- as.integer(delta_start)
 if (length(delta) != selectable_count || any(!delta %in% 0:1)) {
  stop("delta_start must contain one binary state per nonintercept annotation.")
 }
 sweeps <- as.integer(sweeps)
 burn <- as.integer(burn)
 if (sweeps < 1L || burn < 0L || burn >= sweeps) {
  stop("selection sweeps must be positive and burn must be smaller than sweeps.")
 }
 cache <- new.env(parent = emptyenv())
 evaluate <- function(state) {
  key <- .sbayesrc_s_em_model_key(state)
  if (!exists(key, envir = cache, inherits = FALSE)) {
   assign(key, .sbayesrc_s_em_model_laplace(
    A, responsibility, state, pi_a, tau2, intercept_prior_resolved,
    alpha_start, probability_floor, maxit
   ), envir = cache)
  }
  get(key, envir = cache, inherits = FALSE)
 }
 set.seed(as.integer(seed))
 retained <- sweeps - burn
 delta_sum <- numeric(selectable_count)
 alpha_sum <- matrix(0, nrow(alpha_start), ncol(alpha_start))
 model_count <- integer(0)
 best <- evaluate(delta)
 transition <- integer(selectable_count)
 previous <- delta
 for (iteration in seq_len(sweeps)) {
  for (annotation in seq_len(selectable_count)) {
   excluded <- included <- delta
   excluded[annotation] <- 0L
   included[annotation] <- 1L
   model0 <- evaluate(excluded)
   model1 <- evaluate(included)
   log_norm <- .sbayesrc_s_em_log_sum_exp(c(
    model0$log_weight, model1$log_weight
   ))
   probability <- exp(model1$log_weight - log_norm)
   delta[annotation] <- stats::rbinom(1L, 1L, probability)
  }
  current <- evaluate(delta)
  if (current$log_weight > best$log_weight) best <- current
  if (iteration > burn) {
   delta_sum <- delta_sum + delta
   alpha_sum <- alpha_sum + current$alpha_mode
   key <- .sbayesrc_s_em_model_key(delta)
   if (!key %in% names(model_count)) model_count[key] <- 0L
   model_count[key] <- model_count[key] + 1L
   transition <- transition + as.integer(delta != previous)
  }
  previous <- delta
 }
 model_probability <- model_count / sum(model_count)
 list(
  delta_map = best$delta,
  alpha_map = best$alpha_mode,
  annotation_pip_eb = delta_sum / retained,
  alpha_model_average = alpha_sum / retained,
  model_probability = model_probability,
  model_cache_size = length(ls(cache, all.names = TRUE)),
  transition_count = transition,
  max_convergence = max(vapply(
   mget(ls(cache, all.names = TRUE), envir = cache),
   `[[`, integer(1L), "convergence"
  )),
  target = "responsibility_conditioned_laplace_model_distribution"
 )
}

.sbayesrc_s_em_engine <- function(
  state, B, E, A, gamma, alpha_init, delta_init, pi_a, tau2,
  intercept_prior_resolved, inner_function,
  inner_sweeps, inner_burn, final_sweeps, final_burn,
  selection_sweeps, selection_burn,
  damping, tol_alpha, tol_prior, min_outer, max_outer,
  pi_floor, seed, backend, updateB, updateE,
  return_responsibilities, verbose
) {
 alpha <- as.matrix(alpha_init)
 delta <- as.integer(delta_init)
 excluded_rows <- which(delta == 0L) + 1L
 if (length(excluded_rows)) alpha[excluded_rows, ] <- 0
 prior <- .sbayesrc_mcem_component_prior(A, alpha, pi_floor)
 history <- vector("list", max_outer)
 alpha_history <- delta_history <- pip_history <- vector("list", max_outer)
 converged <- FALSE
 last_responsibility <- last_selection <- NULL
 completed <- 0L
 for (outer in seq_len(max_outer)) {
  inner <- inner_function(state, B, E, alpha, inner_sweeps, inner_burn,
                          seed + outer - 1L, TRUE)
  chain <- inner$chains[[1L]][[1L]]
  responsibility <- chain$information_flow$rb_comp_prob
  if (is.null(responsibility)) {
   stop("The SBayesRC-S-EM inner kernel did not return RB responsibilities.")
  }
  selection <- .sbayesrc_s_em_selection_update(
   A, responsibility, delta, alpha, pi_a, tau2,
   intercept_prior_resolved, selection_sweeps, selection_burn,
   seed + 100000L + outer, pi_floor
  )
  delta_new <- selection$delta_map
  alpha_target <- selection$alpha_map
  alpha_new <- .sbayesrc_mcem_damped_update(alpha, alpha_target, damping)
  excluded_rows <- which(delta_new == 0L) + 1L
  if (length(excluded_rows)) alpha_new[excluded_rows, ] <- 0
  prior_new <- .sbayesrc_mcem_component_prior(A, alpha_new, pi_floor)
  max_delta_alpha <- max(abs(alpha_new - alpha))
  max_delta_prior <- max(abs(prior_new - prior))
  delta_changed <- sum(delta_new != delta)
  completed <- outer
  history[[outer]] <- data.frame(
   outer = outer,
   max_delta_alpha = max_delta_alpha,
   max_delta_prior = max_delta_prior,
   delta_changed = delta_changed,
   prior_expected_active = sum(1 - prior_new[, 1L]),
   rb_expected_active = sum(1 - responsibility[, 1L]),
   B = as.numeric(inner$variance$vb[1L, 1L]),
   E = as.numeric(inner$variance$ve[1L, 1L]),
   selection_models_visited = selection$model_cache_size,
   selection_transitions = sum(selection$transition_count),
   selection_max_convergence = selection$max_convergence,
   inner_sweeps = inner_sweeps,
   inner_burn = inner_burn,
   stringsAsFactors = FALSE
  )
  alpha_history[[outer]] <- alpha_new
  delta_history[[outer]] <- delta_new
  pip_history[[outer]] <- selection$annotation_pip_eb
  state <- .sbayesrc_mcem_extract_state(inner, state)
  if (isTRUE(updateB)) B[1L, 1L] <- history[[outer]]$B
  if (isTRUE(updateE)) E[1L, 1L] <- history[[outer]]$E
  alpha <- alpha_new
  delta <- delta_new
  prior <- prior_new
  last_responsibility <- responsibility
  last_selection <- selection
  if (isTRUE(verbose)) {
   message(sprintf(
    paste0("SBayesRC-S-EM backend=%s outer=%d max_dAlpha=%.6g ",
           "max_dPrior=%.6g delta_changes=%d"),
    backend, outer, max_delta_alpha, max_delta_prior, delta_changed
   ))
  }
  if (delta_changed == 0L && .sbayesrc_mcem_has_converged(
      outer, min_outer, max_delta_alpha, max_delta_prior,
      tol_alpha, tol_prior)) {
   converged <- TRUE
   break
  }
 }
 final <- inner_function(
  state, B, E, alpha, final_sweeps, final_burn,
  seed + max_outer + 20000L, isTRUE(return_responsibilities)
 )
 final_state <- .sbayesrc_mcem_extract_state(final, state)
 if (is.matrix(final$marker$b)) {
  final$marker$b[, 1L] <- final_state$b[[1L]]
 } else {
  final$marker$b[[1L]] <- final_state$b[[1L]]
 }
 result <- list(
  method = "SBayesRC-S-EM",
  algorithm = "MCEM-Laplace",
  inference = "mcem",
  target = "observed_data_delta_alpha_MAP_with_conditional_Laplace_EB_PIP",
  backend = backend,
  converged = converged,
  n_outer = completed,
  delta_map = delta,
  alpha_map = alpha,
  alpha_model_average = last_selection$alpha_model_average,
  annotation_pip_eb = last_selection$annotation_pip_eb,
  component_prior = prior,
  history = list(
   summary = do.call(rbind, history[seq_len(completed)]),
   alpha = alpha_history[seq_len(completed)],
   delta = delta_history[seq_len(completed)],
   annotation_pip_eb = pip_history[seq_len(completed)]
  ),
  damping = damping,
  tol_alpha = tol_alpha,
  tol_prior = tol_prior,
  min_outer = min_outer,
  max_outer = max_outer,
  inner_sweeps = inner_sweeps,
  inner_burn = inner_burn,
  selection_sweeps = selection_sweeps,
  selection_burn = selection_burn,
  final_sweeps = final_sweeps,
  final_burn = final_burn,
  pi_A_mode = "fixed",
  pi_A_fixed = pi_a,
  tau2_mode = "fixed",
  tau2_fixed = tau2,
  intercept_prior_mode = "proper_normal",
  mixture_prior_mode = "annotation_stick_intercepts_no_global_Pi_update",
  genomic_hyperparameters = list(
   updateB = isTRUE(updateB), updateE = isTRUE(updateE),
   B_final = final$variance$vb, E_final = final$variance$ve
  ),
  last_estep_responsibilities = if (isTRUE(return_responsibilities)) {
   last_responsibility
  } else NULL,
  final_genomic_responsibilities = if (isTRUE(return_responsibilities)) {
   final$chains[[1L]][[1L]]$information_flow$rb_comp_prob
  } else NULL
 )
 structure(
  list(genomic = final, mcem = result),
  class = c("sblr_sbayesrc_s_em_phase5c", "list")
 )
}

.stblr_mcem_sbayesrc_s_csr <- function(
  wy, ww, yy, b_init, comp_init, r_init, ld_prefix,
  B, E, ssb_prior, sse_prior, A, gamma, alpha_init, delta_init,
  pi_a, tau2, intercept_prior_resolved, n,
  inner_sweeps = 1000L, inner_burn = 300L,
  selection_sweeps = 2000L, selection_burn = 500L,
  final_sweeps = inner_sweeps, final_burn = inner_burn,
  damping = 0.5, tol_alpha = 1e-3, tol_prior = 1e-3,
  min_outer = 3L, max_outer = 50L,
  pi_floor = 1e-12, nub = 4, nue = 4,
  ncores = 1L, seed = 10L, return_responsibilities = TRUE,
  verbose = FALSE
) {
 A <- as.matrix(A)
 alpha <- .sbayesrc_validate_alpha(as.matrix(alpha_init), gamma)
 gamma <- .sbayesrc_validate_gamma(gamma)
 if (length(wy) != 1L || length(ww) != 1L || length(yy) != 1L ||
     length(b_init) != 1L || length(comp_init) != 1L ||
     length(r_init) != 1L || length(n) != 1L ||
     nrow(A) != length(wy[[1L]]) || nrow(alpha) != ncol(A)) {
  stop("CSR SBayesRC-S-EM requires conformable single-trait inputs.")
 }
 hyper <- .sbayesrc_s_em_validate_hyperparameters(pi_a, tau2, ncol(alpha))
 .sbayesrc_mcem_intercept_prior(intercept_prior_resolved, ncol(alpha))
 control <- .sbayesrc_mcem_validate_controls(
  inner_sweeps, inner_burn, final_sweeps, final_burn,
  min_outer, max_outer, ncores, seed, damping, tol_alpha, tol_prior
 )
 state <- list(
  b = lapply(b_init, as.numeric),
  component = lapply(comp_init, as.numeric),
  r = lapply(r_init, as.numeric)
 )
 inner_function <- function(state, B, E, alpha, sweeps, burn, seed, capture) {
  .sbayesrc_mcem_inner_csr(
   wy, ww, yy, state, ld_prefix, B, E, ssb_prior, sse_prior,
   A, gamma, alpha, hyper$tau2, intercept_prior_resolved, n,
   sweeps, burn, seed, control$ncores, pi_floor, nub, nue, 0, capture
  )
 }
 .sbayesrc_s_em_engine(
  state, B, E, A, gamma, alpha, delta_init, hyper$pi_a, hyper$tau2,
  intercept_prior_resolved, inner_function,
  control$inner_sweeps, control$inner_burn,
  control$final_sweeps, control$final_burn,
  selection_sweeps, selection_burn,
  control$damping, control$tol_alpha, control$tol_prior,
  control$min_outer, control$max_outer, pi_floor, control$seed,
  "csr_reference", FALSE, FALSE, return_responsibilities, verbose
 )
}

.stblr_mcem_sbayesrc_s_block_eigen <- function(
  stats, Glist, annotation, block_start,
  B, E, ssb_prior, sse_prior, gamma, alpha_init, delta_init,
  pi_a, tau2, intercept_prior_resolved,
  b_init = NULL, comp_init = NULL,
  representation = "low_rank", eigen_prop = 0.995,
  eigen_filter = "hard_truncate", eigen_tau = 0.01, eigen_eta = 0,
  residual_policy = "gctb_block", block_ve_mode = "allMixVe",
  resam_thresh = 1.1, minimum_ve_ratio = 0.7,
  low_rank_residual_rebuild_every = 100L,
  updateB = FALSE, updateE = FALSE,
  inner_sweeps = 1000L, inner_burn = 300L,
  selection_sweeps = 2000L, selection_burn = 500L,
  final_sweeps = inner_sweeps, final_burn = inner_burn,
  damping = 0.5, tol_alpha = 1e-3, tol_prior = 1e-3,
  min_outer = 3L, max_outer = 50L,
  pi_floor = 1e-12, nub = 4, nue = 4,
  ncores = 1L, seed = 10L, return_responsibilities = TRUE,
  verbose = FALSE
) {
 A <- as.matrix(annotation)
 alpha <- .sbayesrc_validate_alpha(as.matrix(alpha_init), gamma)
 gamma <- .sbayesrc_validate_gamma(gamma)
 marker_count <- nrow(A)
 if (marker_count < 1L || ncol(A) < 2L || nrow(alpha) != ncol(A)) {
  stop("Block SBayesRC-S-EM annotation and alpha dimensions are inconsistent.")
 }
 hyper <- .sbayesrc_s_em_validate_hyperparameters(pi_a, tau2, ncol(alpha))
 intercept <- .sbayesrc_mcem_intercept_prior(
  intercept_prior_resolved, ncol(alpha)
 )
 control <- .sbayesrc_mcem_validate_controls(
  inner_sweeps, inner_burn, final_sweeps, final_burn,
  min_outer, max_outer, ncores, seed, damping, tol_alpha, tol_prior
 )
 if (!identical(residual_policy, "gctb_block")) {
  stop("Phase-5C block MCEM requires residual_policy = 'gctb_block'.")
 }
 resolved_mode <- if (isTRUE(updateE)) "allMixVe" else "fixVe"
 if (!identical(block_ve_mode, resolved_mode)) {
  stop("block_ve_mode must be '", resolved_mode,
       "' for the requested Phase-5C E-update mode.")
 }
 if (is.null(b_init)) b_init <- list(rep(0, marker_count))
 if (is.null(comp_init)) comp_init <- list(rep(0L, marker_count))
 state <- list(
  b = lapply(b_init, as.numeric),
  component = lapply(comp_init, as.numeric),
  r = list(rep(0, marker_count))
 )
 intercept_spec <- list(
  distribution = "normal", mean = intercept$mean,
  sd = 1 / sqrt(intercept$precision)
 )
 ssb_input <- if (is.list(ssb_prior) && length(ssb_prior) == 1L) {
  as.numeric(ssb_prior[[1L]])
 } else ssb_prior
 sse_input <- if (is.list(sse_prior) && length(sse_prior) == 1L) {
  as.numeric(sse_prior[[1L]])
 } else sse_prior
 inner_function <- function(state, B, E, alpha, sweeps, burn, seed, capture) {
  .stblr_csr_sbayesrc_block_eigen(
   stats = stats, Glist = Glist, annotation = A,
   block_start = block_start, representation = representation,
   eigen_prop = eigen_prop, eigen_filter = eigen_filter,
   eigen_tau = eigen_tau, eigen_eta = eigen_eta,
   low_rank_residual_rebuild_every = low_rank_residual_rebuild_every,
   residual_policy = residual_policy, block_ve_mode = resolved_mode,
   resam_thresh = resam_thresh, minimum_ve_ratio = minimum_ve_ratio,
   block_ve_keep_history = FALSE, gamma = gamma, B = B, E = E,
   ssb_prior = ssb_input, sse_prior = sse_input,
   updateAlpha = FALSE, updateB = updateB, updateE = updateE, adjE = 0,
   nit = sweeps, nburn = burn, nthin = 1L, ncores = control$ncores,
   seed = seed, nchains = 1L, keep_chains = TRUE,
   b_init = state$b, comp_init = state$component, use_comp_init = TRUE,
   use_r_init = FALSE, add_intercept = FALSE,
   standardize_annotations = FALSE, center_binary_annotations = FALSE,
   alpha_init = alpha, sigmaSqAlpha_init = hyper$tau2,
   annotation_intercept_prior = intercept_spec,
   pi_floor = pi_floor, nub = nub, nue = nue,
   .diagnostic_updateSigmaSqAlpha = FALSE,
   .information_diagnostics = capture, .return_raw = TRUE
  )
 }
 .sbayesrc_s_em_engine(
  state, as.matrix(B), as.matrix(E), A, gamma, alpha, delta_init,
  hyper$pi_a, hyper$tau2, intercept_prior_resolved, inner_function,
  control$inner_sweeps, control$inner_burn,
  control$final_sweeps, control$final_burn,
  selection_sweeps, selection_burn,
  control$damping, control$tol_alpha, control$tol_prior,
  control$min_outer, control$max_outer, pi_floor, control$seed,
  "block_eigen", updateB, updateE, return_responsibilities, verbose
 )
}
