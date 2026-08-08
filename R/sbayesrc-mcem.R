# Internal Phase-5A MCEM-SBayesRC helpers.  This is a distinct observed-data
# MAP / empirical-Bayes procedure; it is not the joint SBayesRC sampler.

.sbayesrc_mcem_clamp_probability <- function(x, floor = 1e-12) {
 if (!is.numeric(x) || any(!is.finite(x)) || !is.numeric(floor) ||
     length(floor) != 1L || !is.finite(floor) || floor <= 0 || floor >= 0.5) {
  stop("probabilities and floor must be finite, with floor in (0, 0.5).")
 }
 pmin(pmax(x, floor), 1 - floor)
}

.sbayesrc_mcem_component_prior <- function(A, alpha, floor = 1e-12) {
 A <- as.matrix(A)
 alpha <- as.matrix(alpha)
 if (!is.numeric(A) || !is.numeric(alpha) || any(!is.finite(A)) ||
     any(!is.finite(alpha)) || ncol(A) != nrow(alpha) ||
     ncol(alpha) < 1L) {
  stop("A and alpha must be finite, conformable numeric matrices.")
 }
 q <- matrix(
  .sbayesrc_mcem_clamp_probability(
   stats::pnorm(drop(A %*% alpha)), floor = floor
  ),
  nrow = nrow(A), ncol = ncol(alpha)
 )
 probability <- matrix(0, nrow(A), ncol(alpha) + 1L)
 remaining <- rep(1, nrow(A))
 for (stick in seq_len(ncol(alpha))) {
  probability[, stick] <- remaining * (1 - q[, stick])
  remaining <- remaining * q[, stick]
 }
 probability[, ncol(probability)] <- remaining
 probability <- pmax(probability, floor)
 probability <- probability / rowSums(probability)
 colnames(probability) <- paste0("component_", seq_len(ncol(probability)) - 1L)
 probability
}

.sbayesrc_mcem_soft_stick_information <- function(responsibility) {
 responsibility <- as.matrix(responsibility)
 if (!is.numeric(responsibility) || ncol(responsibility) < 2L ||
     any(!is.finite(responsibility)) || any(responsibility < 0) ||
     any(abs(rowSums(responsibility) - 1) > 1e-8)) {
  stop("responsibility must be a finite row-normalized probability matrix.")
 }
 component_count <- ncol(responsibility)
 stick_count <- component_count - 1L
 eligible <- success <- matrix(0, nrow(responsibility), stick_count)
 for (stick in seq_len(stick_count)) {
  eligible[, stick] <- rowSums(
   responsibility[, stick:component_count, drop = FALSE]
  )
  success[, stick] <- rowSums(
   responsibility[, (stick + 1L):component_count, drop = FALSE]
  )
 }
 colnames(eligible) <- colnames(success) <- paste0("stick_", seq_len(stick_count))
 list(
  eligible = eligible,
  success = success,
  rate = ifelse(eligible > 1e-14, success / eligible, 0)
 )
}

.sbayesrc_mcem_damped_update <- function(alpha, target, damping) {
 alpha <- as.matrix(alpha)
 target <- as.matrix(target)
 if (!identical(dim(alpha), dim(target)) || any(!is.finite(alpha)) ||
     any(!is.finite(target)) || !is.numeric(damping) || length(damping) != 1L ||
     !is.finite(damping) || damping <= 0 || damping > 1) {
  stop("alpha, target, and damping do not define a valid damped update.")
 }
 (1 - damping) * alpha + damping * target
}

.sbayesrc_mcem_has_converged <- function(outer, min_outer, delta_alpha,
                                         delta_prior, tol_alpha, tol_prior) {
 outer >= min_outer && delta_alpha < tol_alpha && delta_prior < tol_prior
}

.sbayesrc_mcem_intercept_prior <- function(intercept_prior_resolved,
                                           stick_count) {
 prior <- as.matrix(intercept_prior_resolved)
 if (!is.numeric(prior) || ncol(prior) != stick_count || nrow(prior) < 3L ||
     any(!is.finite(prior[seq_len(3L), , drop = FALSE]))) {
  stop("intercept_prior_resolved has invalid dimensions or values.")
 }
 row_index <- function(name, fallback) {
  if (!is.null(rownames(prior)) && name %in% rownames(prior)) {
   match(name, rownames(prior))
  } else {
   fallback
  }
 }
 type <- prior[row_index("type", 1L), ]
 mean <- prior[row_index("mean", 2L), ]
 precision <- prior[row_index("precision", 3L), ]
 if (any(type != 0) || any(!is.finite(mean)) ||
     any(!is.finite(precision)) || any(precision <= 0)) {
  stop("Phase-5A MCEM requires a proper normal intercept prior for every stick.")
 }
 list(mean = as.numeric(mean), precision = as.numeric(precision))
}

.sbayesrc_mcem_m_step <- function(A, responsibility, alpha_start,
                                  intercept_prior_resolved,
                                  sigmaSqAlpha, probability_floor = 1e-12,
                                  maxit = 500L) {
 A <- as.matrix(A)
 responsibility <- as.matrix(responsibility)
 alpha_start <- as.matrix(alpha_start)
 stick <- .sbayesrc_mcem_soft_stick_information(responsibility)
 stick_count <- ncol(responsibility) - 1L
 annotation_count <- ncol(A)
 if (!is.numeric(A) || any(!is.finite(A)) ||
     nrow(A) != nrow(responsibility) ||
     !identical(dim(alpha_start), c(annotation_count, stick_count)) ||
     any(!is.finite(alpha_start))) {
  stop("A, responsibility, and alpha_start have incompatible dimensions.")
 }
 if (!is.numeric(sigmaSqAlpha) || length(sigmaSqAlpha) != stick_count ||
     any(!is.finite(sigmaSqAlpha)) || any(sigmaSqAlpha <= 0)) {
  stop("sigmaSqAlpha must contain one positive finite variance per stick.")
 }
 if (!is.numeric(maxit) || length(maxit) != 1L || !is.finite(maxit) ||
     maxit < 1 || maxit != as.integer(maxit)) {
  stop("maxit must be a positive integer.")
 }
 intercept <- .sbayesrc_mcem_intercept_prior(
  intercept_prior_resolved, stick_count
 )
 alpha <- matrix(0, annotation_count, stick_count, dimnames = dimnames(alpha_start))
 convergence <- integer(stick_count)
 objective <- numeric(stick_count)
 for (stick_index in seq_len(stick_count)) {
  eligible <- stick$eligible[, stick_index]
  success <- stick$success[, stick_index]
  prior_mean <- c(intercept$mean[stick_index], rep(0, annotation_count - 1L))
  prior_precision <- c(
   intercept$precision[stick_index],
   rep(1 / sigmaSqAlpha[stick_index], annotation_count - 1L)
  )
  fn <- function(parameter) {
   eta <- drop(A %*% parameter)
   probability <- .sbayesrc_mcem_clamp_probability(
    stats::pnorm(eta), floor = probability_floor
   )
   log_likelihood <- sum(
    success * log(probability) +
     (eligible - success) * log1p(-probability)
   )
   log_prior <- -0.5 * sum(
    prior_precision * (parameter - prior_mean)^2
   )
   -(log_likelihood + log_prior)
  }
  gr <- function(parameter) {
   eta <- drop(A %*% parameter)
   probability <- .sbayesrc_mcem_clamp_probability(
    stats::pnorm(eta), floor = probability_floor
   )
   derivative <- stats::dnorm(eta) * (
    success / probability -
     (eligible - success) / (1 - probability)
   )
   likelihood_gradient <- drop(crossprod(A, derivative))
   prior_gradient <- -prior_precision * (parameter - prior_mean)
   -(likelihood_gradient + prior_gradient)
  }
  fit <- stats::optim(
   par = alpha_start[, stick_index], fn = fn, gr = gr, method = "BFGS",
   control = list(maxit = as.integer(maxit), reltol = 1e-12)
  )
  alpha[, stick_index] <- fit$par
  convergence[stick_index] <- fit$convergence
  objective[stick_index] <- fit$value
 }
 list(
  alpha = alpha,
  convergence = convergence,
  objective = objective,
  soft_stick = stick
 )
}

.sbayesrc_mcem_objective <- function(
  A, responsibility, alpha, intercept_prior_resolved,
  sigmaSqAlpha, probability_floor = 1e-12
) {
 A <- as.matrix(A)
 responsibility <- as.matrix(responsibility)
 alpha <- as.matrix(alpha)
 stick <- .sbayesrc_mcem_soft_stick_information(responsibility)
 stick_count <- ncol(responsibility) - 1L
 if (!identical(dim(alpha), c(ncol(A), stick_count))) {
  stop("alpha has incompatible dimensions for the MCEM objective.")
 }
 intercept <- .sbayesrc_mcem_intercept_prior(
  intercept_prior_resolved, stick_count
 )
 q_annotation <- log_prior_alpha <- numeric(stick_count)
 for (index in seq_len(stick_count)) {
  probability <- .sbayesrc_mcem_clamp_probability(
   stats::pnorm(drop(A %*% alpha[, index])), floor = probability_floor
  )
  q_annotation[index] <- sum(
   stick$success[, index] * log(probability) +
    (stick$eligible[, index] - stick$success[, index]) *
     log1p(-probability)
  )
  prior_mean <- c(intercept$mean[index], rep(0, ncol(A) - 1L))
  prior_precision <- c(
   intercept$precision[index], rep(1 / sigmaSqAlpha[index], ncol(A) - 1L)
  )
  log_prior_alpha[index] <- -0.5 * sum(
   prior_precision * (alpha[, index] - prior_mean)^2
  )
 }
 list(
  Q_annotation_by_stick = q_annotation,
  log_prior_alpha_by_stick = log_prior_alpha,
  Q_annotation = sum(q_annotation),
  log_prior_alpha = sum(log_prior_alpha),
  Q_total = sum(q_annotation) + sum(log_prior_alpha)
 )
}

.sbayesrc_mcem_extract_state <- function(raw, initial_state) {
 if (!inherits(raw, "stblr_raw") || !is.list(raw$marker) ||
     !all(c("b", "r", "state") %in% names(raw$marker))) {
  stop("The MCEM inner kernel returned an invalid stblr_raw state.")
 }
 extract_one <- function(x, name) {
  if (is.matrix(x) && ncol(x) == 1L) return(as.numeric(x[, 1L]))
  if (!is.list(x) || length(x) != 1L) {
   stop("Phase-5A MCEM currently supports exactly one trait (", name, ").")
  }
  as.numeric(x[[1L]])
 }
 # The shared CSR aggregate starts its top-level b matrix at b_init before
 # adding the final per-chain state. With one chain, subtracting the known
 # input recovers the final state without touching the native transition.
 final_b <- extract_one(raw$marker$b, "b") - initial_state$b[[1L]]
 list(
  b = list(final_b),
  r = list(extract_one(raw$marker$r, "r")),
  component = list(extract_one(raw$marker$state, "state"))
 )
}

.sbayesrc_mcem_inner_csr <- function(
  wy, ww, yy, state, ld_prefix, B, E, ssb_prior, sse_prior,
  A, gamma, alpha, sigmaSqAlpha, intercept_prior_resolved, n,
  retained_sweeps, burn, seed, ncores, pi_floor, nub, nue, adjE,
  capture_responsibility = TRUE
) {
 dummy_delta <- rep.int(1L, ncol(A) - 1L)
 .st_sbayesrc_selection_csr(
  wy, ww, yy, state$b, state$component, TRUE, state$r, TRUE,
  ld_prefix, B, E, ssb_prior, sse_prior, A, gamma, alpha,
  dummy_delta, 0.5, sigmaSqAlpha, 1, 1, 2, 2, integer(),
  FALSE, FALSE, FALSE, intercept_prior_resolved, pi_floor, nub, nue,
  FALSE, FALSE, adjE, as.integer(n), as.integer(retained_sweeps),
  as.integer(burn), 1L, as.integer(ncores), as.integer(seed), 1L,
  as.integer(seed), FALSE, isTRUE(capture_responsibility)
 )
}

.sbayesrc_mcem_engine <- function(
  state, B, E, A, gamma, alpha_init, sigmaSqAlpha_init,
  intercept_prior_resolved, inner_function,
  inner_sweeps, inner_burn, final_sweeps, final_burn,
  damping, tol_alpha, tol_prior, min_outer, max_outer,
  pi_floor, seed, backend, updateB, updateE,
  return_responsibilities, verbose, objective_diagnostics = FALSE,
  diagnostic_responsibility_iterations = integer()
) {
 alpha <- as.matrix(alpha_init)
 prior <- .sbayesrc_mcem_component_prior(A, alpha, pi_floor)
 alpha_history <- array(
  NA_real_, c(nrow(alpha), ncol(alpha), max_outer + 1L),
  dimnames = list(
   annotation = rownames(alpha), stick = colnames(alpha),
   outer = 0:max_outer
  )
 )
 alpha_history[, , 1L] <- alpha
 history <- vector("list", max_outer)
 converged <- FALSE
 last_responsibility <- NULL
 completed <- 0L
 responsibility_checkpoint <- list()

 for (outer in seq_len(max_outer)) {
  inner <- inner_function(
   state, B, E, alpha, inner_sweeps, inner_burn,
   seed + outer - 1L, TRUE
  )
  chain <- inner$chains[[1L]][[1L]]
  responsibility <- chain$information_flow$rb_comp_prob
  if (is.null(responsibility)) {
   stop("The MCEM inner kernel did not return RB responsibilities.")
  }
  if (outer %in% diagnostic_responsibility_iterations) {
   responsibility_checkpoint[[as.character(outer)]] <- responsibility
  }
  m_step <- .sbayesrc_mcem_m_step(
   A, responsibility, alpha, intercept_prior_resolved,
   sigmaSqAlpha_init, probability_floor = pi_floor
  )
  alpha_target <- m_step$alpha
  alpha_new <- .sbayesrc_mcem_damped_update(alpha, alpha_target, damping)
  prior_new <- .sbayesrc_mcem_component_prior(A, alpha_new, pi_floor)
  max_delta_alpha <- max(abs(alpha_new - alpha))
  max_delta_prior <- max(abs(prior_new - prior))
  hard_probability <- chain$component$prob
  rb_hard <- if (is.matrix(hard_probability) &&
                 identical(dim(hard_probability), dim(responsibility))) {
   abs(responsibility - hard_probability)
  } else {
   NA_real_
  }
  rb_hard_mean <- if (all(is.na(rb_hard))) NA_real_ else
   mean(rb_hard, na.rm = TRUE)
  rb_hard_max <- if (all(is.na(rb_hard))) NA_real_ else
   max(rb_hard, na.rm = TRUE)
  objective <- if (isTRUE(objective_diagnostics)) {
   list(
    current = .sbayesrc_mcem_objective(
     A, responsibility, alpha, intercept_prior_resolved,
     sigmaSqAlpha_init, pi_floor
    ),
    target = .sbayesrc_mcem_objective(
     A, responsibility, alpha_target, intercept_prior_resolved,
     sigmaSqAlpha_init, pi_floor
    ),
    updated = .sbayesrc_mcem_objective(
     A, responsibility, alpha_new, intercept_prior_resolved,
     sigmaSqAlpha_init, pi_floor
    )
   )
  } else NULL
  completed <- outer
  history[[outer]] <- data.frame(
   outer = outer,
   max_delta_alpha = max_delta_alpha,
   max_delta_prior = max_delta_prior,
   prior_expected_active = sum(1 - prior_new[, 1L]),
   rb_expected_active = sum(1 - responsibility[, 1L]),
   B = as.numeric(inner$variance$vb[1L, 1L]),
   E = as.numeric(inner$variance$ve[1L, 1L]),
   rb_hard_mean_abs = rb_hard_mean,
   rb_hard_max_abs = rb_hard_max,
   m_step_max_convergence = max(m_step$convergence),
   m_step_objective = sum(m_step$objective),
   Q_annotation_current = if (is.null(objective)) NA_real_ else
    objective$current$Q_annotation,
   log_prior_alpha_current = if (is.null(objective)) NA_real_ else
    objective$current$log_prior_alpha,
   Q_total_current = if (is.null(objective)) NA_real_ else
    objective$current$Q_total,
   Q_annotation_target = if (is.null(objective)) NA_real_ else
    objective$target$Q_annotation,
   log_prior_alpha_target = if (is.null(objective)) NA_real_ else
    objective$target$log_prior_alpha,
   Q_total_target = if (is.null(objective)) NA_real_ else
    objective$target$Q_total,
   Q_total_updated = if (is.null(objective)) NA_real_ else
    objective$updated$Q_total,
   inner_sweeps = inner_sweeps,
   inner_burn = inner_burn,
   stringsAsFactors = FALSE
  )
  alpha_history[, , outer + 1L] <- alpha_new
  state <- .sbayesrc_mcem_extract_state(inner, state)
  if (isTRUE(updateB)) B[1L, 1L] <- history[[outer]]$B
  if (isTRUE(updateE)) E[1L, 1L] <- history[[outer]]$E
  alpha <- alpha_new
  prior <- prior_new
  last_responsibility <- responsibility
  if (isTRUE(verbose)) {
   message(sprintf(
    "SBayesRC-EM backend=%s outer=%d max_dAlpha=%.6g max_dPrior=%.6g RB_active=%.3f",
    backend, outer, max_delta_alpha, max_delta_prior,
    history[[outer]]$rb_expected_active
   ))
  }
  if (.sbayesrc_mcem_has_converged(
      outer, min_outer, max_delta_alpha, max_delta_prior,
      tol_alpha, tol_prior)) {
   converged <- TRUE
   break
  }
 }

 final <- inner_function(
  state, B, E, alpha, final_sweeps, final_burn,
  seed + max_outer + 10000L, isTRUE(return_responsibilities)
 )
 final_state <- .sbayesrc_mcem_extract_state(final, state)
 if (is.matrix(final$marker$b)) {
  final$marker$b[, 1L] <- final_state$b[[1L]]
 } else {
  final$marker$b[[1L]] <- final_state$b[[1L]]
 }
 history_summary <- do.call(rbind, history[seq_len(completed)])
 history_alpha <- alpha_history[, , seq_len(completed + 1L), drop = FALSE]
 mcem <- list(
  method = "SBayesRC-EM",
  algorithm = "MCEM",
  inference = "mcem",
  target = "observed_data_alpha_MAP_empirical_Bayes",
  backend = backend,
  converged = converged,
  n_outer = completed,
  alpha_map = alpha,
  component_prior = prior,
  history = list(
   summary = history_summary, alpha = history_alpha,
   responsibility_checkpoint = responsibility_checkpoint
  ),
  damping = damping,
  tol_alpha = tol_alpha,
  tol_prior = tol_prior,
  min_outer = min_outer,
  max_outer = max_outer,
  inner_sweeps = inner_sweeps,
  inner_burn = inner_burn,
  final_sweeps = final_sweeps,
  final_burn = final_burn,
  sigmaSqAlpha_mode = "fixed_prior_variance",
  sigmaSqAlpha_fixed = sigmaSqAlpha_init,
  mixture_prior_mode = "annotation_stick_intercepts_no_global_Pi_update",
  objective_diagnostics = isTRUE(objective_diagnostics),
  genomic_hyperparameters = list(
   updateB = isTRUE(updateB), updateE = isTRUE(updateE),
   B_final = final$variance$vb, E_final = final$variance$ve
  ),
  genomic_hyperparameters_fixed = !isTRUE(updateB) && !isTRUE(updateE),
  last_estep_responsibilities = if (isTRUE(return_responsibilities)) {
   last_responsibility
  } else NULL,
  final_genomic_responsibilities = if (isTRUE(return_responsibilities)) {
   final$chains[[1L]][[1L]]$information_flow$rb_comp_prob
  } else NULL,
  # Phase-5A compatibility aliases; the explicit names above are preferred.
  e_step_responsibilities = if (isTRUE(return_responsibilities)) {
   last_responsibility
  } else NULL,
  final_responsibilities = if (isTRUE(return_responsibilities)) {
   final$chains[[1L]][[1L]]$information_flow$rb_comp_prob
  } else NULL
 )
 structure(
  list(genomic = final, mcem = mcem),
  class = c("sblr_sbayesrc_em_phase5b", "list")
 )
}

.sbayesrc_mcem_validate_controls <- function(
  inner_sweeps, inner_burn, final_sweeps, final_burn,
  min_outer, max_outer, ncores, seed, damping, tol_alpha, tol_prior
) {
 integer_scalar <- function(x, name, minimum = 1L) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
      x != as.integer(x) || x < minimum) {
   stop(name, " must be an integer >= ", minimum, ".")
  }
  as.integer(x)
 }
 controls <- list(
  inner_sweeps = integer_scalar(inner_sweeps, "inner_sweeps"),
  inner_burn = integer_scalar(inner_burn, "inner_burn", 0L),
  final_sweeps = integer_scalar(final_sweeps, "final_sweeps"),
  final_burn = integer_scalar(final_burn, "final_burn", 0L),
  min_outer = integer_scalar(min_outer, "min_outer"),
  max_outer = integer_scalar(max_outer, "max_outer"),
  ncores = integer_scalar(ncores, "ncores"),
  seed = integer_scalar(seed, "seed", 0L)
 )
 if (controls$min_outer > controls$max_outer) {
  stop("min_outer cannot exceed max_outer.")
 }
 for (item in list(damping = damping, tol_alpha = tol_alpha,
                   tol_prior = tol_prior)) {
  if (!is.numeric(item) || length(item) != 1L || !is.finite(item) || item <= 0) {
   stop("damping, tol_alpha, and tol_prior must be positive finite scalars.")
  }
 }
 if (damping > 1) stop("damping must be in (0, 1].")
 c(controls, list(
  damping = as.numeric(damping), tol_alpha = as.numeric(tol_alpha),
  tol_prior = as.numeric(tol_prior)
 ))
}

.stblr_mcem_sbayesrc_csr <- function(
  wy, ww, yy, b_init, comp_init, r_init, ld_prefix,
  B, E, ssb_prior, sse_prior, A, gamma, alpha_init,
  sigmaSqAlpha_init, intercept_prior_resolved, n,
  inner_sweeps = 1000L, inner_burn = 300L,
  final_sweeps = inner_sweeps, final_burn = inner_burn,
  damping = 0.5, tol_alpha = 1e-3, tol_prior = 1e-3,
  min_outer = 3L, max_outer = 50L,
  pi_floor = 1e-12, nub = 4, nue = 4, adjE = 0,
  ncores = 1L, seed = 10L, return_responsibilities = TRUE,
  verbose = FALSE, .objective_diagnostics = FALSE,
  .diagnostic_responsibility_iterations = integer()
) {
 if (length(wy) != 1L || length(ww) != 1L || length(b_init) != 1L ||
     length(comp_init) != 1L || length(r_init) != 1L || length(yy) != 1L ||
     length(n) != 1L) {
  stop("Phase-5A MCEM currently supports exactly one trait.")
 }
 A <- as.matrix(A)
 alpha <- as.matrix(alpha_init)
 gamma <- .sbayesrc_validate_gamma(gamma)
 alpha <- .sbayesrc_validate_alpha(alpha, gamma)
 marker_count <- length(wy[[1L]])
 if (!is.numeric(A) || any(!is.finite(A)) || nrow(A) != marker_count ||
     ncol(A) < 2L || nrow(alpha) != ncol(A) ||
     length(ww[[1L]]) != marker_count || length(b_init[[1L]]) != marker_count ||
     length(comp_init[[1L]]) != marker_count || length(r_init[[1L]]) != marker_count) {
  stop("MCEM genomic inputs and annotation dimensions are inconsistent.")
 }
 sigmaSqAlpha_init <- as.numeric(sigmaSqAlpha_init)
 .sbayesrc_mcem_intercept_prior(intercept_prior_resolved, ncol(alpha))
 if (length(sigmaSqAlpha_init) != ncol(alpha) ||
     any(!is.finite(sigmaSqAlpha_init)) || any(sigmaSqAlpha_init <= 0)) {
  stop("sigmaSqAlpha_init must provide a fixed positive variance per stick.")
 }
 control <- .sbayesrc_mcem_validate_controls(
  inner_sweeps, inner_burn, final_sweeps, final_burn,
  min_outer, max_outer, ncores, seed, damping, tol_alpha, tol_prior
 )
 if (!is.numeric(adjE) || length(adjE) != 1L || !is.finite(adjE) ||
     adjE != 0) {
  stop("Phase-5A qualification requires adjE = 0 and fixed genomic variances.")
 }

 state <- list(
  b = lapply(b_init, as.numeric),
  component = lapply(comp_init, as.numeric),
  r = lapply(r_init, as.numeric)
 )
 inner_function <- function(state, B, E, alpha, sweeps, burn, seed, capture) {
  .sbayesrc_mcem_inner_csr(
   wy, ww, yy, state, ld_prefix, B, E, ssb_prior, sse_prior,
   A, gamma, alpha, sigmaSqAlpha_init, intercept_prior_resolved, n,
   sweeps, burn, seed, control$ncores, pi_floor, nub, nue, adjE, capture
  )
 }
 .sbayesrc_mcem_engine(
  state, B, E, A, gamma, alpha, sigmaSqAlpha_init,
  intercept_prior_resolved, inner_function,
  control$inner_sweeps, control$inner_burn,
  control$final_sweeps, control$final_burn,
  control$damping, control$tol_alpha, control$tol_prior,
  control$min_outer, control$max_outer, pi_floor, control$seed,
  "csr_reference", FALSE, FALSE, return_responsibilities, verbose,
  .objective_diagnostics,
  .diagnostic_responsibility_iterations
 )
}

.stblr_mcem_sbayesrc_block_eigen <- function(
  stats, Glist, annotation, block_start,
  B, E, ssb_prior, sse_prior, gamma, alpha_init,
  sigmaSqAlpha_init, intercept_prior_resolved,
  b_init = NULL, comp_init = NULL,
  representation = "low_rank", eigen_prop = 0.995,
  eigen_filter = "hard_truncate", eigen_tau = 0.01, eigen_eta = 0,
  residual_policy = "gctb_block", block_ve_mode = "allMixVe",
  resam_thresh = 1.1, minimum_ve_ratio = 0.7,
  low_rank_residual_rebuild_every = 100L,
  updateB = FALSE, updateE = FALSE,
  inner_sweeps = 1000L, inner_burn = 300L,
  final_sweeps = inner_sweeps, final_burn = inner_burn,
  damping = 0.5, tol_alpha = 1e-3, tol_prior = 1e-3,
  min_outer = 3L, max_outer = 50L,
  pi_floor = 1e-12, nub = 4, nue = 4,
  ncores = 1L, seed = 10L, return_responsibilities = TRUE,
  verbose = FALSE, .objective_diagnostics = FALSE,
  .diagnostic_responsibility_iterations = integer()
) {
 A <- as.matrix(annotation)
 alpha <- .sbayesrc_validate_alpha(as.matrix(alpha_init), gamma)
 gamma <- .sbayesrc_validate_gamma(gamma)
 marker_count <- nrow(A)
 if (marker_count < 1L || ncol(A) < 2L || nrow(alpha) != ncol(A)) {
  stop("Block MCEM annotation and alpha dimensions are inconsistent.")
 }
 sigmaSqAlpha_init <- as.numeric(sigmaSqAlpha_init)
 intercept <- .sbayesrc_mcem_intercept_prior(
  intercept_prior_resolved, ncol(alpha)
 )
 if (length(sigmaSqAlpha_init) != ncol(alpha) ||
     any(!is.finite(sigmaSqAlpha_init)) || any(sigmaSqAlpha_init <= 0)) {
  stop("sigmaSqAlpha_init must provide a fixed positive variance per stick.")
 }
 control <- .sbayesrc_mcem_validate_controls(
  inner_sweeps, inner_burn, final_sweeps, final_burn,
  min_outer, max_outer, ncores, seed, damping, tol_alpha, tol_prior
 )
 if (!identical(residual_policy, "gctb_block")) {
  stop("Phase-5B block MCEM requires residual_policy = 'gctb_block'.")
 }
 resolved_mode <- if (isTRUE(updateE)) "allMixVe" else "fixVe"
 if (!identical(block_ve_mode, resolved_mode)) {
  stop("block_ve_mode must be '", resolved_mode,
       "' for the requested Phase-5B E-update mode.")
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
   alpha_init = alpha, sigmaSqAlpha_init = sigmaSqAlpha_init,
   annotation_intercept_prior = intercept_spec,
   pi_floor = pi_floor, nub = nub, nue = nue,
   .diagnostic_updateSigmaSqAlpha = FALSE,
   .information_diagnostics = capture, .return_raw = TRUE
  )
 }
 .sbayesrc_mcem_engine(
  state, as.matrix(B), as.matrix(E), A, gamma, alpha,
  sigmaSqAlpha_init, intercept_prior_resolved, inner_function,
  control$inner_sweeps, control$inner_burn,
  control$final_sweeps, control$final_burn,
  control$damping, control$tol_alpha, control$tol_prior,
  control$min_outer, control$max_outer, pi_floor, control$seed,
  "block_eigen", updateB, updateE, return_responsibilities, verbose,
  .objective_diagnostics,
  .diagnostic_responsibility_iterations
 )
}
