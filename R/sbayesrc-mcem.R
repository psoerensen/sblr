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
  verbose = FALSE
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
 integer_scalar <- function(x, name, minimum = 1L) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
      x != as.integer(x) || x < minimum) {
   stop(name, " must be an integer >= ", minimum, ".")
  }
  as.integer(x)
 }
 inner_sweeps <- integer_scalar(inner_sweeps, "inner_sweeps")
 inner_burn <- integer_scalar(inner_burn, "inner_burn", 0L)
 final_sweeps <- integer_scalar(final_sweeps, "final_sweeps")
 final_burn <- integer_scalar(final_burn, "final_burn", 0L)
 min_outer <- integer_scalar(min_outer, "min_outer")
 max_outer <- integer_scalar(max_outer, "max_outer")
 ncores <- integer_scalar(ncores, "ncores")
 seed <- integer_scalar(seed, "seed", 0L)
 if (min_outer > max_outer) stop("min_outer cannot exceed max_outer.")
 for (item in list(damping = damping, tol_alpha = tol_alpha,
                   tol_prior = tol_prior)) {
  if (!is.numeric(item) || length(item) != 1L || !is.finite(item) || item <= 0) {
   stop("damping, tol_alpha, and tol_prior must be positive finite scalars.")
  }
 }
 if (damping > 1) stop("damping must be in (0, 1].")
 if (!is.numeric(adjE) || length(adjE) != 1L || !is.finite(adjE) ||
     adjE != 0) {
  stop("Phase-5A qualification requires adjE = 0 and fixed genomic variances.")
 }

 state <- list(
  b = lapply(b_init, as.numeric),
  component = lapply(comp_init, as.numeric),
  r = lapply(r_init, as.numeric)
 )
 prior <- .sbayesrc_mcem_component_prior(A, alpha, pi_floor)
 alpha_history <- array(
  NA_real_, c(nrow(alpha), ncol(alpha), max_outer + 1L),
  dimnames = c(dimnames(alpha), list(outer = 0:max_outer))
 )
 alpha_history[, , 1L] <- alpha
 history <- vector("list", max_outer)
 converged <- FALSE
 last_responsibility <- NULL
 completed <- 0L

 for (outer in seq_len(max_outer)) {
  inner <- .sbayesrc_mcem_inner_csr(
   wy, ww, yy, state, ld_prefix, B, E, ssb_prior, sse_prior,
   A, gamma, alpha, sigmaSqAlpha_init, intercept_prior_resolved, n,
   inner_sweeps, inner_burn, seed + outer - 1L, ncores, pi_floor,
   nub, nue, adjE, TRUE
  )
  chain <- inner$chains[[1L]][[1L]]
  responsibility <- chain$information_flow$rb_comp_prob
  if (is.null(responsibility)) {
   stop("The MCEM inner kernel did not return RB responsibilities.")
  }
  m_step <- .sbayesrc_mcem_m_step(
   A, responsibility, alpha, intercept_prior_resolved,
   sigmaSqAlpha_init, probability_floor = pi_floor
  )
  alpha_target <- m_step$alpha
  alpha_new <- (1 - damping) * alpha + damping * alpha_target
  prior_new <- .sbayesrc_mcem_component_prior(A, alpha_new, pi_floor)
  max_delta_alpha <- max(abs(alpha_new - alpha))
  max_delta_prior <- max(abs(prior_new - prior))
  completed <- outer
  history[[outer]] <- data.frame(
   outer = outer,
   max_delta_alpha = max_delta_alpha,
   max_delta_prior = max_delta_prior,
   prior_expected_active = sum(1 - prior_new[, 1L]),
   rb_expected_active = sum(1 - responsibility[, 1L]),
   m_step_max_convergence = max(m_step$convergence),
   stringsAsFactors = FALSE
  )
  alpha_history[, , outer + 1L] <- alpha_new
  state <- .sbayesrc_mcem_extract_state(inner, state)
  alpha <- alpha_new
  prior <- prior_new
  last_responsibility <- responsibility
  if (isTRUE(verbose)) {
   message(sprintf(
    "MCEM outer=%d max_dAlpha=%.6g max_dPrior=%.6g RB_active=%.3f",
    outer, max_delta_alpha, max_delta_prior,
    history[[outer]]$rb_expected_active
   ))
  }
  if (outer >= min_outer && max_delta_alpha < tol_alpha &&
      max_delta_prior < tol_prior) {
   converged <- TRUE
   break
  }
 }

 final <- .sbayesrc_mcem_inner_csr(
  wy, ww, yy, state, ld_prefix, B, E, ssb_prior, sse_prior,
  A, gamma, alpha, sigmaSqAlpha_init, intercept_prior_resolved, n,
  final_sweeps, final_burn, seed + max_outer + 10000L, ncores,
 pi_floor, nub, nue, adjE, isTRUE(return_responsibilities)
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
  method = "MCEM-SBayesRC",
  target = "observed_data_alpha_MAP_empirical_Bayes",
  converged = converged,
  n_outer = completed,
  alpha_map = alpha,
  component_prior = prior,
  history = list(summary = history_summary, alpha = history_alpha),
  damping = damping,
  tol_alpha = tol_alpha,
  tol_prior = tol_prior,
  min_outer = min_outer,
  max_outer = max_outer,
  inner_sweeps = inner_sweeps,
  inner_burn = inner_burn,
  final_sweeps = final_sweeps,
  final_burn = final_burn,
  sigmaSqAlpha_fixed = sigmaSqAlpha_init,
  genomic_hyperparameters_fixed = TRUE,
  e_step_responsibilities = if (isTRUE(return_responsibilities)) {
   last_responsibility
  } else NULL,
  final_responsibilities = if (isTRUE(return_responsibilities)) {
   final$chains[[1L]][[1L]]$information_flow$rb_comp_prob
  } else NULL
 )
 structure(
  list(genomic = final, mcem = mcem),
  class = c("sblr_mcem_sbayesrc_phase5a", "list")
 )
}
