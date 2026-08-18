# Four-pattern MT-BayesC-Pi reference utilities for T = 2.

mt_pattern_space <- function() {
  out <- as.matrix(expand.grid(delta1 = 0:1, delta2 = 0:1,
                               KEEP.OUT.ATTRS = FALSE))
  storage.mode(out) <- "integer"
  rownames(out) <- apply(out, 1L, paste, collapse = "_")
  out
}

mt_validate_patterns <- function(patterns, probability = NULL) {
  patterns <- as.matrix(patterns)
  if (!is.numeric(patterns) || ncol(patterns) < 1L ||
      any(!patterns %in% 0:1) || anyDuplicated(as.data.frame(patterns)) ||
      sum(rowSums(patterns) == 0L) != 1L) {
    stop("patterns must be unique binary rows with exactly one null row.",
         call. = FALSE)
  }
  storage.mode(patterns) <- "integer"
  if (!is.null(probability)) {
    if (length(probability) != nrow(patterns) ||
        any(!is.finite(probability)) || any(probability <= 0) ||
        abs(sum(probability) - 1) > 1e-10) {
      stop("pattern probabilities must be positive, finite, and sum to one.",
           call. = FALSE)
    }
  }
  patterns
}

mt_rdirichlet <- function(shape) {
  if (!is.numeric(shape) || any(!is.finite(shape)) || any(shape <= 0)) {
    stop("Dirichlet shapes must be finite and strictly positive.", call. = FALSE)
  }
  draw <- rgamma(length(shape), shape = shape)
  if (any(!is.finite(draw)) || sum(draw) <= 0) {
    stop("Dirichlet draw was not finite and positive.", call. = FALSE)
  }
  draw / sum(draw)
}

mt_pattern_conditional <- function(partial, x, Vb, Ve, patterns,
                                   pattern_probability) {
  patterns <- mt_validate_patterns(patterns, pattern_probability)
  T <- ncol(patterns)
  if (!all(dim(partial) == c(length(x), T))) {
    stop("partial residual, marker, and pattern dimensions differ.", call. = FALSE)
  }
  mt_assert_spd(Vb, "Vb")
  mt_assert_spd(Ve, "Ve")
  Ei <- solve(Ve)
  score_full <- as.numeric(Ei %*% as.numeric(crossprod(x, partial)))
  xx <- sum(x^2)
  log_weight <- log(pattern_probability)
  active_mean <- active_covariance <- vector("list", nrow(patterns))

  for (s in seq_len(nrow(patterns))) {
    active <- which(patterns[s, ] == 1L)
    if (!length(active)) next
    prior <- Vb[active, active, drop = FALSE]
    prior_precision <- solve(prior)
    precision <- prior_precision + xx * Ei[active, active, drop = FALSE]
    covariance <- solve(precision)
    score <- score_full[active]
    mean <- as.numeric(covariance %*% score)
    log_weight[s] <- log_weight[s] + 0.5 * (
      mt_logdet_spd(prior_precision) - mt_logdet_spd(precision) +
        sum(score * mean)
    )
    active_mean[[s]] <- mean
    active_covariance[[s]] <- mt_symmetrize(covariance)
  }
  probability <- exp(log_weight - max(log_weight))
  probability <- probability / sum(probability)
  list(probability = probability, log_weight = log_weight,
       active_mean = active_mean, active_covariance = active_covariance)
}

mt_complete_latent <- function(pattern, active_value, Vb) {
  pattern <- as.integer(pattern)
  T <- length(pattern)
  active <- which(pattern == 1L)
  inactive <- which(pattern == 0L)
  if (!length(active)) return(mt_rmvnorm(rep(0, T), Vb))
  if (length(active_value) != length(active) || any(!is.finite(active_value))) {
    stop("active_value must align with the active pattern coordinates.",
         call. = FALSE)
  }
  out <- numeric(T)
  out[active] <- active_value
  if (length(inactive)) {
    Vaa <- Vb[active, active, drop = FALSE]
    Via <- Vb[inactive, active, drop = FALSE]
    mean_i <- as.numeric(Via %*% solve(Vaa, active_value))
    covariance_i <- Vb[inactive, inactive, drop = FALSE] -
      Via %*% solve(Vaa, t(Via))
    out[inactive] <- mt_rmvnorm(mean_i, mt_symmetrize(covariance_i))
  }
  out
}

mt_draw_pattern_effect <- function(pattern_index, conditional, patterns, Vb) {
  pattern <- patterns[pattern_index, ]
  active <- which(pattern == 1L)
  if (!length(active)) {
    beta <- mt_rmvnorm(rep(0, ncol(patterns)), Vb)
    return(list(beta = beta, alpha = numeric(ncol(patterns))))
  }
  active_value <- mt_rmvnorm(conditional$active_mean[[pattern_index]],
                             conditional$active_covariance[[pattern_index]])
  beta <- mt_complete_latent(pattern, active_value, Vb)
  list(beta = beta, alpha = pattern * beta)
}

mt_pattern_configuration_reference <- function(X, Y, Ve, Vb, patterns,
                                               pattern_probability = NULL,
                                               dirichlet_prior = NULL) {
  dat <- mt_validate_data(X, Y, Ve)
  patterns <- mt_validate_patterns(patterns)
  if (ncol(patterns) != dat$T || dat$M > 2L) {
    stop("enumeration requires aligned patterns and M <= 2.", call. = FALSE)
  }
  if (is.null(pattern_probability) == is.null(dirichlet_prior)) {
    stop("supply exactly one of fixed probabilities or a Dirichlet prior.",
         call. = FALSE)
  }
  if (!is.null(pattern_probability)) {
    mt_validate_patterns(patterns, pattern_probability)
  } else if (length(dirichlet_prior) != nrow(patterns) ||
             any(!is.finite(dirichlet_prior)) || any(dirichlet_prior <= 0)) {
    stop("Dirichlet prior must be positive and match the pattern count.",
         call. = FALSE)
  }
  configurations <- mt_all_states(seq_len(nrow(patterns)), dat$M)
  log_weight <- numeric(nrow(configurations))
  conditional_mean <- vector("list", nrow(configurations))
  for (q in seq_len(nrow(configurations))) {
    state <- configurations[q, ]
    masks <- patterns[state, , drop = FALSE]
    active <- which(rowSums(masks) > 0L)
    residual_covariance <- kronecker(Ve, diag(dat$N))
    alpha_mean <- matrix(0, dat$M, dat$T)
    if (length(active)) {
      H <- do.call(cbind, lapply(active, function(j) {
        kronecker(diag(dat$T), dat$X[, j]) %*% diag(masks[j, ])
      }))
      prior <- kronecker(diag(length(active)), Vb)
      marginal <- residual_covariance + H %*% prior %*% t(H)
      precision <- solve(prior) + crossprod(H, solve(residual_covariance, H))
      covariance <- solve(precision)
      mean <- as.numeric(covariance %*%
        crossprod(H, solve(residual_covariance, mt_vec_y(dat$Y))))
      for (a in seq_along(active)) {
        idx <- ((a - 1L) * dat$T + 1L):(a * dat$T)
        alpha_mean[active[a], ] <- masks[active[a], ] * mean[idx]
      }
      log_weight[q] <- mt_log_mvn_zero(mt_vec_y(dat$Y), marginal)
    } else {
      log_weight[q] <- mt_log_mvn_zero(mt_vec_y(dat$Y), residual_covariance)
    }
    counts <- tabulate(state, nbins = nrow(patterns))
    if (!is.null(pattern_probability)) {
      log_weight[q] <- log_weight[q] + sum(counts * log(pattern_probability))
    } else {
      a <- dirichlet_prior
      log_weight[q] <- log_weight[q] + lgamma(sum(a)) -
        lgamma(sum(a) + dat$M) + sum(lgamma(a + counts) - lgamma(a))
    }
    conditional_mean[[q]] <- alpha_mean
  }
  probability <- exp(log_weight - max(log_weight))
  probability <- probability / sum(probability)
  pattern_marginal <- matrix(0, dat$M, nrow(patterns))
  for (q in seq_len(nrow(configurations))) for (j in seq_len(dat$M)) {
    pattern_marginal[j, configurations[q, j]] <-
      pattern_marginal[j, configurations[q, j]] + probability[q]
  }
  alpha_mean <- Reduce(`+`, Map(`*`, conditional_mean, probability))
  pi_mean <- if (is.null(dirichlet_prior)) pattern_probability else {
    Reduce(`+`, lapply(seq_len(nrow(configurations)), function(q) {
      counts <- tabulate(configurations[q, ], nbins = nrow(patterns))
      probability[q] * (dirichlet_prior + counts) /
        (sum(dirichlet_prior) + dat$M)
    }))
  }
  list(configurations = configurations,
       configuration_probability = probability,
       pattern_probability = pattern_marginal,
       trait_pip = pattern_marginal %*% patterns,
       pleiotropic_probability = pattern_marginal[, which(rowSums(patterns) == dat$T)],
       alpha_mean = alpha_mean, fitted_mean = X %*% alpha_mean,
       pi_mean = pi_mean)
}

mt_regional_fixed_reference <- function(X, Y, Ve, patterns, region,
                                        region_levels, Vb_by_region,
                                        Pi_by_region) {
  dat <- mt_validate_data(X, Y, Ve)
  patterns <- mt_validate_patterns(patterns)
  region_index <- match(region, region_levels)
  if (dat$M > 2L || anyNA(region_index) || length(Vb_by_region) != length(region_levels) ||
      !all(dim(Pi_by_region) == c(length(region_levels), nrow(patterns)))) {
    stop("regional exact reference requires M <= 2 and aligned regional inputs.",
         call. = FALSE)
  }
  lapply(Vb_by_region, mt_assert_spd)
  if (any(Pi_by_region <= 0) || any(abs(rowSums(Pi_by_region) - 1) > 1e-10)) {
    stop("regional pattern probabilities must be positive and normalized.", call. = FALSE)
  }
  configurations <- mt_all_states(seq_len(nrow(patterns)), dat$M)
  log_weight <- numeric(nrow(configurations))
  for (q in seq_len(nrow(configurations))) {
    state <- configurations[q, ]
    masks <- patterns[state, , drop = FALSE]
    H <- do.call(cbind, lapply(seq_len(dat$M), function(j) {
      kronecker(diag(dat$T), dat$X[, j]) %*% diag(masks[j, ])
    }))
    prior <- matrix(0, dat$M * dat$T, dat$M * dat$T)
    for (j in seq_len(dat$M)) {
      idx <- ((j - 1L) * dat$T + 1L):(j * dat$T)
      prior[idx, idx] <- Vb_by_region[[region_index[j]]]
    }
    covariance <- kronecker(Ve, diag(dat$N)) + H %*% prior %*% t(H)
    log_weight[q] <- mt_log_mvn_zero(mt_vec_y(Y), covariance) +
      sum(log(Pi_by_region[cbind(region_index, state)]))
  }
  probability <- exp(log_weight - max(log_weight))
  probability <- probability / sum(probability)
  marginal <- matrix(0, dat$M, nrow(patterns))
  for (q in seq_len(nrow(configurations))) for (j in seq_len(dat$M)) {
    marginal[j, configurations[q, j]] <- marginal[j, configurations[q, j]] + probability[q]
  }
  list(configurations = configurations,
       configuration_probability = probability,
       pattern_probability = marginal,
       trait_pip = marginal %*% patterns)
}
