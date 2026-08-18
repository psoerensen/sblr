# Exact state enumeration and low-dimensional numerical covariance references.

mt_vec_y <- function(Y) as.numeric(Y)

mt_state_marginal <- function(X, Y, Ve, state, Vb, gamma = NULL, q = NULL) {
  dat <- mt_validate_data(X, Y, Ve)
  state <- as.integer(state)
  if (length(state) != dat$M) stop("state length must equal M.", call. = FALSE)
  active <- which(state > 0L)
  residual_covariance <- kronecker(Ve, diag(dat$N))
  if (!length(active)) {
    return(list(log_density = mt_log_mvn_zero(mt_vec_y(Y), residual_covariance),
                mean = matrix(0, dat$M, dat$T), covariance = NULL))
  }
  H <- do.call(cbind, lapply(active, function(j) kronecker(diag(dat$T), X[, j])))
  scales <- if (is.null(gamma)) {
    rep(1, length(active))
  } else {
    q <- mt_validate_marker_scale(q, dat$M)
    gamma[state[active] + 1L] * q[active]
  }
  prior <- matrix(0, length(active) * dat$T, length(active) * dat$T)
  for (active_index in seq_along(active)) {
    idx <- ((active_index - 1L) * dat$T + 1L):(active_index * dat$T)
    prior[idx, idx] <- scales[active_index] * Vb
  }
  marginal <- residual_covariance + H %*% prior %*% t(H)
  precision <- solve(prior) + crossprod(H, solve(residual_covariance, H))
  posterior_covariance <- solve(precision)
  posterior_mean <- as.numeric(posterior_covariance %*%
                                 crossprod(H, solve(residual_covariance, mt_vec_y(Y))))
  mean_matrix <- matrix(0, dat$M, dat$T)
  for (active_index in seq_along(active)) {
    idx <- ((active_index - 1L) * dat$T + 1L):(active_index * dat$T)
    mean_matrix[active[active_index], ] <- posterior_mean[idx]
  }
  list(log_density = mt_log_mvn_zero(mt_vec_y(Y), marginal),
       mean = mean_matrix, covariance = posterior_covariance)
}

mt_all_states <- function(values, M) {
  as.matrix(expand.grid(rep(list(values), M), KEEP.OUT.ATTRS = FALSE))
}

mt_exact_bayesc_known_vb <- function(X, Y, Ve, Vb, pi) {
  dat <- mt_validate_data(X, Y, Ve); mt_assert_spd(Vb, "Vb")
  states <- mt_all_states(0:1, dat$M)
  log_weight <- numeric(nrow(states)); conditional_mean <- vector("list", nrow(states))
  for (s in seq_len(nrow(states))) {
    ref <- mt_state_marginal(X, Y, Ve, states[s, ], Vb)
    log_weight[s] <- ref$log_density + sum(stats::dbinom(states[s, ], 1L, pi, log = TRUE))
    conditional_mean[[s]] <- ref$mean
  }
  probability <- exp(log_weight - max(log_weight)); probability <- probability / sum(probability)
  alpha_mean <- Reduce(`+`, Map(`*`, conditional_mean, probability))
  list(states = states, state_probability = probability,
       pip = colSums(states * probability), alpha_mean = alpha_mean,
       fitted_mean = X %*% alpha_mean, log_normalizer = max(log_weight) + log(sum(exp(log_weight - max(log_weight)))))
}

mt_exact_bayesr_known_vb <- function(X, Y, Ve, Vb,
                                     component_probability, gamma, q = NULL) {
  dat <- mt_validate_data(X, Y, Ve); mt_assert_spd(Vb, "Vb")
  q <- mt_validate_marker_scale(q, dat$M)
  if (!is.numeric(gamma) || length(gamma) < 2L || gamma[1L] != 0 ||
      any(!is.finite(gamma)) || any(gamma[-1L] <= 0) ||
      !is.numeric(component_probability) ||
      length(component_probability) != length(gamma) ||
      any(!is.finite(component_probability)) || any(component_probability <= 0)) {
    stop("invalid BayesR probabilities or component multipliers.", call. = FALSE)
  }
  component_probability <- component_probability / sum(component_probability)
  states <- mt_all_states(0:(length(gamma) - 1L), dat$M)
  log_weight <- numeric(nrow(states)); conditional_mean <- vector("list", nrow(states))
  for (s in seq_len(nrow(states))) {
    ref <- mt_state_marginal(X, Y, Ve, states[s, ], Vb, gamma, q)
    log_weight[s] <- ref$log_density + sum(log(component_probability[states[s, ] + 1L]))
    conditional_mean[[s]] <- ref$mean
  }
  probability <- exp(log_weight - max(log_weight)); probability <- probability / sum(probability)
  alpha_mean <- Reduce(`+`, Map(`*`, conditional_mean, probability))
  component_marginal <- matrix(0, dat$M, length(gamma))
  for (j in seq_len(dat$M)) for (k in 0:(length(gamma) - 1L)) {
    component_marginal[j, k + 1L] <- sum(probability[states[, j] == k])
  }
  list(states = states, state_probability = probability,
       component_probability = component_marginal,
       pip = 1 - component_marginal[, 1L], alpha_mean = alpha_mean,
       fitted_mean = X %*% alpha_mean, q = q)
}

mt_vb_grid_reference <- function(X, Y, Ve, pi, prior_df, prior_scale,
                                 grid_n = 27L, span_log_sd = 1.5,
                                 span_z = 2.2) {
  dat <- mt_validate_data(X, Y, Ve)
  if (dat$T != 2L || dat$M > 2L) stop("grid reference requires T=2 and M<=2.", call. = FALSE)
  prior_mean <- mt_iw_mean(prior_df, prior_scale)
  center <- log(sqrt(diag(prior_mean)))
  a <- seq(center[1L] - span_log_sd, center[1L] + span_log_sd, length.out = grid_n)
  b <- seq(center[2L] - span_log_sd, center[2L] + span_log_sd, length.out = grid_n)
  z <- seq(-span_z, span_z, length.out = grid_n)
  states <- mt_all_states(0:1, dat$M)
  sum_weight <- 0; sum_V <- matrix(0, 2L, 2L); sum_V2 <- matrix(0, 2L, 2L)
  records <- vector("list", grid_n^3); log_weight <- numeric(grid_n^3); q <- 0L
  for (ai in seq_along(a)) for (bi in seq_along(b)) for (zi in seq_along(z)) {
    q <- q + 1L; s1 <- exp(a[ai]); s2 <- exp(b[bi]); rho <- tanh(z[zi])
    V <- matrix(c(s1^2, rho * s1 * s2, rho * s1 * s2, s2^2), 2L, 2L)
    state_log <- numeric(nrow(states))
    for (s in seq_len(nrow(states))) {
      state_log[s] <- mt_state_marginal(X, Y, Ve, states[s, ], V)$log_density +
        sum(stats::dbinom(states[s, ], 1L, pi, log = TRUE))
    }
    marginal_log <- max(state_log) + log(sum(exp(state_log - max(state_log))))
    log_jacobian <- log(4) + 3 * a[ai] + 3 * b[bi] + log1p(-rho^2)
    log_weight[q] <- mt_log_iwishart(V, prior_df, prior_scale) + marginal_log + log_jacobian
    records[[q]] <- V
  }
  weight <- exp(log_weight - max(log_weight)); weight <- weight / sum(weight)
  for (i in seq_along(weight)) {
    sum_V <- sum_V + weight[i] * records[[i]]
    sum_V2 <- sum_V2 + weight[i] * records[[i]]^2
  }
  boundary <- array(weight, c(grid_n, grid_n, grid_n))
  boundary_mass <- sum(boundary[c(1L, grid_n), , ]) +
    sum(boundary[, c(1L, grid_n), ]) + sum(boundary[, , c(1L, grid_n)])
  list(mean = sum_V, variance = sum_V2 - sum_V^2,
       boundary_mass_upper_bound = boundary_mass,
       grid = c(log_sd1 = grid_n, log_sd2 = grid_n, fisher_z = grid_n))
}
