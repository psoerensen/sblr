# Exact development references for particle-marginal alpha feasibility.
# These functions are internal and do not alter a production sampler.

.sbayesrc_particle_auxiliary <- function(particles, markers) {
  list(component = matrix(stats::rnorm(particles * markers), particles, markers),
       effect = matrix(stats::rnorm(particles * markers), particles, markers),
       selected_path = stats::rnorm(1L))
}

.sbayesrc_correlate_auxiliary <- function(auxiliary, rho) {
  if (!is.numeric(rho) || length(rho) != 1L || !is.finite(rho) ||
      abs(rho) >= 1)
    stop("rho must be a finite scalar strictly between -1 and 1.", call. = FALSE)
  list(component = rho * auxiliary$component +
         sqrt(1 - rho^2) * matrix(stats::rnorm(length(auxiliary$component)),
                                  nrow(auxiliary$component)),
       effect = rho * auxiliary$effect +
         sqrt(1 - rho^2) * matrix(stats::rnorm(length(auxiliary$effect)),
                                  nrow(auxiliary$effect)),
       selected_path = rho * auxiliary$selected_path +
         sqrt(1 - rho^2) * stats::rnorm(1L))
}

.sbayesrc_sis_block_likelihood <- function(
    Q, w, component_probability, gamma, vb, ve, particles = 16L,
    auxiliary = NULL, retain_paths = TRUE) {
  Q <- as.matrix(Q)
  w <- as.numeric(w)
  component_probability <- as.matrix(component_probability)
  particles <- as.integer(particles)
  markers <- ncol(Q)
  rank <- nrow(Q)
  components <- length(gamma)
  if (particles < 1L || length(w) != rank ||
      nrow(component_probability) != markers ||
      ncol(component_probability) != components ||
      any(component_probability <= 0) ||
      max(abs(rowSums(component_probability) - 1)) > 1e-10 ||
      !is.finite(vb) || vb <= 0 || !is.finite(ve) || ve <= 0)
    stop("Particle-marginal block inputs are invalid.", call. = FALSE)
  if (is.null(auxiliary)) auxiliary <- .sbayesrc_particle_auxiliary(particles, markers)
  if (!identical(dim(auxiliary$component), c(particles, markers)) ||
      !identical(dim(auxiliary$effect), c(particles, markers)) ||
      length(auxiliary$selected_path) != 1L ||
      any(!is.finite(unlist(auxiliary, use.names = FALSE))))
    stop("Particle auxiliary randomness has invalid dimensions or values.",
         call. = FALSE)

  residual <- matrix(w, rank, particles)
  log_weight <- numeric(particles)
  component_path <- if (retain_paths) matrix(0L, particles, markers) else NULL
  beta_path <- if (retain_paths) matrix(0, particles, markers) else NULL
  ess_prefix <- numeric(markers)
  diagonal <- colSums(Q^2)

  for (marker in seq_len(markers)) {
    column <- Q[, marker]
    score <- colSums(residual * column)
    log_component <- matrix(
      rep(log(component_probability[marker, ]), each = particles),
      particles, components)
    conditional_mean <- matrix(0, particles, components)
    conditional_sd <- numeric(components)
    for (component in which(gamma > 0)) {
      prior_variance <- vb * gamma[[component]]
      variance <- 1 / (diagonal[[marker]] / ve + 1 / prior_variance)
      mean <- variance * score / ve
      conditional_mean[, component] <- mean
      conditional_sd[[component]] <- sqrt(variance)
      log_component[, component] <- log_component[, component] + 0.5 *
        (log(variance / prior_variance) + mean^2 / variance)
    }
    row_maximum <- apply(log_component, 1L, max)
    component_mass <- exp(log_component - row_maximum)
    normalizer <- row_maximum + log(rowSums(component_mass))
    probability <- component_mass / rowSums(component_mass)
    uniforms <- stats::pnorm(auxiliary$component[, marker])
    cumulative <- t(apply(probability, 1L, cumsum))
    selected_component <- vapply(seq_len(particles), function(particle)
      which(uniforms[[particle]] <= cumulative[particle, ])[[1L]], integer(1L))
    selected_beta <- numeric(particles)
    active <- gamma[selected_component] > 0
    if (any(active)) {
      row <- which(active)
      selected_beta[row] <- conditional_mean[
        cbind(row, selected_component[row])] +
        conditional_sd[selected_component[row]] * auxiliary$effect[row, marker]
    }
    residual <- residual - tcrossprod(column, selected_beta)
    log_weight <- log_weight + normalizer
    normalized <- exp(log_weight - .sbayesrc_log_sum_exp(log_weight))
    ess_prefix[[marker]] <- 1 / sum(normalized^2)
    if (retain_paths) {
      component_path[, marker] <- selected_component - 1L
      beta_path[, marker] <- selected_beta
    }
  }

  log_estimate <- .sbayesrc_log_sum_exp(log_weight) - log(particles)
  normalized <- exp(log_weight - .sbayesrc_log_sum_exp(log_weight))
  result <- list(log_likelihood = log_estimate,
                 log_path_weight = log_weight,
                 final_ess = 1 / sum(normalized^2),
                 minimum_prefix_ess = min(ess_prefix),
                 median_prefix_ess = stats::median(ess_prefix),
                 auxiliary = auxiliary)
  if (retain_paths) {
    selected_uniform <- stats::pnorm(auxiliary$selected_path)
    selected <- which(selected_uniform <= cumsum(normalized))[[1L]]
    result$selected_path <- selected
    result$component <- component_path[selected, ]
    result$beta <- beta_path[selected, ]
  }
  result
}

.sbayesrc_alpha_log_prior <- function(alpha, intercept_mean, intercept_sd,
                                       sigma_sq_alpha) {
  alpha <- as.matrix(alpha)
  sticks <- ncol(alpha)
  if (length(intercept_mean) != sticks || length(intercept_sd) != sticks ||
      length(sigma_sq_alpha) != sticks || any(intercept_sd <= 0) ||
      any(sigma_sq_alpha <= 0))
    stop("Resolved alpha-prior dimensions are inconsistent.", call. = FALSE)
  sum(stats::dnorm(alpha[1L, ], intercept_mean, intercept_sd, log = TRUE)) +
    sum(vapply(seq_len(sticks), function(stick)
      sum(stats::dnorm(alpha[-1L, stick], 0, sqrt(sigma_sq_alpha[[stick]]),
                       log = TRUE)), numeric(1L)))
}

.sbayesrc_particle_marginal_log_ratio <- function(
    alpha, proposed_alpha, log_likelihood, proposed_log_likelihood,
    intercept_mean, intercept_sd, sigma_sq_alpha,
    log_q_reverse = 0, log_q_forward = 0) {
  .sbayesrc_alpha_log_prior(proposed_alpha, intercept_mean, intercept_sd,
                            sigma_sq_alpha) -
    .sbayesrc_alpha_log_prior(alpha, intercept_mean, intercept_sd,
                              sigma_sq_alpha) +
    proposed_log_likelihood - log_likelihood + log_q_reverse - log_q_forward
}
