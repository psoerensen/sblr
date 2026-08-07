# Development/reference implementation used for posterior-correctness and
# sampler-feasibility audits. Not a production sampler option and not part of
# the supported public API. These helpers implement the exact conditional
# particle-Gibbs block reference used for retained block-eigen SBayesRC audits.

.sbayesrc_log_sum_exp <- function(x) {
  maximum <- max(x)
  maximum + log(sum(exp(x - maximum)))
}

.sbayesrc_particle_marker_parameters <- function(score, diagonal, ve, vb,
                                                  gamma, probability) {
  ncomponent <- length(gamma)
  log_weight <- log(probability)
  mean <- variance <- numeric(ncomponent)
  for (component in seq_len(ncomponent)[gamma > 0]) {
    prior_variance <- vb * gamma[[component]]
    precision <- diagonal / ve + 1 / prior_variance
    variance[[component]] <- 1 / precision
    mean[[component]] <- variance[[component]] * score / ve
    log_weight[[component]] <- log_weight[[component]] + 0.5 * (
      log(variance[[component]] / prior_variance) +
        mean[[component]]^2 / variance[[component]])
  }
  list(log_weight = log_weight, mean = mean, variance = variance,
       log_normalizer = .sbayesrc_log_sum_exp(log_weight))
}

.sbayesrc_multinomial_resample <- function(weight, n) {
  cumulative <- cumsum(weight / sum(weight))
  vapply(stats::runif(n), function(value) which(value <= cumulative)[[1L]],
         integer(1L))
}

.sbayesrc_particle_block_step <- function(
    Q, w, component, beta, component_probability, gamma, vb, ve,
    particles = 16L, resampling_threshold = 0.5,
    retain_diagnostics = FALSE) {
  Q <- as.matrix(Q)
  w <- as.numeric(w)
  component_probability <- as.matrix(component_probability)
  marker_count <- ncol(Q)
  rank <- nrow(Q)
  particles <- as.integer(particles)
  if (length(w) != rank || length(component) != marker_count ||
      length(beta) != marker_count || nrow(component_probability) != marker_count ||
      ncol(component_probability) != length(gamma))
    stop("Particle-block dimensions are inconsistent.", call. = FALSE)
  if (particles < 2L || !is.finite(resampling_threshold) ||
      resampling_threshold <= 0 || resampling_threshold > 1 ||
      !is.finite(vb) || vb <= 0 || !is.finite(ve) || ve <= 0 ||
      any(!is.finite(Q)) || any(!is.finite(w)) ||
      any(!is.finite(component_probability)) ||
      any(component_probability <= 0) ||
      max(abs(rowSums(component_probability) - 1)) > 1e-10)
    stop("Particle-block inputs are invalid.", call. = FALSE)
  if (any(component < 0L | component >= length(gamma)) ||
      any(beta[component == 0L] != 0))
    stop("The retained conditional path is invalid.", call. = FALSE)

  residual <- matrix(w, nrow = rank, ncol = particles)
  component_path <- matrix(0L, particles, marker_count)
  beta_path <- matrix(0, particles, marker_count)
  log_weight <- rep(-log(particles), particles)
  ess <- numeric(marker_count)
  ancestor_diversity <- integer(marker_count)
  resampled <- logical(marker_count)

  for (marker in seq_len(marker_count)) {
    if (marker > 1L) {
      normalized <- exp(log_weight - .sbayesrc_log_sum_exp(log_weight))
      current_ess <- 1 / sum(normalized^2)
      if (current_ess < resampling_threshold * particles) {
        ancestor <- c(1L, .sbayesrc_multinomial_resample(
          normalized, particles - 1L))
        residual <- residual[, ancestor, drop = FALSE]
        component_path <- component_path[ancestor, , drop = FALSE]
        beta_path <- beta_path[ancestor, , drop = FALSE]
        log_weight[] <- -log(particles)
        resampled[[marker]] <- TRUE
        ancestor_diversity[[marker]] <- length(unique(ancestor))
      } else {
        ancestor_diversity[[marker]] <- particles
      }
    } else {
      ancestor_diversity[[marker]] <- particles
    }

    column <- Q[, marker]
    diagonal <- sum(column^2)
    next_log_weight <- numeric(particles)
    for (particle in seq_len(particles)) {
      score <- sum(column * residual[, particle])
      proposal <- .sbayesrc_particle_marker_parameters(
        score, diagonal, ve, vb, gamma,
        component_probability[marker, ])
      if (particle == 1L) {
        selected_component <- component[[marker]] + 1L
        selected_beta <- beta[[marker]]
      } else {
        probability <- exp(proposal$log_weight - proposal$log_normalizer)
        selected_component <- sample.int(length(gamma), 1L, prob = probability)
        selected_beta <- if (gamma[[selected_component]] == 0) 0 else
          stats::rnorm(1L, proposal$mean[[selected_component]],
                       sqrt(proposal$variance[[selected_component]]))
      }
      component_path[particle, marker] <- selected_component - 1L
      beta_path[particle, marker] <- selected_beta
      residual[, particle] <- residual[, particle] - column * selected_beta
      next_log_weight[[particle]] <- log_weight[[particle]] +
        proposal$log_normalizer
    }
    log_weight <- next_log_weight
    normalized <- exp(log_weight - .sbayesrc_log_sum_exp(log_weight))
    ess[[marker]] <- 1 / sum(normalized^2)
  }

  final_probability <- exp(log_weight - .sbayesrc_log_sum_exp(log_weight))
  selected <- sample.int(particles, 1L, prob = final_probability)
  result <- list(component = component_path[selected, ],
                 beta = beta_path[selected, ])
  if (retain_diagnostics) result$diagnostics <- list(
    particle_ess = ess, minimum_ess = min(ess), median_ess = stats::median(ess),
    resampling_count = sum(resampled), ancestor_diversity = ancestor_diversity,
    final_unique_paths = nrow(unique(cbind(component_path, beta_path))),
    selected_reference_path = selected == 1L,
    allocation_changes = sum(component_path[selected, ] != component),
    active_count_jump = sum(component_path[selected, ] > 0L) - sum(component > 0L))
  result
}

.sbayesrc_exact_block_allocation <- function(
    Q, w, component_probability, gamma, vb, ve) {
  Q <- as.matrix(Q)
  marker_count <- ncol(Q)
  states <- as.matrix(expand.grid(rep(list(seq_along(gamma) - 1L), marker_count)))
  log_weight <- numeric(nrow(states))
  for (state_index in seq_len(nrow(states))) {
    state <- states[state_index, ]
    value <- sum(log(component_probability[cbind(seq_len(marker_count), state + 1L)]))
    active <- which(gamma[state + 1L] > 0)
    if (length(active)) {
      prior_variance <- vb * gamma[state[active] + 1L]
      precision <- crossprod(Q[, active, drop = FALSE]) / ve +
        diag(1 / prior_variance, length(active))
      covariance <- solve(precision)
      score <- crossprod(Q[, active, drop = FALSE], w) / ve
      mean <- covariance %*% score
      value <- value + 0.5 * (determinant(covariance, logarithm = TRUE)$modulus -
        sum(log(prior_variance)) + crossprod(mean, precision %*% mean))
    }
    log_weight[[state_index]] <- value
  }
  probability <- exp(log_weight - .sbayesrc_log_sum_exp(log_weight))
  list(states = states, probability = probability,
       log_normalizer = .sbayesrc_log_sum_exp(log_weight),
       pip = colSums((states > 0L) * probability),
       active_count = vapply(0:marker_count, function(value)
         sum(probability[rowSums(states > 0L) == value]), numeric(1L)))
}
