# Independent test-only helpers for the observed-d SBayesRC-S Phase-2
# reference. They extend, but do not modify, the validated Phase-1 helpers.

.sbs2_rtrunc_probit <- function(eta, outcome) {
  eta <- as.numeric(eta)
  outcome <- as.integer(outcome)
  stopifnot(
    length(eta) == length(outcome), all(is.finite(eta)),
    all(outcome %in% c(0L, 1L))
  )
  uniform <- stats::runif(length(eta))
  standard <- numeric(length(eta))
  continued <- outcome == 1L
  if (any(continued)) {
    tail_mass <- stats::pnorm(eta[continued])
    survival <- uniform[continued] * tail_mass
    survival <- pmax(survival, .Machine$double.xmin)
    standard[continued] <- stats::qnorm(survival, lower.tail = FALSE)
  }
  stopped <- !continued
  if (any(stopped)) {
    lower_mass <- stats::pnorm(-eta[stopped])
    probability <- uniform[stopped] * lower_mass
    probability <- pmax(probability, .Machine$double.xmin)
    standard[stopped] <- stats::qnorm(probability)
  }
  latent <- eta + standard
  if (any(!is.finite(latent)) || any(latent[continued] <= 0) ||
      any(latent[stopped] > 0)) {
    stop("Truncated-normal draw violated the finite sign contract")
  }
  latent
}

.sbs2_log_bf_direct <- function(x, residual, tau2) {
  # Direct one-dimensional Gaussian integration, kept algebraically separate
  # from the rank-one log1p formula used by the primary Phase-1 update.
  likelihood_precision <- sum(x * x)
  score <- sum(x * residual)
  posterior_precision <- likelihood_precision + 1 / tau2
  -0.5 * log(tau2) - 0.5 * log(posterior_precision) +
    0.5 * score^2 / posterior_precision
}

.sbs2_validate_hierarchy <- function(outcome, annotation, eligible) {
  annotation <- as.matrix(annotation)
  stopifnot(length(outcome) == length(eligible), length(outcome) > 0L)
  expected <- seq_len(nrow(annotation))
  for (stick in seq_along(outcome)) {
    rows <- as.integer(eligible[[stick]])
    d <- as.integer(outcome[[stick]])
    stopifnot(
      identical(rows, expected), length(d) == length(rows),
      all(d %in% c(0L, 1L)), nrow(annotation[rows, , drop = FALSE]) == length(d)
    )
    if (stick < length(outcome)) expected <- rows[d == 1L]
  }
  invisible(TRUE)
}

.sbs2_draw_coefficients <- function(latent, annotation, selected, tau2,
                                    method = c("primary", "direct"),
                                    intercept_mean = 0,
                                    intercept_variance = 1) {
  method <- match.arg(method)
  posterior <- .sbs_stick_model(
    latent, annotation, selected, tau2, intercept_mean, intercept_variance
  )
  if (method == "primary") return(.sbs_draw_from_stick_posterior(posterior))
  chol_covariance <- chol(posterior$covariance)
  draw <- posterior$mean + drop(t(chol_covariance) %*%
                                  stats::rnorm(length(posterior$mean)))
  list(intercept = draw[1L], slopes = draw[-1L])
}

.sbs2_run_chain <- function(outcome, annotation, eligible, pi_a, tau2,
                            iterations, burn, initial_delta,
                            method = c("primary", "direct"),
                            fixed_delta = NULL, initial_alpha = NULL,
                            initial_intercept = NULL,
                            intercept_prior = NULL) {
  method <- match.arg(method)
  annotation <- as.matrix(annotation)
  .sbs2_validate_hierarchy(outcome, annotation, eligible)
  annotation_count <- ncol(annotation)
  stick_count <- length(outcome)
  if (is.null(intercept_prior)) {
    intercept_prior <- .sbs_intercept_prior(stick_count)
  }
  stopifnot(
    length(tau2) == stick_count, all(is.finite(tau2)), all(tau2 > 0),
    pi_a > 0, pi_a < 1, iterations > burn, burn >= 0
  )
  delta <- as.integer(if (is.null(fixed_delta)) initial_delta else fixed_delta)
  stopifnot(length(delta) == annotation_count, all(delta %in% c(0L, 1L)))
  alpha <- if (is.null(initial_alpha)) {
    matrix(0, annotation_count, stick_count)
  } else {
    matrix(as.numeric(initial_alpha), annotation_count, stick_count)
  }
  alpha[delta == 0L, ] <- 0
  intercept <- if (is.null(initial_intercept)) {
    vapply(seq_len(stick_count), function(stick) {
      d <- outcome[[stick]]
      if (!length(d)) intercept_prior$mean[stick] else
        stats::qnorm((sum(d) + 0.5) / (length(d) + 1))
    }, numeric(1L))
  } else {
    as.numeric(initial_intercept)
  }
  stopifnot(length(intercept) == stick_count)

  retained <- iterations - burn
  delta_draws <- matrix(NA_integer_, retained, annotation_count)
  alpha_draws <- array(NA_real_, c(retained, annotation_count, stick_count))
  intercept_draws <- matrix(NA_real_, retained, stick_count)
  q_sum <- matrix(0, nrow(annotation), stick_count)
  pi_sum <- matrix(0, nrow(annotation), stick_count + 1L)
  kept <- 0L
  for (iteration in seq_len(iterations)) {
    latent <- vector("list", stick_count)
    for (stick in seq_len(stick_count)) {
      rows <- eligible[[stick]]
      eta <- intercept[stick] + drop(annotation[rows, , drop = FALSE] %*%
                                       alpha[, stick])
      latent[[stick]] <- .sbs2_rtrunc_probit(eta, outcome[[stick]])
    }

    if (is.null(fixed_delta)) {
      for (j in seq_len(annotation_count)) {
        log_odds <- stats::qlogis(pi_a)
        moments <- vector("list", stick_count)
        for (stick in seq_len(stick_count)) {
          rows <- eligible[[stick]]
          other <- setdiff(seq_len(annotation_count), j)
          residual <- latent[[stick]] - intercept[stick]
          if (length(other)) {
            residual <- residual - drop(
              annotation[rows, other, drop = FALSE] %*% alpha[other, stick]
            )
          }
          x <- annotation[rows, j]
          precision <- sum(x * x) + 1 / tau2[stick]
          score <- sum(x * residual)
          moments[[stick]] <- c(mean = score / precision,
                                variance = 1 / precision)
          log_odds <- log_odds + if (method == "primary") {
            .sbs_log_bf(x, residual, tau2[stick])
          } else {
            .sbs2_log_bf_direct(x, residual, tau2[stick])
          }
        }
        delta[j] <- stats::rbinom(1L, 1L, stats::plogis(log_odds))
        if (delta[j] == 1L) {
          for (stick in seq_len(stick_count)) {
            alpha[j, stick] <- stats::rnorm(
              1L, moments[[stick]]["mean"], sqrt(moments[[stick]]["variance"])
            )
          }
        } else {
          alpha[j, ] <- 0
        }
      }
    }

    selected <- delta == 1L
    for (stick in seq_len(stick_count)) {
      rows <- eligible[[stick]]
      draw <- .sbs2_draw_coefficients(
        latent[[stick]], annotation[rows, , drop = FALSE], selected,
        tau2[stick], method, intercept_prior$mean[stick],
        intercept_prior$variance[stick]
      )
      intercept[stick] <- draw$intercept
      alpha[, stick] <- 0
      if (any(selected)) alpha[selected, stick] <- draw$slopes
    }

    if (iteration > burn) {
      kept <- kept + 1L
      delta_draws[kept, ] <- delta
      alpha_draws[kept, , ] <- alpha
      intercept_draws[kept, ] <- intercept
      q <- stats::pnorm(sweep(annotation %*% alpha, 2L, intercept, `+`))
      if (any(!is.finite(q)) || any(q < 0) || any(q > 1)) {
        stop("Continuation probabilities violated their contract")
      }
      mixture <- .sbs_component_probability(q)
      if (any(!is.finite(mixture)) ||
          max(abs(rowSums(mixture) - 1)) > 1e-12) {
        stop("Mixture probabilities violated their contract")
      }
      q_sum <- q_sum + q
      pi_sum <- pi_sum + mixture
    }
  }
  list(
    delta_draws = delta_draws,
    alpha_draws = alpha_draws,
    intercept_draws = intercept_draws,
    q_mean = q_sum / retained,
    pi_mean = pi_sum / retained
  )
}

.sbs2_calibrate_intercept <- function(annotation, slopes, target) {
  stats::uniroot(function(intercept) {
    mean(stats::pnorm(intercept + drop(annotation %*% slopes))) - target
  }, c(-10, 10), tol = 1e-12)$root
}

.sbs2_fixture <- function(observations = 420L, seed = 20261001L,
                          signal = c(1, 1, 0), signal_scale = 1,
                          target_continue = c(0.55, 0.55, 0.55)) {
  set.seed(seed)
  annotation <- cbind(
    enriched_binary = stats::rbinom(observations, 1L, 0.18),
    continuous_signal = as.numeric(scale(stats::rnorm(observations))),
    null_annotation = as.numeric(scale(stats::rnorm(observations)))
  )
  base_alpha <- rbind(
    enriched_binary = c(0.75, 0.5, 0.35),
    continuous_signal = c(0.32, 0.24, 0.18),
    null_annotation = c(0, 0, 0)
  )
  true_alpha <- base_alpha * signal * signal_scale
  eligible <- vector("list", 3L)
  latent_true <- vector("list", 3L)
  outcome <- vector("list", 3L)
  intercept <- numeric(3L)
  rows <- seq_len(observations)
  for (stick in seq_len(3L)) {
    eligible[[stick]] <- rows
    design <- annotation[rows, , drop = FALSE]
    intercept[stick] <- .sbs2_calibrate_intercept(
      design, true_alpha[, stick], target_continue[stick]
    )
    eta <- intercept[stick] + drop(design %*% true_alpha[, stick])
    latent_true[[stick]] <- eta + stats::rnorm(length(rows))
    outcome[[stick]] <- as.integer(latent_true[[stick]] > 0)
    if (stick < 3L) rows <- rows[outcome[[stick]] == 1L]
  }
  .sbs2_validate_hierarchy(outcome, annotation, eligible)
  list(
    annotation = annotation, eligible = eligible, outcome = outcome,
    latent_true = latent_true, pi_a = 0.35, tau2 = rep(0.8, 3L),
    intercept_prior = .sbs_intercept_prior(3L),
    true_delta = as.integer(signal != 0), true_alpha = true_alpha,
    true_intercept = intercept
  )
}

.sbs2_rhat <- function(chain_matrices) {
  sample_count <- nrow(chain_matrices[[1L]])
  chain_count <- length(chain_matrices)
  parameter_count <- ncol(chain_matrices[[1L]])
  out <- numeric(parameter_count)
  for (parameter in seq_len(parameter_count)) {
    values <- vapply(chain_matrices, function(x) mean(x[, parameter]), numeric(1L))
    variances <- vapply(chain_matrices, function(x) stats::var(x[, parameter]), numeric(1L))
    within <- mean(variances)
    between <- sample_count * stats::var(values)
    out[parameter] <- if (within > 0) {
      sqrt((((sample_count - 1) / sample_count) * within +
              between / sample_count) / within)
    } else {
      NA_real_
    }
  }
  out
}

.sbs2_switch_diagnostics <- function(draws) {
  draws <- as.matrix(draws)
  out <- lapply(seq_len(ncol(draws)), function(j) {
    difference <- diff(draws[, j])
    runs <- rle(draws[, j])$lengths
    c(zero_to_one = sum(difference == 1L),
      one_to_zero = sum(difference == -1L),
      mean_run_length = mean(runs))
  })
  do.call(rbind, out)
}

.sbs2_standard_continuous_chain <- function(outcome, annotation, eligible,
                                            tau2, iterations, burn,
                                            intercept_prior = NULL) {
  annotation <- as.matrix(annotation)
  .sbs2_validate_hierarchy(outcome, annotation, eligible)
  stick_count <- length(outcome)
  if (is.null(intercept_prior)) {
    intercept_prior <- .sbs_intercept_prior(stick_count)
  }
  coefficient_count <- ncol(annotation) + 1L
  coefficients <- matrix(0, coefficient_count, stick_count)
  retained <- iterations - burn
  draws <- array(NA_real_, c(retained, coefficient_count, stick_count))
  q_sum <- matrix(0, nrow(annotation), stick_count)
  pi_sum <- matrix(0, nrow(annotation), stick_count + 1L)
  kept <- 0L
  for (iteration in seq_len(iterations)) {
    for (stick in seq_len(stick_count)) {
      rows <- eligible[[stick]]
      design <- cbind(Intercept = 1, annotation[rows, , drop = FALSE])
      eta <- drop(design %*% coefficients[, stick])
      latent <- .sbs2_rtrunc_probit(eta, outcome[[stick]])
      posterior <- .sbs_stick_model(
        latent, annotation[rows, , drop = FALSE],
        rep(TRUE, ncol(annotation)), tau2[stick],
        intercept_prior$mean[stick], intercept_prior$variance[stick]
      )
      draw <- .sbs_draw_from_stick_posterior(posterior)
      coefficients[, stick] <- c(draw$intercept, draw$slopes)
    }
    if (iteration > burn) {
      kept <- kept + 1L
      draws[kept, , ] <- coefficients
      q <- stats::pnorm(sweep(
        annotation %*% coefficients[-1L, , drop = FALSE],
        2L, coefficients[1L, ], `+`
      ))
      q_sum <- q_sum + q
      pi_sum <- pi_sum + .sbs_component_probability(q)
    }
  }
  list(
    coefficient_draws = draws,
    q_mean = q_sum / retained,
    pi_mean = pi_sum / retained
  )
}
