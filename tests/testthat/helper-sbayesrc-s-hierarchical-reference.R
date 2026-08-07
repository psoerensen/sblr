# Independent test-only helpers for the SBayesRC-S Phase-3 hierarchy.
# These functions extend the validated Phase-1/2 references without calling
# production SBayesRC code or exposing a package API.

.sbs3_log_beta_binomial_prior <- function(delta, a_pi, b_pi) {
  included <- sum(delta)
  lbeta(a_pi + included, b_pi + length(delta) - included) -
    lbeta(a_pi, b_pi)
}

.sbs3_collapsed_log_odds <- function(delta, annotation, a_pi, b_pi) {
  included_other <- sum(delta) - delta[annotation]
  log(a_pi + included_other) -
    log(b_pi + length(delta) - 1L - included_other)
}

.sbs3_draw_pi <- function(delta, a_pi, b_pi) {
  stats::rbeta(1L, a_pi + sum(delta), b_pi + length(delta) - sum(delta))
}

.sbs3_exact_pi_posterior <- function(z, annotation, eligible, tau2,
                                     a_pi, b_pi) {
  annotation <- as.matrix(annotation)
  states <- .sbs_model_states(ncol(annotation))
  log_weight <- numeric(nrow(states))
  for (model in seq_len(nrow(states))) {
    selected <- states[model, ] == 1L
    log_weight[model] <- .sbs3_log_beta_binomial_prior(
      states[model, ], a_pi, b_pi
    )
    for (stick in seq_along(z)) {
      rows <- eligible[[stick]]
      log_weight[model] <- log_weight[model] + .sbs_stick_model(
        z[[stick]], annotation[rows, , drop = FALSE], selected, tau2[stick]
      )$log_marginal
    }
  }
  model_probability <- exp(log_weight - .sbs_log_sum_exp(log_weight))
  names(model_probability) <- rownames(states)
  annotation_pip <- drop(crossprod(model_probability, states))
  names(annotation_pip) <- colnames(annotation)
  list(
    states = states,
    log_weight = log_weight,
    model_probability = model_probability,
    annotation_pip = annotation_pip,
    expected_included = sum(annotation_pip)
  )
}

.sbs3_log_ig <- function(tau2, shape, scale) {
  ifelse(
    is.finite(tau2) & tau2 > 0,
    shape * log(scale) - lgamma(shape) - (shape + 1) * log(tau2) -
      scale / tau2,
    -Inf
  )
}

.sbs3_draw_tau2 <- function(alpha, delta, shape, scale) {
  selected <- delta == 1L
  posterior_shape <- shape + sum(selected) / 2
  posterior_scale <- scale + 0.5 * sum(alpha[selected]^2)
  1 / stats::rgamma(1L, shape = posterior_shape, rate = posterior_scale)
}

.sbs3_tau_conditional <- function(alpha, delta, shape, scale) {
  selected <- delta == 1L
  list(
    shape = shape + sum(selected) / 2,
    scale = scale + 0.5 * sum(alpha[selected]^2)
  )
}

.sbs3_bfdr <- function(annotation_pip, threshold) {
  selected <- which(annotation_pip >= threshold)
  if (!length(selected)) return(NA_real_)
  sum(1 - annotation_pip[selected]) / length(selected)
}

.sbs3_run_chain <- function(annotation, eligible, iterations, burn,
                            initial_delta, pi_a, tau2,
                            outcome = NULL, fixed_z = NULL,
                            learn_pi = FALSE, collapsed_pi = FALSE,
                            a_pi = 1, b_pi = 1,
                            learn_tau = FALSE, a_tau = 3, b_tau = 1.6,
                            fixed_delta = NULL,
                            initial_alpha = NULL,
                            initial_intercept = NULL,
                            method = c("primary", "direct")) {
  method <- match.arg(method)
  annotation <- as.matrix(annotation)
  annotation_count <- ncol(annotation)
  stick_count <- length(eligible)
  stopifnot(
    iterations > burn, burn >= 0L,
    length(tau2) == stick_count, all(is.finite(tau2)), all(tau2 > 0),
    pi_a > 0, pi_a < 1, a_pi > 0, b_pi > 0,
    a_tau > 0, b_tau > 0,
    xor(is.null(outcome), is.null(fixed_z))
  )
  if (!is.null(outcome)) {
    .sbs2_validate_hierarchy(outcome, annotation, eligible)
  } else {
    stopifnot(length(fixed_z) == stick_count)
    for (stick in seq_len(stick_count)) {
      stopifnot(length(fixed_z[[stick]]) == length(eligible[[stick]]))
    }
  }

  delta <- as.integer(if (is.null(fixed_delta)) initial_delta else fixed_delta)
  stopifnot(length(delta) == annotation_count, all(delta %in% c(0L, 1L)))
  alpha <- if (is.null(initial_alpha)) {
    matrix(0, annotation_count, stick_count)
  } else {
    matrix(as.numeric(initial_alpha), annotation_count, stick_count)
  }
  alpha[delta == 0L, ] <- 0
  intercept <- if (!is.null(initial_intercept)) {
    as.numeric(initial_intercept)
  } else if (!is.null(outcome)) {
    vapply(outcome, function(d) {
      stats::qnorm((sum(d) + 0.5) / (length(d) + 1))
    }, numeric(1L))
  } else {
    vapply(fixed_z, mean, numeric(1L))
  }
  current_pi <- pi_a
  current_tau <- as.numeric(tau2)

  retained <- iterations - burn
  delta_draws <- matrix(NA_integer_, retained, annotation_count)
  alpha_draws <- array(NA_real_, c(retained, annotation_count, stick_count))
  intercept_draws <- matrix(NA_real_, retained, stick_count)
  pi_a_draws <- numeric(retained)
  tau2_draws <- matrix(NA_real_, retained, stick_count)
  included_draws <- integer(retained)
  q_sum <- matrix(0, nrow(annotation), stick_count)
  component_sum <- matrix(0, nrow(annotation), stick_count + 1L)
  kept <- 0L

  for (iteration in seq_len(iterations)) {
    latent <- if (!is.null(fixed_z)) fixed_z else vector("list", stick_count)
    if (is.null(fixed_z)) {
      for (stick in seq_len(stick_count)) {
        rows <- eligible[[stick]]
        eta <- intercept[stick] + drop(
          annotation[rows, , drop = FALSE] %*% alpha[, stick]
        )
        latent[[stick]] <- .sbs2_rtrunc_probit(eta, outcome[[stick]])
      }
    }

    if (is.null(fixed_delta)) {
      for (j in seq_len(annotation_count)) {
        log_odds <- if (learn_pi && collapsed_pi) {
          .sbs3_collapsed_log_odds(delta, j, a_pi, b_pi)
        } else {
          stats::qlogis(current_pi)
        }
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
          precision <- sum(x * x) + 1 / current_tau[stick]
          score <- sum(x * residual)
          moments[[stick]] <- c(
            mean = score / precision,
            variance = 1 / precision
          )
          log_odds <- log_odds + if (method == "primary") {
            .sbs_log_bf(x, residual, current_tau[stick])
          } else {
            .sbs2_log_bf_direct(x, residual, current_tau[stick])
          }
        }
        delta[j] <- stats::rbinom(1L, 1L, stats::plogis(log_odds))
        if (delta[j] == 1L) {
          for (stick in seq_len(stick_count)) {
            alpha[j, stick] <- stats::rnorm(
              1L, moments[[stick]]["mean"],
              sqrt(moments[[stick]]["variance"])
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
        current_tau[stick], method
      )
      intercept[stick] <- draw$intercept
      alpha[, stick] <- 0
      if (any(selected)) alpha[selected, stick] <- draw$slopes
    }

    if (learn_pi) current_pi <- .sbs3_draw_pi(delta, a_pi, b_pi)
    if (learn_tau) {
      for (stick in seq_len(stick_count)) {
        current_tau[stick] <- .sbs3_draw_tau2(
          alpha[, stick], delta, a_tau, b_tau
        )
      }
    }

    if (iteration > burn) {
      kept <- kept + 1L
      delta_draws[kept, ] <- delta
      alpha_draws[kept, , ] <- alpha
      intercept_draws[kept, ] <- intercept
      pi_a_draws[kept] <- current_pi
      tau2_draws[kept, ] <- current_tau
      included_draws[kept] <- sum(delta)
      q <- stats::pnorm(sweep(annotation %*% alpha, 2L, intercept, `+`))
      component <- .sbs_component_probability(q)
      if (any(!is.finite(q)) || any(q < 0) || any(q > 1) ||
          any(!is.finite(component)) ||
          max(abs(rowSums(component) - 1)) > 1e-12) {
        stop("Phase-3 probability contract failed")
      }
      q_sum <- q_sum + q
      component_sum <- component_sum + component
    }
  }

  list(
    delta_draws = delta_draws,
    alpha_draws = alpha_draws,
    intercept_draws = intercept_draws,
    pi_a_draws = pi_a_draws,
    tau2_draws = tau2_draws,
    included_draws = included_draws,
    q_mean = q_sum / retained,
    pi_mean = component_sum / retained
  )
}

.sbs3_pool <- function(chains, annotation_names = NULL) {
  pooled <- .sbs_pool_chains(chains, annotation_names)
  pooled$pi_a <- unlist(lapply(chains, `[[`, "pi_a_draws"), use.names = FALSE)
  pooled$tau2 <- do.call(rbind, lapply(chains, `[[`, "tau2_draws"))
  pooled$included <- unlist(lapply(chains, `[[`, "included_draws"),
                            use.names = FALSE)
  pooled$chain_pip <- do.call(rbind, lapply(
    chains, function(x) colMeans(x$delta_draws)
  ))
  pooled$chain_pi_mean <- vapply(
    chains, function(x) mean(x$pi_a_draws), numeric(1L)
  )
  pooled$chain_tau_mean <- do.call(rbind, lapply(
    chains, function(x) colMeans(x$tau2_draws)
  ))
  pooled$switching <- lapply(chains, function(x) {
    .sbs2_switch_diagnostics(x$delta_draws)
  })
  pooled
}

.sbs3_tau_quadrature <- function(z, annotation, eligible, pi_a,
                                 a_tau, b_tau,
                                 lower = -12, upper = 12) {
  stopifnot(ncol(annotation) == 1L, length(z) == 1L)
  rows <- eligible[[1L]]
  log_m0 <- .sbs_stick_model(
    z[[1L]], annotation[rows, , drop = FALSE], FALSE, 1
  )$log_marginal
  log_kernel <- function(log_tau, moment = 0) {
    tau <- exp(log_tau)
    log_m1 <- .sbs_stick_model(
      z[[1L]], annotation[rows, , drop = FALSE], TRUE, tau
    )$log_marginal
    exp(log_m1 + .sbs3_log_ig(tau, a_tau, b_tau) +
          log_tau + moment * log_tau)
  }
  # Scale both integrals by a deterministic grid maximum.
  grid <- seq(lower, upper, length.out = 1001L)
  log_grid <- vapply(grid, function(u) {
    tau <- exp(u)
    .sbs_stick_model(
      z[[1L]], annotation[rows, , drop = FALSE], TRUE, tau
    )$log_marginal + .sbs3_log_ig(tau, a_tau, b_tau) + u
  }, numeric(1L))
  shift <- max(log_grid)
  integral <- function(power = 0) stats::integrate(
    function(u) {
      tau <- exp(u)
      log_m1 <- .sbs_stick_model(
        z[[1L]], annotation[rows, , drop = FALSE], TRUE, tau
      )$log_marginal
      exp(log_m1 + .sbs3_log_ig(tau, a_tau, b_tau) + u - shift) *
        tau^power
    }, lower, upper, rel.tol = 1e-9, subdivisions = 1000L
  )$value * exp(shift)
  marginal1 <- integral(0)
  posterior_delta <- pi_a * marginal1 /
    ((1 - pi_a) * exp(log_m0) + pi_a * marginal1)
  mean_tau_included <- integral(1) / marginal1
  list(
    pip = posterior_delta,
    mean_tau_given_inclusion = mean_tau_included,
    mean_tau = (1 - posterior_delta) * b_tau / (a_tau - 1) +
      posterior_delta * mean_tau_included,
    marginal_excluded = exp(log_m0),
    marginal_included = marginal1
  )
}

.sbs3_prior_predictive <- function(draws, annotation_count, a_pi, b_pi,
                                   a_tau, b_tau, stick_count = 3L) {
  pi_a <- stats::rbeta(draws, a_pi, b_pi)
  delta <- vapply(pi_a, function(p) {
    sum(stats::rbinom(annotation_count, 1L, p))
  }, integer(1L))
  tau2 <- matrix(
    1 / stats::rgamma(draws * stick_count, shape = a_tau, rate = b_tau),
    draws, stick_count
  )
  list(pi_a = pi_a, included = delta, tau2 = tau2)
}

.sbs3_moderate_fixture <- function(observations = 260L,
                                   annotation_count = 12L,
                                   seed = 20270401L,
                                   hierarchical = FALSE,
                                   pi_true = 0.25,
                                   tau_true = rep(0.8, 3L)) {
  stopifnot(annotation_count >= 10L, length(tau_true) == 3L)
  set.seed(seed)
  annotation <- matrix(
    stats::rnorm(observations * annotation_count),
    observations, annotation_count
  )
  annotation[, 1L] <- stats::rbinom(observations, 1L, 0.2)
  annotation[, 2L] <- as.numeric(scale(annotation[, 2L]))
  annotation[, 3L] <- as.numeric(scale(annotation[, 3L]))
  # Column four is a strong, but not identical, proxy for annotation two.
  annotation[, 4L] <- as.numeric(scale(
    0.88 * annotation[, 2L] + sqrt(1 - 0.88^2) * annotation[, 4L]
  ))
  if (annotation_count > 4L) {
    annotation[, 5:annotation_count] <- apply(
      annotation[, 5:annotation_count, drop = FALSE], 2L,
      function(x) as.numeric(scale(x))
    )
  }
  colnames(annotation) <- c(
    "enriched_binary", "continuous_signal", "secondary_signal",
    "continuous_proxy", paste0("null_", seq_len(annotation_count - 4L))
  )

  if (hierarchical) {
    delta <- stats::rbinom(annotation_count, 1L, pi_true)
    if (!any(delta)) delta[sample.int(annotation_count, 1L)] <- 1L
    alpha <- matrix(0, annotation_count, 3L)
    for (stick in seq_len(3L)) {
      alpha[delta == 1L, stick] <- stats::rnorm(
        sum(delta), 0, sqrt(tau_true[stick])
      )
    }
  } else {
    delta <- c(1L, 1L, 1L, rep(0L, annotation_count - 3L))
    alpha <- matrix(0, annotation_count, 3L)
    alpha[1L, ] <- c(0.8, 0.55, 0.35)
    alpha[2L, ] <- c(0.38, 0.28, 0.18)
    alpha[3L, ] <- c(-0.3, -0.22, -0.15)
  }

  eligible <- vector("list", 3L)
  latent_true <- vector("list", 3L)
  outcome <- vector("list", 3L)
  intercept <- numeric(3L)
  rows <- seq_len(observations)
  for (stick in seq_len(3L)) {
    eligible[[stick]] <- rows
    design <- annotation[rows, , drop = FALSE]
    intercept[stick] <- .sbs2_calibrate_intercept(
      design, alpha[, stick], 0.58
    )
    eta <- intercept[stick] + drop(design %*% alpha[, stick])
    latent_true[[stick]] <- eta + stats::rnorm(length(rows))
    outcome[[stick]] <- as.integer(latent_true[[stick]] > 0)
    if (stick < 3L) rows <- rows[outcome[[stick]] == 1L]
  }
  .sbs2_validate_hierarchy(outcome, annotation, eligible)
  list(
    annotation = annotation, eligible = eligible, outcome = outcome,
    latent_true = latent_true, true_delta = delta, true_alpha = alpha,
    true_intercept = intercept, pi_true = pi_true, tau_true = tau_true
  )
}
