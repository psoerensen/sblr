# Independent test-only mathematics for the fixed-z SBayesRC-S Phase-1
# reference model. These helpers do not call production SBayesRC code.

.sbs_log_sum_exp <- function(x) {
  maximum <- max(x)
  maximum + log(sum(exp(x - maximum)))
}

.sbs_intercept_prior <- function(stick_count,
                                 component_probability = NULL,
                                 mean = "initial_mixture", sd = 1) {
  if (is.null(component_probability)) {
    component_probability <- rep(1 / (stick_count + 1L), stick_count + 1L)
  }
  stopifnot(
    length(component_probability) == stick_count + 1L,
    all(is.finite(component_probability)), all(component_probability >= 0),
    sum(component_probability) > 0, length(sd) %in% c(1L, stick_count),
    all(is.finite(sd)), all(sd > 0)
  )
  probability <- component_probability / sum(component_probability)
  remaining <- rev(cumsum(rev(probability)))
  stick_probability <- remaining[-1L] / remaining[-length(remaining)]
  center <- if (is.character(mean)) {
    stopifnot(identical(mean, "initial_mixture"))
    stats::qnorm(stick_probability)
  } else {
    stopifnot(length(mean) %in% c(1L, stick_count), all(is.finite(mean)))
    rep(as.numeric(mean), length.out = stick_count)
  }
  scale <- rep(as.numeric(sd), length.out = stick_count)
  native <- rbind(
    type = rep(0, stick_count), mean = center, precision = scale^-2
  )
  list(
    mean = center, sd = scale, variance = scale^2,
    precision = scale^-2, component_probability = probability,
    stick_probability = stick_probability, native = native
  )
}

.sbs_component_probability <- function(q) {
  q <- as.matrix(q)
  probability <- matrix(0, nrow(q), ncol(q) + 1L)
  remaining <- rep(1, nrow(q))
  for (stick in seq_len(ncol(q))) {
    probability[, stick] <- remaining * (1 - q[, stick])
    remaining <- remaining * q[, stick]
  }
  probability[, ncol(probability)] <- remaining
  probability
}

.sbs_model_states <- function(annotation_count) {
  states <- as.matrix(expand.grid(rep(list(0:1), annotation_count)))
  storage.mode(states) <- "integer"
  colnames(states) <- paste0("annotation", seq_len(annotation_count))
  rownames(states) <- apply(states, 1L, paste0, collapse = "")
  states
}

.sbs_stick_model <- function(z, annotation, selected, tau2,
                             intercept_mean = 0,
                             intercept_variance = 1) {
  annotation <- as.matrix(annotation)
  selected <- as.logical(selected)
  design <- cbind(
    Intercept = rep(1, nrow(annotation)),
    annotation[, selected, drop = FALSE]
  )
  stopifnot(length(intercept_mean) == 1L, is.finite(intercept_mean),
            length(intercept_variance) == 1L,
            is.finite(intercept_variance), intercept_variance > 0)
  prior_precision <- c(intercept_variance^-1,
                       rep(tau2^-1, sum(selected)))
  prior_mean <- c(intercept_mean, rep(0, sum(selected)))
  precision <- crossprod(design) + diag(prior_precision, nrow = ncol(design))
  rhs <- drop(crossprod(design, z)) + prior_precision * prior_mean
  chol_precision <- chol(precision)
  mean <- backsolve(chol_precision, forwardsolve(t(chol_precision), rhs))
  covariance <- chol2inv(chol_precision)
  log_det_precision <- 2 * sum(log(diag(chol_precision)))
  log_det_prior_covariance <- log(intercept_variance) +
    sum(selected) * log(tau2)
  prior_quadratic <- sum(prior_precision * prior_mean^2)
  log_marginal <- -0.5 * log_det_prior_covariance -
    0.5 * log_det_precision + 0.5 * sum(rhs * mean) -
    0.5 * prior_quadratic
  list(
    log_marginal = log_marginal,
    mean = mean,
    covariance = covariance,
    precision = precision,
    chol_precision = chol_precision,
    selected = which(selected)
  )
}

.sbs_stick_log_marginal_dense <- function(z, annotation, selected, tau2,
                                           intercept_mean,
                                           intercept_variance) {
  annotation <- as.matrix(annotation)
  selected <- as.logical(selected)
  design <- cbind(
    Intercept = rep(1, nrow(annotation)),
    annotation[, selected, drop = FALSE]
  )
  prior_mean <- c(intercept_mean, rep(0, sum(selected)))
  prior_variance <- c(intercept_variance, rep(tau2, sum(selected)))
  covariance <- diag(length(z)) +
    design %*% diag(prior_variance, nrow = length(prior_variance)) %*%
    t(design)
  residual <- z - drop(design %*% prior_mean)
  if (!length(z)) return(0)
  chol_covariance <- chol(covariance)
  solved <- backsolve(
    chol_covariance, forwardsolve(t(chol_covariance), residual)
  )
  -sum(log(diag(chol_covariance))) +
    0.5 * (sum(z^2) - sum(residual * solved))
}

.sbs_intercept_prior_predictive <- function(prior) {
  probability_quantile <- function(p) {
    stats::pnorm(prior$mean + prior$sd * stats::qnorm(p))
  }
  cbind(
    median = probability_quantile(0.5),
    lower_50 = probability_quantile(0.25),
    upper_50 = probability_quantile(0.75),
    lower_95 = probability_quantile(0.025),
    upper_95 = probability_quantile(0.975),
    probability_below_001 = stats::pnorm(
      (stats::qnorm(0.01) - prior$mean) / prior$sd
    ),
    probability_above_099 = stats::pnorm(
      (prior$mean - stats::qnorm(0.99)) / prior$sd
    )
  )
}

.sbs_exact_posterior <- function(z, annotation, eligible, pi_a, tau2,
                                 intercept_prior = NULL) {
  annotation <- as.matrix(annotation)
  annotation_count <- ncol(annotation)
  stick_count <- length(z)
  states <- .sbs_model_states(annotation_count)
  stopifnot(
    length(eligible) == stick_count,
    length(tau2) == stick_count,
    pi_a > 0, pi_a < 1
  )
  if (is.null(intercept_prior)) {
    intercept_prior <- .sbs_intercept_prior(stick_count)
  }
  stopifnot(length(intercept_prior$mean) == stick_count,
            length(intercept_prior$variance) == stick_count)
  log_weight <- numeric(nrow(states))
  posterior <- vector("list", nrow(states))
  for (model in seq_len(nrow(states))) {
    selected <- states[model, ] == 1L
    stick_posterior <- vector("list", stick_count)
    log_weight[model] <- sum(stats::dbinom(states[model, ], 1, pi_a, log = TRUE))
    for (stick in seq_len(stick_count)) {
      rows <- eligible[[stick]]
      stick_posterior[[stick]] <- .sbs_stick_model(
        z[[stick]], annotation[rows, , drop = FALSE], selected, tau2[stick],
        intercept_prior$mean[stick], intercept_prior$variance[stick]
      )
      log_weight[model] <- log_weight[model] +
        stick_posterior[[stick]]$log_marginal
    }
    posterior[[model]] <- stick_posterior
  }
  model_probability <- exp(log_weight - .sbs_log_sum_exp(log_weight))
  names(model_probability) <- rownames(states)
  annotation_pip <- drop(crossprod(model_probability, states))
  names(annotation_pip) <- colnames(annotation)

  alpha_mean <- matrix(0, annotation_count, stick_count,
                       dimnames = list(colnames(annotation), paste0("stick", seq_len(stick_count))))
  alpha_included_numerator <- alpha_mean
  intercept_mean <- numeric(stick_count)
  q_by_model <- array(NA_real_, c(nrow(annotation), stick_count, nrow(states)))
  pi_by_model <- array(NA_real_, c(nrow(annotation), stick_count + 1L, nrow(states)))
  for (model in seq_len(nrow(states))) {
    selected <- states[model, ] == 1L
    for (stick in seq_len(stick_count)) {
      post <- posterior[[model]][[stick]]
      intercept_mean[stick] <- intercept_mean[stick] +
        model_probability[model] * post$mean[1L]
      if (any(selected)) {
        alpha_mean[selected, stick] <- alpha_mean[selected, stick] +
          model_probability[model] * post$mean[-1L]
        alpha_included_numerator[selected, stick] <-
          alpha_included_numerator[selected, stick] +
          model_probability[model] * post$mean[-1L]
      }
      prediction_design <- cbind(1, annotation[, selected, drop = FALSE])
      eta_mean <- drop(prediction_design %*% post$mean)
      eta_variance <- rowSums((prediction_design %*% post$covariance) * prediction_design)
      q_by_model[, stick, model] <- stats::pnorm(
        eta_mean / sqrt(1 + eta_variance)
      )
    }
    pi_by_model[, , model] <- .sbs_component_probability(q_by_model[, , model])
  }
  alpha_mean_given_inclusion <- alpha_included_numerator /
    matrix(annotation_pip, annotation_count, stick_count)
  q_mean <- apply(q_by_model * rep(model_probability, each = nrow(annotation) * stick_count),
                  c(1L, 2L), sum)
  pi_mean <- apply(pi_by_model * rep(model_probability, each = nrow(annotation) * (stick_count + 1L)),
                   c(1L, 2L), sum)
  list(
    states = states,
    model_probability = model_probability,
    annotation_pip = annotation_pip,
    intercept_mean = intercept_mean,
    alpha_mean = alpha_mean,
    alpha_mean_given_inclusion = alpha_mean_given_inclusion,
    q_mean = q_mean,
    pi_mean = pi_mean,
    stick_posterior = posterior,
    log_weight = log_weight
  )
}

.sbs_log_bf <- function(x, residual, tau2) {
  s <- sum(x * x)
  t <- sum(x * residual)
  -0.5 * log1p(tau2 * s) + 0.5 * tau2 * t^2 / (1 + tau2 * s)
}

.sbs_log_bf_direct <- function(x, residual, tau2) {
  covariance <- diag(length(x)) + tau2 * tcrossprod(x)
  chol_covariance <- chol(covariance)
  solved <- backsolve(chol_covariance, forwardsolve(t(chol_covariance), residual))
  -sum(log(diag(chol_covariance))) -
    0.5 * sum(residual * solved) + 0.5 * sum(residual^2)
}

.sbs_draw_stick_coefficients <- function(z, annotation, selected, tau2,
                                         intercept_mean = 0,
                                         intercept_variance = 1) {
  posterior <- .sbs_stick_model(
    z, annotation, selected, tau2, intercept_mean, intercept_variance
  )
  .sbs_draw_from_stick_posterior(posterior)
}

.sbs_draw_from_stick_posterior <- function(posterior) {
  draw <- posterior$mean + backsolve(
    posterior$chol_precision, stats::rnorm(length(posterior$mean))
  )
  list(intercept = draw[1L], slopes = draw[-1L])
}

.sbs_mcmc_chain <- function(z, annotation, eligible, pi_a, tau2,
                            iterations, burn, initial_delta,
                            intercept_prior = NULL) {
  annotation <- as.matrix(annotation)
  annotation_count <- ncol(annotation)
  stick_count <- length(z)
  if (is.null(intercept_prior)) {
    intercept_prior <- .sbs_intercept_prior(stick_count)
  }
  delta <- as.integer(initial_delta)
  alpha <- matrix(0, annotation_count, stick_count)
  intercept <- numeric(stick_count)
  retained <- iterations - burn
  delta_draws <- matrix(NA_integer_, retained, annotation_count)
  alpha_draws <- array(NA_real_, c(retained, annotation_count, stick_count))
  intercept_draws <- matrix(NA_real_, retained, stick_count)
  q_sum <- matrix(0, nrow(annotation), stick_count)
  pi_sum <- matrix(0, nrow(annotation), stick_count + 1L)
  states <- .sbs_model_states(annotation_count)
  posterior_cache <- lapply(seq_len(nrow(states)), function(model) {
    selected <- states[model, ] == 1L
    lapply(seq_len(stick_count), function(stick) {
      rows <- eligible[[stick]]
      .sbs_stick_model(
        z[[stick]], annotation[rows, , drop = FALSE], selected, tau2[stick],
        intercept_prior$mean[stick], intercept_prior$variance[stick]
      )
    })
  })
  kept <- 0L
  for (iteration in seq_len(iterations)) {
    for (j in seq_len(annotation_count)) {
      log_odds <- stats::qlogis(pi_a)
      moments <- vector("list", stick_count)
      for (stick in seq_len(stick_count)) {
        rows <- eligible[[stick]]
        other <- setdiff(seq_len(annotation_count), j)
        residual <- z[[stick]] - intercept[stick]
        if (length(other)) {
          residual <- residual - drop(annotation[rows, other, drop = FALSE] %*%
                                        alpha[other, stick])
        }
        x <- annotation[rows, j]
        s <- sum(x * x)
        t <- sum(x * residual)
        variance <- 1 / (s + tau2[stick]^-1)
        moments[[stick]] <- c(mean = variance * t, variance = variance)
        log_odds <- log_odds + .sbs_log_bf(x, residual, tau2[stick])
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
    selected <- delta == 1L
    model <- match(paste0(delta, collapse = ""), rownames(states))
    for (stick in seq_len(stick_count)) {
      draw <- .sbs_draw_from_stick_posterior(posterior_cache[[model]][[stick]])
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
      q_sum <- q_sum + q
      pi_sum <- pi_sum + .sbs_component_probability(q)
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

.sbs_fixture <- function() {
  set.seed(20260811)
  observations <- 180L
  annotation <- cbind(
    enriched_binary = as.numeric(seq_len(observations) %% 5L == 0L),
    continuous_signal = as.numeric(scale(sin(seq(0, 4 * pi, length.out = observations)))),
    null_annotation = as.numeric(scale(cos(seq(0, 7 * pi, length.out = observations))))
  )
  eligible <- list(
    seq_len(observations),
    seq.int(1L, observations, by = 2L),
    seq.int(1L, observations, by = 4L)
  )
  true_intercept <- c(-0.8, -0.25, 0.1)
  true_alpha <- rbind(
    enriched_binary = c(0.4, 0.3, 0.2),
    continuous_signal = c(0.22, 0.16, 0.1),
    null_annotation = c(0, 0, 0)
  )
  z <- lapply(seq_along(eligible), function(stick) {
    rows <- eligible[[stick]]
    true_intercept[stick] + drop(annotation[rows, ] %*% true_alpha[, stick]) +
      stats::rnorm(length(rows), 0, 1)
  })
  list(
    annotation = annotation, eligible = eligible, z = z,
    pi_a = 0.35, tau2 = c(0.8, 0.8, 0.8),
    intercept_prior = .sbs_intercept_prior(3L),
    true_delta = c(1L, 1L, 0L), true_alpha = true_alpha,
    true_intercept = true_intercept
  )
}

.sbs_pool_chains <- function(chains, annotation_names = NULL) {
  delta <- do.call(rbind, lapply(chains, `[[`, "delta_draws"))
  chain_sizes <- vapply(chains, function(x) dim(x$alpha_draws)[1L], integer(1L))
  alpha <- array(NA_real_, c(sum(chain_sizes), dim(chains[[1L]]$alpha_draws)[2:3]))
  offset <- 0L
  for (chain in seq_along(chains)) {
    rows <- offset + seq_len(chain_sizes[chain])
    alpha[rows, , ] <- chains[[chain]]$alpha_draws
    offset <- offset + chain_sizes[chain]
  }
  intercept <- do.call(rbind, lapply(chains, `[[`, "intercept_draws"))
  if (!is.null(annotation_names)) colnames(delta) <- annotation_names
  model_key <- apply(delta, 1L, paste0, collapse = "")
  model_probability <- table(factor(model_key, levels = rownames(
    .sbs_model_states(ncol(delta))
  ))) / nrow(delta)
  annotation_pip <- colMeans(delta)
  alpha_mean <- apply(alpha, c(2L, 3L), mean)
  alpha_mean_given_inclusion <- alpha_mean
  for (j in seq_len(ncol(delta))) {
    included <- delta[, j] == 1L
    alpha_mean_given_inclusion[j, ] <- apply(
      alpha[included, j, , drop = FALSE], 3L, mean
    )
  }
  list(
    delta = delta,
    alpha = alpha,
    intercept = intercept,
    model_probability = as.numeric(model_probability),
    annotation_pip = annotation_pip,
    alpha_mean = alpha_mean,
    alpha_mean_given_inclusion = alpha_mean_given_inclusion,
    intercept_mean = colMeans(intercept),
    q_mean = Reduce(`+`, lapply(chains, `[[`, "q_mean")) / length(chains),
    pi_mean = Reduce(`+`, lapply(chains, `[[`, "pi_mean")) / length(chains)
  )
}
