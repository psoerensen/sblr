# Independent R-side Phase-4A parity helpers. The C++ hierarchy is always the
# implementation under test; these calculations use the validated R oracle.

.sbs4a_alpha_matrix <- function(intercept, alpha) {
  rbind(Intercept = as.numeric(intercept), as.matrix(alpha))
}

.sbs4a_math_reference <- function(annotation, eligible, latent, alpha,
                                  delta, pi_a, tau2,
                                  a_pi, b_pi, a_tau, b_tau,
                                  intercept_prior = NULL) {
  annotation <- as.matrix(annotation)
  selectable <- ncol(annotation)
  sticks <- length(eligible)
  if (is.null(intercept_prior)) intercept_prior <- .sbs_intercept_prior(sticks)
  s <- t <- log_bf <- conditional_mean <- conditional_variance <-
    matrix(0, selectable, sticks)
  inclusion_logit <- numeric(selectable)
  eta <- vector("list", sticks)
  for (stick in seq_len(sticks)) {
    rows <- eligible[[stick]]
    eta[[stick]] <- alpha[1L, stick] +
      drop(annotation[rows, , drop = FALSE] %*% alpha[-1L, stick])
  }
  for (j in seq_len(selectable)) {
    inclusion_logit[j] <- stats::qlogis(pi_a)
    for (stick in seq_len(sticks)) {
      rows <- eligible[[stick]]
      other <- setdiff(seq_len(selectable), j)
      residual <- latent[[stick]] - alpha[1L, stick]
      if (length(other)) {
        residual <- residual - drop(
          annotation[rows, other, drop = FALSE] %*% alpha[other + 1L, stick]
        )
      }
      x <- annotation[rows, j]
      s[j, stick] <- sum(x^2)
      t[j, stick] <- sum(x * residual)
      conditional_variance[j, stick] <- 1 /
        (s[j, stick] + 1 / tau2[stick])
      conditional_mean[j, stick] <- conditional_variance[j, stick] * t[j, stick]
      log_bf[j, stick] <- -0.5 * log1p(tau2[stick] * s[j, stick]) +
        0.5 * tau2[stick] * t[j, stick]^2 /
        (1 + tau2[stick] * s[j, stick])
      inclusion_logit[j] <- inclusion_logit[j] + log_bf[j, stick]
    }
  }
  beta_parameters <- c(
    a_pi + sum(delta), b_pi + length(delta) - sum(delta)
  )
  ig_parameters <- cbind(
    shape = rep(a_tau + sum(delta) / 2, sticks),
    scale = vapply(seq_len(sticks), function(stick) {
      b_tau + 0.5 * sum(alpha[-1L, stick][delta == 1L]^2)
    }, numeric(1L))
  )
  q <- stats::pnorm(sweep(annotation %*% alpha[-1L, , drop = FALSE],
                          2L, alpha[1L, ], `+`))
  intercept_conditional <- t(vapply(seq_len(sticks), function(stick) {
    rows <- eligible[[stick]]
    residual <- latent[[stick]]
    if (any(delta == 1L)) {
      residual <- residual - drop(
        annotation[rows, delta == 1L, drop = FALSE] %*%
          alpha[which(delta == 1L) + 1L, stick]
      )
    }
    precision <- length(rows) + intercept_prior$precision[stick]
    c(
      mean = (sum(residual) + intercept_prior$precision[stick] *
                intercept_prior$mean[stick]) / precision,
      variance = 1 / precision
    )
  }, numeric(2L)))
  list(
    eta = eta, s = s, t = t, log_bf = log_bf,
    inclusion_logit = inclusion_logit,
    inclusion_probability = stats::plogis(inclusion_logit),
    conditional_mean = conditional_mean,
    conditional_variance = conditional_variance,
    beta_parameters = beta_parameters,
    ig_parameters = ig_parameters,
    intercept_conditional = intercept_conditional,
    q = q,
    component_probability = .sbs_component_probability(q)
  )
}

.sbs4a_cpp_chain <- function(fixture, seed, iterations, burn, initial_delta,
                             fixed_delta = integer(), a_pi = 1, b_pi = 1,
                             a_tau = 3, b_tau = 1.6) {
  .st_bayesrc_selection_hierarchy(
    fixture$annotation, fixture$eligible, fixture$outcome,
    matrix(0, ncol(fixture$annotation) + 1L, length(fixture$eligible)),
    as.integer(initial_delta), 0.35, rep(0.8, length(fixture$eligible)),
    a_pi, b_pi, a_tau, b_tau,
    fixture$intercept_prior$native, 1e-12,
    iterations, burn, seed, as.integer(fixed_delta)
  )
}

.sbs4a_cpp_summary <- function(chains) {
  delta <- do.call(rbind, lapply(chains, `[[`, "delta_draws"))
  alpha <- array(NA_real_, c(
    sum(vapply(chains, function(x) dim(x$alpha_draws)[1L], integer(1L))),
    dim(chains[[1L]]$alpha_draws)[2:3]
  ))
  offset <- 0L
  for (chain in chains) {
    rows <- offset + seq_len(dim(chain$alpha_draws)[1L])
    alpha[rows, , ] <- chain$alpha_draws
    offset <- max(rows)
  }
  list(
    delta = delta,
    alpha = alpha,
    annotation_pip = colMeans(delta),
    alpha_mean = apply(alpha[, -1L, , drop = FALSE], c(2L, 3L), mean),
    intercept_mean = apply(alpha[, 1L, , drop = FALSE], 3L, mean),
    pi_a_mean = mean(unlist(lapply(chains, `[[`, "pi_A_draws"))),
    tau2_mean = colMeans(do.call(rbind, lapply(chains, `[[`, "tau2_draws"))),
    included_mean = mean(unlist(lapply(chains, `[[`, "included_draws"))),
    q_mean = Reduce(`+`, lapply(chains, `[[`, "q_mean")) / length(chains),
    component_probability_mean = Reduce(
      `+`, lapply(chains, `[[`, "component_probability_mean")
    ) / length(chains),
    chain_pip = do.call(rbind, lapply(chains, function(x) {
      colMeans(x$delta_draws)
    }))
  )
}
