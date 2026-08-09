st_logvar_eta_q_internal <- sblr:::st_logvar_eta_q_internal
st_logvar_loglik_bayesc_internal <- sblr:::st_logvar_loglik_bayesc_internal
st_logvar_loglik_bayesr_internal <- sblr:::st_logvar_loglik_bayesr_internal
st_logvar_ess_fixture_internal <- sblr:::st_logvar_ess_fixture_internal

test_that("native log-variance eta and q match deterministic R fixtures", {
  X <- cbind(
    binary = c(-0.5, 0.5, -0.5, 0.5),
    continuous = c(-1.1618950039, -0.387298335, 0.387298335, 1.1618950039)
  )
  theta <- c(0.4, -0.2)
  expected_eta <- drop(X %*% theta)

  got <- st_logvar_eta_q_internal(X, theta)

  expect_equal(drop(got$eta), expected_eta, tolerance = 1e-14)
  expect_equal(drop(got$q), exp(expected_eta), tolerance = 1e-14)
  expect_equal(got$min_log_q, min(expected_eta), tolerance = 1e-14)
  expect_equal(got$max_log_q, max(expected_eta), tolerance = 1e-14)
  expect_equal(exp(mean(log(got$q))), 1, tolerance = 1e-10)
})

test_that("native BayesC-LV and BayesR-LV theta likelihoods match R oracle", {
  X <- cbind(
    binary = c(-0.5, 0.5, -0.5, 0.5),
    continuous = c(-1.1618950039, -0.387298335, 0.387298335, 1.1618950039)
  )
  theta <- c(0.4, -0.2)
  b <- c(0.1, 0, -0.3, 0.2)
  state <- c(1L, 0L, 1L, 1L)
  component <- c(1L, 0L, 2L, 3L)
  gamma <- c(0, 0.01, 0.1, 1)
  vb <- 0.05
  eta <- drop(X %*% theta)

  active_c <- state > 0L
  expected_c <- -0.5 * sum(
    eta[active_c] + b[active_c]^2 / (vb * exp(eta[active_c]))
  )
  active_r <- component > 0L
  expected_r <- -0.5 * sum(
    eta[active_r] + b[active_r]^2 /
      (vb * gamma[component[active_r] + 1L] * exp(eta[active_r]))
  )

  expect_equal(
    st_logvar_loglik_bayesc_internal(theta, X, b, state, vb),
    expected_c,
    tolerance = 1e-13
  )
  expect_equal(
    st_logvar_loglik_bayesr_internal(
      theta, X, b, component, vb, gamma
    ),
    expected_r,
    tolerance = 1e-13
  )
})

test_that("theta zero and fixed q reduce to ordinary marker kernels", {
  score <- 1.7
  wi <- 230
  ve <- 0.8
  vb <- 0.03
  pi_c <- 0.15
  pi_r <- c(0.8, 0.12, 0.06, 0.02)
  gamma <- c(0, 0.01, 0.1, 1)

  bayesc_weights <- function(q) {
    vbi <- vb * q
    denom <- ve + wi * vbi
    c(
      log(1 - pi_c),
      log(pi_c) + 0.5 * log(ve / denom) +
        0.5 * score^2 * vbi / (ve * denom)
    )
  }
  bayesr_weights <- function(q) {
    out <- log(pi_r)
    for (k in 2:length(gamma)) {
      vk <- vb * gamma[k] * q
      denom <- ve + wi * vk
      out[k] <- out[k] + 0.5 * log(ve / denom) +
        0.5 * score^2 * vk / (ve * denom)
    }
    out
  }

  expect_identical(bayesc_weights(exp(0)), bayesc_weights(1))
  expect_identical(bayesr_weights(exp(0)), bayesr_weights(1))

  q <- exp(0.37)
  direct_vb_c <- vb * q
  denom_c <- ve + wi * direct_vb_c
  direct_c <- c(
    log(1 - pi_c),
    log(pi_c) + 0.5 * log(ve / denom_c) +
      0.5 * score^2 * direct_vb_c / (ve * denom_c)
  )
  direct_r <- log(pi_r)
  for (k in 2:length(gamma)) {
    direct_vk <- vb * gamma[k] * q
    denom <- ve + wi * direct_vk
    direct_r[k] <- direct_r[k] + 0.5 * log(ve / denom) +
      0.5 * score^2 * direct_vk / (ve * denom)
  }
  expect_equal(bayesc_weights(q), direct_c, tolerance = 1e-15)
  expect_equal(bayesr_weights(q), direct_r, tolerance = 1e-15)
})

test_that("empty active sets draw theta directly from the frozen prior", {
  X <- cbind(x1 = seq(-1, 1, length.out = 20), x2 = rep(c(-0.5, 0.5), 10))
  fit <- st_logvar_ess_fixture_internal(
    theta = c(0, 0), annotation = X, effect = numeric(20),
    state = integer(20), marker_variance = 0.05,
    theta_prior_sd = 0.7, updates = 6000L, seed = 1729L,
    empty_active_set = TRUE
  )

  expect_equal(colMeans(fit$draws), c(0, 0), tolerance = 0.03)
  expect_equal(apply(fit$draws, 2, sd), c(0.7, 0.7), tolerance = 0.03)
  expect_equal(fit$theta_updates, 6000)
  expect_equal(fit$mean_likelihood_evaluations_per_update, 0)
  expect_equal(fit$max_likelihood_evaluations, 0)
  expect_equal(fit$mean_bracket_contractions, 0)
})

test_that("one-dimensional native ESS agrees with numerical integration", {
  X <- matrix(drop(scale(seq(-2, 2, length.out = 61))), ncol = 1)
  b <- 0.13 * sin(seq_len(nrow(X)) / 4) + 0.04
  state <- rep(1L, nrow(X))
  vb <- 0.05
  prior_sd <- 0.7

  grid <- seq(-2, 2, length.out = 20001)
  log_posterior <- vapply(grid, function(value) {
    eta <- drop(X) * value
    -0.5 * sum(eta + b^2 / (vb * exp(eta))) -
      0.5 * (value / prior_sd)^2
  }, numeric(1))
  weight <- exp(log_posterior - max(log_posterior))
  weight <- weight / sum(weight)
  expected_mean <- sum(grid * weight)
  expected_sd <- sqrt(sum((grid - expected_mean)^2 * weight))

  fit <- st_logvar_ess_fixture_internal(
    theta = 0, annotation = X, effect = b, state = state,
    marker_variance = vb, theta_prior_sd = prior_sd,
    updates = 10000L, seed = 8128L
  )
  kept <- drop(fit$draws[1001:10000, , drop = FALSE])

  expect_equal(mean(kept), expected_mean, tolerance = 0.04)
  expect_equal(sd(kept), expected_sd, tolerance = 0.04)
  expect_equal(fit$theta_updates, 10000)
  expect_gte(fit$mean_likelihood_evaluations_per_update, 1)
  expect_gte(fit$max_likelihood_evaluations, 1)
  expect_true(is.finite(fit$min_log_q))
  expect_true(is.finite(fit$max_log_q))
})

test_that("log-variance numerical guards fail clearly without clamping", {
  expect_error(
    st_logvar_eta_q_internal(matrix(1, 1, 1), 701),
    "outside the positive finite double range"
  )
  expect_error(
    st_logvar_eta_q_internal(matrix(Inf, 1, 1), 0),
    "must be finite"
  )
  expect_error(
    st_logvar_loglik_bayesc_internal(
      0, matrix(0, 1, 1), 0.1, 1L, 0
    ),
    "positive finite"
  )
  expect_error(
    st_logvar_loglik_bayesr_internal(
      0, matrix(0, 1, 1), 0.1, 1L, 0.1, c(0, 0)
    ),
    "positive gamma"
  )
  expect_error(
    st_logvar_ess_fixture_internal(
      0, matrix(0, 1, 1), 0.1, 1L, 0.1,
      theta_prior_sd = -0.7
    ),
    "positive finite prior SD"
  )
})
