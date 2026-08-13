.learned_logistic_softplus <- function(x) {
  out <- numeric(length(x))
  positive <- x > 0
  out[positive] <- x[positive] + log1p(exp(-x[positive]))
  out[!positive] <- log1p(exp(x[!positive]))
  out
}

.learned_logistic_trapezoid <- function(x, y) {
  sum(diff(x) * (head(y, -1L) + tail(y, -1L)) / 2)
}

.learned_logistic_reference <- function(state, offset, prior_a, prior_b,
                                        grid_size = 100001L) {
  z <- seq(-25, 25, length.out = grid_size)
  log_density <-
    -prior_a * .learned_logistic_softplus(-z) -
    prior_b * .learned_logistic_softplus(z)
  for (j in seq_along(state)) {
    marker_logit <- z + offset[j]
    log_density <- log_density + if (state[j] == 1L) {
      -.learned_logistic_softplus(-marker_logit)
    } else {
      -.learned_logistic_softplus(marker_logit)
    }
  }
  density <- exp(log_density - max(log_density))
  density <- density / .learned_logistic_trapezoid(z, density)
  probability <- plogis(z)
  mean <- .learned_logistic_trapezoid(z, probability * density)
  variance <- .learned_logistic_trapezoid(
    z, (probability - mean)^2 * density
  )
  interval_mass <- diff(z) * (head(density, -1L) + tail(density, -1L)) / 2
  cdf <- c(0, cumsum(interval_mass))
  quantiles <- vapply(c(0.25, 0.5, 0.75), function(level) {
    probability[which(cdf >= level)[1L]]
  }, numeric(1))
  list(
    mean = mean,
    variance = variance,
    quantiles = quantiles,
    z = z,
    density = density
  )
}

.learned_logistic_native_draws <- function(state, offset, prior_a, prior_b,
                                            draws, burnin, seed,
                                            pi_init = 0.3) {
  sblr:::stblr_learned_logistic_pi_draws_internal(
    as.integer(state), as.numeric(offset), pi_init,
    prior_a, prior_b, as.integer(draws), as.integer(burnin), as.integer(seed)
  )
}

test_that("R-side learned calibration is unclipped while fixed calibration is unchanged", {
  A <- matrix(c(-1, 0, 1), ncol = 1L)
  learned <- sblr:::.stblr_make_prior_from_annotations(
    A = A, nt = 1L, pi_base = 0.7, beta_pi = 1,
    clip_probability = FALSE
  )
  fixed <- sblr:::.stblr_make_prior_from_annotations(
    A = A, nt = 1L, pi_base = 0.7, beta_pi = 1
  )
  expect_gt(max(learned$pi_marker[[1L]]), 0.5)
  expect_equal(
    learned$pi_marker[[1L]],
    plogis(qlogis(0.7) + c(-1, 0, 1)),
    tolerance = 1e-15
  )
  expect_lte(max(fixed$pi_marker[[1L]]), 0.5)
})

test_that("zero learned-logistic offsets reduce exactly to the Beta conditional", {
  state <- c(1L, 0L, 1L, 0L, 0L, 1L)
  offset <- numeric(length(state))
  prior_a <- 0.7
  prior_b <- 1.3
  shape_a <- prior_a + sum(state)
  shape_b <- prior_b + length(state) - sum(state)
  beta_mean <- shape_a / (shape_a + shape_b)
  beta_variance <- shape_a * shape_b /
    ((shape_a + shape_b)^2 * (shape_a + shape_b + 1))

  reference <- .learned_logistic_reference(state, offset, prior_a, prior_b)
  expect_equal(reference$mean, beta_mean, tolerance = 2e-8)
  expect_equal(reference$variance, beta_variance, tolerance = 2e-8)

  p <- plogis(reference$z)
  beta_z_density <- dbeta(p, shape_a, shape_b) * p * (1 - p)
  beta_z_density <- beta_z_density /
    .learned_logistic_trapezoid(reference$z, beta_z_density)
  expect_lt(max(abs(reference$density - beta_z_density)), 2e-8)

  native <- .learned_logistic_native_draws(
    state, offset, prior_a, prior_b, draws = 50000L, burnin = 0L, seed = 7301L
  )
  native_repeat <- .learned_logistic_native_draws(
    state, offset, prior_a, prior_b, draws = 50000L, burnin = 0L, seed = 7301L
  )
  expect_true(native$conjugate_reduction)
  expect_identical(native$draws, native_repeat$draws)
  expect_lt(abs(mean(native$draws) - beta_mean), 0.004)
  expect_lt(abs(stats::var(native$draws) - beta_variance), 0.002)
})

test_that("centered nonzero offsets use the nonconjugate conditional", {
  state <- c(1L, 0L, 0L, 0L, 0L, 1L)
  offset <- c(-4, -2, -1, 1, 2, 4)
  prior_a <- 0.7
  prior_b <- 1.3
  expect_equal(sum(offset), 0)

  reference <- .learned_logistic_reference(state, offset, prior_a, prior_b)
  old_shape_a <- prior_a + sum(state)
  old_shape_b <- prior_b + length(state) - sum(state)
  old_beta_mean <- old_shape_a / (old_shape_a + old_shape_b)
  expect_gt(abs(reference$mean - old_beta_mean), 0.05)

  native <- .learned_logistic_native_draws(
    state, offset, prior_a, prior_b,
    draws = 40000L, burnin = 2000L, seed = 7302L
  )
  expect_false(native$conjugate_reduction)
  expect_lt(abs(mean(native$draws) - reference$mean), 0.012)
  expect_lt(abs(stats::var(native$draws) - reference$variance), 0.006)
  expect_lt(max(abs(
    unname(stats::quantile(native$draws, c(0.25, 0.5, 0.75))) -
      reference$quantiles
  )), 0.025)
  expect_true(all(is.finite(native$marker_probability)))
  expect_true(all(native$marker_probability > 0 & native$marker_probability < 1))
})

test_that("unequal uncentered offsets agree with numerical quadrature", {
  state <- c(1L, 0L, 1L, 0L, 0L)
  offset <- c(-3, -1, 0.5, 2, 4)
  prior_a <- 1.4
  prior_b <- 0.8
  expect_false(isTRUE(all.equal(sum(offset), 0)))

  reference <- .learned_logistic_reference(state, offset, prior_a, prior_b)
  old_beta_mean <- (prior_a + sum(state)) /
    (prior_a + prior_b + length(state))
  expect_gt(abs(reference$mean - old_beta_mean), 0.05)

  native <- .learned_logistic_native_draws(
    state, offset, prior_a, prior_b,
    draws = 40000L, burnin = 2000L, seed = 7304L
  )
  expect_false(native$conjugate_reduction)
  expect_lt(abs(mean(native$draws) - reference$mean), 0.012)
  expect_lt(abs(stats::var(native$draws) - reference$variance), 0.006)
  expect_lt(max(abs(
    unname(stats::quantile(native$draws, c(0.25, 0.5, 0.75))) -
      reference$quantiles
  )), 0.025)
})

test_that("the removed probability ceiling defined a different target", {
  state <- c(1L, 0L, 1L, 0L, 0L, 1L)
  prior_a <- 0.7
  prior_b <- 1.3
  p <- seq(1e-7, 1 - 1e-7, length.out = 200001L)
  capped <- pmin(p, 0.5)
  log_density <-
    (prior_a - 1) * log(p) + (prior_b - 1) * log1p(-p) +
    sum(state) * log(capped) +
    (length(state) - sum(state)) * log1p(-capped)
  density <- exp(log_density - max(log_density))
  density <- density / .learned_logistic_trapezoid(p, density)
  capped_mean <- .learned_logistic_trapezoid(p, p * density)
  beta_mean <- (prior_a + sum(state)) /
    (prior_a + prior_b + length(state))
  expect_gt(abs(capped_mean - beta_mean), 0.08)
})

test_that("eta and marker probability calculations use stable unclipped logits", {
  A <- matrix(c(1, 0, 0, 1, 1, 1, -1, 2), nrow = 4L, byrow = TRUE)
  eta <- c(0.7, -0.4)
  state <- c(1L, 0L, 1L, 0L)
  base_pi <- 0.3
  sigma_eta <- 1.2
  offset <- as.numeric(A %*% eta)
  offset <- offset - mean(offset)
  marker_logit <- qlogis(base_pi) + offset
  expected <- sum(ifelse(
    state == 1L,
    -.learned_logistic_softplus(-marker_logit),
    -.learned_logistic_softplus(marker_logit)
  )) - sum(eta^2) / (2 * sigma_eta^2)
  observed <- sblr:::stblr_learned_logistic_eta_logpost_internal(
    A, eta, state, base_pi, sigma_eta
  )
  expect_equal(observed, expected, tolerance = 1e-12)

  predictor <- c(-1000, -50, 0, 50, 1000)
  terms <- sblr:::stblr_learned_logistic_probability_terms_internal(predictor)
  expect_true(all(is.finite(terms)))
  expect_true(all(terms[, "probability"] > 0 & terms[, "probability"] < 1))
  expect_equal(
    terms[, "log_probability"],
    -.learned_logistic_softplus(-predictor),
    tolerance = 1e-13
  )
  expect_equal(
    terms[, "log_complement"],
    -.learned_logistic_softplus(predictor),
    tolerance = 1e-13
  )
})

test_that("learned-logistic pi update validates boundaries and inputs", {
  call <- function(state, offset, a = 1, b = 1, pi = 0.3) {
    .learned_logistic_native_draws(
      state, offset, a, b, draws = 2000L, burnin = 100L,
      seed = 7303L, pi_init = pi
    )
  }
  for (state in list(rep(0L, 5L), rep(1L, 5L))) {
    for (shape in c(0.4, 2.5)) {
      result <- call(state, c(-1000, -50, 0, 50, 1000), shape, shape)
      expect_true(all(is.finite(result$draws)))
      expect_true(all(result$draws > 0 & result$draws < 1))
      expect_true(all(result$marker_probability > 0 &
                      result$marker_probability < 1))
    }
  }
  expect_error(call(c(0L, 2L), c(0, 1)), "states must be zero or one")
  expect_error(call(c(0L, 1L), c(0, Inf)), "offsets must be finite")
  expect_error(call(c(0L, 1L), c(0, 1), a = 0), "shapes must be finite and positive")
  expect_error(call(c(0L, 1L), c(0, 1), b = Inf), "shapes must be finite and positive")
  expect_error(call(c(0L, 1L), c(0, 1), pi = 0), "inside \\(0, 1\\)")
  expect_error(
    sblr:::stblr_learned_logistic_eta_logpost_internal(
      matrix(1, 2, 1), Inf, c(0L, 1L), 0.3, 1
    ),
    "predictor became non-finite"
  )
})
