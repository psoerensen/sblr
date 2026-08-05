reference_stick_state <- function(component, nstep) {
  indicator <- outer(as.integer(component), seq_len(nstep) - 1L, `>`)
  eligible <- lapply(seq_len(nstep), function(j) {
    if (j == 1L) seq_along(component) else which(indicator[, j - 1L])
  })
  list(indicator = indicator, eligible = eligible,
    continuation = vapply(seq_len(nstep), function(j)
      sum(indicator[eligible[[j]], j]), integer(1)),
    stopping = vapply(seq_len(nstep), function(j)
      sum(!indicator[eligible[[j]], j]), integer(1)))
}

reference_alpha_conditionals <- function(X, z, alpha, sigma_sq,
                                         intercept_mean, intercept_sd) {
  residual <- z - drop(X %*% alpha)
  out <- matrix(NA_real_, ncol(X), 2L)
  for (k in seq_len(ncol(X))) {
    x <- X[, k]
    diagonal <- sum(x^2)
    likelihood_rhs <- sum(x * residual) + diagonal * alpha[k]
    prior_precision <- if (k == 1L) intercept_sd^-2 else sigma_sq^-1
    prior_mean <- if (k == 1L) intercept_mean else 0
    variance <- 1 / (diagonal + prior_precision)
    out[k, ] <- c(variance *
      (likelihood_rhs + prior_precision * prior_mean), variance)
  }
  out
}

reference_variance_conditional <- function(alpha, a, b) {
  list(df = length(alpha) + a, numerator = sum(alpha^2) + b)
}

reference_prior_mapping <- function(a, b) {
  data.frame(nu0 = a, scale0 = b / a,
    prior_mean = if (a > 2) b / (a - 2) else Inf,
    prior_variance = if (a > 4)
      2 * b^2 / ((a - 2)^2 * (a - 4)) else Inf)
}

reference_component_probability <- function(stick) {
  stick <- as.matrix(stick)
  out <- matrix(0, nrow(stick), ncol(stick) + 1L)
  remaining <- rep(1, nrow(stick))
  for (j in seq_len(ncol(stick))) {
    out[, j] <- remaining * (1 - stick[, j])
    remaining <- remaining * stick[, j]
  }
  out[, ncol(out)] <- remaining
  out
}

test_that("alpha hierarchy stick eligibility follows continuation ordering", {
  component <- c(0L, 1L, 2L, 3L, 0L, 3L)
  state <- reference_stick_state(component, 3L)

  expect_identical(state$eligible, list(1:6, c(2L, 3L, 4L, 6L),
                                        c(3L, 4L, 6L)))
  expect_identical(state$continuation, c(4L, 3L, 2L))
  expect_identical(state$stopping, c(2L, 1L, 1L))
  native <- getFromNamespace(".st_bayesrc_annotation_update", "sblr")(
    cbind(intercept = 1, signal = seq_along(component)), component,
    matrix(0, 2, 3), rep(1, 3),
    sblr:::.sbayesrc_resolve_intercept_prior(rep(0.25, 4))$native,
    2, 2, 91L)
  expect_identical(native$eligible, lengths(state$eligible))
  expect_identical(native$continuation, state$continuation)
})

test_that("fixed-latent alpha conditional moments match independent algebra", {
  X <- cbind(intercept = 1,
             binary = c(0, 1, 0, 1, 1),
             continuous = c(-1, -0.25, 0.5, 1.25, 2))
  z <- c(-0.8, 0.4, -0.2, 1.1, 1.7)
  alpha <- c(-0.3, 0.6, -0.2)
  reference <- reference_alpha_conditionals(
    X, z, alpha, sigma_sq = 1.7, intercept_mean = -0.5,
    intercept_sd = 0.8)
  native <- getFromNamespace(
    ".st_bayesrc_alpha_conditional_moments", "sblr")(
      X, z, alpha, 1.7, -0.5, 0.8)

  expect_equal(unname(native), unname(reference), tolerance = 1e-14)
  expect_true(all(native[, 2L] > 0))
})

test_that("sigmaSqAlpha conditional and prior parameterization are exact", {
  teaching <- reference_prior_mapping(2, 2)
  production <- reference_prior_mapping(4, 4)
  expect_equal(teaching[c("nu0", "scale0")], data.frame(nu0 = 2, scale0 = 1))
  expect_equal(production[c("nu0", "scale0")], data.frame(nu0 = 4, scale0 = 1))
  expect_true(is.infinite(teaching$prior_mean))
  expect_equal(production$prior_mean, 2)
  expect_true(is.infinite(production$prior_variance))

  alpha <- c(-0.5, 0.75, 1.25)
  conditional <- reference_variance_conditional(alpha, 8, 3)
  draws <- getFromNamespace(".st_bayesrc_sigma_sq_alpha_draws", "sblr")(
    alpha, 8, 3, 50000L, 712L)
  probability <- c(0.1, 0.5, 0.9)
  expected <- conditional$numerator /
    stats::qchisq(1 - probability, conditional$df)
  expect_equal(unname(stats::quantile(draws, probability)), expected,
               tolerance = 0.025)
  expect_equal(mean(draws),
    conditional$numerator / (conditional$df - 2), tolerance = 0.02)
})

test_that("independent stick conversion matches package component probabilities", {
  stick <- rbind(c(0.1, 0.6, 0.25), c(0.8, 0.4, 0.9))
  reference <- reference_component_probability(stick)
  for (i in seq_len(nrow(stick))) {
    expect_equal(reference[i, ],
      sblr:::.sbayesrc_stick_to_component_prob(stick[i, ]),
      tolerance = 1e-15)
  }
  expect_equal(rowSums(reference), rep(1, nrow(reference)), tolerance = 1e-15)
  expect_true(all(reference > 0))
})
