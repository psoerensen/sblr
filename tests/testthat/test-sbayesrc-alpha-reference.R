alpha_reference_fixture <- function() {
  observations <- 180L
  annotation <- cbind(
    enriched_binary = rep(c(0, 1), length.out = observations),
    continuous_signal = as.numeric(scale(seq(-2, 2, length.out = observations))),
    null_annotation = as.numeric(scale(cos(seq(0, 5 * pi, length.out = observations))))
  )
  design <- cbind(Intercept = 1, annotation)
  latent <- drop(design %*% c(-0.45, 0.7, -0.25, 0.1)) +
    0.8 * sin(seq(0, 8 * pi, length.out = observations))
  list(annotation = annotation, latent = latent, tau2 = 1.3)
}

test_that("fixed-z scalar alpha Gibbs matches the exact Gaussian target", {
  fixture <- alpha_reference_fixture()
  exact <- .sbayesrc_alpha_exact_posterior(
    fixture$latent, fixture$annotation, fixture$tau2
  )
  set.seed(20260807)
  draws <- .sbayesrc_alpha_scalar_chain(
    fixture$latent, fixture$annotation, fixture$tau2,
    iterations = 8000L, burn = 1000L
  )

  expect_equal(colMeans(draws), exact$mean, tolerance = 0.02)
  expect_equal(stats::cov(draws), exact$covariance, tolerance = 0.006)
})

test_that("blocked alpha reference matches the same exact Gaussian target", {
  fixture <- alpha_reference_fixture()
  exact <- .sbayesrc_alpha_exact_posterior(
    fixture$latent, fixture$annotation, fixture$tau2
  )
  set.seed(20260808)
  draws <- .sbayesrc_alpha_blocked_draws(
    fixture$latent, fixture$annotation, fixture$tau2, draws = 6000L
  )

  expect_equal(colMeans(draws), exact$mean, tolerance = 0.02)
  expect_equal(stats::cov(draws), exact$covariance, tolerance = 0.006)
})

test_that("production scalar conditional moments match the independent formula", {
  fixture <- alpha_reference_fixture()
  alpha <- c(-0.2, 0.5, -0.1, 0.25)
  intercept_mean <- -0.35
  intercept_sd <- 0.9
  reference <- .sbayesrc_alpha_scalar_moments(
    fixture$latent, fixture$annotation, alpha, fixture$tau2,
    intercept_mean = intercept_mean,
    intercept_precision = intercept_sd^-2
  )
  production <- getFromNamespace(
    ".st_bayesrc_alpha_conditional_moments", "sblr"
  )(
    cbind(Intercept = 1, fixture$annotation), fixture$latent, alpha,
    fixture$tau2, intercept_mean, intercept_sd
  )

  expect_equal(unname(production), unname(reference), tolerance = 1e-14)

  flat_reference <- .sbayesrc_alpha_scalar_moments(
    fixture$latent, fixture$annotation, alpha, fixture$tau2
  )
  near_flat_production <- getFromNamespace(
    ".st_bayesrc_alpha_conditional_moments", "sblr"
  )(
    cbind(Intercept = 1, fixture$annotation), fixture$latent, alpha,
    fixture$tau2, 0, 1e8
  )
  expect_equal(
    unname(near_flat_production), unname(flat_reference), tolerance = 1e-13
  )
})

test_that("later sticks use matching eligible outcomes and annotation rows", {
  component <- c(0L, 0L, 1L, 1L, 2L, 2L, 3L, 3L)
  annotation <- cbind(
    Intercept = 1,
    enriched_binary = c(0, 1, 0, 1, 0, 1, 0, 1),
    continuous_signal = seq(-1, 1, length.out = length(component)),
    null_annotation = rep(c(-1, 1), length.out = length(component))
  )
  prior <- sblr:::.sbayesrc_resolve_intercept_prior(rep(0.25, 4))
  update <- getFromNamespace(".st_bayesrc_annotation_update", "sblr")
  production <- update(
    annotation, component, matrix(0, 4, 3), rep(1, 3), prior$native,
    2, 2, 718L
  )

  expect_identical(production$eligible, c(8L, 6L, 4L))
  expect_identical(production$continuation, c(6L, 4L, 2L))
  for (stick in seq_len(3L)) {
    eligible <- if (stick == 1L) {
      seq_along(component)
    } else {
      which(component > stick - 2L)
    }
    expect_equal(sum(annotation[eligible, "Intercept"]^2), length(eligible))
  }
})

test_that("probit sticks map to finite normalized marker probabilities", {
  annotation <- cbind(
    Intercept = 1,
    enriched_binary = c(0, 1, 0, 1, 1),
    continuous_signal = c(-1.2, -0.4, 0, 0.7, 1.3),
    null_annotation = c(0.5, -0.8, 0.1, 0.9, -0.2)
  )
  alpha <- matrix(c(
    -1.1, 0.8, 0.25, 0,
    -0.3, 0.35, -0.15, 0.05,
    0.1, -0.2, 0.1, 0
  ), nrow = 4L, ncol = 3L)
  stick <- stats::pnorm(annotation %*% alpha)
  reference <- .sbayesrc_alpha_reference_component_probability(stick)
  production <- sbayesrc_marker_pi(annotation, alpha)

  expect_true(all(is.finite(stick)))
  expect_true(all(stick >= 0 & stick <= 1))
  expect_true(all(is.finite(production)))
  expect_true(all(production >= 0 & production <= 1))
  expect_equal(rowSums(production), rep(1, nrow(annotation)), tolerance = 1e-15)
  expect_equal(unname(production), unname(reference), tolerance = 1e-15)
})

