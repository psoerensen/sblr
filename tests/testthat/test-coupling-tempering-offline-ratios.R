test_that("offline complete-state exchange ratio matches four target terms", {
  annotation <- cbind(1, c(0, 1, 0, 1))
  baseline <- qnorm(c(0.7, 0.4))
  alpha_a <- matrix(c(0.2, -0.3, 0.4, 0.5), 2, 2)
  alpha_b <- matrix(c(-0.1, 0.7, -0.4, 0.2), 2, 2)
  component_a <- c(0L, 1L, 2L, 0L)
  component_b <- c(1L, 0L, 2L, 1L)

  ell <- .offline_log_allocation_prior
  observed <- .offline_complete_exchange_ratio(
    annotation, baseline, 0, 0.5,
    component_a, alpha_a, component_b, alpha_b
  )
  expected <- ell(annotation, alpha_b, baseline, 0, component_b) +
    ell(annotation, alpha_a, baseline, 0.5, component_a) -
    ell(annotation, alpha_a, baseline, 0, component_a) -
    ell(annotation, alpha_b, baseline, 0.5, component_b)

  expect_equal(observed, expected, tolerance = 1e-13)
  expect_equal(
    observed,
    -.offline_complete_exchange_ratio(
      annotation, baseline, 0, 0.5,
      component_b, alpha_b, component_a, alpha_a
    ),
    tolerance = 1e-13
  )
})

test_that("joint alpha-sigma exchange cancels hierarchy-prior factors", {
  annotation <- cbind(1, c(-1, 0, 1, 2))
  baseline <- qnorm(c(0.8, 0.3))
  alpha_a <- matrix(c(0.3, -0.2, -0.5, 0.4), 2, 2)
  alpha_b <- matrix(c(-0.4, 0.8, 0.6, -0.1), 2, 2)
  sigma_a <- c(0.4, 1.7)
  sigma_b <- c(2.1, 0.6)
  component_a <- c(0L, 1L, 2L, 0L)
  component_b <- c(1L, 0L, 2L, 1L)
  ell <- .offline_log_allocation_prior

  observed <- .offline_alpha_sigma_exchange_ratio(
    annotation, baseline, 0.5, 1,
    component_a, alpha_a, sigma_a,
    component_b, alpha_b, sigma_b,
    prior_df = 2, prior_scale = 2
  )
  expected <- ell(annotation, alpha_b, baseline, 0.5, component_a) +
    ell(annotation, alpha_a, baseline, 1, component_b) -
    ell(annotation, alpha_a, baseline, 0.5, component_a) -
    ell(annotation, alpha_b, baseline, 1, component_b)

  expect_equal(observed, expected, tolerance = 1e-13)
})

test_that("alpha-only exchange includes destination variance densities", {
  annotation <- cbind(1, c(0, 1, 0, 1))
  baseline <- qnorm(c(0.75, 0.35))
  alpha_a <- matrix(c(0.1, -0.5, 0.2, 0.7), 2, 2)
  alpha_b <- matrix(c(-0.3, 1.1, -0.6, 0.4), 2, 2)
  sigma_a <- c(0.25, 2)
  sigma_b <- c(1.5, 0.5)
  component_a <- c(0L, 1L, 2L, 0L)
  component_b <- c(1L, 0L, 2L, 1L)
  ell <- .offline_log_allocation_prior
  h <- .offline_log_nonintercept_prior

  observed <- .offline_alpha_only_exchange_ratio(
    annotation, baseline, 0, 1,
    component_a, alpha_a, sigma_a,
    component_b, alpha_b, sigma_b
  )
  expected <- ell(annotation, alpha_b, baseline, 0, component_a) +
    ell(annotation, alpha_a, baseline, 1, component_b) -
    ell(annotation, alpha_a, baseline, 0, component_a) -
    ell(annotation, alpha_b, baseline, 1, component_b) +
    h(alpha_b, sigma_a) + h(alpha_a, sigma_b) -
    h(alpha_a, sigma_a) - h(alpha_b, sigma_b)

  expect_equal(observed, expected, tolerance = 1e-13)
})

test_that("offline tempered probabilities are finite and normalized", {
  annotation <- cbind(1, seq(-2, 2, length.out = 7))
  alpha <- matrix(c(-1, 0.8, 0.4, -0.6), 2, 2)
  baseline <- qnorm(c(0.9, 0.25))

  for (lambda in c(0, 0.5, 1)) {
    probability <- .offline_tempered_component_probability(
      annotation, alpha, baseline, lambda
    )
    expect_true(all(is.finite(probability)))
    expect_true(all(probability > 0))
    expect_equal(rowSums(probability), rep(1, nrow(annotation)),
                 tolerance = 1e-14)
  }
})
