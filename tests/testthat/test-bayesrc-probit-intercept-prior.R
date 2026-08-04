test_that("component and continuation probabilities round trip in package order", {
  component <- c(0.91, 0.04, 0.03, 0.02)
  stick <- sblr:::.sbayesrc_component_to_stick_prob(component)

  expect_equal(stick, c(0.09, 5 / 9, 0.4), tolerance = 1e-14)
  expect_equal(
    sblr:::.sbayesrc_stick_to_component_prob(stick),
    component,
    tolerance = 1e-14
  )
  expect_equal(
    sblr:::.sbayesrc_component_to_stick_prob(c(1, 0, 0, 0)),
    c(.Machine$double.eps, 0.5, 0.5)
  )
})

test_that("proper intercept prior resolves separately from initialization", {
  component <- c(0.999, rep(0.001 / 3, 3))
  default <- sblr:::.sbayesrc_resolve_intercept_prior(component)
  custom <- sblr:::.sbayesrc_resolve_intercept_prior(
    component,
    list(distribution = "normal", mean = c(-2, 0, 1), sd = c(0.5, 1, 2))
  )

  expect_identical(default$distribution, "normal")
  expect_false(default$legacy_flat)
  expect_equal(unname(default$mean), unname(qnorm(default$stick_probability)))
  expect_equal(unname(default$sd), rep(1, 3))
  expect_equal(unname(custom$mean), c(-2, 0, 1))
  expect_equal(unname(custom$precision), c(4, 1, 0.25))
  expect_equal(unname(custom$native[2, ]), c(-2, 0, 1))

  A <- matrix(1, 4, 1)
  initialized <- make_sbayesrc_alpha_init(
    A, alpha_init = matrix(7, 1, 3), pi_init = 0.001)
  expect_equal(unname(initialized$alpha_init), matrix(7, 1, 3))
  expect_equal(unname(default$mean), unname(qnorm(default$stick_probability)))
})

test_that("intercept prior validates dimensions and legacy migration", {
  component <- c(0.9, 0.05, 0.03, 0.02)
  expect_error(
    sblr:::.sbayesrc_resolve_intercept_prior(component, list(sd = c(1, 2))),
    "scalar or stick-specific"
  )
  expect_error(
    sblr:::.sbayesrc_resolve_intercept_prior(component, list(mean = Inf)),
    "finite"
  )
  expect_error(
    sblr:::.sbayesrc_resolve_intercept_prior(component, intercept_flat = FALSE),
    "historical"
  )
  expect_warning(
    legacy <- sblr:::.sbayesrc_resolve_intercept_prior(
      component, intercept_flat = TRUE),
    "legacy improper"
  )
  expect_true(legacy$legacy_flat)
})

test_that("empty sticks are prior-only Gibbs updates", {
  update <- getFromNamespace(".st_bayesrc_annotation_update", "sblr")
  prior <- sblr:::.sbayesrc_resolve_intercept_prior(
    c(0.8, 0.15, 0.05),
    list(mean = c(-0.5, 0.75), sd = c(0.5, 1.25))
  )
  # Stick two is empty because every component is zero. Stick one is populated.
  out <- update(
    cbind(intercept = rep(1, 8), signal = seq(-1, 1, length.out = 8)),
    rep(0, 8), matrix(0, 2, 2), c(1, 1), prior$native, 2, 2, 19L
  )

  expect_equal(out$eligible, c(8L, 0L))
  expect_equal(out$continuation, c(0L, 0L))
  expect_equal(out$prior_only, c(0L, 1L))
  expect_true(all(is.finite(out$alpha)))
  expect_true(all(is.finite(out$sigmaSqAlpha)))
  expect_true(all(out$sigmaSqAlpha > 0))
})

test_that("proper prior remains finite under both directions of separation", {
  update <- getFromNamespace(".st_bayesrc_annotation_update", "sblr")
  prior <- sblr:::.sbayesrc_resolve_intercept_prior(c(0.9, 0.1))
  run <- function(component) {
    alpha <- matrix(0, 1, 1)
    sigma <- 1
    for (seed in seq_len(250)) {
      out <- update(matrix(1, length(component), 1), component, alpha, sigma,
                    prior$native, 2, 2, seed)
      alpha <- out$alpha
      sigma <- out$sigmaSqAlpha
    }
    c(alpha, sigma)
  }

  expect_true(all(is.finite(run(rep(0, 20)))))
  expect_true(all(is.finite(run(rep(1, 20)))))
})

test_that("intercept-only Gibbs chains match numerical posterior references", {
  update <- getFromNamespace(".st_bayesrc_annotation_update", "sblr")
  prior <- sblr:::.sbayesrc_resolve_intercept_prior(
    c(0.75, 0.25), list(mean = -0.4, sd = 0.8))
  grid <- seq(-5, 5, length.out = 40001)
  dx <- grid[2] - grid[1]
  reference <- function(success, total) {
    log_weight <- success * pnorm(grid, log.p = TRUE) +
      (total - success) * pnorm(grid, lower.tail = FALSE, log.p = TRUE) +
      dnorm(grid, -0.4, 0.8, log = TRUE)
    weight <- exp(log_weight - max(log_weight))
    weight <- weight / (sum(weight) * dx)
    mean <- sum(grid * weight) * dx
    c(mean = mean, variance = sum((grid - mean)^2 * weight) * dx)
  }
  sample_case <- function(success, total, seed_offset) {
    component <- c(rep(1, success), rep(0, total - success))
    alpha <- matrix(0, 1, 1)
    sigma <- 1
    draws <- numeric(1200)
    for (iteration in seq_along(draws)) {
      out <- update(matrix(1, total, 1), component, alpha, sigma,
                    prior$native, 2, 2, seed_offset + iteration)
      alpha <- out$alpha
      sigma <- out$sigmaSqAlpha
      draws[iteration] <- alpha[1, 1]
    }
    draws[-seq_len(200)]
  }

  for (case in list(c(20, 20), c(0, 20), c(19, 20), c(10, 20), c(18, 20))) {
    ref <- reference(case[1], case[2])
    draws <- sample_case(case[1], case[2], 10000 * case[1] + case[2])
    expect_equal(mean(draws), unname(ref["mean"]), tolerance = 0.15)
    expect_equal(var(draws), unname(ref["variance"]), tolerance = 0.15)
  }
})

test_that("legacy flat prior fails for unidentified sticks", {
  update <- getFromNamespace(".st_bayesrc_annotation_update", "sblr")
  expect_warning(
    legacy <- sblr:::.sbayesrc_resolve_intercept_prior(
      c(0.9, 0.1), intercept_flat = TRUE),
    "legacy improper"
  )
  expect_error(
    update(matrix(1, 10, 1), rep(1, 10), matrix(0, 1, 1), 1,
           legacy$native, 2, 2, 1L),
    "complete separation"
  )
})
