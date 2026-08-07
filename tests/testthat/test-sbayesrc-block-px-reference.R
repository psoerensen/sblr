test_that("PX scale ratio matches direct integrated latent density", {
  design <- cbind(1, c(-1, -0.2, 0.4, 1.1, 0.7))
  latent <- c(-1.2, -0.3, 0.2, 1.4, 0.8)
  prior_mean <- c(qnorm(0.2), 0)
  prior_precision <- c(1, 0.7)
  for (log_scale in c(-0.4, -0.1, 0, 0.2, 0.55)) {
    scale <- exp(log_scale)
    direct <- length(latent) * log_scale +
      sblr:::.sbayesrc_px_log_latent_marginal(
        scale * latent, design, prior_mean, prior_precision) -
      sblr:::.sbayesrc_px_log_latent_marginal(
        latent, design, prior_mean, prior_precision)
    expect_equal(
      sblr:::.sbayesrc_px_log_scale_ratio(
        log_scale, latent, design, prior_mean, prior_precision),
      direct, tolerance = 1e-12)
  }
})

test_that("PX scale move satisfies forward-reverse detailed balance", {
  design <- cbind(1, c(0, 1, 0, 1, 1, 0))
  latent <- c(-0.8, 0.5, -1.1, 1.2, 0.3, -0.2)
  prior_mean <- c(-0.9, 0)
  prior_precision <- c(1, 1.4)
  log_scale <- 0.37
  scale <- exp(log_scale)
  forward <- sblr:::.sbayesrc_px_log_scale_ratio(
    log_scale, latent, design, prior_mean, prior_precision)
  reverse <- sblr:::.sbayesrc_px_log_scale_ratio(
    -log_scale, scale * latent, design, prior_mean, prior_precision)
  expect_equal(forward, -reverse, tolerance = 1e-12)

  log_old <- sblr:::.sbayesrc_px_log_latent_marginal(
    latent, design, prior_mean, prior_precision)
  log_new <- sblr:::.sbayesrc_px_log_latent_marginal(
    scale * latent, design, prior_mean, prior_precision)
  forward_flow <- log_old + min(0, forward)
  reverse_flow <- log_new + length(latent) * log_scale + min(0, reverse)
  expect_equal(forward_flow, reverse_flow, tolerance = 1e-12)
})

test_that("blocked PX alpha conditional matches independent Gaussian algebra", {
  design <- cbind(1, c(-1, 0, 0.5, 1))
  latent <- c(-0.7, -0.1, 0.8, 1.3)
  prior_mean <- c(-0.4, 0)
  prior_precision <- c(0.8, 1.7)
  observed <- sblr:::.sbayesrc_px_alpha_conditional(
    latent, design, prior_mean, prior_precision)
  precision <- crossprod(design) + diag(prior_precision)
  expected_covariance <- solve(precision)
  expected_mean <- expected_covariance %*%
    (crossprod(design, latent) + prior_precision * prior_mean)
  expect_equal(observed$mean, drop(expected_mean), tolerance = 1e-13)
  expect_equal(observed$covariance, expected_covariance, tolerance = 1e-13)
})

test_that("PX deterministic references are RNG neutral and reject bad inputs", {
  set.seed(90117)
  seed <- .Random.seed
  design <- cbind(1, c(0, 1, 0))
  expect_true(is.finite(sblr:::.sbayesrc_px_log_scale_ratio(
    0.1, c(-1, 0.4, -0.2), design, c(-0.5, 0), c(1, 1))))
  expect_identical(.Random.seed, seed)
  expect_error(sblr:::.sbayesrc_px_factor(
    c(1, 2), matrix(1, 3, 1), 0, 1), "Invalid")
  expect_equal(sblr:::.sbayesrc_px_log_scale_ratio(
    1000, c(-1, 1, -0.5), design, c(-0.5, 0), c(1, 1)), -Inf)
})

test_that("native PX annotation update is finite, deterministic, and diagnosed", {
  update <- getFromNamespace(".st_bayesrc_annotation_update", "sblr")
  prior <- sblr:::.sbayesrc_resolve_intercept_prior(c(0.8, 0.15, 0.05))
  native <- rbind(
    prior$native,
    update_sigmaSqAlpha = c(1, 1),
    allocation_updates_per_cycle = c(1, 1),
    annotation_updates_per_cycle = c(1, 1),
    coupling_tempering = c(0, 0),
    coupling_swap_every = c(0, 0),
    px_sandwich = c(1, 1),
    px_log_scale_sd = c(0.45, 0.45))
  annotation <- cbind(intercept = 1, signal = seq(-1, 1, length.out = 18))
  component <- rep(0:2, each = 6L)
  args <- list(annotation, component, matrix(0, 2, 2), c(1, 1),
               native, 2, 2, 81723L)
  first <- do.call(update, args)
  second <- do.call(update, args)
  expect_identical(first, second)
  expect_true(all(is.finite(first$alpha)))
  expect_true(all(is.finite(first$sigmaSqAlpha)))
  expect_equal(first$px_attempted, c(1L, 1L))
  expect_true(all(first$px_accepted %in% 0:1))
  expect_true(all(first$px_abs_log_scale >= 0))
  expect_true(all(first$px_alpha_jump > 0))
})
