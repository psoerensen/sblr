test_that("BayesRC truncated-normal latent draws remain valid in extreme tails", {
  sampler <- getFromNamespace(".st_bayesrc_truncated_normal_draws", "sblr")
  location <- c(-40, -8, 0, 8, 40)
  value <- sampler(location, draws = 2000L, seed = 431L)

  expect_true(all(is.finite(value)))
  expect_true(all(value[, seq_along(location), drop = FALSE] > 0))
  expect_true(all(value[, length(location) + seq_along(location), drop = FALSE] < 0))
  expect_identical(value, sampler(location, draws = 2000L, seed = 431L))
})

test_that("BayesRC tail sampler has the intended truncated-normal mean", {
  sampler <- getFromNamespace(".st_bayesrc_truncated_normal_draws", "sblr")
  value <- sampler(c(-8, 8), draws = 20000L, seed = 917L)
  inverse_mills <- exp(stats::dnorm(8, log = TRUE) -
    stats::pnorm(8, lower.tail = FALSE, log.p = TRUE))
  expected_positive <- -8 + inverse_mills

  expect_equal(mean(value[, 1L]), expected_positive, tolerance = 0.006)
  expect_equal(mean(value[, 4L]), -expected_positive, tolerance = 0.006)
})
