test_that("SBayesRC-S flat-intercept model enumeration is normalized", {
  fixture <- .sbs_fixture()
  exact <- .sbs_exact_posterior(
    fixture$z, fixture$annotation, fixture$eligible,
    fixture$pi_a, fixture$tau2
  )

  expect_length(exact$model_probability, 2^ncol(fixture$annotation))
  expect_true(all(is.finite(exact$model_probability)))
  expect_true(all(exact$model_probability >= 0))
  expect_equal(sum(exact$model_probability), 1, tolerance = 1e-14)
  expect_true(all(exact$annotation_pip >= 0 & exact$annotation_pip <= 1))
  expect_true(all(is.finite(exact$q_mean)))
  expect_true(all(exact$q_mean >= 0 & exact$q_mean <= 1))
  expect_true(all(is.finite(exact$pi_mean)))
  expect_true(all(exact$pi_mean >= 0 & exact$pi_mean <= 1))
  expect_equal(rowSums(exact$pi_mean), rep(1, nrow(fixture$annotation)),
               tolerance = 1e-14)
})

test_that("collapsed annotation Bayes factor matches direct Gaussian integration", {
  fixture <- .sbs_fixture()
  rows <- fixture$eligible[[2L]]
  x <- fixture$annotation[rows, "continuous_signal"]
  residual <- fixture$z[[2L]] - fixture$true_intercept[2L] -
    fixture$annotation[rows, "enriched_binary"] * 0.17

  formula <- .sbs_log_bf(x, residual, fixture$tau2[2L])
  direct <- .sbs_log_bf_direct(x, residual, fixture$tau2[2L])
  expect_equal(formula, direct, tolerance = 1e-10)
})

test_that("global model odds agree with integrated model weights", {
  fixture <- .sbs_fixture()
  exact <- .sbs_exact_posterior(
    fixture$z, fixture$annotation, fixture$eligible,
    fixture$pi_a, fixture$tau2
  )
  excluded <- match("000", rownames(exact$states))
  included <- match("100", rownames(exact$states))
  integrated_odds <- stats::qlogis(fixture$pi_a) +
    sum(vapply(seq_along(fixture$z), function(stick) {
      exact$stick_posterior[[included]][[stick]]$log_marginal -
        exact$stick_posterior[[excluded]][[stick]]$log_marginal
    }, numeric(1L)))
  expect_equal(
    exact$log_weight[included] - exact$log_weight[excluded],
    integrated_odds, tolerance = 1e-13
  )
  expect_equal(
    unname(log(
      exact$model_probability[included] / exact$model_probability[excluded]
    )),
    integrated_odds, tolerance = 1e-13
  )
})

test_that("collapsed SBayesRC-S MCMC reproduces the exact fixed-z posterior", {
  fixture <- .sbs_fixture()
  exact <- .sbs_exact_posterior(
    fixture$z, fixture$annotation, fixture$eligible,
    fixture$pi_a, fixture$tau2
  )
  initial <- list(c(0L, 0L, 0L), c(1L, 1L, 1L),
                  c(1L, 0L, 1L), c(0L, 1L, 0L))
  chains <- lapply(seq_along(initial), function(chain) {
    set.seed(20260900L + chain)
    .sbs_mcmc_chain(
      fixture$z, fixture$annotation, fixture$eligible,
      fixture$pi_a, fixture$tau2,
      iterations = 6000L, burn = 1000L, initial_delta = initial[[chain]]
    )
  })
  pooled <- .sbs_pool_chains(chains, colnames(fixture$annotation))

  expect_lte(max(abs(pooled$annotation_pip - exact$annotation_pip)), 0.02)
  expect_lte(0.5 * sum(abs(pooled$model_probability - exact$model_probability)),
             0.03)
  expect_lte(max(abs(pooled$model_probability - exact$model_probability)), 0.02)
  expect_lte(max(abs(pooled$alpha_mean - exact$alpha_mean)), 0.04)
  expect_lte(max(abs(
    pooled$alpha_mean_given_inclusion - exact$alpha_mean_given_inclusion
  )), 0.04)
  expect_lte(max(abs(pooled$intercept_mean - exact$intercept_mean)), 0.04)
  expect_lte(max(abs(pooled$q_mean - exact$q_mean)), 0.015)
  expect_lte(max(abs(pooled$pi_mean - exact$pi_mean)), 0.015)
})

test_that("zero-information annotation retains its prior inclusion probability", {
  fixture <- .sbs_fixture()
  fixture$annotation[, "null_annotation"] <- 0
  exact <- .sbs_exact_posterior(
    fixture$z, fixture$annotation, fixture$eligible,
    fixture$pi_a, fixture$tau2
  )
  expect_equal(
    unname(exact$annotation_pip["null_annotation"]),
    fixture$pi_a, tolerance = 1e-13
  )
  for (stick in seq_along(fixture$z)) {
    rows <- fixture$eligible[[stick]]
    expect_equal(
      .sbs_log_bf(
        fixture$annotation[rows, "null_annotation"],
        fixture$z[[stick]], fixture$tau2[stick]
      ),
      0, tolerance = 1e-15
    )
  }

  set.seed(20260920)
  chain <- .sbs_mcmc_chain(
    fixture$z, fixture$annotation, fixture$eligible,
    fixture$pi_a, fixture$tau2,
    iterations = 8000L, burn = 1000L,
    initial_delta = c(0L, 0L, 0L)
  )
  expect_equal(mean(chain$delta_draws[, 3L]), fixture$pi_a, tolerance = 0.025)
})

test_that("annotation permutations preserve the exact and MCMC posterior", {
  fixture <- .sbs_fixture()
  permutation <- c(3L, 1L, 2L)
  exact <- .sbs_exact_posterior(
    fixture$z, fixture$annotation, fixture$eligible,
    fixture$pi_a, fixture$tau2
  )
  exact_permuted <- .sbs_exact_posterior(
    fixture$z, fixture$annotation[, permutation, drop = FALSE],
    fixture$eligible, fixture$pi_a, fixture$tau2
  )
  expect_equal(
    unname(exact_permuted$annotation_pip),
    unname(exact$annotation_pip[permutation]), tolerance = 1e-13
  )

  set.seed(20260930)
  permuted_chain <- .sbs_mcmc_chain(
    fixture$z, fixture$annotation[, permutation, drop = FALSE],
    fixture$eligible, fixture$pi_a, fixture$tau2,
    iterations = 7000L, burn = 1000L,
    initial_delta = c(1L, 0L, 1L)
  )
  expect_equal(
    colMeans(permuted_chain$delta_draws),
    unname(exact$annotation_pip[permutation]), tolerance = 0.025
  )
})

test_that("duplicate annotations have exchange-symmetric inclusion posteriors", {
  fixture <- .sbs_fixture()
  duplicate <- cbind(
    copy_a = fixture$annotation[, "continuous_signal"],
    copy_b = fixture$annotation[, "continuous_signal"],
    null_annotation = fixture$annotation[, "null_annotation"]
  )
  exact <- .sbs_exact_posterior(
    fixture$z, duplicate, fixture$eligible, fixture$pi_a, fixture$tau2
  )
  expect_equal(unname(exact$annotation_pip["copy_a"]),
               unname(exact$annotation_pip["copy_b"]), tolerance = 1e-14)

  set.seed(20260940)
  chain <- .sbs_mcmc_chain(
    fixture$z, duplicate, fixture$eligible,
    fixture$pi_a, fixture$tau2,
    iterations = 8000L, burn = 1000L,
    initial_delta = c(0L, 1L, 0L)
  )
  expect_equal(mean(chain$delta_draws[, 1L]), mean(chain$delta_draws[, 2L]),
               tolerance = 0.04)
})
