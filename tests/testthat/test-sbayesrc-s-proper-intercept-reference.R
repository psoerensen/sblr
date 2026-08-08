test_that("proper intercept marginal agrees with dense Gaussian integration", {
  fixture <- .sbs_fixture()
  prior <- fixture$intercept_prior
  for (stick in seq_along(fixture$z)) {
    rows <- fixture$eligible[[stick]]
    for (selected in list(c(FALSE, FALSE, FALSE), c(TRUE, FALSE, TRUE))) {
      analytic <- .sbs_stick_model(
        fixture$z[[stick]], fixture$annotation[rows, , drop = FALSE],
        selected, fixture$tau2[stick], prior$mean[stick],
        prior$variance[stick]
      )$log_marginal
      dense <- .sbs_stick_log_marginal_dense(
        fixture$z[[stick]], fixture$annotation[rows, , drop = FALSE],
        selected, fixture$tau2[stick], prior$mean[stick],
        prior$variance[stick]
      )
      expect_equal(analytic, dense, tolerance = 1e-10)
    }
  }
})

test_that("proper intercept prior is mixture-centred with controlled probit tails", {
  prior <- .sbs_intercept_prior(3L, rep(0.25, 4L), sd = 1)
  expect_equal(prior$stick_probability, c(0.75, 2 / 3, 0.5),
               tolerance = 1e-14)
  expect_equal(stats::pnorm(prior$mean), prior$stick_probability,
               tolerance = 1e-14)
  predictive <- .sbs_intercept_prior_predictive(prior)
  expect_true(all(is.finite(predictive)))
  expect_true(all(predictive[, "probability_below_001"] < 0.02))
  expect_true(all(predictive[, "probability_above_099"] < 0.06))
})

test_that("empty sticks reduce exactly to proper coefficient priors", {
  prior <- .sbs_intercept_prior(3L)
  empty_annotation <- matrix(numeric(), 0L, 3L)
  posterior <- .sbs_stick_model(
    numeric(), empty_annotation, c(TRUE, FALSE, TRUE), 0.8,
    prior$mean[3L], prior$variance[3L]
  )
  expect_equal(posterior$mean, c(prior$mean[3L], 0, 0), tolerance = 1e-14)
  expect_equal(diag(posterior$covariance),
               c(prior$variance[3L], 0.8, 0.8), tolerance = 1e-14)
  expect_equal(.sbs_log_bf(numeric(), numeric(), 0.8), 0,
               tolerance = 1e-14)
})

test_that("observed hierarchy supports empty and repopulated later sticks", {
  set.seed(20271011)
  annotation <- cbind(a = stats::rnorm(24), b = stats::rnorm(24))
  prior <- .sbs_intercept_prior(3L)
  empty <- list(0:0, 0:0, integer())
  empty[[1L]] <- rep(0L, 24L)
  empty[[2L]] <- integer()
  eligible_empty <- list(seq_len(24L), integer(), integer())
  empty_chain <- .sbs3_run_chain(
    annotation, eligible_empty, 300L, 50L, c(0L, 1L), 0.35,
    rep(0.8, 3L), outcome = empty, learn_pi = TRUE, learn_tau = TRUE,
    intercept_prior = prior
  )
  expect_true(all(is.finite(c(
    empty_chain$intercept_draws, empty_chain$alpha_draws,
    empty_chain$pi_a_draws, empty_chain$tau2_draws,
    empty_chain$q_mean, empty_chain$pi_mean
  ))))
  expect_equal(rowSums(empty_chain$pi_mean), rep(1, nrow(annotation)),
               tolerance = 1e-12)

  # A controlled following state repopulates both later sticks. The same
  # coherent prior applies without a state-dependent branch.
  outcome_full <- list(
    rep(c(1L, 0L), 12L), rep(c(1L, 0L), 6L), rep(c(0L, 1L), 3L)
  )
  eligible_full <- list(seq_len(24L), seq.int(1L, 23L, 2L),
                        seq.int(1L, 21L, 4L))
  full_chain <- .sbs3_run_chain(
    annotation, eligible_full, 300L, 50L, c(1L, 0L), 0.35,
    rep(0.8, 3L), outcome = outcome_full, learn_pi = TRUE,
    learn_tau = TRUE, intercept_prior = prior
  )
  expect_true(all(is.finite(c(full_chain$intercept_draws,
                              full_chain$q_mean, full_chain$pi_mean))))
})
