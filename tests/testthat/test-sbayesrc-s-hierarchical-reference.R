test_that("learned pi_A beta and collapsed beta-binomial conditionals agree", {
  delta <- c(1L, 0L, 1L, 0L, 0L)
  a_pi <- 1.5
  b_pi <- 4
  annotation <- 2L
  included_other <- sum(delta[-annotation])
  expected <- log(a_pi + included_other) -
    log(b_pi + length(delta) - 1L - included_other)
  expect_equal(
    .sbs3_collapsed_log_odds(delta, annotation, a_pi, b_pi),
    expected, tolerance = 1e-15
  )

  set.seed(20270101)
  draws <- replicate(30000L, .sbs3_draw_pi(delta, a_pi, b_pi))
  expect_equal(
    mean(draws),
    (a_pi + sum(delta)) / (a_pi + b_pi + length(delta)),
    tolerance = 0.006
  )
})

test_that("fixed-z learned-pi_A chains reproduce beta-binomial enumeration", {
  fixture <- .sbs_fixture()
  a_pi <- 1
  b_pi <- 1
  exact <- .sbs3_exact_pi_posterior(
    fixture$z, fixture$annotation, fixture$eligible, fixture$tau2,
    a_pi, b_pi
  )
  expect_equal(sum(exact$model_probability), 1, tolerance = 1e-14)
  expect_equal(exact$expected_included, sum(exact$annotation_pip),
               tolerance = 1e-15)

  initial <- list(c(0L, 0L, 0L), c(1L, 1L, 1L),
                  c(1L, 0L, 1L), c(0L, 1L, 0L))
  run_route <- function(collapsed, seed_offset) lapply(seq_along(initial), function(i) {
    set.seed(seed_offset + i)
    .sbs3_run_chain(
      fixture$annotation, fixture$eligible, 5000L, 1000L,
      initial[[i]], fixture$pi_a, fixture$tau2,
      fixed_z = fixture$z, learn_pi = TRUE, collapsed_pi = collapsed,
      a_pi = a_pi, b_pi = b_pi
    )
  })
  explicit <- .sbs3_pool(run_route(FALSE, 20270110L),
                         colnames(fixture$annotation))
  collapsed <- .sbs3_pool(run_route(TRUE, 20270120L),
                          colnames(fixture$annotation))

  expect_lte(max(abs(explicit$annotation_pip - exact$annotation_pip)), 0.025)
  expect_lte(max(abs(collapsed$annotation_pip - exact$annotation_pip)), 0.025)
  expect_lte(max(abs(explicit$annotation_pip - collapsed$annotation_pip)), 0.03)
  expect_equal(mean(explicit$included), sum(explicit$annotation_pip),
               tolerance = 1e-14)
  expect_equal(mean(collapsed$included), sum(collapsed$annotation_pip),
               tolerance = 1e-14)
})

test_that("learned-pi_A observed-d explicit and collapsed routes agree", {
  fixture <- .sbs2_fixture(observations = 150L, seed = 20270130L)
  set.seed(20270131L)
  explicit <- .sbs3_run_chain(
    fixture$annotation, fixture$eligible, 3500L, 500L,
    c(0L, 0L, 0L), fixture$pi_a, fixture$tau2,
    outcome = fixture$outcome, learn_pi = TRUE,
    a_pi = 1, b_pi = 1
  )
  set.seed(20270132L)
  collapsed <- .sbs3_run_chain(
    fixture$annotation, fixture$eligible, 3500L, 500L,
    c(1L, 1L, 1L), fixture$pi_a, fixture$tau2,
    outcome = fixture$outcome, learn_pi = TRUE, collapsed_pi = TRUE,
    a_pi = 1, b_pi = 1
  )
  expect_lte(max(abs(
    colMeans(explicit$delta_draws) - colMeans(collapsed$delta_draws)
  )), 0.05)
  expect_lte(max(abs(explicit$q_mean - collapsed$q_mean)), 0.035)
  expect_lte(max(abs(explicit$pi_mean - collapsed$pi_mean)), 0.035)
})

test_that("tau2 inverse-gamma conditional and empty-model reversion are exact", {
  alpha <- c(0.4, 0, -0.7, 0.2)
  delta <- c(1L, 0L, 1L, 1L)
  conditional <- .sbs3_tau_conditional(alpha, delta, 3, 1.6)
  expect_equal(conditional$shape, 4.5)
  expect_equal(conditional$scale, 1.6 + 0.5 * sum(alpha[delta == 1L]^2))

  set.seed(20270201)
  draws <- replicate(40000L, .sbs3_draw_tau2(alpha, delta, 3, 1.6))
  expect_equal(
    mean(draws), conditional$scale / (conditional$shape - 1),
    tolerance = 0.01
  )
  expected_quantile <- 1 / stats::qgamma(
    c(0.975, 0.5, 0.025), conditional$shape, rate = conditional$scale
  )
  expect_equal(unname(stats::quantile(draws, c(0.025, 0.5, 0.975))),
               expected_quantile, tolerance = 0.025)

  empty <- .sbs3_tau_conditional(rep(0, 4L), rep(0L, 4L), 3, 1.6)
  expect_identical(empty, list(shape = 3, scale = 1.6))
})

test_that("fixed-z learned tau2 chain agrees with one-dimensional quadrature", {
  fixture <- .sbs_fixture()
  annotation <- fixture$annotation[, "continuous_signal", drop = FALSE]
  eligible <- fixture$eligible[1L]
  z <- fixture$z[1L]
  pi_a <- 0.35
  oracle <- .sbs3_tau_quadrature(
    z, annotation, eligible, pi_a, a_tau = 3, b_tau = 1.6
  )

  set.seed(20270210)
  chain <- .sbs3_run_chain(
    annotation, eligible, 18000L, 3000L, 0L, pi_a, 0.8,
    fixed_z = z, learn_tau = TRUE, a_tau = 3, b_tau = 1.6
  )
  expect_equal(mean(chain$delta_draws), oracle$pip, tolerance = 0.03)
  expect_equal(mean(chain$tau2_draws), oracle$mean_tau, tolerance = 0.09)
})

test_that("all-excluded learned tau2 follows its prior and probabilities normalize", {
  fixture <- .sbs2_fixture(observations = 120L, seed = 20270220L)
  set.seed(20270221L)
  chain <- .sbs3_run_chain(
    fixture$annotation, fixture$eligible, 5000L, 500L,
    rep(0L, 3L), fixture$pi_a, fixture$tau2,
    outcome = fixture$outcome, learn_tau = TRUE,
    a_tau = 3, b_tau = 1.6, fixed_delta = rep(0L, 3L)
  )
  expect_identical(unique(as.numeric(chain$alpha_draws)), 0)
  expect_equal(colMeans(chain$tau2_draws), rep(0.8, 3L), tolerance = 0.04)
  expect_equal(rowSums(chain$pi_mean), rep(1, nrow(fixture$annotation)),
               tolerance = 1e-12)
})

test_that("prior predictive inclusion count and BFDR reference are correct", {
  set.seed(20270230L)
  prior <- .sbs3_prior_predictive(
    50000L, annotation_count = 20L, a_pi = 1, b_pi = 9,
    a_tau = 3, b_tau = 1.6
  )
  expect_equal(mean(prior$included), 2, tolerance = 0.04)
  expect_equal(mean(prior$pi_a), 0.1, tolerance = 0.003)
  expect_equal(colMeans(prior$tau2), rep(0.8, 3L), tolerance = 0.025)
  expect_equal(.sbs3_bfdr(c(0.95, 0.8, 0.2), 0.8), 0.125)
  expect_true(is.na(.sbs3_bfdr(c(0.2, 0.3), 0.8)))
})

test_that("joint learned pi_A and tau2 explicit and collapsed routes agree", {
  fixture <- .sbs2_fixture(observations = 170L, seed = 20270301L)
  set.seed(20270302L)
  explicit <- .sbs3_run_chain(
    fixture$annotation, fixture$eligible, 5000L, 750L,
    c(0L, 0L, 0L), fixture$pi_a, fixture$tau2,
    outcome = fixture$outcome,
    learn_pi = TRUE, learn_tau = TRUE,
    a_pi = 1, b_pi = 1, a_tau = 3, b_tau = 1.6
  )
  set.seed(20270303L)
  collapsed <- .sbs3_run_chain(
    fixture$annotation, fixture$eligible, 5000L, 750L,
    c(1L, 1L, 1L), fixture$pi_a, fixture$tau2,
    outcome = fixture$outcome,
    learn_pi = TRUE, collapsed_pi = TRUE, learn_tau = TRUE,
    a_pi = 1, b_pi = 1, a_tau = 3, b_tau = 1.6
  )
  expect_lte(max(abs(
    colMeans(explicit$delta_draws) - colMeans(collapsed$delta_draws)
  )), 0.06)
  expect_lte(max(abs(explicit$q_mean - collapsed$q_mean)), 0.04)
  expect_lte(max(abs(explicit$pi_mean - collapsed$pi_mean)), 0.04)
  expect_lte(max(abs(
    colMeans(explicit$tau2_draws) - colMeans(collapsed$tau2_draws)
  )), 0.12)
  expect_equal(rowSums(explicit$pi_mean),
               rep(1, nrow(fixture$annotation)), tolerance = 1e-12)
})

test_that("full hierarchy all-excluded state has exact conditional behavior", {
  fixture <- .sbs2_fixture(observations = 100L, seed = 20270310L)
  set.seed(20270311L)
  excluded <- .sbs3_run_chain(
    fixture$annotation, fixture$eligible, 6000L, 1000L,
    rep(0L, 3L), fixture$pi_a, fixture$tau2,
    outcome = fixture$outcome,
    learn_pi = TRUE, learn_tau = TRUE,
    a_pi = 1, b_pi = 9, a_tau = 3, b_tau = 1.6,
    fixed_delta = rep(0L, 3L)
  )
  expect_identical(unique(as.numeric(excluded$alpha_draws)), 0)
  expect_lt(abs(mean(excluded$pi_a_draws) - 1 / 13), 0.003)
  expect_equal(colMeans(excluded$tau2_draws), rep(0.8, 3L), tolerance = 0.04)
  expect_identical(unique(excluded$included_draws), 0L)
})

test_that("learned hierarchy preserves permutation and duplicate symmetry", {
  fixture <- .sbs_fixture()
  exact <- .sbs3_exact_pi_posterior(
    fixture$z, fixture$annotation, fixture$eligible, fixture$tau2, 1, 9
  )
  permutation <- c(3L, 1L, 2L)
  permuted <- .sbs3_exact_pi_posterior(
    fixture$z, fixture$annotation[, permutation, drop = FALSE],
    fixture$eligible, fixture$tau2, 1, 9
  )
  expect_equal(unname(permuted$annotation_pip),
               unname(exact$annotation_pip[permutation]), tolerance = 1e-13)

  duplicate <- cbind(
    copy_a = fixture$annotation[, "continuous_signal"],
    copy_b = fixture$annotation[, "continuous_signal"],
    null = fixture$annotation[, "null_annotation"]
  )
  duplicate_exact <- .sbs3_exact_pi_posterior(
    fixture$z, duplicate, fixture$eligible, fixture$tau2, 1, 9
  )
  expect_equal(unname(duplicate_exact$annotation_pip[1L]),
               unname(duplicate_exact$annotation_pip[2L]), tolerance = 1e-14)

  zero <- rep(0, nrow(fixture$annotation))
  for (stick in seq_along(fixture$eligible)) {
    expect_equal(.sbs_log_bf(
      zero[fixture$eligible[[stick]]], fixture$z[[stick]], fixture$tau2[stick]
    ), 0, tolerance = 1e-15)
  }
})
