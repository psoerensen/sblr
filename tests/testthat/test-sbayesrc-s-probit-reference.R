test_that("Albert-Chib inverse-CDF draws satisfy finite truncation contracts", {
  eta <- rep(c(-6, -2, 0, 2, 6), each = 200L)
  set.seed(20261010)
  positive <- .sbs2_rtrunc_probit(eta, rep(1L, length(eta)))
  set.seed(20261011)
  negative <- .sbs2_rtrunc_probit(eta, rep(0L, length(eta)))

  expect_true(all(is.finite(positive)))
  expect_true(all(is.finite(negative)))
  expect_true(all(positive > 0))
  expect_true(all(negative <= 0))
})

test_that("observed-stick eligibility uses identical nested rows", {
  fixture <- .sbs2_fixture(observations = 180L)
  expect_silent(.sbs2_validate_hierarchy(
    fixture$outcome, fixture$annotation, fixture$eligible
  ))
  expect_identical(fixture$eligible[[1L]], seq_len(nrow(fixture$annotation)))
  expect_identical(
    fixture$eligible[[2L]],
    fixture$eligible[[1L]][fixture$outcome[[1L]] == 1L]
  )
  expect_identical(
    fixture$eligible[[3L]],
    fixture$eligible[[2L]][fixture$outcome[[2L]] == 1L]
  )
  broken <- fixture$eligible
  broken[[2L]] <- rev(broken[[2L]])
  expect_error(.sbs2_validate_hierarchy(
    fixture$outcome, fixture$annotation, broken
  ))
})

test_that("direct and rank-one collapsed Bayes factors agree", {
  fixture <- .sbs2_fixture(observations = 160L)
  rows <- fixture$eligible[[2L]]
  x <- fixture$annotation[rows, 2L]
  residual <- fixture$latent_true[[2L]] - fixture$true_intercept[2L]
  expect_equal(
    .sbs_log_bf(x, residual, fixture$tau2[2L]),
    .sbs2_log_bf_direct(x, residual, fixture$tau2[2L]),
    tolerance = 1e-13
  )
})

test_that("all-included SBayesRC-S reduces to continuous-alpha SBayesRC", {
  fixture <- .sbs2_fixture(observations = 150L)
  set.seed(20261020)
  selection <- .sbs2_run_chain(
    fixture$outcome, fixture$annotation, fixture$eligible,
    fixture$pi_a, fixture$tau2, iterations = 2500L, burn = 500L,
    initial_delta = rep(1L, 3L), fixed_delta = rep(1L, 3L)
  )
  set.seed(20261021)
  continuous <- .sbs2_standard_continuous_chain(
    fixture$outcome, fixture$annotation, fixture$eligible, fixture$tau2,
    iterations = 2500L, burn = 500L
  )
  selection_coefficients <- array(NA_real_, c(2000L, 4L, 3L))
  selection_coefficients[, 1L, ] <- selection$intercept_draws
  selection_coefficients[, -1L, ] <- selection$alpha_draws

  selection_mean <- apply(selection_coefficients, c(2L, 3L), mean)
  continuous_mean <- apply(continuous$coefficient_draws, c(2L, 3L), mean)
  selection_sd <- apply(selection_coefficients, c(2L, 3L), stats::sd)
  continuous_sd <- apply(continuous$coefficient_draws, c(2L, 3L), stats::sd)
  expect_lte(max(abs(selection_mean - continuous_mean)), 0.09)
  expect_lte(max(abs(selection_sd - continuous_sd)), 0.08)
  expect_lte(max(abs(selection$q_mean - continuous$q_mean)), 0.03)
  expect_lte(max(abs(selection$pi_mean - continuous$pi_mean)), 0.03)
})

test_that("all-excluded SBayesRC-S is an intercept-only hierarchy", {
  fixture <- .sbs2_fixture(observations = 180L)
  set.seed(20261030)
  excluded <- .sbs2_run_chain(
    fixture$outcome, fixture$annotation, fixture$eligible,
    fixture$pi_a, fixture$tau2, iterations = 1500L, burn = 300L,
    initial_delta = rep(0L, 3L), fixed_delta = rep(0L, 3L)
  )
  expect_identical(unique(as.numeric(excluded$alpha_draws)), 0)
  expected_q <- vapply(seq_len(3L), function(stick) {
    mean(stats::pnorm(excluded$intercept_draws[, stick]))
  }, numeric(1L))
  expect_equal(colMeans(excluded$q_mean), expected_q, tolerance = 1e-12)
  expect_equal(rowSums(excluded$pi_mean), rep(1, nrow(fixture$annotation)),
               tolerance = 1e-12)
})

test_that("zero observed-d annotation retains prior inclusion probability", {
  fixture <- .sbs2_fixture(observations = 150L)
  fixture$annotation[, 3L] <- 0
  set.seed(20261040)
  chain <- .sbs2_run_chain(
    fixture$outcome, fixture$annotation, fixture$eligible,
    fixture$pi_a, fixture$tau2, iterations = 4000L, burn = 500L,
    initial_delta = rep(0L, 3L)
  )
  expect_lte(abs(mean(chain$delta_draws[, 3L]) - fixture$pi_a), 0.04)
})

test_that("observed-d annotation permutations preserve posterior summaries", {
  fixture <- .sbs2_fixture(observations = 150L)
  permutation <- c(3L, 1L, 2L)
  set.seed(20261050)
  original <- .sbs2_run_chain(
    fixture$outcome, fixture$annotation, fixture$eligible,
    fixture$pi_a, fixture$tau2, iterations = 3500L, burn = 500L,
    initial_delta = c(0L, 1L, 0L)
  )
  set.seed(20261051)
  permuted <- .sbs2_run_chain(
    fixture$outcome, fixture$annotation[, permutation, drop = FALSE],
    fixture$eligible, fixture$pi_a, fixture$tau2,
    iterations = 3500L, burn = 500L,
    initial_delta = c(0L, 0L, 1L)
  )
  expect_lte(max(abs(
    colMeans(permuted$delta_draws) -
      colMeans(original$delta_draws)[permutation]
  )), 0.06)
  permuted_alpha <- apply(permuted$alpha_draws, c(2L, 3L), mean)
  original_alpha <- apply(original$alpha_draws, c(2L, 3L), mean)
  expect_lte(max(abs(
    permuted_alpha - original_alpha[permutation, , drop = FALSE]
  )), 0.10)
})

test_that("primary and direct observed-d reference samplers agree", {
  fixture <- .sbs2_fixture(observations = 160L)
  set.seed(20261060)
  primary <- .sbs2_run_chain(
    fixture$outcome, fixture$annotation, fixture$eligible,
    fixture$pi_a, fixture$tau2, iterations = 4000L, burn = 500L,
    initial_delta = c(0L, 0L, 0L), method = "primary"
  )
  set.seed(20261061)
  direct <- .sbs2_run_chain(
    fixture$outcome, fixture$annotation, fixture$eligible,
    fixture$pi_a, fixture$tau2, iterations = 4000L, burn = 500L,
    initial_delta = c(1L, 1L, 1L), method = "direct"
  )
  expect_lte(max(abs(
    colMeans(primary$delta_draws) - colMeans(direct$delta_draws)
  )), 0.05)
  expect_lte(max(abs(primary$q_mean - direct$q_mean)), 0.03)
  expect_lte(max(abs(primary$pi_mean - direct$pi_mean)), 0.03)
  expect_equal(rowSums(primary$pi_mean), rep(1, nrow(fixture$annotation)),
               tolerance = 1e-12)
})
