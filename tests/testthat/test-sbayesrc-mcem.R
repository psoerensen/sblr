test_that("MCEM probability mapping and soft sticks match the R oracle", {
  set.seed(20271270L)
  A <- cbind(Intercept = 1, binary = rbinom(30, 1, 0.3),
             continuous = rnorm(30))
  alpha <- matrix(c(-0.7, 0.4, -0.2, -0.3, 0.2, 0.1), 3L, 2L)
  probability <- sblr:::.sbayesrc_mcem_component_prior(A, alpha)
  reference <- .mcem_ref_component_prior(A, alpha)
  expect_equal(unname(probability), reference, tolerance = 1e-14)
  expect_equal(rowSums(probability), rep(1, nrow(A)), tolerance = 1e-14)

  responsibility <- matrix(runif(30 * 3), 30L, 3L)
  responsibility <- responsibility / rowSums(responsibility)
  observed <- sblr:::.sbayesrc_mcem_soft_stick_information(responsibility)
  expected <- .mcem_ref_soft_stick(responsibility)
  expect_equal(unname(observed$eligible), expected$eligible, tolerance = 1e-15)
  expect_equal(unname(observed$success), expected$success, tolerance = 1e-15)
  expect_true(all(observed$success <= observed$eligible + 1e-15))
})

test_that("MCEM proper-prior probit M-step matches the independent R oracle", {
  set.seed(20271271L)
  A <- cbind(Intercept = 1, binary = rbinom(80, 1, 0.35),
             continuous = rnorm(80))
  responsibility <- matrix(rexp(80 * 4), 80L, 4L)
  responsibility <- responsibility / rowSums(responsibility)
  component_probability <- c(0.75, 0.15, 0.07, 0.03)
  intercept <- sblr:::.sbayesrc_resolve_intercept_prior(
    component_probability, list(mean = "initial_mixture", sd = c(0.8, 1, 1.2))
  )
  sigma <- c(0.7, 1.1, 1.4)
  starts <- list(
    matrix(0, 3L, 3L),
    matrix(rnorm(9, sd = 0.4), 3L, 3L)
  )
  for (start in starts) {
    observed <- sblr:::.sbayesrc_mcem_m_step(
      A, responsibility, start, intercept$native, sigma
    )
    expected <- .mcem_ref_m_step(
      A, responsibility, start, intercept$mean, intercept$sd, sqrt(sigma)
    )
    expect_equal(observed$alpha, expected, tolerance = 2e-8)
    expect_true(all(observed$convergence == 0L))
  }
})

test_that("MCEM inner RB capture is RNG-neutral and uses fixed alpha priors", {
  fixture <- .sbs4b_fixture(36L, 20271272L)
  off <- .sbs4b_run(
    fixture, 20271273L, 100L, 20L, updateB = FALSE, updateE = FALSE,
    update_hierarchy = FALSE, selection_enabled = FALSE,
    information_diagnostics = FALSE
  )
  on <- .sbs4b_run(
    fixture, 20271273L, 100L, 20L, updateB = FALSE, updateE = FALSE,
    update_hierarchy = FALSE, selection_enabled = FALSE,
    information_diagnostics = TRUE
  )
  off_chain <- off$chains[[1L]][[1L]]
  on_chain <- on$chains[[1L]][[1L]]
  for (field in c("marker", "trace", "pi", "component", "annotation",
                  "convergence_trace", "selection")) {
    expect_identical(on_chain[[field]], off_chain[[field]], info = field)
  }
  expected_prior <- sblr:::.sbayesrc_mcem_component_prior(
    fixture$A, fixture$alpha_init
  )
  expect_equal(on_chain$information_flow$prior_comp_prob, unname(expected_prior),
               tolerance = 1e-12)
  expect_equal(rowSums(on_chain$information_flow$rb_comp_prob),
               rep(1, 36L), tolerance = 1e-12)
})

test_that("MCEM reproduces the orthogonal observed-data MAP target", {
  fixture <- .mcem_exact_fixture(FALSE, 6L, 20271274L)
  on.exit(.mcem_cleanup_csr(fixture$prefix), add = TRUE)
  exact <- .mcem_exact_target(fixture)
  fit <- .mcem_run_fixture(fixture, 20271275L, 500L, 150L, 12L)
  expect_identical(fit$mcem$method, "MCEM-SBayesRC")
  expect_identical(fit$mcem$target, "observed_data_alpha_MAP_empirical_Bayes")
  expect_true(isTRUE(fit$mcem$genomic_hyperparameters_fixed))
  expect_lt(sqrt(mean((fit$mcem$alpha_map - exact$alpha)^2)), 0.12)
  expect_equal(rowSums(fit$mcem$component_prior), rep(1, 6L), tolerance = 1e-12)
  expect_s3_class(fit$genomic, "stblr_raw")
})

test_that("MCEM approaches the tiny correlated-LD exact target", {
  fixture <- .mcem_exact_fixture(TRUE, 6L, 20271276L)
  on.exit(.mcem_cleanup_csr(fixture$prefix), add = TRUE)
  exact <- .mcem_exact_target(fixture)
  fit <- .mcem_run_fixture(fixture, 20271277L, 1500L, 400L, 30L)
  alpha_rmse <- sqrt(mean((fit$mcem$alpha_map - exact$alpha)^2))
  prior_rmse <- sqrt(mean((
    fit$mcem$component_prior - .mcem_ref_component_prior(fixture$A, exact$alpha)
  )^2))
  expect_lt(alpha_rmse, 0.08)
  expect_lt(prior_rmse, 0.03)
  expect_equal(rowSums(fit$mcem$final_responsibilities), rep(1, 6L),
               tolerance = 1e-12)
})
