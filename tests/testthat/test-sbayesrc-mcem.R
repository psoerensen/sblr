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

test_that("MCEM damping and convergence contracts are explicit", {
  alpha <- matrix(c(-1, 0, 1, 2), 2L, 2L)
  target <- matrix(c(1, 2, 3, 4), 2L, 2L)
  expect_equal(
    sblr:::.sbayesrc_mcem_damped_update(alpha, target, 0.5),
    (alpha + target) / 2
  )
  expect_equal(
    sblr:::.sbayesrc_mcem_damped_update(alpha, target, 1), target
  )
  expect_false(sblr:::.sbayesrc_mcem_has_converged(
    2, 3, 0, 0, 1e-3, 1e-3
  ))
  expect_false(sblr:::.sbayesrc_mcem_has_converged(
    3, 3, 2e-3, 0, 1e-3, 1e-3
  ))
  expect_true(sblr:::.sbayesrc_mcem_has_converged(
    3, 3, 5e-4, 5e-4, 1e-3, 1e-3
  ))
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
  expect_identical(fit$mcem$method, "SBayesRC-EM")
  expect_identical(fit$mcem$algorithm, "MCEM")
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

test_that("block-eigen RB collection is RNG neutral", {
 fixture <- .mcem_block_fixture(marker_count = 8L, sample_count = 60L)
 on.exit(.mcem_cleanup_block_fixture(fixture), add = TRUE)
 prior <- list(
  distribution = "normal", mean = fixture$intercept_prior["mean", ],
  sd = 1 / sqrt(fixture$intercept_prior["precision", ])
 )
 common <- list(
  stats = fixture$stats, Glist = fixture$Glist,
  annotation = fixture$A, block_start = 1L,
  representation = "low_rank", eigen_prop = 0.999999,
  residual_policy = "gctb_block", block_ve_mode = "fixVe",
  gamma = fixture$gamma, B = fixture$B, E = fixture$E,
  ssb_prior = fixture$ssb_prior[[1L]], sse_prior = fixture$sse_prior[[1L]],
  updateAlpha = FALSE, updateB = FALSE, updateE = FALSE,
  nit = 40L, nburn = 15L, ncores = 1L, seed = 771L,
  keep_chains = TRUE, comp_init = list(rep(0L, nrow(fixture$A))),
  use_comp_init = TRUE, add_intercept = FALSE,
  standardize_annotations = FALSE,
  alpha_init = matrix(c(-0.5, -0.4, 0.1, -0.1), 2L, 2L),
  sigmaSqAlpha_init = c(1, 1), annotation_intercept_prior = prior,
  .diagnostic_updateSigmaSqAlpha = FALSE, .return_raw = TRUE
 )
 off <- do.call(sblr:::.stblr_csr_sbayesrc_block_eigen,
                c(common, list(.information_diagnostics = FALSE)))
 on <- do.call(sblr:::.stblr_csr_sbayesrc_block_eigen,
               c(common, list(.information_diagnostics = TRUE)))
 expect_identical(on$marker, off$marker)
 expect_identical(on$trace, off$trace)
 expect_identical(on$variance, off$variance)
 expect_identical(on$pi, off$pi)
 expect_identical(on$annotation, off$annotation)
 expect_identical(on$component, off$component)
 rb <- on$chains[[1L]][[1L]]$information_flow$rb_comp_prob
 expect_equal(rowSums(rb), rep(1, nrow(rb)), tolerance = 1e-12)
 expect_identical(on$meta$residual_policy, "gctb_block")
 expect_identical(on$meta$block_ve_mode, "fixVe")
})

test_that("fixed-parameter block-eigen MCEM agrees with CSR reference", {
 fixture <- .mcem_block_fixture(marker_count = 12L, sample_count = 80L)
 on.exit(.mcem_cleanup_block_fixture(fixture), add = TRUE)
 starts <- list(
  baseline = matrix(c(-0.6, -0.5, 0, 0), 2L, 2L),
  dispersed = matrix(c(0.2, 0.1, -0.6, 0.5), 2L, 2L)
 )
 run_csr <- function(alpha_start, seed) {
  sblr:::.stblr_mcem_sbayesrc_csr(
   wy = fixture$stats$wy, ww = fixture$stats$ww, yy = fixture$stats$yy,
   b_init = list(rep(0, nrow(fixture$A))),
   comp_init = list(rep(0L, nrow(fixture$A))),
   r_init = fixture$stats$wy, ld_prefix = fixture$prefix,
   B = fixture$B, E = fixture$E,
   ssb_prior = fixture$ssb_prior, sse_prior = fixture$sse_prior,
   A = fixture$A, gamma = fixture$gamma, alpha_init = alpha_start,
   sigmaSqAlpha_init = c(1, 1),
   intercept_prior_resolved = fixture$intercept_prior,
   n = fixture$stats$n, inner_sweeps = 800L, inner_burn = 300L,
   final_sweeps = 900L, final_burn = 350L, max_outer = 35L,
   seed = seed, ncores = 1L
  )
 }
 run_block <- function(alpha_start, seed) {
  sblr:::.stblr_mcem_sbayesrc_block_eigen(
   stats = fixture$stats, Glist = fixture$Glist,
   annotation = fixture$A, block_start = 1L,
   B = fixture$B, E = fixture$E,
   ssb_prior = fixture$ssb_prior, sse_prior = fixture$sse_prior,
   gamma = fixture$gamma, alpha_init = alpha_start,
   sigmaSqAlpha_init = c(1, 1),
   intercept_prior_resolved = fixture$intercept_prior,
   representation = "low_rank", eigen_prop = 0.999999,
   residual_policy = "gctb_block", block_ve_mode = "fixVe",
   inner_sweeps = 800L, inner_burn = 300L,
   final_sweeps = 900L, final_burn = 350L, max_outer = 35L,
   seed = seed, ncores = 1L
  )
 }
 csr <- run_csr(starts$baseline, 881L)
 block <- run_block(starts$baseline, 882L)
 block_dispersed <- run_block(starts$dispersed, 883L)
 expect_lt(max(abs(csr$mcem$alpha_map - block$mcem$alpha_map)), 0.12)
 expect_lt(max(abs(csr$mcem$component_prior - block$mcem$component_prior)), 0.04)
 expect_lt(max(abs(csr$mcem$last_estep_responsibilities -
                   block$mcem$last_estep_responsibilities)), 0.08)
 expect_lt(abs(sum(1 - csr$mcem$component_prior[, 1L]) -
               sum(1 - block$mcem$component_prior[, 1L])), 0.4)
 expect_lt(max(abs(block$mcem$alpha_map - block_dispersed$mcem$alpha_map)), 0.12)
  expect_identical(block$mcem$backend, "block_eigen")
  expect_identical(block$mcem$method, "SBayesRC-EM")
  expect_identical(block$mcem$algorithm, "MCEM")
  expect_identical(block$mcem$sigmaSqAlpha_mode, "fixed_prior_variance")
  expect_identical(
   block$mcem$mixture_prior_mode,
   "annotation_stick_intercepts_no_global_Pi_update"
  )
  expect_equal(
   block$genomic$annotation$alpha_final[[1L]], block$mcem$alpha_map,
   tolerance = 1e-12
  )
  expect_equal(
   rowSums(block$mcem$final_genomic_responsibilities),
   rep(1, nrow(block$mcem$component_prior)), tolerance = 1e-12
  )
  expect_false(identical(
   block$mcem$last_estep_responsibilities,
   block$mcem$final_genomic_responsibilities
  ))
  expect_identical(block$genomic$meta$residual_policy, "gctb_block")
 expect_identical(block$genomic$meta$block_ve_mode, "fixVe")
})

test_that("learned block residual variance remains RNG neutral and stable", {
 fixture <- .mcem_block_fixture(seed = 9190L, marker_count = 12L,
                                sample_count = 90L)
 on.exit(.mcem_cleanup_block_fixture(fixture), add = TRUE)
 prior <- list(
  distribution = "normal", mean = fixture$intercept_prior["mean", ],
  sd = 1 / sqrt(fixture$intercept_prior["precision", ])
 )
 raw_args <- list(
  stats = fixture$stats, Glist = fixture$Glist,
  annotation = fixture$A, block_start = 1L,
  representation = "low_rank", eigen_prop = 0.999999,
  residual_policy = "gctb_block", block_ve_mode = "allMixVe",
  gamma = fixture$gamma, B = fixture$B, E = fixture$E,
  ssb_prior = fixture$ssb_prior[[1L]], sse_prior = fixture$sse_prior[[1L]],
  updateAlpha = FALSE, updateB = FALSE, updateE = TRUE,
  nit = 60L, nburn = 20L, ncores = 1L, seed = 772L,
  keep_chains = TRUE, comp_init = list(rep(0L, nrow(fixture$A))),
  use_comp_init = TRUE, add_intercept = FALSE,
  standardize_annotations = FALSE,
  alpha_init = matrix(c(-0.5, -0.4, 0.1, -0.1), 2L, 2L),
  sigmaSqAlpha_init = c(1, 1), annotation_intercept_prior = prior,
  .diagnostic_updateSigmaSqAlpha = FALSE, .return_raw = TRUE
 )
 off <- do.call(sblr:::.stblr_csr_sbayesrc_block_eigen,
                c(raw_args, list(.information_diagnostics = FALSE)))
 on <- do.call(sblr:::.stblr_csr_sbayesrc_block_eigen,
               c(raw_args, list(.information_diagnostics = TRUE)))
 expect_identical(on$marker, off$marker)
 expect_identical(on$trace, off$trace)
 expect_identical(on$variance, off$variance)
 expect_identical(on$diagnostics$block_residual,
                  off$diagnostics$block_residual)

 baseline <- matrix(c(-0.6, -0.5, 0, 0), 2L, 2L)
 dispersed <- matrix(c(0.2, 0.1, -0.6, 0.5), 2L, 2L)
 fixed <- .mcem_run_block(fixture, baseline, 991L)
 learned <- .mcem_run_block(fixture, baseline, 992L, updateE = TRUE)
 learned_dispersed <- .mcem_run_block(
  fixture, dispersed, 993L, updateE = TRUE
 )
 expect_true(all(is.finite(learned$mcem$history$summary$E)))
 expect_true(all(learned$mcem$history$summary$E > 0))
 expect_lt(max(abs(learned$mcem$alpha_map -
                   learned_dispersed$mcem$alpha_map)), 0.15)
 expect_lt(max(abs(fixed$mcem$alpha_map - learned$mcem$alpha_map)), 0.35)
 expect_true(isTRUE(learned$mcem$genomic_hyperparameters$updateE))
 expect_identical(learned$genomic$meta$block_ve_mode, "allMixVe")
})

test_that("learned B and joint B/E produce stable SBayesRC-EM solutions", {
 fixture <- .mcem_block_fixture(seed = 9195L, marker_count = 12L,
                                sample_count = 90L)
 on.exit(.mcem_cleanup_block_fixture(fixture), add = TRUE)
 baseline <- matrix(c(-0.6, -0.5, 0, 0), 2L, 2L)
 dispersed <- matrix(c(0.2, 0.1, -0.6, 0.5), 2L, 2L)
 learned_b <- .mcem_run_block(
  fixture, baseline, 1101L, updateB = TRUE, updateE = FALSE
 )
 learned_b_dispersed <- .mcem_run_block(
  fixture, dispersed, 1102L, updateB = TRUE, updateE = FALSE
 )
 learned_both <- .mcem_run_block(
  fixture, baseline, 1103L, updateB = TRUE, updateE = TRUE
 )
 learned_both_dispersed <- .mcem_run_block(
  fixture, dispersed, 1104L, updateB = TRUE, updateE = TRUE
 )
 for (fit in list(learned_b, learned_b_dispersed,
                  learned_both, learned_both_dispersed)) {
  expect_true(all(is.finite(fit$mcem$history$summary$B)))
  expect_true(all(is.finite(fit$mcem$history$summary$E)))
  expect_true(all(fit$mcem$history$summary$B > 0))
  expect_true(all(fit$mcem$history$summary$E > 0))
  expect_true(all(is.finite(fit$mcem$alpha_map)))
  expect_equal(rowSums(fit$mcem$component_prior),
               rep(1, nrow(fit$mcem$component_prior)), tolerance = 1e-10)
 }
 expect_lt(max(abs(learned_b$mcem$alpha_map -
                   learned_b_dispersed$mcem$alpha_map)), 0.18)
 expect_lt(max(abs(learned_both$mcem$alpha_map -
                   learned_both_dispersed$mcem$alpha_map)), 0.18)
 expect_true(isTRUE(learned_b$mcem$genomic_hyperparameters$updateB))
 expect_false(isTRUE(learned_b$mcem$genomic_hyperparameters$updateE))
 expect_true(isTRUE(learned_both$mcem$genomic_hyperparameters$updateB))
 expect_true(isTRUE(learned_both$mcem$genomic_hyperparameters$updateE))
})
