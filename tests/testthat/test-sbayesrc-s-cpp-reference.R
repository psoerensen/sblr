test_that("Phase-4A deterministic C++ hierarchy mathematics matches R", {
  fixture <- .sbs2_fixture(observations = 180L, seed = 20270801L)
  alpha <- .sbs4a_alpha_matrix(
    fixture$true_intercept,
    fixture$true_alpha + matrix(seq(-0.03, 0.03, length.out = 9L), 3L, 3L)
  )
  delta <- c(1L, 1L, 0L)
  reference <- .sbs4a_math_reference(
    fixture$annotation, fixture$eligible, fixture$latent_true,
    alpha, delta, 0.35, fixture$tau2, 1, 9, 3, 1.6,
    fixture$intercept_prior
  )
  native <- .st_bayesrc_selection_math(
    fixture$annotation, fixture$eligible, fixture$latent_true,
    alpha, delta, 0.35, fixture$tau2, 1, 9, 3, 1.6,
    fixture$intercept_prior$native, 1e-12
  )
  for (field in c(
    "s", "t", "log_bf", "inclusion_logit", "inclusion_probability",
    "conditional_mean", "conditional_variance", "beta_parameters",
    "ig_parameters", "intercept_conditional", "q", "component_probability"
  )) {
    expect_equal(as.numeric(native[[field]]), as.numeric(reference[[field]]),
                 tolerance = 1e-12, info = field)
  }
  for (stick in seq_along(reference$eta)) {
    expect_equal(as.numeric(native$eta[[stick]]),
                 as.numeric(reference$eta[[stick]]), tolerance = 1e-12)
  }
})

test_that("Phase-4A C++ and R observed-d hierarchies have posterior parity", {
  fixture <- .sbs2_fixture(observations = 180L, seed = 20270810L)
  initial <- list(c(0L, 0L, 0L), c(1L, 1L, 1L))
  cpp <- lapply(seq_along(initial), function(i) {
    .sbs4a_cpp_chain(fixture, 20270810L + i, 4500L, 750L, initial[[i]])
  })
  r <- lapply(seq_along(initial), function(i) {
    set.seed(20270820L + i)
    .sbs3_run_chain(
      fixture$annotation, fixture$eligible, 4500L, 750L, initial[[i]],
      0.35, fixture$tau2, outcome = fixture$outcome,
      learn_pi = TRUE, learn_tau = TRUE,
      a_pi = 1, b_pi = 1, a_tau = 3, b_tau = 1.6
    )
  })
  cpp_summary <- .sbs4a_cpp_summary(cpp)
  r_summary <- .sbs3_pool(r, colnames(fixture$annotation))
  expect_lte(max(abs(
    cpp_summary$annotation_pip - r_summary$annotation_pip
  )), 0.03)
  expect_lte(max(abs(cpp_summary$q_mean - r_summary$q_mean)), 0.02)
  expect_lte(max(abs(
    cpp_summary$component_probability_mean - r_summary$pi_mean
  )), 0.02)
  expect_lte(abs(cpp_summary$pi_a_mean - mean(r_summary$pi_a)), 0.05)
  expect_lte(max(abs(
    cpp_summary$tau2_mean - colMeans(r_summary$tau2)
  )), 0.10)
  alpha_sd <- apply(r_summary$alpha, c(2L, 3L), stats::sd)
  expect_lte(max(abs(cpp_summary$alpha_mean - r_summary$alpha_mean) /
                   pmax(alpha_sd, 0.05)), 0.15)
})

test_that("Phase-4A structural selection guards hold", {
  fixture <- .sbs2_fixture(observations = 140L, seed = 20270830L)
  fixture$annotation[, 3L] <- 0
  alpha <- .sbs4a_alpha_matrix(fixture$true_intercept, fixture$true_alpha)
  math <- .st_bayesrc_selection_math(
    fixture$annotation, fixture$eligible, fixture$latent_true,
    alpha, c(1L, 1L, 0L), 0.35, fixture$tau2,
    1, 9, 3, 1.6, fixture$intercept_prior$native, 1e-12
  )
  expect_equal(math$log_bf[3L, ], rep(0, 3L), tolerance = 1e-14)
  expect_equal(rowSums(math$component_probability),
               rep(1, nrow(fixture$annotation)), tolerance = 1e-12)

  excluded <- .sbs4a_cpp_chain(
    fixture, 20270831L, 5000L, 750L, rep(0L, 3L), rep(0L, 3L),
    a_pi = 1, b_pi = 9
  )
  expect_identical(unique(as.numeric(excluded$alpha_draws[, -1L, ])), 0)
  expect_true(all(is.finite(excluded$alpha_draws[, 1L, ])))
  expect_equal(colMeans(excluded$tau2_draws), rep(0.8, 3L), tolerance = 0.05)
  expect_equal(rowSums(excluded$component_probability_mean),
               rep(1, nrow(fixture$annotation)), tolerance = 1e-12)
})

test_that("Phase-4A moderate-J C++ hierarchy is finite and exchange coherent", {
  fixture <- .sbs3_moderate_fixture(
    observations = 180L, annotation_count = 12L, seed = 20270840L
  )
  chains <- list(
    .sbs4a_cpp_chain(fixture, 20270841L, 3000L, 500L, rep(0L, 12L),
                     a_pi = 1, b_pi = 9),
    .sbs4a_cpp_chain(fixture, 20270842L, 3000L, 500L, rep(1L, 12L),
                     a_pi = 1, b_pi = 9)
  )
  summary <- .sbs4a_cpp_summary(chains)
  expect_true(all(is.finite(c(
    summary$annotation_pip, summary$alpha_mean, summary$pi_a_mean,
    summary$tau2_mean, summary$q_mean, summary$component_probability_mean
  ))))
  expect_equal(rowSums(summary$component_probability_mean),
               rep(1, nrow(fixture$annotation)), tolerance = 1e-12)
  expect_lte(max(abs(summary$chain_pip[1L, ] - summary$chain_pip[2L, ])),
             0.20)
  expect_gt(mean(summary$annotation_pip[1:3]),
            mean(summary$annotation_pip[5:12]))
})

test_that("Phase-4A C++ empty sticks use the proper intercept prior", {
  fixture <- .sbs2_fixture(observations = 80L, seed = 20271090L)
  fixture$eligible <- list(seq_len(80L), integer(), integer())
  fixture$outcome <- list(rep(0L, 80L), integer(), integer())
  chain <- .sbs4a_cpp_chain(
    fixture, 20271091L, 5000L, 500L, rep(0L, 3L), rep(0L, 3L)
  )
  later <- chain$alpha_draws[, 1L, 2:3, drop = FALSE]
  later_matrix <- matrix(later, nrow = dim(later)[1L], ncol = 2L)
  expect_lte(max(abs(colMeans(later_matrix) -
                       fixture$intercept_prior$mean[2:3])), 0.06)
  expect_lte(max(abs(apply(later_matrix, 2L, stats::var) -
                       fixture$intercept_prior$variance[2:3])), 0.08)
  expect_true(all(is.finite(c(chain$q_mean,
                              chain$component_probability_mean))))
  expect_equal(rowSums(chain$component_probability_mean),
               rep(1, nrow(fixture$annotation)), tolerance = 1e-12)
})
