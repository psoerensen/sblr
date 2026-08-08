test_that("SBayesRC-S-EM Laplace model probabilities match quadrature", {
 fixture <- .sbs_em_ref_fixture()
 exact <- .sbs_em_ref_exact(
  fixture$A, fixture$responsibility, fixture$pi_a, fixture$tau2,
  fixture$intercept_prior, order = 17L
 )
 selection <- .sbayesrc_s_em_selection_update(
  fixture$A, fixture$responsibility, c(0L, 0L), fixture$alpha_start,
  fixture$pi_a, fixture$tau2, fixture$intercept_prior,
  sweeps = 6000L, burn = 1000L, seed = 7701L
 )
 expect_lt(max(abs(
  selection$annotation_pip_eb - exact$annotation_pip_eb
 )), 0.035)
 expect_true(all(selection$annotation_pip_eb >= 0 &
                 selection$annotation_pip_eb <= 1))
 expect_gt(selection$annotation_pip_eb[1L], selection$annotation_pip_eb[2L])
 expect_equal(selection$delta_map, c(1L, 0L))
 expect_equal(selection$alpha_map[3L, ], c(0, 0))
})

test_that("SBayesRC-S-EM selection is shared, proper, and start robust", {
 fixture <- .sbs_em_ref_fixture(empty_last_stick = TRUE)
 first <- .sbayesrc_s_em_selection_update(
  fixture$A, fixture$responsibility, c(0L, 0L), fixture$alpha_start,
  fixture$pi_a, fixture$tau2, fixture$intercept_prior,
  sweeps = 4000L, burn = 500L, seed = 7702L
 )
 second <- .sbayesrc_s_em_selection_update(
  fixture$A, fixture$responsibility, c(1L, 1L), fixture$alpha_start + 0.2,
  fixture$pi_a, fixture$tau2, fixture$intercept_prior,
  sweeps = 4000L, burn = 500L, seed = 7703L
 )
 expect_lt(max(abs(first$annotation_pip_eb - second$annotation_pip_eb)), 0.05)
 expect_true(all(is.finite(first$alpha_map)))
 expect_length(first$delta_map, ncol(fixture$A) - 1L)
 empty_model <- .sbayesrc_s_em_model_laplace(
  fixture$A, fixture$responsibility, c(0L, 0L), fixture$pi_a,
  fixture$tau2, fixture$intercept_prior, fixture$alpha_start
 )
 expect_equal(
  unname(empty_model$sticks[[2L]]$mode[1L]),
  unname(fixture$intercept_prior["mean", 2L]), tolerance = 1e-6
 )
})

test_that("SBayesRC-S-EM correlated annotations retain uncertainty", {
 fixture <- .sbs_em_ref_fixture(correlated = TRUE)
 exact <- .sbs_em_ref_exact(
  fixture$A, fixture$responsibility, fixture$pi_a, fixture$tau2,
  fixture$intercept_prior, order = 15L
 )
 expect_true(all(exact$annotation_pip_eb > 0.05))
 expect_lt(abs(diff(exact$annotation_pip_eb)), 0.40)
 expect_equal(sum(exact$model_probability), 1, tolerance = 1e-12)
})

test_that("SBayesRC-S-EM result semantics distinguish EB and joint PIPs", {
 fixture <- .sbs_em_ref_fixture()
 selection <- .sbayesrc_s_em_selection_update(
  fixture$A, fixture$responsibility, fixture$delta, fixture$alpha_start,
  fixture$pi_a, fixture$tau2, fixture$intercept_prior,
  sweeps = 1000L, burn = 200L, seed = 7704L
 )
 expect_identical(
  selection$target,
  "responsibility_conditioned_laplace_model_distribution"
 )
 expect_false("annotation_pip" %in% names(selection))
})

test_that("SBayesRC-S-EM CSR result has explicit MAP and EB semantics", {
 fixture <- .mcem_exact_fixture(correlated = TRUE)
 on.exit(.mcem_cleanup_csr(fixture$prefix), add = TRUE)
 fit <- .sbs_em_run_csr(fixture, delta_start = 0L, sweeps = 120L,
                        burn = 40L, outer = 3L)
 expect_identical(fit$mcem$method, "SBayesRC-S-EM")
 expect_identical(fit$mcem$algorithm, "MCEM-Laplace")
 expect_identical(fit$mcem$pi_A_mode, "fixed")
 expect_identical(fit$mcem$tau2_mode, "fixed")
 expect_length(fit$mcem$annotation_pip_eb, 1L)
 expect_false("annotation_pip" %in% names(fit$mcem))
 expect_equal(
  unname(fit$mcem$alpha_map[2L, ] * (1L - fit$mcem$delta_map)), c(0, 0)
 )
 expect_false(identical(
  fit$mcem$last_estep_responsibilities,
  fit$mcem$final_genomic_responsibilities
 ))
})

test_that("SBayesRC-S-EM block route retains learned B and E semantics", {
 fixture <- .mcem_block_fixture(marker_count = 8L, sample_count = 55L)
 on.exit(.mcem_cleanup_block_fixture(fixture), add = TRUE)
 fit <- .sbs_em_run_block(
  fixture, delta_start = 1L, updateB = TRUE, updateE = TRUE,
  sweeps = 100L, burn = 35L, outer = 3L
 )
 expect_identical(fit$mcem$backend, "block_eigen")
 expect_true(fit$mcem$genomic_hyperparameters$updateB)
 expect_true(fit$mcem$genomic_hyperparameters$updateE)
 expect_true(all(is.finite(c(
  fit$mcem$genomic_hyperparameters$B_final,
  fit$mcem$genomic_hyperparameters$E_final,
  fit$mcem$annotation_pip_eb
 ))))
 expect_equal(rowSums(fit$mcem$component_prior), rep(1, nrow(fixture$A)),
              tolerance = 1e-10)
})
