test_that("internal CSR SBayesRC-S backend returns coherent genomic state", {
  fixture <- .sbs4b_fixture(48L, 20270901L)
  fit <- .sbs4b_run(fixture, 20270902L, 180L, 40L)
  annotation <- fit$chains[[1L]][[1L]]$annotation
  empty_diagnostics <- annotation$empty_stick_diagnostics
  expect_identical(fit$meta$model, "sbayesrc_selection")
  expect_identical(fit$meta$backend, "csr_sbayesrc_selection_internal")
  expect_length(annotation$annotation_pip, 3L)
  expect_true(all(annotation$annotation_pip >= 0 &
                  annotation$annotation_pip <= 1))
  expect_equal(dim(annotation$alpha), c(4L, 3L))
  expect_equal(dim(annotation$alpha_mean_given_inclusion), c(3L, 3L))
  expect_true(all(is.finite(c(
    annotation$alpha, annotation$annotation_pi_A,
    annotation$annotation_tau2, annotation$annotation_included_mean
  ))))
  expect_equal(rowSums(fit$component$prob[[1L]]),
               rep(1, nrow(fixture$A)), tolerance = 1e-12)
  expect_equal(dim(fit$marker$bm), c(nrow(fixture$A), 1L))
  expect_equal(dim(fit$marker$dm), c(nrow(fixture$A), 1L))
})

test_that("genomic SBayesRC-S fixed-delta structural bridges hold", {
  fixture <- .sbs4b_fixture(60L, 20270910L)
  excluded <- .sbs4b_run(
    fixture, 20270911L, 300L, 75L, rep(0L, 3L), initial_delta = rep(0L, 3L)
  )
  excluded_annotation <- excluded$chains[[1L]][[1L]]$annotation
  expect_identical(as.numeric(excluded_annotation$annotation_pip), rep(0, 3L))
  expect_identical(as.numeric(excluded_annotation$alpha[-1L, ]), rep(0, 9L))
  expect_true(all(is.finite(excluded_annotation$alpha[1L, ])))

  included <- .sbs4b_run(
    fixture, 20270912L, 300L, 50L, rep(1L, 3L),
    initial_delta = rep(1L, 3L), update_hierarchy = FALSE
  )
  standard <- .sbs4b_run_standard(
    fixture, 20270912L, 300L, 50L, update_hierarchy = FALSE
  )
  selection_alpha <- included$chains[[1L]][[1L]]$annotation$alpha
  standard_alpha <- standard$chains[[1L]][[1L]]$annotation$alpha
  expect_identical(
    as.numeric(included$chains[[1L]][[1L]]$annotation$annotation_pip),
    rep(1, 3L)
  )
  expect_identical(selection_alpha, standard_alpha)
  expect_identical(included$marker, standard$marker)
  expect_identical(included$trace, standard$trace)
  expect_identical(included$component, standard$component)
})

test_that("disabled selection policy preserves standard SBayesRC RNG", {
  fixture <- .sbs4b_fixture(40L, 20270920L)
  first <- .sbs4b_run_standard(
    fixture, 20270921L, 80L, 20L, update_hierarchy = FALSE
  )
  second <- .sbs4b_run_standard(
    fixture, 20270921L, 80L, 20L, update_hierarchy = FALSE
  )
  expect_identical(first$marker, second$marker)
  expect_identical(first$trace, second$trace)
  expect_identical(first$annotation, second$annotation)
  expect_identical(first$component, second$component)
})

test_that("genomic SBayesRC-S supports an empty proper-intercept stick", {
  fixture <- .sbs4b_fixture(40L, 20270930L)
  fixture$comp_init <- list(rep(0L, 40L))
  fixture$b_init <- list(rep(0, 40L))
  fixture$r_init <- fixture$wy
  fit <- .sbs4b_run(fixture, 20270931L, 120L, 20L)
  annotation <- fit$chains[[1L]][[1L]]$annotation
  empty_diagnostics <- annotation$empty_stick_diagnostics
  expect_true(all(is.finite(c(
    annotation$alpha, annotation$annotation_pip,
    annotation$annotation_pi_A, annotation$annotation_tau2
  ))))
  expect_equal(rowSums(fit$component$prob[[1L]]),
               rep(1, nrow(fixture$A)), tolerance = 1e-12)
  expect_equal(dim(empty_diagnostics), c(3L, 4L))
  expect_true(all(empty_diagnostics[2:3, 1L] >= 1))
  expect_true(all(empty_diagnostics[2:3, 2L] >= 1))
  expect_true(all(empty_diagnostics[2:3, 3L] >= 1))
})

test_that("Phase-4C internal controls retain compact hierarchy traces", {
  fixture <- .sbs4b_fixture(36L, 20270940L)
  fit <- .sbs4b_run(
    fixture, 20270941L, 100L, 20L,
    update_pi_A = FALSE, update_tau2 = FALSE,
    hierarchy_sweeps = 2L, genomic_sweeps = 1L
  )
  trace <- fit$chains[[1L]][[1L]]$convergence_trace
  expect_equal(dim(trace$annotation_delta), c(100L, 3L))
  expect_length(trace$annotation_pi_A, 100L)
  expect_length(trace$annotation_included_count, 100L)
  expect_identical(as.numeric(trace$annotation_pi_A), rep(0.35, 100L))
  expect_equal(trace$sigmaSqAlpha, matrix(0.8, 100L, 3L))
  expect_identical(
    as.numeric(trace$annotation_included_count),
    as.numeric(rowSums(trace$annotation_delta))
  )
})
