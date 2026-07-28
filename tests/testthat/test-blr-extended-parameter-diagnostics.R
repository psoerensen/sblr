test_that("extended controls fail closed and memory guards do not truncate", {
  expect_error(.blr_convergence_controls(
    "none", list(selected_markers = 1L), 2L), "cannot be used")
  expect_error(.blr_convergence_controls(
    "core", list(extended_groups = "probability"), 2L),
    "requires convergence")
  expect_error(.blr_convergence_controls(
    "extended", list(selected_markers = "all"), 2L), "shortcuts")
  memory <- .blr_extended_trace_memory(4L, 1000L, 100000L, 0L, TRUE)
  controls <- .blr_convergence_controls(
    "extended", list(max_trace_gb = 1e-6, warn = FALSE), 4L)
  expect_error(.blr_enforce_trace_guard(
    memory, controls, c(chains = 4L, draws = 1000L, covariance = 100000L), 0L),
    "allow_large_traces")
})

test_that("MT BayesRC annotation diagnostics retain column-major identities", {
  x <- .mt_bayesrc_fixture(); on.exit(.mt_bayesrc_cleanup(x), add = TRUE)
  common <- .mt_bayesrc_common(); common$annotations <- x$annotations
  common$updateAlpha <- TRUE; common$nchains <- 2L; common$ncores <- 1L
  common$convergence <- "extended"
  common$convergence_control <- list(
    warn = FALSE, extended_groups = "annotations", keep_traces = TRUE)
  fit <- do.call(mtblr_csr, c(list(
    stats = x$stats, ld_prefix = x$prefix, ld_metadata = x$ld_metadata),
    common))
  annotations <- fit$convergence$summary[
    fit$convergence$summary$tier == 2L, , drop = FALSE]
  expected <- as.vector(outer(
    fit$model_parameters$annotations$processed_annotation_names,
    fit$model_parameters$annotations$stick_names,
    function(annotation, stick) paste0("alpha[", annotation, ",", stick, "]")))
  expect_identical(annotations$quantity[seq_along(expected)], expected)
  expect_true(any(annotations$is_intercept))
  expect_true(all(grepl("^sigmaSqAlpha", tail(
    annotations$quantity, length(fit$model_parameters$annotations$stick_names)))))
})

test_that("sampled ST selection-S uses chain-private unthinned states", {
  chains <- list(T1 = list(
    chain1 = list(selection_s = seq(-1, 0, length.out = 10)),
    chain2 = list(selection_s = seq(-.8, .2, length.out = 10))))
  controls <- .blr_convergence_controls(
    "extended", list(warn = FALSE, extended_groups = "selection_s"), 2L)
  bundle <- .blr_st_extended_bundle(
    chains, "T1", "sbayesc", "csr", 8L, 2L,
    list(input = list(prior_kernel = "bayesc",
                      estimate_selection_s = TRUE)), controls)
  expect_identical(bundle$quantities$group, "selection_s")
  expect_equal(bundle$values[, 1L, 1L], chains$T1$chain1$selection_s[3:10])
  expect_equal(bundle$values[, 2L, 1L], chains$T1$chain2$selection_s[3:10])
})
