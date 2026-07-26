test_that("Phase 17V exposes exact convergence formals and controls", {
  f <- formals(mtblr_bed)
  expect_identical(names(f)[39:42], c(
    "convergence", "convergence_control", "memory_warning_gb", "verbose"))
  expect_identical(f$convergence, quote(c("auto", "none", "core")))
  expect_null(f$convergence_control)
  defaults <- sblr:::.mtblr_bed_convergence_controls("auto", NULL, 2L)
  expect_identical(defaults[c("warn", "rhat_threshold",
    "ess_per_chain_threshold", "mcse_mean_over_sd_threshold",
    "keep_traces")], list(warn = TRUE, rhat_threshold = 1.01,
      ess_per_chain_threshold = 100,
      mcse_mean_over_sd_threshold = 0.05, keep_traces = FALSE))
  expect_true(defaults$trace_route_required)
  expect_false(sblr:::.mtblr_bed_convergence_controls(
    "auto", NULL, 1L)$trace_route_required)
  expect_true(sblr:::.mtblr_bed_convergence_controls(
    "core", list(keep_traces = TRUE), 1L)$trace_route_required)
})

test_that("Phase 17V convergence controls fail closed", {
  resolver <- sblr:::.mtblr_bed_convergence_controls
  invalid <- list(
    1, data.frame(warn = TRUE), list(TRUE),
    structure(list(TRUE, FALSE), names = c("warn", "warn")),
    list(unknown = TRUE), list(warn = NA), list(warn = 1),
    list(warn = c(TRUE, FALSE)), list(keep_traces = NA),
    list(keep_traces = 1), list(keep_traces = c(TRUE, FALSE)),
    list(rhat_threshold = Inf), list(rhat_threshold = 0),
    list(ess_per_chain_threshold = -1),
    list(mcse_mean_over_sd_threshold = c(.05, .1)))
  for (value in invalid) expect_error(resolver("auto", value, 2L))
  expect_error(resolver("none", list(keep_traces = TRUE), 2L),
               "cannot retain")
  expect_silent(resolver("core", list(warn = FALSE), 2L))
  expect_silent(resolver("auto", list(rhat_threshold = 1.02), 2L))
})

test_that("Phase 17V constructors validate disabled and quiet states", {
  control <- sblr:::.mtblr_convergence_control()
  none <- sblr:::.mtblr_convergence_not_requested(
    c("T1", "T2"), TRUE, TRUE, 1L, 8L, control)
  expect_false(none$requested)
  expect_false(none$computed)
  expect_identical(none$scope, "none")
  expect_identical(none$overall_status, "not_requested")
  expect_identical(nrow(none$summary), 0L)
  unavailable <- sblr:::.mtblr_convergence_unavailable(
    c("T1", "T2"), FALSE, FALSE, 1L, 8L, control)
  expect_true(unavailable$requested)
  expect_false(unavailable$computed)
  expect_identical(nrow(unavailable$summary), 6L)
  expect_true(all(unavailable$summary$status[
    unavailable$summary$group %in% c("B_diag", "E_diag")] == "not_updated"))
  expect_true(all(unavailable$summary$status[
    unavailable$summary$group == "G_diag"] == "unavailable_single_chain"))
  expect_true(all(!unavailable$summary$rhat_available))
})

test_that("Phase 17V public routing and Tier 1 output are exact", {
  case <- phase17p_case()
  on.exit(phase17p_cleanup(case), add = TRUE)
  none_args <- phase17v_public_args(case, convergence = "none")
  none <- do.call(mtblr_bed, none_args)
  expect_identical(none$input$convergence_trace_route,
                   "mtblr_bed_chains_internal")
  expect_identical(none$convergence$overall_status, "not_requested")
  expect_null(none$convergence_traces)
  auto1 <- do.call(mtblr_bed, phase17v_public_args(case))
  expect_identical(auto1$input$convergence_trace_route,
                   "mtblr_bed_chains_internal")
  expect_identical(auto1$convergence$overall_status, "unavailable")
  expect_false(auto1$input$convergence_warning_emitted)
  args <- phase17v_public_args(case, nchains = 2L, convergence = "core",
                              keep_traces = TRUE)
  fit <- do.call(mtblr_bed, args)
  expect_identical(fit$input$convergence_trace_route,
                   "mtblr_bed_convergence_trace_internal")
  expect_identical(nrow(fit$convergence$summary), 3L * ncol(case$Y))
  expect_identical(dim(fit$convergence_traces$values),
                   c(8L, 2L, 3L * ncol(case$Y)))
  expected <- phase17v_internal_diagnostic(
    args, colnames(case$Y), keep_traces = TRUE)
  expect_identical(fit$convergence$summary,
                   expected$raw$diagnostics$convergence$summary)
  expect_identical(fit$convergence$overview,
                   expected$raw$diagnostics$convergence$overview)
  expect_identical(fit$convergence_traces, expected$convergence_traces)
  raw <- phase17s_internal(none_args)
  expect_equal(phase17s_fit_numerics(none), phase17s_raw_numerics(raw),
               tolerance = 1e-12)
})

test_that("Phase 17V public warnings follow mode and suppression policy", {
  case <- phase17p_case()
  on.exit(phase17p_cleanup(case), add = TRUE)
  expect_silent(do.call(mtblr_bed, phase17v_public_args(case)))
  expect_silent(do.call(mtblr_bed, phase17v_public_args(
    case, convergence = "none", warn = TRUE)))
  expect_warning(do.call(mtblr_bed, phase17v_public_args(
    case, convergence = "core", warn = TRUE)),
    "diagnostics are unavailable")
  expect_silent(do.call(mtblr_bed, phase17v_public_args(
    case, convergence = "core", warn = FALSE)))
  memory_args <- phase17v_public_args(case)
  memory_args$memory_warning_gb <- 1e-12
  expect_warning(do.call(mtblr_bed, memory_args),
    "convergence=auto.*convergence trace capture=FALSE")
})

test_that("Phase 17V warning policy is aggregated and suppressible", {
  flagged <- phase17v_warning_fixture(flagged = TRUE)
  message <- sblr:::.mtblr_convergence_warning_messages(flagged, "auto")
  expect_length(message, 1L)
  expect_match(message, "advisory requires review")
  expect_match(message, "max R-hat")
  expect_match(message, "fit\\$convergence")
  partial <- phase17v_warning_fixture("computed_partial", FALSE, 2L, 5L)
  expect_length(sblr:::.mtblr_convergence_warning_messages(
    partial, "core"), 1L)
  unavailable <- sblr:::.mtblr_convergence_unavailable(
    "T1", TRUE, TRUE, 1L, 8L, sblr:::.mtblr_convergence_control())
  expect_match(sblr:::.mtblr_convergence_warning_messages(
    unavailable, "core"), "unavailable")
  expect_length(sblr:::.mtblr_convergence_warning_messages(
    unavailable, "auto"), 0L)
})

test_that("Phase 17V memory accounting separates diagnostic storage", {
  controls <- sblr:::.mtblr_bed_convergence_controls("core", NULL, 4L)
  memory <- sblr:::.mtblr_bed_convergence_memory(
    "core", controls, 4L, 100L, 2L)
  expect_identical(memory$trace_capture_bytes, 8 * 4 * 100 * 3 * 2)
  expect_gt(memory$maximum_workspace_bytes, 0)
  kept_controls <- sblr:::.mtblr_bed_convergence_controls(
    "core", list(keep_traces = TRUE), 4L)
  kept <- sblr:::.mtblr_bed_convergence_memory(
    "core", kept_controls, 4L, 100L, 2L)
  expect_identical(kept$retained_trace_bytes, kept$trace_capture_bytes)
  none <- sblr:::.mtblr_bed_convergence_memory(
    "none", sblr:::.mtblr_bed_convergence_controls("none", NULL, 4L),
    4L, 100L, 2L)
  expect_identical(none$estimated_total_bytes, 0)
  quiet <- sblr:::.mtblr_bed_convergence_memory(
    "auto", sblr:::.mtblr_bed_convergence_controls("auto", NULL, 1L),
    1L, 100L, 2L)
  expect_identical(quiet$trace_capture_bytes, 0)
  expect_gt(quiet$summary_output_bytes, 0)
})

test_that("Phase 17V protects native and unrelated public surfaces", {
  expect_false("convergence" %in% names(formals(mtblr_csr)))
  expect_false("convergence" %in% names(formals(mtblr_block_eigen)))
  expect_false("convergence" %in% names(formals(stblr_bed)))
  expect_false("convergence" %in% names(formals(sblr)))
  root <- blr_repo_path()
  skip_if(is.null(root), "source checkout unavailable")
  protected <- c(
    "src/blr_mt_bed_convergence_types.h",
    "src/blr_mt_bed_convergence_trace_impl.h",
    "src/blr_mt_bed_core_impl.h", "src/blr_mt_bed_chains_types.h",
    "src/blr_mt_bed_chains_execution_impl.h",
    "src/blr_mt_bed_chains_aggregate_impl.h", "src/mtblr.cpp",
    "R/RcppExports.R", "src/RcppExports.cpp", "DESCRIPTION", "NAMESPACE")
  expect_identical(system2("git", c("diff", "--name-only", "--", protected),
                           stdout = TRUE), character())
})
