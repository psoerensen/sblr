test_that("Phase 17T fixes post-burn, split, and rank contracts", {
  traces <- matrix(1:48, 4, 12)
  post <- phase17t_postburn(traces, 4)
  expect_identical(post, traces[, 5:12, drop = FALSE])
  changed <- traces; changed[, 1:4] <- -999
  expect_equal(phase17t_postburn(changed, 4), post)
  split <- phase17t_split_chains(matrix(1:36, 4, 9))
  expect_equal(dim(split), c(8, 4))
  expect_false(any(split == matrix(1:36, 4, 9)[, 5]))
  tied <- phase17t_rank_normalize(c(1, 1, 2, 3))
  expected <- qnorm((c(1.5, 1.5, 3, 4) - 3 / 8) / (4 + 1 / 4))
  expect_equal(tied, expected, tolerance = 1e-15)
})

test_that("Phase 17T R-hat contract detects location and scale failures", {
  x <- phase17t_fixtures()
  expect_lt(phase17t_rhat(x$well_mixed), 1.01)
  expect_gt(phase17t_rank_rhat(x$shifted), 1.01)
  expect_gt(phase17t_folded_rhat(x$scales), 1.01)
  expect_equal(phase17t_rhat(x$scales),
               max(phase17t_rank_rhat(x$scales),
                   phase17t_folded_rhat(x$scales)))
  expect_gt(phase17t_rhat(x$drift), 1.01)
  expect_true(all(is.finite(phase17t_rank_normalize(x$tied))))
  expect_true(all(is.finite(phase17t_rank_normalize(x$binary))))
})

test_that("Phase 17T ESS and MCSE oracles distinguish mixing behavior", {
  x <- phase17t_fixtures()
  expect_lt(phase17t_ess_bulk(x$positive_ar),
            phase17t_ess_bulk(x$well_mixed))
  expect_true(is.finite(phase17t_ess_bulk(x$negative_ar)))
  expect_lt(phase17t_ess_tail(x$poor_tail),
            phase17t_ess_tail(x$well_mixed))
  expect_true(is.finite(phase17t_ess_mean(x$well_mixed)))
  mcse <- phase17t_mcse_mean(x$well_mixed)
  expect_equal(unname(mcse["mcse_mean"]),
               unname(mcse["posterior_sd"] / sqrt(mcse["ess_mean"])))
  expect_equal(unname(mcse["mcse_mean_over_sd"]),
               unname(1 / sqrt(mcse["ess_mean"])))
})

test_that("Phase 17T statuses are failure-closed", {
  x <- phase17t_fixtures()
  expect_identical(phase17t_status(x$four_chains, FALSE), "not_updated")
  expect_identical(phase17t_status(x$constant), "constant")
  expect_identical(phase17t_status(x$nonfinite), "nonfinite")
  expect_identical(phase17t_status(x$one_chain), "unavailable_single_chain")
  expect_identical(phase17t_status(x$short), "insufficient_draws")
  expect_identical(phase17t_status(x$two_chains),
                   "computed_fewer_than_four_chains")
  expect_identical(phase17t_status(x$four_chains), "computed")
  expect_identical(phase17t_status(x$one_constant), "constant_chain_mismatch")
  expect_gt(phase17t_rhat(x$one_constant), 1.01)
  expect_true(all(is.na(phase17t_mcse_mean(x$constant)[c(1, 3, 4)])))
})

test_that("Phase 17T thresholds and overview are advisory and aggregated", {
  flags <- phase17t_flags(1.02, 350, 399, 0.06, 4)
  expect_true(all(flags))
  expect_false(any(phase17t_flags(1.01, 400, 400, 0.05, 4)))
  tab <- data.frame(quantity = c("B_diag[T1]", "G_diag[T1]"),
    status = c("computed", "constant"), rhat = c(1.02, NA),
    ess_bulk = c(300, NA), ess_tail = c(250, NA),
    mcse_mean_over_sd = c(.06, NA), nchains = 4,
    rhat_flag = c(TRUE, FALSE), ess_bulk_flag = c(TRUE, FALSE),
    ess_tail_flag = c(TRUE, FALSE), mcse_flag = c(TRUE, FALSE))
  overview <- phase17t_overview(tab)
  expect_identical(overview$overall_status, "warning")
  expect_equal(overview$n_flagged, 1)
  expect_false(overview$fewer_than_four_chains)
})

test_that("Phase 17T memory formulas count only diagnostic storage", {
  got <- phase17t_memory(4, 1000, 5, 16, 10)
  expect_equal(unname(got["tier1_trace_bytes"]), 8 * 4 * 1000 * 3 * 5)
  expect_equal(unname(got["covariance_trace_bytes"]), 8 * 4 * 1000 * 3 * 15)
  expect_equal(unname(got["probability_trace_bytes"]), 8 * 4 * 1000 * 2)
  expect_equal(unname(got["full_pi_trace_bytes"]), 8 * 4 * 1000 * 16)
  expect_equal(unname(got["selected_b_trace_bytes"]), 8 * 4 * 1000 * 10 * 5)
  expect_equal(unname(got["selected_d_trace_bytes"]), 4 * 4 * 1000 * 10 * 5)
  expect_equal(unname(got["per_quantity_workspace_bytes"]), 8 * 4 * 1000)
})

test_that("Phase 17T tiers and retention are explicitly separate", {
  scope <- phase17t_scope_contract()
  expect_identical(scope$tier1, c("B_diag", "G_diag", "E_diag"))
  expect_setequal(scope$tier2_requires_new_traces,
                  c("B_lower", "G_lower", "E_lower", "pi_null", "pi_active"))
  expect_identical(scope$tier3_opt_in,
                   c("selected_marker_b", "selected_marker_d"))
  expect_false(scope$all_markers_default)
  expect_false(scope$full_pi_default)
  expect_false(scope$keep_chains_required)
  expect_false(scope$diagnostic_thinning)
  expect_identical(scope$trace_retention, "independent_bundle")
})

test_that("Phase 17T contract remains protected after public activation", {
  expect_true("convergence" %in% names(formals(mtblr_bed)))
  fit_formals <- names(formals(mtblr_bed))
  expect_identical(tail(fit_formals, 8),
    c("nchains", "ncores", "chain_seeds", "keep_chains",
      "convergence", "convergence_control", "memory_warning_gb", "verbose"))
  root <- blr_repo_path()
  skip_if(is.null(root), "source checkout unavailable")
  source <- paste(readLines(file.path(root, "R", "mtblr-bed.R"), warn = FALSE),
                  collapse = "\n")
  native <- paste(readLines(file.path(root, "src", "mtblr.cpp"), warn = FALSE),
                  collapse = "\n")
  expect_true(grepl("diagnostics\\$convergence", source))
  expect_true(grepl("MtBedConvergenceTraceBundle", native, fixed = TRUE))
  expect_false(grepl("\\.mtblr_convergence_ess <-", source))
  expect_true(grepl("iterationwise_chain_mean", source, fixed = TRUE))
})
