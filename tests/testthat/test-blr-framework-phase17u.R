test_that("Phase 17U matches pinned posterior 1.6.1 diagnostic values", {
  fixtures <- phase17t_fixtures()
  well <- phase17u_fixture_metrics(fixtures$well_mixed)
  expect_equal(well$rhat, 0.974679434480896, tolerance = 1e-14)
  expect_equal(well$ess_bulk, 54.458829511151, tolerance = 1e-12)
  expect_equal(well$ess_q05, 197.402597402597, tolerance = 1e-11)
  expect_equal(well$ess_q95, 197.402597402597, tolerance = 1e-11)
  expect_equal(well$ess_tail, 197.402597402597, tolerance = 1e-11)
  expect_equal(well$ess_mean, 43.8040345821326, tolerance = 1e-12)
  expect_equal(well$mcse_mean, 0.183995239970637, tolerance = 1e-13)
  expect_equal(well$rhat, phase17t_rhat(fixtures$well_mixed),
               tolerance = 1e-12)
  expect_equal(well$ess_bulk, phase17t_ess_bulk(fixtures$well_mixed),
               tolerance = 1e-8)
  expect_equal(well$ess_tail, phase17t_ess_tail(fixtures$well_mixed),
               tolerance = 1e-8)
  expect_equal(well$ess_mean, phase17t_ess_mean(fixtures$well_mixed),
               tolerance = 1e-8)

  poor <- phase17u_fixture_metrics(fixtures$poor_tail)
  expect_equal(poor$ess_tail, 21.6910453085979, tolerance = 1e-12)
  expect_lt(poor$ess_tail, well$ess_tail)
  expect_gt(phase17u_fixture_metrics(fixtures$shifted)$rhat, 1.01)
  expect_gt(phase17u_fixture_metrics(fixtures$scales)$rhat_folded, 1.01)
  expect_gt(phase17u_fixture_metrics(fixtures$drift)$rhat, 1.01)
})

test_that("Phase 17U is numerically compatible with posterior when present", {
  skip_if_not_installed("posterior")
  for (chains in phase17t_fixtures()[c(
    "well_mixed", "shifted", "scales", "positive_ar", "poor_tail",
    "tied", "binary")]) {
    got <- phase17u_fixture_metrics(chains)
    x <- t(chains)
    expect_equal(got$rhat, posterior::rhat(x), tolerance = 1e-12)
    expect_equal(got$ess_bulk,
                 suppressWarnings(posterior::ess_bulk(x)),
                 tolerance = 1e-8)
    expect_equal(got$ess_tail,
                 suppressWarnings(posterior::ess_tail(x)),
                 tolerance = 1e-8)
    expect_equal(got$ess_mean,
                 suppressWarnings(posterior::ess_mean(x)),
                 tolerance = 1e-8)
    expect_equal(got$mcse_mean,
                 suppressWarnings(posterior::mcse_mean(x)),
                 tolerance = 1e-8)
  }
})

test_that("Phase 17U separates metric availability and statuses", {
  fixtures <- phase17t_fixtures()
  one <- phase17u_fixture_metrics(fixtures$one_chain)
  expect_identical(one$status, "unavailable_single_chain")
  expect_false(any(unlist(one[c(
    "rhat_available", "ess_bulk_available", "ess_tail_available",
    "ess_mean_available", "mcse_mean_available")])))

  four <- phase17u_fixture_metrics(fixtures$well_mixed[, 1:4])
  five <- phase17u_fixture_metrics(fixtures$well_mixed[, 1:5])
  six <- phase17u_fixture_metrics(fixtures$well_mixed[, 1:6])
  expect_identical(four$status, "computed_partial")
  expect_true(four$rhat_available)
  expect_false(four$ess_bulk_available)
  expect_identical(five$status, "computed_partial")
  expect_true(five$rhat_available)
  expect_false(five$ess_mean_available)
  expect_true(six$rhat_available)
  expect_true(six$ess_bulk_available)
  expect_true(six$mcse_mean_available)

  expect_identical(phase17u_fixture_metrics(
    fixtures$constant)$status, "constant")
  mismatch <- phase17u_fixture_metrics(fixtures$one_constant)
  expect_identical(mismatch$status, "constant_chain_mismatch")
  expect_true(mismatch$rhat_flag)
  expect_identical(phase17u_fixture_metrics(
    fixtures$nonfinite)$status, "nonfinite")
  not_updated <- sblr:::.blr_convergence_scalar(
    t(fixtures$well_mixed), updated = FALSE)
  expect_identical(not_updated$status, "not_updated")
  expect_true(all(is.na(unlist(not_updated[c(
    "rhat", "ess_bulk", "ess_tail", "ess_mean", "mcse_mean")]))))
})

test_that("Phase 17U rank, split, tail, and ESS stabilization are explicit", {
  odd <- matrix(seq_len(20), 5, 4)
  split <- sblr:::.blr_convergence_split_chains(odd)
  expect_identical(dim(split), c(2L, 8L))
  expect_false(any(split == odd[3L, 1L]))
  tied <- sblr:::.blr_convergence_rank_normalize(
    matrix(c(1, 1, 2, 3), 2))
  expected <- qnorm((c(1.5, 1.5, 3, 4) - 3 / 8) / (4 + 1 / 4))
  expect_equal(as.vector(tied), expected, tolerance = 1e-15)

  antithetic <- phase17u_fixture_metrics(
    phase17t_fixtures()$well_mixed)
  expect_gt(antithetic$ess_tail, length(
    phase17t_fixtures()$well_mixed))
  expect_lte(antithetic$ess_tail,
             length(phase17t_fixtures()$well_mixed) *
               log10(length(phase17t_fixtures()$well_mixed)))
  expect_true(is.logical(antithetic$ess_stability_bound_applied))
})

test_that("Phase 17U trace route preserves raw and extracts post-burn draws", {
  case <- phase17o_case(nt = 2L, updates = TRUE)
  case$nit <- 8L
  case$nburn <- 3L
  on.exit(phase17o_cleanup(case), add = TRUE)
  ordinary <- phase17r_call(case, nchains = 2L, ncores = 1L,
                            keep_chains = TRUE)
  native <- phase17u_native_call(case, nchains = 2L, ncores = 1L,
                                 keep_chains = TRUE)
  expect_identical(phase17r_without_timing(native$raw),
                   phase17r_without_timing(ordinary))
  bundle <- native$trace_bundle
  expect_silent(sblr:::.blr_validate_convergence_trace_bundle(
    bundle, nt = 2L, updateB = TRUE, updateE = TRUE))
  expect_identical(dim(bundle$values), c(8L, 2L, 10L))
  expect_identical(as.character(bundle$quantities$group),
                   rep(c("vbs", "vgs", "ves", "vle", "vld"), each = 2L))
  expect_identical(as.integer(bundle$quantities$trait_index),
                   rep(1:2, 5L))
  source_trace <- ordinary$chains$chain1$trace$vbs
  expect_identical(bundle$values[, 1L, 1L],
                   source_trace[4:11, 1L])
  expect_identical(bundle$values[, 2L, 10L],
                   ordinary$chains$chain2$trace$vld[4:11, 2L])
  expect_equal(bundle$values[, , 9:10],
               bundle$values[, , 3:4] - bundle$values[, , 7:8],
               tolerance = 1e-12)
})

test_that("Phase 17U diagnostics are independent of compact chains and workers", {
  case <- phase17o_case(nt = 2L, updates = TRUE)
  case$nit <- 8L
  on.exit(phase17o_cleanup(case), add = TRUE)
  dropped <- phase17u_native_call(case, nchains = 2L, ncores = 1L,
                                  keep_chains = FALSE)
  retained <- phase17u_native_call(case, nchains = 2L, ncores = 1L,
                                   keep_chains = TRUE)
  expect_identical(dropped$trace_bundle, retained$trace_bundle)
  d1 <- phase17u_diagnose(dropped, colnames(case$Y), TRUE, TRUE)
  d2 <- phase17u_diagnose(retained, colnames(case$Y), TRUE, TRUE)
  expect_identical(d1$raw$diagnostics$convergence,
                   d2$raw$diagnostics$convergence)
  expect_null(d1$convergence_traces)
  kept <- phase17u_diagnose(
    dropped, colnames(case$Y), TRUE, TRUE, keep_traces = TRUE)
  expect_identical(dim(kept$convergence_traces$values), c(8L, 2L, 10L))
  expect_identical(dimnames(kept$convergence_traces$values)[[2L]],
                   c("chain1", "chain2"))

  if (isTRUE(dropped$raw$diagnostics$mt_bed$openmp_available)) {
    parallel <- phase17u_native_call(case, nchains = 2L, ncores = 2L)
    expect_identical(dropped$trace_bundle, parallel$trace_bundle)
    expect_identical(
      phase17r_without_timing(dropped$raw),
      phase17r_without_timing(parallel$raw))
  }
})

test_that("Phase 17U update flags and nthin do not change Tier 1 ownership", {
  case <- phase17o_case(nt = 1L, updates = FALSE)
  case$nit <- 8L
  on.exit(phase17o_cleanup(case), add = TRUE)
  first <- phase17u_native_call(case, nchains = 2L)
  case$nthin <- 2L
  second <- phase17u_native_call(case, nchains = 2L)
  expect_identical(first$trace_bundle, second$trace_bundle)
  expect_identical(first$trace_bundle$quantities$updated,
                   c(FALSE, TRUE, FALSE, TRUE, TRUE))
  diagnosed <- phase17u_diagnose(first, "T1", FALSE, FALSE)
  tab <- diagnosed$raw$diagnostics$convergence$summary
  expect_identical(tab$status[c(1L, 3L)],
                   c("not_updated", "not_updated"))
  expect_true(all(tab$status[c(2L, 4L, 5L)] %in% c(
    "computed", "computed_fewer_than_four_chains",
    "computed_partial", "constant", "constant_chain_mismatch")))
  expect_true(all(is.na(tab$rhat[c(1L, 3L)])))
})

test_that("Phase 17U actual scope covers chain, trait, and covariance modes", {
  specifications <- list(
    list(nt = 1L, mode = "diagonal", nchains = 1L,
         updateB = TRUE, updateE = TRUE),
    list(nt = 2L, mode = "full", nchains = 2L,
         updateB = TRUE, updateE = FALSE),
    list(nt = 3L, mode = "diagonal", nchains = 4L,
         updateB = FALSE, updateE = TRUE)
  )
  for (spec in specifications) {
    case <- phase17o_case(
      nt = spec$nt, residual_covariance = spec$mode, updates = TRUE)
    case$nit <- 8L
    case$updateB <- spec$updateB
    case$updateE <- spec$updateE
    on.exit(phase17o_cleanup(case), add = TRUE)
    native <- phase17u_native_call(
      case, nchains = spec$nchains, ncores = 1L)
    diagnosed <- phase17u_diagnose(
      native, colnames(case$Y), spec$updateB, spec$updateE)
    convergence <- diagnosed$raw$diagnostics$convergence
    expect_identical(nrow(convergence$summary), 5L * spec$nt)
    expect_identical(convergence$nchains, spec$nchains)
    expect_identical(convergence$postburn_draws_per_chain, 8L)
    expect_true(all(convergence$summary$status[
      convergence$summary$group == "vbs"] ==
        if (spec$updateB) convergence$summary$status[
          convergence$summary$group == "vbs"] else "not_updated"))
    expect_true(all(convergence$summary$status[
      convergence$summary$group == "ves"] ==
        if (spec$updateE) convergence$summary$status[
          convergence$summary$group == "ves"] else "not_updated"))
  }
})

test_that("Phase 17U diagnostics are seed reproducible and fit local", {
  case <- phase17o_case(nt = 2L, residual_covariance = "full",
                        updates = TRUE)
  case$nit <- 8L
  on.exit(phase17o_cleanup(case), add = TRUE)
  args <- phase17u_args(
    case, nchains = 2L, ncores = 1L,
    chain_seeds = c(-1L, 19L), keep_chains = FALSE)
  first <- do.call(sblr:::mtblr_bed_convergence_trace_internal, args)
  second <- do.call(sblr:::mtblr_bed_convergence_trace_internal, args)
  expect_identical(first$trace_bundle, second$trace_bundle)
  first_diagnostic <- phase17u_diagnose(
    first, colnames(case$Y), TRUE, TRUE)
  second_diagnostic <- phase17u_diagnose(
    second, colnames(case$Y), TRUE, TRUE)
  expect_identical(first_diagnostic$raw$diagnostics$convergence,
                   second_diagnostic$raw$diagnostics$convergence)

  dense_case <- phase17o_case()
  on.exit(phase17o_cleanup(dense_case), add = TRUE)
  invisible(do.call(sblr:::mtblr, phase17o_dense_args(dense_case)))
  after_summary <- do.call(
    sblr:::mtblr_bed_convergence_trace_internal, args)
  expect_identical(first$trace_bundle, after_summary$trace_bundle)

  scalar <- phase17p_case()
  on.exit(phase17p_cleanup(scalar), add = TRUE)
  invisible(do.call(stblr_bed, c(
    list(y = scalar$Y[, 1L], Glist = scalar$Glist, rows = scalar$rows),
    list(nit = 1L, nburn = 0L, nthin = 1L, ncores = 1L,
         updateB = FALSE, updateE = FALSE, updatePi = FALSE))))
  after_scalar <- do.call(
    sblr:::mtblr_bed_convergence_trace_internal, args)
  expect_identical(first$trace_bundle, after_scalar$trace_bundle)
})

test_that("Phase 17U internal diagnostics are fresh-process reproducible", {
  skip_if_not_installed("callr")
  case <- phase17o_case(nt = 2L, updates = TRUE)
  case$nit <- 8L
  on.exit(phase17o_cleanup(case), add = TRUE)
  args <- phase17u_args(case, nchains = 2L, ncores = 1L)
  expected <- do.call(sblr:::mtblr_bed_convergence_trace_internal, args)
  root <- blr_repo_path()
  fresh <- callr::r(function(args, root) {
    if (!is.null(root)) {
      pkgload::load_all(root, compile = FALSE, quiet = TRUE)
    } else {
      library(sblr)
    }
    do.call(getFromNamespace(
      "mtblr_bed_convergence_trace_internal", "sblr"), args)
  }, list(args = args, root = root))
  expect_identical(expected$trace_bundle, fresh$trace_bundle)
  expect_identical(phase17r_without_timing(expected$raw),
                   phase17r_without_timing(fresh$raw))
})

test_that("Phase 17U validators and memory accounting are failure closed", {
  fixture <- phase17u_bundle_from_chains(
    phase17t_fixtures()$well_mixed)
  result <- sblr:::.blr_convergence_tier1(fixture, "T1")
  expect_silent(sblr:::.blr_validate_convergence_result(result))
  bad <- fixture
  bad$values <- bad$values[-1]
  expect_error(sblr:::.blr_validate_convergence_trace_bundle(bad),
               "dimensions")
  bad <- fixture
  bad$values[1L, 1L, 1L] <- Inf
  expect_error(sblr:::.blr_validate_convergence_trace_bundle(bad),
               "quantities")
  memory <- sblr:::.blr_convergence_memory_estimate(4, 1000, 5, TRUE)
  expect_identical(memory$trace_capture_bytes, 8 * 4 * 1000 * 5 * 5)
  expect_identical(memory$retained_trace_bytes,
                   memory$trace_capture_bytes)
  expect_false(memory$measured_rss)
  expect_false(memory$measured_peak_rss)
  expect_lt(sblr:::.blr_convergence_memory_estimate(
    4, 1000, 5, FALSE)$estimated_total_bytes,
    memory$estimated_total_bytes)
})

test_that("Phase 17U engine remains protected after public activation", {
  expect_true("convergence" %in% names(formals(mtblr_bed)))
  case <- phase17p_case()
  on.exit(phase17p_cleanup(case), add = TRUE)
  fit <- do.call(mtblr_bed, phase17p_public_args(case))
  expect_true("convergence" %in% names(fit))
  expect_true("convergence_traces" %in% names(fit))
  root <- blr_repo_path()
  skip_if(is.null(root), "source checkout unavailable")
  public <- paste(readLines(file.path(root, "R", "mtblr-bed.R"),
                           warn = FALSE), collapse = "\n")
  expect_true(grepl("mtblr_bed_convergence_trace_internal", public,
                    fixed = TRUE))
  expect_true(grepl("mtblr_bed_chains_internal", public, fixed = TRUE))
})
