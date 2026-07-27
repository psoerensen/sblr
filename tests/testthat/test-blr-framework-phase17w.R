test_that("Phase 17W fixes strict-lower covariance ordering and names", {
  expect_equal(nrow(phase17w_strict_lower(1)), 0)
  expect_identical(unname(phase17w_strict_lower(2)), matrix(c(2L, 1L), 1L))
  expect_identical(unname(phase17w_strict_lower(3)),
    rbind(c(2L, 1L), c(3L, 1L), c(3L, 2L)))
  expect_equal(nrow(phase17w_strict_lower(5)), 10)

  got <- phase17w_covariance_descriptors(c("zeta", "alpha", "mu"))
  expect_identical(got$group, rep(c("B_cov", "G_cov", "E_cov"), each = 3L))
  expect_identical(got$quantity[1:3],
    c("B[alpha,zeta]", "B[mu,zeta]", "B[mu,alpha]"))
  expect_identical(got$trait[1:3], c("alpha", "mu", "mu"))
  expect_identical(got$trait2[1:3], c("zeta", "zeta", "alpha"))
  expect_false(any(got$row == got$column))
})

test_that("Phase 17W covariance update and structural semantics are exact", {
  fixed <- phase17w_covariance_descriptors(c("T1", "T2"),
    updateB = FALSE, updateE = FALSE)
  expect_identical(fixed$status, c("not_updated", "eligible", "not_updated"))
  expect_identical(fixed$derived, c(FALSE, TRUE, FALSE))
  expect_identical(fixed$captured, c(FALSE, TRUE, FALSE))

  diagonal <- phase17w_covariance_descriptors(c("T1", "T2"),
    residual = "diagonal", updateE = TRUE)
  expect_identical(diagonal$status, c("eligible", "eligible", "structural_zero"))
  expect_false(diagonal$captured[3])
  expect_true(diagonal$structural[3])
})

test_that("Phase 17W probability mass and pattern selection deduplicate", {
  models <- rbind(c(0L, 0L), c(1L, 0L), c(1L, 1L))
  names <- c("null", "trait1", "both")
  expect_identical(phase17w_null_index(models), 1L)
  mass <- phase17w_probability_plan(models, names, "mass")
  expect_identical(mass$physical_model_indices, 1L)
  expect_identical(mass$primary_mass_quantity, "pi_active")
  expect_identical(mass$complement, "pi_null")
  expect_identical(mass$diagnostic_key, "pi_mass:null_active")

  expect_identical(phase17w_resolve_patterns(models, names, c("both", "null")),
                   c(3L, 1L))
  expect_identical(phase17w_resolve_patterns(models, names, c(2L, 1L)),
                   c(2L, 1L))
  expect_error(phase17w_resolve_patterns(models, names, c("null", "null")),
               "unique")
  expect_error(phase17w_resolve_patterns(models, names, "missing"), "unknown")
  expect_error(phase17w_probability_plan(models, names, "selected"), "requires")
  expect_error(phase17w_probability_plan(models, names, "all", 1L), "only")
  expect_identical(phase17w_probability_plan(models, names, "all")$
                     physical_model_indices, 1:3)
  expect_false(phase17w_probability_plan(models, names, "mass", updatePi = FALSE)$captured)
})

test_that("Phase 17W marker selection uses final public marker order", {
  ids <- c("m3", "m1", "m2")
  expect_identical(phase17w_resolve_markers(ids,
    marker_ids_request = c("m2", "m3")), c(3L, 1L))
  expect_identical(phase17w_resolve_markers(ids,
    marker_indices = c(2L, 1L)), c(2L, 1L))
  expect_error(phase17w_resolve_markers(ids, "unknown"), "unknown")
  expect_error(phase17w_resolve_markers(c("m1", "m1"), "m1"), "ambiguous")
  expect_error(phase17w_resolve_markers(ids, c("m1", "m1")), "invalid")
  expect_error(phase17w_resolve_markers(ids, marker_indices = c(1L, 1L)), "invalid")
  expect_error(phase17w_resolve_markers(ids, marker_indices = 4L), "invalid")
  expect_error(phase17w_resolve_markers(ids, "m1", 1L), "choose")
})

test_that("Phase 17W memory formulas reuse overlap and exclude shared data", {
  got <- phase17w_extended_memory(4, 1000, 5, residual = "full",
    probability = "mass", P = 16, K = 10, keep_traces = TRUE)
  qoff <- 10
  expect_equal(unname(got["B_cov_bytes"]), 8 * 4 * 1000 * qoff)
  expect_equal(unname(got["G_cov_bytes"]), 8 * 4 * 1000 * qoff)
  expect_equal(unname(got["E_cov_bytes"]), 8 * 4 * 1000 * qoff)
  expect_equal(unname(got["pi_mass_bytes"]), 8 * 4 * 1000)
  expect_equal(unname(got["probability_unique_bytes"]),
               unname(got["pi_mass_bytes"]))
  expect_equal(unname(got["selected_b_bytes"]), 8 * 4 * 1000 * 10 * 5)
  expect_equal(unname(got["selected_d_bytes"]), 4 * 4 * 1000 * 10 * 5)
  expect_equal(unname(got["retained_trace_bytes"]),
               unname(got["captured_trace_bytes"]))
  expect_false(any(grepl("genotype|phenotype", names(got))))
  expect_equal(unname(got["workspace_bytes"]), 8 * 4 * 1000 * 8)

  diagonal <- phase17w_extended_memory(2, 100, 3, updateB = FALSE,
    updateE = TRUE, residual = "diagonal", probability = "none",
    marker_quantities = character())
  expect_equal(unname(diagonal[c("B_cov_bytes", "E_cov_bytes")]), c(0, 0))
  expect_gt(unname(diagonal["G_cov_bytes"]), 0)
  all_pi <- phase17w_extended_memory(2, 100, 3, probability = "all", P = 16,
                                    K = 0, marker_quantities = character())
  expect_equal(unname(all_pi["probability_unique_bytes"]), 8 * 2 * 100 * 16)
})

test_that("Phase 17W quantity counts and large-request policy are failure safe", {
  expect_identical(phase17w_quantity_counts(3, 2, 4),
    c(tier1 = 9L, covariance = 9L, probability = 2L,
      marker_b = 12L, marker_d = 12L))
  expect_equal(phase17w_resolved_limit_gb(8), 2)
  expect_equal(phase17w_resolved_limit_gb(1), .25)
  below <- phase17w_large_request(.1 * 1024^3)
  above <- phase17w_large_request(3 * 1024^3)
  expect_identical(below$decision, "allow")
  expect_identical(above$decision, "require_explicit_override")
  expect_identical(phase17w_large_request(3 * 1024^3,
    allow_large_traces = TRUE)$decision, "allow")
})

test_that("Phase 17W future controls reject contradictions", {
  expect_identical(phase17w_validate_future_control()$probability, "mass")
  expect_error(phase17w_validate_future_control("core", list(probability = "mass")),
               "require")
  expect_error(phase17w_validate_future_control(extended = list(unknown = 1)),
               "unknown")
  expect_error(phase17w_validate_future_control(extended = structure(
    list("mass", "all"), names = c("probability", "probability"))), "invalid")
  expect_error(phase17w_validate_future_control(extended = list(
    probability = "selected")), "requires")
  expect_error(phase17w_validate_future_control(extended = list(
    probability = "mass", probability_models = "0_0")), "conflicts")
  expect_error(phase17w_validate_future_control(extended = list(
    marker_ids = "m1", marker_indices = 1L)), "conflict")
  expect_error(phase17w_validate_future_control(extended = list(
    marker_quantities = "b")), "require")
  expect_error(phase17w_validate_future_control(extended = list(
    marker_ids = "m1", marker_quantities = character())), "invalid")
  expect_error(phase17w_validate_future_control(extended = list(
    allow_large_traces = 1)), "logical")
  expect_error(phase17w_validate_future_control(extended = list(
    trace_limit_gb = 0)), "positive")
})

test_that("Phase 17W scalar applicability keeps binary and zero-inflated caveats", {
  binary <- rbind(c(rep(0, 9), 1, rep(0, 9), 1),
                  c(rep(0, 8), 1, 0, rep(0, 8), 1, 0),
                  c(rep(0, 7), 1, 0, 0, rep(0, 7), 1, 0, 0),
                  c(rep(0, 6), 1, 0, 0, 0, rep(0, 6), 1, 0, 0, 0))
  got <- sblr:::.blr_convergence_scalar(t(binary))
  expect_true(got$status %in% c("computed", "computed_partial"))
  expect_true(is.finite(got$ess_mean))
  expect_true(is.finite(got$mcse_mean))

  constant <- sblr:::.blr_convergence_scalar(matrix(0, 20, 4))
  expect_identical(constant$status, "constant")
  expect_true(is.na(constant$rhat))
  mismatch <- binary; mismatch[1, ] <- 0
  expect_identical(sblr:::.blr_convergence_scalar(t(mismatch))$status,
                   "constant_chain_mismatch")
  zero_inflated <- binary * matrix(rep(seq_len(20), each = 4), 4, 20)
  expect_true(sblr:::.blr_convergence_scalar(t(zero_inflated))$status %in%
                c("computed", "computed_partial"))
})

test_that("Phase 17W warnings group details without row-wise emission", {
  tab <- data.frame(group = c("B_cov", "B_cov", "marker_d"),
    status = c("computed", "computed", "constant_chain_mismatch"),
    rhat_flag = c(TRUE, FALSE, FALSE), ess_bulk_flag = FALSE,
    ess_tail_flag = c(FALSE, TRUE, FALSE), mcse_flag = FALSE)
  got <- phase17w_warning_groups(tab)
  expect_equal(as.integer(got), c(2, 1))
  expect_identical(names(got), c("B_cov", "marker_d"))
})

test_that("Phase 17W is contract-only and protects production", {
  root <- blr_repo_path()
  skip_if(is.na(root), "source architecture assertion requires a checkout")
  protected <- c("R/mtblr-bed.R", "R/mtblr-convergence.R", "R/mtblr-csr.R",
    "R/RcppExports.R", "src/RcppExports.cpp", "src/mtblr.cpp",
    "src/blr_mt_bed_core_impl.h", "src/blr_mt_bed_types.h",
    "src/blr_mt_bed_chains_types.h", "src/blr_mt_bed_chains_execution_impl.h",
    "src/blr_mt_bed_chains_aggregate_impl.h",
    "src/blr_mt_bed_convergence_types.h",
    "src/blr_mt_bed_convergence_trace_impl.h", "DESCRIPTION", "NAMESPACE",
    "man/mtblr_bed.Rd")
  changed <- system2("git", c("diff", "--name-only", "--", protected),
                     stdout = TRUE)
  expect_length(changed, 0)
  public <- readLines(file.path(root, "R", "mtblr-bed.R"), warn = FALSE)
  expect_false(any(grepl('convergence = c\\("auto", "none", "core", "extended"',
                         public)))
  native <- paste(readLines(file.path(root, "src", "mtblr.cpp"), warn = FALSE),
                  collapse = "\n")
  expect_false(grepl("extended_convergence", native, fixed = TRUE))
  expect_identical(deparse(formals(sblr::mtblr_bed)$convergence),
                   "c(\"auto\", \"none\", \"core\")")
})
