test_that("public model semantics separate data level from effect scale", {
  cases <- list(
    list("bayesc", "packed_bed", "bayesc", "individual", "unit"),
    list("sbayesc", "csr", "bayesc", "summary_statistics", "unit"),
    list("bayesr", "packed_bed", "bayesr", "individual", "component"),
    list("sbayesr", "block_eigen", "bayesr", "summary_statistics", "component"),
    list("bayesrc", "packed_bed", "bayesrc", "individual", "component"),
    list("sbayesrc", "csr", "bayesrc", "summary_statistics", "component"))
  for (case in cases) {
    resolved <- sblr:::.blr_model_semantics(case[[1L]], case[[2L]])
    expect_identical(resolved$prior_kernel, case[[3L]])
    expect_identical(resolved$data_level, case[[4L]])
    expect_identical(resolved$effect_scale_policy, case[[5L]])
    expect_false(resolved$maf_effect_s_active)
    expect_identical(resolved$model_semantics_version, 2L)
  }
  expect_identical(
    sblr:::.blr_model_semantics("sbayesr", "csr", maf_effect_s = 0)$effect_scale_policy,
    "component_maf_s")
  expect_identical(
    sblr:::.blr_model_semantics("sbayesc", "csr", maf_effect_s = -1)$effect_scale_policy,
    "maf_s")
})

test_that("operator model matrices are failure-closed before preparation", {
  expect_identical(eval(formals(stblr_csr)$method), c("sbayesc", "sbayesr"))
  expect_identical(eval(formals(stblr_block_eigen)$method),
                   c("sbayesc", "sbayesr", "sbayesrc"))
  expect_identical(eval(formals(stblr_bed)$method),
                   c("bayesc", "bayesr", "bayesrc"))
  expect_identical(formals(mtblr_csr)$method, "sbayesc")
  expect_identical(formals(mtblr_block_eigen)$method, "sbayesc")
  expect_identical(formals(mtblr_bed)$method, "bayesc")
  expect_error(stblr_csr(list(), method = "bayesr"),
               "s.*prefix denotes summary statistics")
  expect_error(stblr_block_eigen(list(), NULL, 1L, method = "bayesr"),
               "s.*prefix denotes summary statistics")
  expect_error(stblr_bed(numeric(), list(), method = "sbayesr"),
               "s.*prefix denotes summary statistics")
  expect_error(mtblr_csr(list(), method = "bayesr"),
               "s.*prefix denotes summary statistics")
  expect_error(mtblr_block_eigen(list(), NULL, 1L, method = "bayesr"),
               "s.*prefix denotes summary statistics")
  expect_error(mtblr_bed(numeric(), list(bedfiles = "x"), method = "sbayesr"),
               "s.*prefix denotes summary statistics")
})

test_that("maf_effect_s NULL and minus one are numerically equal but semantic distinct", {
  patterns <- sblr:::.mtblr_models(matrix(c(0L, 1L), 2L, 1L),
                                    c(.8, .2), .2, 1L)
  base <- sblr:::.mtblr_bayesr_spec(
    "bayesr", patterns, NULL, 2L, c(0, 1))
  unit <- sblr:::.mtblr_bayesr_spec(
    "bayesr", patterns, c(.2, .4), 2L, c(0, 1), maf_effect_s = -1)
  expect_equal(base$marker_scale, unit$marker_scale, tolerance = 0)
  expect_null(base$maf_effect_s)
  expect_identical(unit$maf_effect_s, -1)
  expect_false(base$model_parameters$mixture$maf_effect_s_active)
  expect_true(unit$model_parameters$mixture$maf_effect_s_active)
})

test_that("selection MAF provenance is explicit and fallback is opt-in", {
  explicit <- sblr:::.mtblr_resolve_effect_maf(
    c(.1, .2), TRUE, 2L)
  expect_identical(explicit$effect_maf_source, "explicit_effect_maf")
  expect_false(explicit$effect_maf_fallback_used)
  reference <- data.frame(marker_id = c("a", "b"),
                          allele_frequency = c(.1, .2))
  expect_error(sblr:::.mtblr_resolve_effect_maf(
    NULL, TRUE, 2L, reference_marker_metadata = reference),
    "fallback is disabled")
  fallback <- sblr:::.mtblr_resolve_effect_maf(
    NULL, TRUE, 2L, reference_marker_metadata = reference,
    allow_reference_maf_for_maf_effect_s = TRUE)
  expect_identical(fallback$effect_maf_source,
                   "reference_panel_frequency")
  expect_true(fallback$effect_maf_fallback_used)
})

test_that("ambiguous old fits are not silently reinterpreted", {
  old <- list(input = list(model = "sbayesr"))
  expect_error(sblr:::.blr_validate_fit_model_semantics(old),
               "cannot be silently reinterpreted")
  current <- list(input = list(
    model_semantics_version = 2L,
    model_semantics = "s_prefix_means_summary_statistics"))
  expect_invisible(sblr:::.blr_validate_fit_model_semantics(current))
})

test_that("annotation policies retain independent summary model labels", {
  for (policy in c("prior", "learned", "group")) {
    fit <- list(input = list())
    standardized <- sblr:::.standardize_stblr_annotation_fit(fit, policy)
    expect_identical(standardized$input$model, "sbayesc")
  }
  mixture <- sblr:::.standardize_stblr_annotation_fit(
    list(input = list()), "sbayesrc")
  expect_identical(mixture$input$model, "sbayesrc")
})
