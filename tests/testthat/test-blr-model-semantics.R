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
  expect_identical(formals(mtblr_bed)$method, "bayesc")
  expect_error(stblr_csr(list(), method = "bayesr"),
               "s.*prefix denotes summary statistics")
  expect_error(stblr_block_eigen(list(), NULL, 1L, method = "bayesr"),
               "s.*prefix denotes summary statistics")
  expect_error(stblr_bed(numeric(), list(), method = "sbayesr"),
               "s.*prefix denotes summary statistics")
  expect_identical(formals(mtblr_csr)$method, "bayesc")
  expect_identical(formals(mtblr_block_eigen)$method, "bayesc")
  expect_error(mtblr_bed(numeric(), list(bedfiles = "x"),
                         method = "sbayesr"),
               "method")
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
