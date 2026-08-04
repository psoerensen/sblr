test_that("SBayesRC invalid scales stop with failure-only diagnostic context", {
  path <- tempfile(fileext = ".tsv")
  on.exit(unlink(path), add = TRUE)

  expect_error(sblr:::.st_sbayesrc_invalid_scale_diagnostic_fixture(path),
  "invalid projected residual scale at iteration 1")
  expect_true(file.exists(path))
  diagnostic <- readLines(path, warn = FALSE)
  expect_true(any(grepl("^meta\\toperator\\tcsr$", diagnostic)))
  expect_true(any(grepl("^meta\\titeration\\t1$", diagnostic)))
  expect_true(any(grepl("^meta\\tmaintained_scale\\t-", diagnostic)))
  expect_true(any(grepl("^meta\\trebuilt_scale\\t-", diagnostic)))
  expect_true(any(grepl("^meta\\tquadratic_scale\\t-", diagnostic)))
  expect_true(any(startsWith(diagnostic, "marker\tindex\teffect")))
})
