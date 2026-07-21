source(blr_fixture_path("blr_phase17b_mt_default",
  "blr-phase17b-mt-default-reference.R"))

test_that("Phase 17B historical multivariate fixtures remain immutable", {
  paths <- file.path("tests/testthat/fixtures/blr_phase17b_mt_default",
    sprintf("config-%d.rds", 1:3))
  blr_expect_fixture_md5(paths, c(
    "2b820e1cdd9e731f1f0ffea442ef4e53",
    "d6d8abec35168a088a9accab87b3c6d0",
    "48fe8040d52b1d23c6e7437d632ebebf"))
  for (path in paths) {
    ref <- readRDS(blr_fixture_path("blr_phase17b_mt_default", basename(path)))
    expect_identical(length(ref$raw), 20L)
    expect_identical(ref$metadata$reference_mode,
      "structure_exact_numeric_tolerance")
    expect_identical(ref$metadata$numeric_tolerance, 1e-12)
  }
})
test_that("Phase 17B preserves the historical fixed-B defect evidence", {
  ref <- readRDS(blr_fixture_path("blr_phase17b_mt_default", "config-2.rds"))
  config <- phase17b_mt_config(2L)
  dimnames(config$vb) <- dimnames(ref$fit$vb)
  expect_false(isTRUE(all.equal(ref$fit$vb, config$vb, tolerance = 0)))
  expect_match(blr_source_text("docs/dev/blr_framework_phase17b_report.md"),
    "updateB = FALSE does not keep B fixed", fixed = TRUE)
})

test_that("Phase 17B preserves the historical denominator evidence", {
  ref <- readRDS(blr_fixture_path("blr_phase17b_mt_default", "config-1.rds"))
  config <- phase17b_mt_config(1L)
  accumulated <- sum(0:(config$nit + config$nburn - 1L) > config$nburn)
  expect_identical(accumulated, config$nit - 1L)
  expect_equal(sum(unname(ref$fit$pim)), accumulated / config$nit,
    tolerance = 1e-12)
  expect_false(isTRUE(all.equal(sum(unname(ref$fit$pim)), 1,
    tolerance = 1e-12)))
})
