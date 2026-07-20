source(file.path(blr_test_root,
  "tests/testthat/fixtures/blr_phase17b_mt_default",
  "blr-phase17b-mt-default-reference.R"))
source(file.path(blr_test_root,
  "tests/testthat/fixtures/blr_phase17c_mt_default_corrected",
  "blr-phase17c-mt-default-corrected-reference.R"))

phase17c_reference <- function(id) readRDS(file.path(blr_test_root,
  "tests/testthat/fixtures/blr_phase17c_mt_default_corrected",
  sprintf("config-%d.rds", id)))

test_that("Phase 17C uniquely owns corrected raw and formatted references", {
  for (id in 1:3) {
    ref <- phase17c_reference(id)
    expect_equal(phase17c_mt_capture(id, FALSE), ref$raw, tolerance = 1e-12)
    expect_equal(phase17c_mt_capture(id, TRUE), ref$fit, tolerance = 1e-12)
    expect_identical(length(ref$raw), 20L)
    expect_identical(ref$metadata$reference_mode,
      "structure_exact_numeric_tolerance")
    expect_identical(ref$metadata$numeric_tolerance, 1e-12)
  }
})
test_that("Phase 17C owns fixed update-control behavior", {
  config <- phase17c_mt_config(2L)
  fit <- phase17c_mt_capture(2L, TRUE)
  expect_identical(unname(fit$vb), unname(config$vb))
  expect_identical(unname(fit$ve), unname(config$ve))
  expect_identical(unname(fit$pi), c(.8, rep(.2 / 3, 3)))
  expect_true(all(fit$covb == 0))
  expect_true(all(fit$cove == 0))
  expect_true(all(fit$pim == 0))
  expect_true(all(is.finite(unlist(fit))))
})

test_that("Phase 17C owns retained-count and normalization behavior", {
  for (id in 1:3) {
    config <- phase17c_mt_config(id)
    expected <- phase17c_mt_retained_counts(config)
    metadata <- phase17c_reference(id)$metadata
    for (field in names(expected))
      expect_identical(metadata[[field]], expected[[field]])
    expect_identical(expected$marker_retained_count,
      as.integer(ceiling(config$nit / config$nthin)))
  }
  for (id in c(1L, 3L)) {
    fit <- phase17c_mt_capture(id, TRUE)
    config <- phase17c_mt_config(id)
    retained <- (config$nburn + 1L):(config$nburn + config$nit)
    expect_equal(sum(unname(fit$pim)), 1, tolerance = 1e-12)
    expect_equal(diag(fit$covb), colMeans(fit$vbs[retained, , drop = FALSE]),
      tolerance = 1e-12)
    expect_equal(diag(fit$covg), colMeans(fit$vgs[retained, , drop = FALSE]),
      tolerance = 1e-12)
    expect_equal(diag(fit$cove), colMeans(fit$ves[retained, , drop = FALSE]),
      tolerance = 1e-12)
  }
})

test_that("Phase 17C owns corrected scientific identities", {
  for (id in 1:3) {
    fit <- phase17c_mt_capture(id, TRUE)
    expect_identical(dim(fit$bm), c(4L, 2L))
    expect_identical(rownames(fit$bm), paste0("M", 1:4))
    expect_identical(colnames(fit$bm), c("TraitA", "TraitB"))
    expect_true(all(is.finite(unlist(fit))))
    expect_true(all(fit$dm >= 0 & fit$dm <= 1))
    expect_true(all(fit$d %in% 0:1))
    for (field in c("covb", "covg", "cove", "vb", "vg", "ve")) {
      expect_equal(fit[[field]], t(fit[[field]]), tolerance = 1e-12)
      expect_true(min(eigen(fit[[field]], symmetric = TRUE,
        only.values = TRUE)$values) > -1e-8)
    }
  }
})

test_that("Phase 17C owns current same-process determinism", {
  a <- phase17c_mt_capture(1L, TRUE)
  expect_equal(phase17c_mt_capture(1L, TRUE), a, tolerance = 1e-12)
  invisible(phase17c_mt_capture(2L, TRUE))
  expect_equal(phase17c_mt_capture(1L, TRUE), a, tolerance = 1e-12)
})
