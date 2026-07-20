test_that("representative canonical BLR routes are fresh-process reproducible", {
  skip_if_not(identical(Sys.getenv("SBLR_RUN_EXTENDED_REPRODUCIBILITY"), "true"))
  skip_if_not_installed("callr")
  observed <- callr::r(function(root) {
    setwd(root); pkgload::load_all(".", compile = FALSE, quiet = TRUE)
    source("tests/testthat/fixtures/blr-phase2-reference.R")
    source("tests/testthat/fixtures/blr-phase10b-scheduled-reference.R")
    source("tests/testthat/fixtures/blr-phase11b-bed-bayesc-reference.R")
    source("tests/testthat/fixtures/blr_phase17b_mt_default/blr-phase17b-mt-default-reference.R")
    source("tests/testthat/fixtures/blr_phase17c_mt_default_corrected/blr-phase17c-mt-default-corrected-reference.R")
    scalar <- phase2_reference_objects(phase2_reference_configurations$one_trait_one_chain_one_core)
    scheduled <- phase10b_run(phase10b_configs$skip_two_one)
    bed <- phase11b_capture("multichain", 1L, 2L, 71L)
    mt <- phase17c_mt_capture(1L, TRUE)
    old <- Sys.getenv("OMP_NUM_THREADS", unset = NA_character_)
    on.exit(if (is.na(old)) Sys.unsetenv("OMP_NUM_THREADS") else Sys.setenv(OMP_NUM_THREADS = old), add = TRUE)
    Sys.setenv(OMP_NUM_THREADS = "1"); mt_one <- phase17c_mt_capture(1L, TRUE)
    Sys.setenv(OMP_NUM_THREADS = "2"); mt_two <- phase17c_mt_capture(1L, TRUE)
    list(scalar=scalar, scheduled=scheduled, bed=bed, mt=mt, mt_one=mt_one, mt_two=mt_two)
  }, list(root = blr_test_root))
  oldwd <- setwd(blr_test_root)
  on.exit(setwd(oldwd), add = TRUE)
  source("tests/testthat/fixtures/blr-phase2-reference.R", local = TRUE)
  source("tests/testthat/fixtures/blr-phase10b-scheduled-reference.R", local = TRUE)
  source("tests/testthat/fixtures/blr-phase11b-bed-bayesc-reference.R", local = TRUE)
  source("tests/testthat/fixtures/blr_phase17b_mt_default/blr-phase17b-mt-default-reference.R", local = TRUE)
  source("tests/testthat/fixtures/blr_phase17c_mt_default_corrected/blr-phase17c-mt-default-corrected-reference.R", local = TRUE)
  normalize_scalar <- function(x) {
    if (!is.list(x)) return(x)
    for (name in names(x)) {
      if (name %in% c("seconds_mean", "seconds_max")) x[[name]][] <- 0
      else if (name == "ld_prefix") x[[name]] <- "<fixture>"
      else x[[name]] <- normalize_scalar(x[[name]])
    }
    x
  }
  expect_equal(normalize_scalar(observed$scalar), normalize_scalar(
    phase2_reference_objects(
      phase2_reference_configurations$one_trait_one_chain_one_core)),
    tolerance=1e-12)
  expect_equal(observed$scheduled, phase10b_run(phase10b_configs$skip_two_one), tolerance=1e-12)
  expect_equal(observed$bed, phase11b_capture("multichain", 1L, 2L, 71L), tolerance=1e-12)
  expect_equal(observed$mt, phase17c_mt_capture(1L, TRUE), tolerance=1e-12)
  expect_equal(observed$mt_one, observed$mt_two, tolerance=1e-12)
})
