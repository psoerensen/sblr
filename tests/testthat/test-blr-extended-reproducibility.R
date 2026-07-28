test_that("representative canonical BLR routes are fresh-process reproducible", {
  skip_if_not(identical(Sys.getenv("SBLR_RUN_EXTENDED_REPRODUCIBILITY"), "true"))
  skip_if_not_installed("callr")
  root <- blr_repo_path()
  observed <- callr::r(function(root) {
    setwd(root); pkgload::load_all(".", compile = FALSE, quiet = TRUE)
    source("tests/testthat/fixtures/st-bayesc-csr-reference.R")
    source("tests/testthat/fixtures/st-bayesc-scheduled-csr-reference.R")
    source("tests/testthat/fixtures/st-bayesc-bed-reference.R")
    source("tests/testthat/fixtures/mt-bayesc-reference.R")
    scalar <- st_bayesc_csr_reference_objects(st_bayesc_csr_reference_configurations$one_trait_one_chain_one_core)
    scheduled <- st_bayesc_scheduled_run(st_bayesc_scheduled_configs$skip_two_one)
    bed <- st_bayesc_bed_capture("multichain", 1L, 2L, 71L)
    mt <- mt_bayesc_reference_capture(1L, TRUE)
    old <- Sys.getenv("OMP_NUM_THREADS", unset = NA_character_)
    on.exit(if (is.na(old)) Sys.unsetenv("OMP_NUM_THREADS") else Sys.setenv(OMP_NUM_THREADS = old), add = TRUE)
    Sys.setenv(OMP_NUM_THREADS = "1"); mt_one <- mt_bayesc_reference_capture(1L, TRUE)
    Sys.setenv(OMP_NUM_THREADS = "2"); mt_two <- mt_bayesc_reference_capture(1L, TRUE)
    list(scalar=scalar, scheduled=scheduled, bed=bed, mt=mt, mt_one=mt_one, mt_two=mt_two)
  }, list(root = root))
  oldwd <- setwd(root)
  on.exit(setwd(oldwd), add = TRUE)
  source("tests/testthat/fixtures/st-bayesc-csr-reference.R", local = TRUE)
  source("tests/testthat/fixtures/st-bayesc-scheduled-csr-reference.R", local = TRUE)
  source("tests/testthat/fixtures/st-bayesc-bed-reference.R", local = TRUE)
  source("tests/testthat/fixtures/mt-bayesc-reference.R", local = TRUE)
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
    st_bayesc_csr_reference_objects(
      st_bayesc_csr_reference_configurations$one_trait_one_chain_one_core)),
    tolerance=1e-12)
  expect_equal(observed$scheduled, st_bayesc_scheduled_run(st_bayesc_scheduled_configs$skip_two_one), tolerance=1e-12)
  expect_equal(observed$bed, st_bayesc_bed_capture("multichain", 1L, 2L, 71L), tolerance=1e-12)
  expect_equal(observed$mt, mt_bayesc_reference_capture(1L, TRUE), tolerance=1e-12)
  expect_equal(observed$mt_one, observed$mt_two, tolerance=1e-12)
})
