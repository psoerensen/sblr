phase9b2_path <- function(...) {
  path <- file.path(...)
  if (file.exists(path)) path else file.path("..", "..", ...)
}

source(phase9b2_path(
  "tests", "testthat", "fixtures", "blr-phase9a-annotation-reference.R"
))

test_that("Phase 9B2 exposes an explicit fixed-prior callable boundary", {
  types <- paste(readLines(phase9b2_path(
    "src", "blr_csr_prior_bayesc_types.h"
  ), warn = FALSE), collapse = "\n")
  core <- paste(readLines(phase9b2_path(
    "src", "blr_csr_prior_bayesc_core_impl.h"
  ), warn = FALSE), collapse = "\n")
  source_text <- paste(readLines(phase9b2_path(
    "src", "st_cpg_omp_csr_prior.cpp"
  ), warn = FALSE), collapse = "\n")

  expect_match(types, "struct CsrPriorBayesCExecutionContext", fixed = TRUE)
  expect_match(types, "struct CsrPriorBayesCExecutionResult", fixed = TRUE)
  expect_match(core, "run_csr_prior_bayesc(", fixed = TRUE)
  expect_match(source_text, "return sblr::core::run_csr_prior_bayesc(context);", fixed = TRUE)
  expect_false(grepl("Rcpp::|SEXP", paste(types, core), perl = TRUE))
})

test_that("Phase 9B2 retains one active numerical marker loop", {
  files <- vapply(c(
    "st_cpg_omp_csr_prior.cpp", "blr_csr_prior_bayesc_core_impl.h"
  ), function(name) phase9b2_path("src", name), character(1))
  marker <- "for (int isort = 0; isort < m; ++isort)"
  active <- vapply(files, function(path) {
    lines <- readLines(path, warn = FALSE)
    sum(grepl(marker, lines, fixed = TRUE) & !grepl("^\\s*//", lines))
  }, integer(1))
  expect_equal(sum(active), 1L)
  expect_equal(unname(active[basename(files) == "blr_csr_prior_bayesc_core_impl.h"]), 1L)
})

test_that("Phase 9B2 fixed-prior frozen references remain exact", {
  for (name in names(phase9a_configs$prior)) {
    reference <- readRDS(phase9b2_path(
      "tests", "testthat", "fixtures", "blr_phase9a_prior",
      paste0(name, ".rds")
    ))
    config <- phase9a_configs$prior[[name]]
    expect_identical(
      phase9a_normalize(phase9a_run("prior", config, TRUE)), reference$raw,
      info = paste(name, "raw")
    )
    expect_identical(
      phase9a_normalize(phase9a_run("prior", config, FALSE)), reference$fit,
      info = paste(name, "formatted")
    )
  }
})

test_that("Phase 9B2 leaves wrapper-level multichain aggregation in place", {
  source_text <- paste(readLines(phase9b2_path(
    "src", "st_cpg_omp_csr_prior.cpp"
  ), warn = FALSE), collapse = "\n")
  expect_match(source_text, "for (int chain = 0; chain < nchains; ++chain)", fixed = TRUE)
  expect_match(source_text, "stblr_csr_prior_bayesc_result_to_raw", fixed = TRUE)
  expect_match(source_text, "cpg_prior_chains_raw_v1", fixed = TRUE)
})
