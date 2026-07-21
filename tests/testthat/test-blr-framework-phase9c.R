phase9c_path <- function(...) {
  parts <- c(...)
  if (identical(parts[1:3], c("tests", "testthat", "fixtures")))
    do.call(blr_fixture_path, as.list(parts[-(1:3)])) else
    do.call(blr_repo_path, as.list(parts))
}

source(phase9c_path(
  "tests", "testthat", "fixtures", "blr-phase9a-annotation-reference.R"
))

phase9c_active <- function(lines, text) {
  sum(grepl(text, lines, fixed = TRUE) & !grepl("^\\s*//", lines))
}

test_that("Phase 9C canonical architecture is singular and binding safe", {
  source_lines <- readLines(phase9c_path("src", "st_cpg_omp_csr_prior.cpp"), warn = FALSE)
  core_lines <- readLines(phase9c_path("src", "blr_csr_prior_bayesc_core_impl.h"), warn = FALSE)
  types_lines <- readLines(phase9c_path("src", "blr_csr_prior_bayesc_types.h"), warn = FALSE)
  source_text <- paste(source_lines, collapse = "\n")
  core_text <- paste(core_lines, collapse = "\n")
  types_text <- paste(types_lines, collapse = "\n")

  expect_equal(phase9c_active(core_lines, "run_csr_prior_bayesc("), 1L)
  expect_equal(phase9c_active(core_lines, "for (int isort = 0; isort < m; ++isort)"), 1L)
  expect_equal(phase9c_active(source_lines, "static Rcpp::List stblr_csr_prior_bayesc_result_to_raw("), 1L)
  expect_equal(phase9c_active(source_lines, "for (int chain = 0; chain < nchains; ++chain)"), 1L)
  expect_equal(phase9c_active(source_lines, "#include \"blr_csr_prior_bayesc_core_impl.h\""), 1L)
  expect_match(core_text, "#ifndef SBLR_BLR_CSR_PRIOR_BAYESC_CORE_IMPL_H", fixed = TRUE)
  expect_match(core_text, "inline CsrPriorBayesCExecutionResult run_csr_prior_bayesc", fixed = TRUE)
  expect_false(grepl("Rcpp::|SEXP|PyObject|pybind11", paste(core_text, types_text), perl = TRUE))
  expect_false(grepl("use_old|use_new|old_path|new_path|fallback", source_text, ignore.case = TRUE))
  expect_match(types_text, "const arma::mat* marker_probability", fixed = TRUE)
  expect_match(types_text, "const arma::mat* marker_multiplier", fixed = TRUE)
  expect_match(types_text, "const void* ld_storage", fixed = TRUE)
  expect_match(core_text, "std::mt19937 gen_t", fixed = TRUE)
})

test_that("Phase 9C permanent fixed-prior references remain exact", {
  for (name in names(phase9a_configs$prior)) {
    reference <- readRDS(phase9c_path(
      "tests", "testthat", "fixtures", "blr_phase9a_prior", paste0(name, ".rds")
    ))
    config <- phase9a_configs$prior[[name]]
    expect_equal(phase9a_normalize(phase9a_run("prior", config, TRUE)), reference$raw, tolerance=1e-12,
                     info = paste(name, "raw"))
    expect_equal(phase9a_normalize(phase9a_run("prior", config, FALSE)), reference$fit, tolerance=1e-12,
                     info = paste(name, "formatted"))
  }
})

test_that("Phase 9C canonical fixed-prior execution is reproducible", {
  comparable <- function(value) {
    value <- phase9a_normalize(value)
    value$input$ncores <- 0L
    value
  }
  config <- phase9a_configs$prior$fixed_chains
  config$ncores <- 1L
  one <- comparable(phase9a_run("prior", config, FALSE))
  expect_identical(comparable(phase9a_run("prior", config, FALSE)), one)
  config$ncores <- 2L
  two <- comparable(phase9a_run("prior", config, FALSE))
  expect_identical(two, one)
  expect_identical(comparable(phase9a_run("prior", config, FALSE)), two)
  config$ncores <- 1L
  expect_identical(comparable(phase9a_run("prior", config, FALSE)), one)
  invisible(phase9a_run("annotation", phase9a_configs$annotation$annot_fixed, FALSE))
  expect_identical(comparable(phase9a_run("prior", config, FALSE)), one)
})

test_that("Phase 9C preserves the current trait-dimension rejection", {
  config <- phase9a_configs$prior$fixed_one
  config$nt <- 2L
  expect_error(phase9a_run("prior", config, FALSE))
})

test_that("Phase 9C protects canonical and adjacent implementations", {
  protected <- c(
    "src/blr_csr_bayesc_types.h" = "4d0eb5380007195a8d34e7b2e081dec4",
    "src/blr_csr_bayesc_core_impl.h" = "c548157cc9e5804272e714983bdcb798",
    "NAMESPACE" = "a1f389e8ea9ab5abef440767a11b8378"
  )
  actual <- unname(tools::md5sum(vapply(names(protected), phase9c_path, character(1))))
  expect_identical(actual, unname(protected))
})
