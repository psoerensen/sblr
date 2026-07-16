phase9b3_path <- function(...) {
  path <- file.path(...)
  if (file.exists(path)) path else file.path("..", "..", ...)
}

source(phase9b3_path(
  "tests", "testthat", "fixtures", "blr-phase9a-annotation-reference.R"
))

phase9b3_active <- function(lines, text) {
  sum(grepl(text, lines, fixed = TRUE) & !grepl("^\\s*//", lines))
}

test_that("Phase 9B3 has one core, marker loop, converter, and aggregation", {
  source_lines <- readLines(phase9b3_path("src", "st_cpg_omp_csr_prior.cpp"), warn = FALSE)
  core_lines <- readLines(phase9b3_path("src", "blr_csr_prior_bayesc_core_impl.h"), warn = FALSE)
  types <- paste(readLines(phase9b3_path("src", "blr_csr_prior_bayesc_types.h"), warn = FALSE), collapse = "\n")
  source_text <- paste(source_lines, collapse = "\n")
  core_text <- paste(core_lines, collapse = "\n")

  expect_equal(phase9b3_active(core_lines, "run_csr_prior_bayesc("), 1L)
  expect_equal(phase9b3_active(core_lines, "for (int isort = 0; isort < m; ++isort)"), 1L)
  expect_equal(phase9b3_active(source_lines, "static Rcpp::List stblr_csr_prior_bayesc_result_to_raw("), 1L)
  expect_equal(phase9b3_active(source_lines, "for (int chain = 0; chain < nchains; ++chain)"), 1L)
  expect_equal(phase9b3_active(source_lines, "#include \"blr_csr_prior_bayesc_core_impl.h\""), 1L)
  expect_false(grepl("legacy|fallback|use_old|use_new|pybind11", paste(core_text, types), ignore.case = TRUE))
  expect_false(grepl("Rcpp::|SEXP|PyObject", paste(core_text, types), perl = TRUE))
  expect_match(source_text, "return sblr::core::run_csr_prior_bayesc(context);", fixed = TRUE)
  expect_match(types, "const arma::mat* marker_probability", fixed = TRUE)
  expect_match(types, "const arma::mat* marker_multiplier", fixed = TRUE)
})

test_that("Phase 9B3 fixed-prior references remain exact", {
  for (name in names(phase9a_configs$prior)) {
    reference <- readRDS(phase9b3_path(
      "tests", "testthat", "fixtures", "blr_phase9a_prior", paste0(name, ".rds")
    ))
    config <- phase9a_configs$prior[[name]]
    expect_identical(phase9a_normalize(phase9a_run("prior", config, TRUE)), reference$raw,
                     info = paste(name, "raw"))
    expect_identical(phase9a_normalize(phase9a_run("prior", config, FALSE)), reference$fit,
                     info = paste(name, "formatted"))
  }
})

test_that("Phase 9B3 reproducibility is exact across scheduling and call order", {
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

test_that("Phase 9B3 preserves unsupported multiple-trait behavior", {
  config <- phase9a_configs$prior$fixed_one
  config$nt <- 2L
  expect_error(phase9a_run("prior", config, FALSE))
})

test_that("Phase 9B3 keeps public and protected surfaces unchanged", {
  expect_identical(
    unname(tools::md5sum(phase9b3_path("NAMESPACE"))),
    "f5b6ee37a3972aa436357bdc8f602f4e"
  )
  protected <- c(
    "src/st_cpg_omp_csr_group.cpp" = "f9701762d2e0245a40e996c89a4addb2",
    "src/st_cpg_omp_csr_annot.cpp" = "baaf3a0919ba97c78401066f7ac7d6f3"
  )
  actual <- unname(tools::md5sum(vapply(names(protected), phase9b3_path, character(1))))
  expect_identical(actual, unname(protected))
})
