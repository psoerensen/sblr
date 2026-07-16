phase9d3_path <- function(...) {
  path <- file.path(...)
  if (file.exists(path)) path else file.path("..", "..", ...)
}

source(phase9d3_path(
  "tests", "testthat", "fixtures", "blr-phase9a-annotation-reference.R"
))

phase9d3_active <- function(lines, text) {
  sum(grepl(text, lines, fixed = TRUE) & !grepl("^\\s*//", lines))
}

test_that("Phase 9D3 migrated group architecture is singular", {
  source_lines <- readLines(phase9d3_path("src", "st_cpg_omp_csr_group.cpp"), warn = FALSE)
  core_lines <- readLines(phase9d3_path("src", "blr_csr_group_bayesc_core_impl.h"), warn = FALSE)
  type_lines <- readLines(phase9d3_path("src", "blr_csr_group_bayesc_types.h"), warn = FALSE)
  source_text <- paste(source_lines, collapse = "\n")
  core_text <- paste(core_lines, collapse = "\n")
  type_text <- paste(type_lines, collapse = "\n")

  expect_equal(phase9d3_active(core_lines, "run_csr_group_bayesc("), 1L)
  expect_equal(phase9d3_active(core_lines, "for (int isort = 0; isort < m; ++isort)"), 1L)
  expect_equal(phase9d3_active(source_lines, "stblr_csr_group_bayesc_result_to_raw("), 3L)
  expect_equal(phase9d3_active(source_lines, "static Rcpp::List stblr_csr_group_bayesc_result_to_raw("), 1L)
  expect_equal(phase9d3_active(source_lines, "for (int chain = 0; chain < nchains; ++chain)"), 1L)
  expect_equal(phase9d3_active(source_lines, "run_csr_group_bayesc(context)"), 1L)
  expect_match(source_text, "const sblr::core::CsrGroupBayesCExecutionResult& execution_result", fixed = TRUE)
  expect_match(source_text, "context.marker_group=&group", fixed = TRUE)
  expect_match(type_text, "GroupBayesCPolicyView group_policy", fixed = TRUE)
  expect_match(type_text, "const arma::Row<int>* marker_group", fixed = TRUE)
  expect_match(core_text, "sampleGroupVbMultipliers_ST_csr_group(", fixed = TRUE)
  expect_match(core_text, "normalize_group_vb", fixed = TRUE)
  expect_false(grepl("Rcpp|SEXP|RObject|NumericVector|NumericMatrix|Nullable|pybind11|Python.h",
                     paste(core_text, type_text), perl = TRUE))
  expect_false(grepl("use_old|use_new|old_path|new_path|fallback", source_text,
                     ignore.case = TRUE))
})

test_that("Phase 9D3 permanent group references remain exact", {
  expect_length(phase9a_configs$group, 3L)
  for (name in names(phase9a_configs$group)) {
    reference <- readRDS(phase9d3_path(
      "tests", "testthat", "fixtures", "blr_phase9a_group", paste0(name, ".rds")
    ))
    config <- phase9a_configs$group[[name]]
    expect_identical(phase9a_normalize(phase9a_run("group", config, TRUE)), reference$raw,
                     info = paste(name, "raw"))
    expect_identical(phase9a_normalize(phase9a_run("group", config, FALSE)), reference$fit,
                     info = paste(name, "formatted"))
  }
})

test_that("Phase 9D3 group route remains reproducible", {
  comparable <- function(value) {
    value <- phase9a_normalize(value)
    value$input$ncores <- 0L
    value
  }
  config <- phase9a_configs$group$group_chains
  config$ncores <- 1L
  one <- comparable(phase9a_run("group", config, FALSE))
  expect_identical(comparable(phase9a_run("group", config, FALSE)), one)
  config$ncores <- 2L
  two <- comparable(phase9a_run("group", config, FALSE))
  expect_identical(two, one)
  expect_identical(comparable(phase9a_run("group", config, FALSE)), two)
  config$ncores <- 1L
  expect_identical(comparable(phase9a_run("group", config, FALSE)), one)
  invisible(phase9a_run("annotation", phase9a_configs$annotation$annot_fixed, FALSE))
  expect_identical(comparable(phase9a_run("group", config, FALSE)), one)
})

test_that("Phase 9D3 preserves policy coverage and restrictions", {
  expect_identical(phase9a_configs$group$group_one$normalize, TRUE)
  expect_identical(phase9a_configs$group$group_chains$normalize, FALSE)
  expect_identical(phase9a_configs$group$group_explicit$seeds, c(41L, 42L))
  inputs <- phase9a_inputs(1L)
  expect_identical(unname(inputs$group), c("coding", "background", "coding", "background"))
  contracts <- paste(readLines(
    phase9d3_path("src", "blr_csr_annotation_bayesc_types.h"), warn = FALSE
  ), collapse = "\n")
  expect_match(contracts, "group.marker_group must be zero-based", fixed = TRUE)
  expect_match(contracts, "group policy does not permit empty groups", fixed = TRUE)
  expect_match(contracts, "group value dimensions are inconsistent", fixed = TRUE)
})
