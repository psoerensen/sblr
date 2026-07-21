phase9e_path <- function(...) {
  parts <- c(...)
  if (identical(parts[1:3], c("tests", "testthat", "fixtures")))
    do.call(blr_fixture_path, as.list(parts[-(1:3)])) else
    do.call(blr_repo_path, as.list(parts))
}

source(phase9e_path(
  "tests", "testthat", "fixtures", "blr-phase9a-annotation-reference.R"
))

phase9e_active <- function(lines, text) {
  sum(grepl(text, lines, fixed = TRUE) & !grepl("^\\s*//", lines))
}

test_that("canonical group BayesC architecture is singular and build-safe", {
  source_lines <- readLines(phase9e_path("src", "st_cpg_omp_csr_group.cpp"), warn = FALSE)
  core_lines <- readLines(phase9e_path("src", "blr_csr_group_bayesc_core_impl.h"), warn = FALSE)
  type_lines <- readLines(phase9e_path("src", "blr_csr_group_bayesc_types.h"), warn = FALSE)
  source_text <- paste(source_lines, collapse = "\n")
  core_text <- paste(core_lines, collapse = "\n")
  type_text <- paste(type_lines, collapse = "\n")
  all_sources <- list.files(phase9e_path("src"), pattern = "\\.(cpp|h)$", full.names = TRUE)
  includes <- vapply(all_sources, function(path) {
    any(grepl('#include "blr_csr_group_bayesc_core_impl.h"', readLines(path, warn = FALSE), fixed = TRUE))
  }, logical(1))

  expect_equal(phase9e_active(core_lines, "run_csr_group_bayesc("), 1L)
  expect_equal(phase9e_active(core_lines, "for (int isort = 0; isort < m; ++isort)"), 1L)
  expect_equal(phase9e_active(source_lines, "static Rcpp::List stblr_csr_group_bayesc_result_to_raw("), 1L)
  expect_equal(phase9e_active(source_lines, "for (int chain = 0; chain < nchains; ++chain)"), 1L)
  expect_equal(basename(all_sources[includes]), "st_cpg_omp_csr_group.cpp")
  expect_match(core_text, "#ifndef SBLR_BLR_CSR_GROUP_BAYESC_CORE_IMPL_H", fixed = TRUE)
  expect_match(core_text, "SBLR_CSR_GROUP_BAYESC_CORE_IMPL_TRANSLATION_UNIT", fixed = TRUE)
  expect_match(core_text, "Binding-neutral group BayesC implementation detail", fixed = TRUE)
  expect_match(source_text, "run_csr_group_bayesc(context)", fixed = TRUE)
  expect_match(source_text, "context.marker_group=&group", fixed = TRUE)
  expect_match(type_text, "const arma::Row<int>* marker_group", fixed = TRUE)
  expect_match(type_text, "GroupBayesCPolicyView group_policy", fixed = TRUE)
  expect_false(grepl("Rcpp|SEXP|RObject|NumericVector|NumericMatrix|Nullable|pybind11|Python.h",
                     paste(core_text, type_text), perl = TRUE))
  expect_false(grepl("use_old|use_new|old_path|new_path|fallback", source_text,
                     ignore.case = TRUE))

  marker_start <- grep("for (int isort = 0; isort < m; ++isort)", core_lines, fixed = TRUE)
  marker_window <- paste(core_lines[marker_start:min(length(core_lines), marker_start + 100L)],
                         collapse = "\n")
  expect_false(grepl("new |malloc|calloc|realloc|push_back|resize|std::vector<|std::map|unordered_map",
                     marker_window, perl = TRUE))
  expect_match(core_text, "normalize_group_vb", fixed = TRUE)
})

test_that("canonical group raw and formatted references remain exact", {
  expect_length(phase9a_configs$group, 3L)
  for (name in names(phase9a_configs$group)) {
    reference <- readRDS(phase9e_path(
      "tests", "testthat", "fixtures", "blr_phase9a_group", paste0(name, ".rds")
    ))
    config <- phase9a_configs$group[[name]]
    expect_equal(phase9a_normalize(phase9a_run("group", config, TRUE)), reference$raw, tolerance=1e-12,
                     info = paste(name, "raw"))
    expect_equal(phase9a_normalize(phase9a_run("group", config, FALSE)), reference$fit, tolerance=1e-12,
                     info = paste(name, "formatted"))
  }
})

test_that("canonical group route remains exactly reproducible", {
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

test_that("canonical group policy and unsupported cases stay protected", {
  expect_identical(phase9a_configs$group$group_one$normalize, TRUE)
  expect_identical(phase9a_configs$group$group_chains$normalize, FALSE)
  expect_identical(phase9a_configs$group$group_explicit$seeds, c(41L, 42L))
  inputs <- phase9a_inputs(1L)
  expect_identical(unname(inputs$group), c("coding", "background", "coding", "background"))
  contracts <- paste(readLines(
    phase9e_path("src", "blr_csr_annotation_bayesc_types.h"), warn = FALSE
  ), collapse = "\n")
  expect_match(contracts, "group.marker_group must be zero-based", fixed = TRUE)
  expect_match(contracts, "group policy does not permit empty groups", fixed = TRUE)
  expect_match(contracts, "group value dimensions are inconsistent", fixed = TRUE)
})

test_that("canonical group public and protected boundaries remain frozen", {
  protected <- c(
    "NAMESPACE" = "a1f389e8ea9ab5abef440767a11b8378",
    "src/st_cpg_omp_csr_prior.cpp" = "cce51072da6ddc3c18d58ab3b1f3c6df",
    "docs/dev/stblr_raw_schema.md" = "82ac9ba4b7d8edc6f3e16ee3a26d8466"
  )
  actual <- unname(tools::md5sum(vapply(names(protected), phase9e_path, character(1))))
  expect_identical(actual, unname(protected))
})
