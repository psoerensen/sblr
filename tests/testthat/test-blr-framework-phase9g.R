phase9g_path <- function(...) {
  parts <- c(...)
  if (identical(parts[1:3], c("tests", "testthat", "fixtures")))
    do.call(blr_fixture_path, as.list(parts[-(1:3)])) else
    do.call(blr_repo_path, as.list(parts))
}

source(phase9g_path(
  "tests", "testthat", "fixtures", "blr-phase9a-annotation-reference.R"
))

phase9g_active <- function(lines, text) {
  sum(grepl(text, lines, fixed = TRUE) & !grepl("^\\s*//", lines))
}

test_that("Phase 9G permanently protects the canonical architecture", {
  source_lines <- readLines(phase9g_path("src", "st_cpg_omp_csr_annot.cpp"), warn = FALSE)
  core_lines <- readLines(phase9g_path(
    "src", "blr_csr_learned_annotation_bayesc_core_impl.h"
  ), warn = FALSE)
  type_lines <- readLines(phase9g_path(
    "src", "blr_csr_learned_annotation_bayesc_types.h"
  ), warn = FALSE)
  source_text <- paste(source_lines, collapse = "\n")
  core_text <- paste(core_lines, collapse = "\n")
  type_text <- paste(type_lines, collapse = "\n")
  all_sources <- list.files(phase9g_path("src"), pattern = "\\.(cpp|h)$", full.names = TRUE)
  includes <- vapply(all_sources, function(path) {
    any(grepl('#include "blr_csr_learned_annotation_bayesc_core_impl.h"',
              readLines(path, warn = FALSE), fixed = TRUE))
  }, logical(1))

  expect_match(core_lines[[1]], "#ifndef SBLR_BLR_CSR_LEARNED_ANNOTATION_BAYESC_CORE_IMPL_H", fixed = TRUE)
  expect_equal(phase9g_active(core_lines, "run_csr_learned_annotation_bayesc("), 1L)
  expect_equal(phase9g_active(core_lines, "for (int isort = 0; isort < m; ++isort)"), 1L)
  expect_equal(phase9g_active(source_lines,
                              "stblr_csr_learned_annotation_bayesc_result_to_raw("), 3L)
  expect_equal(phase9g_active(source_lines,
                              "for (int chain = 0; chain < nchains; ++chain)"), 1L)
  expect_equal(basename(all_sources[includes]), "st_cpg_omp_csr_annot.cpp")
  expect_match(type_text, "LearnedAnnotationBayesCPolicyView annotation_policy", fixed = TRUE)
  expect_match(type_text, "const arma::mat* annotation", fixed = TRUE)
  expect_match(type_text, "const void* ld_storage", fixed = TRUE)
  expect_match(type_text, "Trait-local sampler, coefficient, RNG, and accumulator state", fixed = TRUE)
  expect_match(source_text, 'annotation_policy.layout="column_major_marker_by_annotation"', fixed = TRUE)
  expect_match(source_text, "annotation_policy.includes_intercept=false", fixed = TRUE)
  expect_false(grepl("use_old|use_new|old_path|new_path|fallback", source_text,
                     ignore.case = TRUE))
  expect_false(grepl("Rcpp|SEXP|RObject|NumericVector|NumericMatrix|Nullable|pybind11|Python.h",
                     paste(core_text, type_text), perl = TRUE))
  expect_false(grepl("new |malloc|calloc|realloc|resize\\(",
                     paste(core_lines[240:330], collapse = "\n"), perl = TRUE))
})

test_that("Phase 9G permanently protects the learned policy", {
  core_text <- paste(readLines(phase9g_path(
    "src", "blr_csr_learned_annotation_bayesc_core_impl.h"
  ), warn = FALSE), collapse = "\n")
  adapter_lines <- readLines(phase9g_path("src", "st_cpg_omp_csr_annot.cpp"), warn = FALSE)
  adapter_start <- grep("^Rcpp::List stblr_cpg_omp_csr_annot\\(", adapter_lines)
  adapter_end <- grep("^// // \\[\\[Rcpp::depends", adapter_lines)
  adapter_text <- paste(adapter_lines[adapter_start:(adapter_end[[1]] - 1L)], collapse = "\n")

  expect_match(core_text, "make_pi_from_annotation(A, eta_pi_t, pi_t[1], pi_min, pi_max)", fixed = TRUE)
  expect_match(core_text, "make_vb_multiplier_from_annotation(", fixed = TRUE)
  expect_match(core_text, "rw_update_eta_pi(", fixed = TRUE)
  expect_match(core_text, "rw_update_eta_vb(", fixed = TRUE)
  expect_match(core_text, "(it + 1) % annot_update_every == 0", fixed = TRUE)
  expect_match(core_text, "std::mt19937 gen_t", fixed = TRUE)
  expect_false(grepl("make_pi_from_annotation\\(|make_vb_multiplier_from_annotation\\(|rw_update_eta_",
                     adapter_text, perl = TRUE))
})

test_that("Phase 9G permanent learned-annotation references remain exact", {
  expect_length(phase9a_configs$annotation, 3L)
  for (name in names(phase9a_configs$annotation)) {
    reference <- readRDS(phase9g_path(
      "tests", "testthat", "fixtures", "blr_phase9a_annotation", paste0(name, ".rds")
    ))
    config <- phase9a_configs$annotation[[name]]
    tol <- if (name == "annot_learned") 1e-8 else 1e-12
    expect_equal(phase9a_normalize(phase9a_run("annotation", config, TRUE)), reference$raw, tolerance=tol,
                     info = paste(name, "raw"))
    expect_equal(phase9a_normalize(phase9a_run("annotation", config, FALSE)), reference$fit, tolerance=tol,
                     info = paste(name, "formatted"))
  }
})

test_that("Phase 9G canonical route remains reproducible", {
  comparable <- function(value) {
    value <- phase9a_normalize(value)
    value$input$ncores <- 0L
    value
  }
  config <- phase9a_configs$annotation$annot_learned
  config$ncores <- 1L
  one <- comparable(phase9a_run("annotation", config, FALSE))
  expect_identical(comparable(phase9a_run("annotation", config, FALSE)), one)
  config$ncores <- 2L
  two <- comparable(phase9a_run("annotation", config, FALSE))
  expect_identical(two, one)
  expect_identical(comparable(phase9a_run("annotation", config, FALSE)), two)
  config$ncores <- 1L
  expect_identical(comparable(phase9a_run("annotation", config, FALSE)), one)
  invisible(phase9a_run("group", phase9a_configs$group$group_one, FALSE))
  expect_identical(comparable(phase9a_run("annotation", config, FALSE)), one)
})

test_that("Phase 9G protects public and unrelated backends", {
  protected <- c(
    "src/st_cpg_omp_csr_prior.cpp" = "cce51072da6ddc3c18d58ab3b1f3c6df",
    "src/st_cpg_omp_csr_group.cpp" = "87e923f7f8ee6420e39d9f041263d11b",
    "NAMESPACE" = "1aae574d7dc2a324d4460e3477639f9a",
    "docs/dev/stblr_raw_schema.md" = "82ac9ba4b7d8edc6f3e16ee3a26d8466"
  )
  actual <- unname(tools::md5sum(vapply(names(protected), phase9g_path, character(1))))
  expect_identical(actual, unname(protected))
})
