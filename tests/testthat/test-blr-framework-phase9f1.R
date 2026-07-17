phase9f1_path <- function(...) {
  path <- file.path(...)
  if (file.exists(path)) path else file.path("..", "..", ...)
}

source(phase9f1_path(
  "tests", "testthat", "fixtures", "blr-phase9a-annotation-reference.R"
))

phase9f1_active <- function(lines, text) {
  sum(grepl(text, lines, fixed = TRUE) & !grepl("^\\s*//", lines))
}

test_that("Phase 9F1 learned-annotation execution block is singular", {
  source_lines <- readLines(phase9f1_path("src", "st_cpg_omp_csr_annot.cpp"), warn = FALSE)
  core_lines <- readLines(phase9f1_path(
    "src", "blr_csr_learned_annotation_bayesc_core_impl.h"
  ), warn = FALSE)
  source_text <- paste(source_lines, collapse = "\n")
  core_text <- paste(core_lines, collapse = "\n")
  all_sources <- list.files(phase9f1_path("src"), pattern = "\\.(cpp|h)$", full.names = TRUE)
  includes <- vapply(all_sources, function(path) {
    any(grepl('#include "blr_csr_learned_annotation_bayesc_core_impl.h"',
              readLines(path, warn = FALSE), fixed = TRUE))
  }, logical(1))

  expect_equal(phase9f1_active(core_lines, "for (int isort = 0; isort < m; ++isort)"), 1L)
  expect_equal(phase9f1_active(source_lines, "for (int isort = 0; isort < m; ++isort)"), 0L)
  expect_equal(basename(all_sources[includes]), "st_cpg_omp_csr_annot.cpp")
  expect_match(source_text, "SBLR_CSR_LEARNED_ANNOTATION_BAYESC_CORE_IMPL_TRANSLATION_UNIT", fixed = TRUE)
  expect_match(core_text, "#ifndef SBLR_BLR_CSR_LEARNED_ANNOTATION_BAYESC_CORE_IMPL_H", fixed = TRUE)
  expect_match(core_text, "implementation detail", fixed = TRUE)
  expect_false(grepl("use_old|use_new|old_path|new_path|fallback", source_text,
                     ignore.case = TRUE))
})

test_that("Phase 9F1 preserves learned-annotation policy statements", {
  core_text <- paste(readLines(phase9f1_path(
    "src", "blr_csr_learned_annotation_bayesc_core_impl.h"
  ), warn = FALSE), collapse = "\n")
  source_text <- paste(readLines(phase9f1_path(
    "src", "st_cpg_omp_csr_annot.cpp"
  ), warn = FALSE), collapse = "\n")

  expect_match(core_text, "make_pi_from_annotation(A, eta_pi_t, pi_t[1], pi_min, pi_max)", fixed = TRUE)
  expect_match(core_text, "make_vb_multiplier_from_annotation(", fixed = TRUE)
  expect_match(core_text, "rw_update_eta_pi(", fixed = TRUE)
  expect_match(core_text, "rw_update_eta_vb(", fixed = TRUE)
  expect_match(core_text, "(it + 1) % annot_update_every == 0", fixed = TRUE)
  expect_match(core_text, "vb_multiplier_min", fixed = TRUE)
  expect_match(core_text, "vb_multiplier_max", fixed = TRUE)
  expect_match(core_text, "eta_pi_accepted", fixed = TRUE)
  expect_match(core_text, "eta_pi_proposed", fixed = TRUE)
  expect_match(core_text, "eta_vb_accepted", fixed = TRUE)
  expect_match(core_text, "eta_vb_proposed", fixed = TRUE)
  expect_match(source_text, "stblr_csr_learned_annotation_bayesc_result_to_raw(", fixed = TRUE)
  expect_equal(length(gregexpr("for (int chain = 0; chain < nchains; ++chain)",
                              source_text, fixed = TRUE)[[1]]), 1L)
})

test_that("Phase 9F1 learned-annotation references remain exact", {
  expect_length(phase9a_configs$annotation, 3L)
  for (name in names(phase9a_configs$annotation)) {
    reference <- readRDS(phase9f1_path(
      "tests", "testthat", "fixtures", "blr_phase9a_annotation", paste0(name, ".rds")
    ))
    config <- phase9a_configs$annotation[[name]]
    expect_identical(phase9a_normalize(phase9a_run("annotation", config, TRUE)), reference$raw,
                     info = paste(name, "raw"))
    expect_identical(phase9a_normalize(phase9a_run("annotation", config, FALSE)), reference$fit,
                     info = paste(name, "formatted"))
  }
})

test_that("Phase 9F1 learned-annotation route remains reproducible", {
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

test_that("Phase 9F1 protects canonical and unrelated backends", {
  protected <- c(
    "src/st_cpg_omp_csr.cpp" = "92dafc0266d5a0e72aea000224154cef",
    "src/st_cpg_omp_csr_prior.cpp" = "cce51072da6ddc3c18d58ab3b1f3c6df",
    "src/st_cpg_omp_csr_group.cpp" = "87e923f7f8ee6420e39d9f041263d11b",
    "src/blr_csr_group_bayesc_core_impl.h" = "00d30ba51bb9fa7d4ef2e71f709b741f",
    "src/st_cpg_omp_csr_bayesr.cpp" = "0a005f9d5a19037285fd4869fdc4dcf0",
    "src/st_sbayesrc_omp_csr.cpp" = "8c1b03d8f5b93e6831ccbed856c77ead",
    "src/st_block_eigen.cpp" = "49f0a62c9fe235967a264b0f8de144a7",
    "src/st_cpg_omp_csr_scheduled.cpp" = "fdac03befb742f4f6fa7c22ccbbbc920",
    "src/mt_cpg_omp_csr.cpp" = "aec85896b5c30db3014efaeb5e3c3a96",
    "NAMESPACE" = "f5b6ee37a3972aa436357bdc8f602f4e",
    "src/RcppExports.cpp" = "8e4fd28d6658c32900792c316eb639fd",
    "R/RcppExports.R" = "67243ea7c7dd31154b2d1e2e5f9637e7",
    "docs/dev/stblr_raw_schema.md" = "82ac9ba4b7d8edc6f3e16ee3a26d8466"
  )
  actual <- unname(tools::md5sum(vapply(names(protected), phase9f1_path, character(1))))
  expect_identical(actual, unname(protected))
})
