phase9d2_path <- function(...) {
  path <- file.path(...)
  if (file.exists(path)) path else file.path("..", "..", ...)
}

source(phase9d2_path(
  "tests", "testthat", "fixtures", "blr-phase9a-annotation-reference.R"
))

phase9d2_active <- function(lines, text) {
  sum(grepl(text, lines, fixed = TRUE) & !grepl("^\\s*//", lines))
}

test_that("Phase 9D2 typed group execution boundary is singular and binding neutral", {
  source_lines <- readLines(phase9d2_path("src", "st_cpg_omp_csr_group.cpp"), warn = FALSE)
  core_lines <- readLines(phase9d2_path("src", "blr_csr_group_bayesc_core_impl.h"), warn = FALSE)
  type_lines <- readLines(phase9d2_path("src", "blr_csr_group_bayesc_types.h"), warn = FALSE)
  source_text <- paste(source_lines, collapse = "\n")
  core_text <- paste(core_lines, collapse = "\n")
  type_text <- paste(type_lines, collapse = "\n")

  expect_equal(phase9d2_active(core_lines, "run_csr_group_bayesc("), 1L)
  expect_equal(phase9d2_active(core_lines, "for (int isort = 0; isort < m; ++isort)"), 1L)
  expect_equal(phase9d2_active(source_lines, "run_csr_group_bayesc(context)"), 1L)
  expect_equal(phase9d2_active(source_lines, "for (int chain = 0; chain < nchains; ++chain)"), 1L)
  expect_equal(phase9d2_active(source_lines, "static Rcpp::List stblr_csr_group_bayesc_result_to_raw("), 1L)
  expect_equal(phase9d2_active(source_lines, "#include \"blr_csr_group_bayesc_core_impl.h\""), 1L)
  expect_match(type_text, "struct CsrGroupBayesCExecutionContext", fixed = TRUE)
  expect_match(type_text, "struct CsrGroupBayesCExecutionResult", fixed = TRUE)
  expect_match(type_text, "GroupBayesCPolicyView group_policy", fixed = TRUE)
  expect_match(type_text, "const arma::Row<int>* marker_group", fixed = TRUE)
  expect_match(type_text, "const void* ld_storage", fixed = TRUE)
  expect_match(type_text, "validate_group_policy(policy, data.marker_count, 1)", fixed = TRUE)
  expect_match(source_text, "context.marker_group=&group", fixed = TRUE)
  expect_match(source_text, "context.group_policy.zero_based_index", fixed = FALSE)
  expect_match(core_text, "const arma::Row<int>& group=*context.marker_group", fixed = TRUE)
  expect_match(core_text, "std::mt19937 gen_t", fixed = TRUE)
  expect_false(grepl("Rcpp|SEXP|RObject|NumericVector|NumericMatrix|List|Nullable|pybind11|Python.h",
                     paste(core_text, type_text), perl = TRUE))
  expect_false(grepl("use_old|use_new|old_path|new_path|fallback", source_text,
                     ignore.case = TRUE))

  cpp <- list.files(phase9d2_path("src"), pattern = "\\.(cpp|h)$", full.names = TRUE)
  includes <- sum(vapply(cpp, function(path) {
    any(grepl("#include \"blr_csr_group_bayesc_core_impl.h\"",
              readLines(path, warn = FALSE), fixed = TRUE))
  }, logical(1)))
  expect_equal(includes, 1L)
})

test_that("Phase 9D2 permanent group references remain exact", {
  expect_length(phase9a_configs$group, 3L)
  for (name in names(phase9a_configs$group)) {
    reference <- readRDS(phase9d2_path(
      "tests", "testthat", "fixtures", "blr_phase9a_group", paste0(name, ".rds")
    ))
    config <- phase9a_configs$group[[name]]
    expect_identical(phase9a_normalize(phase9a_run("group", config, TRUE)), reference$raw,
                     info = paste(name, "raw"))
    expect_identical(phase9a_normalize(phase9a_run("group", config, FALSE)), reference$fit,
                     info = paste(name, "formatted"))
  }
})

test_that("Phase 9D2 typed group core is reproducible across core order", {
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

test_that("Phase 9D2 retains group policy and unsupported-case contracts", {
  expect_identical(phase9a_configs$group$group_one$normalize, TRUE)
  expect_identical(phase9a_configs$group$group_chains$normalize, FALSE)
  expect_identical(phase9a_configs$group$group_explicit$seeds, c(41L, 42L))
  inputs <- phase9a_inputs(1L)
  expect_identical(unname(inputs$group), c("coding", "background", "coding", "background"))

  type_text <- paste(readLines(
    phase9d2_path("src", "blr_csr_annotation_bayesc_types.h"), warn = FALSE
  ), collapse = "\n")
  expect_match(type_text, "group.marker_group must be zero-based", fixed = TRUE)
  expect_match(type_text, "group policy does not permit empty groups", fixed = TRUE)
  expect_match(type_text, "group value dimensions are inconsistent", fixed = TRUE)
})
