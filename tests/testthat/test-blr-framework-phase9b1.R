phase9b1_path <- function(...) {
  path <- file.path(...)
  if (file.exists(path)) path else file.path("..", "..", ...)
}

source(phase9b1_path(
  "tests", "testthat", "fixtures", "blr-phase9a-annotation-reference.R"
))

phase9b1_active <- function(lines, text) {
  sum(grepl(text, lines, fixed = TRUE) & !grepl("^\\s*//", lines))
}

test_that("Phase 9B1 has one active fixed-prior BayesC execution block", {
  source_lines <- readLines(
    phase9b1_path("src", "st_cpg_omp_csr_prior.cpp"), warn = FALSE
  )
  core_lines <- readLines(
    phase9b1_path("src", "blr_csr_prior_bayesc_core_impl.h"), warn = FALSE
  )

  marker_loop <- "for (int isort = 0; isort < m; ++isort)"
  mcmc_loop <- "for (int it = 0; it < nit + nburn; ++it)"
  include_line <- "#include \"blr_csr_prior_bayesc_core_impl.h\""

  expect_equal(phase9b1_active(source_lines, marker_loop), 0L)
  expect_equal(phase9b1_active(core_lines, marker_loop), 1L)
  expect_equal(phase9b1_active(source_lines, mcmc_loop), 0L)
  # The second identical iteration header copies traces into the result after
  # sampling; only the first contains the unique marker traversal above.
  expect_equal(phase9b1_active(core_lines, mcmc_loop), 2L)
  expect_equal(phase9b1_active(source_lines, include_line), 1L)
  expect_match(
    paste(core_lines, collapse = "\n"),
    "SBLR_BLR_CSR_PRIOR_BAYESC_CORE_IMPL_H", fixed = TRUE
  )
})

test_that("Phase 9B1 preserves fixed-prior accesses and inline conversion", {
  source_text <- paste(readLines(
    phase9b1_path("src", "st_cpg_omp_csr_prior.cpp"), warn = FALSE
  ), collapse = "\n")
  core_text <- paste(readLines(
    phase9b1_path("src", "blr_csr_prior_bayesc_core_impl.h"), warn = FALSE
  ), collapse = "\n")

  expect_match(core_text, "use_pi_marker ? pi_marker_t(iu) : pi_t[1]", fixed = TRUE)
  expect_match(core_text, "use_vb_multiplier ? vb_multiplier_t(iu) : 1.0", fixed = TRUE)
  expect_match(core_text, "samplePi_ST_prior", fixed = TRUE)
  expect_match(core_text, "sampleB_ST_csr_prior", fixed = TRUE)
  expect_match(source_text, "static Rcpp::List cpg_prior_raw_v1(", fixed = TRUE)
  expect_match(source_text, "return cpg_prior_raw_v1(", fixed = TRUE)
  expect_false(grepl("getenv", paste(source_text, core_text), fixed = TRUE))
  expect_false(grepl("new_path", paste(source_text, core_text), fixed = TRUE))
  expect_false(grepl("legacy_path", paste(source_text, core_text), fixed = TRUE))
})

test_that("Phase 9B1 implementation header has one intended inclusion site", {
  cpp_paths <- list.files(
    phase9b1_path("src"), pattern = "\\.(cpp|h)$", full.names = TRUE
  )
  occurrences <- vapply(cpp_paths, function(path) {
    sum(grepl(
      "#include \"blr_csr_prior_bayesc_core_impl.h\"",
      readLines(path, warn = FALSE), fixed = TRUE
    ))
  }, integer(1))
  expect_equal(sum(occurrences), 1L)
  expect_equal(
    basename(cpp_paths[occurrences == 1L]), "st_cpg_omp_csr_prior.cpp"
  )
})

test_that("Phase 9B1 fixed-prior frozen references remain exact", {
  for (name in names(phase9a_configs$prior)) {
    reference <- readRDS(phase9b1_path(
      "tests", "testthat", "fixtures", "blr_phase9a_prior",
      paste0(name, ".rds")
    ))
    config <- phase9a_configs$prior[[name]]
    expect_identical(
      phase9a_normalize(phase9a_run("prior", config, TRUE)),
      reference$raw, info = paste(name, "raw")
    )
    expect_identical(
      phase9a_normalize(phase9a_run("prior", config, FALSE)),
      reference$fit, info = paste(name, "formatted")
    )
  }
})

test_that("Phase 9B1 fixed-prior calls remain reproducible", {
  comparable <- function(value) {
    value <- phase9a_normalize(value)
    value$input$ncores <- 0L
    value
  }
  config <- phase9a_configs$prior$fixed_chains

  first <- phase9a_run("prior", config, FALSE)
  second <- phase9a_run("prior", config, FALSE)
  expect_identical(phase9a_normalize(first), phase9a_normalize(second))

  config$ncores <- 1L
  one <- comparable(phase9a_run("prior", config, FALSE))
  config$ncores <- 2L
  two <- comparable(phase9a_run("prior", config, FALSE))
  expect_identical(two, one)
  expect_identical(comparable(phase9a_run("prior", config, FALSE)), two)
  config$ncores <- 1L
  expect_identical(comparable(phase9a_run("prior", config, FALSE)), one)

  invisible(phase9a_run(
    "annotation", phase9a_configs$annotation$annot_fixed, FALSE
  ))
  expect_identical(comparable(phase9a_run("prior", config, FALSE)), one)
})
