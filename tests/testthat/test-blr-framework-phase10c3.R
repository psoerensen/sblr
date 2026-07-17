phase10c3_root <- normalizePath(file.path(testthat::test_path(), "..", ".."),
  winslash = "/", mustWork = TRUE)
source(file.path(testthat::test_path(), "fixtures",
  "blr-phase10b-scheduled-reference.R"))

phase10c3_text <- function(path) paste(readLines(file.path(phase10c3_root, path),
  warn = FALSE), collapse = "\n")
phase10c3_comparable <- function(x) {
  if (is.list(x) && !is.null(x$input)) x$input$ncores <- 0L
  if (is.list(x)) for (nm in intersect(names(x), c("seconds_mean", "seconds_max"))) x[[nm]][] <- 0
  x
}

test_that("Phase 10C3 has one core, scheduler, converter, aggregation, and RNG path", {
  source <- phase10c3_text("src/st_cpg_omp_csr_scheduled.cpp")
  core <- phase10c3_text("src/blr_csr_scheduled_bayesc_core_impl.h")
  types <- phase10c3_text("src/blr_csr_scheduled_bayesc_types.h")
  expect_equal(length(gregexpr("run_csr_scheduled_bayesc(", core, fixed = TRUE)[[1]]), 1L)
  expect_equal(length(gregexpr("for (int it = 0; it < total_it; ++it)", core,
    fixed = TRUE)[[1]]), 1L)
  expect_equal(length(gregexpr("ScheduledChainRng chain_rng(task_seed)", core,
    fixed = TRUE)[[1]]), 1L)
  expect_equal(length(gregexpr("stblr_csr_scheduled_bayesc_result_to_raw(", source,
    fixed = TRUE)[[1]]), 2L) # definition and sole call
  expect_match(types, "struct CsrScheduledBayesCExecutionContext", fixed = TRUE)
  expect_match(types, "struct CsrScheduledBayesCExecutionResult", fixed = TRUE)
  expect_match(source, "return stblr_csr_scheduled_bayesc_result_to_raw(", fixed = TRUE)
  expect_false(grepl("static thread_local|old_path|new_path|execution_selector|rng_selector|fallback",
    paste(source, core, types)))
  expect_false(grepl("Rcpp|SEXP|pybind11|Python.h", paste(core, types)))
  expect_false(grepl("const arma::mat& bm_mat=execution_result", source, fixed = TRUE))
})

test_that("scheduled implementation header inclusion remains restricted", {
  files <- list.files(file.path(phase10c3_root, "src"), recursive = TRUE, full.names = TRUE)
  files <- files[grepl("\\.(cpp|cc|cxx|h|hpp)$", files)]
  hits <- vapply(files, function(path) any(grepl(
    '#include "blr_csr_scheduled_bayesc_core_impl.h"',
    readLines(path, warn = FALSE), fixed = TRUE)), logical(1))
  expect_identical(basename(files[hits]), "st_cpg_omp_csr_scheduled.cpp")
})

test_that("Phase 10B corrected raw and formatted references remain exact", {
  expect_length(phase10b_configs, 3L)
  for (nm in names(phase10b_configs)) {
    ref <- readRDS(file.path(testthat::test_path(), "fixtures",
      "blr_phase10b_scheduled_csr", paste0(nm, ".rds")))
    observed <- phase10b_run(phase10b_configs[[nm]])
    expect_identical(observed$raw, ref$raw, info = paste(nm, "raw"))
    expect_identical(observed$fit, ref$fit, info = paste(nm, "formatted"))
    expect_identical(ref$metadata$rng_ownership_version, "scheduled_chain_rng_v1")
  }
})

test_that("converter closure preserves fit-local and worker-independent results", {
  a <- phase10b_configs$skip_two_one
  b <- a; b$seeds <- c(5101L, 5102L); b$full <- 2L; b$base <- 3L
  first <- phase10b_run(a)$fit
  expect_identical(phase10b_run(a)$fit, first)
  invisible(phase10b_run(b))
  expect_identical(phase10b_run(a)$fit, first)
  run_core <- function(k) { a$ncores <- k; phase10c3_comparable(phase10b_run(a)$fit) }
  one <- run_core(1L)
  expect_identical(run_core(2L), one)
  expect_identical(run_core(2L), one)
  expect_identical(run_core(1L), one)
  dense <- phase10b_configs$dense_one
  dense_first <- phase10b_run(dense)$fit
  invisible(sblr::stblr_csr(stats = phase10a_stats(), ld_prefix = phase10a_prefix(),
    scheduled = FALSE, pi_init = .35, pi_prior_mean = .35, pi_prior_strength = 3,
    updateB = FALSE, updateE = FALSE, updatePi = FALSE, nit = 5L, nburn = 1L,
    seed = 444L))
  expect_identical(phase10b_run(dense)$fit, dense_first)
  invisible(phase10b_run(a))
  expect_identical(phase10b_run(dense)$fit, dense_first)
})

test_that("dense scheduling remains distinct and unsupported models remain unsupported", {
  scheduled <- phase10b_run(phase10b_configs$dense_one)$fit
  ordinary <- sblr::stblr_csr(stats = phase10a_stats(), ld_prefix = phase10a_prefix(),
    scheduled = FALSE, pi_init = .35, pi_prior_mean = .35, pi_prior_strength = 3,
    updateB = FALSE, updateE = FALSE, updatePi = FALSE, nit = 8L, nburn = 2L,
    nthin = 1L, seed = 1001L, nchains = 1L, ncores = 1L, updateLDswap = FALSE)
  expect_false(identical(phase10b_normalize(scheduled), phase10b_normalize(ordinary)))
  expect_false("model" %in% names(formals(sblr::stblr_csr)))
})

test_that("Phase 10C3 protects canonical and unrelated native backends", {
  protected <- c(
    "src/st_cpg_omp_csr.cpp" = "92dafc0266d5a0e72aea000224154cef",
    "src/st_cpg_omp_csr_bayesr.cpp" = "0a005f9d5a19037285fd4869fdc4dcf0",
    "src/st_sbayesrc_omp_csr.cpp" = "8c1b03d8f5b93e6831ccbed856c77ead",
    "src/st_cpg_omp_csr_prior.cpp" = "cce51072da6ddc3c18d58ab3b1f3c6df",
    "src/st_cpg_omp_csr_group.cpp" = "87e923f7f8ee6420e39d9f041263d11b",
    "src/st_cpg_omp_csr_annot.cpp" = "59bd49f048d116d0fe61d73d79bd4693",
    "src/st_cpg_omp_individual_scheduled.cpp" = "0d726fe3faf5deec887381c1458ab6b6",
    "src/st_cpg_omp_individual_scheduled_chains.cpp" = "f58fbefcffb183b9d54a96b398321dfb",
    "src/st_block_eigen.cpp" = "49f0a62c9fe235967a264b0f8de144a7",
    "src/mt_cpg_omp_csr.cpp" = "aec85896b5c30db3014efaeb5e3c3a96",
    "NAMESPACE" = "f5b6ee37a3972aa436357bdc8f602f4e")
  expect_identical(unname(tools::md5sum(file.path(phase10c3_root, names(protected)))),
    unname(protected))
})
