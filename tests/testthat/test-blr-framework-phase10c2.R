phase10c2_root <- normalizePath(file.path(testthat::test_path(), "..", ".."),
  winslash = "/", mustWork = TRUE)
source(file.path(testthat::test_path(), "fixtures",
  "blr-phase10b-scheduled-reference.R"))

phase10c2_text <- function(path) {
  paste(readLines(file.path(phase10c2_root, path), warn = FALSE), collapse = "\n")
}

phase10c2_comparable <- function(x) {
  if (is.list(x) && !is.null(x$input)) x$input$ncores <- 0L
  if (is.list(x) && !is.null(x$meta)) x$meta$ncores <- 0L
  if (is.list(x)) {
    for (nm in intersect(names(x), c("seconds_mean", "seconds_max"))) x[[nm]][] <- 0
  }
  x
}

test_that("Phase 10C2 typed scheduled boundary is singular and binding neutral", {
  source <- phase10c2_text("src/st_cpg_omp_csr_scheduled.cpp")
  core <- phase10c2_text("src/blr_csr_scheduled_bayesc_core_impl.h")
  types <- phase10c2_text("src/blr_csr_scheduled_bayesc_types.h")
  scheduled <- phase10c2_text("src/blr_scheduled_execution_types.h")

  expect_equal(length(gregexpr("run_csr_scheduled_bayesc(", core, fixed = TRUE)[[1]]), 1L)
  expect_match(source, "run_csr_scheduled_bayesc(execution_context)", fixed = TRUE)
  expect_match(types, "struct CsrScheduledBayesCExecutionContext", fixed = TRUE)
  expect_match(types, "struct CsrScheduledBayesCExecutionResult", fixed = TRUE)
  expect_match(types, "const Operator& ld", fixed = TRUE)
  expect_match(types, "const ScheduledExecutionControl& scheduled", fixed = TRUE)
  expect_match(scheduled, "struct ScheduledSweepControl", fixed = TRUE)
  expect_match(scheduled, "struct NullSkipControl", fixed = TRUE)
  expect_match(scheduled, "struct CandidateControl", fixed = TRUE)
  expect_match(scheduled, "struct NeighborWakeupControl", fixed = TRUE)
  expect_false(grepl("Rcpp|SEXP|RObject|NumericVector|NumericMatrix|pybind11|Python.h",
    paste(core, types)))
  expect_false(grepl("old_path|new_path|route_selector|execution_selector|fallback",
    paste(source, core, types)))
})

test_that("scheduler loop, traversal, and chain RNG remain in one callable core", {
  source <- phase10c2_text("src/st_cpg_omp_csr_scheduled.cpp")
  core <- phase10c2_text("src/blr_csr_scheduled_bayesc_core_impl.h")
  native <- paste(source, core, sep = "\n")

  expect_equal(length(gregexpr("for (int it = 0; it < total_it; ++it)", core,
    fixed = TRUE)[[1]]), 1L)
  expect_equal(length(gregexpr("ScheduledChainRng chain_rng(task_seed)", core,
    fixed = TRUE)[[1]]), 1L)
  expect_false(grepl("static thread_local", native, fixed = TRUE))
  expect_false(grepl("thread_local", native, fixed = TRUE))
  expect_match(core, "#pragma omp parallel for num_threads(nthreads) schedule(static)",
    fixed = TRUE)
  expect_match(core, "for (int marker : active_list)", fixed = TRUE)
  expect_match(core, "for (int marker : candidate_list)", fixed = TRUE)
  expect_match(core, "const std::vector<int>& due", fixed = TRUE)
  expect_match(core, "auto wakeup_neighbors =", fixed = TRUE)
  expect_match(core, "last_updated[static_cast<std::size_t>(marker)] == it",
    fixed = TRUE)
  expect_false(grepl("ScheduledChainRng chain_rng", source, fixed = TRUE))
})

test_that("borrowed inputs and chain-owned mutable state are explicit", {
  types <- phase10c2_text("src/blr_csr_scheduled_bayesc_types.h")
  core <- phase10c2_text("src/blr_csr_scheduled_bayesc_core_impl.h")

  for (needle in c("const arma::mat& wy", "const arma::mat& ww",
      "const std::vector<int>& marker_order",
      "const std::vector<int>& sample_sizes",
      "const std::vector<std::vector<double>>& initial_d",
      "const std::vector<std::vector<double>>& initial_r")) {
    expect_match(types, needle, fixed = TRUE)
  }
  expect_match(core, "arma::rowvec b_t(m, arma::fill::zeros)", fixed = TRUE)
  expect_match(core, "std::vector<int> scheduled_at", fixed = TRUE)
  expect_match(core, "std::vector<int> candidate_list", fixed = TRUE)
  expect_match(core, "ScheduledChainRng chain_rng(task_seed)", fixed = TRUE)
  expect_false(grepl("context.ld =", core, fixed = TRUE))
})

test_that("implementation inclusion and binding converter boundaries are permanent", {
  source <- phase10c2_text("src/st_cpg_omp_csr_scheduled.cpp")
  core <- phase10c2_text("src/blr_csr_scheduled_bayesc_core_impl.h")
  src <- list.files(file.path(phase10c2_root, "src"), recursive = TRUE,
    full.names = TRUE)
  src <- src[grepl("\\.(cpp|cc|cxx|h|hpp)$", src)]
  hits <- vapply(src, function(path) {
    any(grepl('#include "blr_csr_scheduled_bayesc_core_impl.h"',
      readLines(path, warn = FALSE), fixed = TRUE))
  }, logical(1))

  expect_identical(basename(src[hits]), "st_cpg_omp_csr_scheduled.cpp")
  expect_match(source, "auto execution_result=sblr::core::run_csr_scheduled_bayesc",
    fixed = TRUE)
  expect_match(source, "stblr_csr_scheduled_bayesc_result_to_raw(", fixed = TRUE)
  expect_match(source, "Rcpp::List marker=Rcpp::List::create(", fixed = TRUE)
  expect_match(source, 'Rcpp::Named("schema")=Rcpp::List::create(', fixed = TRUE)
  expect_false(grepl("Rcpp::List::create", core, fixed = TRUE))
})

test_that("Phase 10B corrected references remain exact behind typed execution", {
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

test_that("typed execution preserves same-process and core-order reproducibility", {
  a <- phase10b_configs$skip_two_one
  b <- a
  b$seeds <- c(5101L, 5102L)
  b$full <- 2L
  b$base <- 3L
  first <- phase10b_run(a)$fit
  expect_identical(phase10b_run(a)$fit, first)
  invisible(phase10b_run(b))
  expect_identical(phase10b_run(a)$fit, first)

  run_core <- function(k) {
    a$ncores <- k
    phase10c2_comparable(phase10b_run(a)$fit)
  }
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

test_that("fresh-process artifacts and dense-reduction classification remain unchanged", {
  for (nm in names(phase10b_configs)) {
    ref <- readRDS(file.path(testthat::test_path(), "fixtures",
      "blr_phase10b_scheduled_csr", paste0(nm, ".rds")))
    expect_identical(ref$metadata$reference_mode, "fresh R process")
  }
  scheduled <- phase10b_run(phase10b_configs$dense_one)$fit
  ordinary <- sblr::stblr_csr(stats = phase10a_stats(), ld_prefix = phase10a_prefix(),
    scheduled = FALSE, pi_init = .35, pi_prior_mean = .35, pi_prior_strength = 3,
    updateB = FALSE, updateE = FALSE, updatePi = FALSE, nit = 8L, nburn = 2L,
    nthin = 1L, seed = 1001L, nchains = 1L, ncores = 1L, updateLDswap = FALSE)
  expect_false(identical(phase10b_normalize(scheduled), phase10b_normalize(ordinary)))
})

test_that("Phase 10C2 protects canonical and unrelated native backends", {
  protected <- c(
    "src/st_cpg_omp_csr.cpp" = "92dafc0266d5a0e72aea000224154cef",
    "src/st_cpg_omp_csr_bayesr.cpp" = "0a005f9d5a19037285fd4869fdc4dcf0",
    "src/st_sbayesrc_omp_csr.cpp" = "8c1b03d8f5b93e6831ccbed856c77ead",
    "src/st_cpg_omp_csr_prior.cpp" = "cce51072da6ddc3c18d58ab3b1f3c6df",
    "src/st_cpg_omp_csr_group.cpp" = "87e923f7f8ee6420e39d9f041263d11b",
    "src/st_cpg_omp_csr_annot.cpp" = "59bd49f048d116d0fe61d73d79bd4693",
    "src/st_cpg_omp_individual_scheduled.cpp" = "ead3070c97cfcc733fe2144fe67d6b72",
    "src/st_cpg_omp_individual_scheduled_chains.cpp" = "947ecab35bd7444101c5d23d956202c0",
    "src/st_block_eigen.cpp" = "49f0a62c9fe235967a264b0f8de144a7",
    "src/mt_cpg_omp_csr.cpp" = "aec85896b5c30db3014efaeb5e3c3a96",
    "NAMESPACE" = "f5b6ee37a3972aa436357bdc8f602f4e"
  )
  actual <- unname(tools::md5sum(file.path(phase10c2_root, names(protected))))
  expect_identical(actual, unname(protected))
})
