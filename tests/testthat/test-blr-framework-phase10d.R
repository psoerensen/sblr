source(file.path(testthat::test_path(), "fixtures",
  "blr-phase10b-scheduled-reference.R"))

phase10d_text <- function(path) paste(readLines(blr_repo_path(path),
  warn = FALSE), collapse = "\n")
phase10d_comparable <- function(x) {
  if (is.list(x) && !is.null(x$input)) x$input$ncores <- 0L
  if (is.list(x)) for (nm in intersect(names(x), c("seconds_mean", "seconds_max"))) x[[nm]][] <- 0
  x
}

test_that("canonical scheduled BayesC architecture is singular and build safe", {
  source <- phase10d_text("src/st_cpg_omp_csr_scheduled.cpp")
  core <- phase10d_text("src/blr_csr_scheduled_bayesc_core_impl.h")
  types <- phase10d_text("src/blr_csr_scheduled_bayesc_types.h")
  scheduled <- phase10d_text("src/blr_scheduled_execution_types.h")
  expect_match(core, "#ifndef SBLR_BLR_CSR_SCHEDULED_BAYESC_CORE_IMPL_H", fixed = TRUE)
  expect_match(core, "Implementation detail: included only by st_cpg_omp_csr_scheduled.cpp", fixed = TRUE)
  expect_equal(source_match_count("run_csr_scheduled_bayesc(", core, fixed = TRUE), 1L)
  expect_equal(source_match_count("for (int it = 0; it < total_it; ++it)", core, fixed = TRUE), 1L)
  expect_equal(source_match_count("ScheduledChainRng chain_rng(task_seed)", core, fixed = TRUE), 1L)
  expect_equal(source_match_count("stblr_csr_scheduled_bayesc_result_to_raw(", source,
    fixed = TRUE), 2L)
  expect_match(types, "struct CsrScheduledBayesCExecutionContext", fixed = TRUE)
  expect_match(types, "struct CsrScheduledBayesCExecutionResult", fixed = TRUE)
  expect_match(types, "const Operator& ld", fixed = TRUE)
  expect_match(types, "const ScheduledExecutionControl& scheduled", fixed = TRUE)
  expect_match(scheduled, "struct ScheduledChainRng", fixed = TRUE)
  expect_false(grepl("Rcpp|SEXP|RObject|pybind11|Python.h", paste(core, types)))
  expect_false(grepl("static thread_local|old_path|new_path|execution_selector|rng_selector|fallback",
    paste(source, core, types)))
})

test_that("implementation inclusion is restricted to the binding source", {
  files <- list.files(blr_repo_path("src"), recursive = TRUE, full.names = TRUE)
  files <- files[grepl("\\.(cpp|cc|cxx|h|hpp)$", files)]
  hits <- vapply(files, function(path) any(grepl(
    '#include "blr_csr_scheduled_bayesc_core_impl.h"',
    readLines(path, warn = FALSE), fixed = TRUE)), logical(1))
  expect_identical(basename(files[hits]), "st_cpg_omp_csr_scheduled.cpp")
})

test_that("canonical corrected raw and formatted references remain exact", {
  expect_length(phase10b_configs, 3L)
  for (nm in names(phase10b_configs)) {
    ref <- readRDS(file.path(testthat::test_path(), "fixtures",
      "blr_phase10b_scheduled_csr", paste0(nm, ".rds")))
    observed <- phase10b_run(phase10b_configs[[nm]])
    expect_equal(observed$raw, ref$raw, tolerance=1e-12, info = paste(nm, "raw"))
    expect_equal(observed$fit, ref$fit, tolerance=1e-12, info = paste(nm, "formatted"))
    expect_identical(ref$metadata$rng_ownership_version, "scheduled_chain_rng_v1")
    expect_identical(ref$metadata$reference_mode, "fresh R process")
  }
})

test_that("canonical path remains fit-local and worker independent", {
  a <- phase10b_configs$skip_two_one
  b <- a; b$seeds <- c(6101L, 6102L); b$full <- 2L; b$base <- 3L
  first <- phase10b_run(a)$fit
  expect_identical(phase10b_run(a)$fit, first)
  invisible(phase10b_run(b))
  expect_identical(phase10b_run(a)$fit, first)
  run_core <- function(k) { a$ncores <- k; phase10d_comparable(phase10b_run(a)$fit) }
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

test_that("scheduler semantics and unsupported boundaries remain explicit", {
  core <- phase10d_text("src/blr_csr_scheduled_bayesc_core_impl.h")
  for (needle in c("((it % full_sweep_every) == 0)", "for (int marker : active_list)",
      "for (int marker : candidate_list)", "const std::vector<int>& due",
      "auto wakeup_neighbors =", "adaptive_skip_length_csr_scheduled(")) {
    expect_match(core, needle, fixed = TRUE)
  }
  expect_false("model" %in% names(formals(sblr::stblr_csr)))
  expect_match(phase10d_text("src/st_cpg_omp_csr_scheduled.cpp"),
    "keep_chains is not yet supported for scheduled CSR", fixed = TRUE)
})

test_that("dense scheduled execution retains its documented nonidentity", {
  scheduled <- phase10b_run(phase10b_configs$dense_one)$fit
  ordinary <- sblr::stblr_csr(stats = phase10a_stats(), ld_prefix = phase10a_prefix(),
    scheduled = FALSE, pi_init = .35, pi_prior_mean = .35, pi_prior_strength = 3,
    updateB = FALSE, updateE = FALSE, updatePi = FALSE, nit = 8L, nburn = 2L,
    nthin = 1L, seed = 1001L, nchains = 1L, ncores = 1L, updateLDswap = FALSE)
  expect_false(identical(phase10b_normalize(scheduled), phase10b_normalize(ordinary)))
})

test_that("Phase 10D protects canonical and unrelated native backends", {
  protected <- c(
    "src/st_cpg_omp_csr_prior.cpp" = "cce51072da6ddc3c18d58ab3b1f3c6df",
    "src/st_cpg_omp_csr_group.cpp" = "87e923f7f8ee6420e39d9f041263d11b",
    "src/st_cpg_omp_csr_annot.cpp" = "59bd49f048d116d0fe61d73d79bd4693",
    "src/st_cpg_omp_individual_scheduled.cpp" = "0d726fe3faf5deec887381c1458ab6b6",
    "NAMESPACE" = "1aae574d7dc2a324d4460e3477639f9a")
  expect_identical(unname(tools::md5sum(vapply(names(protected), blr_repo_path, character(1)))),
    unname(protected))
})
