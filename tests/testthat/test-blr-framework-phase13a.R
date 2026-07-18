phase13a_root <- normalizePath(file.path(testthat::test_path(), "..", ".."),
  winslash = "/", mustWork = TRUE)
source(file.path(testthat::test_path(), "fixtures", "blr-phase13a-bed-bayesr-reference.R"))
phase13a_text <- function(path) paste(readLines(file.path(phase13a_root, path),
  warn = FALSE), collapse = "\n")

test_that("Phase 13A BayesR route and future per-chain seam are discoverable", {
  public <- phase13a_text("R/sparse_ld_bed_helper.R")
  native <- phase13a_text("src/stblr_cpg_omp_bed_marker_scheduled_chains_bayesr.cpp")
  core <- phase13a_text("src/blr_bed_bayesr_core_impl.h")
  expect_match(public, 'method = c("bayesc", "bayesr", "bayesrc")', fixed = TRUE)
  expect_match(public, ".fit_stblr_bed_bayesr", fixed = TRUE)
  expect_source_count("stblr_cpg_omp_bed_marker_scheduled_chains_bayesr(", native, 1L)
  expect_source_count("BedBayesRChainExecutionResult run_bed_bayesr_chain(", core, 1L)
  expect_match(native, "#pragma omp parallel for num_threads(nthreads) schedule(static)", fixed = TRUE)
})

test_that("Phase 13A binding-neutral audit contracts encode validated semantics", {
  types <- phase13a_text("src/blr_bed_bayesr_audit_types.h")
  for (symbol in c("BedBayesRComponentSpec", "BedBayesRSchedulerSpec",
      "BedBayesRExecutionAuditSpec", "BedBayesROwnershipVocabulary",
      "BedBayesRChainResultVocabulary", "BedBayesRExecutionResultVocabulary"))
    expect_source_count(paste0("struct ", symbol), types, 1L)
  for (error in c("null component index must be zero", "probabilities must sum to one",
      "Dirichlet parameters must be finite and positive", "chain and core counts must be positive"))
    expect_match(types, error, fixed = TRUE)
  expect_false(grepl("Rcpp|SEXP|Python.h|pybind11", types))
})

test_that("BayesR RNG and scheduler ownership are logical-chain local", {
  native <- paste(
    phase13a_text("src/stblr_cpg_omp_bed_marker_scheduled_chains_bayesr.cpp"),
    phase13a_text("src/blr_bed_bayesr_core_impl.h"), sep = "\n")
  for (needle in c("std::mt19937 gen_t(chain_seed)",
      "std::uniform_real_distribution<double> runif(0.0, 1.0)",
      "std::normal_distribution<double> norm01(0.0, 1.0)",
      "std::uniform_int_distribution<int> jitter_dist",
      "std::gamma_distribution<double> rgamma(shape, 1.0)"))
    expect_match(native, needle, fixed = TRUE)
  expect_false(grepl("static thread_local|thread_local", native))
  for (state in c("scheduled_at", "last_updated", "candidate_list", "active_list",
      "last_interesting")) expect_match(native, state, fixed = TRUE)
  expect_match(native, "for (int marker : active_list)", fixed = TRUE)
  expect_match(native, "for (int marker : candidate_list)", fixed = TRUE)
  expect_match(native, "const std::vector<int>& due", fixed = TRUE)
})

test_that("Phase 13A frozen raw and formatted references are exact", {
  configs <- list(one_chain_one_core = c(1L, 1L, 71L),
    two_chains_one_core = c(1L, 2L, 73L),
    two_chains_two_cores = c(2L, 2L, 73L))
  for (name in names(configs)) {
    z <- configs[[name]]
    ref <- readRDS(file.path(testthat::test_path(), "fixtures",
      "blr_phase13a_bed_bayesr", paste0(name, ".rds")))
    observed <- phase13a_capture(z[1], z[2], z[3])
    expect_identical(phase13a_normalize(observed$raw), phase13a_normalize(ref$raw))
    expect_identical(phase13a_normalize(observed$fit), phase13a_normalize(ref$fit))
    expect_identical(ref$metadata$rng_ownership, "logical-chain-owned")
  }
})

test_that("BayesR component identities match current schema semantics", {
  x <- phase13a_capture(1L, 2L, 73L)
  expect_equal(unname(rowSums(x$fit$pi)), rep(1, nrow(x$fit$pi)), tolerance = 1e-12)
  expect_equal(unname(rowSums(x$fit$pim)), rep(1, nrow(x$fit$pim)), tolerance = 1e-12)
  cp <- x$fit$comp_prob[[1L]]
  expect_equal(unname(rowSums(cp)), rep(1, nrow(cp)), tolerance = 1e-12)
  expect_equal(x$fit$dm[, 1L], 1 - cp[, 1L], tolerance = 1e-12)
  expect_identical(unname(x$fit$mixture_var), c(0, .01, .1, 1))
  expect_true(all(x$fit$component[, 1L] >= 0 & x$fit$component[, 1L] <= 3))
})

test_that("BayesR reproducibility is call-order, core, and chain-count independent", {
  a <- phase13a_normalize(phase13a_capture(1L, 2L, 73L))
  expect_identical(phase13a_normalize(phase13a_capture(1L, 2L, 73L)), a)
  phase13a_capture(1L, 1L, 79L)
  expect_identical(phase13a_normalize(phase13a_capture(1L, 2L, 73L)), a)
  expect_identical(phase13a_normalize(phase13a_capture(2L, 2L, 73L)), a)
  expect_identical(phase13a_normalize(phase13a_capture(2L, 2L, 73L)), a)
  expect_identical(phase13a_normalize(phase13a_capture(1L, 2L, 73L)), a)
  phase11a_capture("bayesc", 1L, 1L, 81L)
  expect_identical(phase13a_normalize(phase13a_capture(1L, 2L, 73L)), a)
  phase11a_capture("bayesrc", 1L, 1L, 83L)
  expect_identical(phase13a_normalize(phase13a_capture(1L, 2L, 73L)), a)
})

test_that("BayesR reductions and nonreductions are explicit", {
  fixed <- phase13a_capture(1L, 1L, 71L, updatePi = FALSE)
  expect_identical(unname(fixed$fit$pi[1L, ]), c(.95, .03, .015, .005))
  expect_equal(unname(fixed$fit$pim[1L, ]), c(.95, .03, .015, .005), tolerance = 1e-15)
  dense_a <- phase13a_normalize(phase13a_capture(1L, 1L, 71L,
    full_sweep_every = 0L, null_skip_base = 50L))
  dense_b <- phase13a_normalize(phase13a_capture(1L, 1L, 71L,
    full_sweep_every = 0L, null_skip_base = 1L))
  expect_false(identical(dense_a, dense_b)) # initial scheduler jitter consumes RNG
  bayesc <- phase13a_normalize(phase11a_capture("bayesc", 1L, 1L, 71L))
  expect_false(identical(dense_a$fit$bm, bayesc$fit$bm))
})

test_that("BayesR production and protected sources remain byte-identical", {
  protected <- c(
    "src/stblr_cpg_omp_bed_marker_scheduled_chains_bayesrc.cpp" = "5904c60b32165a7ae73bfc9d6c0f920c",
    "src/st_cpg_omp_individual_scheduled_chains.cpp" = "43c71b13d8259a95f88d8a95498b213b",
    "src/st_cpg_omp_csr_bayesr.cpp" = "0a005f9d5a19037285fd4869fdc4dcf0",
    "src/st_block_eigen.cpp" = "49f0a62c9fe235967a264b0f8de144a7",
    "src/mt_cpg_omp_csr.cpp" = "aec85896b5c30db3014efaeb5e3c3a96",
    "NAMESPACE" = "f5b6ee37a3972aa436357bdc8f602f4e")
  expect_identical(unname(tools::md5sum(file.path(phase13a_root, names(protected)))),
    unname(protected))
})

test_that("Phase 13A fresh-process references are exact", {
  skip_if_not(identical(Sys.getenv("SBLR_RUN_PHASE13A_FRESH"), "true"))
  skip_if_not_installed("callr")
  for (name in c("one_chain_one_core", "two_chains_two_cores")) {
    ref <- readRDS(file.path(testthat::test_path(), "fixtures",
      "blr_phase13a_bed_bayesr", paste0(name, ".rds")))
    z <- if (name == "one_chain_one_core") c(1L, 1L, 71L) else c(2L, 2L, 73L)
    observed <- callr::r(function(root, z) {
      setwd(root); pkgload::load_all(".", compile = FALSE, quiet = TRUE)
      source(file.path("tests", "testthat", "fixtures", "blr-phase13a-bed-bayesr-reference.R"))
      phase13a_capture(z[1], z[2], z[3])
    }, list(root = phase13a_root, z = z))
    expect_identical(phase13a_normalize(observed$raw), phase13a_normalize(ref$raw))
    expect_identical(phase13a_normalize(observed$fit), phase13a_normalize(ref$fit))
  }
})
