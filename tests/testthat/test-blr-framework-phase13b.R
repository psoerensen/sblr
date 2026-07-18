phase13b_root <- normalizePath(file.path(testthat::test_path(), "..", ".."),
  winslash = "/", mustWork = TRUE)
source(file.path(testthat::test_path(), "fixtures", "blr-phase13a-bed-bayesr-reference.R"))
phase13b_text <- function(path) paste(readLines(file.path(phase13b_root, path),
  warn = FALSE), collapse = "\n")

test_that("Phase 13B mechanically isolates one BayesR chain implementation", {
  adapter <- phase13b_text("src/stblr_cpg_omp_bed_marker_scheduled_chains_bayesr.cpp")
  core <- phase13b_text("src/blr_bed_bayesr_core_impl.h")
  expect_source_count('#include "blr_bed_bayesr_core_impl.h"', adapter, 1L)
  expect_source_count("static ChainResultBayesR run_one_bayesr_chain(", core, 1L)
  expect_source_count("for (int it = 0; it < total_it; ++it)", core, 1L)
  expect_source_count("run_one_bayesr_chain(", adapter, 1L)
  expect_source_count("#pragma omp parallel for num_threads(nthreads) schedule(static)", adapter, 1L)
  expect_false(grepl("old_path|new_path|use_legacy|execution_selector", paste(adapter, core)))
})

test_that("Phase 13B header is guarded and included only by the intended source", {
  core <- phase13b_text("src/blr_bed_bayesr_core_impl.h")
  expect_match(core, "#ifndef SBLR_BLR_BED_BAYESR_CORE_IMPL_H", fixed = TRUE)
  expect_match(core, "#define SBLR_BLR_BED_BAYESR_CORE_IMPL_H", fixed = TRUE)
  expect_match(core, "Implementation detail", fixed = TRUE)
  cpp <- list.files(file.path(phase13b_root, "src"), pattern = "\\.(cpp|h)$",
    full.names = TRUE)
  consumers <- cpp[vapply(cpp, function(path) any(grepl(
    '#include "blr_bed_bayesr_core_impl.h"', readLines(path, warn = FALSE), fixed = TRUE)),
    logical(1))]
  expect_identical(basename(consumers), "stblr_cpg_omp_bed_marker_scheduled_chains_bayesr.cpp")
  expect_false(grepl("Rcpp::List::create|Rcpp::NumericMatrix|fopen|fseek|fread", core))
  expect_false(grepl("Python.h|pybind11", core))
})

test_that("Phase 13B preserves component, scheduler, and RNG mechanics", {
  core <- phase13b_text("src/blr_bed_bayesr_core_impl.h")
  for (needle in c("logp[0] = std::log", "vb * c[", "runif(gen)",
      "if (u < cumsum)", "if (d_new > 0)", "std::mt19937 gen_t(chain_seed)",
      "std::uniform_real_distribution<double> runif(0.0, 1.0)",
      "std::normal_distribution<double> norm01(0.0, 1.0)",
      "std::uniform_int_distribution<int> jitter_dist", "std::gamma_distribution<double>",
      "for (int marker : active_list)", "for (int marker : candidate_list)",
      "const std::vector<int>& due", "last_updated", "scheduled_at"))
    expect_match(core, needle, fixed = TRUE)
  expect_false(grepl("static thread_local|thread_local|std::discrete_distribution", core))
})

test_that("Phase 13B keeps decoding, dispatch, aggregation, and conversion in adapter", {
  adapter <- phase13b_text("src/stblr_cpg_omp_bed_marker_scheduled_chains_bayesr.cpp")
  core <- phase13b_text("src/blr_bed_bayesr_core_impl.h")
  for (needle in c("br_read_bed_blocked", "#pragma omp parallel for",
      "// Aggregate across chains", "// Build named raw schema v1", "Rcpp::List marker"))
    expect_match(adapter, needle, fixed = TRUE)
  expect_false(grepl("br_read_bed_blocked|Rcpp::List marker|Aggregate across chains", core))
})

test_that("Phase 13B Phase 13A references and identities remain exact", {
  configs <- list(one_chain_one_core = c(1L, 1L, 71L),
    two_chains_one_core = c(1L, 2L, 73L), two_chains_two_cores = c(2L, 2L, 73L))
  for (name in names(configs)) {
    z <- configs[[name]]
    ref <- readRDS(file.path(testthat::test_path(), "fixtures",
      "blr_phase13a_bed_bayesr", paste0(name, ".rds")))
    observed <- phase13a_capture(z[1], z[2], z[3])
    expect_identical(phase13a_normalize(observed$raw), phase13a_normalize(ref$raw))
    expect_identical(phase13a_normalize(observed$fit), phase13a_normalize(ref$fit))
  }
  x <- phase13a_capture(1L, 2L, 73L)$fit
  expect_equal(unname(rowSums(x$comp_prob[[1L]])), rep(1, nrow(x$comp_prob[[1L]])), tolerance = 1e-12)
  expect_equal(x$dm[, 1L], 1 - x$comp_prob[[1L]][, 1L], tolerance = 1e-12)
})

test_that("Phase 13B protected files remain unchanged", {
  protected <- c(
    "src/stblr_cpg_omp_bed_marker_scheduled_chains_bayesrc.cpp" = "5904c60b32165a7ae73bfc9d6c0f920c",
    "src/st_cpg_omp_individual_scheduled_chains.cpp" = "43c71b13d8259a95f88d8a95498b213b",
    "src/st_cpg_omp_csr_bayesr.cpp" = "0a005f9d5a19037285fd4869fdc4dcf0",
    "src/st_block_eigen.cpp" = "49f0a62c9fe235967a264b0f8de144a7",
    "src/mt_cpg_omp_csr.cpp" = "aec85896b5c30db3014efaeb5e3c3a96",
    "R/RcppExports.R" = "9d13ea00b326c7e0cd606194d13a8bca",
    "src/RcppExports.cpp" = "b4859db0f6308fa7e38051ddcf32d245",
    "NAMESPACE" = "f5b6ee37a3972aa436357bdc8f602f4e")
  expect_identical(unname(tools::md5sum(file.path(phase13b_root, names(protected)))),
    unname(protected))
})
