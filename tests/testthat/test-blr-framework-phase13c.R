phase13c_root <- normalizePath(file.path(testthat::test_path(), "..", ".."),
  winslash = "/", mustWork = TRUE)
source(file.path(testthat::test_path(), "fixtures", "blr-phase13a-bed-bayesr-reference.R"))
phase13c_text <- function(path) paste(readLines(file.path(phase13c_root, path),
  warn = FALSE), collapse = "\n")

test_that("Phase 13C typed BayesR chain boundary is singular and binding-neutral", {
  types <- phase13c_text("src/blr_bed_bayesr_types.h")
  family_types <- phase13c_text("src/blr_bed_family_types.h")
  core <- phase13c_text("src/blr_bed_bayesr_core_impl.h")
  adapter <- phase13c_text("src/stblr_cpg_omp_bed_marker_scheduled_chains_bayesr.cpp")
  for (symbol in c("BedBayesRComponentSpec",
      "BedBayesRSchedulerControl", "BedBayesRChainExecutionContext",
      "BedBayesRChainExecutionResult", "BedBayesRProgressEvent"))
    expect_source_count(paste0("struct ", symbol), types, 1L)
  expect_source_count("using BedBayesRPackedGenotypeView", types, 1L)
  expect_source_count("struct BedPackedGenotypeView", family_types, 1L)
  expect_source_count("BedBayesRChainExecutionResult run_bed_bayesr_chain(", core, 1L)
  expect_source_count("for (int it = 0; it < total_it; ++it)", core, 1L)
  expect_source_count("run_bed_bayesr_chain(context)", adapter, 1L)
  expect_source_count("#pragma omp parallel for num_threads(nthreads) schedule(static)", adapter, 1L)
  expect_false(grepl("run_one_bayesr_chain|old_path|new_path|execution_selector", paste(core, adapter)))
  expect_false(grepl("Rcpp|SEXP|Python.h|pybind11", paste(types, core)))
})

test_that("Phase 13C preserves RNG, inverse-CDF, scheduler, and genotype boundaries", {
  types <- phase13c_text("src/blr_bed_bayesr_types.h")
  family_types <- phase13c_text("src/blr_bed_family_types.h")
  core <- phase13c_text("src/blr_bed_bayesr_core_impl.h")
  adapter <- phase13c_text("src/stblr_cpg_omp_bed_marker_scheduled_chains_bayesr.cpp")
  for (needle in c("std::mt19937 gen_t(chain_seed)",
      "std::uniform_real_distribution<double> runif", "std::normal_distribution<double> norm01",
      "std::uniform_int_distribution<int> jitter_dist", "if (u < cumsum)",
      "for (int marker : active_list)", "for (int marker : candidate_list)",
      "const std::vector<int>& due", "progress_events.push_back"))
    expect_match(core, needle, fixed = TRUE)
  expect_false(grepl("static thread_local|thread_local|std::discrete_distribution", core))
  expect_false(grepl("fopen|fseek|fread|ifstream|Rcpp::Rcout", core))
  expect_match(family_types, "const PackedGenotype& storage", fixed = TRUE)
  expect_match(adapter, "br_read_bed_blocked", fixed = TRUE)
  expect_match(adapter, "Rcpp::Rcout", fixed = TRUE)
  expect_match(adapter, "aggregate_bed_bayesr_results", fixed = TRUE)
  expect_match(adapter, "stblr_bed_bayesr_result_to_raw", fixed = TRUE)
})

test_that("Phase 13C preserves all Phase 13A frozen references exactly", {
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
})

test_that("Phase 13C progress is adapter-controlled and numerically inert", {
  quiet <- phase13a_normalize(phase13a_capture(1L, 1L, 71L, progress_every = 0))
  visible <- capture.output(x <- phase13a_capture(1L, 1L, 71L, progress_every = 1))
  expect_match(paste(visible, collapse = "\n"), "progress chain 0, trait 0", fixed = TRUE)
  x$fit$input$progress_every <- quiet$fit$input$progress_every
  expect_identical(phase13a_normalize(x), quiet)
})
