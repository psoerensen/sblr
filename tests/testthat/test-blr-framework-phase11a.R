phase11a_root <- normalizePath(file.path(testthat::test_path(), "..", ".."),
  winslash = "/", mustWork = TRUE)
source(file.path(testthat::test_path(), "fixtures", "blr-phase11a-bed-reference.R"))
phase11a_text <- function(path) paste(readLines(file.path(phase11a_root, path),
  warn = FALSE), collapse = "\n")

test_that("Phase 11A route and model inventory remains discoverable", {
  exports <- phase11a_text("R/RcppExports.R")
  public <- phase11a_text("R/sparse_ld_bed_helper.R")
  for (symbol in c("stblr_cpg_omp_bed_marker_sparse",
      "stblr_cpg_omp_bed_marker_scheduled",
      "stblr_cpg_omp_bed_marker_scheduled_chains",
      "stblr_cpg_omp_bed_marker_scheduled_chains_bayesr",
      "stblr_cpg_omp_bed_marker_scheduled_chains_bayesrc"))
    expect_match(exports, symbol, fixed = TRUE)
  expect_match(public, "stblr_bed <- function(", fixed = TRUE)
  expect_match(public, 'method = c("bayesc", "bayesr", "bayesrc")', fixed = TRUE)
  expect_match(public, "this backend always uses unscheduled full sweeps", fixed = TRUE)
})

test_that("binding-neutral BED audit contracts encode current distinctions", {
  types <- phase11a_text("src/blr_genotype_backend_audit_types.h")
  expect_match(types, "struct PackedBedSourceView", fixed = TRUE)
  expect_match(types, "struct AdaptiveBedSchedulerControl", fixed = TRUE)
  expect_match(types, "struct BedExecutionAuditContract", fixed = TRUE)
  expect_match(types, "struct BedExecutionResultVocabulary", fixed = TRUE)
  expect_match(types, "enum class SchedulerKind", fixed = TRUE)
  expect_match(types, "full_sweep_only", fixed = TRUE)
  expect_match(types, "worker_thread_persistent", fixed = TRUE)
  for (field_error in c("bed_paths must not be empty", "packed BED must be SNP-major",
      "full_sweep_every must be non-negative", "explicit chain seeds must match chains",
      "BayesRC annotations must have one row per marker and positive columns"))
    expect_match(types, field_error, fixed = TRUE)
  expect_false(grepl("Rcpp|SEXP|pybind11|Python.h", types))
})

test_that("RNG ownership risks and safe backends are localized", {
  single <- phase11a_text("src/st_cpg_omp_individual_scheduled.cpp")
  chains <- phase11a_text("src/blr_bed_scheduled_bayesc_core_impl.h")
  rng <- phase11a_text("src/blr_bed_scheduled_bayesc_rng.h")
  bayesr <- paste(phase11a_text("src/stblr_cpg_omp_bed_marker_scheduled_chains_bayesr.cpp"),
    phase11a_text("src/blr_bed_bayesr_core_impl.h"), sep = "\n")
  bayesrc <- paste(
    phase11a_text("src/stblr_cpg_omp_bed_marker_scheduled_chains_bayesrc.cpp"),
    phase11a_text("src/blr_bed_bayesrc_core_impl.h"), sep = "\n")
  active <- paste(grep("^\\s*//", strsplit(paste(single, chains), "\\n")[[1L]],
                       invert = TRUE, value = TRUE), collapse = "\n")
  expect_false(grepl("static thread_local", active, fixed = TRUE))
  expect_match(rng, "BedScheduledBayesCChainRng", fixed = TRUE)
  expect_match(rng, "std::normal_distribution<double>", fixed = TRUE)
  expect_match(rng, "std::uniform_real_distribution<double>", fixed = TRUE)
  expect_match(bayesr, "std::normal_distribution<double> norm01(0.0, 1.0)", fixed = TRUE)
  expect_match(bayesrc, "std::normal_distribution<double> norm01(0.0, 1.0)", fixed = TRUE)
  diagnostic <- sblr:::blr_phase10a_distribution_cache_diagnostic_cpp(1701L, 2L)
  expect_true(isTRUE(diagnostic$cached_state_survives_engine_reseed))
})

test_that("fresh-process BED references are frozen by backend", {
  for (model in c("bayesc", "bayesr", "bayesrc")) {
    ref <- readRDS(file.path(testthat::test_path(), "fixtures", "blr_phase11a",
      paste0(model, ".rds")))
    expect_identical(ref$metadata$reference_mode, "fresh R process")
    expect_identical(ref$metadata$model, model)
    expect_s3_class(ref$raw, "stblr_raw_v1")
    expect_true(is.list(ref$fit))
    expect_true(all(c("bm", "dm", "input") %in% names(ref$fit)))
  }
})

test_that("safe BED BayesR and BayesRC references reproduce exactly", {
  for (model in c("bayesr", "bayesrc")) {
    ref <- readRDS(file.path(testthat::test_path(), "fixtures", "blr_phase11a",
      paste0(model, ".rds")))
    observed <- phase11a_capture(model, 1L, 1L, 71L)
    expect_identical(phase11a_normalize(observed$raw), phase11a_normalize(ref$raw))
    expect_identical(phase11a_normalize(observed$fit), phase11a_normalize(ref$fit))
  }
})

test_that("BED BayesC correction removes the Phase 11A worker assignment risk", {
  a <- phase11a_capture("bayesc", 1L, 2L, 71L)
  repeat_a <- phase11a_capture("bayesc", 1L, 2L, 71L)
  expect_identical(phase11a_normalize(repeat_a), phase11a_normalize(a))
  two_core <- phase11a_capture("bayesc", 2L, 2L, 71L)
  expect_identical(phase11a_normalize(two_core), phase11a_normalize(a))
})

test_that("safe BED models are core-order independent after metadata normalization", {
  for (model in c("bayesr", "bayesrc")) {
    nchains <- if (model == "bayesrc") 2L else 1L
    one <- phase11a_normalize(phase11a_capture(model, 1L, nchains, 71L))
    expect_identical(phase11a_normalize(phase11a_capture(model, 2L, nchains, 71L)), one)
    expect_identical(phase11a_normalize(phase11a_capture(model, 2L, nchains, 71L)), one)
    expect_identical(phase11a_normalize(phase11a_capture(model, 1L, nchains, 71L)), one)
  }
})

test_that("scheduler semantics are similar only for BayesC and BayesR", {
  chains <- phase11a_text("src/blr_bed_scheduled_bayesc_core_impl.h")
  bayesr <- paste(phase11a_text("src/stblr_cpg_omp_bed_marker_scheduled_chains_bayesr.cpp"),
    phase11a_text("src/blr_bed_bayesr_core_impl.h"), sep = "\n")
  bayesrc <- paste(
    phase11a_text("src/stblr_cpg_omp_bed_marker_scheduled_chains_bayesrc.cpp"),
    phase11a_text("src/blr_bed_bayesrc_core_impl.h"), sep = "\n")
  for (needle in c("full_sweep_every", "null_skip_base", "candidate_lifetime",
      "scheduled_at", "candidate_list")) {
    expect_match(chains, needle, fixed = TRUE)
    expect_match(bayesr, needle, fixed = TRUE)
  }
  expect_false(grepl("candidate_list|scheduled_at", bayesrc))
  expect_match(bayesrc, "Exact full sweep: every marker is visited once", fixed = TRUE)
})

test_that("Phase 11A leaves production and protected sources unchanged", {
  protected <- c(
    "src/st_cpg_omp_individual.cpp" = "667a0445503ef9f6b23dbab1e0114b4d",
    "src/st_cpg_omp_individual_scheduled.cpp" = "0d726fe3faf5deec887381c1458ab6b6",
    "src/stblr_cpg_omp_bed_marker_scheduled_chains_bayesrc.cpp" = "72d4a9fa0a7cd51071328c2d62d0192b",
    "src/st_cpg_omp_csr.cpp" = "92dafc0266d5a0e72aea000224154cef",
    "src/st_block_eigen.cpp" = "49f0a62c9fe235967a264b0f8de144a7",
    "src/mt_cpg_omp_csr.cpp" = "aec85896b5c30db3014efaeb5e3c3a96",
    "NAMESPACE" = "f5b6ee37a3972aa436357bdc8f602f4e")
  expect_identical(unname(tools::md5sum(file.path(phase11a_root, names(protected)))),
    unname(protected))
  migrated <- phase11a_text("src/st_cpg_omp_individual_scheduled_chains.cpp")
  expect_match(migrated, "run_bed_scheduled_bayesc_chain(context)", fixed = TRUE)
  expect_match(migrated, "aggregate_bed_scheduled_bayesc_results", fixed = TRUE)
  expect_match(migrated, "stblr_bed_scheduled_bayesc_result_to_raw", fixed = TRUE)
})

test_that("fresh-process references can be checked explicitly", {
  skip_if_not(identical(Sys.getenv("SBLR_RUN_PHASE11A_FRESH"), "true"))
  skip_if_not_installed("callr")
  for (model in c("bayesr", "bayesrc")) {
    observed <- callr::r(function(root, model) {
      setwd(root); pkgload::load_all(".", compile = FALSE, quiet = TRUE)
      source(file.path("tests", "testthat", "fixtures", "blr-phase11a-bed-reference.R"))
      phase11a_capture(model, 1L, if (model == "bayesc") 2L else 1L, 71L)
    }, list(root = phase11a_root, model = model))
    ref <- readRDS(file.path(testthat::test_path(), "fixtures", "blr_phase11a",
      paste0(model, ".rds")))
    expect_identical(phase11a_normalize(observed$raw), phase11a_normalize(ref$raw))
    expect_identical(phase11a_normalize(observed$fit), phase11a_normalize(ref$fit))
  }
})
