phase11d_root <- normalizePath(file.path(testthat::test_path(), "..", ".."),
  winslash = "/", mustWork = TRUE)
owd <- setwd(phase11d_root); on.exit(setwd(owd), add = TRUE)
source(file.path("tests", "testthat", "fixtures", "blr-phase11b-bed-bayesc-reference.R"))
read11d <- function(path) paste(readLines(path, warn = FALSE), collapse = "\n")
active11d <- function(x) paste(grep("^\\s*//", strsplit(x, "\n")[[1]],
  invert = TRUE, value = TRUE), collapse = "\n")
count11d <- function(pattern, x) {
  hit <- gregexpr(pattern, x, perl = TRUE)[[1]]
  if (identical(hit, -1L)) 0L else length(hit)
}

test_that("public BED BayesC has one closed typed execution path", {
  src <- active11d(read11d("src/st_cpg_omp_individual_scheduled_chains.cpp"))
  core <- active11d(read11d("src/blr_bed_scheduled_bayesc_core_impl.h"))
  agg <- active11d(read11d("src/blr_bed_scheduled_bayesc_aggregate_impl.h"))
  types <- active11d(read11d("src/blr_bed_scheduled_bayesc_types.h"))
  expect_identical(count11d("struct BedScheduledBayesCChainExecutionContext", types), 1L)
  expect_identical(count11d("struct BedScheduledBayesCChainExecutionResult", types), 1L)
  expect_identical(count11d("struct BedScheduledBayesCExecutionResult", types), 1L)
  expect_identical(count11d("run_bed_scheduled_bayesc_chain\\s*\\(", core), 1L)
  expect_identical(count11d("aggregate_bed_scheduled_bayesc_results\\s*\\(", agg), 1L)
  expect_identical(count11d("stblr_bed_scheduled_bayesc_result_to_raw\\s*\\(", src), 2L)
  expect_identical(count11d("for \\(int it = 0; it < total_it; \\+\\+it\\)", core), 1L)
  expect_identical(count11d("for \\(int job = 0; job < njobs; \\+\\+job\\)", src), 3L)
  expect_match(src, "#pragma omp parallel for num_threads(nthreads) schedule(static)",
    fixed = TRUE)
  expect_identical(count11d("BedScheduledBayesCChainRng chain_rng", core), 1L)
  expect_match(src, "aggregate_bed_scheduled_bayesc_results", fixed = TRUE)
  expect_match(src, "stblr_bed_scheduled_bayesc_result_to_raw", fixed = TRUE)
})

test_that("numerical and aggregation headers remain binding neutral", {
  core <- active11d(read11d("src/blr_bed_scheduled_bayesc_core_impl.h"))
  agg <- active11d(read11d("src/blr_bed_scheduled_bayesc_aggregate_impl.h"))
  types <- active11d(read11d("src/blr_bed_scheduled_bayesc_types.h"))
  expect_false(grepl("Rcpp|SEXP|Python|pybind11", paste(core, agg, types)))
  expect_false(grepl("mt19937|ChainRng|uniform_real_distribution|normal_distribution", agg))
  expect_false(grepl("fopen|fread|fseek|bed_files|FILE\\s*\\*", paste(core, agg)))
  expect_false(grepl("static thread_local std::(normal|uniform)", paste(core, agg)))
  expect_match(types, "const PackedGenotype& storage", fixed = TRUE)
  expect_false(grepl("std::vector<uint8_t>.*genotype|std::vector<std::uint8_t>.*genotype",
    paste(core, agg)))
})

test_that("aggregation owns numerical summaries and converter owns R schema", {
  src <- active11d(read11d("src/st_cpg_omp_individual_scheduled_chains.cpp"))
  agg <- active11d(read11d("src/blr_bed_scheduled_bayesc_aggregate_impl.h"))
  for (needle in c("bm_diff", "dm_diff", "inv_chains", "mean_total_log_cpo",
      "mean_seconds", "max_seconds")) expect_match(agg, needle, fixed = TRUE)
  for (needle in c("Rcpp::List marker", "stblr_raw_v1", "R_NilValue"))
    expect_match(src, needle, fixed = TRUE)
  expect_false(grepl("bm_diff|dm_diff|inv_chains", src))
  expect_false(grepl("Rcpp|R_NilValue|stblr_raw_v1", agg))
  expect_false(grepl("old_path|new_path|aggregation_selector|rng_selector|legacy_rng",
    paste(src, agg)))
})

test_that("Phase 11B corrected raw and formatted references remain exact", {
  specs <- list(single_1x1 = c("single", 1L, 1L),
    multichain_2x1 = c("multichain", 2L, 1L),
    multichain_2x2 = c("multichain", 2L, 2L))
  for (nm in names(specs)) {
    z <- specs[[nm]]
    ref <- readRDS(file.path("tests", "testthat", "fixtures",
      "blr_phase11b_bed_bayesc", paste0(nm, ".rds")))
    got <- phase11b_capture(z[[1]], as.integer(z[[3]]), as.integer(z[[2]]), 71L)
    expect_identical(phase11a_normalize(got$raw), phase11a_normalize(ref$raw))
    expect_identical(phase11a_normalize(got$fit), phase11a_normalize(ref$fit))
  }
})

test_that("same-process core-order and intervening fits remain exact", {
  a <- phase11a_normalize(phase11b_capture("multichain", 1L, 2L, 71L))
  expect_identical(phase11a_normalize(phase11b_capture("multichain", 1L, 2L, 71L)), a)
  invisible(phase11b_capture("multichain", 1L, 1L, 99L))
  expect_identical(phase11a_normalize(phase11b_capture("multichain", 1L, 2L, 71L)), a)
  expect_identical(phase11a_normalize(phase11b_capture("multichain", 2L, 2L, 71L)), a)
  expect_identical(phase11a_normalize(phase11b_capture("multichain", 2L, 2L, 71L)), a)
  expect_identical(phase11a_normalize(phase11b_capture("multichain", 1L, 2L, 71L)), a)
  invisible(phase11a_capture("bayesr", 1L, 1L, 91L))
  expect_identical(phase11a_normalize(phase11b_capture("multichain", 1L, 2L, 71L)), a)
  invisible(phase11a_capture("bayesrc", 1L, 1L, 92L))
  expect_identical(phase11a_normalize(phase11b_capture("multichain", 1L, 2L, 71L)), a)
})

test_that("route nonidentity and protected sources are permanent", {
  single <- phase11a_normalize(phase11b_capture("single", 1L, 1L, 71L))
  multi <- phase11a_normalize(phase11b_capture("multichain", 1L, 1L, 71L))
  expect_false(identical(single$raw$marker$bm, multi$raw$marker$bm))
  protected <- c(
    "src/st_cpg_omp_individual_scheduled.cpp" = "0d726fe3faf5deec887381c1458ab6b6",
    "src/st_cpg_omp_individual.cpp" = "667a0445503ef9f6b23dbab1e0114b4d",
    "src/blr_bed_scheduled_bayesc_rng.h" = "002468fa8afd7d0c491f61ea4324f982",
    "src/stblr_cpg_omp_bed_marker_scheduled_chains_bayesrc.cpp" = "9812450ae103f5026d1632b2bcc31e95",
    "src/st_block_eigen.cpp" = "49f0a62c9fe235967a264b0f8de144a7",
    "src/mt_cpg_omp_csr.cpp" = "aec85896b5c30db3014efaeb5e3c3a96",
    "R/RcppExports.R" = "9d13ea00b326c7e0cd606194d13a8bca",
    "NAMESPACE" = "f5b6ee37a3972aa436357bdc8f602f4e")
  expect_identical(unname(tools::md5sum(names(protected))), unname(protected))
})


test_that("canonical implementation headers are guarded and singly included", {
  core <- read11d("src/blr_bed_scheduled_bayesc_core_impl.h")
  agg <- read11d("src/blr_bed_scheduled_bayesc_aggregate_impl.h")
  expect_match(core, "#ifndef SBLR_BLR_BED_SCHEDULED_BAYESC_CORE_IMPL_H", fixed = TRUE)
  expect_match(agg, "#ifndef SBLR_BLR_BED_SCHEDULED_BAYESC_AGGREGATE_IMPL_H", fixed = TRUE)
  sources <- list.files("src", pattern = "\\.(cpp|h)$", full.names = TRUE)
  core_users <- sources[vapply(sources, function(x)
    grepl('#include "blr_bed_scheduled_bayesc_core_impl.h"', read11d(x), fixed = TRUE), logical(1))]
  agg_users <- sources[vapply(sources, function(x)
    grepl('#include "blr_bed_scheduled_bayesc_aggregate_impl.h"', read11d(x), fixed = TRUE), logical(1))]
  expect_identical(basename(core_users), "st_cpg_omp_individual_scheduled_chains.cpp")
  expect_identical(basename(agg_users), "st_cpg_omp_individual_scheduled_chains.cpp")
})

test_that("documentation marks only the public packed-BED BayesC route canonical", {
  plan <- read11d("docs/dev/blr_framework_implementation_plan.md")
  matrix <- read11d("docs/dev/blr_model_capability_matrix.md")
  expect_match(plan, "Phase 11D public scheduled packed-BED BayesC canonicalization", fixed = TRUE)
  expect_match(matrix, "Canonical corrected logical-chain RNG", fixed = TRUE)
  expect_match(matrix, "experimental single-chain route remains unchanged and noncanonical", fixed = TRUE)
})
