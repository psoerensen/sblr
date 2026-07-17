phase11c2_root <- normalizePath(file.path(testthat::test_path(), "..", ".."),
  winslash = "/", mustWork = TRUE)
owd <- setwd(phase11c2_root); on.exit(setwd(owd), add = TRUE)
source(file.path("tests", "testthat", "fixtures", "blr-phase11b-bed-bayesc-reference.R"))
read11c2 <- function(path) paste(readLines(path, warn = FALSE), collapse = "\n")
active11c2 <- function(x) paste(grep("^\\s*//", strsplit(x, "\n")[[1]],
  invert = TRUE, value = TRUE), collapse = "\n")
count11c2 <- function(pattern, x) {
  hit <- gregexpr(pattern, x, perl = TRUE)[[1]]
  if (identical(hit, -1L)) 0L else length(hit)
}

test_that("public BED BayesC uses one typed per-chain execution boundary", {
  src <- active11c2(read11c2("src/st_cpg_omp_individual_scheduled_chains.cpp"))
  core <- active11c2(read11c2("src/blr_bed_scheduled_bayesc_core_impl.h"))
  types <- active11c2(read11c2("src/blr_bed_scheduled_bayesc_types.h"))
  expect_identical(count11c2("run_bed_scheduled_bayesc_chain\\s*\\(", core), 1L)
  expect_identical(count11c2("for \\(int it = 0; it < total_it; \\+\\+it\\)", core), 1L)
  expect_identical(count11c2("BedScheduledBayesCChainRng chain_rng", core), 1L)
  expect_identical(count11c2("struct BedScheduledBayesCChainExecutionContext", types), 1L)
  expect_identical(count11c2("struct BedScheduledBayesCChainExecutionResult", types), 1L)
  expect_identical(count11c2("struct BedPackedGenotypeView", types), 1L)
  expect_match(src, "BedScheduledBayesCChainExecutionContext<FastPackedBedMatrix> context", fixed = TRUE)
  expect_match(src, "run_bed_scheduled_bayesc_chain(context)", fixed = TRUE)
  expect_false(grepl("run_one_scheduled_bed_chain", paste(src, core)))
  expect_false(grepl("old_path|new_path|legacy_rng|rng_selector", paste(src, core)))
})

test_that("typed numerical core is binding neutral and I/O free", {
  core <- active11c2(read11c2("src/blr_bed_scheduled_bayesc_core_impl.h"))
  types <- active11c2(read11c2("src/blr_bed_scheduled_bayesc_types.h"))
  expect_false(grepl("Rcpp|SEXP|Python|pybind11", paste(core, types)))
  expect_false(grepl("fread|fseek|fopen|bed_files|std::FILE|FILE\\s*\\*", paste(core, types)))
  expect_false(grepl("static thread_local std::(normal|uniform)", core))
  expect_match(types, "const PackedGenotype& storage", fixed = TRUE)
  expect_match(types, "const std::uint8_t* packed_markers", fixed = TRUE)
  expect_false(grepl("std::vector<uint8_t>.*genotype|std::vector<std::uint8_t>.*genotype", core))
})

test_that("scheduler contracts and mutable ownership remain chain-local", {
  types <- read11c2("src/blr_bed_scheduled_bayesc_types.h")
  core <- read11c2("src/blr_bed_scheduled_bayesc_core_impl.h")
  for (needle in c("ScheduledSweepControl", "NullSkipControl", "CandidateControl"))
    expect_match(types, needle, fixed = TRUE)
  expect_false(grepl("NeighborWakeupControl", types, fixed = TRUE))
  for (needle in c("scheduled_at", "last_updated", "is_candidate",
      "candidate_list", "active_list", "BedScheduledBayesCChainRng chain_rng"))
    expect_match(core, needle, fixed = TRUE)
  expect_match(core, "active_list", fixed = TRUE)
  expect_match(core, "candidate_list", fixed = TRUE)
  expect_match(core, "const std::vector<int>& due", fixed = TRUE)
})

test_that("dispatch and decoding remain in adapter after migration closure", {
  src <- read11c2("src/st_cpg_omp_individual_scheduled_chains.cpp")
  core <- read11c2("src/blr_bed_scheduled_bayesc_core_impl.h")
  aggregate <- read11c2("src/blr_bed_scheduled_bayesc_aggregate_impl.h")
  for (needle in c("#pragma omp parallel for", "job_results", "std::fread",
      "std::fseek", "stblr_bed_scheduled_bayesc_result_to_raw"))
    expect_match(src, needle, fixed = TRUE)
  expect_false(grepl("#pragma omp parallel for|job_results|Rcpp::List|stblr_raw_v1", core))
  expect_match(aggregate, "aggregate_bed_scheduled_bayesc_results", fixed = TRUE)
  expect_false(grepl("Rcpp|SEXP|mt19937|ChainRng", aggregate))
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

test_that("unselected and protected sources remain unchanged", {
  protected <- c(
    "src/st_cpg_omp_individual_scheduled.cpp" = "0d726fe3faf5deec887381c1458ab6b6",
    "src/st_cpg_omp_individual.cpp" = "667a0445503ef9f6b23dbab1e0114b4d",
    "src/blr_bed_scheduled_bayesc_rng.h" = "002468fa8afd7d0c491f61ea4324f982",
    "src/stblr_cpg_omp_bed_marker_scheduled_chains_bayesr.cpp" = "85a5e45e03c59ce62654496a2f076fe9",
    "src/stblr_cpg_omp_bed_marker_scheduled_chains_bayesrc.cpp" = "5904c60b32165a7ae73bfc9d6c0f920c",
    "src/st_cpg_omp_csr_scheduled.cpp" = "fa2148492bdee4a5a363f7ecdf67c789",
    "src/st_block_eigen.cpp" = "49f0a62c9fe235967a264b0f8de144a7",
    "src/mt_cpg_omp_csr.cpp" = "aec85896b5c30db3014efaeb5e3c3a96",
    "R/RcppExports.R" = "9d13ea00b326c7e0cd606194d13a8bca",
    "NAMESPACE" = "f5b6ee37a3972aa436357bdc8f602f4e")
  expect_identical(unname(tools::md5sum(names(protected))), unname(protected))
})
