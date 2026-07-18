phase11c1_root <- normalizePath(file.path(testthat::test_path(), "..", ".."),
  winslash = "/", mustWork = TRUE)
owd <- setwd(phase11c1_root); on.exit(setwd(owd), add = TRUE)
source(file.path("tests", "testthat", "fixtures", "blr-phase11b-bed-bayesc-reference.R"))
read11c1 <- function(path) paste(readLines(path, warn = FALSE), collapse = "\n")
active11c1 <- function(x) paste(grep("^\\s*//", strsplit(x, "\n")[[1]],
  invert = TRUE, value = TRUE), collapse = "\n")
count11c1 <- function(pattern, x) {
  hit <- gregexpr(pattern, x, perl = TRUE)[[1]]
  if (identical(hit, -1L)) 0L else length(hit)
}

test_that("public scheduled packed-BED BayesC uses one extracted numerical body", {
  src <- active11c1(read11c1("src/st_cpg_omp_individual_scheduled_chains.cpp"))
  core <- active11c1(read11c1("src/blr_bed_scheduled_bayesc_core_impl.h"))
  public <- read11c1("R/sparse_ld_bed_helper.R")
  expect_match(public, "stblr_cpg_omp_bed_marker_scheduled_chains", fixed = TRUE)
  expect_identical(count11c1("run_bed_scheduled_bayesc_chain\\s*\\(", core), 1L)
  expect_identical(count11c1("for \\(int it = 0; it < total_it; \\+\\+it\\)", core), 1L)
  expect_match(src, '#include "blr_bed_scheduled_bayesc_core_impl.h"', fixed = TRUE)
  expect_false(grepl("for \\(int it = 0; it < total_it", src))
  expect_match(core, "BedScheduledBayesCChainRng chain_rng", fixed = TRUE)
  expect_false(grepl("static thread_local std::(normal|uniform)", core))
  expect_false(grepl("Rcpp::List|SEXP|pybind11|Python.h", core))
})

test_that("implementation header is guarded and included by one translation unit", {
  core <- read11c1("src/blr_bed_scheduled_bayesc_core_impl.h")
  expect_match(core, "#ifndef SBLR_BLR_BED_SCHEDULED_BAYESC_CORE_IMPL_H", fixed = TRUE)
  sources <- list.files("src", pattern = "\\.(cpp|h)$", full.names = TRUE)
  containing <- sources[vapply(sources, function(x)
    grepl('#include "blr_bed_scheduled_bayesc_core_impl.h"', read11c1(x), fixed = TRUE),
    logical(1))]
  expect_identical(basename(containing), "st_cpg_omp_individual_scheduled_chains.cpp")
})

test_that("BED decoding and R result construction remain in the binding source", {
  src <- read11c1("src/st_cpg_omp_individual_scheduled_chains.cpp")
  core <- read11c1("src/blr_bed_scheduled_bayesc_core_impl.h")
  for (needle in c("std::fread", "std::fseek", "Rcpp::List marker", "stblr_raw_v1"))
    expect_match(src, needle, fixed = TRUE)
  expect_false(grepl("std::fread|std::fseek|Rcpp::List marker|stblr_raw_v1", core))
})

test_that("Phase 11B corrected references remain exact after extraction", {
  specs <- list(single_1x1 = c("single", 1L, 1L),
    multichain_2x1 = c("multichain", 2L, 1L),
    multichain_2x2 = c("multichain", 2L, 2L))
  for (nm in names(specs)) {
    z <- specs[[nm]]; ref <- readRDS(file.path("tests", "testthat", "fixtures",
      "blr_phase11b_bed_bayesc", paste0(nm, ".rds")))
    got <- phase11b_capture(z[[1]], as.integer(z[[3]]), as.integer(z[[2]]), 71L)
    expect_identical(phase11a_normalize(got$raw), phase11a_normalize(ref$raw))
    expect_identical(phase11a_normalize(got$fit), phase11a_normalize(ref$fit))
  }
})

test_that("unselected and protected production sources remain unchanged", {
  protected <- c(
    "src/st_cpg_omp_individual_scheduled.cpp" = "0d726fe3faf5deec887381c1458ab6b6",
    "src/st_cpg_omp_individual.cpp" = "667a0445503ef9f6b23dbab1e0114b4d",
    "src/stblr_cpg_omp_bed_marker_scheduled_chains_bayesrc.cpp" = "9ef7d514895f80b8561de831798f2701",
    "src/st_block_eigen.cpp" = "49f0a62c9fe235967a264b0f8de144a7",
    "src/mt_cpg_omp_csr.cpp" = "aec85896b5c30db3014efaeb5e3c3a96",
    "R/RcppExports.R" = "9d13ea00b326c7e0cd606194d13a8bca",
    "NAMESPACE" = "f5b6ee37a3972aa436357bdc8f602f4e")
  expect_identical(unname(tools::md5sum(names(protected))), unname(protected))
})
