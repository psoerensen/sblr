phase16a_root <- normalizePath(file.path(testthat::test_path(), "..", ".."),
  winslash = "/", mustWork = TRUE)
phase16a_text <- function(path) paste(readLines(file.path(phase16a_root, path),
  warn = FALSE), collapse = "\n")

test_that("every packed-BED BayesC route has one explicit disposition", {
  report <- phase16a_text("docs/dev/blr_framework_phase16a_report.md")
  expect_match(report, "scheduled multichain.*canonical public", perl = TRUE)
  expect_match(report, "scheduled single-chain.*retain as explicitly experimental", perl = TRUE)
  expect_match(report, "sparse marker scheduler.*retain as explicitly experimental", perl = TRUE)
  expect_match(report, "commented historical duplicates.*historical audit artifact", perl = TRUE)
})

test_that("canonical and experimental route reachability is explicit", {
  r <- phase16a_text("R/sparse_ld_bed_helper.R")
  exports <- phase16a_text("R/RcppExports.R")
  expect_match(r, "stblr_cpg_omp_bed_marker_scheduled_chains", fixed = TRUE)
  expect_match(r, "stblr_cpg_omp_bed_marker_scheduled", fixed = TRUE)
  expect_match(r, "stblr_cpg_omp_bed_marker_sparse", fixed = TRUE)
  expect_match(r, "Experimental lower-level choices", fixed = TRUE)
  expect_match(exports, "stblr_cpg_omp_bed_marker_scheduled <- function", fixed = TRUE)
  expect_match(exports, "stblr_cpg_omp_bed_marker_sparse <- function", fixed = TRUE)
})

test_that("experimental implementations remain distinct and noncanonical", {
  single <- phase16a_text("src/st_cpg_omp_individual_scheduled.cpp")
  sparse <- phase16a_text("src/st_cpg_omp_individual.cpp")
  canonical <- phase16a_text("src/st_cpg_omp_individual_scheduled_chains.cpp")
  expect_match(single, "stblr_cpg_omp_bed_marker_scheduled(", fixed = TRUE)
  expect_match(sparse, "stblr_cpg_omp_bed_marker_sparse(", fixed = TRUE)
  expect_match(sparse, "null_update_prob", fixed = TRUE)
  expect_match(canonical, "run_bed_scheduled_bayesc_chain", fixed = TRUE)
  expect_match(canonical, "aggregate_bed_scheduled_bayesc_results", fixed = TRUE)
  expect_false(grepl("run_bed_scheduled_bayesc_chain", single, fixed = TRUE))
  expect_false(grepl("run_bed_scheduled_bayesc_chain", sparse, fixed = TRUE))
})

test_that("experimental RNG and genotype ownership remain bounded", {
  for (path in c("src/st_cpg_omp_individual_scheduled.cpp",
                 "src/st_cpg_omp_individual.cpp")) {
    x <- phase16a_text(path)
    active <- paste(grep("^//", strsplit(x, "\n", fixed = TRUE)[[1]],
      value = TRUE, invert = TRUE), collapse = "\n")
    expect_false(grepl("static thread_local std::(normal|uniform)", active))
    expect_match(active, "PackedBedMatrix", fixed = TRUE)
  }
})

test_that("single-chain experimental reference remains exact", {
  old <- setwd(phase16a_root); on.exit(setwd(old), add = TRUE)
  source("tests/testthat/fixtures/blr-phase11b-bed-bayesc-reference.R")
  ref <- readRDS("tests/testthat/fixtures/blr_phase11b_bed_bayesc/single_1x1.rds")
  got <- phase11b_capture("single", 1L, 1L, 71L)
  expect_identical(phase11a_normalize(got$raw), phase11a_normalize(ref$raw))
  expect_identical(phase11a_normalize(got$fit), phase11a_normalize(ref$fit))
})

test_that("sparse experimental route is same-process deterministic", {
  old <- setwd(phase16a_root); on.exit(setwd(old), add = TRUE)
  source("tests/testthat/fixtures/blr-phase11a-bed-reference.R")
  x <- phase11a_fixture()
  run <- function() sblr::stblr_bed_marker(x$Glist, x$y, backend = "sparse",
    pi_init = .5, pi_prior_mean = .5, pi_prior_strength = 4,
    nit = 6L, nburn = 2L, nthin = 1L, seed = 71L, ncores = 1L,
    updateB = FALSE, updateE = FALSE, rebuild_every = 2L,
    full_sweep_every = 2L, null_update_prob = .5)
  expect_identical(phase11a_normalize(run()), phase11a_normalize(run()))
})

test_that("current support policy and schemas remain coherent", {
  matrix <- phase16a_text("docs/dev/blr_model_capability_matrix.md")
  schema <- phase16a_text("docs/dev/stblr_raw_schema.md")
  expect_match(matrix, "Explicitly experimental", fixed = TRUE)
  expect_match(matrix, "scheduled single-chain and sparse", fixed = TRUE)
  expect_match(schema, "schema `class = \"stblr_raw\"", fixed = TRUE)
  expect_true("stblr_bed_marker" %in% getNamespaceExports("sblr"))
  expect_true(all(c("bayesc", "bayesr", "bayesrc") %in%
    eval(formals(sblr::stblr_bed)$method)))
})

test_that("protected canonical numerical sources are unchanged", {
  protected <- c(
    "src/blr_bed_scheduled_bayesc_core_impl.h" = "723cee003504c1fdcd075b965cb63d83",
    "src/blr_bed_bayesr_core_impl.h" = "afe77e26d2cf2b8e3d64088221b33e14",
    "src/blr_bed_bayesrc_core_impl.h" = "82365cf3f1f5306c57b980f59b4d83d3",
    "src/st_block_eigen.cpp" = "49f0a62c9fe235967a264b0f8de144a7",
    "NAMESPACE" = "f5b6ee37a3972aa436357bdc8f602f4e")
  expect_identical(unname(tools::md5sum(file.path(phase16a_root,
    names(protected)))), unname(protected))
})
