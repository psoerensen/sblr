phase11b_root <- normalizePath(file.path(testthat::test_path(), "..", ".."),
  winslash = "/", mustWork = TRUE)
owd <- setwd(phase11b_root); on.exit(setwd(owd), add = TRUE)
source(file.path("tests", "testthat", "fixtures", "blr-phase11b-bed-bayesc-reference.R"))
phase11b_text <- function(path) paste(readLines(path, warn = FALSE), collapse = "\n")

test_that("scheduled BED BayesC RNG state is logical-chain owned", {
  rng <- phase11b_text("src/blr_bed_scheduled_bayesc_rng.h")
  single <- phase11b_text("src/st_cpg_omp_individual_scheduled.cpp")
  multi <- phase11b_text("src/st_cpg_omp_individual_scheduled_chains.cpp")
  multi_core <- phase11b_text("src/blr_bed_scheduled_bayesc_core_impl.h")
  expect_match(rng, "struct BedScheduledBayesCChainRng", fixed = TRUE)
  expect_match(rng, "std::mt19937 engine", fixed = TRUE)
  expect_match(rng, "std::normal_distribution<double> normal", fixed = TRUE)
  expect_match(rng, "std::uniform_real_distribution<double> uniform", fixed = TRUE)
  expect_match(multi, "#include \"blr_bed_scheduled_bayesc_core_impl.h\"",
    fixed = TRUE)
  for (x in list(single, paste(multi, multi_core, sep = "\n"))) {
    active <- paste(grep("^//", strsplit(x, "\n", fixed = TRUE)[[1]],
      value = TRUE, invert = TRUE), collapse = "\n")
    expect_false(grepl("static thread_local std::(normal|uniform)", active))
    expect_match(active, "BedScheduledBayesCChainRng chain_rng", fixed = TRUE)
    expect_match(active, "rng.uniform(rng.engine)", fixed = TRUE)
    expect_match(active, "rng.normal(rng.engine)", fixed = TRUE)
  }
})

test_that("post-correction raw and formatted references are exact", {
  specs <- list(single_1x1 = c("single", 1L, 1L),
    multichain_2x1 = c("multichain", 2L, 1L),
    multichain_2x2 = c("multichain", 2L, 2L))
  for (nm in names(specs)) {
    z <- specs[[nm]]; ref <- readRDS(file.path("tests", "testthat", "fixtures",
      "blr_phase11b_bed_bayesc", paste0(nm, ".rds")))
    got <- phase11b_capture(z[[1]], as.integer(z[[3]]), as.integer(z[[2]]), 71L)
    expect_identical(ref$metadata$rng_ownership, "bed_scheduled_bayesc_chain_rng_v1")
    expect_identical(phase11a_normalize(got$raw), phase11a_normalize(ref$raw))
    expect_identical(phase11a_normalize(got$fit), phase11a_normalize(ref$fit))
  }
})

test_that("same-process and intervening-fit reproducibility is exact", {
  a <- phase11b_capture("multichain", 1L, 2L, 71L)
  expect_identical(phase11a_normalize(phase11b_capture("multichain", 1L, 2L, 71L)),
    phase11a_normalize(a))
  invisible(phase11b_capture("multichain", 1L, 1L, 99L))
  expect_identical(phase11a_normalize(phase11b_capture("multichain", 1L, 2L, 71L)),
    phase11a_normalize(a))
  invisible(phase11a_capture("bayesr", 1L, 1L, 91L))
  expect_identical(phase11a_normalize(phase11b_capture("multichain", 1L, 2L, 71L)),
    phase11a_normalize(a))
  invisible(phase11a_capture("bayesrc", 1L, 1L, 92L))
  expect_identical(phase11a_normalize(phase11b_capture("multichain", 1L, 2L, 71L)),
    phase11a_normalize(a))
})

test_that("multichain BayesC is core-order and chain-count independent", {
  a <- phase11a_normalize(phase11b_capture("multichain", 1L, 2L, 71L))
  expect_identical(phase11a_normalize(phase11b_capture("multichain", 2L, 2L, 71L)), a)
  expect_identical(phase11a_normalize(phase11b_capture("multichain", 2L, 2L, 71L)), a)
  expect_identical(phase11a_normalize(phase11b_capture("multichain", 1L, 2L, 71L)), a)
  one <- phase11a_normalize(phase11b_capture("multichain", 1L, 1L, 71L))
  invisible(phase11b_capture("multichain", 2L, 2L, 73L))
  expect_identical(phase11a_normalize(phase11b_capture("multichain", 1L, 1L, 71L)), one)
})

test_that("single and multichain-one-chain routes have documented seed nonidentity", {
  single <- phase11a_normalize(phase11b_capture("single", 1L, 1L, 71L))
  multi <- phase11a_normalize(phase11b_capture("multichain", 1L, 1L, 71L))
  expect_false(identical(single$raw$marker$bm, multi$raw$marker$bm))
})

test_that("Phase 11A references and protected sources remain unchanged", {
  for (model in c("bayesc", "bayesr", "bayesrc"))
    expect_true(file.exists(file.path("tests", "testthat", "fixtures",
      "blr_phase11a", paste0(model, ".rds"))))
  protected <- c(
    "src/stblr_cpg_omp_bed_marker_scheduled_chains_bayesrc.cpp" = "1f09f1e420c20a395c26867ad1d912d2",
    "src/st_cpg_omp_individual.cpp" = "667a0445503ef9f6b23dbab1e0114b4d",
    "NAMESPACE" = "f5b6ee37a3972aa436357bdc8f602f4e")
  expect_identical(unname(tools::md5sum(names(protected))), unname(protected))
})

test_that("fresh-process corrected references can be checked explicitly", {
  skip_if_not(identical(Sys.getenv("SBLR_RUN_PHASE11B_FRESH"), "true"))
  skip_if_not_installed("callr")
  observed <- callr::r(function(root) {
    setwd(root); pkgload::load_all(".", compile = FALSE, quiet = TRUE)
    source(file.path("tests", "testthat", "fixtures", "blr-phase11b-bed-bayesc-reference.R"))
    phase11b_capture("multichain", 2L, 2L, 71L)
  }, list(root = phase11b_root))
  ref <- readRDS(file.path("tests", "testthat", "fixtures",
    "blr_phase11b_bed_bayesc", "multichain_2x2.rds"))
  expect_identical(phase11a_normalize(observed$raw), phase11a_normalize(ref$raw))
  expect_identical(phase11a_normalize(observed$fit), phase11a_normalize(ref$fit))
})
