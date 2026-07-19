phase17b_root <- normalizePath(file.path(testthat::test_path(), "..", ".."),
  winslash = "/", mustWork = TRUE)
phase17b_text <- function(path) paste(readLines(file.path(phase17b_root, path),
  warn = FALSE), collapse = "\n")
source(file.path(phase17b_root, "tests/testthat/fixtures/blr_phase17b_mt_default",
  "blr-phase17b-mt-default-reference.R"))

test_that("authoritative public route is singular and unchanged", {
  route <- phase17b_text("R/interface_mtblr.R")
  call <- 'fit <- .Call("_sblr_mtblr"'
  expect_match(route, 'if(algorithm=="default")', fixed = TRUE)
  expect_source_count(call, route, 1L)
  expect_identical(source_match_count("missing Phase 17B route", route), 0L)
  expect_identical(source_match_count(call, paste(route, route)), 2L)
  expect_match(route, "seed <- sample.int(.Machine$integer.max, 1)", fixed = TRUE)
  expect_match(route, 'names(fit) <- c("bm","dm","wy","r","b","d","o"',
    fixed = TRUE)
})

test_that("three references have exact structure and narrow numeric tolerance", {
  for (id in 1:3) {
    ref <- readRDS(file.path(phase17b_root,
      "tests/testthat/fixtures/blr_phase17b_mt_default",
      sprintf("config-%d.rds", id)))
    expect_equal(phase17b_mt_capture(id, FALSE), ref$raw, tolerance = 1e-12,
      info = paste("raw config", id))
    expect_equal(phase17b_mt_capture(id, TRUE), ref$fit, tolerance = 1e-12,
      info = paste("formatted config", id))
    expect_identical(ref$metadata$reference_mode,
      "structure_exact_numeric_tolerance")
    expect_identical(ref$metadata$numeric_tolerance, 1e-12)
    expect_true(isTRUE(ref$metadata$structure_exact))
  }
})

test_that("same-process and intervening multivariate calls are exact", {
  a1 <- phase17b_mt_capture(1, TRUE)
  a2 <- phase17b_mt_capture(1, TRUE)
  expect_equal(a1, a2, tolerance = 1e-12)
  invisible(phase17b_mt_capture(2, TRUE))
  expect_equal(a1, phase17b_mt_capture(1, TRUE), tolerance = 1e-12)
  old <- Sys.getenv("OMP_NUM_THREADS", unset = NA_character_)
  on.exit(if (is.na(old)) Sys.unsetenv("OMP_NUM_THREADS") else
    Sys.setenv(OMP_NUM_THREADS = old), add = TRUE)
  Sys.setenv(OMP_NUM_THREADS = "1")
  one <- phase17b_mt_capture(1, TRUE)
  Sys.setenv(OMP_NUM_THREADS = "2")
  two <- phase17b_mt_capture(1, TRUE)
  expect_equal(one, two, tolerance = 1e-12)
})

test_that("intervening canonical CSR and packed-BED fits have no effect", {
  a <- phase17b_mt_capture(1, TRUE)
  source(file.path(phase17b_root, "tests/testthat/fixtures/blr-phase5a-bayesr-reference.R"),
    local = environment())
  invisible(phase5a_bayesr_run(phase5a_bayesr_configs$one_chain, raw = FALSE))
  expect_equal(phase17b_mt_capture(1, TRUE), a, tolerance = 1e-12)
  source(file.path(phase17b_root, "tests/testthat/fixtures/blr-phase13a-bed-bayesr-reference.R"),
    local = environment())
  invisible(phase13a_capture(ncores = 1L, nchains = 1L, seed = 71L))
  expect_equal(phase17b_mt_capture(1, TRUE), a, tolerance = 1e-12)
})

test_that("fresh-process public default reference is exact", {
  skip_if_not(identical(Sys.getenv("SBLR_RUN_PHASE17B_FRESH"), "true"))
  skip_if_not_installed("callr")
  observed <- callr::r(function(root) {
    setwd(root); pkgload::load_all(".", compile = FALSE, quiet = TRUE)
    source("tests/testthat/fixtures/blr_phase17b_mt_default/blr-phase17b-mt-default-reference.R")
    phase17b_mt_capture(1L, TRUE)
  }, list(root = phase17b_root))
  expected <- readRDS(file.path(phase17b_root,
    "tests/testthat/fixtures/blr_phase17b_mt_default/config-1.rds"))$fit
  expect_equal(observed, expected, tolerance = 1e-12)
})

test_that("formatted scientific identities and ordering are stable", {
  for (id in 1:3) {
    fit <- phase17b_mt_capture(id, TRUE)
    expect_identical(dim(fit$bm), c(4L, 2L))
    expect_identical(rownames(fit$bm), paste0("M", 1:4))
    expect_identical(colnames(fit$bm), c("TraitA", "TraitB"))
    expect_true(all(is.finite(unlist(fit))))
    expect_true(all(fit$dm >= 0 & fit$dm <= 1))
    expect_true(all(fit$d %in% 0:1))
    for (field in c("covb", "covg", "cove", "vb", "vg", "ve")) {
      x <- fit[[field]]
      expect_equal(x, t(x), tolerance = 1e-12)
      expect_true(min(eigen(x, symmetric = TRUE, only.values = TRUE)$values) > -1e-8)
      expect_identical(rownames(x), c("TraitA", "TraitB"))
    }
    expect_identical(nrow(fit$vbs), phase17b_mt_config(id)$nit +
      phase17b_mt_config(id)$nburn)
  }
})

test_that("update-control audit records the public updateB defect", {
  fit <- phase17b_mt_capture(2, TRUE)
  config <- phase17b_mt_config(2)
  dimnames(config$vb) <- dimnames(fit$vb)
  dimnames(config$ve) <- dimnames(fit$ve)
  expect_false(isTRUE(all.equal(fit$vb, config$vb, tolerance = 0)))
  expect_identical(fit$ve, config$ve)
  expect_identical(unname(fit$pi), c(.8, rep(.2 / 3, 3)))
  report <- phase17b_text("docs/dev/blr_framework_phase17b_report.md")
  expect_match(report, "updateB = FALSE does not keep B fixed", fixed = TRUE)
})

test_that("legacy posterior denominator is frozen evidence for Phase 17C", {
  config <- phase17b_mt_config(1L)
  fit <- phase17b_mt_capture(1L, formatted = TRUE)
  iteration_index <- 0:(config$nit + config$nburn - 1L)
  n_accumulated <- sum(iteration_index > config$nburn)
  expect_identical(n_accumulated, config$nit - 1L)
  expect_equal(sum(unname(fit$pim)), n_accumulated / config$nit,
    tolerance = 1e-12)
  expect_false(isTRUE(all.equal(sum(unname(fit$pim)), 1,
    tolerance = 1e-12)))
  core <- phase17b_text("src/mtblr.cpp")
  public <- substr(core, regexpr("mtblr(", core, fixed = TRUE)[1],
    regexpr("// [[Rcpp::export]]\nstd::vector<std::vector<std::vector<double>>>  mtblr_hybrid",
      core, fixed = TRUE)[1] - 1L)
  expect_source_count("if(it>nburn) pis[k] = pis[k] + pi[k];", public, 1L)
  expect_source_count("result[17][t][i] = pis[i]/nit;", public, 1L)
  succeed("Frozen legacy evidence only; this is not the desired future policy")
})

test_that("fast and extended CI include Phase 17B at the intended level", {
  fast <- phase17b_text(".github/workflows/blr-framework.yml")
  extended <- phase17b_text(".github/workflows/blr-framework-extended.yml")
  expect_match(fast, "blr-framework-phase(10|11|12|17b)", fixed = TRUE)
  expect_match(extended, 'SBLR_RUN_PHASE17B_FRESH: "true"', fixed = TRUE)
  expect_match(extended, "devtools::test('.')", fixed = TRUE)
  expect_false(grepl("blr_phase17b_mt_default_audit.R", fast, fixed = TRUE))
  expect_false(grepl("blr_phase17b_mt_default_audit.R", extended, fixed = TRUE))
})

test_that("RNG ownership and update order are structurally frozen", {
  core <- phase17b_text("src/mtblr.cpp")
  expect_equal(length(gregexpr("std::mt19937 gen(seed);", core, fixed = TRUE)[[1]]), 3L)
  public <- substr(core, regexpr("mtblr(", core, fixed = TRUE)[1],
    regexpr("// [[Rcpp::export]]\nstd::vector<std::vector<std::vector<double>>>  mtblr_hybrid",
      core, fixed = TRUE)[1] - 1L)
  expect_match(public, "sampleBset", fixed = TRUE)
  expect_match(public, "sampleBetaCPG_Mt_latent", fixed = TRUE)
  expect_match(public, "samplePi(cmodel, pi, gen)", fixed = TRUE)
  expect_match(public, "sampleB(nt, m, nub, B", fixed = TRUE)
  expect_match(public, "computeG(nt, m", fixed = TRUE)
  expect_match(public, "sampleE(nt, m, nue, E", fixed = TRUE)
  expect_false(grepl("omp_get_thread_num", public, fixed = TRUE))
  expect_false(grepl("static std::mt19937", public, fixed = TRUE))
  expect_false(grepl("thread_local", public, fixed = TRUE))
})

test_that("worker-sensitive CPG OpenMP risk remains explicit", {
  omp <- phase17b_text("src/mt_cpg_omp.cpp")
  report <- phase17b_text("docs/dev/blr_framework_phase17b_report.md")
  expect_match(omp, "seed + 100000 * it + omp_get_thread_num()", fixed = TRUE)
  expect_match(report, "root$1$1", fixed = TRUE)
  expect_match(report, "0.0614042", fixed = TRUE)
  expect_match(report, "0.05991059", fixed = TRUE)
  expect_match(report, "correct RNG before retention", fixed = TRUE)
})

test_that("every implementation and typed audit vocabulary is classified", {
  report <- phase17b_text("docs/dev/blr_framework_phase17b_report.md")
  for (x in c("mtblr()", "mtblr_cpg()", "mtblr_cpg_arma()",
      "mtblr_cpg_omp()", "mtblr_cpg_omp_csr()", "mtblr_eigen()",
      "mtblr_hybrid()")) expect_match(report, x, fixed = TRUE)
  audit <- phase17b_text("src/blr_mt_default_audit_types.h")
  for (x in c("MtDataSpec", "MtModelSpec", "MtCovariancePriorSpec",
      "MtExecutionAuditSpec", "MtOwnershipAuditSpec",
      "MtChainResultVocabulary", "MtExecutionResultVocabulary"))
    expect_match(audit, x, fixed = TRUE)
  expect_false(grepl("Rcpp|SEXP|mt19937|#include.*mtblr", audit))
})

test_that("production and public boundaries remain byte-identical", {
  protected <- c(
    "src/mtblr.cpp"="419472a9d17afbf39edfcafb98bba459",
    "src/mt_cpg.cpp"="49a2c308b127de69cfe3bdf9df2be227",
    "src/mt_cpg_arma.cpp"="f911293210e4a29017f64a92769ec814",
    "src/mt_cpg_omp.cpp"="4c2e24988bd3151674be3c8982a36118",
    "src/mt_cpg_omp_csr.cpp"="aec85896b5c30db3014efaeb5e3c3a96",
    "src/st_block_eigen.cpp"="49f0a62c9fe235967a264b0f8de144a7",
    "src/st_block_eigen.h"="bec3bc1e41841ab77747e34dc9818574",
    "src/st_cpg_omp_csr.cpp"="92dafc0266d5a0e72aea000224154cef",
    "src/st_cpg_omp_csr_bayesr.cpp"="0a005f9d5a19037285fd4869fdc4dcf0",
    "src/st_sbayesrc_omp_csr.cpp"="8c1b03d8f5b93e6831ccbed856c77ead",
    "R/RcppExports.R"="9d13ea00b326c7e0cd606194d13a8bca",
    "src/RcppExports.cpp"="b4859db0f6308fa7e38051ddcf32d245",
    "NAMESPACE"="f5b6ee37a3972aa436357bdc8f602f4e")
  actual <- unname(tools::md5sum(file.path(phase17b_root, names(protected))))
  expect_identical(actual, unname(protected))
})
