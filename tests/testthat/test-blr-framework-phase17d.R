phase17d_root <- normalizePath(file.path(testthat::test_path(), "..", ".."),
  winslash = "/", mustWork = TRUE)
phase17d_text <- function(path) paste(readLines(file.path(phase17d_root, path),
  warn = FALSE), collapse = "\n")
source(file.path(phase17d_root,
  "tests/testthat/fixtures/blr_phase17b_mt_default",
  "blr-phase17b-mt-default-reference.R"))
source(file.path(phase17d_root,
  "tests/testthat/fixtures/blr_phase17c_mt_default_corrected",
  "blr-phase17c-mt-default-corrected-reference.R"))

phase17d_reference <- function(id) readRDS(file.path(phase17d_root,
  "tests/testthat/fixtures/blr_phase17c_mt_default_corrected",
  sprintf("config-%d.rds", id)))

test_that("one guarded Phase 17D boundary has one typed successor", {
  adapter <- phase17d_text("src/mtblr.cpp")
  core <- phase17d_text("src/blr_mt_default_core_impl.h")
  expect_source_count("std::vector<std::vector<std::vector<double>>>  mtblr(",
    adapter, 1L)
  expect_source_count('#include "blr_mt_default_core_impl.h"', adapter, 1L)
  users <- list.files(file.path(phase17d_root, "src"), recursive = TRUE,
    full.names = TRUE)
  users <- users[grepl("\\.(cpp|h)$", users)]
  include_users <- users[vapply(users, function(path)
    source_match_count('#include "blr_mt_default_core_impl.h"',
      paste(readLines(path, warn = FALSE), collapse = "\n")) > 0L, logical(1))]
  expect_identical(basename(include_users), "mtblr.cpp")
  expect_match(core, "#ifndef SBLR_BLR_MT_DEFAULT_CORE_IMPL_H", fixed = TRUE)
  expect_source_count("for ( int it = 0; it < execution.nit+execution.nburn; it++)",
    core, 1L)
  public_start <- regexpr("std::vector<std::vector<std::vector<double>>>  mtblr(",
    adapter, fixed = TRUE)[1]
  hybrid_start <- regexpr("std::vector<std::vector<std::vector<double>>>  mtblr_hybrid(",
    adapter, fixed = TRUE)[1]
  public <- substr(adapter, public_start, hybrid_start - 1L)
  expect_source_count("for ( int it = 0; it < nit+nburn; it++)", public, 0L)
  expect_source_count("run_mt_default_core(", public, 1L)
  expect_match(public, "// Summarize results", fixed = TRUE)
  expect_match(public, "result.resize(20);", fixed = TRUE)
  expect_match(public, "return result;", fixed = TRUE)
})

test_that("Phase 17D mechanical body evidence remains documented", {
  header <- readLines(file.path(phase17d_root,
    "src/blr_mt_default_core_impl.h"), warn = FALSE)
  first_body <- match(" // Define local variables", header)
  guard <- match(" if (marker_retained_count <= 0.0) {", header)
  last_body <- guard + match(" }", header[(guard + 1L):length(header)])
  body <- header[first_body:last_body]
  expect_identical(length(body), 240L)
  report <- phase17d_text("docs/dev/blr_framework_phase17d_report.md")
  expect_match(report, "MECHANICAL_LINES=240", fixed = TRUE)
  expect_match(report, "IDENTICAL=TRUE", fixed = TRUE)
  adapter <- phase17d_text("src/mtblr.cpp")
  expect_match(adapter, " // Summarize results", fixed = TRUE)
  expect_match(adapter, "result.resize(20);", fixed = TRUE)
})

test_that("corrected execution contracts remain singular", {
  core <- phase17d_text("src/blr_mt_default_core_impl.h")
  expect_source_count("std::mt19937 gen(execution.seed);", core, 1L)
  expect_source_forbidden(core, c("omp_get_thread_num", "static std::mt19937",
    "thread_local", "Rcpp", "SEXP", "pybind11", "ifstream", "ofstream"))
  expect_source_count("if (execution.updateB) {", core, 1L)
  expect_source_count("sampleBset(nt, m, nub, B", core, 1L)
  expect_source_count("sampleB_latent(nt, m, nub, B", core, 1L)
  expect_source_count("if(execution.updateB && method==4)", core, 1L)
  expect_source_count("it >= execution.nburn", core, 10L)
  expect_source_count("(it - execution.nburn) % execution.nthin", core, 2L)
  for (name in c("marker_retained_count", "covb_retained_count",
      "covg_retained_count", "cove_retained_count", "pi_retained_count"))
    expect_true(source_match_count(name, core) >= 3L)
  expect_match(core, "marker_retained_count <= 0.0", fixed = TRUE)
  expect_source_forbidden(core, c("aggregate_mt", "result_to_raw"))
})

test_that("Phase 17C corrected references remain exact", {
  for (id in 1:3) {
    ref <- phase17d_reference(id)
    expect_equal(phase17c_mt_capture(id, FALSE), ref$raw, tolerance = 1e-12)
    expect_equal(phase17c_mt_capture(id, TRUE), ref$fit, tolerance = 1e-12)
    expect_identical(length(ref$raw), 20L)
  }
})

test_that("fixed controls, counts, and scientific identities remain exact", {
  fixed <- phase17c_mt_capture(2L, TRUE)
  config <- phase17c_mt_config(2L)
  expect_identical(unname(fixed$vb), unname(config$vb))
  expect_identical(unname(fixed$ve), unname(config$ve))
  expect_true(all(fixed$covb == 0) && all(fixed$cove == 0) &&
    all(fixed$pim == 0))
  for (id in c(1L, 3L)) {
    fit <- phase17c_mt_capture(id, TRUE)
    expect_equal(sum(unname(fit$pim)), 1, tolerance = 1e-12)
    expect_true(all(is.finite(unlist(fit))))
    expect_true(all(fit$d %in% 0:1) && all(fit$dm >= 0 & fit$dm <= 1))
    for (field in c("covb", "covg", "cove", "vb", "vg", "ve")) {
      expect_equal(fit[[field]], t(fit[[field]]), tolerance = 1e-12)
      expect_true(min(eigen(fit[[field]], symmetric = TRUE,
        only.values = TRUE)$values) > -1e-8)
    }
  }
})

test_that("same process and thread environment remain reproducible", {
  a <- phase17c_mt_capture(1L, TRUE)
  expect_equal(phase17c_mt_capture(1L, TRUE), a, tolerance = 1e-12)
  invisible(phase17c_mt_capture(2L, TRUE))
  expect_equal(phase17c_mt_capture(1L, TRUE), a, tolerance = 1e-12)
  old <- Sys.getenv("OMP_NUM_THREADS", unset = NA_character_)
  on.exit(if (is.na(old)) Sys.unsetenv("OMP_NUM_THREADS") else
    Sys.setenv(OMP_NUM_THREADS = old), add = TRUE)
  Sys.setenv(OMP_NUM_THREADS = "1"); one <- phase17c_mt_capture(1L, TRUE)
  Sys.setenv(OMP_NUM_THREADS = "2")
  expect_equal(phase17c_mt_capture(1L, TRUE), one, tolerance = 1e-12)
})

test_that("fresh process matches corrected formatted reference", {
  skip_if_not(identical(Sys.getenv("SBLR_RUN_PHASE17D_FRESH"), "true"))
  skip_if_not_installed("callr")
  observed <- callr::r(function(root) {
    setwd(root); pkgload::load_all(".", compile = FALSE, quiet = TRUE)
    source(paste0("tests/testthat/fixtures/blr_phase17b_mt_default/",
      "blr-phase17b-mt-default-reference.R"))
    source(paste0("tests/testthat/fixtures/blr_phase17c_mt_default_corrected/",
      "blr-phase17c-mt-default-corrected-reference.R"))
    phase17c_mt_capture(1L, TRUE)
  }, list(root = phase17d_root))
  expect_equal(observed, phase17d_reference(1L)$fit, tolerance = 1e-12)
})

test_that("fast and extended CI cover Phase 17D at the intended level", {
  fast <- phase17d_text(".github/workflows/blr-framework.yml")
  extended <- phase17d_text(".github/workflows/blr-framework-extended.yml")
  expect_match(fast, "blr-framework-phase(10|11|12|17b|17c|17d|17e)",
    fixed = TRUE)
  expect_match(extended, 'SBLR_RUN_PHASE17D_FRESH: "true"', fixed = TRUE)
  expect_false(grepl("blr_phase17d_mt_default_extraction.R", fast,
    fixed = TRUE))
})

test_that("historical and corrected fixture files remain immutable", {
  files <- c(
    "tests/testthat/fixtures/blr_phase17b_mt_default/config-1.rds",
    "tests/testthat/fixtures/blr_phase17b_mt_default/config-2.rds",
    "tests/testthat/fixtures/blr_phase17b_mt_default/config-3.rds",
    "tests/testthat/fixtures/blr_phase17c_mt_default_corrected/config-1.rds",
    "tests/testthat/fixtures/blr_phase17c_mt_default_corrected/config-2.rds",
    "tests/testthat/fixtures/blr_phase17c_mt_default_corrected/config-3.rds")
  expected <- c("2b820e1cdd9e731f1f0ffea442ef4e53",
    "d6d8abec35168a088a9accab87b3c6d0",
    "48fe8040d52b1d23c6e7437d632ebebf",
    "d64c29b872546006c0dfb0303403abce",
    "0d699c71e02348b39167bc1695b87b5e",
    "087ac6e06a68763a9ac558c5b21b8f3a")
  expect_identical(unname(tools::md5sum(file.path(phase17d_root, files))),
    expected)
})

test_that("public and protected production boundaries remain byte-identical", {
  protected <- c(
    "R/interface_mtblr.R" = "18a0adea26a0495da597c7c59b5c2c1c",
    "R/RcppExports.R" = "9d13ea00b326c7e0cd606194d13a8bca",
    "src/RcppExports.cpp" = "b4859db0f6308fa7e38051ddcf32d245",
    "NAMESPACE" = "f5b6ee37a3972aa436357bdc8f602f4e",
    "src/mt_cpg.cpp" = "49a2c308b127de69cfe3bdf9df2be227",
    "src/mt_cpg_arma.cpp" = "f911293210e4a29017f64a92769ec814",
    "src/mt_cpg_omp.cpp" = "4c2e24988bd3151674be3c8982a36118",
    "src/mt_cpg_omp_csr.cpp" = "aec85896b5c30db3014efaeb5e3c3a96",
    "src/st_block_eigen.cpp" = "49f0a62c9fe235967a264b0f8de144a7",
    "src/st_block_eigen.h" = "bec3bc1e41841ab77747e34dc9818574",
    "src/st_cpg_omp_csr.cpp" = "92dafc0266d5a0e72aea000224154cef",
    "src/st_cpg_omp_csr_bayesr.cpp" = "0a005f9d5a19037285fd4869fdc4dcf0",
    "src/st_sbayesrc_omp_csr.cpp" = "8c1b03d8f5b93e6831ccbed856c77ead")
  expect_identical(unname(tools::md5sum(file.path(phase17d_root,
    names(protected)))), unname(protected))
  expect_match(phase17d_text("src/mt_cpg_omp.cpp"),
    "seed + 100000 * it + omp_get_thread_num()", fixed = TRUE)
})
