phase17e_root <- normalizePath(file.path(testthat::test_path(), "..", ".."),
  winslash = "/", mustWork = TRUE)
phase17e_text <- function(path) paste(readLines(file.path(phase17e_root, path),
  warn = FALSE), collapse = "\n")
source(file.path(phase17e_root,
  "tests/testthat/fixtures/blr_phase17b_mt_default",
  "blr-phase17b-mt-default-reference.R"))
source(file.path(phase17e_root,
  "tests/testthat/fixtures/blr_phase17c_mt_default_corrected",
  "blr-phase17c-mt-default-corrected-reference.R"))

phase17e_reference <- function(id) readRDS(file.path(phase17e_root,
  "tests/testthat/fixtures/blr_phase17c_mt_default_corrected",
  sprintf("config-%d.rds", id)))

phase17e_public_source <- function() {
  adapter <- phase17e_text("src/mtblr.cpp")
  first <- regexpr("std::vector<std::vector<std::vector<double>>>  mtblr(",
    adapter, fixed = TRUE)[1]
  last <- regexpr("std::vector<std::vector<std::vector<double>>>  mtblr_hybrid(",
    adapter, fixed = TRUE)[1]
  substr(adapter, first, last - 1L)
}

test_that("all typed binding-neutral contracts are active", {
  types <- phase17e_text("src/blr_mt_default_types.h")
  for (name in c("MtDefaultDataView", "MtDefaultModelSpec",
      "MtDefaultCovariancePriorView", "MtDefaultExecutionSpec",
      "MtDefaultInitialState", "MtDefaultCoreResult"))
    expect_source_count(paste("struct", name), types, 1L)
  expect_source_forbidden(types, c("Rcpp", "SEXP", "RObject",
    "NumericVector", "NumericMatrix", "DataFrame", "CharacterVector"))
  for (field in c("wy", "ww", "yy", "XXvalues", "XXindices", "n",
      "models", "sets", "ssb_prior", "sse_prior"))
    expect_true(source_match_count(paste0("& ", field), types) >= 1L)
  expect_match(types, "std::vector<std::vector<double>> b;", fixed = TRUE)
  expect_match(types, "arma::mat B;", fixed = TRUE)
  expect_match(types, "arma::mat E;", fixed = TRUE)
  expect_match(types, "std::vector<double> pi;", fixed = TRUE)
})

test_that("one callable core and one thin public adapter are active", {
  core <- phase17e_text("src/blr_mt_default_core_impl.h")
  adapter <- phase17e_text("src/mtblr.cpp")
  public <- phase17e_public_source()
  expect_source_count("inline MtDefaultCoreResult run_mt_default_core(", core, 1L)
  expect_source_count("for ( int it = 0; it < execution.nit+execution.nburn; it++)",
    core, 1L)
  expect_source_count("run_mt_default_core(", public, 1L)
  expect_source_count("for ( int it = 0;", public, 0L)
  expect_source_count("std::mt19937", public, 0L)
  expect_source_forbidden(public, c("sampleBset(", "sampleB_latent(",
    "sampleBetaCPG_Mt_latent(", "samplePi(", "sampleB(", "computeG(",
    "sampleE("))
  expect_source_count('#include "blr_mt_default_core_impl.h"', adapter, 1L)
  expect_match(public, "result.resize(20);", fixed = TRUE)
  expect_match(public, "return result;", fixed = TRUE)
})

test_that("typed ownership and result vocabulary are explicit", {
  types <- phase17e_text("src/blr_mt_default_types.h")
  adapter <- phase17e_public_source()
  expect_match(adapter, "std::move(b), std::move(B), std::move(E), std::move(pi)",
    fixed = TRUE)
  expect_match(adapter, "std::move(initial_state)", fixed = TRUE)
  for (field in c("nt", "m", "nmodels", "marker_retained_count",
      "covb_retained_count", "covg_retained_count", "cove_retained_count",
      "pi_retained_count", "bm", "dm", "r", "b", "d", "order", "vbs",
      "vgs", "ves", "cvbm", "cvgm", "cvem", "B", "G", "E", "pi",
      "pis", "pistrait", "pismarker"))
    expect_true(source_match_count(field, types) >= 1L)
  expect_source_forbidden(types, c("dimnames", "schema", "R_NilValue",
    "stblr_raw", "result_to_raw", "aggregate_mt", "operator<"))
})

test_that("Phase 17D numerical statement and operation order are preserved", {
  header <- readLines(file.path(phase17e_root,
    "src/blr_mt_default_core_impl.h"), warn = FALSE)
  first <- match(" // Define local variables", header)
  guard <- match(" if (marker_retained_count <= 0.0) {", header)
  last <- guard + match(" }", header[(guard + 1L):length(header)])
  body <- header[first:last]
  expect_identical(length(body), 240L)
  normalized <- gsub("execution\\.nit", "nit", body)
  normalized <- gsub("execution\\.nburn", "nburn", normalized)
  normalized <- gsub("execution\\.nthin", "nthin", normalized)
  normalized <- gsub("execution\\.seed", "seed", normalized)
  normalized <- gsub("execution\\.updateB", "updateB", normalized)
  normalized <- gsub("execution\\.updateE", "updateE", normalized)
  normalized <- gsub("execution\\.updatePi", "updatePi", normalized)
  path <- tempfile(); on.exit(unlink(path), add = TRUE)
  writeChar(paste(normalized, collapse = "\n"), path, eos = NULL,
    useBytes = TRUE)
  expect_identical(unname(tools::md5sum(path)),
    "7e8ea9e4812ce57a701416f8896a97cc")
  expect_identical(sum(grepl(";\\s*(//.*)?$", normalized)), 92L)
})

test_that("corrected guards, counts, and RNG remain in the core", {
  core <- phase17e_text("src/blr_mt_default_core_impl.h")
  expect_source_count("std::mt19937 gen(execution.seed);", core, 1L)
  expect_source_forbidden(core, c("omp_get_thread_num", "static std::mt19937",
    "thread_local", "Rcpp", "SEXP", "R::", "arma::arma_rng", "ifstream",
    "ofstream"))
  expect_source_count("if (execution.updateB) {", core, 1L)
  expect_source_count("if(execution.updateB && method==4)", core, 1L)
  expect_source_count("it >= execution.nburn", core, 10L)
  expect_source_count("(it - execution.nburn) % execution.nthin", core, 2L)
  for (name in c("marker_retained_count", "covb_retained_count",
      "covg_retained_count", "cove_retained_count", "pi_retained_count")) {
    expect_true(source_match_count(name, core) >= 4L)
    expect_true(source_match_count(name,
      phase17e_text("src/blr_mt_default_types.h")) >= 1L)
  }
  expect_match(core, "marker_retained_count <= 0.0", fixed = TRUE)
  expect_source_forbidden(core, c("aggregate_mt", "result_to_raw",
    "GenotypeOperator", "DataOperator"))
})

test_that("Phase 17C raw and formatted references remain exact", {
  for (id in 1:3) {
    ref <- phase17e_reference(id)
    expect_equal(phase17c_mt_capture(id, FALSE), ref$raw, tolerance = 1e-12)
    expect_equal(phase17c_mt_capture(id, TRUE), ref$fit, tolerance = 1e-12)
    expect_identical(length(ref$raw), 20L)
  }
})

test_that("corrected identities and same-process reproducibility remain exact", {
  fixed <- phase17c_mt_capture(2L, TRUE)
  config <- phase17c_mt_config(2L)
  expect_identical(unname(fixed$vb), unname(config$vb))
  expect_identical(unname(fixed$ve), unname(config$ve))
  expect_true(all(fixed$covb == 0) && all(fixed$cove == 0) &&
    all(fixed$pim == 0))
  a <- phase17c_mt_capture(1L, TRUE)
  expect_equal(phase17c_mt_capture(1L, TRUE), a, tolerance = 1e-12)
  invisible(phase17c_mt_capture(2L, TRUE))
  expect_equal(phase17c_mt_capture(1L, TRUE), a, tolerance = 1e-12)
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

test_that("fresh process matches the corrected formatted reference", {
  skip_if_not(identical(Sys.getenv("SBLR_RUN_PHASE17E_FRESH"), "true"))
  skip_if_not_installed("callr")
  observed <- callr::r(function(root) {
    setwd(root)
    pkgload::load_all(".", compile = FALSE, quiet = TRUE)
    source(paste0("tests/testthat/fixtures/blr_phase17b_mt_default/",
      "blr-phase17b-mt-default-reference.R"))
    source(paste0("tests/testthat/fixtures/blr_phase17c_mt_default_corrected/",
      "blr-phase17c-mt-default-corrected-reference.R"))
    phase17c_mt_capture(1L, TRUE)
  }, list(root = phase17e_root))
  expect_equal(observed, phase17e_reference(1L)$fit, tolerance = 1e-12)
})

test_that("fixtures and protected production boundaries remain immutable", {
  fixture_files <- c(
    "tests/testthat/fixtures/blr_phase17b_mt_default/config-1.rds",
    "tests/testthat/fixtures/blr_phase17b_mt_default/config-2.rds",
    "tests/testthat/fixtures/blr_phase17b_mt_default/config-3.rds",
    "tests/testthat/fixtures/blr_phase17c_mt_default_corrected/config-1.rds",
    "tests/testthat/fixtures/blr_phase17c_mt_default_corrected/config-2.rds",
    "tests/testthat/fixtures/blr_phase17c_mt_default_corrected/config-3.rds")
  fixture_md5 <- c("2b820e1cdd9e731f1f0ffea442ef4e53",
    "d6d8abec35168a088a9accab87b3c6d0", "48fe8040d52b1d23c6e7437d632ebebf",
    "d64c29b872546006c0dfb0303403abce", "0d699c71e02348b39167bc1695b87b5e",
    "087ac6e06a68763a9ac558c5b21b8f3a")
  expect_identical(unname(tools::md5sum(file.path(phase17e_root,
    fixture_files))), fixture_md5)
  protected <- c(`R/interface_mtblr.R` = "18a0adea26a0495da597c7c59b5c2c1c",
    `R/RcppExports.R` = "9d13ea00b326c7e0cd606194d13a8bca",
    `src/RcppExports.cpp` = "b4859db0f6308fa7e38051ddcf32d245",
    NAMESPACE = "f5b6ee37a3972aa436357bdc8f602f4e",
    `src/mt_cpg.cpp` = "49a2c308b127de69cfe3bdf9df2be227",
    `src/mt_cpg_arma.cpp` = "f911293210e4a29017f64a92769ec814",
    `src/mt_cpg_omp.cpp` = "4c2e24988bd3151674be3c8982a36118",
    `src/mt_cpg_omp_csr.cpp` = "aec85896b5c30db3014efaeb5e3c3a96",
    `src/st_block_eigen.cpp` = "49f0a62c9fe235967a264b0f8de144a7",
    `src/st_block_eigen.h` = "bec3bc1e41841ab77747e34dc9818574",
    `src/st_cpg_omp_csr.cpp` = "92dafc0266d5a0e72aea000224154cef",
    `src/st_cpg_omp_csr_bayesr.cpp` = "0a005f9d5a19037285fd4869fdc4dcf0",
    `src/st_sbayesrc_omp_csr.cpp` = "8c1b03d8f5b93e6831ccbed856c77ead")
  expect_identical(unname(tools::md5sum(file.path(phase17e_root,
    names(protected)))), unname(protected))
  expect_match(phase17e_text("src/mt_cpg_omp.cpp"),
    "seed + 100000 * it + omp_get_thread_num()", fixed = TRUE)
})

test_that("CI covers ordinary and fresh Phase 17E tests", {
  fast <- phase17e_text(".github/workflows/blr-framework.yml")
  extended <- phase17e_text(".github/workflows/blr-framework-extended.yml")
  expect_match(fast,
    "blr-framework-phase(10|11|12|17b|17c|17d|17e)", fixed = TRUE)
  expect_match(extended, 'SBLR_RUN_PHASE17E_FRESH: "true"', fixed = TRUE)
  expect_false(grepl("blr_phase17e_mt_default_typed_boundary.R", fast,
    fixed = TRUE))
})
