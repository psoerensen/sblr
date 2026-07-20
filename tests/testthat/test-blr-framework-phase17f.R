phase17f_root <- normalizePath(file.path(testthat::test_path(), "..", ".."),
  winslash = "/", mustWork = TRUE)
phase17f_text <- function(path) paste(readLines(file.path(phase17f_root, path),
  warn = FALSE), collapse = "\n")
source(file.path(phase17f_root,
  "tests/testthat/fixtures/blr_phase17b_mt_default",
  "blr-phase17b-mt-default-reference.R"))
source(file.path(phase17f_root,
  "tests/testthat/fixtures/blr_phase17c_mt_default_corrected",
  "blr-phase17c-mt-default-corrected-reference.R"))
phase17f_reference <- function(id) readRDS(file.path(phase17f_root,
  "tests/testthat/fixtures/blr_phase17c_mt_default_corrected",
  sprintf("config-%d.rds", id)))
phase17f_public_source <- function() {
  adapter <- phase17f_text("src/mtblr.cpp")
  first <- regexpr("std::vector<std::vector<std::vector<double>>>  mtblr(",
    adapter, fixed = TRUE)[1]
  last <- regexpr("std::vector<std::vector<std::vector<double>>>  mtblr_hybrid(",
    adapter, fixed = TRUE)[1]
  substr(adapter, first, last - 1L)
}

test_that("one typed binding-neutral finalization boundary is active", {
  types <- phase17f_text("src/blr_mt_default_types.h")
  finalizer <- phase17f_text("src/blr_mt_default_finalize_impl.h")
  adapter <- phase17f_text("src/mtblr.cpp")
  public <- phase17f_public_source()
  expect_source_count("struct MtDefaultFinalResult", types, 1L)
  expect_source_count("inline MtDefaultFinalResult finalize_mt_default_result(",
    finalizer, 1L)
  expect_source_count("finalize_mt_default_result(", public, 1L)
  expect_source_count('#include "blr_mt_default_finalize_impl.h"', adapter, 1L)
  expect_source_forbidden(paste(types, finalizer), c("Rcpp", "SEXP", "RObject",
    "NumericVector", "NumericMatrix", "dimnames", "schema_version"))
  expect_source_forbidden(finalizer, c("std::mt19937", "sampleBset(",
    "sampleB_latent(", "sampleBetaCPG_Mt_latent(", "samplePi(", "sampleB(",
    "computeG(", "sampleE(", "result.resize(20)", "Rcpp", "SEXP"))
})

test_that("all posterior arithmetic is outside the legacy adapter", {
  finalizer <- phase17f_text("src/blr_mt_default_finalize_impl.h")
  public <- phase17f_public_source()
  expect_source_count("/result.marker_retained_count", finalizer, 2L)
  expect_source_count("/ result.covb_retained_count", finalizer, 1L)
  expect_source_count("/ result.covg_retained_count", finalizer, 1L)
  expect_source_count("/ result.cove_retained_count", finalizer, 1L)
  expect_source_count("/ result.pi_retained_count", finalizer, 1L)
  for (name in c("covb_retained_count > 0.0", "covg_retained_count > 0.0",
      "cove_retained_count > 0.0", "pi_retained_count > 0.0"))
    expect_source_count(name, finalizer, 1L)
  expect_source_forbidden(public, c("/marker_retained_count",
    "covb_retained_count >", "covg_retained_count >",
    "cove_retained_count >", "pi_retained_count >", "/nit"))
})

test_that("core, final, and positional responsibilities remain singular", {
  core <- phase17f_text("src/blr_mt_default_core_impl.h")
  public <- phase17f_public_source()
  expect_source_count("for ( int it = 0; it < execution.nit+execution.nburn; it++)",
    core, 1L)
  expect_source_count("std::mt19937 gen(execution.seed);", core, 1L)
  expect_source_count("run_mt_default_core(", public, 1L)
  expect_source_count("finalize_mt_default_result(", public, 1L)
  expect_source_count("result.resize(20);", public, 1L)
  expect_source_count("for ( int it = 0;", public, 0L)
  expect_source_count("std::mt19937", public, 0L)
  expect_source_forbidden(public, c("sampleBset(", "samplePi(", "sampleE("))
  for (position in 0:19)
    expect_true(source_match_count(sprintf("result[%d]", position), public) >= 2L)
})

test_that("final result vocabulary is numerical and explicit", {
  types <- phase17f_text("src/blr_mt_default_types.h")
  for (field in c("marker_retained_count", "covb_retained_count",
      "covg_retained_count", "cove_retained_count", "pi_retained_count",
      "bm", "dm", "r", "b", "d", "marker_order", "vbs", "vgs", "ves",
      "covb", "covg", "cove", "vb", "vg", "ve", "pi_final", "pi_mean",
      "pitrait", "pimarker"))
    expect_true(source_match_count(field, types) >= 1L)
  expect_source_forbidden(types, c("marker_names", "trait_names", "class_name",
    "R_NilValue", "result_to_raw", "chain_results"))
})

test_that("finalization arithmetic and disabled policies preserve Phase 17E", {
  finalizer <- phase17f_text("src/blr_mt_default_finalize_impl.h")
  expect_identical(sum(grepl("/.*retained_count", strsplit(finalizer, "\n")[[1]])),
    6L)
  expect_source_count("for (int t=0; t < result.nt; t++)", finalizer, 1L)
  expect_source_count("for (int i=0; i < result.m; i++)", finalizer, 1L)
  expect_source_count("for (int t1=0; t1 < result.nt; t1++)", finalizer, 1L)
  expect_source_count("for (int t2=0; t2 < result.nt; t2++)", finalizer, 1L)
  expect_source_count(": 0.0;", finalizer, 4L)
  expect_match(finalizer, "result.bm[t][i]", fixed = TRUE)
  expect_match(finalizer, "result.dm[t][i]", fixed = TRUE)
})

test_that("shared naming and future operator requirements are explicit", {
  naming <- phase17f_text("docs/dev/stblr_backend_naming.md")
  inventory <- phase17f_text("docs/dev/stblr_backend_computation_inventory.md")
  expect_match(naming, "Phase 17F shared ST/MT naming plan", fixed = TRUE)
  for (term in c("sample size", "marker IDs", "trait IDs",
      "summary cross-products", "posterior marker means", "pi_final",
      "pi_mean", "pi_trace", "sample-overlap policy"))
    expect_match(naming, term, fixed = TRUE)
  expect_match(naming, "SparseLdCsrView", fixed = TRUE)
  expect_match(naming, "one view per trait/study", fixed = TRUE)
  expect_match(naming, "Fully independent CSR structures", fixed = TRUE)
  expect_match(naming, "no artificial union pattern", fixed = TRUE)
  expect_match(naming, "canonical marker ID and effect-allele orientation",
    fixed = TRUE)
  expect_match(naming, "known-overlap", fixed = TRUE)
  expect_match(inventory, "zero-based", fixed = TRUE)
  expect_match(inventory, "one operator per trait", fixed = TRUE)
  expect_match(inventory, "block-eigen", fixed = TRUE)
  expect_source_forbidden(paste(
    phase17f_text("src/blr_mt_default_types.h"),
    phase17f_text("src/blr_mt_default_finalize_impl.h")),
    c("row_ptr", "column_index", "eigenvectors", "LdOperator",
      "SparseLdCsrView", "operator_bundle"))
})

test_that("Phase 17C raw and formatted references remain exact", {
  for (id in 1:3) {
    ref <- phase17f_reference(id)
    expect_equal(phase17c_mt_capture(id, FALSE), ref$raw, tolerance = 1e-12)
    expect_equal(phase17c_mt_capture(id, TRUE), ref$fit, tolerance = 1e-12)
    expect_identical(length(ref$raw), 20L)
  }
})

test_that("corrected scientific identities and reproducibility remain exact", {
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
  skip_if_not(identical(Sys.getenv("SBLR_RUN_PHASE17F_FRESH"), "true"))
  skip_if_not_installed("callr")
  observed <- callr::r(function(root) {
    setwd(root)
    pkgload::load_all(".", compile = FALSE, quiet = TRUE)
    source(paste0("tests/testthat/fixtures/blr_phase17b_mt_default/",
      "blr-phase17b-mt-default-reference.R"))
    source(paste0("tests/testthat/fixtures/blr_phase17c_mt_default_corrected/",
      "blr-phase17c-mt-default-corrected-reference.R"))
    phase17c_mt_capture(1L, TRUE)
  }, list(root = phase17f_root))
  expect_equal(observed, phase17f_reference(1L)$fit, tolerance = 1e-12)
})

test_that("fixtures and protected production boundaries remain immutable", {
  fixtures <- c(
    "tests/testthat/fixtures/blr_phase17b_mt_default/config-1.rds",
    "tests/testthat/fixtures/blr_phase17b_mt_default/config-2.rds",
    "tests/testthat/fixtures/blr_phase17b_mt_default/config-3.rds",
    "tests/testthat/fixtures/blr_phase17c_mt_default_corrected/config-1.rds",
    "tests/testthat/fixtures/blr_phase17c_mt_default_corrected/config-2.rds",
    "tests/testthat/fixtures/blr_phase17c_mt_default_corrected/config-3.rds")
  expected <- c("2b820e1cdd9e731f1f0ffea442ef4e53",
    "d6d8abec35168a088a9accab87b3c6d0", "48fe8040d52b1d23c6e7437d632ebebf",
    "d64c29b872546006c0dfb0303403abce", "0d699c71e02348b39167bc1695b87b5e",
    "087ac6e06a68763a9ac558c5b21b8f3a")
  expect_identical(unname(tools::md5sum(file.path(phase17f_root, fixtures))),
    expected)
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
  expect_identical(unname(tools::md5sum(file.path(phase17f_root,
    names(protected)))), unname(protected))
  expect_match(phase17f_text("src/mt_cpg_omp.cpp"),
    "seed + 100000 * it + omp_get_thread_num()", fixed = TRUE)
})

test_that("CI covers ordinary and fresh Phase 17F tests", {
  fast <- phase17f_text(".github/workflows/blr-framework.yml")
  extended <- phase17f_text(".github/workflows/blr-framework-extended.yml")
  expect_match(fast,
    "blr-framework-phase(10|11|12|17b|17c|17d|17e|17f)", fixed = TRUE)
  expect_match(extended, 'SBLR_RUN_PHASE17F_FRESH: "true"', fixed = TRUE)
  expect_false(grepl("blr_phase17f_mt_default_finalization.R", fast,
    fixed = TRUE))
})
