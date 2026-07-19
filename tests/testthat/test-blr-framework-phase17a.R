phase17a_root <- normalizePath(file.path(testthat::test_path(), "..", ".."),
  winslash = "/", mustWork = TRUE)
phase17a_text <- function(path) paste(readLines(file.path(phase17a_root, path),
  warn = FALSE), collapse = "\n")

test_that("all remaining non-packed backend routes are classified", {
  old <- setwd(phase17a_root); on.exit(setwd(old), add = TRUE)
  env <- new.env(parent = globalenv())
  sys.source("tools/audit/blr_phase17a_backend_inventory.R", env)
  expect_equal(nrow(env$inventory), 17L)
  expect_false(anyNA(env$inventory$disposition))
  expect_true(all(env$inventory$generated_wrapper))
  expect_setequal(unique(env$inventory$family),
    c("scalar CSR", "block-eigen", "multivariate"))
})

test_that("public and internal reachability is explicit", {
  report <- phase17a_text("docs/dev/blr_framework_phase17a_report.md")
  for (token in c("supported public", "internal research", "native-only",
                  "supported public legacy")) expect_match(report, token, fixed = TRUE)
  ns <- getNamespaceExports("sblr")
  expect_true(all(c("stblr_csr", "stblr_csr_annot", "sblr") %in% ns))
  expect_false("mtblr_cpg_omp_csr" %in% ns)
  expect_false(".stblr_csr_bayesc_block_eigen" %in% ns)
})

test_that("capability, maturity, RNG, and reference inventories agree", {
  plan <- phase17a_text("docs/dev/blr_framework_implementation_plan.md")
  matrix <- phase17a_text("docs/dev/blr_model_capability_matrix.md")
  for (x in list(plan, matrix)) {
    expect_match(x, "Phase 17A", fixed = TRUE)
    expect_match(x, "legacy multivariate", ignore.case = TRUE)
    expect_match(x, "block-eigen", fixed = TRUE)
  }
  report <- phase17a_text("docs/dev/blr_framework_phase17a_report.md")
  for (token in c("canonical architecture", "audited but noncanonical",
      "legacy architecture", "logical-chain safe", "worker-sensitive risk",
      "strong permanent references", "smoke test only", "no deterministic reference"))
    expect_match(report, token, fixed = TRUE)
})

test_that("block-eigen routes have an explicit experimental disposition", {
  r <- paste(phase17a_text("R/sparse_ld_bed_helper.R"),
             phase17a_text("R/stblr-csr-sbayesrc.R"))
  for (symbol in c(".stblr_csr_bayesc_block_eigen",
      ".stblr_csr_bayesr_block_eigen", ".stblr_csr_sbayesrc_block_eigen"))
    expect_match(r, symbol, fixed = TRUE)
  report <- phase17a_text("docs/dev/blr_framework_phase17a_report.md")
  expect_match(report, "retain experimental (P2)", fixed = TRUE)
  expect_match(report, "CSR plus fit-local BED-derived dense blocks", fixed = TRUE)
})

test_that("every multivariate implementation has a disposition", {
  report <- phase17a_text("docs/dev/blr_framework_phase17a_report.md")
  for (symbol in c("mtblr()", "mtblr_cpg()", "mtblr_cpg_arma()",
      "mtblr_cpg_omp()", "mtblr_eigen()", "mtblr_hybrid()",
      "mtblr_cpg_omp_csr()")) expect_match(report, symbol, fixed = TRUE)
  omp <- phase17a_text("src/mt_cpg_omp.cpp")
  expect_match(omp, "omp_get_thread_num()", fixed = TRUE)
  expect_match(omp, "seed + 100000 * it + omp_get_thread_num()", fixed = TRUE)
})

test_that("public routing and legacy schema boundary remain unchanged", {
  route <- phase17a_text("R/interface_mtblr.R")
  for (algorithm in c("default", "cpg", "cpg_arma", "cpg_omp", "eigen"))
    expect_match(route, paste0('algorithm=="', algorithm, '"'), fixed = TRUE)
  expect_match(route, 'names(fit) <- c("bm","dm","wy","r"', fixed = TRUE)
  expect_match(route, "seed <- sample.int(.Machine$integer.max, 1)", fixed = TRUE)
})

test_that("production numerical sources and public boundaries are unchanged", {
  protected <- c(
    "src/st_cpg_omp_csr.cpp"="92dafc0266d5a0e72aea000224154cef",
    "src/st_cpg_omp_csr_scheduled.cpp"="abeabf03db69e3358fb4850c0a432db2",
    "src/st_cpg_omp_csr_bayesr.cpp"="0a005f9d5a19037285fd4869fdc4dcf0",
    "src/st_sbayesrc_omp_csr.cpp"="8c1b03d8f5b93e6831ccbed856c77ead",
    "src/st_cpg_omp_csr_prior.cpp"="cce51072da6ddc3c18d58ab3b1f3c6df",
    "src/st_cpg_omp_csr_group.cpp"="87e923f7f8ee6420e39d9f041263d11b",
    "src/st_cpg_omp_csr_annot.cpp"="59bd49f048d116d0fe61d73d79bd4693",
    "src/st_block_eigen.cpp"="49f0a62c9fe235967a264b0f8de144a7",
    "src/st_block_eigen.h"="bec3bc1e41841ab77747e34dc9818574",
    "src/mtblr.cpp"="419472a9d17afbf39edfcafb98bba459",
    "src/mt_cpg.cpp"="49a2c308b127de69cfe3bdf9df2be227",
    "src/mt_cpg_arma.cpp"="f911293210e4a29017f64a92769ec814",
    "src/mt_cpg_omp.cpp"="4c2e24988bd3151674be3c8982a36118",
    "src/mt_cpg_omp_csr.cpp"="aec85896b5c30db3014efaeb5e3c3a96",
    "R/RcppExports.R"="9d13ea00b326c7e0cd606194d13a8bca",
    "src/RcppExports.cpp"="b4859db0f6308fa7e38051ddcf32d245",
    "NAMESPACE"="f5b6ee37a3972aa436357bdc8f602f4e")
  actual <- unname(tools::md5sum(file.path(phase17a_root, names(protected))))
  expect_identical(actual, unname(protected))
})
