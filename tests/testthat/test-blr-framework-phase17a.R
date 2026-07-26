phase17a_text <- function(path) paste(readLines(blr_repo_path(path),
  warn = FALSE), collapse = "\n")

test_that("all remaining non-packed backend routes are classified", {
  old <- setwd(blr_repo_path()); on.exit(setwd(old), add = TRUE)
  env <- new.env(parent = globalenv())
  sys.source("tools/audit/blr_phase17a_backend_inventory.R", env)
  expect_equal(nrow(env$inventory), 13L)
  expect_false(anyNA(env$inventory$disposition))
  expect_true(all(env$inventory$generated_wrapper))
  expect_equal(nrow(env$historical_disposition), 4L)
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
  inventory <- phase17a_text("tools/audit/blr_phase17a_backend_inventory.R")
  expect_match(inventory, "historical_disposition", fixed = TRUE)
  expect_match(inventory, "retired and removed", fixed = TRUE)
})

test_that("public routing retains one supported legacy schema boundary", {
  route <- phase17a_text("R/interface_mtblr.R")
  expect_match(route, '!identical(algorithm, "default")', fixed = TRUE)
  expect_equal(source_match_count('.Call("_sblr_mtblr"', route, fixed = TRUE), 1L)
  for (algorithm in c("_sblr_mtblr_cpg", "_sblr_mtblr_cpg_arma",
      "_sblr_mtblr_cpg_omp", "_sblr_mtblr_eigen"))
    expect_false(grepl(algorithm, route, fixed = TRUE))
  expect_match(route, 'names(fit) <- c("bm","dm","wy","r"', fixed = TRUE)
  expect_match(route, "seed <- sample.int(.Machine$integer.max, 1)", fixed = TRUE)
})

test_that("canonical scalar and block-eigen sources remain protected", {
  protected <- c(
    "src/st_cpg_omp_csr_scheduled.cpp"="abeabf03db69e3358fb4850c0a432db2",
    "src/st_cpg_omp_csr_prior.cpp"="cce51072da6ddc3c18d58ab3b1f3c6df",
    "src/st_cpg_omp_csr_group.cpp"="87e923f7f8ee6420e39d9f041263d11b",
    "src/st_cpg_omp_csr_annot.cpp"="59bd49f048d116d0fe61d73d79bd4693",
    "NAMESPACE"="7519d0b7f23694a1ac78c1110bbf6e0b")
  actual <- unname(tools::md5sum(vapply(names(protected), blr_repo_path, character(1))))
  expect_identical(actual, unname(protected))
})
