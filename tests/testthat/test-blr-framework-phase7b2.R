phase7b2_path <- function(...) {
 p <- file.path(...)
 if (file.exists(p)) p else file.path("..", "..", ...)
}
source(phase7b2_path("tests", "testthat", "fixtures", "blr-phase7a-sbayesrc-reference.R"))

test_that("Phase 7B2 exposes one explicit operator-templated callable core", {
 src <- paste(readLines(phase7b2_path("src", "st_sbayesrc_omp_csr.cpp"), warn = FALSE), collapse = "\n")
 core <- paste(readLines(phase7b2_path("src", "blr_csr_sbayesrc_core_impl.h"), warn = FALSE), collapse = "\n")
 expect_match(core, "template <class Operator>\nCsrSBayesRCExecutionResult run_csr_sbayesrc(", fixed = TRUE)
 expect_match(core, "CsrSBayesRCExecutionContext<Operator>& context", fixed = TRUE)
 expect_match(src, "CsrSBayesRCExecutionContext<Operator> execution_context", fixed = TRUE)
 expect_match(src, "run_csr_sbayesrc(execution_context)", fixed = TRUE)
 expect_match(src, "SBayesRCOperatorContext<CsrOperator>", fixed = TRUE)
 expect_match(src, "SBayesRCOperatorContext<BlockEigenOperator>", fixed = TRUE)
 expect_equal(source_match_count("for (int it = 0; it < nit + nburn; ++it)", core, fixed = TRUE), 1L)
 expect_false(grepl("#include \"blr_csr_sbayesrc_core_impl.h\"\n\n const bool return_chain_summaries", src, fixed = TRUE))
 expect_false(any(vapply(c("getenv", "new_path", "old_path", "legacy_route"),
                         function(token) grepl(token, src, fixed = TRUE), logical(1))))
})

test_that("Phase 7B2 core is binding neutral and ownership is explicit", {
 core <- paste(readLines(phase7b2_path("src", "blr_csr_sbayesrc_core_impl.h"), warn = FALSE), collapse = "\n")
 for (token in c("Rcpp", "RcppArmadillo", "SEXP", "RObject", "NumericVector",
                  "NumericMatrix", "Rcpp::stop", "Rcpp::Rcout", "pybind11", "Python.h")) {
  expect_false(grepl(token, core, fixed = TRUE), info = token)
 }
 for (token in c("const Operator& op", "const arma::mat& A", "const arma::vec& gamma",
                  "const SBayesRCLDLDFriends& ld_swap_friends",
                  "storage_outlives_execution", "std::mt19937 gen_t(task_seed)")) {
  expect_match(core, token, fixed = TRUE)
 }
 expect_equal(source_match_count("#include \"blr_csr_sbayesrc_core_impl.h\"",
                              paste(readLines(phase7b2_path("src", "st_sbayesrc_omp_csr.cpp"), warn = FALSE), collapse = "\n"),
                              fixed = TRUE), 1L)
})

test_that("Phase 7B2 frozen raw and formatted references remain exact", {
 for (nm in names(phase7a_sbayesrc_configs)) {
  ref <- readRDS(phase7b2_path("tests", "testthat", "fixtures", "blr_phase7a_sbayesrc", paste0(nm, ".rds")))
  cfg <- phase7a_sbayesrc_configs[[nm]]
  expect_identical(phase7a_sbayesrc_normalize(phase7a_sbayesrc_run(cfg, TRUE)), ref$raw,
                   info = paste(nm, "raw"))
  expect_identical(phase7a_sbayesrc_normalize(phase7a_sbayesrc_run(cfg, FALSE)), ref$fit,
                   info = paste(nm, "fit"))
 }
})

test_that("Phase 7B2 reproducibility and alpha modes remain exact", {
 comparable <- function(x) {
  x <- phase7a_sbayesrc_normalize(x)
  x$input$ncores <- 0L
  x
 }
 cfg <- phase7a_sbayesrc_configs$learned_two_cores
 cfg$ncores <- 1L
 one <- comparable(phase7a_sbayesrc_run(cfg, FALSE))
 expect_identical(comparable(phase7a_sbayesrc_run(cfg, FALSE)), one)
 cfg$ncores <- 2L
 two <- comparable(phase7a_sbayesrc_run(cfg, FALSE))
 expect_identical(two, one)
 expect_identical(comparable(phase7a_sbayesrc_run(cfg, FALSE)), two)
 cfg$ncores <- 1L
 expect_identical(comparable(phase7a_sbayesrc_run(cfg, FALSE)), one)

 fixed <- phase7a_sbayesrc_configs$fixed_one_chain
 expect_identical(phase7a_sbayesrc_normalize(phase7a_sbayesrc_run(fixed, FALSE)),
                  readRDS(phase7b2_path("tests", "testthat", "fixtures", "blr_phase7a_sbayesrc", "fixed_one_chain.rds"))$fit)
 kept <- phase7a_sbayesrc_configs$learned_explicit_keep
 expect_identical(phase7a_sbayesrc_normalize(phase7a_sbayesrc_run(kept, FALSE)),
                  readRDS(phase7b2_path("tests", "testthat", "fixtures", "blr_phase7a_sbayesrc", "learned_explicit_keep.rds"))$fit)
})

test_that("Phase 7B2 keeps one named converter and protected source files", {
 src <- paste(readLines(phase7b2_path("src", "st_sbayesrc_omp_csr.cpp"), warn = FALSE), collapse = "\n")
 expect_equal(source_match_count("static Rcpp::List stblr_csr_sbayesrc_result_to_raw(", src, fixed = TRUE), 2L)
 expect_equal(source_match_count("return stblr_csr_sbayesrc_result_to_raw(execution_result, binding_metadata);", src, fixed = TRUE), 1L)
 expect_match(src, "Rcpp::List comp_prob_out", fixed = TRUE)
 expect_match(src, "Rcpp::List marker = Rcpp::List::create", fixed = TRUE)
 paths <- c("src/st_cpg_omp_csr.cpp", "src/blr_csr_bayesc_types.h", "src/blr_csr_bayesc_core_impl.h",
            "src/st_cpg_omp_csr_bayesr.cpp", "src/blr_csr_bayesr_types.h", "src/blr_csr_bayesr_core_impl.h",
            "src/st_block_eigen.cpp", "src/st_block_eigen.h",
            "src/stblr_cpg_omp_bed_marker_scheduled_chains_bayesrc.cpp",
             "NAMESPACE")
 expected <- c("92dafc0266d5a0e72aea000224154cef", "e5975c311c69fe536db57dd21f01334f", "f7c617cbfc172639c1f8aea1bd8b1876",
               "0a005f9d5a19037285fd4869fdc4dcf0", "bf1d4b73065207ca361c7abdab3cb253", "4dac6bef2df917613df8e1a827640303",
               "49f0a62c9fe235967a264b0f8de144a7", "bec3bc1e41841ab77747e34dc9818574", "21ee1d04ae816644c4d918202b29d515",
               "f5b6ee37a3972aa436357bdc8f602f4e")
 expect_identical(unname(tools::md5sum(vapply(paths, phase7b2_path, character(1)))), expected)
})
