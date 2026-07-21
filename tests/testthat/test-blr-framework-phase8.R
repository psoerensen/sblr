phase8_path <- function(...) {
 parts <- c(...)
 if (identical(parts[1:3], c("tests", "testthat", "fixtures")))
   do.call(blr_fixture_path, as.list(parts[-(1:3)])) else
   do.call(blr_repo_path, as.list(parts))
}
source(phase8_path("tests", "testthat", "fixtures", "blr-phase7a-sbayesrc-reference.R"))

test_that("Phase 8 canonical SBayesRC architecture has one core and converter", {
 src <- paste(readLines(phase8_path("src", "st_sbayesrc_omp_csr.cpp"), warn = FALSE), collapse = "\n")
 core <- paste(readLines(phase8_path("src", "blr_csr_sbayesrc_core_impl.h"), warn = FALSE), collapse = "\n")
 all_cpp <- paste(vapply(list.files(phase8_path("src"), pattern = "\\.(cpp|h)$", full.names = TRUE),
                         function(p) paste(readLines(p, warn = FALSE), collapse = "\n"), character(1)), collapse = "\n")
 expect_equal(source_match_count("CsrSBayesRCExecutionResult run_csr_sbayesrc(", core, fixed = TRUE), 1L)
 expect_equal(source_match_count("for (int it = 0; it < nit + nburn; ++it)", core, fixed = TRUE), 1L)
 expect_equal(source_match_count("static Rcpp::List stblr_csr_sbayesrc_result_to_raw(", src, fixed = TRUE), 2L)
 expect_equal(source_match_count("return stblr_csr_sbayesrc_result_to_raw(execution_result, binding_metadata);", src, fixed = TRUE), 1L)
 expect_equal(source_match_count("#include \"blr_csr_sbayesrc_core_impl.h\"", all_cpp, fixed = TRUE), 1L)
 expect_match(src, "CsrSBayesRCExecutionContext<Operator> execution_context", fixed = TRUE)
 expect_match(src, "SBayesRCOperatorContext<CsrOperator>", fixed = TRUE)
 expect_match(src, "SBayesRCOperatorContext<BlockEigenOperator>", fixed = TRUE)
 expect_false(any(vapply(c("getenv", "old_path", "new_path", "legacy_route", "route_selector", "fallback_route"),
                         function(x) grepl(x, src, fixed = TRUE), logical(1))))
})

test_that("Phase 8 canonical core is binding neutral with explicit ownership", {
 core <- paste(readLines(phase8_path("src", "blr_csr_sbayesrc_core_impl.h"), warn = FALSE), collapse = "\n")
 for (token in c("Rcpp", "RcppArmadillo", "SEXP", "RObject", "NumericVector", "NumericMatrix",
                  "Rcpp::stop", "Rcpp::Rcout", "pybind11", "Python.h"))
  expect_false(grepl(token, core, fixed = TRUE), info = token)
 for (token in c("const Operator& op", "const arma::mat& A", "const arma::vec& gamma",
                  "const SBayesRCLDLDFriends& ld_swap_friends", "storage_outlives_execution",
                  "std::mt19937 gen_t(task_seed)", "arma::mat b_task(ntasks, m"))
  expect_match(core, token, fixed = TRUE)
 expect_false(grepl("std::vector<CsrSBayesRCDataView>", core, fixed = TRUE))
 expect_false(grepl("std::vector<SBayesRCAnnotationDesignView>", core, fixed = TRUE))
})

test_that("Phase 8 permanent raw and formatted references remain exact", {
 for (nm in names(phase7a_sbayesrc_configs)) {
  ref <- readRDS(phase8_path("tests", "testthat", "fixtures", "blr_phase7a_sbayesrc", paste0(nm, ".rds")))
  cfg <- phase7a_sbayesrc_configs[[nm]]
  expect_equal(phase7a_sbayesrc_normalize(phase7a_sbayesrc_run(cfg, TRUE)), ref$raw, tolerance=1e-12, info = paste(nm, "raw"))
  expect_equal(phase7a_sbayesrc_normalize(phase7a_sbayesrc_run(cfg, FALSE)), ref$fit, tolerance=1e-12, info = paste(nm, "fit"))
 }
})

test_that("Phase 8 reproducibility is exact across scheduling and intervening fits", {
 comparable <- function(x) { x <- phase7a_sbayesrc_normalize(x); x$input$ncores <- 0L; x }
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
 source(phase8_path("tests", "testthat", "fixtures", "blr-phase5a-bayesr-reference.R"), local = TRUE)
 invisible(phase5a_bayesr_run(phase5a_bayesr_configs$one_chain, FALSE))
 expect_identical(comparable(phase7a_sbayesrc_run(cfg, FALSE)), one)
})

test_that("Phase 8 annotation, alpha, probability, and chain contracts remain exact", {
 for (nm in c("fixed_one_chain", "learned_explicit_keep", "multiple_traits", "explicit_scales")) {
  ref <- readRDS(phase8_path("tests", "testthat", "fixtures", "blr_phase7a_sbayesrc", paste0(nm, ".rds")))
  fit <- phase7a_sbayesrc_normalize(phase7a_sbayesrc_run(phase7a_sbayesrc_configs[[nm]], FALSE))
  expect_equal(fit, ref$fit, tolerance=1e-12, info = nm)
  expect_equal(fit$comp_prob, ref$fit$comp_prob, tolerance=1e-12, info = paste(nm, "component probabilities"))
  expect_identical(fit$input$annotation, ref$fit$input$annotation, info = paste(nm, "annotation"))
 }
})

test_that("Phase 8 protected implementations and public namespace are unchanged", {
 paths <- c("src/blr_csr_bayesc_types.h", "src/blr_csr_bayesc_core_impl.h",
            "src/blr_csr_bayesr_types.h", "src/blr_csr_bayesr_core_impl.h",
            "src/stblr_cpg_omp_bed_marker_scheduled_chains_bayesrc.cpp",
             "NAMESPACE")
 expected <- c("4d0eb5380007195a8d34e7b2e081dec4", "c548157cc9e5804272e714983bdcb798",
               "bf1d4b73065207ca361c7abdab3cb253", "4dac6bef2df917613df8e1a827640303",
               "72d4a9fa0a7cd51071328c2d62d0192b",
               "a1f389e8ea9ab5abef440767a11b8378")
 expect_identical(unname(tools::md5sum(vapply(paths, phase8_path, character(1)))), expected)
})
