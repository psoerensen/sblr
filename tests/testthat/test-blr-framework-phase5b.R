phase5b_path <- function(...) { p<-file.path(...); if(file.exists(p)) p else file.path("..","..",...) }
source(phase5b_path("tests","testthat","fixtures","blr-phase5a-bayesr-reference.R"))

test_that("Phase 5B preserves all frozen BayesR raw and formatted results", {
  for(nm in names(phase5a_bayesr_configs)) {
    cfg<-phase5a_bayesr_configs[[nm]]
    ref<-readRDS(phase5b_path("tests","testthat","fixtures","blr_phase5a_bayesr",paste0(nm,".rds")))
    expect_identical(phase5a_bayesr_normalize(phase5a_bayesr_run(cfg,TRUE)),ref$raw,info=paste(nm,"raw"))
    expect_identical(phase5a_bayesr_normalize(phase5a_bayesr_run(cfg,FALSE)),ref$fit,info=paste(nm,"fit"))
  }
})

test_that("Phase 5B has one core loop and one binding converter", {
  binding<-paste(readLines(phase5b_path("src","st_cpg_omp_csr_bayesr.cpp"),warn=FALSE),collapse="\n")
  core<-paste(readLines(phase5b_path("src","blr_csr_bayesr_core_impl.h"),warn=FALSE),collapse="\n")
  expect_equal(length(gregexpr("for (int it = 0; it < trace_len",paste(binding,core),fixed=TRUE)[[1L]]),1L)
  expect_equal(length(gregexpr("csr_bayesr_result_to_raw",binding,fixed=TRUE)[[1L]]),2L)
  expect_match(binding,"run_bayesr_execution(execution_context)",fixed=TRUE)
  expect_false(grepl("getenv",paste(binding,core),fixed=TRUE))
  for(token in c("Rcpp","SEXP","RObject","NumericVector","NumericMatrix","pybind11","Python.h"))
    expect_false(grepl(token,core,fixed=TRUE),info=token)
})
