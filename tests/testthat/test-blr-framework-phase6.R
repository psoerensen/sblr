phase6_path<-function(...){parts<-c(...);if(identical(parts[1:3],c("tests","testthat","fixtures")))do.call(blr_fixture_path,as.list(parts[-(1:3)])) else do.call(blr_repo_path,as.list(parts))}
source(phase6_path("tests","testthat","fixtures","blr-phase5a-bayesr-reference.R"))

test_that("canonical CSR BayesR retains permanent exact references",{
 for(nm in names(phase5a_bayesr_configs)){
  cfg<-phase5a_bayesr_configs[[nm]]
  ref<-readRDS(phase6_path("tests","testthat","fixtures","blr_phase5a_bayesr",paste0(nm,".rds")))
  expect_equal(phase5a_bayesr_raw_science(phase5a_bayesr_normalize(phase5a_bayesr_run(cfg,TRUE))),
    phase5a_bayesr_raw_science(ref$raw),tolerance=1e-12,info=paste(nm,"raw scientific fields"))
  expect_equal(phase5a_bayesr_fit_science(phase5a_bayesr_normalize(phase5a_bayesr_run(cfg,FALSE))),
    phase5a_bayesr_fit_science(ref$fit),tolerance=1e-12,info=paste(nm,"fit scientific fields"))
 }
})

test_that("canonical CSR BayesR has one guarded core and converter",{
 binding<-paste(readLines(phase6_path("src","st_cpg_omp_csr_bayesr.cpp"),warn=FALSE),collapse="\n")
 core<-paste(readLines(phase6_path("src","blr_csr_bayesr_core_impl.h"),warn=FALSE),collapse="\n")
 all<-paste(binding,core)
 expect_equal(source_match_count("for (int it = 0; it < trace_len",all,fixed=TRUE),1L)
 expect_equal(source_match_count("stblr_csr_bayesr_result_to_raw",binding,fixed=TRUE),2L)
 expect_match(binding,"run_csr_bayesr(execution_context)",fixed=TRUE)
 expect_match(core,"template <class Operator>",fixed=TRUE)
 expect_match(core,"SBLR_CSR_BAYESR_CORE_IMPL_TRANSLATION_UNIT",fixed=TRUE)
 expect_equal(length(grep("blr_csr_bayesr_core_impl.h",list.files(phase6_path("src"),full.names=TRUE),value=TRUE)),1L)
 expect_false(grepl("getenv",all,fixed=TRUE))
 for(token in c("Rcpp","SEXP","RObject","NumericVector","NumericMatrix","pybind11","Python.h"))expect_false(grepl(token,core,fixed=TRUE),info=token)
 expect_match(core,"std::mt19937 gen_t(task_seed)",fixed=TRUE)
 expect_match(core,"const arma::rowvec& ww_t = op.diag()",fixed=TRUE)
})
