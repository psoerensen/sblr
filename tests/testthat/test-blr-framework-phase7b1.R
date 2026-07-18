phase7b1_path<-function(...) { p<-file.path(...); if(file.exists(p)) p else file.path("..","..",...) }
source(phase7b1_path("tests","testthat","fixtures","blr-phase7a-sbayesrc-reference.R"))

test_that("Phase 7B1 has one active shared SBayesRC MCMC block", {
 src.lines<-readLines(phase7b1_path("src","st_sbayesrc_omp_csr.cpp"),warn=FALSE)
 core.lines<-readLines(phase7b1_path("src","blr_csr_sbayesrc_core_impl.h"),warn=FALSE)
 active<-function(x,pattern) sum(grepl(pattern,x,fixed=TRUE)&!grepl("^\\s*//",x))
 expect_equal(active(src.lines,"for (int it = 0; it < nit + nburn; ++it)"),0L)
 expect_equal(active(core.lines,"for (int it = 0; it < nit + nburn; ++it)"),1L)
 expect_equal(active(src.lines,"#include \"blr_csr_sbayesrc_core_impl.h\""),1L)
 expect_equal(active(core.lines,"CsrSBayesRCExecutionResult run_csr_sbayesrc("),1L)
 expect_match(paste(core.lines,collapse="\n"),"SBLR_BLR_CSR_SBAYESRC_CORE_IMPL_H",fixed=TRUE)
 expect_false(any(grepl("getenv",c(src.lines,core.lines),fixed=TRUE)))
})

test_that("Phase 7B1 keeps operator sharing and the permanent converter", {
 src<-paste(readLines(phase7b1_path("src","st_sbayesrc_omp_csr.cpp"),warn=FALSE),collapse="\n")
 expect_match(src,"SBayesRCOperatorContext<CsrOperator>",fixed=TRUE)
 expect_match(src,"SBayesRCOperatorContext<BlockEigenOperator>",fixed=TRUE)
 expect_equal(source_match_count("#include \"blr_csr_sbayesrc_core_impl.h\"",src,fixed=TRUE),1L)
 expect_match(src,"run_csr_sbayesrc(execution_context)",fixed=TRUE)
 expect_match(src,"stblr_csr_sbayesrc_result_to_raw",fixed=TRUE)
 expect_match(src,"Rcpp::List comp_prob_out",fixed=TRUE)
 expect_match(src,"Rcpp::List marker = Rcpp::List::create",fixed=TRUE)
})

test_that("Phase 7B1 frozen references remain exact", {
 for(nm in names(phase7a_sbayesrc_configs)) {
  ref<-readRDS(phase7b1_path("tests","testthat","fixtures","blr_phase7a_sbayesrc",paste0(nm,".rds"))); cfg<-phase7a_sbayesrc_configs[[nm]]
  expect_identical(phase7a_sbayesrc_normalize(phase7a_sbayesrc_run(cfg,TRUE)),ref$raw,info=paste(nm,"raw"))
  expect_identical(phase7a_sbayesrc_normalize(phase7a_sbayesrc_run(cfg,FALSE)),ref$fit,info=paste(nm,"fit"))
 }
})
