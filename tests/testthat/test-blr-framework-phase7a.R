phase7a_path<-function(...) { p<-file.path(...); if(file.exists(p)) p else file.path("..","..",...) }
source(phase7a_path("tests","testthat","fixtures","blr-phase7a-sbayesrc-reference.R"))

phase7a_contract<-function() list(
 data=list(marker_count=4,trait_count=2,sample_size=c(80L,80L),marker_order=paste0("m",1:4),trait_order=c("T1","T2"),shared_read_only=TRUE,per_chain_payload=FALSE,storage_outlives_execution=TRUE),
 annotation=list(marker_count=4,annotation_count=3,values=matrix(c(1,1,1,1,0,1,0,1,1,0,0,1),4,3),annotation_order=c("intercept","coding","qtl"),layout="column_major",includes_intercept=TRUE,standardized=FALSE,centered_binary=FALSE,shared_read_only=TRUE,per_chain_payload=FALSE,storage_outlives_execution=TRUE),
 component=list(scales=c(0,.1,1),null_component=0,scale_interpretation="variance_multiplier"),
 alpha=list(annotation_count=3,step_count=2,initial_values=rep(0,6),initial_variance=c(1,1),intercept_flat=TRUE,update=TRUE,variance_prior_a=2,variance_prior_b=2,update_every=2L),
 probability=list(transformation="probit_stick_breaking",reference_component=0,stick_order=c(1,2),probability_floor=1e-12,floor_then_normalize=TRUE),
 controls=list(iterations=10L,burnin=2L,thinning=1L,chains=2L,cores=2L,seed=73L,chain_seeds=c(701L,702L),keep_chains=TRUE,update_ld_swap=FALSE,ld_swap_probability=.05,ld_swap_r2=.8,ld_swap_max_friends=50L,ld_swap_moves=1L),output=list(keep_chains=TRUE,diagnostics=TRUE))

test_that("Phase 7A SBayesRC contract round-trips exactly without execution", {
 x<-phase7a_contract(); y<-blr_phase7a_validate_sbayesrc_contract_cpp(x)
 expect_true(y$validated); expect_false(y$invokes_sampler)
 for(nm in names(x)) expect_identical(y[[nm]],x[[nm]],info=nm)
})

test_that("Phase 7A validation rejects ownership, annotation, component, alpha, and execution defects", {
 bad<-function(group,field,value) { x<-phase7a_contract(); x[[group]][[field]]<-value; x }
 cases<-list(
  list("data","marker_count",0,"marker_count"),list("data","trait_count",0,"trait_count"),
  list("data","shared_read_only",FALSE,"read-only"),list("data","per_chain_payload",TRUE,"per-chain"),list("data","storage_outlives_execution",FALSE,"outlives"),
  list("annotation","marker_count",3,"dimensions"),list("annotation","annotation_count",0,"dimensions"),list("annotation","annotation_order",c("a"),"annotation_order"),
  list("annotation","values",matrix(1,3,3),"dimensions"),list("annotation","values",matrix(c(Inf,rep(0,11)),4,3),"finite"),list("annotation","layout","row_major","column_major"),
  list("annotation","shared_read_only",FALSE,"read-only"),list("annotation","per_chain_payload",TRUE,"per-chain"),
  list("component","scales",c(.1,1),"null"),list("component","null_component",1,"null"),list("component","scales",c(0,-.1,1),"scales"),
  list("alpha","step_count",3,"dimensions"),list("alpha","initial_values",rep(0,5),"dimensions"),list("alpha","initial_values",c(Inf,rep(0,5)),"finite"),
  list("alpha","initial_variance",c(1,-1),"positive"),list("alpha","variance_prior_a",0,"priors"),list("alpha","update_every",0L,"update_every"),
  list("probability","transformation","softmax","probit"),list("probability","stick_order",c(2,1),"stick_order"),list("probability","probability_floor",1,"probability_floor"),
  list("controls","iterations",0L,"MCMC"),list("controls","chains",0L,"counts"),list("controls","cores",0L,"counts"),list("controls","chain_seeds",701L,"chain_seeds"),list("controls","ld_swap_probability",2,"LD-swap") )
 for(z in cases) expect_error(blr_phase7a_validate_sbayesrc_contract_cpp(bad(z[[1]],z[[2]],z[[3]])),z[[4]],info=paste(z[[1]],z[[2]]))
})

test_that("Phase 7A typed source is binding neutral and complete", {
 z<-paste(readLines(phase7a_path("src","blr_csr_sbayesrc_types.h"),warn=FALSE),collapse="\n")
 for(token in c("Rcpp","RcppArmadillo","SEXP","RObject","NumericVector","NumericMatrix","Nullable","Rcpp::stop","Rcpp::Rcout","R::rnorm","R::rchisq","R::pnorm","R::qnorm","arma::randn","arma::randu","pybind11","Python.h")) expect_false(grepl(token,z,fixed=TRUE),info=token)
 for(token in c("CsrSBayesRCDataView","SBayesRCAnnotationDesignView","SBayesRCAlphaSpec","SBayesRCProbabilityPolicy","CsrSBayesRCExecutionInput","CsrSBayesRCResult")) expect_match(z,token,fixed=TRUE)
})

test_that("Phase 7A frozen production SBayesRC raw and formatted references remain exact", {
 for(nm in names(phase7a_sbayesrc_configs)) {
  ref<-readRDS(phase7a_path("tests","testthat","fixtures","blr_phase7a_sbayesrc",paste0(nm,".rds"))); cfg<-phase7a_sbayesrc_configs[[nm]]
  expect_identical(phase7a_sbayesrc_normalize(phase7a_sbayesrc_run(cfg,TRUE)),ref$raw,info=paste(nm,"raw"))
  expect_identical(phase7a_sbayesrc_normalize(phase7a_sbayesrc_run(cfg,FALSE)),ref$fit,info=paste(nm,"fit"))
 }
})

test_that("Phase 7A SBayesRC reproducibility is exact across cores and call order", {
 comparable<-function(x) { x<-phase7a_sbayesrc_normalize(x); x$input$ncores<-0L; x }
 cfg<-phase7a_sbayesrc_configs$learned_two_cores; cfg$ncores<-1L
 one<-comparable(phase7a_sbayesrc_run(cfg,FALSE))
 repeat_one<-comparable(phase7a_sbayesrc_run(cfg,FALSE))
 cfg$ncores<-2L; two<-comparable(phase7a_sbayesrc_run(cfg,FALSE))
 expect_identical(repeat_one,one)
 expect_identical(two,one)
 expect_identical(comparable(phase7a_sbayesrc_run(cfg,FALSE)),two)
 cfg$ncores<-1L
 expect_identical(comparable(phase7a_sbayesrc_run(cfg,FALSE)),one)
 source(phase7a_path("tests","testthat","fixtures","blr-phase5a-bayesr-reference.R"),local=TRUE)
 invisible(phase5a_bayesr_run(phase5a_bayesr_configs$one_chain,FALSE))
 expect_identical(comparable(phase7a_sbayesrc_run(cfg,FALSE)),one)
})

test_that("Phase 7A protects extracted production structure and backend hashes", {
 src<-paste(readLines(phase7a_path("src","st_sbayesrc_omp_csr.cpp"),warn=FALSE),collapse="\n")
 core<-paste(readLines(phase7a_path("src","blr_csr_sbayesrc_core_impl.h"),warn=FALSE),collapse="\n")
 expect_match(src,"#include \"blr_csr_sbayesrc_core_impl.h\"",fixed=TRUE)
 expect_match(core,"for (int it = 0; it < nit + nburn; ++it)",fixed=TRUE)
 paths<-c("src/st_cpg_omp_csr.cpp","src/blr_csr_bayesc_types.h","src/blr_csr_bayesc_core_impl.h","src/st_cpg_omp_csr_bayesr.cpp","src/blr_csr_bayesr_types.h","src/blr_csr_bayesr_core_impl.h","src/st_block_eigen.cpp","src/st_block_eigen.h","src/stblr_cpg_omp_bed_marker_scheduled_chains_bayesrc.cpp","src/st_cpg_omp_csr_annot.cpp","NAMESPACE")
 expected<-c("92dafc0266d5a0e72aea000224154cef","e5975c311c69fe536db57dd21f01334f","f7c617cbfc172639c1f8aea1bd8b1876","0a005f9d5a19037285fd4869fdc4dcf0","bf1d4b73065207ca361c7abdab3cb253","4dac6bef2df917613df8e1a827640303","49f0a62c9fe235967a264b0f8de144a7","bec3bc1e41841ab77747e34dc9818574","5904c60b32165a7ae73bfc9d6c0f920c","baaf3a0919ba97c78401066f7ac7d6f3","f5b6ee37a3972aa436357bdc8f602f4e")
 expect_identical(unname(tools::md5sum(vapply(paths,phase7a_path,character(1)))),expected)
})
