phase9a_path<-function(...){parts<-c(...);if(identical(parts[1:3],c("tests","testthat","fixtures")))do.call(blr_fixture_path,as.list(parts[-(1:3)])) else do.call(blr_repo_path,as.list(parts))}
source(phase9a_path("tests","testthat","fixtures","blr-phase9a-annotation-reference.R"))
phase9a_common_spec<-function(){list(data=list(marker_count=4L,trait_count=1L,sample_size=50L,shared_read_only=TRUE,per_chain_payload=FALSE,storage_outlives_execution=TRUE),controls=list(iterations=6L,burnin=2L,thinning=1L,chains=2L,cores=2L,seed=9L,chain_seeds=c(11L,12L),keep_chains=TRUE,update_marker_variance=FALSE,update_residual_variance=FALSE,update_global_probability=FALSE,update_ld_swap=FALSE,ld_swap_probability=.05,ld_swap_r2=.8,ld_swap_max_friends=50L,ld_swap_moves=1L),output=list(keep_chains=TRUE,diagnostics=TRUE))}

test_that("Phase 9A fixed-prior contract validates and round-trips",{
 s<-phase9a_common_spec();s$policy<-list(marker_probability=c(.2,.3,.4,.5),marker_multiplier=c(.7,1,1.4,2),use_marker_probability=TRUE,use_marker_multiplier=TRUE,shared_read_only=TRUE,per_chain_payload=FALSE,storage_outlives_execution=TRUE)
 z<-sblr:::blr_phase9a_validate_fixed_prior_bayesc_cpp(s);expect_identical(z$spec,s);expect_true(z$validated);expect_false(z$invokes_sampler);expect_identical(z$policy,"fixed_marker_prior")
 bad<-s;bad$policy$marker_probability<-c(.2,.3);expect_error(sblr:::blr_phase9a_validate_fixed_prior_bayesc_cpp(bad),"length")
 bad<-s;bad$policy$marker_multiplier[1]<-0;expect_error(sblr:::blr_phase9a_validate_fixed_prior_bayesc_cpp(bad),"positive")
})
test_that("Phase 9A group contract validates mapping, order, and normalization",{
 s<-phase9a_common_spec();s$policy<-list(marker_group=c(0L,1L,0L,1L),group_count=2L,group_order=c("coding","background"),initial_probability=c(.35,.2),initial_multiplier=c(1.4,.7),prior_a=c(1,1),prior_b=c(1,1),update_probability=TRUE,update_multiplier=TRUE,normalize_multiplier=TRUE,zero_based_index=TRUE,shared_read_only=TRUE,per_chain_payload=FALSE,storage_outlives_execution=TRUE)
 z<-sblr:::blr_phase9a_validate_group_bayesc_cpp(s);expect_identical(z$spec,s);expect_false(z$invokes_sampler)
 bad<-s;bad$policy$marker_group[1]<-2L;expect_error(sblr:::blr_phase9a_validate_group_bayesc_cpp(bad),"index")
 bad<-s;bad$policy$marker_group<-rep(0L,4);expect_error(sblr:::blr_phase9a_validate_group_bayesc_cpp(bad),"empty")
})
test_that("Phase 9A learned-annotation contract validates exact link policy",{
 s<-phase9a_common_spec();A<-cbind(intercept=1,coding=c(1,0,1,0),qtl=c(0,1,1,0));s$policy<-list(marker_count=4L,annotation_count=3L,trait_count=1L,annotation=A,annotation_order=colnames(A),layout="column_major_marker_by_annotation",includes_intercept=TRUE,eta_probability_init=rep(0,3),eta_multiplier_init=rep(0,3),learn_probability=TRUE,learn_multiplier=TRUE,probability_prior_sd=1,multiplier_prior_sd=1,probability_proposal_sd=.02,multiplier_proposal_sd=.02,update_every=2L,probability_min=.01,probability_max=.9,multiplier_min=.1,multiplier_max=10,probability_link="centered_logit_offset",multiplier_link="centered_exponential",shared_read_only=TRUE,per_chain_payload=FALSE,storage_outlives_execution=TRUE)
 z<-sblr:::blr_phase9a_validate_learned_annotation_bayesc_cpp(s);expect_identical(z$spec,s);expect_false(z$invokes_sampler)
 bad<-s;bad$policy$update_every<-0L;expect_error(sblr:::blr_phase9a_validate_learned_annotation_bayesc_cpp(bad),"controls")
 bad<-s;bad$policy$probability_link<-"softmax";expect_error(sblr:::blr_phase9a_validate_learned_annotation_bayesc_cpp(bad),"link")
 bad<-s;bad$policy$shared_read_only<-FALSE;expect_error(sblr:::blr_phase9a_validate_learned_annotation_bayesc_cpp(bad),"borrowed")
})
test_that("Phase 9A common validation rejects invalid execution ownership and seeds",{
 s<-phase9a_common_spec();s$policy<-list(marker_probability=c(.2,.3,.4,.5),marker_multiplier=rep(1,4),use_marker_probability=TRUE,use_marker_multiplier=TRUE,shared_read_only=TRUE,per_chain_payload=FALSE,storage_outlives_execution=TRUE)
 bad<-s;bad$data$marker_count<-0L;expect_error(sblr:::blr_phase9a_validate_fixed_prior_bayesc_cpp(bad),"marker_count")
 bad<-s;bad$controls$chain_seeds<-1L;expect_error(sblr:::blr_phase9a_validate_fixed_prior_bayesc_cpp(bad),"chain_seeds")
 bad<-s;bad$data$per_chain_payload<-TRUE;expect_error(sblr:::blr_phase9a_validate_fixed_prior_bayesc_cpp(bad),"borrowed")
})
test_that("Phase 9A binding-neutral contracts contain no binding APIs",{
 x<-paste(readLines(phase9a_path("src","blr_csr_annotation_bayesc_types.h"),warn=FALSE),collapse="\n")
 for(tok in c("Rcpp","RcppArmadillo","SEXP","RObject","NumericVector","NumericMatrix","Nullable","Rcpp::stop","Rcpp::Rcout","R::rnorm","R::rchisq","arma::randn","arma::randu","pybind11","Python.h"))expect_false(grepl(tok,x,fixed=TRUE),info=tok)
})
test_that("Phase 9A all production raw and formatted references are exact",{
 for(backend in names(phase9a_configs))for(nm in names(phase9a_configs[[backend]])){ref<-readRDS(phase9a_path("tests","testthat","fixtures",paste0("blr_phase9a_",backend),paste0(nm,".rds")));cfg<-phase9a_configs[[backend]][[nm]];tol<-if(backend=="annotation"&&nm=="annot_learned")1e-8 else 1e-12;expect_equal(phase9a_normalize(phase9a_run(backend,cfg,TRUE)),ref$raw,tolerance=tol,info=paste(backend,nm,"raw"));expect_equal(phase9a_normalize(phase9a_run(backend,cfg,FALSE)),ref$fit,tolerance=tol,info=paste(backend,nm,"fit"))}
})
test_that("Phase 9A reproducibility is exact for each supported policy",{
 for(backend in names(phase9a_configs)){cfg<-phase9a_configs[[backend]][[2]];a<-phase9a_normalize(phase9a_run(backend,cfg,FALSE));b<-phase9a_normalize(phase9a_run(backend,cfg,FALSE));expect_identical(a,b,info=backend)}
})
test_that("Phase 9A core ordering and intervening annotation fits are exact",{
 comparable<-function(x){x<-phase9a_normalize(x);x$input$ncores<-0L;x}
 for(backend in names(phase9a_configs)){
  cfg<-phase9a_configs[[backend]][[2]];cfg$ncores<-1L;one<-comparable(phase9a_run(backend,cfg,FALSE))
  cfg$ncores<-2L;two<-comparable(phase9a_run(backend,cfg,FALSE));expect_identical(two,one,info=paste(backend,"1/2 cores"));expect_identical(comparable(phase9a_run(backend,cfg,FALSE)),two)
  cfg$ncores<-1L;expect_identical(comparable(phase9a_run(backend,cfg,FALSE)),one)
  other<-setdiff(names(phase9a_configs),backend)[1];invisible(phase9a_run(other,phase9a_configs[[other]][[1]],FALSE));expect_identical(comparable(phase9a_run(backend,cfg,FALSE)),one,info=paste(backend,"intervening"))
 }
})
test_that("Phase 9A learned-annotation production route remains structurally protected",{
 source_lines <- readLines(phase9a_path("src/st_cpg_omp_csr_annot.cpp"), warn = FALSE)
 core_lines <- readLines(phase9a_path("src/blr_csr_learned_annotation_bayesc_core_impl.h"), warn = FALSE)
 expect_equal(sum(grepl("for (int isort = 0; isort < m; ++isort)", core_lines, fixed = TRUE)), 1L)
 expect_true(any(grepl('#include "blr_csr_learned_annotation_bayesc_core_impl.h"', source_lines, fixed = TRUE)))
 expect_true(any(grepl("make_pi_from_annotation(A, eta_pi_t", core_lines, fixed = TRUE)))
 expect_true(any(grepl("make_vb_multiplier_from_annotation(", core_lines, fixed = TRUE)))
})
