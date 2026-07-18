phase14a_root <- normalizePath(file.path(testthat::test_path(),"..",".."),winslash="/",mustWork=TRUE)
source(file.path(testthat::test_path(),"fixtures","blr-phase14a-bed-bayesrc-reference.R"))
phase14a_text <- function(path) paste(readLines(file.path(phase14a_root,path),warn=FALSE),collapse="\n")

test_that("Phase 14A route and full-sweep production seam are discoverable",{
 r <- phase14a_text("R/sparse_ld_bed_helper.R"); n <- phase14a_text("R/stblr-bed-bayesrc-internal.R")
 cpp <- paste(phase14a_text("src/stblr_cpg_omp_bed_marker_scheduled_chains_bayesrc.cpp"),
              phase14a_text("src/blr_bed_bayesrc_core_impl.h"),sep="\n")
 expect_match(r,'method = c("bayesc", "bayesr", "bayesrc")',fixed=TRUE)
 expect_match(r,".stblr_bed_bayesrc_native(",fixed=TRUE)
 expect_match(n,".stblr_align_bed_bayesrc_annotations",fixed=TRUE)
 expect_source_count("run_bed_bayesrc_chain(",cpp,2L)
 expect_source_count("for (int it = 0; it < total_it; ++it)",cpp,1L)
 expect_match(cpp,"Exact full sweep: every marker is visited once",fixed=TRUE)
 expect_false(grepl("candidate_list|scheduled_at|due_buckets|null_skip",cpp))
})

test_that("Phase 14A binding-neutral audit contracts encode the model",{
 h <- phase14a_text("src/blr_bed_bayesrc_audit_types.h")
 for(x in c("BedBayesRCComponentSpec","BedBayesRCAnnotationSpec","BedBayesRCCoefficientPriorSpec",
             "BedBayesRCExecutionAuditSpec","BedBayesRCOwnershipAuditSpec",
             "BedBayesRCChainResultVocabulary","BedBayesRCExecutionResultVocabulary"))
  expect_source_count(paste0("struct ",x),h,1L)
 for(x in c("null_component must be zero","annotation coefficient dimensions mismatch",
             "packed-BED BayesRC is full-sweep only","invalid BayesRC ownership contract"))
  expect_match(h,x,fixed=TRUE)
 expect_false(grepl("Rcpp|SEXP|Python.h|pybind11|FILE|fopen",h))
})

test_that("Phase 14A annotation alignment preserves production policy",{
 ids <- c("rs1","rs2")
 x <- data.frame(marker_id=c("extra","rs2","rs1"),binary=c(0,1,0),continuous=c(3,2,1))
 z <- sblr:::.stblr_align_bed_bayesrc_annotations(x,ids,add_intercept=TRUE,
  standardize_annotations=FALSE)
 expect_identical(rownames(z$A),ids); expect_identical(colnames(z$A),c("Intercept","binary","continuous"))
 expect_equal(unname(z$A[,"binary"]),c(0,1)); expect_identical(z$alignment$unused_annotation_rows,1L)
 expect_true(z$preprocessing$intercept_added)
 expect_error(sblr:::.stblr_align_bed_bayesrc_annotations(transform(x,marker_id=c("rs1","rs1","x")),ids),"unique")
 expect_error(sblr:::.stblr_align_bed_bayesrc_annotations(x[x$marker_id!="rs2",],ids),"missing")
 expect_error(sblr:::.stblr_align_bed_bayesrc_annotations(cbind(intercept=rep(1,2),constant=rep(2,2)),ids),"positive variance")
 expect_error(sblr:::.stblr_align_bed_bayesrc_annotations(matrix(c(0,NA),2),ids),"non-finite")
})

test_that("Phase 14A probit-stick and annotation identities hold",{
 target <- c(.95,.03,.015,.005)
 alpha <- sblr:::.bayesr_pi_to_probit_stick_intercepts(target)
 p <- pnorm(alpha); remaining <- 1; got <- numeric(4)
 for(k in 1:3){got[k] <- remaining*(1-p[k]); remaining <- remaining*p[k]}; got[4] <- remaining
 expect_equal(got,target,tolerance=1e-12); expect_equal(sum(got),1,tolerance=1e-12)
 z <- phase14a_capture(1,1,141,TRUE,TRUE)$raw; cp <- z$component$prob[[1L]]
 expect_true(all(is.finite(cp)&cp>=0)); expect_equal(rowSums(cp),rep(1,nrow(cp)),tolerance=1e-12)
 expect_equal(z$marker$dm[,1],1-cp[,1],tolerance=1e-12)
 expect_true(all(z$marker$state[,1]>=0&z$marker$state[,1]<ncol(cp)))
 expect_equal(dim(z$annotation$alpha_mean[[1L]]),c(3L,3L))
 expect_equal(z$pi$final[1,],colMeans(z$annotation$marker_prior_final[[1L]]),tolerance=1e-12)
})

test_that("Phase 14A frozen raw and formatted references are exact",{
 cfg <- list(one_chain_one_core=c(1,1,141,1,0),two_chains_one_core=c(1,2,143,1,1),two_chains_two_cores=c(2,2,143,1,1))
 for(nm in names(cfg)){z<-cfg[[nm]]; ref<-readRDS(file.path(testthat::test_path(),"fixtures","blr_phase14a_bed_bayesrc",paste0(nm,".rds")))
  observed<-phase14a_normalize(phase14a_capture(z[1],z[2],z[3],as.logical(z[4]),as.logical(z[5])))
  expect_identical(observed$raw,ref$raw); expect_identical(observed$fit,ref$fit)
 }
})

test_that("Phase 14A same-process, core-order and intervening fits are exact",{
 a <- phase14a_normalize(phase14a_capture(1,2,143,TRUE,TRUE))
 expect_identical(phase14a_normalize(phase14a_capture(1,2,143,TRUE,TRUE)),a)
 phase14a_capture(1,1,149,FALSE,FALSE)
 expect_identical(phase14a_normalize(phase14a_capture(1,2,143,TRUE,TRUE)),a)
 expect_identical(phase14a_normalize(phase14a_capture(2,2,143,TRUE,TRUE)),a)
 phase11a_capture("bayesr",1,1,71); phase11a_capture("bayesc",1,1,71)
 expect_identical(phase14a_normalize(phase14a_capture(1,2,143,TRUE,TRUE)),a)
 one <- phase14a_normalize(phase14a_capture(1,1,143,TRUE,TRUE))
 phase14a_capture(1,2,151,TRUE,TRUE)
 expect_identical(phase14a_normalize(phase14a_capture(1,1,143,TRUE,TRUE)),one)
})

test_that("Phase 14A fixed intercept reduction remains BayesR-exact",{
 x<-phase11a_fixture(); target<-c(.95,.03,.015,.005)
 fixed<-phase14a_capture(1,1,141,FALSE,FALSE)$raw
 expect_equal(fixed$pi$final[1,],target,tolerance=1e-12)
 expect_equal(unname(fixed$annotation$alpha_final[[1]]),unname(sblr:::.bayesr_pi_to_probit_stick_intercepts(target)),tolerance=1e-12)
 expect_match(phase14a_text("tests/testthat/test-individual-bayesrc.R"),"reduces to fixed-pi BayesR",fixed=TRUE)
 expect_error(sblr::stblr_bed(y=x$y,Glist=x$Glist,method="bayesrc",annotation=NULL),"required")
})

test_that("Phase 14A production and protected sources remain frozen",{
 protected<-c("src/stblr_cpg_omp_bed_marker_scheduled_chains_bayesrc.cpp"="9812450ae103f5026d1632b2bcc31e95",
 "src/blr_bed_bayesrc_core_impl.h"="82365cf3f1f5306c57b980f59b4d83d3",
 "src/st_bayesrc_annotation_prior.h"="1e7072512f4246fc2a36e79de655d8c5","src/blr_bed_bayesr_core_impl.h"="afe77e26d2cf2b8e3d64088221b33e14",
 "src/blr_bed_scheduled_bayesc_core_impl.h"="723cee003504c1fdcd075b965cb63d83","src/blr_csr_sbayesrc_core_impl.h"="d06ec2a530e8c914201ee22b6be65739",
 "R/RcppExports.R"="9d13ea00b326c7e0cd606194d13a8bca","src/RcppExports.cpp"="b4859db0f6308fa7e38051ddcf32d245","NAMESPACE"="f5b6ee37a3972aa436357bdc8f602f4e")
 expect_identical(unname(tools::md5sum(file.path(phase14a_root,names(protected)))),unname(protected))
})

test_that("Phase 14A fresh-process references can be checked explicitly",{
 skip_if_not(identical(Sys.getenv("SBLR_RUN_PHASE14A_FRESH"),"true")); skip_if_not_installed("callr")
 observed<-callr::r(function(root){setwd(root);pkgload::load_all(".",compile=FALSE,quiet=TRUE);library(testthat);source("tests/testthat/fixtures/blr-phase14a-bed-bayesrc-reference.R");phase14a_normalize(phase14a_capture(2,2,143,TRUE,TRUE))},list(root=phase14a_root))
 ref<-readRDS(file.path(testthat::test_path(),"fixtures","blr_phase14a_bed_bayesrc","two_chains_two_cores.rds"))
 expect_identical(observed$raw,ref$raw);expect_identical(observed$fit,ref$fit)
})
