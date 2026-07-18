phase14c_root<-normalizePath(file.path(testthat::test_path(),"..",".."),winslash="/",mustWork=TRUE)
source(file.path(testthat::test_path(),"fixtures","blr-phase14a-bed-bayesrc-reference.R"))
phase14c_text<-function(path) paste(readLines(file.path(phase14c_root,path),warn=FALSE),collapse="\n")

test_that("Phase 14C activates one typed chain boundary",{
 types<-phase14c_text("src/blr_bed_bayesrc_types.h")
 core<-phase14c_text("src/blr_bed_bayesrc_core_impl.h")
 adapter<-phase14c_text("src/stblr_cpg_omp_bed_marker_scheduled_chains_bayesrc.cpp")
 for(x in c("BedBayesRCComponentSpec","BedBayesRCCoefficientPriorSpec",
  "BedBayesRCChainExecutionResult")) expect_source_count(paste0("struct ",x),types,1L)
 for(x in c("BedBayesRCAnnotationSpec","BedBayesRCPackedGenotypeView",
  "BedBayesRCChainExecutionContext")) expect_source_count(paste0("struct ",x),types,1L)
 expect_source_count("run_bed_bayesrc_chain(",paste(core,adapter),2L)
 expect_source_count("for (int it = 0; it < total_it; ++it)",core,1L)
 expect_source_count("#pragma omp parallel for num_threads(ncores) schedule(static)",adapter,1L)
 expect_false(grepl("run_one_bayesrc_chain|ChainResultBayesRC|old_path|new_path|fallback|execution_selector|rng_selector",paste(core,adapter,types)))
})

test_that("Phase 14C numerical headers are binding neutral",{
 types<-phase14c_text("src/blr_bed_bayesrc_types.h")
 core<-phase14c_text("src/blr_bed_bayesrc_core_impl.h")
 probability<-phase14c_text("src/blr_normal_probability.h")
 expect_false(grepl("Rcpp|SEXP|RObject|NumericVector|NumericMatrix|CharacterVector|DataFrame|pybind11|Python",paste(types,core)))
 expect_false(grepl("R::pnorm|R::qnorm",core,fixed=FALSE))
 expect_source_count("struct StandardNormalProbability",probability,1L)
 expect_source_count("static inline double cdf(double x)",probability,1L)
 expect_source_count("static inline double quantile(double p)",probability,1L)
 expect_match(probability,"R::pnorm(x, 0.0, 1.0, 1, 0)",fixed=TRUE)
 expect_match(probability,"R::qnorm(p, 0.0, 1.0, 1, 0)",fixed=TRUE)
})

test_that("Phase 14C preserves probability latent alpha and full-sweep operations",{
 core<-phase14c_text("src/blr_bed_bayesrc_core_impl.h")
 helper<-phase14c_text("src/st_bayesrc_annotation_prior.h")
 for(x in c("std::max(val, pi_floor)","snp_pi.row","/= s",
  "st_bayesrc_sample_truncated_normal_std","(ci > j) ? 1 : 0",
  "std::normal_distribution<double> norm(mean, sd)",
  "std::chi_squared_distribution<double> rchisq(df)")) expect_match(helper,x,fixed=TRUE)
 for(x in c("std::mt19937 gen(static_cast<unsigned int>(context.chain_seed))",
  "std::uniform_real_distribution<double> runif(0.0, 1.0)",
  "std::normal_distribution<double> norm01(0.0, 1.0)",
  "for (int marker : marker_order)","const double vbk = vb * gamma",
  "component_new = ncomponent - 1","if (component_new > 0)")) expect_match(core,x,fixed=TRUE)
 expect_false(grepl("thread_local|omp_get_thread_num|full_sweep_every|null_skip|candidate_list|active_list|due_bucket|scheduler",paste(types<-phase14c_text("src/blr_bed_bayesrc_types.h"),core)))
})

test_that("Phase 14C typed chain boundary remains active below aggregation",{
 types<-phase14c_text("src/blr_bed_bayesrc_types.h")
 core<-phase14c_text("src/blr_bed_bayesrc_core_impl.h")
 adapter<-phase14c_text("src/stblr_cpg_omp_bed_marker_scheduled_chains_bayesrc.cpp")
 expect_match(types,"const PackedGenotype& storage",fixed=TRUE)
 expect_match(types,"const AnnotationMatrix& matrix",fixed=TRUE)
 expect_match(adapter,"br_read_bed_blocked",fixed=TRUE)
 expect_match(adapter,"BedBayesRCChainExecutionContext<",fixed=TRUE)
 expect_match(adapter,"aggregate_bed_bayesrc_results(jobs,aggregation_context)",fixed=TRUE)
 expect_match(adapter,"stblr_bed_bayesrc_result_to_raw(result,metadata)",fixed=TRUE)
 expect_false(grepl("fopen|fseek|fread|bed_files|Rcpp::List|R_NilValue|marker_id|factor|add_intercept",core))
})

test_that("Phase 14C references and reproducibility remain exact",{
 cfg<-list(one_chain_one_core=c(1,1,141,1,0),two_chains_one_core=c(1,2,143,1,1),two_chains_two_cores=c(2,2,143,1,1))
 for(nm in names(cfg)){z<-cfg[[nm]];ref<-readRDS(file.path(testthat::test_path(),"fixtures","blr_phase14a_bed_bayesrc",paste0(nm,".rds")))
  observed<-phase14a_normalize(phase14a_capture(z[1],z[2],z[3],as.logical(z[4]),as.logical(z[5])))
  expect_identical(observed$raw,ref$raw);expect_identical(observed$fit,ref$fit)}
 a<-phase14a_normalize(phase14a_capture(1,2,143,TRUE,TRUE))
 expect_identical(phase14a_normalize(phase14a_capture(1,2,143,TRUE,TRUE)),a)
 phase14a_capture(1,1,149,FALSE,FALSE)
 expect_identical(phase14a_normalize(phase14a_capture(2,2,143,TRUE,TRUE)),a)
 phase11a_capture("bayesr",1,1,71);phase11a_capture("bayesc",1,1,71)
 expect_identical(phase14a_normalize(phase14a_capture(1,2,143,TRUE,TRUE)),a)
})

test_that("Phase 14C identities and reductions remain exact",{
 z<-phase14a_capture(1,1,141,FALSE,FALSE)$raw;cp<-z$component$prob[[1L]]
 expect_true(all(is.finite(cp)&cp>=0));expect_equal(rowSums(cp),rep(1,nrow(cp)),tolerance=1e-12)
 expect_equal(z$marker$dm[,1],1-cp[,1],tolerance=1e-12)
 expect_true(all(z$marker$state[,1]>=0&z$marker$state[,1]<ncol(cp)))
 expect_equal(dim(z$annotation$alpha_final[[1L]]),c(1L,3L))
 expect_equal(z$pi$final[1,],c(.95,.03,.015,.005),tolerance=1e-12)
 expect_equal(unname(z$annotation$alpha_final[[1]]),unname(sblr:::.bayesr_pi_to_probit_stick_intercepts(c(.95,.03,.015,.005))),tolerance=1e-12)
 expect_match(phase14c_text("tests/testthat/test-individual-bayesrc.R"),"reduces to fixed-pi BayesR",fixed=TRUE)
})

test_that("Phase 14C protects public and adjacent backend boundaries",{
 protected<-c("src/blr_bed_bayesr_core_impl.h"="afe77e26d2cf2b8e3d64088221b33e14",
  "src/blr_bed_scheduled_bayesc_core_impl.h"="723cee003504c1fdcd075b965cb63d83",
  "src/blr_csr_sbayesrc_core_impl.h"="d06ec2a530e8c914201ee22b6be65739",
  "src/st_block_eigen.cpp"="49f0a62c9fe235967a264b0f8de144a7",
  "src/mt_cpg_omp_csr.cpp"="aec85896b5c30db3014efaeb5e3c3a96",
  "R/RcppExports.R"="9d13ea00b326c7e0cd606194d13a8bca",
  "src/RcppExports.cpp"="b4859db0f6308fa7e38051ddcf32d245",
  "NAMESPACE"="f5b6ee37a3972aa436357bdc8f602f4e")
 expect_identical(unname(tools::md5sum(file.path(phase14c_root,names(protected)))),unname(protected))
})

test_that("Phase 14C fresh-process references can be checked explicitly",{
 skip_if_not(identical(Sys.getenv("SBLR_RUN_PHASE14A_FRESH"),"true"));skip_if_not_installed("callr")
 observed<-callr::r(function(root){setwd(root);pkgload::load_all(".",compile=FALSE,quiet=TRUE);library(testthat);source("tests/testthat/fixtures/blr-phase14a-bed-bayesrc-reference.R");phase14a_normalize(phase14a_capture(2,2,143,TRUE,TRUE))},list(root=phase14c_root))
 ref<-readRDS(file.path(testthat::test_path(),"fixtures","blr_phase14a_bed_bayesrc","two_chains_two_cores.rds"))
 expect_identical(observed$raw,ref$raw);expect_identical(observed$fit,ref$fit)
})
