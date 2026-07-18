phase14d_root<-normalizePath(file.path(testthat::test_path(),"..",".."),winslash="/",mustWork=TRUE)
source(file.path(testthat::test_path(),"fixtures","blr-phase14a-bed-bayesrc-reference.R"))
phase14d_text<-function(path) paste(readLines(file.path(phase14d_root,path),warn=FALSE),collapse="\n")

test_that("Phase 14D has one typed aggregate and one aggregation path",{
 types<-phase14d_text("src/blr_bed_bayesrc_types.h")
 aggregate<-phase14d_text("src/blr_bed_bayesrc_aggregate_impl.h")
 adapter<-phase14d_text("src/stblr_cpg_omp_bed_marker_scheduled_chains_bayesrc.cpp")
 expect_source_count("struct BedBayesRCExecutionResult",types,1L)
 expect_source_count("struct BedBayesRCAggregationContext",types,1L)
 expect_source_count("aggregate_bed_bayesrc_results(",paste(aggregate,adapter),2L)
 expect_source_count("stblr_bed_bayesrc_result_to_raw(",adapter,2L)
 expect_source_count("#pragma omp parallel for num_threads(ncores) schedule(static)",adapter,1L)
 expect_match(adapter,"aggregate_bed_bayesrc_results(jobs,aggregation_context)",fixed=TRUE)
 expect_match(adapter,"stblr_bed_bayesrc_result_to_raw(result,metadata)",fixed=TRUE)
 expect_false(grepl("aggregation_selector|conversion_selector|old_path|new_path|fallback",paste(types,aggregate,adapter)))
})

test_that("Phase 14D aggregation is binding neutral and singular",{
 aggregate<-phase14d_text("src/blr_bed_bayesrc_aggregate_impl.h")
 expect_false(grepl("Rcpp|SEXP|RObject|NumericVector|NumericMatrix|CharacterVector|DataFrame|pybind11|Python",aggregate))
 expect_false(grepl("R::pnorm|R::qnorm|std::mt19937|uniform_real_distribution|normal_distribution",aggregate))
 expect_false(grepl("fopen|fseek|fread|ifstream|bed_files|br_read_bed",aggregate))
 for(x in c("out.bm.col(t)+=z.bm.t()","out.comp_prob[t]+=z.comp_prob",
  "out.alpha_final[t]+=z.annot_alpha_final","out.log_cpo(t)+=z.log_cpo",
  "out.component_counts.row(t)=arma::sum(out.comp_prob[t],0)")) expect_match(aggregate,x,fixed=TRUE)
})

test_that("Phase 14D final marker priors use the canonical probability helper",{
 aggregate<-phase14d_text("src/blr_bed_bayesrc_aggregate_impl.h")
 helper<-phase14d_text("src/st_bayesrc_annotation_prior.h")
 expect_source_count("st_bayesrc_compute_snp_pi(",aggregate,1L)
 expect_source_count("st_bayesrc_compute_snp_pi(",helper,1L)
 expect_match(aggregate,"context.annotation,z.annot_alpha_final,context.pi_floor",fixed=TRUE)
 for(x in c("StandardNormalProbability::cdf","std::max(val, pi_floor)","prod_prev *= pk","/= s"))
  expect_match(helper,x,fixed=TRUE)
})

test_that("Phase 14D keeps dispatch and optional genotype diagnostics in the adapter",{
 adapter<-phase14d_text("src/stblr_cpg_omp_bed_marker_scheduled_chains_bayesrc.cpp")
 aggregate<-phase14d_text("src/blr_bed_bayesrc_aggregate_impl.h")
 for(x in c("br_read_bed_blocked","BedBayesRCChainExecutionContext<","br_xb(G,maps,order,result.b.col(t).t())",
  "br_dot_residual","BedBayesRCBindingMetadata")) expect_match(adapter,x,fixed=TRUE)
 expect_false(grepl("br_xb|br_dot_residual|FastPackedBedMatrixBR|marker_maps",aggregate))
})

test_that("Phase 14D references remain exactly equal",{
 cfg<-list(one_chain_one_core=c(1,1,141,1,0),two_chains_one_core=c(1,2,143,1,1),two_chains_two_cores=c(2,2,143,1,1))
 for(nm in names(cfg)){z<-cfg[[nm]];ref<-readRDS(file.path(testthat::test_path(),"fixtures","blr_phase14a_bed_bayesrc",paste0(nm,".rds")))
  observed<-phase14a_normalize(phase14a_capture(z[1],z[2],z[3],as.logical(z[4]),as.logical(z[5])))
  expect_identical(observed$raw,ref$raw);expect_identical(observed$fit,ref$fit)}
})

test_that("Phase 14D reproducibility identities and reduction remain exact",{
 a<-phase14a_normalize(phase14a_capture(1,2,143,TRUE,TRUE))
 expect_identical(phase14a_normalize(phase14a_capture(1,2,143,TRUE,TRUE)),a)
 phase14a_capture(1,1,149,FALSE,FALSE)
 expect_identical(phase14a_normalize(phase14a_capture(2,2,143,TRUE,TRUE)),a)
 phase11a_capture("bayesr",1,1,71);phase11a_capture("bayesc",1,1,71)
 expect_identical(phase14a_normalize(phase14a_capture(1,2,143,TRUE,TRUE)),a)
 z<-phase14a_capture(1,1,141,FALSE,FALSE)$raw;cp<-z$component$prob[[1L]]
 expect_equal(rowSums(cp),rep(1,nrow(cp)),tolerance=1e-12)
 expect_equal(z$marker$dm[,1],1-cp[,1],tolerance=1e-12)
 expect_true(all(z$marker$state[,1]>=0&z$marker$state[,1]<ncol(cp)))
 expect_equal(z$pi$final[1,],c(.95,.03,.015,.005),tolerance=1e-12)
 expect_match(phase14d_text("tests/testthat/test-individual-bayesrc.R"),"reduces to fixed-pi BayesR",fixed=TRUE)
})

test_that("Phase 14D protects adjacent backends and generated interfaces",{
 protected<-c("src/blr_bed_bayesr_core_impl.h"="afe77e26d2cf2b8e3d64088221b33e14",
  "src/blr_bed_scheduled_bayesc_core_impl.h"="723cee003504c1fdcd075b965cb63d83",
  "src/blr_csr_sbayesrc_core_impl.h"="d06ec2a530e8c914201ee22b6be65739",
  "src/st_block_eigen.cpp"="49f0a62c9fe235967a264b0f8de144a7",
  "src/mt_cpg_omp_csr.cpp"="aec85896b5c30db3014efaeb5e3c3a96",
  "R/RcppExports.R"="9d13ea00b326c7e0cd606194d13a8bca",
  "src/RcppExports.cpp"="b4859db0f6308fa7e38051ddcf32d245","NAMESPACE"="f5b6ee37a3972aa436357bdc8f602f4e")
 expect_identical(unname(tools::md5sum(file.path(phase14d_root,names(protected)))),unname(protected))
})

test_that("Phase 14D fresh-process reference is opt-in",{
 skip_if_not(identical(Sys.getenv("SBLR_RUN_PHASE14A_FRESH"),"true"));skip_if_not_installed("callr")
 observed<-callr::r(function(root){setwd(root);pkgload::load_all(".",compile=FALSE,quiet=TRUE);library(testthat);source("tests/testthat/fixtures/blr-phase14a-bed-bayesrc-reference.R");phase14a_normalize(phase14a_capture(2,2,143,TRUE,TRUE))},list(root=phase14d_root))
 ref<-readRDS(file.path(testthat::test_path(),"fixtures","blr_phase14a_bed_bayesrc","two_chains_two_cores.rds"))
 expect_identical(observed$raw,ref$raw);expect_identical(observed$fit,ref$fit)
})
