phase14e_root<-normalizePath(file.path(testthat::test_path(),"..",".."),winslash="/",mustWork=TRUE)
source(file.path(testthat::test_path(),"fixtures","blr-phase14a-bed-bayesrc-reference.R"))
phase14e_text<-function(path) paste(readLines(file.path(phase14e_root,path),warn=FALSE),collapse="\n")

test_that("canonical packed-BED BayesRC has one permanent architecture",{
 types<-phase14e_text("src/blr_bed_bayesrc_types.h");core<-phase14e_text("src/blr_bed_bayesrc_core_impl.h")
 aggregate<-phase14e_text("src/blr_bed_bayesrc_aggregate_impl.h");adapter<-phase14e_text("src/stblr_cpg_omp_bed_marker_scheduled_chains_bayesrc.cpp")
 for(x in c("BedBayesRCComponentSpec","BedBayesRCCoefficientPriorSpec","BedBayesRCChainExecutionResult","BedBayesRCAggregationContext","BedBayesRCExecutionResult"))
  expect_source_count(paste0("struct ",x),types,1L)
 for(x in c("BedBayesRCAnnotationSpec","BedBayesRCPackedGenotypeView","BedBayesRCChainExecutionContext"))
  expect_source_count(paste0("struct ",x),types,1L)
 expect_source_count("run_bed_bayesrc_chain(",paste(core,adapter),2L)
 expect_source_count("aggregate_bed_bayesrc_results(",paste(aggregate,adapter),2L)
 expect_source_count("stblr_bed_bayesrc_result_to_raw(",adapter,2L)
 expect_source_count("for (int it = 0; it < total_it; ++it)",core,1L)
 expect_source_count("#pragma omp parallel for num_threads(ncores) schedule(static)",adapter,1L)
 expect_false(grepl("run_one_bayesrc_chain|ChainResultBayesRC|old_path|new_path|fallback|selector|aggregate2|inline_result",paste(types,core,aggregate,adapter)))
})

test_that("canonical numerical boundaries are binding neutral and full-sweep",{
 types<-phase14e_text("src/blr_bed_bayesrc_types.h");core<-phase14e_text("src/blr_bed_bayesrc_core_impl.h")
 aggregate<-phase14e_text("src/blr_bed_bayesrc_aggregate_impl.h");probability<-phase14e_text("src/blr_normal_probability.h")
 expect_false(grepl("Rcpp|SEXP|RObject|NumericVector|NumericMatrix|DataFrame|pybind11|Python",paste(types,core,aggregate)))
 expect_false(grepl("R::pnorm|R::qnorm",paste(core,aggregate)))
 expect_match(probability,"R::pnorm(x, 0.0, 1.0, 1, 0)",fixed=TRUE)
 expect_match(probability,"R::qnorm(p, 0.0, 1.0, 1, 0)",fixed=TRUE)
 expect_match(core,"for (int marker : marker_order)",fixed=TRUE)
 expect_false(grepl("full_sweep_every|null_skip|candidate_list|active_list|due_bucket|scheduler",paste(types,core,aggregate)))
 expect_false(grepl("fopen|fseek|fread|ifstream|bed_files",paste(core,aggregate)))
})

test_that("canonical probability alpha RNG and ownership policies remain",{
 types<-phase14e_text("src/blr_bed_bayesrc_types.h");core<-phase14e_text("src/blr_bed_bayesrc_core_impl.h")
 aggregate<-phase14e_text("src/blr_bed_bayesrc_aggregate_impl.h");helper<-phase14e_text("src/st_bayesrc_annotation_prior.h")
 for(x in c("const PackedGenotype& storage","const AnnotationMatrix& matrix")) expect_match(types,x,fixed=TRUE)
 for(x in c("std::mt19937 gen(static_cast<unsigned int>(context.chain_seed))","std::uniform_real_distribution<double> runif(0.0, 1.0)","std::normal_distribution<double> norm01(0.0, 1.0)")) expect_match(core,x,fixed=TRUE)
 for(x in c("std::max(val, pi_floor)","prod_prev *= pk","snp_pi.row","/= s","st_bayesrc_sample_truncated_normal_std","std::normal_distribution<double> norm(mean, sd)","std::chi_squared_distribution<double> rchisq(df)")) expect_match(helper,x,fixed=TRUE)
 expect_source_count("st_bayesrc_compute_snp_pi(",aggregate,1L)
 expect_false(grepl("thread_local|omp_get_thread_num|std::mt19937|uniform_real_distribution|normal_distribution",aggregate))
})

test_that("Phase 14A fixtures are permanent canonical references",{
 cfg<-list(one_chain_one_core=c(1,1,141,1,0),two_chains_one_core=c(1,2,143,1,1),two_chains_two_cores=c(2,2,143,1,1))
 for(nm in names(cfg)){z<-cfg[[nm]];ref<-readRDS(file.path(testthat::test_path(),"fixtures","blr_phase14a_bed_bayesrc",paste0(nm,".rds")))
  observed<-phase14a_normalize(phase14a_capture(z[1],z[2],z[3],as.logical(z[4]),as.logical(z[5])))
  expect_identical(observed$raw,ref$raw);expect_identical(observed$fit,ref$fit)}
})

test_that("canonical reproducibility identities and reduction remain exact",{
 a<-phase14a_normalize(phase14a_capture(1,2,143,TRUE,TRUE))
 expect_identical(phase14a_normalize(phase14a_capture(1,2,143,TRUE,TRUE)),a)
 phase14a_capture(1,1,149,FALSE,FALSE)
 expect_identical(phase14a_normalize(phase14a_capture(2,2,143,TRUE,TRUE)),a)
 phase11a_capture("bayesr",1,1,71);phase11a_capture("bayesc",1,1,71)
 expect_identical(phase14a_normalize(phase14a_capture(1,2,143,TRUE,TRUE)),a)
 z<-phase14a_capture(1,1,141,FALSE,FALSE)$raw;cp<-z$component$prob[[1L]]
 expect_true(all(is.finite(cp)&cp>=0));expect_equal(rowSums(cp),rep(1,nrow(cp)),tolerance=1e-12)
 expect_equal(z$marker$dm[,1],1-cp[,1],tolerance=1e-12)
 expect_true(all(z$marker$state[,1]>=0&z$marker$state[,1]<ncol(cp)))
 expect_equal(dim(z$annotation$alpha_final[[1L]]),c(1L,3L))
 expect_equal(z$pi$final[1,],c(.95,.03,.015,.005),tolerance=1e-12)
 expect_match(phase14e_text("tests/testthat/test-individual-bayesrc.R"),"reduces to fixed-pi BayesR",fixed=TRUE)
})

test_that("canonical route and unsupported policies remain explicit",{
 adapter<-phase14e_text("src/stblr_cpg_omp_bed_marker_scheduled_chains_bayesrc.cpp")
 expect_match(adapter,"aggregate_bed_bayesrc_results(jobs,aggregation_context)",fixed=TRUE)
 expect_match(adapter,"br_xb(G,maps,order,result.b.col(t).t())",fixed=TRUE)
 expect_match(adapter,"stblr_bed_bayesrc_result_to_raw(result,metadata)",fixed=TRUE)
 rsrc<-phase14e_text("R/sparse_ld_bed_helper.R")
 expect_match(rsrc,'method == "bayesrc"',fixed=TRUE)
 expect_false(grepl("chain_seeds|full_sweep_every|null_skip_base|candidate_threshold",adapter))
})

test_that("canonical BayesRC protects other backends and interfaces",{
 protected<-c("src/blr_bed_bayesr_core_impl.h"="afe77e26d2cf2b8e3d64088221b33e14",
  "src/blr_bed_scheduled_bayesc_core_impl.h"="723cee003504c1fdcd075b965cb63d83",
  "src/blr_csr_sbayesrc_core_impl.h"="d06ec2a530e8c914201ee22b6be65739",
  "src/st_block_eigen.cpp"="49f0a62c9fe235967a264b0f8de144a7","NAMESPACE"="ab1479ce78ea20b39bf8b94f9bc0aa62")
 expect_identical(unname(tools::md5sum(file.path(phase14e_root,names(protected)))),unname(protected))
})
