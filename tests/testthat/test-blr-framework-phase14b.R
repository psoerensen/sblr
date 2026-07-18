phase14b_root <- normalizePath(file.path(testthat::test_path(),"..",".."),winslash="/",mustWork=TRUE)
source(file.path(testthat::test_path(),"fixtures","blr-phase14a-bed-bayesrc-reference.R"))
phase14b_text <- function(path) paste(readLines(file.path(phase14b_root,path),warn=FALSE),collapse="\n")

test_that("Phase 14B has one guarded extracted chain implementation",{
 adapter<-phase14b_text("src/stblr_cpg_omp_bed_marker_scheduled_chains_bayesrc.cpp")
 core<-phase14b_text("src/blr_bed_bayesrc_core_impl.h")
 expect_match(core,"#ifndef SBLR_BLR_BED_BAYESRC_CORE_IMPL_H",fixed=TRUE)
 expect_source_count('#include "blr_bed_bayesrc_core_impl.h"',adapter,1L)
 expect_source_count("static BedBayesRCChainExecutionResult run_bed_bayesrc_chain(",core,1L)
 expect_source_count("for (int it = 0; it < total_it; ++it)",core,1L)
 expect_source_count("run_bed_bayesrc_chain(",paste(adapter,core),2L)
 expect_source_count("#pragma omp parallel for num_threads(ncores) schedule(static)",adapter,1L)
 expect_false(grepl("run_one_bayesrc_chain|for \\(int it",adapter))
 expect_false(grepl("old_path|new_path|fallback|execution_selector|rng_selector",paste(adapter,core)))
})

test_that("Phase 14B preserves full sweeps, component sampling and RNG",{
 core<-phase14b_text("src/blr_bed_bayesrc_core_impl.h")
 for(x in c("std::mt19937 gen(static_cast<unsigned int>(context.chain_seed))",
  "std::uniform_real_distribution<double> runif(0.0, 1.0)","std::normal_distribution<double> norm01(0.0, 1.0)",
  "for (int marker : marker_order)","const double vbk = vb * gamma","component_new = ncomponent - 1",
  "if (component_new > 0)","component_t(ju) > 0 ? 1.0 : 0.0")) expect_match(core,x,fixed=TRUE)
 expect_false(grepl("thread_local|omp_get_thread_num|candidate_list|scheduled_at|null_skip|due_bucket|full_sweep_every",core))
})

test_that("Phase 14B preserves the shared probit, latent and alpha boundary",{
 helper<-phase14b_text("src/st_bayesrc_annotation_prior.h")
 core<-phase14b_text("src/blr_bed_bayesrc_core_impl.h")
 for(x in c("StandardNormalProbability::cdf","StandardNormalProbability::quantile","std::max(val, pi_floor)","snp_pi.row","/= s",
  "st_bayesrc_sample_truncated_normal_std","(ci > j) ? 1 : 0","std::normal_distribution<double> norm(mean, sd)",
  "std::chi_squared_distribution<double> rchisq(df)")) expect_match(helper,x,fixed=TRUE)
 expect_source_count("st_bayesrc_compute_snp_pi(annotation, annot_alpha, pi_floor)",core,2L)
 expect_source_count("st_bayesrc_update_annotation_prior(",core,1L)
 expect_false(grepl("factor|marker_id|rownames|colnames|add_intercept|data.frame",core))
})

test_that("Phase 14B leaves decoding, dispatch, aggregation and conversion adapter-owned",{
 adapter<-phase14b_text("src/stblr_cpg_omp_bed_marker_scheduled_chains_bayesrc.cpp")
 core<-phase14b_text("src/blr_bed_bayesrc_core_impl.h")
 for(x in c("br_read_bed_blocked","std::vector<sblr::core::BedBayesRCChainExecutionResult> jobs","const sblr::core::BedBayesRCChainExecutionResult& z = jobs[ch * nt + t]",
  "Rcpp::List raw = Rcpp::List::create","marker_prior_final[t] += st_bayesrc_compute_snp_pi")) expect_match(adapter,x,fixed=TRUE)
 expect_false(grepl("Rcpp::List|br_read_bed_blocked|fopen|fseek|fread|bed_files|R_NilValue",core))
})

test_that("Phase 14B references and reproducibility remain exact",{
 cfg<-list(one_chain_one_core=c(1,1,141,1,0),two_chains_one_core=c(1,2,143,1,1),two_chains_two_cores=c(2,2,143,1,1))
 for(nm in names(cfg)){z<-cfg[[nm]];ref<-readRDS(file.path(testthat::test_path(),"fixtures","blr_phase14a_bed_bayesrc",paste0(nm,".rds")))
  observed<-phase14a_normalize(phase14a_capture(z[1],z[2],z[3],as.logical(z[4]),as.logical(z[5])))
  expect_identical(observed$raw,ref$raw);expect_identical(observed$fit,ref$fit)}
 a<-phase14a_normalize(phase14a_capture(1,2,143,TRUE,TRUE))
 expect_identical(phase14a_normalize(phase14a_capture(2,2,143,TRUE,TRUE)),a)
 phase11a_capture("bayesr",1,1,71);phase11a_capture("bayesc",1,1,71)
 expect_identical(phase14a_normalize(phase14a_capture(1,2,143,TRUE,TRUE)),a)
})

test_that("Phase 14B identities and fixed-alpha reduction remain protected",{
 z<-phase14a_capture(1,1,141,FALSE,FALSE)$raw;cp<-z$component$prob[[1L]]
 expect_equal(rowSums(cp),rep(1,nrow(cp)),tolerance=1e-12)
 expect_equal(z$marker$dm[,1],1-cp[,1],tolerance=1e-12)
 expect_true(all(z$marker$state[,1]>=0&z$marker$state[,1]<ncol(cp)))
 target<-c(.95,.03,.015,.005)
 expect_equal(z$pi$final[1,],target,tolerance=1e-12)
 expect_equal(unname(z$annotation$alpha_final[[1]]),unname(sblr:::.bayesr_pi_to_probit_stick_intercepts(target)),tolerance=1e-12)
 expect_match(phase14b_text("tests/testthat/test-individual-bayesrc.R"),"reduces to fixed-pi BayesR",fixed=TRUE)
})

test_that("Phase 14B protects unrelated backends and generated interfaces",{
 protected<-c("src/blr_bed_bayesr_core_impl.h"="afe77e26d2cf2b8e3d64088221b33e14",
 "src/blr_bed_scheduled_bayesc_core_impl.h"="723cee003504c1fdcd075b965cb63d83",
 "src/blr_csr_sbayesrc_core_impl.h"="d06ec2a530e8c914201ee22b6be65739",
 "src/st_bayesrc_annotation_prior.h"="1e7072512f4246fc2a36e79de655d8c5",
 "R/RcppExports.R"="9d13ea00b326c7e0cd606194d13a8bca","src/RcppExports.cpp"="b4859db0f6308fa7e38051ddcf32d245","NAMESPACE"="f5b6ee37a3972aa436357bdc8f602f4e")
 expect_identical(unname(tools::md5sum(file.path(phase14b_root,names(protected)))),unname(protected))
})
