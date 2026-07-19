phase15a_root<-normalizePath(file.path(testthat::test_path(),"..",".."),winslash="/",mustWork=TRUE)
phase15a_text<-function(path) paste(readLines(file.path(phase15a_root,path),warn=FALSE),collapse="\n")

test_that("canonical packed-BED families retain one closed architecture each",{
 models<-list(
  bayesc=c("src/blr_bed_scheduled_bayesc_types.h","src/blr_bed_scheduled_bayesc_core_impl.h","src/blr_bed_scheduled_bayesc_aggregate_impl.h","src/st_cpg_omp_individual_scheduled_chains.cpp"),
  bayesr=c("src/blr_bed_bayesr_types.h","src/blr_bed_bayesr_core_impl.h","src/blr_bed_bayesr_aggregate_impl.h","src/stblr_cpg_omp_bed_marker_scheduled_chains_bayesr.cpp"),
  bayesrc=c("src/blr_bed_bayesrc_types.h","src/blr_bed_bayesrc_core_impl.h","src/blr_bed_bayesrc_aggregate_impl.h","src/stblr_cpg_omp_bed_marker_scheduled_chains_bayesrc.cpp"))
 needles<-list(
  bayesc=c("BedScheduledBayesCChainExecutionContext","BedScheduledBayesCChainExecutionResult","run_bed_scheduled_bayesc_chain","BedScheduledBayesCExecutionResult","aggregate_bed_scheduled_bayesc_results","stblr_bed_scheduled_bayesc_result_to_raw"),
  bayesr=c("BedBayesRChainExecutionContext","BedBayesRChainExecutionResult","run_bed_bayesr_chain","BedBayesRExecutionResult","aggregate_bed_bayesr_results","stblr_bed_bayesr_result_to_raw"),
  bayesrc=c("BedBayesRCChainExecutionContext","BedBayesRCChainExecutionResult","run_bed_bayesrc_chain","BedBayesRCExecutionResult","aggregate_bed_bayesrc_results","stblr_bed_bayesrc_result_to_raw"))
 for(nm in names(models)){x<-paste(vapply(models[[nm]],phase15a_text,character(1)),collapse="\n");for(s in needles[[nm]])expect_match(x,s,fixed=TRUE);expect_match(x,"schedule(static)",fixed=TRUE)}
})

test_that("task mapping and logical-chain seed policy are exactly equivalent",{
 adapters<-vapply(c("src/st_cpg_omp_individual_scheduled_chains.cpp","src/stblr_cpg_omp_bed_marker_scheduled_chains_bayesr.cpp","src/stblr_cpg_omp_bed_marker_scheduled_chains_bayesrc.cpp"),phase15a_text,character(1))
 common<-phase15a_text("src/blr_bed_family_types.h")
 expect_match(common,"job / trait_count",fixed=TRUE);expect_match(common,"job % trait_count",fixed=TRUE)
 expect_match(common,"1000003",fixed=TRUE);expect_match(common,"9176",fixed=TRUE)
 for(x in adapters){expect_match(x,"make_bed_family_task_index",fixed=TRUE);expect_match(x,"resolve_bed_family_logical_chain_seed",fixed=TRUE);expect_match(x,"schedule(static)",fixed=TRUE)}
 expect_false(any(grepl("omp_get_thread_num.*seed|seed.*omp_get_thread_num",adapters)))
})

test_that("genotype ownership is common but representation sharing is bounded",{
 common<-phase15a_text("src/blr_bed_family_types.h");ctype<-phase15a_text("src/blr_bed_scheduled_bayesc_types.h");rtype<-phase15a_text("src/blr_bed_bayesr_types.h");rctype<-phase15a_text("src/blr_bed_bayesrc_types.h")
 for(x in c(common,rctype))expect_match(x,"const PackedGenotype& storage",fixed=TRUE)
 for(s in c("packed_markers","packed_size","bytes_per_marker","stride"))expect_match(common,s,fixed=TRUE)
 expect_match(ctype,"blr_bed_family_types.h",fixed=TRUE)
 expect_match(rtype,"using BedBayesRPackedGenotypeView",fixed=TRUE)
 expect_false(grepl("packed_markers|packed_size|stride",rctype))
 cores<-paste(vapply(c("src/blr_bed_scheduled_bayesc_core_impl.h","src/blr_bed_bayesr_core_impl.h","src/blr_bed_bayesrc_core_impl.h"),phase15a_text,character(1)),collapse="\n")
 expect_false(grepl("fopen|fseek|fread|ifstream|bed_files|memory_map",cores))
})

test_that("reuse classifications are explicit and model-specific semantics stay separate",{
 report<-phase15a_text("docs/dev/blr_framework_phase15a_report.md")
 for(x in c("reuse unchanged","reuse through narrow common wrapper","reuse helper only","reuse convention only","retain model-specific","unsafe to consolidate"))expect_match(report,x,fixed=TRUE)
 ccore<-phase15a_text("src/blr_bed_scheduled_bayesc_core_impl.h");rcore<-phase15a_text("src/blr_bed_bayesr_core_impl.h");rccore<-phase15a_text("src/blr_bed_bayesrc_core_impl.h")
 expect_match(ccore,"candidate",fixed=TRUE);expect_match(rcore,"component",fixed=TRUE)
 expect_match(rccore,"st_bayesrc_compute_snp_pi",fixed=TRUE);expect_match(rccore,"for (int marker : marker_order)",fixed=TRUE)
 expect_false(grepl("null_skip|candidate_list|due_bucket",rccore))
})

test_that("Phase 15B changes only audited packed-BED family infrastructure",{
 protected<-c("src/blr_bed_scheduled_bayesc_types.h"="3d7673842ca20656551621ef2e6336b3","src/blr_bed_scheduled_bayesc_core_impl.h"="723cee003504c1fdcd075b965cb63d83","src/blr_bed_scheduled_bayesc_aggregate_impl.h"="c3023c3af5e62d83599af1fab1aa9fa3","src/st_cpg_omp_individual_scheduled_chains.cpp"="aae0fd8bc851cc96f946309eae81defe","src/blr_bed_bayesr_types.h"="67a07af56e9331f8d89cd99843304877","src/blr_bed_bayesr_core_impl.h"="afe77e26d2cf2b8e3d64088221b33e14","src/blr_bed_bayesr_aggregate_impl.h"="8d0749f385319fbbf93e55317747642b","src/stblr_cpg_omp_bed_marker_scheduled_chains_bayesr.cpp"="6f9f73843a7a4d67d796931890f62f1b","src/blr_bed_bayesrc_types.h"="cee9c94df0672674865ec97fca212636","src/blr_bed_bayesrc_core_impl.h"="82365cf3f1f5306c57b980f59b4d83d3","src/blr_bed_bayesrc_aggregate_impl.h"="9ef26e7e1d587149f8801f7fedfed162","src/stblr_cpg_omp_bed_marker_scheduled_chains_bayesrc.cpp"="72d4a9fa0a7cd51071328c2d62d0192b","R/RcppExports.R"="9d13ea00b326c7e0cd606194d13a8bca","src/RcppExports.cpp"="b4859db0f6308fa7e38051ddcf32d245","NAMESPACE"="f5b6ee37a3972aa436357bdc8f602f4e")
 expect_identical(unname(tools::md5sum(file.path(phase15a_root,names(protected)))),unname(protected))
})

test_that("public schema conventions remain common without merging converters",{
 adapters<-vapply(c("src/st_cpg_omp_individual_scheduled_chains.cpp","src/stblr_cpg_omp_bed_marker_scheduled_chains_bayesr.cpp","src/stblr_cpg_omp_bed_marker_scheduled_chains_bayesrc.cpp"),phase15a_text,character(1))
 for(x in adapters){expect_match(x,"stblr_raw_v1",fixed=TRUE);expect_match(x,"R_NilValue",fixed=TRUE)}
 expect_true(length(unique(grep("result_to_raw",adapters,value=TRUE)))==3L)
})
