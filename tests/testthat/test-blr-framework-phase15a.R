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
 for(x in adapters){expect_match(x,"job\\s*/\\s*nt");expect_match(x,"job\\s*%\\s*nt");expect_match(x,"1000003",fixed=TRUE);expect_match(x,"9176",fixed=TRUE);expect_match(x,"schedule(static)",fixed=TRUE)}
 expect_false(any(grepl("omp_get_thread_num.*seed|seed.*omp_get_thread_num",adapters)))
})

test_that("genotype ownership is common but representation sharing is bounded",{
 ctype<-phase15a_text("src/blr_bed_scheduled_bayesc_types.h");rtype<-phase15a_text("src/blr_bed_bayesr_types.h");rctype<-phase15a_text("src/blr_bed_bayesrc_types.h")
 for(x in c(ctype,rtype,rctype))expect_match(x,"const PackedGenotype& storage",fixed=TRUE)
 for(x in c(ctype,rtype))for(s in c("packed_markers","packed_size","bytes_per_marker","stride"))expect_match(x,s,fixed=TRUE)
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

test_that("canonical production sources are byte-identical to the Phase 15A baseline",{
 protected<-c("src/blr_bed_scheduled_bayesc_types.h"="40000fea66bb8ea95183151e83c62c87","src/blr_bed_scheduled_bayesc_core_impl.h"="723cee003504c1fdcd075b965cb63d83","src/blr_bed_scheduled_bayesc_aggregate_impl.h"="c3023c3af5e62d83599af1fab1aa9fa3","src/st_cpg_omp_individual_scheduled_chains.cpp"="43c71b13d8259a95f88d8a95498b213b","src/blr_bed_bayesr_types.h"="137ce5391f885df14fc738655e5a065e","src/blr_bed_bayesr_core_impl.h"="afe77e26d2cf2b8e3d64088221b33e14","src/blr_bed_bayesr_aggregate_impl.h"="8d0749f385319fbbf93e55317747642b","src/stblr_cpg_omp_bed_marker_scheduled_chains_bayesr.cpp"="259b28b3d649dbc3698e0152f67d5fe8","src/blr_bed_bayesrc_types.h"="8d9201efad3832a562da1f14bc2e9d06","src/blr_bed_bayesrc_core_impl.h"="82365cf3f1f5306c57b980f59b4d83d3","src/blr_bed_bayesrc_aggregate_impl.h"="9ef26e7e1d587149f8801f7fedfed162","src/stblr_cpg_omp_bed_marker_scheduled_chains_bayesrc.cpp"="9ef7d514895f80b8561de831798f2701","R/RcppExports.R"="9d13ea00b326c7e0cd606194d13a8bca","src/RcppExports.cpp"="b4859db0f6308fa7e38051ddcf32d245","NAMESPACE"="f5b6ee37a3972aa436357bdc8f602f4e")
 expect_identical(unname(tools::md5sum(file.path(phase15a_root,names(protected)))),unname(protected))
})

test_that("public schema conventions remain common without merging converters",{
 adapters<-vapply(c("src/st_cpg_omp_individual_scheduled_chains.cpp","src/stblr_cpg_omp_bed_marker_scheduled_chains_bayesr.cpp","src/stblr_cpg_omp_bed_marker_scheduled_chains_bayesrc.cpp"),phase15a_text,character(1))
 for(x in adapters){expect_match(x,"stblr_raw_v1",fixed=TRUE);expect_match(x,"R_NilValue",fixed=TRUE)}
 expect_true(length(unique(grep("result_to_raw",adapters,value=TRUE)))==3L)
})
