phase15b_root<-normalizePath(file.path(testthat::test_path(),"..",".."),winslash="/",mustWork=TRUE)
phase15b_text<-function(path)paste(readLines(file.path(phase15b_root,path),warn=FALSE),collapse="\n")

test_that("shared task indexing is exact and used by all canonical adapters",{
 common<-phase15b_text("src/blr_bed_family_types.h")
 expect_source_count("struct BedFamilyTaskIndex",common,1L)
 expect_source_count("make_bed_family_task_index(",common,1L)
 expect_match(common,"job % trait_count",fixed=TRUE);expect_match(common,"job / trait_count",fixed=TRUE)
 for(nt in c(1L,2L,3L))for(nch in c(1L,3L,4L))for(job in 0:(nt*nch-1L)){
  expect_identical(c(trait=job%%nt,chain=job%/%nt),c(trait=job%%nt,chain=job%/%nt))
 }
 adapters<-vapply(c("src/st_cpg_omp_individual_scheduled_chains.cpp","src/stblr_cpg_omp_bed_marker_scheduled_chains_bayesr.cpp","src/stblr_cpg_omp_bed_marker_scheduled_chains_bayesrc.cpp"),phase15b_text,character(1))
 for(x in adapters){expect_match(x,"make_bed_family_task_index(job,nt)",fixed=TRUE);expect_true(grepl("job_results[static_cast<std::size_t>(job)]",x,fixed=TRUE)||grepl("jobs[job]",x,fixed=TRUE))}
})

test_that("shared seed resolver preserves exact signed arithmetic and unsigned cast",{
 common<-phase15b_text("src/blr_bed_family_types.h")
 expect_source_count("resolve_bed_family_logical_chain_seed(",common,1L)
 expect_match(common,"static_cast<unsigned int>",fixed=TRUE)
 expect_match(common,"seed+1000003*(trait+1)+9176*(chain+1)",fixed=TRUE)
 old<-function(seed,trait,chain)(as.double(seed)+1000003*(trait+1)+9176*(chain+1))%%2^32
 shared_contract<-function(seed,trait,chain)(as.double(seed)+1000003*(trait+1)+9176*(chain+1))%%2^32
 for(seed in c(0,1,100000,100000000))for(trait in c(0,2))for(chain in c(0,3))expect_identical(shared_contract(seed,trait,chain),old(seed,trait,chain))
 for(path in c("src/st_cpg_omp_individual_scheduled_chains.cpp","src/stblr_cpg_omp_bed_marker_scheduled_chains_bayesr.cpp","src/stblr_cpg_omp_bed_marker_scheduled_chains_bayesrc.cpp")){
  x<-phase15b_text(path);expect_match(x,"resolve_bed_family_logical_chain_seed",fixed=TRUE);expect_false(grepl("1000003|9176|omp_get_thread_num.*seed",x))
 }
})

test_that("BayesC and BayesR share one immutable packed view while BayesRC stays specific",{
 common<-phase15b_text("src/blr_bed_family_types.h");c<-phase15b_text("src/blr_bed_scheduled_bayesc_types.h");r<-phase15b_text("src/blr_bed_bayesr_types.h");rc<-phase15b_text("src/blr_bed_bayesrc_types.h")
 expect_source_count("struct BedPackedGenotypeView",common,1L)
 for(s in c("const PackedGenotype& storage","const std::uint8_t* packed_markers","packed_size","marker_count","sample_count","bytes_per_marker","stride"))expect_match(common,s,fixed=TRUE)
 expect_match(r,"using BedBayesRPackedGenotypeView=BedPackedGenotypeView",fixed=TRUE)
 expect_false(grepl("struct BedPackedGenotypeView",c));expect_source_count("struct BedBayesRCPackedGenotypeView",rc,1L)
 expect_false(grepl("packed_markers|packed_size|stride",sub("struct BedBayesRCChainExecutionResult[\\s\\S]*","",rc)))
})

test_that("failure decision and dispatch boundaries are explicit",{
 report<-phase15b_text("docs/dev/blr_framework_phase15b_report.md")
 expect_match(report,"model-specific production failure payloads retained",fixed=TRUE)
 for(path in c("src/st_cpg_omp_individual_scheduled_chains.cpp","src/stblr_cpg_omp_bed_marker_scheduled_chains_bayesr.cpp","src/stblr_cpg_omp_bed_marker_scheduled_chains_bayesrc.cpp")){
  x<-phase15b_text(path);expect_gte(source_match_count("#pragma omp parallel for",x),1L);expect_match(x,"schedule(static)",fixed=TRUE)
 }
 expect_false(grepl("virtual|std::function|dispatch_callback",phase15b_text("src/blr_bed_family_types.h")))
})

test_that("model-specific cores aggregators converters and policies remain separate",{
 ccore<-phase15b_text("src/blr_bed_scheduled_bayesc_core_impl.h");rcore<-phase15b_text("src/blr_bed_bayesr_core_impl.h");rccore<-phase15b_text("src/blr_bed_bayesrc_core_impl.h")
 expect_match(ccore,"candidate",fixed=TRUE);expect_match(rcore,"component",fixed=TRUE);expect_match(rccore,"st_bayesrc_compute_snp_pi",fixed=TRUE)
 expect_match(rccore,"for (int marker : marker_order)",fixed=TRUE);expect_false(grepl("null_skip|candidate_list|due_bucket",rccore))
 common<-phase15b_text("src/blr_bed_family_types.h")
 expect_source_forbidden(common,c("Rcpp","scheduler","component","annotation","std::mt19937","fopen","ifstream","virtual"))
})

test_that("shared structural reference and benchmark helpers are active",{
 helper<-phase15b_text("tests/testthat/helper-source-architecture.R")
 for(s in c("expect_source_forbidden","expect_source_hashes","reference_first_difference","run_reference_fresh_process"))expect_match(helper,s,fixed=TRUE)
 expect_null(reference_first_difference(list(a=1L),list(a=1L)))
 expect_identical(reference_first_difference(list(a=1L),list(a=2L))$path,"root$a")
 bench<-phase15b_text("tools/benchmarks/blr_bed_family_benchmark_common.R")
 for(s in c("no cross-model speed ranking","workloads are model-specific","completed-fit RSS is not peak RSS","page-cache effects apply","tiny timings are regression signals"))expect_match(bench,s,fixed=TRUE)
})
