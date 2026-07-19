phase13d_root <- normalizePath(file.path(testthat::test_path(), "..", ".."), winslash="/", mustWork=TRUE)
source(file.path(testthat::test_path(), "fixtures", "blr-phase13a-bed-bayesr-reference.R"))
phase13d_text <- function(path) paste(readLines(file.path(phase13d_root,path),warn=FALSE),collapse="\n")

test_that("Phase 13D has one typed aggregation and one named converter", {
  types <- phase13d_text("src/blr_bed_bayesr_types.h")
  family_types <- phase13d_text("src/blr_bed_family_types.h")
  core <- phase13d_text("src/blr_bed_bayesr_core_impl.h")
  aggregate <- phase13d_text("src/blr_bed_bayesr_aggregate_impl.h")
  adapter <- phase13d_text("src/stblr_cpg_omp_bed_marker_scheduled_chains_bayesr.cpp")
  for (symbol in c("BedBayesRChainExecutionContext", "BedBayesRChainExecutionResult",
      "BedBayesRAggregationContext", "BedBayesRExecutionResult"))
    expect_source_count(paste0("struct ",symbol),types,1L)
  expect_source_count("BedBayesRChainExecutionResult run_bed_bayesr_chain(",core,1L)
  expect_source_count("BedBayesRExecutionResult aggregate_bed_bayesr_results(",aggregate,1L)
  expect_source_count("static Rcpp::List stblr_bed_bayesr_result_to_raw(",adapter,2L)
  expect_source_count("aggregate_bed_bayesr_results(job_results,aggregation_context)",adapter,1L)
  expect_source_count("stblr_bed_bayesr_result_to_raw(result,binding_metadata)",adapter,1L)
  expect_source_count("#pragma omp parallel for num_threads(nthreads) schedule(static)",adapter,1L)
  expect_source_count("for (int it = 0; it < total_it; ++it)",core,1L)
  expect_false(grepl("Rcpp|SEXP|Python.h|pybind11|fopen|fseek|fread",aggregate))
  expect_false(grepl("old_path|new_path|aggregation_selector|rng_selector|fallback",paste(core,aggregate,adapter)))
})

test_that("Phase 13D aggregation permanently protects legacy formulas", {
  aggregate <- phase13d_text("src/blr_bed_bayesr_aggregate_impl.h")
  for (needle in c("ch*nt+t", "out.bm*=inv", "out.dm*=inv",
      "out.component_probability", "out.final_pi*=inv", "out.mean_pi*=inv",
      "out.log_cpo*=inv", "out.seconds_max", "nchains-1",
      "arma::sqrt", "std::min", "std::max"))
    expect_match(aggregate,needle,fixed=TRUE)
})

test_that("Phase 13D preserves core, progress, component and genotype boundaries", {
  core <- phase13d_text("src/blr_bed_bayesr_core_impl.h")
  types <- phase13d_text("src/blr_bed_bayesr_types.h")
  family_types <- phase13d_text("src/blr_bed_family_types.h")
  adapter <- phase13d_text("src/stblr_cpg_omp_bed_marker_scheduled_chains_bayesr.cpp")
  for (needle in c("std::mt19937 gen_t(chain_seed)","if (u < cumsum)",
      "std::uniform_int_distribution<int> jitter_dist","for (int marker : active_list)",
      "for (int marker : candidate_list)","const std::vector<int>& due"))
    expect_match(core,needle,fixed=TRUE)
  expect_false(grepl("Rcpp|SEXP|std::discrete_distribution|thread_local",core))
  expect_match(family_types,"const PackedGenotype& storage",fixed=TRUE)
  expect_match(adapter,"for (const sblr::core::BedBayesRProgressEvent& event",fixed=TRUE)
  expect_match(adapter,"br_read_bed_blocked",fixed=TRUE)
})

test_that("Phase 13D raw and formatted references remain exact", {
  configs <- list(one_chain_one_core=c(1L,1L,71L),two_chains_one_core=c(1L,2L,73L),
    two_chains_two_cores=c(2L,2L,73L))
  for (name in names(configs)) {
    z <- configs[[name]]
    ref <- readRDS(file.path(testthat::test_path(),"fixtures","blr_phase13a_bed_bayesr",paste0(name,".rds")))
    observed <- phase13a_capture(z[1],z[2],z[3])
    expect_identical(phase13a_normalize(observed$raw),phase13a_normalize(ref$raw))
    expect_identical(phase13a_normalize(observed$fit),phase13a_normalize(ref$fit))
  }
})

test_that("Phase 13D identities and call-order reproducibility remain exact", {
  a <- phase13a_normalize(phase13a_capture(1L,2L,73L))
  expect_identical(phase13a_normalize(phase13a_capture(1L,2L,73L)),a)
  phase13a_capture(1L,1L,79L)
  expect_identical(phase13a_normalize(phase13a_capture(1L,2L,73L)),a)
  expect_identical(phase13a_normalize(phase13a_capture(2L,2L,73L)),a)
  x <- phase13a_capture(1L,2L,73L)$fit
  expect_equal(unname(rowSums(x$pi)),rep(1,nrow(x$pi)),tolerance=1e-12)
  expect_equal(unname(rowSums(x$pim)),rep(1,nrow(x$pim)),tolerance=1e-12)
  expect_equal(x$dm[,1L],1-x$comp_prob[[1L]][,1L],tolerance=1e-12)
})
