source(file.path(testthat::test_path(),"fixtures","blr-phase13a-bed-bayesr-reference.R"))
phase13e_text <- function(path) paste(readLines(blr_repo_path(path),warn=FALSE),collapse="\n")
phase13e_pick <- function(x, current, historical = current) {
  if (!is.null(x[[current]])) x[[current]] else x[[historical]]
}
phase13e_raw_science <- function(x) x[c("marker", "trace", "variance", "pi", "component")]
phase13e_fit_science <- function(x) list(
  bm=x$bm, dm=x$dm, wy=x$wy, r=x$r, b=x$b, d=x$d,
  vbs=x$vbs, vgs=x$vgs, ves=x$ves, vle=x$vle, vld=x$vld,
  component=x$component, dm_component_mean=x$dm_component_mean,
  mixture_var=x$mixture_var,
  component_probabilities=phase13e_pick(x,"component_probabilities","comp_prob"),
  pi_final=as.numeric(phase13e_pick(x,"pi_final","final_pi")),
  pi_mean=as.numeric(phase13e_pick(x,"pi_mean","mean_pi")),
  cov_b_mean=phase13e_pick(x,"cov_b_mean","covb"),
  cov_g_mean=phase13e_pick(x,"cov_g_mean","covg"),
  cov_e_mean=phase13e_pick(x,"cov_e_mean","cove"),
  cov_b_final=phase13e_pick(x,"cov_b_final","vb"),
  cov_g_final=phase13e_pick(x,"cov_g_final","vg"),
  cov_e_final=phase13e_pick(x,"cov_e_final","ve")
)

test_that("Phase 13E permanently protects the canonical BayesR architecture", {
  types <- phase13e_text("src/blr_bed_bayesr_types.h")
  family_types <- phase13e_text("src/blr_bed_family_types.h")
  core <- phase13e_text("src/blr_bed_bayesr_core_impl.h")
  aggregate <- phase13e_text("src/blr_bed_bayesr_aggregate_impl.h")
  adapter <- phase13e_text("src/stblr_cpg_omp_bed_marker_scheduled_chains_bayesr.cpp")
  for (symbol in c("BedBayesRComponentSpec",
      "BedBayesRChainExecutionContext","BedBayesRChainExecutionResult",
      "BedBayesRExecutionResult")) expect_source_count(paste0("struct ",symbol),types,1L)
  expect_source_count("using BedBayesRPackedGenotypeView",types,1L)
  expect_source_count("struct BedPackedGenotypeView",family_types,1L)
  expect_source_count("BedBayesRChainExecutionResult run_bed_bayesr_chain(",core,1L)
  expect_source_count("for (int it = 0; it < total_it; ++it)",core,1L)
  expect_source_count("BedBayesRExecutionResult aggregate_bed_bayesr_results(",aggregate,1L)
  expect_source_count("static Rcpp::List stblr_bed_bayesr_result_to_raw(",adapter,2L)
  expect_source_count("#pragma omp parallel for num_threads(nthreads) schedule(static)",adapter,1L)
  expect_source_count("for (const sblr::core::BedBayesRProgressEvent& event",adapter,1L)
  expect_false(grepl("Rcpp|SEXP|Python.h|pybind11|fopen|fseek|fread",paste(core,aggregate,types)))
  expect_false(grepl("thread_local|std::discrete_distribution|old_path|new_path|aggregation_selector|rng_selector",paste(core,aggregate,adapter)))
  expect_match(family_types,"const PackedGenotype& storage",fixed=TRUE)
  expect_match(adapter,"br_read_bed_blocked",fixed=TRUE)
})

test_that("Phase 13E preserves scheduler, sampling and aggregation semantics", {
  core <- phase13e_text("src/blr_bed_bayesr_core_impl.h")
  aggregate <- phase13e_text("src/blr_bed_bayesr_aggregate_impl.h")
  for (needle in c("std::mt19937 gen_t(chain_seed)","std::uniform_int_distribution<int> jitter_dist",
      "if (u < cumsum)","for (int marker : active_list)","for (int marker : candidate_list)",
      "const std::vector<int>& due","progress_events.push_back")) expect_match(core,needle,fixed=TRUE)
  for (needle in c("ch*nt+t","out.bm*=inv","out.dm*=inv","out.final_pi*=inv",
      "out.mean_pi*=inv","nchains-1","arma::sqrt")) expect_match(aggregate,needle,fixed=TRUE)
})

test_that("Phase 13E canonical raw and formatted fixtures are exact", {
  configs <- list(one_chain_one_core=c(1L,1L,71L),two_chains_one_core=c(1L,2L,73L),
    two_chains_two_cores=c(2L,2L,73L))
  for (name in names(configs)) {
    z <- configs[[name]]; ref <- readRDS(file.path(testthat::test_path(),"fixtures",
      "blr_phase13a_bed_bayesr",paste0(name,".rds")))
    observed <- phase13a_capture(z[1],z[2],z[3])
    expect_equal(phase13e_raw_science(phase13a_normalize(observed$raw)),
      phase13e_raw_science(phase13a_normalize(ref$raw)),tolerance=1e-12)
    expect_equal(phase13e_fit_science(phase13a_normalize(observed$fit)),
      phase13e_fit_science(phase13a_normalize(ref$fit)),tolerance=1e-12)
  }
})

test_that("Phase 13E canonical reproducibility and identities remain exact", {
  a <- phase13a_normalize(phase13a_capture(1L,2L,73L))
  expect_identical(phase13a_normalize(phase13a_capture(1L,2L,73L)),a)
  phase13a_capture(1L,1L,79L)
  expect_identical(phase13a_normalize(phase13a_capture(1L,2L,73L)),a)
  expect_identical(phase13a_normalize(phase13a_capture(2L,2L,73L)),a)
  x <- phase13a_capture(1L,2L,73L)$fit
  expect_equal(sum(x$pi_final),1,tolerance=1e-12)
  expect_equal(sum(x$pi_mean),1,tolerance=1e-12)
  expect_equal(x$dm[,1L],1-x$component_probabilities[[1L]][,1L],tolerance=1e-12)
  expect_true(all(x$component[,1L]>=0 & x$component[,1L]<=3))
})

test_that("Phase 13E documentation identifies BayesR as canonical", {
  plan <- phase13e_text("docs/dev/blr_framework_implementation_plan.md")
  matrix <- phase13e_text("docs/dev/blr_model_capability_matrix.md")
  expect_match(plan,"Phase 13E designates the public packed-BED BayesR route canonical",fixed=TRUE)
  expect_match(matrix,"Canonical typed component/genotype/context/core",fixed=TRUE)
  expect_match(plan,"fixtures and the Phase 13D runtime/completed-fit-RSS/I/O measurements",fixed=TRUE)
})
