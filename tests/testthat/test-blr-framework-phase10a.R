phase10a_root <- normalizePath(file.path(testthat::test_path(), "..", ".."), winslash="/", mustWork=TRUE)
source(file.path(testthat::test_path(),"fixtures","blr-phase10a-scheduled-reference.R"))

phase10a_spec <- function(m=6L) list(
  execution=list(marker_count=m,trait_count=1L,iterations=8L,burnin=2L,
    thinning=1L,chains=2L,cores=2L,seed=1001L,chain_seeds=c(11L,12L),
    keep_chains=FALSE),
  rng_ownership=list(engine_owner="chain",distribution_owner="chain",
    lifetime="one_chain_execution",worker_thread_owner="none",
    fit_persistent_distribution_state=FALSE),
  sweep=list(full_sweep_every=3L,iteration_zero_is_full=TRUE),
  skip=list(null_skip_base=2L,null_skip_max=7L,burnin_only=FALSE,
    growth_rule="probability_adaptive"),
  candidate=list(threshold=.1,lifetime=2L),
  neighbor=list(enabled=TRUE,difference_threshold=.02,maximum_neighbors=2L,
    friend_marker_count=m,shared_read_only=TRUE,storage_outlives_execution=TRUE),
  state=list(scheduled_at=rep(-1L,m),last_updated=rep(-1L,m),
    candidate=rep(FALSE,m),in_candidate_list=rep(FALSE,m),
    in_active_list=rep(FALSE,m),last_interesting=rep(-1000000000L,m))
)

test_that("Phase 10A scheduled contracts round trip without sampler or RNG", {
  spec <- phase10a_spec()
  out <- sblr:::blr_phase10a_validate_scheduled_execution_cpp(spec)
  expect_identical(out$schema,"blr_scheduled_execution_contract_v1")
  expect_identical(out$spec,spec)
  expect_true(out$validated); expect_false(out$invokes_sampler); expect_false(out$consumes_rng)
  bad <- spec; bad$sweep$full_sweep_every <- 0L
  expect_error(sblr:::blr_phase10a_validate_scheduled_execution_cpp(bad),"full_sweep_every")
  bad <- spec; bad$skip$null_skip_base <- 8L
  expect_error(sblr:::blr_phase10a_validate_scheduled_execution_cpp(bad),"null_skip_base")
  bad <- spec; bad$candidate$threshold <- 2
  expect_error(sblr:::blr_phase10a_validate_scheduled_execution_cpp(bad),"candidate_threshold")
  bad <- spec; bad$neighbor$shared_read_only <- FALSE
  expect_error(sblr:::blr_phase10a_validate_scheduled_execution_cpp(bad),"borrowed immutable")
  bad <- spec; bad$execution$chain_seeds <- 1L
  expect_error(sblr:::blr_phase10a_validate_scheduled_execution_cpp(bad),"chain_seeds")
  bad <- spec; bad$state$scheduled_at <- integer(5)
  expect_error(sblr:::blr_phase10a_validate_scheduled_execution_cpp(bad),"state dimensions")
})

test_that("persistent normal distribution retains cached state after engine reseed", {
  out <- sblr:::blr_phase10a_distribution_cache_diagnostic_cpp(913L,2L)
  expect_true(out$cached_state_survives_engine_reseed)
  expect_identical(out$first,out$after_distribution_reset)
  expect_false(identical(out$after_engine_reseed_without_distribution_reset,
                         out$after_distribution_reset))
})

test_that("Phase 10A audit artifact and corrected production ownership are explicit", {
  scheduled <- readLines(file.path(phase10a_root,"src","st_cpg_omp_csr_scheduled.cpp"),warn=FALSE)
  scheduled_core <- readLines(file.path(phase10a_root,"src","blr_csr_scheduled_bayesc_core_impl.h"),warn=FALSE)
  scheduled_native <- paste(c(scheduled, scheduled_core),collapse="\n")
  route <- readLines(file.path(phase10a_root,"R","sparse_ld_bed_helper.R"),warn=FALSE)
  expect_false(grepl("static thread_local",scheduled_native,fixed=TRUE))
  expect_match(scheduled_native,"ScheduledChainRng chain_rng",fixed=TRUE)
  expect_match(scheduled_native,"std::chi_squared_distribution<double>",fixed=TRUE)
  expect_match(scheduled_native,"std::gamma_distribution<double>",fixed=TRUE)
  expect_match(paste(route,collapse="\n"),"scheduled CSR BayesR is not currently implemented",fixed=TRUE)
  expect_match(paste(route,collapse="\n"),"stblr_cpg_omp_csr_scheduled",fixed=TRUE)
})

test_that("Phase 10A defective fresh-process references remain immutable audit artifacts", {
  for (nm in names(phase10a_configs)) {
    reference <- readRDS(file.path(testthat::test_path(),"fixtures","blr_phase10a_scheduled",paste0(nm,".rds")))
    expect_identical(reference$metadata$reference_mode,"fresh R process")
    expect_identical(reference$metadata$starting_commit,phase10a_starting_commit)
    expect_identical(reference$raw$schema$class,"stblr_raw")
    expect_identical(reference$raw$schema$version,1L)
  }
})

test_that("same-process call order defect documented by Phase 10A is corrected", {
  a <- phase10a_configs$skip_two_one
  b <- a; b$seeds <- c(2101L,2102L); b$full <- 2L; b$base <- 3L
  first <- phase10a_run(a)$fit
  repeated <- phase10a_run(a)$fit
  expect_identical(repeated,first)
  invisible(phase10a_run(b))
  after_b <- phase10a_run(a)$fit
  expect_identical(after_b,first)
  expect_identical(after_b$bm,first$bm)
})

test_that("core order is independent of worker-thread assignment", {
  cfg <- phase10a_configs$skip_two_one
  run_core <- function(k) { cfg$ncores <- k; phase10a_run(cfg)$fit }
  one <- run_core(1L); two <- run_core(2L); two_again <- run_core(2L); one_again <- run_core(1L)
  comparable <- function(x) { x$input$ncores <- 0L; x }
  expect_identical(comparable(one),comparable(two))
  expect_identical(comparable(two),comparable(two_again))
  expect_identical(comparable(one),comparable(one_again))
})

test_that("dense scheduled controls are audited rather than assumed canonical", {
  z <- phase10a_stats(); prefix <- phase10a_prefix()
  common <- list(stats=z,ld_prefix=prefix,pi_init=.35,pi_prior_mean=.35,
    pi_prior_strength=3,updateB=FALSE,updateE=FALSE,updatePi=FALSE,
    nit=8L,nburn=2L,nthin=1L,seed=333L,nchains=1L,ncores=1L,
    updateLDswap=FALSE)
  scheduled <- do.call(sblr::stblr_csr,c(common,list(scheduled=TRUE,
    full_sweep_every=1L,null_skip_base=1L,null_skip_max=1L,
    candidate_threshold=0,candidate_lifetime=0L,wakeup_ld_neighbors=FALSE)))
  ordinary <- do.call(sblr::stblr_csr,c(common,list(scheduled=FALSE)))
  expect_false(identical(phase10a_normalize(scheduled),phase10a_normalize(ordinary)))
})

test_that("Phase 10A protects the extracted scheduled route and unrelated sources", {
  scheduled <- paste(readLines(file.path(phase10a_root,"src","st_cpg_omp_csr_scheduled.cpp"),
    warn=FALSE),collapse="\n")
  expect_match(scheduled,"#include \"blr_csr_scheduled_bayesc_core_impl.h\"",fixed=TRUE)
  protected <- c(
    "src/st_cpg_omp_csr.cpp"="92dafc0266d5a0e72aea000224154cef",
    "src/st_cpg_omp_csr_bayesr.cpp"="0a005f9d5a19037285fd4869fdc4dcf0",
    "src/st_sbayesrc_omp_csr.cpp"="8c1b03d8f5b93e6831ccbed856c77ead",
    "src/st_cpg_omp_csr_prior.cpp"="cce51072da6ddc3c18d58ab3b1f3c6df",
    "src/st_cpg_omp_csr_group.cpp"="87e923f7f8ee6420e39d9f041263d11b",
    "src/st_cpg_omp_csr_annot.cpp"="59bd49f048d116d0fe61d73d79bd4693",
    "src/st_cpg_omp_individual_scheduled.cpp"="0d726fe3faf5deec887381c1458ab6b6",
    "src/st_cpg_omp_individual_scheduled_chains.cpp"="f58fbefcffb183b9d54a96b398321dfb",
    "NAMESPACE"="f5b6ee37a3972aa436357bdc8f602f4e")
  actual <- unname(tools::md5sum(file.path(phase10a_root,names(protected))))
  expect_identical(actual,unname(protected))
})
