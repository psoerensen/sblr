phase10b_root <- normalizePath(file.path(testthat::test_path(),"..",".."),winslash="/",mustWork=TRUE)
source(file.path(testthat::test_path(),"fixtures","blr-phase10b-scheduled-reference.R"))

phase10b_comparable <- function(x) {
  if (is.list(x) && !is.null(x$input)) x$input$ncores <- 0L
  if (is.list(x) && !is.null(x$meta)) x$meta$ncores <- 0L
  if (is.list(x)) for(nm in intersect(names(x),c("seconds_mean","seconds_max"))) x[[nm]][] <- 0
  x
}

test_that("Phase 10B chain RNG reconstructs exactly and isolates odd other-chain draws", {
  out <- sblr:::blr_phase10b_chain_rng_diagnostic_cpp(717L,9L)
  expect_identical(out$first,out$identical_seed)
  expect_identical(out$first,out$after_odd_other_chain)
  expect_identical(out$owner,"chain")
  expect_identical(out$lifetime,"one_chain_execution")
  expect_identical(out$worker_thread_owner,"none")
})

test_that("scheduled production owns distributions by logical chain", {
  source <- paste(readLines(file.path(phase10b_root,"src","st_cpg_omp_csr_scheduled.cpp"),warn=FALSE),collapse="\n")
  core <- paste(readLines(file.path(phase10b_root,"src","blr_csr_scheduled_bayesc_core_impl.h"),warn=FALSE),collapse="\n")
  native <- paste(source,core,sep="\n")
  types <- paste(readLines(file.path(phase10b_root,"src","blr_scheduled_execution_types.h"),warn=FALSE),collapse="\n")
  expect_false(grepl("static thread_local",native,fixed=TRUE))
  expect_false(grepl("thread_local",native,fixed=TRUE))
  expect_match(core,"ScheduledChainRng chain_rng(task_seed)",fixed=TRUE)
  expect_match(source,"rng.normal(rng.engine)",fixed=TRUE)
  expect_match(source,"rng.uniform(rng.engine)",fixed=TRUE)
  expect_match(types,"std::mt19937 engine",fixed=TRUE)
  expect_match(types,"std::normal_distribution<double> normal",fixed=TRUE)
  expect_match(types,"std::uniform_real_distribution<double> uniform",fixed=TRUE)
  expect_false(grepl("old_rng|new_rng|rng_selector|fallback",native))
})

test_that("post-correction raw and formatted references are exact", {
  for(nm in names(phase10b_configs)) {
    ref <- readRDS(file.path(testthat::test_path(),"fixtures","blr_phase10b_scheduled_csr",paste0(nm,".rds")))
    observed <- phase10b_run(phase10b_configs[[nm]])
    expect_identical(observed$raw,ref$raw,info=paste(nm,"raw"))
    expect_identical(observed$fit,ref$fit,info=paste(nm,"fit"))
    expect_identical(ref$metadata$rng_ownership_version,"scheduled_chain_rng_v1")
  }
})

test_that("same-process fits are independent of intervening scheduled and unscheduled work", {
  a <- phase10b_configs$skip_two_one
  b <- a; b$seeds <- c(3101L,3102L); b$full <- 2L; b$base <- 3L
  first <- phase10b_run(a)$fit
  expect_identical(phase10b_run(a)$fit,first)
  invisible(phase10b_run(b))
  expect_identical(phase10b_run(a)$fit,first)
  dense <- phase10b_configs$dense_one
  one <- phase10b_run(dense)$fit
  invisible(sblr::stblr_csr(stats=phase10a_stats(),ld_prefix=phase10a_prefix(),
    scheduled=FALSE,pi_init=.35,pi_prior_mean=.35,pi_prior_strength=3,
    updateB=FALSE,updateE=FALSE,updatePi=FALSE,nit=5L,nburn=1L,seed=444L))
  expect_identical(phase10b_run(dense)$fit,one)
  invisible(phase10b_run(a))
  expect_identical(phase10b_run(dense)$fit,one)
})

test_that("explicit seeds are independent of 1,2,2,1 core assignment", {
  cfg <- phase10b_configs$skip_two_one
  run_core <- function(k) { cfg$ncores <- k; phase10b_comparable(phase10b_run(cfg)$fit) }
  one <- run_core(1L); two <- run_core(2L); two_again <- run_core(2L); one_again <- run_core(1L)
  expect_identical(two,one)
  expect_identical(two_again,one)
  expect_identical(one_again,one)
})

test_that("fresh-process results equal reused-process results", {
  if (!identical(Sys.getenv("SBLR_RUN_FRESH_PROCESS_TESTS"),"true")) {
    for(nm in names(phase10b_configs)) {
      ref <- readRDS(file.path(testthat::test_path(),"fixtures","blr_phase10b_scheduled_csr",paste0(nm,".rds")))
      expect_identical(ref$metadata$reference_mode,"fresh R process")
    }
    return(invisible())
  }
  for(nm in names(phase10b_configs)) {
    fresh <- callr::r(function(root,name) {
      setwd(root); pkgload::load_all(".",compile=FALSE,quiet=TRUE)
      source(file.path("tests","testthat","fixtures","blr-phase10b-scheduled-reference.R"))
      phase10b_run(phase10b_configs[[name]])
    },list(root=phase10b_root,name=nm))
    reused <- phase10b_run(phase10b_configs[[nm]])
    expect_identical(fresh,reused,info=nm)
  }
})

test_that("scheduler and public limitations remain unchanged", {
  expect_error(sblr::stblr_csr(stats=phase10a_stats(),ld_prefix=phase10a_prefix(),
    scheduled=TRUE,keep_chains=TRUE),"keep_chains")
  expect_error(sblr::stblr_csr(stats=list(),method="bayesr",scheduled=TRUE),
               "scheduled CSR BayesR is not currently implemented")
  scheduled <- phase10b_run(phase10b_configs$dense_one)$fit
  ordinary <- sblr::stblr_csr(stats=phase10a_stats(),ld_prefix=phase10a_prefix(),
    scheduled=FALSE,pi_init=.35,pi_prior_mean=.35,pi_prior_strength=3,
    updateB=FALSE,updateE=FALSE,updatePi=FALSE,nit=8L,nburn=2L,nthin=1L,
    seed=1001L,nchains=1L,ncores=1L,updateLDswap=FALSE)
  expect_false(identical(phase10b_normalize(scheduled),phase10b_normalize(ordinary)))
  expect_identical(scheduled$input$backend,"csr_scheduled_bayesc")
  expect_identical(scheduled$input$scheduled,TRUE)
})

test_that("Phase 10B protects unrelated production sources and NAMESPACE", {
  protected <- c(
    "src/st_cpg_omp_csr.cpp"="92dafc0266d5a0e72aea000224154cef",
    "src/st_cpg_omp_csr_bayesr.cpp"="0a005f9d5a19037285fd4869fdc4dcf0",
    "src/st_sbayesrc_omp_csr.cpp"="8c1b03d8f5b93e6831ccbed856c77ead",
    "src/st_cpg_omp_csr_prior.cpp"="cce51072da6ddc3c18d58ab3b1f3c6df",
    "src/st_cpg_omp_csr_group.cpp"="87e923f7f8ee6420e39d9f041263d11b",
    "src/st_cpg_omp_csr_annot.cpp"="59bd49f048d116d0fe61d73d79bd4693",
    "src/st_cpg_omp_individual_scheduled.cpp"="0d726fe3faf5deec887381c1458ab6b6",
    "src/st_cpg_omp_individual_scheduled_chains.cpp"="43e8ff78b759656c5e897da186c1a548",
    "src/st_block_eigen.cpp"="49f0a62c9fe235967a264b0f8de144a7",
    "src/mt_cpg_omp_csr.cpp"="aec85896b5c30db3014efaeb5e3c3a96",
    "NAMESPACE"="f5b6ee37a3972aa436357bdc8f602f4e")
  expect_identical(unname(tools::md5sum(file.path(phase10b_root,names(protected)))),unname(protected))
})
