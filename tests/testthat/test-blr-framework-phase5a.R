phase5a_path <- function(...) {
  p <- file.path(...); if (file.exists(p)) return(p)
  file.path("..","..",...)
}
source(phase5a_path("tests","testthat","fixtures","blr-phase5a-bayesr-reference.R"))

phase5a_contract <- function() list(
  data=list(marker_count=4,trait_count=2,sample_size=c(100L,100L),
    marker_order=paste0("m",1:4),trait_order=c("T1","T2"),
    shared_read_only=TRUE,per_chain_payload=FALSE,storage_outlives_execution=TRUE),
  component=list(scales=c(0,.01,.1,1),initial_probability=c(.7,.1,.1,.1),
    dirichlet_prior=c(14,2,2,2),null_component=0,update_probability=TRUE,
    scale_interpretation="variance_multiplier"),
  controls=list(iterations=14L,burnin=4L,thinning=1L,chains=2L,cores=2L,seed=31L,
    chain_seeds=c(401L,402L),keep_chains=TRUE,update_marker_variance=TRUE,
    update_residual_variance=FALSE,update_ld_swap=FALSE,ld_swap_probability=.05,
    ld_swap_r2=.8,ld_swap_max_friends=50L,ld_swap_moves=1L),
  output=list(keep_chains=TRUE))

test_that("Phase 5A BayesR typed contract round-trips exactly", {
  x <- phase5a_contract(); y <- blr_phase5a_validate_bayesr_contract_cpp(x)
  expect_true(y$validated); expect_false(y$invokes_sampler)
  expect_identical(y$data,x$data); expect_identical(y$component,x$component)
  expect_identical(y$controls,x$controls); expect_identical(y$output,x$output)
  expect_identical(y$data$marker_order,paste0("m",1:4))
  expect_identical(y$data$trait_order,c("T1","T2"))
  expect_identical(y$component$scales,c(0,.01,.1,1))
  expect_identical(y$component$null_component,0)
})

test_that("Phase 5A BayesR contract rejects invalid dimensions and ownership", {
  mutate <- function(path,value) { x<-phase5a_contract(); x[[path[1]]][[path[2]]]<-value; x }
  expect_error(blr_phase5a_validate_bayesr_contract_cpp(mutate(c("data","marker_count"),0)),"marker_count")
  expect_error(blr_phase5a_validate_bayesr_contract_cpp(mutate(c("data","trait_count"),0)),"trait_count")
  expect_error(blr_phase5a_validate_bayesr_contract_cpp(mutate(c("data","shared_read_only"),FALSE)),"read-only")
  expect_error(blr_phase5a_validate_bayesr_contract_cpp(mutate(c("data","per_chain_payload"),TRUE)),"per-chain")
  expect_error(blr_phase5a_validate_bayesr_contract_cpp(mutate(c("data","storage_outlives_execution"),FALSE)),"outlive")
  expect_error(blr_phase5a_validate_bayesr_contract_cpp(mutate(c("component","scales"),0)),"at least two")
  expect_error(blr_phase5a_validate_bayesr_contract_cpp(mutate(c("component","null_component"),9)),"out of range")
  expect_error(blr_phase5a_validate_bayesr_contract_cpp(mutate(c("component","scales"),c(0,-1,.1,1))),"positive")
  expect_error(blr_phase5a_validate_bayesr_contract_cpp(mutate(c("component","initial_probability"),c(.5,.5))),"dimensions")
  expect_error(blr_phase5a_validate_bayesr_contract_cpp(mutate(c("component","initial_probability"),c(.7,-.1,.2,.2))),"non-negative")
  expect_error(blr_phase5a_validate_bayesr_contract_cpp(mutate(c("controls","iterations"),0L)),"iterations")
  expect_error(blr_phase5a_validate_bayesr_contract_cpp(mutate(c("controls","chains"),0L)),"chains")
  expect_error(blr_phase5a_validate_bayesr_contract_cpp(mutate(c("controls","cores"),0L)),"cores")
  expect_error(blr_phase5a_validate_bayesr_contract_cpp(mutate(c("controls","chain_seeds"),401L)),"chain_seeds")
})

test_that("Phase 5A typed header is binding neutral and encodes result vocabulary", {
  z <- paste(readLines(phase5a_path("src","blr_csr_bayesr_types.h"),warn=FALSE),collapse="\n")
  banned <- c("Rcpp","RcppArmadillo","SEXP","RObject","NumericVector","NumericMatrix",
    "Nullable","Rcpp::stop","Rcpp::Rcout","R::rnorm","R::rchisq","R::pnorm","R::qnorm",
    "arma::randn","arma::randu","pybind11","Python.h")
  for (token in banned) expect_false(grepl(token,z,fixed=TRUE),info=token)
  expect_match(z,"struct CsrBayesRDataView",fixed=TRUE)
  expect_match(z,"struct CsrBayesRExecutionInput",fixed=TRUE)
  expect_match(z,"struct CsrBayesRResult",fixed=TRUE)
  expect_match(z,"per_chain_payload",fixed=TRUE)
})

test_that("Phase 5A frozen production BayesR raw and fit references remain exact", {
  for (nm in names(phase5a_bayesr_configs)) {
    ref <- readRDS(phase5a_path("tests","testthat","fixtures","blr_phase5a_bayesr",paste0(nm,".rds")))
    cfg <- phase5a_bayesr_configs[[nm]]
    expect_identical(phase5a_bayesr_normalize(phase5a_bayesr_run(cfg,TRUE)),ref$raw,info=paste(nm,"raw"))
    expect_identical(phase5a_bayesr_normalize(phase5a_bayesr_run(cfg,FALSE)),ref$fit,info=paste(nm,"fit"))
    expect_identical(ref$metadata$schema_name,"stblr_raw")
    expect_identical(ref$metadata$schema_version,1L)
  }
})

test_that("Phase 5A protects BayesC, block-eigen, and namespace sources", {
  paths <- c("R/sparse_ld_bed_helper.R","docs/dev/stblr_raw_schema.md",
    "src/st_cpg_omp_csr.cpp","src/blr_csr_bayesc_types.h",
    "src/blr_csr_bayesc_core_impl.h","src/st_block_eigen.cpp","src/st_block_eigen.h","NAMESPACE")
  expected <- c("26c5c894058434deea25e3242dd56d4a","82ac9ba4b7d8edc6f3e16ee3a26d8466",
    "92dafc0266d5a0e72aea000224154cef","e5975c311c69fe536db57dd21f01334f",
    "f7c617cbfc172639c1f8aea1bd8b1876","49f0a62c9fe235967a264b0f8de144a7",
    "bec3bc1e41841ab77747e34dc9818574","f5b6ee37a3972aa436357bdc8f602f4e")
  expect_identical(unname(tools::md5sum(vapply(paths,phase5a_path,character(1)))),expected)
  src <- paste(readLines(phase5a_path("src","st_cpg_omp_csr_bayesr.cpp"),warn=FALSE),collapse="\n")
  core <- paste(readLines(phase5a_path("src","blr_csr_bayesr_core_impl.h"),warn=FALSE),collapse="\n")
  report <- paste(readLines(phase5a_path("docs","dev","blr_framework_phase5a_report.md"),warn=FALSE),collapse="\n")
  expect_match(report,"Approved Phase 5B seam",fixed=TRUE)
  expect_false(grepl("blr_phase5a_validate_bayesr_contract_cpp",src,fixed=TRUE))
  expect_equal(length(gregexpr("for (int it = 0; it < trace_len",paste(src,core),fixed=TRUE)[[1L]]),1L)
  expect_match(src,"run_bayesr_execution(execution_context)",fixed=TRUE)
  expect_match(src,"csr_bayesr_result_to_raw",fixed=TRUE)
  expect_false(grepl("getenv",paste(src,core),fixed=TRUE))
})
