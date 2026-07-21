test_that("MT block-eigen contracts are binding-neutral and share the MT core", {
  types <- blr_source_text("src/blr_mt_block_eigen_types.h")
  access <- blr_source_text("src/blr_mt_ld_access.h")
  core <- blr_source_text("src/blr_mt_default_core_impl.h")
  native <- blr_source_text("src/mtblr.cpp")
  expect_match(types, "struct MtBlockEigenBundleView", fixed=TRUE)
  expect_match(types, "std::vector<sblr::core::BlockEigenView> trait_operator", fixed=TRUE)
  expect_match(types, "struct MtBlockEigenDataView", fixed=TRUE)
  expect_false(grepl("Rcpp|SEXP|shared_ptr", types))
  expect_match(access, "trait_operator[static_cast<std::size_t>(trait)]", fixed=TRUE)
  expect_match(core, "run_mt_block_eigen_core", fixed=TRUE)
  expect_equal(source_match_count("for ( int it =", core, fixed=TRUE), 1L)
  expect_equal(source_match_count("sampleBetaCPG_Mt_latent(i", core, fixed=TRUE), 1L)
  expect_match(native, "make_mt_default_legacy_result", fixed=TRUE)
})

test_that("shared block-eigen filters reduce exactly to dense MT", {
  phase17l_compare(phase17l_case(updates=FALSE))
  phase17l_compare(phase17l_case(updates=TRUE),1e-10)
  phase17l_compare(phase17l_case(filters=rep("ridge_fixed",2)),1e-12)
  phase17l_compare(phase17l_case(filters=rep("ridge_lw",2)),1e-12)
})

test_that("initial effects, sets, and three traits retain dense execution", {
  phase17l_compare(phase17l_case(nonzero=TRUE,updates=TRUE),1e-10)
  phase17l_compare(phase17l_case(multiple_sets=TRUE,updates=TRUE),1e-10)
  phase17l_compare(phase17l_case(nt=3L,filters=rep("hard_truncate",3),updates=TRUE),1e-10)
})

test_that("trait-specific boundaries and mixed filters use matching wy", {
  x <- phase17l_case(nt=3L,shared=FALSE,
    blocks=list(c(0L,2L),c(0L,1L,3L),c(0L,3L)),
    filters=c("hard_truncate","ridge_fixed","ridge_lw"),updates=TRUE)
  result <- phase17l_compare(x,1e-10)
  expect_false(identical(x$transformed[[1]],x$block$wy[[1]]))
  expect_identical(x$transformed[[2]],x$block$wy[[2]])
  expect_identical(x$transformed[[3]],x$block$wy[[3]])
  for (t in seq_len(3L)) {
    residual <- x$transformed[[t]]-as.numeric(x$matrices[[t]] %*% result[[5]][[t]])
    expect_equal(result[[4]][[t]],residual,tolerance=1e-10)
  }
})

test_that("one trait and differing analysis sample sizes are supported", {
  phase17l_compare(phase17l_case(nt=1L,filters="ridge_fixed"),1e-12)
  x <- phase17l_case(); x$dense$n <- x$block$n <- c(31L,73L)
  phase17l_compare(x,1e-12)
})

test_that("exactly representable block operators reduce to canonical CSR", {
  phase17l_compare_csr(phase17l_case(
    filters=rep("ridge_fixed",2),blocks=rep(list(0:3),2)),1e-12)
  phase17l_compare_csr(phase17l_case(
    filters=rep("ridge_fixed",2),blocks=rep(list(c(0L,2L)),2)),1e-10)
})

test_that("descriptor validation rejects malformed execution before MCMC", {
  x <- phase17l_case()
  x$block$operator_descriptors <- list(x$block$operator_descriptors[[1]],
    x$block$operator_descriptors[[1]],x$block$operator_descriptors[[1]])
  expect_error(do.call(sblr:::mtblr_block_eigen_internal,x$block),"length one or trait count")
  x <- phase17l_case(); x$block$operator_descriptors[[1]]$af <- c(.2,.3)
  expect_error(do.call(sblr:::mtblr_block_eigen_internal,x$block),"marker count")
  for (field in c("bed_files","n_bed","cls","af","block_start","eigen_filter","eigen_tau","eigen_eta")) {
    x <- phase17l_case(); x$block$operator_descriptors[[1]][field] <- NULL
    expect_error(do.call(sblr:::mtblr_block_eigen_internal,x$block),"lacks")
  }
  x <- phase17l_case(); x$block$operator_descriptors[[1]]$eigen_filter <- "bad"
  expect_error(do.call(sblr:::mtblr_block_eigen_internal,x$block),"eigen_filter")
  x <- phase17l_case(); x$block$method <- 5L
  expect_error(do.call(sblr:::mtblr_block_eigen_internal,x$block),"method = 4")
})

test_that("internal registration does not create a public block-eigen route", {
  expect_true(exists("mtblr_block_eigen_internal",asNamespace("sblr"),inherits=FALSE))
  expect_false("mtblr_block_eigen_internal" %in% getNamespaceExports("sblr"))
  expect_false("mtblr_block_eigen" %in% getNamespaceExports("sblr"))
  expect_false("block_start" %in% names(formals(sblr::mtblr_csr)))
  native <- blr_source_text("src/mtblr.cpp")
  expect_false(grepl("run_mt_block_eigen_core[\\s\\S]*mtblr_eigen\\(", native))
})
