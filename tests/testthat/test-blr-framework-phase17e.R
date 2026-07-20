test_that("Phase 17E owns typed binding-neutral core contracts", {
  types <- blr_source_text("src/blr_mt_default_types.h")
  core <- blr_source_text("src/blr_mt_default_core_impl.h")
  for (name in c("MtDefaultDataView", "MtDefaultModelSpec",
      "MtDefaultCovariancePriorView", "MtDefaultExecutionSpec",
      "MtDefaultInitialState", "MtDefaultCoreResult"))
    expect_source_count(paste("struct", name), types, 1L)
  expect_source_forbidden(paste(types, core), c("Rcpp", "SEXP", "RObject",
    "NumericVector", "NumericMatrix", "R::"))
  expect_source_count("inline MtDefaultCoreResult run_mt_default_core(", core, 1L)
})

test_that("Phase 17E owns core execution and RNG boundaries", {
  core <- blr_source_text("src/blr_mt_default_core_impl.h")
  adapter <- blr_mt_public_source()
  expect_source_count("std::mt19937 gen(execution.seed);", core, 1L)
  expect_source_count("for ( int it = 0; it < execution.nit+execution.nburn; it++)",
    core, 1L)
  expect_source_forbidden(core, c("omp_get_thread_num", "static std::mt19937",
    "thread_local"))
  expect_source_count("std::mt19937", adapter, 0L)
  expect_source_count("run_mt_default_core(", adapter, 1L)
})

test_that("Phase 17E owns explicit borrowed and mutable-state ownership", {
  types <- blr_source_text("src/blr_mt_default_types.h")
  adapter <- blr_source_text("src/mtblr.cpp")
  expect_match(types, "const std::vector<std::vector<double>>& wy", fixed = TRUE)
  expect_match(types, "std::vector<std::vector<double>> b;", fixed = TRUE)
  expect_match(types, "arma::mat B;", fixed = TRUE)
  expect_match(types, "arma::mat E;", fixed = TRUE)
  expect_match(types, "std::vector<double> pi;", fixed = TRUE)
  expect_match(adapter,
    "std::move(b), std::move(B), std::move(E), std::move(pi)", fixed = TRUE)
})
