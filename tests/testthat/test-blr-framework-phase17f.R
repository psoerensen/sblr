test_that("Phase 17F owns the typed binding-neutral finalizer", {
  types <- blr_source_text("src/blr_mt_default_types.h")
  finalizer <- blr_source_text("src/blr_mt_default_finalize_impl.h")
  adapter <- blr_source_text("src/mtblr.cpp")
  expect_source_count("struct MtDefaultFinalResult", types, 1L)
  expect_source_count("inline MtDefaultFinalResult finalize_mt_default_result(",
    finalizer, 1L)
  # Dense, single-chain CSR/block-eigen, Phase 17O BED, and the two aligned
  # Phase 18 native multichain operator adapters share the same finalizer.
  expect_source_count("finalize_mt_default_result(", adapter, 6L)
  expect_source_forbidden(paste(types, finalizer), c("Rcpp", "SEXP", "RObject",
    "NumericVector", "NumericMatrix", "schema_version"))
  expect_source_forbidden(finalizer, c("std::mt19937", "sampleBset(",
    "samplePi(", "sampleE(", "result.resize(20)"))
})

test_that("Phase 17F owns posterior arithmetic and positional separation", {
  finalizer <- blr_source_text("src/blr_mt_default_finalize_impl.h")
  adapter <- blr_mt_public_source()
  legacy <- blr_source_text("src/blr_mt_default_legacy_adapter.h")
  # Phase 19 adds component posterior normalization to the six Phase 17F
  # retained-count divisions without changing their ownership.
  expect_identical(sum(grepl("/.*retained_count",
    strsplit(finalizer, "\n")[[1]])), 7L)
  expect_source_count("MtDefaultLegacyResult result(22);", legacy, 1L)
  expect_source_forbidden(adapter, c("/marker_retained_count",
    "covb_retained_count >", "covg_retained_count >",
    "cove_retained_count >", "pi_retained_count >"))
  for (position in 0:21)
    expect_true(source_match_count(sprintf("result[%d]", position), legacy) >= 1L)
})

test_that("Phase 17F owns shared naming and future operator requirements", {
  naming <- blr_source_text("docs/dev/stblr_backend_naming.md")
  for (term in c("Phase 17F shared ST/MT naming plan", "SparseLdCsrView",
      "one view per trait/study", "Fully independent CSR structures",
      "no artificial union pattern",
      "canonical marker ID and effect-allele orientation",
      "sample-overlap policy"))
    expect_match(naming, term, fixed = TRUE)
  expect_source_forbidden(paste(
    blr_source_text("src/blr_mt_default_types.h"),
    blr_source_text("src/blr_mt_default_finalize_impl.h")),
    c("row_ptr", "column_index", "eigenvectors", "operator_bundle"))
})
