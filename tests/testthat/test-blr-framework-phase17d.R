test_that("Phase 17D extraction history has one active typed successor", {
  core <- blr_source_text("src/blr_mt_default_core_impl.h")
  adapter <- blr_mt_public_source()
  expect_source_count("inline MtDefaultCoreResult run_mt_default_core(", core, 1L)
  expect_source_count("run_mt_default_core(", adapter, 1L)
  expect_source_count("for ( int it = 0;", adapter, 0L)
  expect_match(blr_source_text("docs/dev/blr_framework_phase17d_report.md"),
    "MECHANICAL_LINES=240", fixed = TRUE)
})
