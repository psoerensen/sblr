test_that("permanent BLR contract ownership and tiers are documented", {
  ownership <- blr_source_text("docs/dev/blr_test_contract_ownership.md")
  for (term in c("ordinary CSR BayesC", "ordinary CSR BayesR",
      "ordinary CSR SBayesRC", "fixed-prior CSR BayesC",
      "group CSR BayesC", "learned-annotation CSR BayesC",
      "scheduled CSR BayesC", "scheduled packed-BED BayesC",
      "packed-BED BayesR", "packed-BED BayesRC",
      "corrected dense MT BayesC", "Tier 1", "Tier 2", "Tier 3"))
    expect_match(ownership, term, fixed = TRUE)
})

test_that("extended reproducibility has one stable opt-in owner", {
  extended <- blr_source_text("tests/testthat/test-blr-extended-reproducibility.R")
  workflow <- blr_source_text(".github/workflows/blr-framework-extended.yml")
  expect_source_count("SBLR_RUN_EXTENDED_REPRODUCIBILITY", extended, 1L)
  expect_source_count("SBLR_RUN_EXTENDED_REPRODUCIBILITY", workflow, 1L)
  expect_source_forbidden(workflow, c("SBLR_RUN_PHASE11A_FRESH",
    "SBLR_RUN_PHASE11B_FRESH", "SBLR_RUN_PHASE17B_FRESH",
    "SBLR_RUN_PHASE17C_FRESH", "SBLR_RUN_PHASE17D_FRESH",
    "SBLR_RUN_PHASE17E_FRESH", "SBLR_RUN_PHASE17F_FRESH"))
  expect_match(workflow, "SBLR_RUN_PEAK_RSS", fixed = TRUE)
})

test_that("historical and corrected MT fixture hashes remain permanent", {
  paths <- c(
    file.path("tests/testthat/fixtures/blr_phase17b_mt_default",
      sprintf("config-%d.rds", 1:3)),
    file.path("tests/testthat/fixtures/blr_phase17c_mt_default_corrected",
      sprintf("config-%d.rds", 1:3)))
  blr_expect_fixture_md5(paths, c(
    "2b820e1cdd9e731f1f0ffea442ef4e53",
    "d6d8abec35168a088a9accab87b3c6d0",
    "48fe8040d52b1d23c6e7437d632ebebf",
    "d64c29b872546006c0dfb0303403abce",
    "0d699c71e02348b39167bc1695b87b5e",
    "087ac6e06a68763a9ac558c5b21b8f3a"))
})

