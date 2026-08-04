test_that("CSR operator contracts distinguish full, hard-sparse, and unknown", {
  full <- sblr:::.stblr_csr_operator_contract(list(sparseLD = list(
    source = "make_sparse_ld", prefix = "full", r2_threshold = 0,
    max_distance_variants = 0L, max_distance_bp = 0, block_size = 64L)))
  sparse <- sblr:::.stblr_csr_operator_contract(list(sparseLD = list(
    source = "make_sparse_ld", prefix = "sparse", r2_threshold = 0.001,
    max_distance_variants = 1000L, max_distance_bp = 0,
    block_size = 1024L)))
  unknown <- sblr:::.stblr_csr_operator_contract(NULL, "external")

  expect_identical(full$operator_role, "summary_statistics_reference")
  expect_identical(full$operator_contract, "exact_full_csr_v1")
  expect_false(full$operator_approximate)
  expect_identical(sparse$operator_role, "approximate_summary_statistics")
  expect_identical(sparse$operator_contract, "hard_sparse_csr_v1")
  expect_true(sparse$operator_approximate)
  expect_match(sparse$limitation, "semidefiniteness is not sufficient")
  expect_identical(unknown$operator_role, "supplied_csr_unclassified")
  expect_true(is.na(unknown$operator_approximate))
  expect_identical(unknown$fidelity_evidence, "unavailable")
})
