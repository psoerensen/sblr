test_that("ST BayesC captures one active probability without its complement", {
  chains <- list(T1 = list(
    chain1 = list(pis = seq(.1, .5, length.out = 10)),
    chain2 = list(pis = seq(.2, .6, length.out = 10))))
  controls <- .blr_convergence_controls(
    "extended", list(warn = FALSE, extended_groups = "probability"), 2L)
  bundle <- .blr_st_extended_bundle(
    chains, "T1", "sbayesc", "csr", 8L, 2L,
    list(input = list(prior_kernel = "bayesc", updatePi = TRUE)), controls)
  expect_identical(bundle$quantities$group, "pi_active")
  expect_length(bundle$quantities$diagnostic_key, 1L)
  expect_equal(bundle$values[, 1L, 1L], chains$T1$chain1$pis[3:10])
})
