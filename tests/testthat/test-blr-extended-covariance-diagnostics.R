test_that("strict-lower descriptor order is column major", {
  controls <- .blr_convergence_controls(
    "extended", list(warn = FALSE, extended_groups = "covariance"), 2L)
  plan <- .blr_mtblr_extended_plan(
    controls, paste0("m", 1:3), paste0("T", 1:5), "sbayesc",
    list(patterns = list(names = c("null", "active"))), list(),
    TRUE, TRUE, TRUE, "diagonal", 2L, 8L)
  expect_equal(plan$counts[["covariance"]], 30L)
})
