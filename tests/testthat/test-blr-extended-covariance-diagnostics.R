test_that("MT covariance diagnostics use strict lower order and structural E", {
  fixture <- blr_unified_fixture()
  on.exit(blr_unified_cleanup(fixture), add = TRUE)
  y <- cbind(as.numeric(fixture$y),
             as.numeric(fixture$y) * 0.8 + seq_along(fixture$y) / 100)
  colnames(y) <- c("T1", "T2")
  fit <- mtblr_bed(
    y = y, Glist = fixture$Glist, method = "bayesc",
    residual_covariance = "diagonal", updateB = FALSE, updateE = FALSE,
    nit = 8L, nburn = 2L, nchains = 2L, ncores = 1L,
    convergence = "extended",
    convergence_control = list(
      warn = FALSE, extended_groups = "covariance", keep_traces = TRUE))
  covariance <- fit$convergence$summary[
    fit$convergence$summary$tier == 2L, , drop = FALSE]
  expect_identical(covariance$quantity,
    c("cov_b[T2,T1]", "cov_g[T2,T1]", "cov_e[T2,T1]"))
  expect_identical(covariance$trait, rep("T1", 3L))
  expect_identical(covariance$trait2, rep("T2", 3L))
  expect_identical(covariance$status,
    c("not_updated", covariance$status[[2L]], "structural_zero"))
  expect_false(covariance$captured[[3L]])
  expect_true(all(fit$convergence_traces$values[, , 13L] == 0))
})

test_that("strict-lower descriptor order is column major", {
  controls <- .blr_convergence_controls(
    "extended", list(warn = FALSE, extended_groups = "covariance"), 2L)
  plan <- .blr_mtblr_extended_plan(
    controls, paste0("m", 1:3), paste0("T", 1:5), "sbayesc",
    list(patterns = list(names = c("null", "active"))), list(),
    TRUE, TRUE, TRUE, "diagonal", 2L, 8L)
  expect_equal(plan$counts[["covariance"]], 30L)
})
