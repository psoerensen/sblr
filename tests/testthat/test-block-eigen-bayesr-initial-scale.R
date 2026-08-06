test_that("block-eigen SBayesR compact allocation traces equal the oracle", {
  fixture <- aggregate_block_fixture()
  on.exit(unlink(fixture$path), add = TRUE)
  fit <- stblr_block_eigen(
    stats = fixture$stats, Glist = fixture$Glist,
    block_start = c(1L, 3L), method = "sbayesr",
    representation = "low_rank", eigen_prop = .999999,
    mixture_var = c(0, .01, .1, 1), updateE = FALSE,
    nit = 7L, nburn = 2L, nchains = 2L, ncores = 1L,
    keep_chains = TRUE, seed = 808L, chain_seeds = c(1808L, 2808L),
    convergence = "extended", convergence_control = list(
      warn = FALSE, extended_groups = "probability", keep_traces = TRUE,
      aggregate_component_states = TRUE, selected_markers = fixture$ids,
      selected_marker_quantities = "component"))
  expect_native_aggregate_oracle(fit, 4L, 4L)
  expect_identical(fit$input$operator_representation, "low_rank")
})

test_that("omitted cross-block quadratic can make global projected SSE invalid", {
  dosage <- rbind(
    c(0, 1, 2, 0, 1, 2, 1, 0),
    c(2, 0, 1, 2, 0, 1, 2, 1),
    c(0, 1, 2, 0, 1, 2, 1, 0),
    c(2, 0, 1, 2, 0, 1, 2, 1))
  path <- tempfile(fileext = ".bed")
  on.exit(unlink(path), add = TRUE)
  write_aggregate_bayesr_bed(path, dosage)
  af <- rowMeans(dosage) / 2
  X <- t((dosage - 2 * af) / sqrt(2 * af * (1 - af)))
  beta <- c(1, 0, 1, 0)
  y <- drop(X %*% beta)
  score <- drop(crossprod(X, y))
  contract <- sblr:::stblr_block_low_rank_contract_internal(
    path, ncol(dosage), list(seq_len(nrow(dosage))), NULL, af,
    c(0L, 2L), matrix(score, nrow = 1L), beta, .999999,
    yy = sum(y^2))
  direct_sse <- sum((y - X %*% beta)^2)
  expect_equal(direct_sse, 0, tolerance = 1e-12)
  expect_lt(contract$residual_sse, -1)
  expect_equal(
    contract$residual_sse,
    sum(y^2) - 2 * sum(beta * score) + contract$quadratic_form,
    tolerance = 2e-5)
  expect_gt(sum((X %*% beta)^2) - contract$quadratic_form, 1)
})
