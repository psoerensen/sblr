test_that("official block Ve boundaries are exact", {
  expect_false(sblr:::.stblr_block_ve_decision(
    "fixVe", 0L, 1, 2, 1.1)$selected)
  expect_true(sblr:::.stblr_block_ve_decision(
    "samVe", 0L, 0, 0, 1.1)$selected)
  below <- sblr:::.stblr_block_ve_decision(
    "allMixVe", 0L, 1, 1.1, 1.1)
  above <- sblr:::.stblr_block_ve_decision(
    "allMixVe", 0L, 1, 1.1000001, 1.1)
  expect_false(below$selected)
  expect_true(above$selected)
  expect_false(sblr:::.stblr_block_ve_decision(
    "allMixVe", 0L, 1e-9, 1, 1.1)$selected)
  expect_false(sblr:::.stblr_block_ve_decision(
    "allMixVe", 0L, 1, 1e-9, 1.1)$selected)
  expect_true(sblr:::.stblr_block_ve_decision(
    "allMixVe", 0L, 1e-8, 2e-8, 1.1)$selected)

  expect_false(sblr:::.stblr_block_ve_decision(
    "mixVe", 49L, 1, 2, 1.1)$selected)
  at_trial <- sblr:::.stblr_block_ve_decision(
    "mixVe", 50L, 1, 2, 1.1)
  expect_true(at_trial$selected)
  expect_true(at_trial$permanent_selected)
  expect_true(sblr:::.stblr_block_ve_decision(
    "mixVe", 51L, 1, 0, 1.1, TRUE)$selected)

  expect_false(sblr:::.stblr_block_ve_ratio_accepted(0.7, 0.7))
  expect_true(sblr:::.stblr_block_ve_ratio_accepted(
    0.7 + .Machine$double.eps, 0.7))

  pars <- sblr:::.stblr_block_ve_draw_parameters(13, 2.5, 4, 7L)
  expect_equal(pars$prior_scale, 1.25)
  expect_equal(pars$scale, 18)
  expect_equal(pars$degrees_freedom, 11)
})

test_that("retained factor scale maps to official block units", {
  fixture <- aggregate_block_fixture()
  on.exit(unlink(fixture$path), add = TRUE)
  beta <- c(0.2, -0.1, 0.05, 0)
  contract <- sblr:::stblr_block_low_rank_contract_internal(
    fixture$path, fixture$Glist$n, list(seq_along(fixture$ids)), NULL,
    fixture$Glist$maf[[1L]], c(0L, 2L),
    matrix(fixture$stats$wy[[1L]], nrow = 1L), beta, .999999,
    yy = fixture$stats$yy[[1L]])

  for (block in seq_along(contract$factor)) {
    q_s <- contract$factor[[block]]
    start <- c(1L, 3L)[[block]]
    size <- ncol(q_s)
    idx <- seq.int(start, length.out = size)
    r_s <- contract$residual[
      contract$residual_offset[[block]] + seq_len(nrow(q_s))]
    q_official <- q_s / sqrt(fixture$Glist$n)
    r_official <- r_s / sqrt(fixture$Glist$n)
    d_s <- colSums(q_s^2)
    d_official <- colSums(q_official^2)
    expect_equal(
      fixture$Glist$n * sum(r_official^2),
      contract$block_residual_norm_squared[[block]], tolerance = 1e-10)
    expect_equal(
      sum((q_official %*% beta[idx])^2),
      contract$block_genetic_variance[[block]], tolerance = 1e-10)
    expect_equal(
      sum(beta[idx]^2), contract$block_effect_ss[[block]], tolerance = 1e-12)
    expect_equal(
      as.numeric(crossprod(q_s, r_s) + d_s * beta[idx]),
      fixture$Glist$n * as.numeric(
        crossprod(q_official, r_official) + d_official * beta[idx]),
      tolerance = 1e-10)
    ve <- 1.7
    prior_precision <- seq_along(idx) + 2
    expect_equal(
      d_s / ve + prior_precision,
      fixture$Glist$n * d_official / ve + prior_precision,
      tolerance = 1e-10)
  }

  c_y <- 1.4
  q_s <- contract$factor[[1L]]
  r_s <- contract$residual[seq_len(nrow(q_s))]
  beta_s <- beta[seq_len(ncol(q_s))]
  q_official <- q_s / sqrt(fixture$Glist$n)
  r_official <- r_s / (c_y * sqrt(fixture$Glist$n))
  beta_official <- beta_s / c_y
  expect_equal(sum(r_s^2), c_y^2 * fixture$Glist$n * sum(r_official^2),
               tolerance = 1e-10)
  expect_equal(sum((q_s %*% beta_s)^2) / fixture$Glist$n,
               c_y^2 * sum((q_official %*% beta_official)^2),
               tolerance = 1e-10)
})

test_that("public retained SBayesR defaults to fixed-compatible block Ve when updateE is false", {
  fixture <- aggregate_block_fixture()
  on.exit(unlink(fixture$path), add = TRUE)
  fit <- stblr_block_eigen(
    fixture$stats, fixture$Glist, c(1L, 3L), method = "sbayesr",
    representation = "low_rank", eigen_prop = .999999,
    updateE = FALSE, nit = 5L, nburn = 2L, seed = 818L,
    convergence = "none")
  expect_identical(fit$input$residual_policy, "gctb_block")
  expect_identical(fit$input$block_ve_mode, "fixVe")
  expect_true(all(is.finite(fit$ves)))
  expect_equal(as.numeric(fit$ves), rep(fit$input$phenotype_variance, 7L))
  expect_identical(
    fit$input$heritability_definition,
    "sum of block genetic variances divided by phenotype variance")
  expect_true(is.list(fit$block_ve))
  expect_true(all(fit$block_ve$resampled_per_chain_block == 0))
})

test_that("block Ve diagnostics do not consume sampler RNG", {
  fixture <- aggregate_block_fixture()
  on.exit(unlink(fixture$path), add = TRUE)
  base <- list(
    stats = fixture$stats, Glist = fixture$Glist,
    block_start = c(1L, 3L), method = "sbayesr",
    representation = "low_rank", eigen_prop = .999999,
    residual_policy = "fixed_block", nit = 5L, nburn = 2L,
    seed = 819L, convergence = "none")
  compact <- do.call(stblr_block_eigen, c(base, list(block_ve_keep_history = FALSE)))
  extended <- do.call(stblr_block_eigen, c(base, list(block_ve_keep_history = TRUE)))
  for (field in c("b", "d", "bm", "dm", "vbs", "vgs", "ves", "pis"))
    expect_identical(compact[[field]], extended[[field]])
  expect_null(compact$block_ve$history)
  expect_true(is.list(extended$block_ve$history))
})

test_that("legacy projected policy remains explicit and rejects ambiguous controls", {
  fixture <- aggregate_block_fixture()
  on.exit(unlink(fixture$path), add = TRUE)
  expect_error(stblr_block_eigen(
    fixture$stats, fixture$Glist, c(1L, 3L), method = "sbayesr",
    representation = "low_rank", adjE = 0.9, nit = 2L, nburn = 1L),
    "adjE is only available")
  expect_error(stblr_block_eigen(
    fixture$stats, fixture$Glist, c(1L, 3L), method = "sbayesr",
    representation = "dense_reconstructed", residual_policy = "gctb_block",
    low_rank_residual_rebuild_every = 0L, nit = 2L, nburn = 1L),
    "require representation")
})
