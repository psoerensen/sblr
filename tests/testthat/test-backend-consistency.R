check_stblr_consistency <- getFromNamespace("check_stblr_consistency", "sblr")

make_backend_consistency_fit <- function(nchains = 1L, summaries = FALSE) {
  markers <- paste0("m", 1:3)
  dm <- matrix(c(0.10, 0.20, 0.30), ncol = 1,
               dimnames = list(markers, "trait1"))
  bm <- matrix(c(0.01, -0.02, 0.03), ncol = 1,
               dimnames = list(markers, "trait1"))
  fit <- list(dm = dm, bm = bm, input = list(nchains = nchains))
  if (summaries) {
    fit$dm_chain_mean_sd <- matrix(if (nchains == 1L) rep(0, 3) else c(0.01, 0.02, 0.03),
                        ncol = 1, dimnames = dimnames(dm))
    fit$dm_chain_mean_min <- dm - if (nchains == 1L) 0 else 0.01
    fit$dm_chain_mean_max <- dm + if (nchains == 1L) 0 else 0.01
    fit$bm_chain_mean_sd <- matrix(if (nchains == 1L) rep(0, 3) else c(0.001, 0.002, 0.003),
                        ncol = 1, dimnames = dimnames(bm))
    fit$bm_chain_mean_min <- bm - if (nchains == 1L) 0 else 0.001
    fit$bm_chain_mean_max <- bm + if (nchains == 1L) 0 else 0.001
  }
  fit
}

test_that("backend consistency checker accepts minimal single-chain fit", {
  chk <- check_stblr_consistency(
    make_backend_consistency_fit(nchains = 1L),
    verbose = FALSE
  )

  expect_s3_class(chk, "stblr_backend_check")
  expect_true(chk$ok)
  expect_equal(chk$metadata$nchains, 1L)
  expect_false(chk$metadata$has_chain_summaries)
})

test_that("backend consistency checker accepts valid multi-chain summaries", {
  chk <- check_stblr_consistency(
    make_backend_consistency_fit(nchains = 2L, summaries = TRUE),
    verbose = FALSE
  )

  expect_true(chk$ok)
  expect_true(chk$metadata$has_chain_summaries)
})

test_that("backend consistency checker handles missing multi-chain summaries", {
  fit <- make_backend_consistency_fit(nchains = 2L, summaries = FALSE)

  chk_default <- check_stblr_consistency(fit, verbose = FALSE)
  expect_false(chk_default$ok)
  expect_true("chain_summaries.required" %in% chk_default$checks$check)

  chk_optional <- check_stblr_consistency(
    fit,
    require_chain_summaries = FALSE,
    verbose = FALSE
  )
  expect_true(chk_optional$ok)
})

test_that("backend consistency checker catches invalid dimensions", {
  fit <- make_backend_consistency_fit()
  fit$bm <- fit$bm[1:2, , drop = FALSE]

  chk <- check_stblr_consistency(fit, verbose = FALSE)

  expect_false(chk$ok)
  expect_false(chk$checks$ok[match("dimensions.dm_bm", chk$checks$check)])
})

test_that("backend consistency checker catches invalid min and max summaries", {
  fit <- make_backend_consistency_fit(nchains = 2L, summaries = TRUE)
  fit$dm_chain_mean_min[1, 1] <- fit$dm[1, 1] + 0.1
  fit$bm_chain_mean_max[2, 1] <- fit$bm[2, 1] - 0.1

  chk <- check_stblr_consistency(fit, verbose = FALSE)

  expect_false(chk$ok)
  expect_false(chk$checks$ok[match("chain_summary.dm.bounds", chk$checks$check)])
  expect_false(chk$checks$ok[match("chain_summary.bm.bounds", chk$checks$check)])
})

test_that("backend consistency checker validates single-chain summaries", {
  chk <- check_stblr_consistency(
    make_backend_consistency_fit(nchains = 1L, summaries = TRUE),
    require_chain_summaries = TRUE,
    verbose = FALSE
  )

  expect_true(chk$ok)
  expect_true(chk$checks$ok[match("chain_summary.single.dm_chain_mean_sd_zero", chk$checks$check)])
  expect_true(chk$checks$ok[match("chain_summary.single.bm_equal", chk$checks$check)])
})

test_that("backend consistency checker validates compact chain summaries", {
  fit <- make_backend_consistency_fit(nchains = 2L, summaries = TRUE)
  fit$dm <- matrix(c(0.15, 0.25, 0.35), ncol = 1,
                   dimnames = dimnames(fit$dm))
  fit$bm <- matrix(c(0.015, -0.025, 0.035), ncol = 1,
                   dimnames = dimnames(fit$bm))
  fit$chains <- list(
    task1 = list(trait_index = 1L, chain_index = 1L,
                 dm = c(0.10, 0.20, 0.30), bm = c(0.01, -0.02, 0.03)),
    task2 = list(trait_index = 1L, chain_index = 2L,
                 dm = c(0.20, 0.30, 0.40), bm = c(0.02, -0.03, 0.04))
  )
  fit$dm_chain_mean_min <- pmin(fit$chains$task1$dm, fit$chains$task2$dm)
  fit$dm_chain_mean_max <- pmax(fit$chains$task1$dm, fit$chains$task2$dm)
  fit$bm_chain_mean_min <- pmin(fit$chains$task1$bm, fit$chains$task2$bm)
  fit$bm_chain_mean_max <- pmax(fit$chains$task1$bm, fit$chains$task2$bm)
  for (nm in c("dm_chain_mean_min", "dm_chain_mean_max",
               "bm_chain_mean_min", "bm_chain_mean_max")) {
    fit[[nm]] <- matrix(fit[[nm]], ncol = 1, dimnames = dimnames(fit$dm))
  }

  chk <- check_stblr_consistency(
    fit,
    require_chains = TRUE,
    verbose = FALSE
  )

  expect_true(chk$ok)
  expect_true(chk$checks$ok[match("chains.trait1.dm_mean", chk$checks$check)])
  expect_true(chk$checks$ok[match("chains.trait1.bm_mean", chk$checks$check)])
})

test_that("backend consistency checker validates LD-swap diagnostics", {
  fit <- make_backend_consistency_fit()
  fit$diagnostics$ld_swap <- data.frame(
    attempted = 10,
    accepted = 4,
    acceptance_rate = 0.4
  )

  chk <- check_stblr_consistency(
    fit,
    require_ld_swap = TRUE,
    verbose = FALSE
  )

  expect_true(chk$ok)
  expect_true(chk$checks$ok[match("ld_swap.counts", chk$checks$check)])
})

test_that("backend consistency checker catches invalid LD-swap diagnostics", {
  fit <- make_backend_consistency_fit()
  fit$diagnostics$ld_swap <- data.frame(
    attempted = 2,
    accepted = 3,
    acceptance_rate = 1.5
  )

  chk <- check_stblr_consistency(
    fit,
    require_ld_swap = TRUE,
    verbose = FALSE
  )

  expect_false(chk$ok)
  expect_false(chk$checks$ok[match("ld_swap.counts", chk$checks$check)])
  expect_false(chk$checks$ok[match("ld_swap.acceptance_rate", chk$checks$check)])
})
