test_that("shared scalar convergence engine retains the canonical mathematics", {
  x <- cbind(
    c(-2, -1, 0, 1, 2, -1, 0, 1),
    c(-1.8, -.8, .2, 1.2, 2.2, -.8, .2, 1.2),
    c(-2.1, -.9, .1, .9, 1.9, -1.1, .1, .9),
    c(-1.9, -1.1, -.1, 1.1, 2.1, -.9, -.1, 1.1))
  scalar <- sblr:::.blr_convergence_scalar(x)
  expect_true(all(c("rhat", "rhat_rank", "rhat_folded", "ess_bulk",
                    "ess_tail", "ess_mean", "posterior_sd", "mcse_mean") %in%
                    names(scalar)))
  expect_equal(scalar$rhat,
               max(scalar$rhat_rank, scalar$rhat_folded), tolerance = 0)
  expect_equal(scalar$mcse_mean,
               scalar$posterior_sd / sqrt(scalar$ess_mean), tolerance = 1e-12)
})

test_that("ST adapters use unpooled post-burn chain traces", {
  trace <- function(offset) list(
    vbs = 1:8 + offset, vgs = 2:9 + offset, ves = 3:10 + offset,
    vle = 4:11 + offset, vld = 5:12 + offset)
  chains <- list(T1 = list(trace(0), trace(.25)))
  bundle <- sblr:::.blr_st_convergence_bundle(
    chains, "T1", "bayesc", "csr", nit = 6L, nburn = 2L)
  expect_identical(dim(bundle$values), c(6L, 2L, 5L))
  expect_equal(bundle$values[, 1, 1], 3:8)
  expect_equal(bundle$values[, 2, 1], 3:8 + .25)
  expect_identical(bundle$schema$class, "blr_convergence_trace_bundle")
})

test_that("scheduled ST CSR convergence uses task-private post-burn traces", {
  fixture <- blr_unified_scheduled_csr_fixture(nt = 2L)
  on.exit(blr_unified_cleanup_csr(fixture), add = TRUE)
  call <- function(nthin, keep_chains, ncores = 1L) {
    stblr_csr(
      stats = fixture$stats, ld_prefix = fixture$prefix,
      method = "sbayesc", scheduled = TRUE,
      pi_init = .35, pi_prior_mean = .35, pi_prior_strength = 3,
      updateB = FALSE, updateE = FALSE, updatePi = FALSE,
      nit = 8L, nburn = 2L, nthin = nthin, seed = 1001L,
      nchains = 2L, ncores = ncores, chain_seeds = c(1101L, 1102L),
      keep_chains = keep_chains, convergence = "core",
      convergence_control = list(warn = FALSE, keep_traces = TRUE),
      full_sweep_every = 3L, null_skip_base = 2L, null_skip_max = 7L,
      candidate_threshold = .12, candidate_lifetime = 2L,
      skip_nulls_burnin_only = FALSE, wakeup_ld_neighbors = TRUE,
      wakeup_diff_threshold = .02, wakeup_max_neighbors = 2L,
      updateLDswap = FALSE)
  }
  base <- call(1L, FALSE)
  retained <- call(1L, TRUE)
  thinned <- call(2L, FALSE)
  parallel <- call(1L, FALSE, 2L)
  expect_identical(base$convergence$summary, retained$convergence$summary)
  expect_identical(base$convergence$overview, retained$convergence$overview)
  expect_identical(base$convergence_traces, retained$convergence_traces)
  expect_identical(base$convergence_traces, thinned$convergence_traces)
  expect_identical(base$convergence$summary, thinned$convergence$summary)
  expect_identical(base$convergence_traces, parallel$convergence_traces)
  expect_identical(base$convergence$summary, parallel$convergence$summary)
  expect_null(base$chains)
  expect_length(retained$chains, 4L)
  expect_identical(dim(base$convergence_traces$values), c(8L, 2L, 10L))
  expect_identical(base$convergence_traces$postburn_draws_per_chain, 8L)
  expect_false(base$convergence$trace_retention$burnin_included)
  expect_false(base$convergence$trace_retention$additional_thinning)
})

test_that("MT packed-BED exposes chain-private five-trace convergence", {
  case <- blr_bed_public_case(nt = 2L)
  on.exit(blr_bed_public_cleanup(case), add = TRUE)
  fit <- do.call(mtblr_bed, blr_public_convergence_public_args(
    case, nchains = 2L, ncores = 1L, convergence = "core",
    warn = FALSE, keep_traces = TRUE, nit = 8L, nburn = 2L))
  expect_equal(fit$vld, fit$vgs - fit$vle, tolerance = 1e-12)
  expect_identical(unique(fit$convergence$summary$group),
                   c("vbs", "vgs", "ves", "vle", "vld"))
  expect_identical(dim(fit$convergence_traces$values),
                   c(8L, 2L, 10L))
  expect_false(fit$convergence$trace_retention$burnin_included)
  expect_false(fit$convergence$trace_retention$additional_thinning)
})
