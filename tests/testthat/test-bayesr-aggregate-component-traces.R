test_that("BED BayesR compact allocation traces equal the full-state oracle", {
  fixture <- aggregate_bayesr_fixture()
  on.exit(unlink(fixture$path), add = TRUE)
  fit <- stblr_bed(
    y = fixture$y, Glist = fixture$Glist, method = "bayesr",
    mixture_var = c(0, .01, .1, 1), nit = 9L, nburn = 3L,
    nchains = 2L, ncores = 1L, keep_chains = TRUE, updateE = FALSE,
    seed = 804L, chain_seeds = c(1804L, 2804L),
    convergence = "extended", convergence_control = list(
      warn = FALSE, extended_groups = "probability", keep_traces = TRUE,
      aggregate_component_states = TRUE, selected_markers = fixture$ids,
      selected_marker_quantities = "component"))
  expect_identical(dim(fit$component_count_trace), c(9L, 2L, 4L))
  expect_identical(dim(fit$realized_active_count_trace), c(9L, 2L, 1L))
  expect_identical(dim(fit$stick_eligible_count_trace), c(9L, 2L, 3L))
  expect_native_aggregate_oracle(fit, 4L, 4L)
})

test_that("BED BayesR compact tracing is RNG-neutral", {
  fixture <- aggregate_bayesr_fixture()
  on.exit(unlink(fixture$path), add = TRUE)
  common <- list(
    y = fixture$y, Glist = fixture$Glist, method = "bayesr",
    mixture_var = c(0, .01, .1, 1), nit = 7L, nburn = 2L,
    nchains = 2L, ncores = 1L, keep_chains = TRUE, updateE = FALSE,
    seed = 805L, chain_seeds = c(1805L, 2805L), convergence = "extended")
  run <- function(enabled) do.call(stblr_bed, c(common, list(
    convergence_control = list(
      warn = FALSE, extended_groups = "probability", keep_traces = TRUE,
      aggregate_component_states = enabled))))
  off <- run(FALSE)
  on <- run(TRUE)
  for (field in c("bm", "dm", "b", "component", "vbs", "vgs", "ves",
                  "vle", "vld", "pi_final", "pi_mean",
                  "component_probabilities"))
    expect_equal(on[[field]], off[[field]], tolerance = 0, info = field)
  expect_null(off$component_count_trace)
})
