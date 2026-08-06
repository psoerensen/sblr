test_that("BED BayesRC compact traces use the BayesR stick orientation", {
  fixture <- aggregate_bayesr_fixture()
  on.exit(unlink(fixture$path), add = TRUE)
  annotation <- cbind(intercept = 1, enriched = c(0, 1, 0, 1))
  rownames(annotation) <- fixture$ids
  fit <- stblr_bed(
    y = fixture$y, Glist = fixture$Glist, method = "bayesrc",
    annotation = annotation, mixture_var = c(0, .01, .1, 1),
    add_intercept = FALSE, updateAlpha = FALSE, updateE = FALSE,
    nit = 8L, nburn = 2L, nchains = 2L, ncores = 1L,
    keep_chains = TRUE, seed = 806L, chain_seeds = c(1806L, 2806L),
    convergence = "extended", convergence_control = list(
      warn = FALSE, extended_groups = c("annotations", "probability"),
      keep_traces = TRUE, aggregate_component_states = TRUE,
      selected_markers = fixture$ids,
      selected_marker_quantities = "component"))
  expect_native_aggregate_oracle(fit, 4L, 4L)
  expect_true(all(fit$stick_eligible_count_trace[, , 1L] == 4L))
  expect_identical(
    as.integer(fit$stick_continue_count_trace[, , 1L]),
    as.integer(fit$realized_active_count_trace))
})

test_that("pooled output cannot discard compact trace chain boundaries", {
  fixture <- aggregate_bayesr_fixture()
  on.exit(unlink(fixture$path), add = TRUE)
  expect_error(stblr_bed(
    y = fixture$y, Glist = fixture$Glist, method = "bayesr",
    mixture_var = c(0, .01, .1, 1), nit = 6L, nburn = 4L, nthin = 2L,
    nchains = 1L, ncores = 1L, keep_chains = FALSE, updateE = FALSE,
    seed = 807L, convergence = "extended", convergence_control = list(
      warn = FALSE, extended_groups = "probability", keep_traces = TRUE,
      aggregate_component_states = TRUE, selected_markers = fixture$ids,
      selected_marker_quantities = "component")),
    "requires keep_chains = TRUE")
})
