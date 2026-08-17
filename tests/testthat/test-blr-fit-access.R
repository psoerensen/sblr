fit_access_fixture <- function(family = c("st", "mt")) {
  family <- match.arg(family)
  markers <- paste0("m", 1:3)
  traits <- if (family == "st") "trait1" else c("trait1", "trait2")
  chains <- "chain1"
  draws <- c("draw1", "draw2")
  realised <- array(
    seq_len(2L * length(markers) * length(traits)) / 100,
    dim = c(2L, 1L, length(markers), length(traits)),
    dimnames = list(draws, chains, markers, traits))
  covariance <- array(
    rep(diag(length(traits)), each = 2L),
    dim = c(2L, 1L, length(traits), length(traits)),
    dimnames = list(draws, chains, traits, traits))
  final_effect <- realised[2L, , , , drop = FALSE]
  dim(final_effect) <- c(1L, length(markers), length(traits))
  dimnames(final_effect) <- list(chains, markers, traits)
  final_covariance <- covariance[2L, , , , drop = FALSE]
  dim(final_covariance) <- c(1L, length(traits), length(traits))
  dimnames(final_covariance) <- list(chains, traits, traits)
  raw <- structure(list(
    schema = list(version = 2L),
    input = list(
      schema = list(
        seed_contract_version = 1L,
        retention_contract_version = 1L,
        execution_contract_version = 1L),
      data = list(
        providers = NULL, operator_resources = NULL,
        provider_maps = NULL, likelihood_regime = NULL),
      mcmc = list(
        task_seeds = stats::setNames(123, chains),
        retained_transition_indices = c(2L, 4L)),
      compute = list(scheduler_version = 1L)),
    draws = list(
      realised_effects = realised, latent_effects = realised,
      marker_covariance = covariance, residual_covariance = NULL,
      convergence = NULL),
    final = list(
      realised_effects = final_effect, latent_effects = final_effect,
      marker_covariance = final_covariance, residual_covariance = NULL),
    diagnostics = list(
      convergence = NULL,
      workers = list(
        requested_cores = 1L, configured_workers = 1L,
        actual_team_size = 1L, task_worker_ids = 0L,
        scheduler_version = 1L, logical_task_order = "chain:0")),
    provenance = list(
      operator_resources = NULL, marker_alignment = NULL,
      seed_contract_version = 1L,
      task_seeds = stats::setNames(123, chains))),
    class = c("blr_raw", "list"))
  covariance_mean <- diag(length(traits))
  dimnames(covariance_mean) <- list(traits, traits)
  fit <- structure(list(
    bm = apply(realised, c(3L, 4L), mean),
    dm = matrix(.5, length(markers), length(traits),
                dimnames = list(markers, traits)),
    cov_b_mean = covariance_mean,
    cov_e_mean = NULL, predictions = NULL,
    diagnostics = raw$diagnostics,
    input = list(
      logical_task_ids = "chain:0",
      task_seeds_resolved = stats::setNames(123, chains),
      retained_transition_indices = c(2L, 4L),
      convergence_iteration_indices = 1:4,
      seed_contract_version = 1L,
      retention_contract_version = 1L,
      scheduler_version = 1L),
    data = list(providers = NULL, operator_resources = NULL),
    convergence = NULL, convergence_traces = NULL),
    class = c(if (family == "st") "stblr_fit" else "mtblr_fit",
              "blr_fit", "list"))
  attr(fit, "blr_raw") <- raw
  fit
}

test_that("canonical extraction preserves ST and MT axes", {
  st <- fit_access_fixture("st")
  mt <- fit_access_fixture("mt")
  expect_identical(dim(extract_posterior(st, "pips")), c(3L, 1L))
  expect_identical(
    dim(extract_posterior(st, "realised_effects", "retained")),
    c(2L, 1L, 3L, 1L))
  expect_identical(
    dim(extract_posterior(mt, "realised_effects", "retained",
                          traits = "trait2", chains = "chain1")),
    c(2L, 1L, 3L, 1L))
  expect_null(extract_posterior(mt, "residual_covariance"))
})

test_that("canonical retained summaries use stored post-burn draws", {
  fit <- fit_access_fixture("mt")
  summary <- summarise_posterior(
    fit, quantity = "effect_covariance", traits = "trait1")
  expect_named(summary, c("trait_row", "trait_col", "n", "mean", "sd",
                          "median", "q_lower", "q_upper"))
  expect_identical(summary$n, 2)
  expect_equal(summary$mean, 1)
})

test_that("diagnostics have one canonical extraction route", {
  fit <- fit_access_fixture("st")
  diagnostics <- extract_diagnostics(fit)
  expect_named(diagnostics, c(
    "sampler", "execution", "convergence", "providers", "provenance"))
  expect_identical(diagnostics$sampler, fit$diagnostics)
  expect_identical(diagnostics$execution$task_seeds,
                   stats::setNames(123, "chain1"))
  expect_identical(diagnostics$execution$retained_transition_indices,
                   c(2L, 4L))
  expect_identical(extract_diagnostics(fit, "execution"),
                   diagnostics$execution)
  expect_error(extract_diagnostics(fit, "missing"),
               "Unknown diagnostic namespace")
})
