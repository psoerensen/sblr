test_that("MT BayesR probability diagnostics use named marginals and opt-in joint states", {
  fixture <- blr_unified_fixture()
  on.exit(blr_unified_cleanup(fixture), add = TRUE)
  prefix <- blr_unified_exact_ld_prefix(fixture$dosage)
  on.exit(blr_unified_cleanup_prefix(prefix), add = TRUE)
  metadata <- fixture$stats$marker_metadata
  metadata$effect_allele <- "A"; metadata$other_allele <- "C"
  stats <- fixture$stats; stats$marker_metadata <- metadata
  fit <- mtblr_csr(
    stats = stats, ld_prefix = prefix,
    ld_metadata = list(prefix = prefix, marker_ids = stats$marker_names,
      marker_metadata = metadata, scale = "standardized_genotype",
      source = "make_summary_stats"), method = "sbayesr",
    mixture_var = c(0, .1, 1), models = matrix(c(0L, 1L), 2L, 1L),
    joint_pi = c(.7, .15, .15), joint_pi_prior = rep(1, 3),
    vb = matrix(.1), ve = matrix(.5), updateB = FALSE, updateE = FALSE,
    nit = 8L, nburn = 2L, nchains = 2L, ncores = 1L,
    convergence = "extended", convergence_control = list(
      warn = FALSE, extended_groups = "probability",
      full_probability_states = TRUE, keep_traces = TRUE))
  probability <- fit$convergence$summary[
    fit$convergence$summary$tier == 2L, , drop = FALSE]
  expect_identical(probability$group,
    c(rep("component_pi", 3L), "pattern_pi", rep("joint_pi", 3L)))
  expect_identical(probability$component_name[1:3],
    paste0("component_", 0:2))
  expect_identical(probability$pattern_name[5:7],
    fit$model_parameters$mixture$joint_state_names)
  values <- fit$convergence_traces$values[, , probability$quantity, drop = FALSE]
  expect_true(all(is.finite(values)))
})

test_that("BayesRC probability scope excludes marker-average component priors", {
  x <- .mt_bayesrc_fixture(); on.exit(.mt_bayesrc_cleanup(x), add = TRUE)
  common <- .mt_bayesrc_common(); common$annotations <- x$annotations
  common$nchains <- 2L; common$ncores <- 1L
  common$convergence <- "extended"
  common$convergence_control <- list(
    warn = FALSE, extended_groups = "probability")
  fit <- do.call(mtblr_csr, c(list(
    stats = x$stats, ld_prefix = x$prefix, ld_metadata = x$ld_metadata),
    common))
  expect_true(any(fit$convergence$summary$group == "pattern_pi"))
  expect_false(any(fit$convergence$summary$group == "component_pi"))
})

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

test_that("MT binary BayesC captures one active probability without its complement", {
  fixture <- blr_unified_fixture()
  on.exit(blr_unified_cleanup(fixture), add = TRUE)
  prefix <- blr_unified_exact_ld_prefix(fixture$dosage)
  on.exit(blr_unified_cleanup_prefix(prefix), add = TRUE)
  stats <- fixture$stats
  stats$marker_metadata$effect_allele <- "A"
  stats$marker_metadata$other_allele <- "C"
  metadata <- list(prefix = prefix, marker_ids = stats$marker_names,
    marker_metadata = stats$marker_metadata,
    scale = "standardized_genotype", source = "make_summary_stats")
  fit <- mtblr_csr(
    stats, ld_prefix = prefix, ld_metadata = metadata,
    models = matrix(c(0L, 1L), 2L, 1L), pimodels = c(.8, .2),
    updateB = FALSE, updateE = FALSE, updatePi = TRUE,
    nit = 8L, nburn = 2L, nchains = 2L, ncores = 1L,
    convergence = "extended", convergence_control = list(
      warn = FALSE, extended_groups = "probability", keep_traces = TRUE))
  probability <- fit$convergence$summary[
    fit$convergence$summary$tier == 2L, , drop = FALSE]
  expect_identical(probability$group, "pi_active")
  expect_length(probability$diagnostic_key, 1L)
  expect_identical(dim(fit$convergence_traces$values)[[3L]], 6L)
})
