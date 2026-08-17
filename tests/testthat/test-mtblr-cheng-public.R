phase5b_public_fixture <- function(trait_count) {
  fixture <- blr_unified_fixture()
  base <- fixture$y[, 1L]
  phenotype <- vapply(seq_len(trait_count), function(trait) {
    (0.62 - 0.08 * trait) * base +
      0.015 * trait * (seq_along(base) - mean(seq_along(base)))
  }, numeric(length(base)))
  colnames(phenotype) <- paste0("trait", seq_len(trait_count))
  rownames(phenotype) <- rownames(fixture$y)
  marker_covariance <- outer(
    seq_len(trait_count), seq_len(trait_count),
    function(i, j) 0.04^(abs(i - j) + 1L))
  diag(marker_covariance) <- seq(0.42, 0.30, length.out = trait_count)
  residual_covariance <- outer(
    seq_len(trait_count), seq_len(trait_count),
    function(i, j) 0.08 / (1 + abs(i - j)))
  diag(residual_covariance) <- seq(1.05, 0.78, length.out = trait_count)
  dimnames(marker_covariance) <- dimnames(residual_covariance) <-
    list(colnames(phenotype), colnames(phenotype))
  fixture$phenotype <- phenotype
  fixture$Vb <- marker_covariance
  fixture$Ve <- residual_covariance
  fixture$Psi_b <- diag(0.55, trait_count)
  fixture$Psi_e <- diag(0.65, trait_count)
  dimnames(fixture$Psi_b) <- dimnames(fixture$Psi_e) <-
    list(colnames(phenotype), colnames(phenotype))
  fixture
}

phase5b_public_fit <- function(fixture, sampled = FALSE, chains = 1L,
                               cores = 1L, seed = 51,
                               chain_seeds = NULL, method = "bayesc",
                               component_scales = NULL) {
  trait_count <- ncol(fixture$phenotype)
  mtblr_bed(
    y = fixture$phenotype, Glist = fixture$Glist,
    trait_ids = colnames(fixture$phenotype),
    sample_ids = rownames(fixture$phenotype),
    method = method, component_scales = component_scales,
    initial_scale_probability = if (identical(method, "bayesr"))
      rep(1 / length(component_scales), length(component_scales)) else NULL,
    scale_dirichlet_prior = if (identical(method, "bayesr"))
      rep(1, length(component_scales)) else NULL,
    marker_multipliers = if (identical(method, "bayesr"))
      rep(1, length(fixture$Glist$rsids[[1L]])) else NULL,
    residual_covariance_policy = if (sampled) "sampled_full" else
      "fixed_full",
    fixed_residual_covariance = if (sampled) NULL else fixture$Ve,
    initial_marker_covariance = fixture$Vb,
    marker_covariance_prior_degrees_of_freedom = trait_count + 0.5,
    marker_covariance_prior_scale = fixture$Psi_b,
    initial_residual_covariance = if (sampled) fixture$Ve else NULL,
    residual_covariance_prior_degrees_of_freedom = if (sampled) {
      trait_count - 0.5
    } else NULL,
    residual_covariance_prior_scale = if (sampled) fixture$Psi_e else NULL,
    nburn = 1L, nit = 3L, nthin = 1L, seed = seed,
    nchains = chains, ncores = cores, chain_seeds = chain_seeds,
    keep_traces = TRUE,
    memory_limit_bytes = Inf)
}

phase5b_internal_raw <- function(fixture, sampled = FALSE, chains = 1L,
                                 cores = 1L, seed = 51,
                                 chain_seeds = NULL, method = "bayesc",
                                 component_scales = NULL) {
  trait_count <- ncol(fixture$phenotype)
  sblr:::.blr_cheng_mt_bayesc_bed_qualification(
    y = fixture$phenotype, Glist = fixture$Glist,
    fixed_residual_covariance = if (sampled) NULL else fixture$Ve,
    initial_marker_covariance = fixture$Vb,
    marker_covariance_prior_df = trait_count + 0.5,
    marker_covariance_prior_scale = fixture$Psi_b,
    method = method, component_scales = component_scales,
    initial_scale_probability = if (identical(method, "bayesr"))
      rep(1 / length(component_scales), length(component_scales)) else NULL,
    scale_dirichlet_prior = if (identical(method, "bayesr"))
      rep(1, length(component_scales)) else NULL,
    marker_multipliers = if (identical(method, "bayesr"))
      rep(1, length(fixture$Glist$rsids[[1L]])) else NULL,
    burn_in_iterations = 1L, sampling_iterations = 3L,
    thin_interval = 1L, chains = chains, cores = cores, seed = seed,
    chain_seeds = chain_seeds, keep_traces = TRUE,
    residual_covariance_policy = if (sampled) "sampled_full" else
      "fixed_full",
    initial_residual_covariance = if (sampled) fixture$Ve else NULL,
    residual_covariance_prior_df = if (sampled) trait_count - 0.5 else NULL,
    residual_covariance_prior_scale = if (sampled) fixture$Psi_e else NULL,
    memory_limit_bytes = Inf)
}

test_that("public BED pattern-by-scale BayesR is neutral to qualification", {
  fixture <- phase5b_public_fixture(2L)
  on.exit(blr_unified_cleanup(fixture), add = TRUE)
  scales <- c(.1, .5, 1)
  public <- attr(phase5b_public_fit(
    fixture, seed = 2751L, method = "bayesr", component_scales = scales),
    "blr_raw", exact = TRUE)
  internal <- phase5b_internal_raw(
    fixture, seed = 2751L, method = "bayesr", component_scales = scales)
  expect_identical(public$posterior, internal$posterior)
  expect_identical(public$draws, internal$draws)
  expect_identical(public$final, internal$final)
  expect_identical(public$input$mcmc$task_seeds,
                   internal$input$mcmc$task_seeds)
  expect_silent(sblr:::validate_blr_raw_v2(public))
})

test_that("public corrected BED route is neutral to the qualified sampler", {
  for (case in list(c(2L, FALSE), c(3L, TRUE))) {
    local({
      fixture <- phase5b_public_fixture(case[[1L]])
      on.exit(blr_unified_cleanup(fixture), add = TRUE)
      chain_seeds <- if (case[[1L]] == 2L) -1 else NULL
      fit <- phase5b_public_fit(
        fixture, sampled = case[[2L]], chain_seeds = chain_seeds)
      public <- attr(fit, "blr_raw", exact = TRUE)
      internal <- phase5b_internal_raw(
        fixture, sampled = case[[2L]], chain_seeds = chain_seeds)
      expect_identical(public$posterior, internal$posterior)
      expect_identical(public$draws, internal$draws)
      expect_identical(public$final, internal$final)
      expect_identical(public$derived, internal$derived)
      expect_identical(public$input$mcmc$task_seeds,
                       internal$input$mcmc$task_seeds)
      expect_identical(public$input$mcmc$retained_transition_indices,
                       internal$input$mcmc$retained_transition_indices)
      expect_identical(
        public$diagnostics$qualification$current_legacy_mt_route_used,
        FALSE)
    })
  }
})

test_that("public formatted output is a validated raw-v2 view", {
  fixture <- phase5b_public_fixture(2L)
  on.exit(blr_unified_cleanup(fixture), add = TRUE)
  fit <- phase5b_public_fit(fixture)
  raw <- attr(fit, "blr_raw", exact = TRUE)
  expect_silent(sblr:::validate_blr_raw_v2(raw))
  expect_s3_class(fit, "mtblr_fit")
  expect_s3_class(fit, "blr_fit")
  expect_identical(fit$bm, raw$posterior$realised_effect_mean)
  expect_identical(fit$dm, raw$posterior$pips)
  expect_identical(fit$activity_pattern_probabilities,
                   raw$posterior$activity_pattern_probabilities)
  expect_identical(fit$pleiotropic_probabilities,
                   raw$posterior$pleiotropic_probabilities)
  expect_identical(fit$marker_covariance_draws,
                   raw$draws$marker_covariance)
  expect_null(fit$residual_covariance_draws)
  expect_null(fit$cov_e_mean)
  expect_identical(fit$final_residual_covariance,
                   raw$final$residual_covariance)
  expect_identical(raw$diagnostics$qualification$status,
                   "publicly_supported")
})

test_that("four-trait public chains are serial-parallel identical", {
  fixture <- phase5b_public_fixture(4L)
  on.exit(blr_unified_cleanup(fixture), add = TRUE)
  serial <- attr(phase5b_public_fit(
    fixture, sampled = TRUE, chains = 2L, cores = 1L, seed = 79),
    "blr_raw", exact = TRUE)
  parallel <- attr(phase5b_public_fit(
    fixture, sampled = TRUE, chains = 2L, cores = 2L, seed = 79),
    "blr_raw", exact = TRUE)
  expect_identical(serial$posterior, parallel$posterior)
  expect_identical(serial$draws, parallel$draws)
  expect_identical(serial$final, parallel$final)
  expect_identical(serial$input$mcmc$task_seeds,
                   parallel$input$mcmc$task_seeds)
  if (isTRUE(parallel$diagnostics$workers$openmp_available) &&
      parallel$diagnostics$workers$actual_team_size > 1L) {
    expect_identical(parallel$diagnostics$workers$actual_team_size, 2L)
    expect_gt(length(unique(parallel$diagnostics$workers$task_worker_ids)), 1L)
  }
})

test_that("public summary MT operators expose only corrected Cheng interfaces", {
  expect_true("mtblr_csr" %in% getNamespaceExports("sblr"))
  expect_true("mtblr_block_eigen" %in% getNamespaceExports("sblr"))
  expect_identical(formals(mtblr_csr)$method, "bayesc")
  expect_identical(formals(mtblr_block_eigen)$method, "bayesc")
  expect_false(exists("mtblr_csr_chains_raw_internal",
                      envir = asNamespace("sblr"), inherits = FALSE))
  expect_false(exists("mtblr_block_eigen_chains_raw_internal",
                      envir = asNamespace("sblr"), inherits = FALSE))
})

test_that("public covariance inputs must follow declared trait order", {
  fixture <- phase5b_public_fixture(3L)
  on.exit(blr_unified_cleanup(fixture), add = TRUE)
  bad <- fixture$Vb[c(2L, 1L, 3L), c(2L, 1L, 3L)]
  expect_error(
    mtblr_bed(
      y = fixture$phenotype, Glist = fixture$Glist,
      trait_ids = colnames(fixture$phenotype),
      sample_ids = rownames(fixture$phenotype),
      fixed_residual_covariance = fixture$Ve,
      initial_marker_covariance = bad,
      marker_covariance_prior_degrees_of_freedom = 3.5,
      marker_covariance_prior_scale = fixture$Psi_b,
      nburn = 0L, nit = 1L),
    "initial_marker_covariance dimnames must follow the declared trait order")
})
