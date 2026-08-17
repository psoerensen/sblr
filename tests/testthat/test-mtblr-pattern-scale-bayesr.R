phase7_alleles <- function(marker_ids) {
  data.frame(
    marker_id = marker_ids,
    effect = rep(c("A", "C"), length.out = length(marker_ids)),
    other = rep(c("G", "T"), length.out = length(marker_ids)),
    coding = "effect_allele_count", stringsAsFactors = FALSE)
}

phase7_summary_fixture <- function(T = 3L, representation = "csr") {
  markers <- paste0("m", 1:4)
  traits <- paste0("trait", seq_len(T))
  alleles <- phase7_alleles(markers)
  map <- sblr:::.blr_new_global_marker_map(markers, alleles)
  cross_product <- matrix(c(
    4, .25, 0, 0,
    .25, 3, 0, 0,
    0, 0, 2.5, .125,
    0, 0, .125, 2
  ), 4L, 4L, byrow = TRUE, dimnames = list(markers, markers))
  resource <- switch(representation,
    csr = sblr:::.blr_csr_from_dense_resource(
      "resource", cross_product, markers, alleles),
    full = sblr:::.blr_block_eigen_resource(
      "resource", cross_product, markers, alleles,
      blocks = list(1:2, 3:4)),
    retained = sblr:::.blr_block_eigen_resource(
      "resource", cross_product, markers, alleles,
      blocks = list(1:2, 3:4), retained_ranks = c(1L, 1L)))
  providers <- lapply(seq_len(T), function(trait) {
    id <- paste0("provider", trait)
    sblr:::.blr_new_likelihood_provider(
      provider_id = id, trait_ids = traits[trait],
      operator_resource_id = resource$resource_id,
      local_to_global = stats::setNames(seq_along(markers), markers),
      sufficient_statistics = list(
        score = matrix(
          c(.5, -.25, .35, .15) + .0625 * (trait - 1L), ncol = 1L,
          dimnames = list(marker = markers, trait = traits[trait])),
        residual_scale = .75 + .25 * trait),
      sample_size = stats::setNames(80 + 10 * trait, traits[trait]),
      likelihood_regime = "independent_summary",
      residual_contract = "fixed_provider_residual_scale",
      population = paste0("population", trait),
      effect_scale = "standardized", overlap_group = NULL,
      provenance = list(source = "Phase 7 fixture"))
  })
  names(providers) <- vapply(providers, `[[`, character(1), "provider_id")
  list(
    collection = sblr:::.blr_new_provider_collection(
      map, list(resource = resource), providers,
      likelihood_regime = "independent_summary",
      analysis_mode = "joint_multitrait"),
    traits = traits, markers = markers, cross_product = cross_product,
    resource = resource)
}

phase7_summary_fit <- function(
    fixture, method = "bayesr", representation = NULL,
    scales = c(.1, 1), chains = 1L, cores = 1L, seed = 2701,
    update = TRUE) {
  T <- length(fixture$traits)
  covariance <- diag(seq(.35, .35 + .05 * (T - 1L), length.out = T))
  dimnames(covariance) <- list(fixture$traits, fixture$traits)
  sblr:::.blr_cheng_mt_bayesc_summary_qualification(
    collection = fixture$collection, trait_ids = fixture$traits,
    initial_marker_covariance = covariance,
    marker_covariance_prior_df = T + .5,
    marker_covariance_prior_scale = covariance,
    update_marker_covariance = update,
    update_activity_pattern_probability = update,
    burn_in_iterations = 3L, sampling_iterations = 6L,
    thin_interval = 2L, chains = chains, cores = cores, seed = seed,
    keep_traces = TRUE, memory_limit_bytes = Inf, method = method,
    component_scales = if (identical(method, "bayesr")) scales else NULL,
    initial_scale_probability = if (identical(method, "bayesr"))
      rep(1 / length(scales), length(scales)) else NULL,
    scale_dirichlet_prior = if (identical(method, "bayesr"))
      rep(1, length(scales)) else NULL,
    marker_multipliers = if (identical(method, "bayesr"))
      rep(1, length(fixture$markers)) else NULL)
}

phase7_common_scientific <- function(raw) {
  list(
    posterior = raw$posterior[c(
      "realised_effect_mean", "latent_effect_mean", "trait_pip",
      "activity_pattern_probabilities", "pleiotropic_probabilities",
      "joint_probability_parameter_mean", "marker_covariance_mean")],
    draws = raw$draws[c(
      "joint_states", "trait_activity", "latent_effects", "realised_effects",
      "joint_probability_parameters", "marker_covariance")],
    final = raw$final[c(
      "joint_states", "latent_effects", "realised_effects",
      "joint_probability_parameters", "marker_covariance")],
    retained = raw$input$mcmc$retained_transition_indices,
    convergence = raw$convergence,
    task_seeds = raw$input$mcmc$task_seeds)
}

test_that("Phase 7 pattern-by-scale conditional matches the independent reference", {
  reference <- new.env(parent = baseenv())
  sys.source(test_path("..", "research", "mtblr_covariance",
                       "mtblr_pattern_scale_reference.R"), reference)
  patterns <- as.matrix(expand.grid(rep(list(0:1), 3L)))
  patterns <- matrix(as.integer(t(patterns)), ncol = 3L, byrow = TRUE)
  h <- c(.35, -.22, .41)
  likelihood <- matrix(c(
    2.1, .3, -.15, .3, 1.7, .2, -.15, .2, 1.4
  ), 3L, 3L, byrow = TRUE)
  Vb <- matrix(c(
    .8, .18, -.1, .18, .65, .12, -.1, .12, .7
  ), 3L, 3L, byrow = TRUE)
  pi <- c(.3, rep(.1, 7L))
  scales <- c(.05, .4, 1.2)
  omega <- c(.25, .5, .25)
  q <- 1.35
  expected <- reference$mtblr_pattern_scale_conditional(
    h, likelihood, Vb, patterns, pi, scales, omega, q)
  actual <- sblr:::mtblr_phase7_pattern_scale_contract_internal(
    h, likelihood, Vb, pi, patterns, scales, omega, q)
  expect_equal(actual$probability, expected$probability, tolerance = 2e-13)
  expect_equal(actual$pattern_index + 1L, expected$states$pattern)
  expect_equal(actual$scale_index[-1L] + 1L, expected$states$scale[-1L])
  for (state in 2:length(actual$probability)) {
    root <- sqrt(scales[expected$states$scale[state]] * q)
    expect_equal(as.numeric(actual$base_active_mean[[state]]),
                 expected$active_mean[[state]] / root, tolerance = 2e-13)
    expect_equal(actual$base_active_covariance[[state]],
                 expected$active_covariance[[state]] / root^2,
                 tolerance = 2e-13)
  }
})

test_that("one positive unit scale reduces exactly to Cheng BayesC", {
  fixture <- phase7_summary_fixture(3L, "csr")
  bayesc <- phase7_summary_fit(fixture, method = "bayesc", seed = 2711)
  bayesr <- phase7_summary_fit(
    fixture, method = "bayesr", scales = 1, seed = 2711)
  expect_identical(phase7_common_scientific(bayesr),
                   phase7_common_scientific(bayesc))
  expect_true(all(bayesr$draws$component_assignments >= -1L))
  expect_identical(dim(
    bayesr$posterior$joint_component_assignment_probabilities),
                   c(length(fixture$markers), 1L))
})

test_that("scaled completed effects remove gamma q exactly once in Vb", {
  fixture <- phase7_summary_fixture(3L, "csr")
  scales <- c(.1, .7, 1.4)
  q <- c(.8, 1.1, .65, 1.25)
  T <- length(fixture$traits)
  covariance <- diag(c(.35, .4, .45))
  dimnames(covariance) <- list(fixture$traits, fixture$traits)
  raw <- sblr:::.blr_cheng_mt_bayesc_summary_qualification(
    fixture$collection, fixture$traits, covariance, T + .5, covariance,
    update_marker_covariance = FALSE,
    update_activity_pattern_probability = FALSE,
    burn_in_iterations = 2L, sampling_iterations = 4L,
    thin_interval = 2L, chains = 1L, cores = 1L, seed = 2721,
    keep_traces = TRUE, memory_limit_bytes = Inf, method = "bayesr",
    component_scales = scales,
    initial_scale_probability = rep(1 / 3, 3L),
    scale_dirichlet_prior = rep(1, 3L), marker_multipliers = q)
  theta <- raw$final$latent_effects[1L, , , drop = TRUE]
  component <- raw$final$component_assignments[1L, , drop = TRUE]
  active <- which(component >= 0L)
  expected <- matrix(0, T, T)
  for (marker in active) {
    expected <- expected + tcrossprod(theta[marker, ]) /
      (scales[component[marker] + 1L] * q[marker])
  }
  expect_equal(raw$diagnostics$qualification$covariance_updates[[1L]]$statistic,
               expected, tolerance = 2e-13)
})

test_that("BayesR shares CSR and block-eigen targets and scheduler behavior", {
  csr <- phase7_summary_fit(
    phase7_summary_fixture(3L, "csr"), seed = 2731, update = FALSE)
  full <- phase7_summary_fit(
    phase7_summary_fixture(3L, "full"), seed = 2731, update = FALSE)
  expect_equal(phase7_common_scientific(csr),
               phase7_common_scientific(full), tolerance = 2e-14)
  retained_fixture <- phase7_summary_fixture(3L, "retained")
  reconstructed <- sblr:::.blr_operator_matrix(retained_fixture$resource)
  csr_resource <- sblr:::.blr_csr_from_dense_resource(
    "resource", reconstructed, retained_fixture$markers,
    phase7_alleles(retained_fixture$markers))
  providers <- lapply(retained_fixture$collection$providers, function(x) {
    sblr:::.blr_new_likelihood_provider(
      x$provider_id, x$trait_ids, csr_resource$resource_id,
      x$local_to_global, x$sufficient_statistics, x$sample_size,
      x$likelihood_regime, x$residual_contract, x$population,
      x$effect_scale, x$overlap_group, x$provenance)
  })
  names(providers) <- vapply(providers, `[[`, character(1), "provider_id")
  csr_reconstruction <- retained_fixture
  csr_reconstruction$collection <- sblr:::.blr_new_provider_collection(
    retained_fixture$collection$global_marker_map,
    list(resource = csr_resource), providers,
    likelihood_regime = "independent_summary",
    analysis_mode = "joint_multitrait")
  retained <- phase7_summary_fit(retained_fixture, seed = 2732, update = FALSE)
  reconstructed_fit <- phase7_summary_fit(
    csr_reconstruction, seed = 2732, update = FALSE)
  expect_silent(sblr:::validate_blr_raw_v2(retained))
  expect_identical(retained$draws$joint_states,
                   reconstructed_fit$draws$joint_states)
  expect_identical(retained$draws$component_assignments,
                   reconstructed_fit$draws$component_assignments)
  expect_equal(phase7_common_scientific(retained),
               phase7_common_scientific(reconstructed_fit),
               tolerance = 1e-7)
  fixture <- phase7_summary_fixture(3L, "csr")
  serial <- phase7_summary_fit(fixture, chains = 2L, cores = 1L, seed = 2733)
  parallel <- phase7_summary_fit(fixture, chains = 2L, cores = 2L, seed = 2733)
  expect_identical(phase7_common_scientific(serial),
                   phase7_common_scientific(parallel))
})

test_that("BayesR BED and summary routes reduce under diagonal residual covariance", {
  fixture <- blr_unified_fixture()
  on.exit(blr_unified_cleanup(fixture), add = TRUE)
  markers <- fixture$Glist$rsids[[1L]]
  traits <- c("trait1", "trait2")
  phenotype <- cbind(
    trait1 = fixture$y[, 1L],
    trait2 = .65 * fixture$y[, 1L] +
      c(-.07, .03, .05, -.01, .04, -.06, .02, 0))
  rownames(phenotype) <- rownames(fixture$y)
  frequency <- rowMeans(fixture$dosage) / 2
  standardized <- t((fixture$dosage - 2 * frequency) /
    sqrt(2 * frequency * (1 - frequency)))
  cross_product <- crossprod(standardized)
  dimnames(cross_product) <- list(markers, markers)
  score <- crossprod(standardized, phenotype)
  alleles <- phase7_alleles(markers)
  map <- sblr:::.blr_new_global_marker_map(markers, alleles)
  resource <- sblr:::.blr_csr_from_dense_resource(
    "resource", cross_product, markers, alleles)
  residual_scale <- c(.8, 1.15)
  providers <- lapply(seq_along(traits), function(index) {
    sblr:::.blr_new_likelihood_provider(
      provider_id = paste0("provider", index), trait_ids = traits[index],
      operator_resource_id = "resource",
      local_to_global = stats::setNames(seq_along(markers), markers),
      sufficient_statistics = list(
        score = matrix(score[, index], ncol = 1L,
                       dimnames = list(marker = markers,
                                       trait = traits[index])),
        residual_scale = residual_scale[index]),
      sample_size = stats::setNames(nrow(phenotype), traits[index]),
      likelihood_regime = "independent_summary",
      residual_contract = "fixed_provider_residual_scale",
      population = paste0("population", index), effect_scale = "standardized",
      overlap_group = NULL, provenance = list(source = "Phase 7 reduction"))
  })
  names(providers) <- paste0("provider", seq_along(traits))
  collection <- sblr:::.blr_new_provider_collection(
    map, list(resource = resource), providers,
    likelihood_regime = "independent_summary",
    analysis_mode = "joint_multitrait")
  Vb <- diag(c(.4, .35)); dimnames(Vb) <- list(traits, traits)
  Ve <- diag(residual_scale); dimnames(Ve) <- list(traits, traits)
  common <- list(
    initial_marker_covariance = Vb, marker_covariance_prior_df = 2.5,
    marker_covariance_prior_scale = Vb, update_marker_covariance = FALSE,
    update_activity_pattern_probability = FALSE,
    burn_in_iterations = 2L, sampling_iterations = 5L,
    thin_interval = 1L, chains = 1L, cores = 1L, seed = 2734L,
    keep_traces = TRUE, memory_limit_bytes = Inf, method = "bayesr",
    component_scales = c(.1, 1), initial_scale_probability = c(.4, .6),
    scale_dirichlet_prior = c(1, 1), marker_multipliers = rep(1, 3L))
  bed <- do.call(sblr:::.blr_cheng_mt_bayesc_bed_qualification, c(list(
    y = phenotype, Glist = fixture$Glist,
    fixed_residual_covariance = Ve), common))
  summary <- do.call(sblr:::.blr_cheng_mt_bayesc_summary_qualification, c(list(
    collection = collection, trait_ids = traits), common))
  expect_identical(bed$draws$joint_states, summary$draws$joint_states)
  expect_identical(bed$draws$component_assignments,
                   summary$draws$component_assignments)
  expect_equal(bed$draws$realised_effects, summary$draws$realised_effects,
               tolerance = 2e-9)
  expect_equal(bed$posterior$activity_pattern_probabilities,
               summary$posterior$activity_pattern_probabilities,
               tolerance = 2e-9)
})

test_that("BayesR raw v2 uses marginal scale output without joint expansion", {
  raw <- phase7_summary_fit(phase7_summary_fixture(3L, "csr"), seed = 2741)
  expect_silent(sblr:::validate_blr_raw_v2(raw))
  expect_identical(dim(raw$draws$component_assignments), c(3L, 1L, 4L))
  expect_identical(dim(raw$draws$joint_component_probability_parameters),
                   c(3L, 1L, 2L))
  expect_identical(
    dim(raw$posterior$joint_component_assignment_probabilities), c(4L, 2L))
  expect_false(any(vapply(raw$draws, function(x) {
    is.array(x) && length(dim(x)) >= 2L &&
      all(c(8L, 2L) %in% dim(x))
  }, logical(1))))
  expect_null(raw$draws$residual_covariance)
  expect_null(raw$posterior$predictions)
})

test_that("Phase 7 raw v2 enforces component sub-probability identity", {
  raw <- phase7_summary_fit(phase7_summary_fixture(3L, "csr"), seed = 2742)
  expect_silent(sblr:::validate_blr_raw_v2(raw))
  patterns <- raw$input$prior$probability$activity_patterns
  null_id <- rownames(patterns)[which(rowSums(patterns) == 0L)]
  nonnull <- 1 - raw$posterior$activity_pattern_probabilities[, null_id]
  marker <- which(nonnull > 0)[[1L]]
  invalid <- raw
  invalid$posterior$joint_component_assignment_probabilities[marker, ] <- 0
  expect_error(
    sblr:::validate_blr_raw_v2(invalid),
    "sub-probabilities must sum to one minus the declared null-pattern")
})

test_that("Phase 7 compact pattern and scale diagnostics are contracted", {
  raw <- phase7_summary_fit(phase7_summary_fixture(3L, "csr"), seed = 2743)
  expect_silent(sblr:::validate_blr_raw_v2(raw))

  missing_pattern <- raw
  missing_pattern$diagnostics$qualification["pattern_occupancy_counts"] <- NULL
  expect_error(
    sblr:::validate_blr_raw_v2(missing_pattern),
    "missing required compact transition fields")

  missing_scale <- raw
  missing_scale$diagnostics$qualification["scale_occupancy_counts"] <- NULL
  expect_error(
    sblr:::validate_blr_raw_v2(missing_scale),
    "missing required compact transition fields")

  inconsistent <- raw
  inconsistent$diagnostics$qualification$scale_occupancy_counts[[1L]][1L] <-
    inconsistent$diagnostics$qualification$scale_occupancy_counts[[1L]][1L] + 1
  expect_error(
    sblr:::validate_blr_raw_v2(inconsistent),
    "Scale occupancy totals must equal the non-null marker-update events")
})

test_that("Phase 7 memory includes conditional state workspaces", {
  patterns <- sblr:::.blr_phase4a_patterns(paste0("trait", 1:3))
  base <- sblr:::.blr_phase5a_memory_estimate(
    marker_count = 1000L, trait_count = 3L, observation_count = 0L,
    chains = 2L, retained_draws = 5L, convergence_count = 10L,
    keep_traces = TRUE, sampled_residual = FALSE,
    memory_limit_bytes = Inf, enforce = FALSE)
  adjusted <- sblr:::.blr_phase7_memory_adjust(
    base, 1000L, 3L, 4L, 2L, 5L, 10L, TRUE, Inf, enforce = FALSE,
    activity_patterns = patterns, concurrent_chains = 2L)
  expect_gt(adjusted$estimated_peak_incremental_bytes,
            base$estimated_peak_incremental_bytes)
  expect_identical(adjusted$components[["phase7_marker_component_summary"]],
                   8 * 1000 * 4)
  active <- rowSums(patterns)[-1L]
  states <- (nrow(patterns) - 1L) * 4L + 1L
  expect_equal(
    adjusted$components[["phase7_pattern_scale_state_tables"]],
    2 * states * 28, tolerance = 0)
  expect_equal(
    adjusted$components[["phase7_pattern_scale_active_numeric"]],
    2 * 4 * sum(active + active^2) * 8, tolerance = 0)
  expect_equal(
    adjusted$dimensions$pattern_scale_nonnull_states,
    (nrow(patterns) - 1L) * 4L, tolerance = 0)
})

test_that("Phase 7 rejects infeasible conditional workspaces before entry", {
  patterns <- sblr:::.blr_phase4a_patterns(paste0("trait", 1:12))
  base <- sblr:::.blr_phase5a_memory_estimate(
    marker_count = 1L, trait_count = 12L, observation_count = 0L,
    chains = 1L, retained_draws = 1L, convergence_count = 1L,
    keep_traces = TRUE, sampled_residual = FALSE,
    memory_limit_bytes = Inf, enforce = FALSE)
  large <- sblr:::.blr_phase7_memory_adjust(
    base, 1L, 12L, 1000L, 1L, 1L, 1L, TRUE, 256 * 1024^2,
    enforce = FALSE, activity_patterns = patterns, concurrent_chains = 1L)
  expect_true(large$limit_exceeded)
  expect_gt(
    large$components[["phase7_pattern_scale_active_numeric"]],
    256 * 1024^2)
  expect_error(
    sblr:::.blr_phase7_memory_adjust(
      base, 1L, 12L, 1000L, 1L, 1L, 1L, TRUE, 256 * 1024^2,
      enforce = TRUE, activity_patterns = patterns, concurrent_chains = 1L),
    "before provider construction and sampling.*T = 12, K = 1000.*state count = 4095001")

  fixture <- blr_unified_fixture()
  on.exit(blr_unified_cleanup(fixture), add = TRUE)
  fixture$Glist$rsids[[1L]] <- fixture$Glist$rsids[[1L]][1L]
  fixture$Glist$rsidsLD[[1L]] <- fixture$Glist$rsids[[1L]]
  fixture$Glist$af[[1L]] <- fixture$Glist$af[[1L]][1L]
  phenotype <- matrix(
    rep(fixture$y[, 1L], 12L), ncol = 12L,
    dimnames = list(rownames(fixture$y), paste0("trait", 1:12)))
  provider_entries <- 0L
  native_entries <- 0L
  testthat::with_mocked_bindings(
    expect_error(
      sblr:::.blr_cheng_mt_bayesc_bed_qualification(
        y = phenotype, Glist = fixture$Glist, chr = 1L, cls = list(1L),
        fixed_residual_covariance = diag(12),
        initial_marker_covariance = diag(12),
        marker_covariance_prior_df = 12,
        marker_covariance_prior_scale = diag(12),
        update_marker_covariance = FALSE,
        update_activity_pattern_probability = FALSE,
        burn_in_iterations = 0L, sampling_iterations = 1L,
        memory_limit_bytes = 256 * 1024^2, method = "bayesr",
        component_scales = seq_len(1000L) / 1000,
        initial_scale_probability = rep(1 / 1000, 1000L),
        scale_dirichlet_prior = rep(1, 1000L), marker_multipliers = 1),
      "before provider construction and sampling"),
    .blr_phase4a_bed_collection = function(...) {
      provider_entries <<- provider_entries + 1L
      stop("PROVIDER_ENTERED", call. = FALSE)
    },
    mtblr_phase4a_cheng_bed_internal = function(...) {
      native_entries <<- native_entries + 1L
      stop("NATIVE_ENTERED", call. = FALSE)
    },
    .package = "sblr")
  expect_identical(provider_entries, 0L)
  expect_identical(native_entries, 0L)
})

test_that("Phase 7 native workspace guard matches R dimensions", {
  fixture <- phase7_summary_fixture(3L, "csr")
  patterns <- sblr:::.blr_phase4a_patterns(fixture$traits)
  base <- sblr:::.blr_phase6a_memory_estimate(
    fixture$collection, 3L, 1L, 3L, 6L, TRUE, Inf, enforce = FALSE)
  adjusted <- sblr:::.blr_phase7_memory_adjust(
    base, 4L, 3L, 2L, 1L, 3L, 6L, TRUE, Inf, enforce = FALSE,
    activity_patterns = patterns, concurrent_chains = 1L)
  workspace_fields <- c(
    "phase7_pattern_scale_state_tables",
    "phase7_pattern_scale_active_containers",
    "phase7_pattern_scale_active_numeric",
    "phase7_pattern_scale_kernel_containers")
  expected_workspace <- sum(adjusted$components[workspace_fields])
  original_adjust <- sblr:::.blr_phase7_memory_adjust
  testthat::with_mocked_bindings(
    expect_error(
      phase7_summary_fit(fixture, scales = c(.1, 1), seed = 2744),
      paste0("pattern_scale_conditional_workspace=", expected_workspace)),
    .blr_phase7_memory_adjust = function(
        memory, marker_count, trait_count, scale_count, chains,
        retained_draws, convergence_count, keep_traces, memory_limit_bytes,
        enforce = TRUE, activity_patterns = NULL,
        concurrent_chains = chains, workspace_enabled = scale_count > 0L) {
      value <- original_adjust(
        memory, marker_count, trait_count, scale_count, chains,
        retained_draws, convergence_count, keep_traces, Inf,
        enforce = FALSE, activity_patterns = activity_patterns,
        concurrent_chains = concurrent_chains,
        workspace_enabled = workspace_enabled)
      value$limit_bytes <- 1
      value$limit_exceeded <- FALSE
      value
    },
    .package = "sblr")
})
