general_t_fixture <- function(trait_count) {
  fixture <- blr_unified_fixture()
  base <- fixture$y[, 1L]
  offsets <- seq(-0.20, 0.15, length.out = length(base))
  phenotype <- vapply(seq_len(trait_count), function(trait) {
    (0.55 - 0.12 * trait) * base +
      ((-1)^trait) * (0.04 + 0.015 * trait) * offsets +
      0.01 * trait * seq_along(base)
  }, numeric(length(base)))
  colnames(phenotype) <- paste0("trait", seq_len(trait_count))
  rownames(phenotype) <- rownames(fixture$y)
  marker_covariance <- outer(seq_len(trait_count), seq_len(trait_count),
                             function(i, j) 0.06^(abs(i - j) + 1L))
  diag(marker_covariance) <- seq(0.45, 0.30, length.out = trait_count)
  residual_covariance <- outer(seq_len(trait_count), seq_len(trait_count),
                               function(i, j) 0.10 / (1 + abs(i - j)))
  diag(residual_covariance) <- seq(1.05, 0.75, length.out = trait_count)
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

general_t_fit <- function(fixture, sampled_residual = FALSE,
                          chains = 1L, cores = 1L, seed = 37,
                          sampling_iterations = 3L, ...) {
  trait_count <- ncol(fixture$phenotype)
  sblr:::.blr_cheng_mt_bayesc_bed_qualification(
    y = fixture$phenotype, Glist = fixture$Glist,
    fixed_residual_covariance = if (sampled_residual) NULL else fixture$Ve,
    initial_marker_covariance = fixture$Vb,
    marker_covariance_prior_df = trait_count + 0.5,
    marker_covariance_prior_scale = fixture$Psi_b,
    burn_in_iterations = 1L, sampling_iterations = sampling_iterations,
    thin_interval = 1L, chains = chains, cores = cores, seed = seed,
    residual_covariance_policy = if (sampled_residual) {
      "sampled_full"
    } else "fixed_full",
    initial_residual_covariance = if (sampled_residual) fixture$Ve else NULL,
    residual_covariance_prior_df = if (sampled_residual) {
      trait_count - 0.5
    } else NULL,
    residual_covariance_prior_scale = if (sampled_residual) {
      fixture$Psi_e
    } else NULL, ...)
}

test_that("general Cheng patterns are canonical, complete, and bounded", {
  expected_two <- rbind(
    `0_0` = c(0L, 0L), `1_0` = c(1L, 0L),
    `0_1` = c(0L, 1L), `1_1` = c(1L, 1L))
  colnames(expected_two) <- c("a", "b")
  expect_identical(sblr:::.blr_phase4a_patterns(c("a", "b")), expected_two)

  for (trait_count in 2:4) {
    patterns <- sblr:::.blr_phase4a_patterns(
      paste0("trait", seq_len(trait_count)))
    expect_identical(nrow(patterns), as.integer(2^trait_count))
    expect_identical(nrow(unique(patterns)), nrow(patterns))
    expect_identical(sum(rowSums(patterns) == 0L), 1L)
    expect_identical(sum(rowSums(patterns) == trait_count), 1L)
  }
  expect_error(sblr:::.blr_phase4a_patterns(paste0("t", 1:13)),
               "T in \\[2, 12\\].*8192 patterns")
})

test_that("Phase 5A memory preflight accounts for mandatory exponential output", {
  estimate <- function(marker_count, trait_count, chains = 1,
                       retained_draws = 1, convergence_count = 1) {
    sblr:::.blr_phase5a_memory_estimate(
      marker_count = marker_count, trait_count = trait_count,
      observation_count = 1000, chains = chains,
      retained_draws = retained_draws,
      convergence_count = convergence_count,
      sampled_residual = TRUE, keep_traces = TRUE,
      memory_limit_bytes = Inf)
  }
  cases <- list(
    A = list(M = 100000, T = 6, bytes = 51200000),
    B = list(M = 500000, T = 8, bytes = 1024000000),
    C = list(M = 500000, T = 12, bytes = 16384000000))
  for (case in cases) {
    observed <- estimate(case$M, case$T)
    expect_equal(
      unname(observed$components[
        "marker_pattern_activity_probability_matrix"]),
      case$bytes, tolerance = 0)
    expect_equal(
      unname(observed$components[
        "marker_pattern_joint_state_probability_matrix"]),
      case$bytes, tolerance = 0)
    expect_gt(observed$estimated_peak_incremental_bytes, case$bytes)
    expect_identical(observed$components[[
      "dense_transition_diagnostics"]], 0)
    expect_identical(observed$dense_transition_diagnostics, FALSE)
  }

  one_chain <- estimate(1000, 6, chains = 1, retained_draws = 2)
  two_chains <- estimate(1000, 6, chains = 2, retained_draws = 2)
  more_draws <- estimate(1000, 6, chains = 1, retained_draws = 4)
  expect_gt(two_chains$estimated_peak_incremental_bytes,
            one_chain$estimated_peak_incremental_bytes)
  expect_gt(two_chains$components[["native_retained_output"]],
            one_chain$components[["native_retained_output"]])
  expect_gt(more_draws$components[["native_retained_output"]],
            one_chain$components[["native_retained_output"]])
  expect_equal(
    one_chain$components[["compact_transition_diagnostics"]],
    (64 * 8) + 8, tolerance = 0)

  expect_error(
    sblr:::.blr_phase5a_memory_estimate(
      marker_count = 500000, trait_count = 8,
      observation_count = 1000, chains = 1, retained_draws = 1,
      convergence_count = 1, sampled_residual = FALSE,
      keep_traces = FALSE, memory_limit_bytes = 256 * 1024^2),
    "failed before sampling.*M = 500000.*T = 8.*K = 256")
  expect_error(
    sblr:::.blr_phase5a_memory_estimate(
      marker_count = 500000, trait_count = 12,
      observation_count = 1000, chains = 1, retained_draws = 1,
      convergence_count = 1, sampled_residual = FALSE,
      keep_traces = FALSE, memory_limit_bytes = 256 * 1024^2),
    "failed before sampling.*K = 4096")
  expect_error(
    sblr:::.blr_phase5a_memory_estimate(
      marker_count = 2^41 + 1, trait_count = 12,
      observation_count = 1, chains = 1, retained_draws = 1,
      convergence_count = 1, sampled_residual = FALSE,
      keep_traces = FALSE, memory_limit_bytes = Inf),
               "overflowed component 'marker_pattern_count'.*before sampling")
})

test_that("packed-BED allocation accounting matches native aligned strides", {
  sample_counts <- c(1, 2, 3, 4, 5, 255, 256, 257, 511, 512, 513,
                     200000)
  expected_stride <- c(rep(64, 7), 128, 128, 128, 192, 50048)
  for (index in seq_along(sample_counts)) {
    selected <- sample_counts[[index]]
    source <- selected + 11
    r_value <- sblr:::.blr_phase5a_packed_bed_allocation(
      marker_count = 3, selected_sample_count = selected,
      source_sample_count = source, selected_rows_used = TRUE)
    native <- sblr:::mtblr_phase5a_packed_bed_allocation_internal(
      selected, 3, source, TRUE)
    expect_identical(r_value, native)
    expect_equal(r_value$aligned_stride_bytes,
                 expected_stride[[index]], tolerance = 0)
    expect_equal(r_value$owner_bytes,
                 3 * expected_stride[[index]], tolerance = 0)
    expect_equal(r_value$source_row_buffer_bytes,
                 ceiling(source / 4), tolerance = 0)
  }

  all_samples <- sblr:::.blr_phase5a_packed_bed_allocation(
    marker_count = 7, selected_sample_count = 513,
    source_sample_count = 513, selected_rows_used = FALSE)
  expect_equal(all_samples$owner_bytes, 7 * 192, tolerance = 0)
  expect_equal(all_samples$source_row_buffer_bytes, 0, tolerance = 0)
  expect_identical(
    all_samples,
    sblr:::mtblr_phase5a_packed_bed_allocation_internal(
      513, 7, 513, FALSE))

  one_chain <- sblr:::.blr_phase5a_memory_estimate(
    marker_count = 1000, trait_count = 3, observation_count = 257,
    source_sample_count = 513, selected_rows_used = TRUE,
    chains = 1, retained_draws = 1, convergence_count = 1,
    sampled_residual = FALSE, keep_traces = FALSE,
    memory_limit_bytes = Inf)
  two_chains <- sblr:::.blr_phase5a_memory_estimate(
    marker_count = 1000, trait_count = 3, observation_count = 257,
    source_sample_count = 513, selected_rows_used = TRUE,
    chains = 2, retained_draws = 1, convergence_count = 1,
    sampled_residual = FALSE, keep_traces = FALSE,
    memory_limit_bytes = Inf)
  expect_equal(one_chain$components[["packed_bed_owner"]], 128000,
               tolerance = 0)
  expect_equal(one_chain$components[["packed_bed_source_row_buffer"]], 129,
               tolerance = 0)
  expect_equal(two_chains$components[["packed_bed_owner"]],
               one_chain$components[["packed_bed_owner"]], tolerance = 0)
  expect_error(
    sblr:::.blr_phase5a_packed_bed_allocation(
      marker_count = 2^53, selected_sample_count = 2^53,
      source_sample_count = 2^53, selected_rows_used = FALSE),
    "overflowed component 'packed_bed_owner'")
  expect_error(
    sblr:::mtblr_phase5a_packed_bed_allocation_internal(
      2^53, 2^53, 2^53, FALSE),
    "packed_bed_owner")
})

test_that("packed-BED counterexample fails before provider and native entry", {
  allocation <- sblr:::.blr_phase5a_memory_estimate(
    marker_count = 10000, trait_count = 2,
    observation_count = 200000, source_sample_count = 200000,
    selected_rows_used = FALSE, chains = 1, retained_draws = 1,
    convergence_count = 1, sampled_residual = FALSE,
    keep_traces = FALSE, memory_limit_bytes = Inf)
  expect_equal(allocation$components[["packed_bed_owner"]],
               500480000, tolerance = 0)
  expect_gt(allocation$estimated_peak_incremental_bytes, 256 * 1024^2)

  fixture <- blr_unified_fixture()
  on.exit(blr_unified_cleanup(fixture), add = TRUE)
  marker_count <- 10000L
  source_count <- 200000L
  fixture$Glist$n <- source_count
  fixture$Glist$ids <- NULL
  fixture$Glist$rsids[[1L]] <- paste0("m", seq_len(marker_count))
  fixture$Glist$rsidsLD[[1L]] <- fixture$Glist$rsids[[1L]]
  fixture$Glist$af[[1L]] <- rep(0.25, marker_count)
  phenotype <- matrix(0, source_count, 2L,
                      dimnames = list(NULL, c("trait1", "trait2")))
  provider_entries <- 0L
  native_entries <- 0L
  testthat::with_mocked_bindings(
    expect_error(
      sblr:::.blr_cheng_mt_bayesc_bed_qualification(
        y = phenotype, Glist = fixture$Glist, chr = 1L,
        cls = list(seq_len(marker_count)),
        fixed_residual_covariance = diag(2),
        initial_marker_covariance = diag(2),
        marker_covariance_prior_df = 2,
        marker_covariance_prior_scale = diag(2),
        burn_in_iterations = 0L, sampling_iterations = 1L,
        update_marker_covariance = FALSE,
        update_activity_pattern_probability = FALSE),
      "memory preflight failed before sampling.*packed_bed_owner"),
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

test_that("Phase 5A memory rejection occurs before provider and native entry", {
  fixture <- general_t_fixture(3L)
  on.exit(blr_unified_cleanup(fixture), add = TRUE)
  provider_entries <- 0L
  native_entries <- 0L
  testthat::with_mocked_bindings(
    expect_error(
      general_t_fit(fixture, memory_limit_bytes = 1),
      "memory preflight failed before sampling"),
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

test_that("native packed-BED guard runs before owner construction", {
  fixture <- general_t_fixture(3L)
  on.exit(blr_unified_cleanup(fixture), add = TRUE)
  estimator <- sblr:::.blr_phase5a_memory_estimate
  testthat::with_mocked_bindings(
    expect_error(
      general_t_fit(fixture, sampling_iterations = 1L,
                    memory_limit_bytes = 1),
      "Packed BED allocation guard failed before construction.*packed_bed_owner"),
    .blr_phase5a_memory_estimate = function(...) {
      arguments <- list(...)
      arguments$memory_limit_bytes <- Inf
      arguments$enforce <- FALSE
      value <- do.call(estimator, arguments)
      value$limit_bytes <- 1
      value
    },
    .package = "sblr")
})

test_that("memory preflight and compact validation are zero-RNG", {
  for (trait_count in 2:3) {
    fixture <- general_t_fixture(trait_count)
    on.exit(blr_unified_cleanup(fixture), add = TRUE)
    for (sampled in c(FALSE, TRUE)) {
      finite <- general_t_fit(
        fixture, sampled_residual = sampled,
        sampling_iterations = 2L, memory_limit_bytes = 256 * 1024^2)
      unlimited <- general_t_fit(
        fixture, sampled_residual = sampled,
        sampling_iterations = 2L, memory_limit_bytes = Inf)
      expect_identical(finite$posterior, unlimited$posterior)
      expect_identical(finite$draws, unlimited$draws)
      expect_identical(finite$final, unlimited$final)
      expect_identical(finite$derived, unlimited$derived)
      expect_identical(finite$input$mcmc$task_seeds,
                       unlimited$input$mcmc$task_seeds)
      expect_identical(
        finite$diagnostics$qualification$pattern_occupancy_counts,
        unlimited$diagnostics$qualification$pattern_occupancy_counts)
      expect_identical(
        finite$diagnostics$qualification$pattern_change_counts,
        unlimited$diagnostics$qualification$pattern_change_counts)
    }
  }
})

test_that("T=3 marker conditionals match an independent exact reference", {
  reference <- new.env(parent = globalenv())
  sys.source(test_path("..", "research", "mtblr_covariance",
                       "mtblr_general_t_reference.R"), envir = reference)
  trait_ids <- c("a", "b", "c")
  patterns <- reference$mtblr_general_t_patterns(trait_ids)
  score <- c(0.48, -0.31, 0.22)
  Vb <- matrix(c(.55, .14, -.06, .14, .42, .09,
                 -.06, .09, .38), 3L)
  Ve <- matrix(c(1.20, .28, -.12, .28, .95, .19,
                 -.12, .19, .82), 3L)
  probability <- c(.31, .08, .09, .10, .07, .11, .08, .16)
  expected <- reference$mtblr_general_t_marker_reference(
    score, 1.85, Vb, Ve, probability, patterns)
  observed <- sblr:::mtblr_phase4a_pattern_contract_internal(
    score, 1.85, Vb, Ve, probability, patterns)

  expect_identical(observed$pattern_order, rownames(patterns))
  expect_equal(observed$log_weight, expected$log_weight, tolerance = 1e-12)
  expect_equal(observed$probability, expected$probability, tolerance = 1e-13)
  for (state in 2:nrow(patterns)) {
    expect_equal(as.numeric(observed$active_mean[[state]]),
                 expected$active_mean[[state]], tolerance = 1e-12)
    expect_equal(observed$active_covariance[[state]],
                 expected$active_covariance[[state]], tolerance = 1e-12)
  }
  for (id in c("1_0_0", "1_1_0", "1_0_1")) {
    state <- match(id, rownames(patterns))
    active <- which(patterns[state, ] == 1L)
    active_value <- seq_along(active) / 10
    completion <- reference$mtblr_general_t_completion_reference(
      Vb, patterns[state, ], active_value)
    expect_equal(observed$completion_mean_coefficient[[state]],
                 completion$coefficient, tolerance = 1e-12)
    expect_equal(observed$completion_covariance[[state]],
                 completion$covariance, tolerance = 1e-12)
  }
  expect_null(observed$active_mean[[1L]])
})

test_that("T=3 fixed and sampled chains obey dynamic raw and covariance contracts", {
  reference <- new.env(parent = globalenv())
  sys.source(test_path("..", "research", "mtblr_covariance",
                       "mtblr_general_t_reference.R"), envir = reference)
  fixture <- general_t_fixture(3L)
  on.exit(blr_unified_cleanup(fixture), add = TRUE)
  fixed <- general_t_fit(fixture, sampled_residual = FALSE)
  sampled <- general_t_fit(fixture, sampled_residual = TRUE)

  expect_true(sblr:::validate_blr_raw_v2(fixed))
  expect_true(sblr:::validate_blr_raw_v2(sampled))
  expect_identical(dim(fixed$draws$realised_effects), c(3L, 1L, 3L, 3L))
  expect_identical(dim(fixed$draws$activity_pattern_parameters),
                   c(3L, 1L, 8L))
  expect_identical(dim(sampled$draws$residual_covariance),
                   c(3L, 1L, 3L, 3L))
  expect_null(fixed$draws$residual_covariance)
  expect_null(fixed$diagnostics$qualification$transition_counts)
  expect_identical(
    names(fixed$diagnostics$qualification$pattern_occupancy_counts[[1L]]),
    rownames(fixed$input$prior$probability$activity_patterns))
  expect_equal(
    sum(fixed$diagnostics$qualification$pattern_occupancy_counts[[1L]]),
    3 * (1 + 3), tolerance = 0)
  expect_true(is.numeric(
    fixed$diagnostics$qualification$pattern_change_counts))
  expect_identical(
    fixed$diagnostics$memory$contract,
    "phase5a_peak_incremental_fit_allocation_v1")
  expect_equal(
    fixed$diagnostics$memory$estimated_peak_incremental_bytes,
    fixed$input$output$memory_estimate_bytes, tolerance = 0)

  patterns <- fixed$input$prior$probability$activity_patterns
  expect_equal(unname(fixed$posterior$pips),
               unname(fixed$posterior$activity_pattern_probabilities %*%
                        patterns), tolerance = 1e-15, ignore_attr = TRUE)
  expect_equal(as.numeric(fixed$posterior$pleiotropic_probabilities),
               fixed$posterior$activity_pattern_probabilities[, "1_1_1"],
               tolerance = 1e-15, ignore_attr = TRUE)

  frequency <- rowMeans(fixture$dosage) / 2
  X <- t((fixture$dosage - 2 * frequency) /
           sqrt(2 * frequency * (1 - frequency)))
  last <- dim(sampled$draws$realised_effects)[1L]
  realised <- sampled$draws$realised_effects[last, 1L, , , drop = TRUE]
  latent <- sampled$draws$latent_effects[last, 1L, , , drop = TRUE]
  states <- sampled$draws$joint_states[last, 1L, , drop = TRUE]
  conditional <- reference$mtblr_general_t_covariance_conditionals(
    fixture$phenotype, X, realised, latent, states,
    3.5, fixture$Psi_b, 2.5, fixture$Psi_e)
  diagnostics <- sampled$diagnostics$qualification
  expect_equal(diagnostics$covariance_updates[[1L]]$statistic,
               conditional$marker_scale - fixture$Psi_b,
               tolerance = 1e-12, ignore_attr = TRUE)
  expect_equal(diagnostics$residual_covariance_updates[[1L]]$statistic,
               crossprod(conditional$residual), tolerance = 1e-12,
               ignore_attr = TRUE)
  expect_equal(diagnostics$covariance_updates[[1L]]$degrees_of_freedom,
               conditional$marker_degrees_of_freedom)
  expect_equal(
    diagnostics$residual_covariance_updates[[1L]]$degrees_of_freedom,
    conditional$residual_degrees_of_freedom)

  controlled <- general_t_fit(
    fixture, update_marker_covariance = FALSE,
    update_activity_pattern_probability = FALSE)
  for (draw in seq_len(dim(controlled$draws$marker_covariance)[1L])) {
    expect_equal(controlled$draws$marker_covariance[draw, 1L, , ],
                 fixture$Vb, tolerance = 0, ignore_attr = TRUE)
  }
  expect_equal(
    controlled$draws$activity_pattern_parameters,
    array(rep(1 / 8, 3L * 8L), dim = c(3L, 1L, 8L)),
    tolerance = 0, ignore_attr = TRUE)
})

test_that("T=4 runs both residual policies with complete dynamic axes", {
  fixture <- general_t_fixture(4L)
  on.exit(blr_unified_cleanup(fixture), add = TRUE)
  fixed <- general_t_fit(fixture, sampled_residual = FALSE,
                         sampling_iterations = 2L)
  sampled <- general_t_fit(fixture, sampled_residual = TRUE,
                           sampling_iterations = 2L)
  patterns <- fixed$input$prior$probability$activity_patterns
  expect_identical(dim(patterns), c(16L, 4L))
  expect_identical(rownames(patterns)[c(1L, 16L)],
                   c("0_0_0_0", "1_1_1_1"))
  expect_identical(
    names(fixed$input$prior$probability$activity_pattern_dirichlet),
    rownames(patterns))
  expect_identical(dim(fixed$draws$traitwise_activity), c(2L, 1L, 3L, 4L))
  expect_identical(dim(sampled$draws$marker_covariance), c(2L, 1L, 4L, 4L))
  expect_identical(dim(sampled$draws$residual_covariance), c(2L, 1L, 4L, 4L))
  expect_equal(unname(fixed$posterior$pips),
               unname(fixed$posterior$activity_pattern_probabilities %*%
                        patterns), tolerance = 1e-15, ignore_attr = TRUE)
  expect_true(sblr:::validate_blr_raw_v2(fixed))
  expect_true(sblr:::validate_blr_raw_v2(sampled))
  for (draw in seq_len(dim(sampled$draws$residual_covariance)[1L])) {
    expect_silent(chol(sampled$draws$residual_covariance[draw, 1L, , ]))
  }

  serial <- general_t_fit(fixture, sampled_residual = TRUE,
                          chains = 2L, cores = 1L,
                          sampling_iterations = 2L)
  parallel <- general_t_fit(fixture, sampled_residual = TRUE,
                            chains = 2L, cores = 2L,
                            sampling_iterations = 2L)
  expect_identical(serial$posterior, parallel$posterior)
  expect_identical(serial$draws, parallel$draws)
  expect_identical(serial$final, parallel$final)
  expect_identical(serial$derived, parallel$derived)
  expect_identical(serial$input$mcmc$task_seeds,
                   parallel$input$mcmc$task_seeds)
  capability <- sblr:::sparseLD_thread_info(2L)
  if (isTRUE(capability$openmp) &&
      capability$actual_threads_requested_region > 1L) {
    expect_gt(parallel$diagnostics$workers$actual_team_size, 1L)
    expect_gt(length(unique(parallel$diagnostics$workers$task_worker_ids)),
              1L)
  }
})

test_that("general Cheng conditionals are equivariant to declared trait order", {
  reference <- new.env(parent = globalenv())
  sys.source(test_path("..", "research", "mtblr_covariance",
                       "mtblr_general_t_reference.R"), envir = reference)
  traits <- c("a", "b", "c")
  patterns <- reference$mtblr_general_t_patterns(traits)
  score <- c(.4, -.2, .3)
  Vb <- matrix(c(.5,.1,.04,.1,.4,.07,.04,.07,.35), 3L)
  Ve <- matrix(c(1,.2,.1,.2,.9,.12,.1,.12,.8), 3L)
  probability <- seq_len(8L); probability <- probability / sum(probability)
  original <- sblr:::mtblr_phase4a_pattern_contract_internal(
    score, 2.1, Vb, Ve, probability, patterns)

  permutation <- c(3L, 1L, 2L)
  permuted_patterns <- reference$mtblr_general_t_patterns(traits[permutation])
  old_row_for_new <- apply(permuted_patterns, 1L, function(bits) {
    original_bits <- integer(3L)
    original_bits[permutation] <- bits
    match(paste(original_bits, collapse = "_"), rownames(patterns))
  })
  permuted <- sblr:::mtblr_phase4a_pattern_contract_internal(
    score[permutation], 2.1, Vb[permutation, permutation],
    Ve[permutation, permutation], probability[old_row_for_new],
    permuted_patterns)
  restored <- numeric(8L)
  restored[old_row_for_new] <- permuted$probability
  expect_equal(restored, original$probability, tolerance = 1e-12)
  for (new_state in 2:nrow(permuted_patterns)) {
    old_state <- old_row_for_new[[new_state]]
    old_active <- which(patterns[old_state, ] == 1L)
    new_active <- which(permuted_patterns[new_state, ] == 1L)
    restored_mean <- rep(NA_real_, 3L)
    restored_mean[permutation[new_active]] <-
      permuted$active_mean[[new_state]]
    expected_mean <- rep(NA_real_, 3L)
    expected_mean[old_active] <- original$active_mean[[old_state]]
    expect_equal(restored_mean, expected_mean, tolerance = 1e-12)

    restored_covariance <- matrix(NA_real_, 3L, 3L)
    restored_covariance[permutation[new_active], permutation[new_active]] <-
      permuted$active_covariance[[new_state]]
    expected_covariance <- matrix(NA_real_, 3L, 3L)
    expected_covariance[old_active, old_active] <-
      original$active_covariance[[old_state]]
    expect_equal(restored_covariance, expected_covariance, tolerance = 1e-12)
  }
})

test_that("general-T metadata is validated without compatibility bypass", {
  fixture <- general_t_fixture(3L)
  on.exit(blr_unified_cleanup(fixture), add = TRUE)
  raw <- general_t_fit(fixture, sampling_iterations = 2L)
  invalid <- raw
  invalid$input$prior$probability$activity_patterns[2L, ] <-
    invalid$input$prior$probability$activity_patterns[1L, ]
  for (compatibility in c("phase0-v2", "phase1-r-v2")) {
    invalid$schema$compatibility_id <- compatibility
    expect_error(sblr:::validate_blr_raw_v2(invalid),
                 "complete unique binary pattern matrix")
  }
  reordered <- raw
  reordered$input$prior$probability$activity_patterns <-
    reordered$input$prior$probability$activity_patterns[
      c(1L, 3L, 2L, 4:8), , drop = FALSE]
  expect_error(sblr:::validate_blr_raw_v2(reordered),
               "declared state and trait order|canonically")
  expect_error(sblr:::.blr_phase4a_patterns(paste0("trait", 1:13)),
               "8192 patterns")
})

test_that("general-T compact transition diagnostics are validated exactly", {
  fixture <- general_t_fixture(3L)
  on.exit(blr_unified_cleanup(fixture), add = TRUE)
  raw <- general_t_fit(
    fixture, chains = 2L, cores = 1L, sampling_iterations = 2L)
  validate <- sblr:::validate_blr_raw_v2
  expect_true(validate(raw))
  qualification <- raw$diagnostics$qualification
  patterns <- rownames(raw$input$prior$probability$activity_patterns)
  chains <- c("chain1", "chain2")
  expect_true("transition_counts" %in% names(qualification))
  expect_null(qualification$transition_counts)
  expect_identical(names(qualification$pattern_occupancy_counts), chains)
  expect_identical(names(qualification$pattern_change_counts), chains)
  expect_true(all(vapply(
    qualification$pattern_occupancy_counts,
    function(x) identical(names(x), patterns), logical(1))))

  expect_invalid <- function(mutator, message) {
    candidate <- raw
    candidate <- mutator(candidate)
    expect_error(validate(candidate), message)
  }
  expect_invalid(function(x) {
    x$diagnostics$qualification["transition_counts"] <- NULL
    x
  }, "missing required compact transition fields")
  expect_invalid(function(x) {
    x$diagnostics$qualification <- c(
      x$diagnostics$qualification,
      list(transition_counts = matrix(0, 8L, 8L)))
    x
  }, "compact transition field must occur exactly once")
  for (fabricated in list(matrix(0, 8L, 8L), numeric(8L), list())) {
    expect_invalid(function(x) {
      x$diagnostics$qualification$transition_counts <- fabricated
      x
    }, "transition_counts must be present with value NULL")
  }
  expect_invalid(function(x) {
    x$diagnostics$qualification["pattern_occupancy_counts"] <- NULL
    x
  }, "missing required compact transition fields")
  expect_invalid(function(x) {
    x$diagnostics$qualification$pattern_occupancy_counts[[1L]] <-
      x$diagnostics$qualification$pattern_occupancy_counts[[1L]][-1L]
    x
  }, "pattern_occupancy_counts.*exact declared IDs")
  expect_invalid(function(x) {
    names(x$diagnostics$qualification$pattern_occupancy_counts[[1L]]) <-
      rev(patterns)
    x
  }, "pattern_occupancy_counts.*exact declared IDs")
  expect_invalid(function(x) {
    names(x$diagnostics$qualification$pattern_occupancy_counts) <- rev(chains)
    x
  }, "exact declared chain axis")
  for (bad_count in list(-1, Inf, 0.5)) {
    expect_invalid(function(x) {
      x$diagnostics$qualification$pattern_occupancy_counts[[1L]][1L] <-
        bad_count
      x
    }, "finite nonnegative integer-valued counts")
  }
  expect_invalid(function(x) {
    x$diagnostics$qualification$pattern_change_counts <-
      x$diagnostics$qualification$pattern_change_counts[-1L]
    x
  }, "pattern_change_counts.*exact declared IDs")
  for (bad_count in list(-1, NaN, 0.25)) {
    expect_invalid(function(x) {
      x$diagnostics$qualification$pattern_change_counts[[1L]] <- bad_count
      x
    }, "finite nonnegative integer-valued counts")
  }
  expect_invalid(function(x) {
    x$diagnostics$qualification$pattern_occupancy_counts[[1L]][1L] <-
      x$diagnostics$qualification$pattern_occupancy_counts[[1L]][1L] + 1
    x
  }, "occupancy totals")
  expect_invalid(function(x) {
    x$diagnostics$qualification$pattern_change_counts[[1L]] <- 10
    x
  }, "cannot exceed")

  for (compatibility in c("phase0-v2", "phase1-r-v2")) {
    candidate <- raw
    candidate$schema$compatibility_id <- compatibility
    candidate$diagnostics$qualification$transition_counts <- matrix(0, 8L, 8L)
    expect_error(validate(candidate),
                 "transition_counts must be present with value NULL")
  }
  expect_invalid(function(x) {
    x$input$compute$operator_numerical_controls[
      "transition_diagnostic_policy"] <- NULL
    x
  }, "requires transition_diagnostic_policy")
})

test_that("general-T qualification remains internal and public MT is unchanged", {
  expect_false(".blr_cheng_mt_bayesc_bed_qualification" %in%
                 getNamespaceExports("sblr"))
  expect_false(any(grepl("phase4a|cheng_mt_bayesc_bed_qualification",
                         deparse(body(sblr::mtblr_bed)), fixed = FALSE)))
  expect_false(any(grepl("phase4a|cheng_mt_bayesc_bed_qualification",
                         deparse(body(sblr::mtblr_csr)), fixed = FALSE)))
  expect_false(any(grepl("phase4a|cheng_mt_bayesc_bed_qualification",
                         deparse(body(sblr::mtblr_block_eigen)), fixed = FALSE)))
})
