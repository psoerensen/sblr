# Permanent qualification fixture for the corrected two-trait Cheng BED slice.
phase4a_fixture <- function() {
  fixture <- blr_unified_fixture()
  phenotype <- cbind(
    trait1 = fixture$y[, 1L],
    trait2 = 0.35 * fixture$y[, 1L] +
      c(0.20, -0.10, 0.30, -0.20, 0.10, 0.05, -0.15, 0.25))
  rownames(phenotype) <- rownames(fixture$y)
  list(
    bed = fixture$bed,
    dosage = fixture$dosage,
    Glist = fixture$Glist,
    phenotype = phenotype,
    Ve = matrix(c(1.00, 0.25, 0.25, 0.80), 2L),
    Vb = matrix(c(0.35, 0.10, 0.10, 0.30), 2L),
    Psi = matrix(c(0.50, 0.08, 0.08, 0.45), 2L))
}

phase4b_fit <- function(fixture, ..., cores = 1L, chains = 1L,
                        keep_traces = TRUE) {
  sblr:::.blr_phase4a_cheng_mt_bed(
    y = fixture$phenotype, Glist = fixture$Glist,
    fixed_residual_covariance = NULL,
    initial_marker_covariance = fixture$Vb,
    marker_covariance_prior_df = 4.5,
    marker_covariance_prior_scale = fixture$Psi,
    burn_in_iterations = 2L, sampling_iterations = 4L,
    thin_interval = 2L, chains = chains, cores = cores,
    seed = 17, keep_traces = keep_traces,
    residual_covariance_policy = "sampled_full",
    initial_residual_covariance = fixture$Ve,
    residual_covariance_prior_df = 1.5,
    residual_covariance_prior_scale =
      matrix(c(.60, .12, .12, .55), 2L), ...)
}

phase4a_fit <- function(fixture, ..., cores = 1L, chains = 1L,
                        keep_traces = TRUE) {
  sblr:::.blr_phase4a_cheng_mt_bed(
    y = fixture$phenotype, Glist = fixture$Glist,
    fixed_residual_covariance = fixture$Ve,
    initial_marker_covariance = fixture$Vb,
    marker_covariance_prior_df = 4.5,
    marker_covariance_prior_scale = fixture$Psi,
    burn_in_iterations = 2L, sampling_iterations = 4L,
    thin_interval = 2L, chains = chains, cores = cores,
    seed = 17, keep_traces = keep_traces, ...)
}

phase4a_pattern_oracle <- function(score, marker_sum_squares, Vb, Ve,
                                   probability) {
  patterns <- rbind(c(0L, 0L), c(1L, 0L), c(0L, 1L), c(1L, 1L))
  omega <- solve(Ve)
  h <- as.numeric(omega %*% score)
  log_weight <- log(probability)
  means <- covariances <- vector("list", 4L)
  for (state in 2:4) {
    active <- which(patterns[state, ] == 1L)
    prior <- Vb[active, active, drop = FALSE]
    precision <- solve(prior) +
      marker_sum_squares * omega[active, active, drop = FALSE]
    covariance <- solve(precision)
    mean <- as.numeric(covariance %*% h[active])
    log_weight[state] <- log_weight[state] -
      0.5 * as.numeric(determinant(prior, logarithm = TRUE)$modulus) -
      0.5 * as.numeric(determinant(precision, logarithm = TRUE)$modulus) +
      0.5 * sum(h[active] * mean)
    means[[state]] <- mean
    covariances[[state]] <- covariance
  }
  weight <- exp(log_weight - max(log_weight))
  list(probability = weight / sum(weight), mean = means,
       covariance = covariances)
}

phase4a_relabel_pattern_axes <- function(raw, state_ids) {
  raw$input$model$state_space <- state_ids
  raw$model$state_space <- state_ids
  for (group in c("posterior", "draws", "final")) {
    for (field in names(raw[[group]])) {
      value <- raw[[group]][[field]]
      if (!is.array(value)) next
      axes <- attr(value, "dim_axis_names")
      for (axis in c("joint_state", "activity_pattern")) {
        location <- match(axis, axes)
        if (!is.na(location)) dimnames(value)[[location]] <- state_ids
      }
      raw[[group]][[field]] <- value
    }
  }
  raw
}

test_that("Phase 4a marker contract matches an independent four-pattern oracle", {
  score <- c(0.55, -0.30)
  c_j <- 1.7
  Vb <- matrix(c(0.45, 0.16, 0.16, 0.38), 2L)
  Ve <- matrix(c(1.20, 0.31, 0.31, 0.85), 2L)
  probability <- c(0.52, 0.13, 0.11, 0.24)
  expected <- phase4a_pattern_oracle(score, c_j, Vb, Ve, probability)
  observed <- sblr:::mtblr_phase4a_pattern_contract_internal(
    score, c_j, Vb, Ve, probability)

  expect_identical(observed$pattern_order, c("0_0", "1_0", "0_1", "1_1"))
  expect_equal(observed$probability, expected$probability, tolerance = 1e-13)
  for (state in 2:4) {
    expect_equal(as.numeric(observed$active_mean[[state]]),
                 expected$mean[[state]],
                 tolerance = 1e-13)
    expect_equal(observed$active_covariance[[state]],
                 expected$covariance[[state]], tolerance = 1e-13)
  }
  expect_equal(observed$completion_coefficient,
               c(Vb[2, 1] / Vb[1, 1], Vb[1, 2] / Vb[2, 2]),
               tolerance = 1e-15)
  expect_equal(observed$completion_variance,
               c(Vb[2, 2] - Vb[2, 1]^2 / Vb[1, 1],
                 Vb[1, 1] - Vb[1, 2]^2 / Vb[2, 2]),
               tolerance = 1e-15)

  diagonal <- sblr:::mtblr_phase4a_pattern_contract_internal(
    score, c_j, Vb, diag(diag(Ve)), probability)
  expect_false(isTRUE(all.equal(observed$probability, diagonal$probability,
                                tolerance = 1e-12)))
})

test_that("Phase 4a one-marker Monte Carlo frequencies match exact weights", {
  fixture <- phase4a_fixture()
  on.exit(unlink(fixture$bed), add = TRUE)
  fixture$Glist$rsids <- fixture$Glist$rsidsLD <- list("m1")
  fixture$Glist$chr <- list(1L)
  fixture$Glist$pos <- list(100)
  fixture$Glist$af <- fixture$Glist$af[1L]
  fixture$Glist$maf <- fixture$Glist$maf[1L]
  probability <- c(.52, .13, .11, .24)
  p <- fixture$Glist$af[[1L]][1L]
  x <- (c(0, 1, 2, 0, 1, 2, 1, 0) - 2 * p) /
    sqrt(2 * p * (1 - p))
  score <- as.numeric(crossprod(x, fixture$phenotype))
  expected <- phase4a_pattern_oracle(
    score, sum(x^2), fixture$Vb, fixture$Ve, probability)$probability
  fit <- sblr:::.blr_phase4a_cheng_mt_bed(
    fixture$phenotype, fixture$Glist, fixture$Ve, fixture$Vb,
    4.5, fixture$Psi,
    initial_activity_pattern_probability = probability,
    update_marker_covariance = FALSE,
    update_activity_pattern_probability = FALSE,
    burn_in_iterations = 0L, sampling_iterations = 12000L,
    thin_interval = 1L, chains = 1L, cores = 1L, seed = 104)
  observed <- drop(fit$posterior$activity_pattern_probabilities[1L, ])
  expect_equal(unname(observed), expected, tolerance = 0.018)
})

test_that("Phase 4a learned pattern moments match one-marker enumeration", {
  fixture <- phase4a_fixture()
  on.exit(unlink(fixture$bed), add = TRUE)
  fixture$Glist$rsids <- fixture$Glist$rsidsLD <- list("m1")
  fixture$Glist$chr <- list(1L)
  fixture$Glist$pos <- list(100)
  fixture$Glist$af <- fixture$Glist$af[1L]
  fixture$Glist$maf <- fixture$Glist$maf[1L]
  shape <- c(4.0, 1.2, 0.9, 2.1)
  prior_probability <- shape / sum(shape)
  p <- fixture$Glist$af[[1L]][1L]
  x <- (c(0, 1, 2, 0, 1, 2, 1, 0) - 2 * p) /
    sqrt(2 * p * (1 - p))
  state_probability <- phase4a_pattern_oracle(
    as.numeric(crossprod(x, fixture$phenotype)), sum(x^2),
    fixture$Vb, fixture$Ve, prior_probability)$probability
  expected_pi <- (shape + state_probability) / (sum(shape) + 1)
  fit <- sblr:::.blr_phase4a_cheng_mt_bed(
    fixture$phenotype, fixture$Glist, fixture$Ve, fixture$Vb,
    4.5, fixture$Psi,
    initial_activity_pattern_probability = prior_probability,
    activity_pattern_dirichlet_prior = shape,
    update_marker_covariance = FALSE,
    update_activity_pattern_probability = TRUE,
    burn_in_iterations = 1000L, sampling_iterations = 16000L,
    thin_interval = 1L, chains = 1L, cores = 1L, seed = 37)
  observed_pi <- colMeans(fit$draws$activity_pattern_parameters[, 1L, ])
  expect_equal(unname(observed_pi), expected_pi, tolerance = 0.018)
  expect_equal(
    unname(fit$posterior$activity_pattern_probabilities[1L, ]),
    state_probability, tolerance = 0.022)
})

test_that("Phase 4a produces coherent qualification-only raw-v2 state", {
  fixture <- phase4a_fixture()
  on.exit(unlink(fixture$bed), add = TRUE)
  fit <- phase4a_fit(fixture, chains = 2L)

  expect_true(sblr:::is_blr_raw_v2(fit))
  expect_true(fit$input$prior$marker_covariance$sampled)
  expect_equal(fit$input$prior$marker_covariance$degrees_of_freedom, 4.5)
  expect_null(fit$input$prior$marker_covariance$fixed_value)
  expect_equal(fit$diagnostics$qualification$initial_marker_covariance,
               fixture$Vb, ignore_attr = TRUE)
  expect_identical(fit$model$analysis_mode, "joint_multitrait")
  expect_length(fit$input$data$operator_resources, 1L)
  expect_length(fit$input$data$providers, 1L)
  provider <- fit$input$data$providers[[1L]]
  resource <- fit$input$data$operator_resources[[1L]]
  expect_identical(provider$trait_ids, c("trait1", "trait2"))
  expect_identical(provider$operator_resource_id, resource$resource_id)
  expect_identical(provider$likelihood_regime, "common_sample")
  expect_identical(fit$schema$seed_contract_version, 1L)
  expect_identical(fit$schema$retention_contract_version, 1L)
  expect_identical(fit$input$mcmc$retained_transition_indices, c(2L, 4L))
  expect_identical(dim(fit$draws$realised_effects), c(2L, 2L, 3L, 2L))
  expect_identical(dim(fit$draws$joint_states), c(2L, 2L, 3L))
  expect_identical(dim(fit$draws$marker_covariance), c(2L, 2L, 2L, 2L))
  expect_identical(dim(fit$draws$activity_pattern_parameters),
                   c(2L, 2L, 4L))
  expect_identical(dim(fit$final$marker_covariance), c(2L, 2L, 2L))
  expect_identical(dim(fit$derived$predictions), c(2L, 2L, 8L, 2L))
  expect_identical(fit$diagnostics$qualification$status,
                   "qualification_only")
  expect_false(fit$diagnostics$qualification$current_legacy_mt_route_used)

  patterns <- rbind(c(0, 0), c(1, 0), c(0, 1), c(1, 1))
  state <- fit$draws$joint_states
  for (draw in seq_len(dim(state)[1L])) {
    for (chain in seq_len(dim(state)[2L])) {
      for (marker in seq_len(dim(state)[3L])) {
        index <- state[draw, chain, marker] + 1L
        latent <- fit$draws$latent_effects[draw, chain, marker, ]
        realised <- fit$draws$realised_effects[draw, chain, marker, ]
        if (index == 1L) {
          expect_true(all(is.na(latent)))
          expect_identical(as.numeric(realised), c(0, 0))
        } else {
          expect_true(all(is.finite(latent)))
          expect_equal(realised, patterns[index, ] * latent,
                       tolerance = 1e-14)
        }
      }
    }
  }
  expect_equal(unname(rowSums(fit$posterior$activity_pattern_probabilities)),
               rep(1, 3L), tolerance = 1e-15)
  expect_equal(unname(apply(fit$draws$activity_pattern_parameters,
                            c(1L, 2L), sum)),
               matrix(1, 2L, 2L), tolerance = 1e-14)
  for (draw in seq_len(dim(fit$draws$marker_covariance)[1L])) {
    for (chain in seq_len(dim(fit$draws$marker_covariance)[2L])) {
      value <- fit$draws$marker_covariance[draw, chain, , ]
      expect_equal(value, t(value), tolerance = 1e-14)
      expect_silent(chol(value))
    }
  }
})

test_that("Phase 4a pleiotropic probabilities obey the joint-pattern contract", {
  fixture <- phase4a_fixture()
  on.exit(unlink(fixture$bed), add = TRUE)
  fit <- phase4a_fit(fixture)
  validate <- sblr:::validate_blr_raw_v2

  expect_true(validate(fit))
  expect_identical(
    as.numeric(fit$posterior$pleiotropic_probabilities),
    as.numeric(fit$posterior$activity_pattern_probabilities[, "1_1"]))

  missing <- fit
  missing$posterior$pleiotropic_probabilities <- NULL
  expect_error(validate(missing),
               "pleiotropic_probabilities is required")

  invalid_value <- function(value) {
    candidate <- fit
    candidate$posterior$pleiotropic_probabilities[1L] <- value
    candidate
  }
  expect_error(validate(invalid_value(1.01)), "must lie in \\[0, 1\\]")
  expect_error(validate(invalid_value(-0.01)), "must lie in \\[0, 1\\]")
  for (value in list(NA_real_, NaN, Inf, -Inf)) {
    expect_error(validate(invalid_value(value)), "finite numeric")
  }
  logical_value <- fit
  logical_value$posterior$pleiotropic_probabilities <- sblr:::.blr_make_array(
    rep(TRUE, length(fit$input$data$global_markers)),
    list(marker = fit$input$data$global_markers))
  expect_error(validate(logical_value), "finite numeric")

  wrong_axis <- fit
  wrong_axis$posterior$pleiotropic_probabilities <- sblr:::.blr_make_array(
    c(.1, .2), list(marker = fit$input$data$global_markers[1:2]))
  expect_error(validate(wrong_axis), "Raw axis IDs")

  inconsistent <- fit
  current <- inconsistent$posterior$pleiotropic_probabilities[1L]
  inconsistent$posterior$pleiotropic_probabilities[1L] <-
    current + if (current <= .5) .25 else -.25
  expect_error(validate(inconsistent), "must equal.*1_1")

  invalid_pattern <- phase4a_relabel_pattern_axes(
    fit, c("0_0", "1_0", "0_1", "1_x"))
  expect_error(validate(invalid_pattern),
               "identifiable binary activity-pattern IDs")
  missing_pattern <- fit
  missing_pattern$input$model$state_space <- c("0_0", "1_0", "0_1")
  missing_pattern$model$state_space <- c("0_0", "1_0", "0_1")
  expect_error(validate(missing_pattern),
               "exactly one identifiable pleiotropic")

  bad_labels <- fit
  dimnames(bad_labels$posterior$activity_pattern_probabilities)[[2L]] <-
    c("0_0", "1_0", "0_1", "0_1")
  expect_error(validate(bad_labels), "Raw axis IDs")

  zero <- fit
  zero$posterior$activity_pattern_probabilities[] <- 0
  zero$posterior$activity_pattern_probabilities[, "0_0"] <- 1
  zero$posterior$pleiotropic_probabilities[] <- 0
  expect_true(validate(zero))

  one <- fit
  one$posterior$activity_pattern_probabilities[] <- 0
  one$posterior$activity_pattern_probabilities[, "1_1"] <- 1
  one$posterior$pleiotropic_probabilities[] <- 1
  expect_true(validate(one))

  for (compatibility_id in c("phase0-v2", "phase1-r-v2")) {
    bypass <- invalid_value(2)
    bypass$schema$compatibility_id <- compatibility_id
    expect_error(validate(bypass), "must lie in \\[0, 1\\]")
  }
})

test_that("Phase 4a fixed controls and diagnostics are zero-RNG", {
  fixture <- phase4a_fixture()
  on.exit(unlink(fixture$bed), add = TRUE)
  arguments <- list(
    y = fixture$phenotype, Glist = fixture$Glist,
    fixed_residual_covariance = fixture$Ve,
    initial_marker_covariance = fixture$Vb,
    marker_covariance_prior_df = 4.5,
    marker_covariance_prior_scale = fixture$Psi,
    initial_activity_pattern_probability = c(.55, .12, .11, .22),
    update_marker_covariance = FALSE,
    update_activity_pattern_probability = FALSE,
    burn_in_iterations = 1L, sampling_iterations = 3L,
    thin_interval = 1L, chains = 2L, cores = 1L, seed = 0)
  with_trace <- do.call(sblr:::.blr_phase4a_cheng_mt_bed,
                        c(arguments, list(keep_traces = TRUE)))
  without_trace <- do.call(sblr:::.blr_phase4a_cheng_mt_bed,
                           c(arguments, list(keep_traces = FALSE)))

  expect_identical(with_trace$draws$realised_effects,
                   without_trace$draws$realised_effects)
  expect_identical(with_trace$draws$joint_states,
                   without_trace$draws$joint_states)
  expect_identical(with_trace$final$realised_effects,
                   without_trace$final$realised_effects)
  expect_true(is.null(without_trace$draws$convergence))
  expect_true(is.null(without_trace$diagnostics$convergence))
  expect_false(with_trace$input$prior$marker_covariance$sampled)
  expect_null(with_trace$input$prior$marker_covariance$degrees_of_freedom)
  expect_null(with_trace$input$prior$marker_covariance$scale)
  expect_equal(with_trace$input$prior$marker_covariance$fixed_value,
               fixture$Vb, ignore_attr = TRUE)
  for (draw in 1:3) for (chain in 1:2) {
    expect_equal(with_trace$draws$marker_covariance[draw, chain, , ],
                 fixture$Vb, tolerance = 0, ignore_attr = TRUE)
  }
  expected_probability <- c(.55, .12, .11, .22)
  expect_equal(unname(drop(
                 with_trace$draws$activity_pattern_parameters[1, 1, ])),
               expected_probability, tolerance = 0)
})

test_that("Phase 4a chain scheduling is deterministic across workers", {
  fixture <- phase4a_fixture()
  on.exit(unlink(fixture$bed), add = TRUE)
  serial <- phase4a_fit(fixture, chains = 2L, cores = 1L)
  repeat_serial <- phase4a_fit(fixture, chains = 2L, cores = 1L)
  parallel <- phase4a_fit(fixture, chains = 2L, cores = 2L)
  scientific_fields <- c(
    "realised_effects", "latent_effects", "joint_states",
    "activity_pattern_parameters", "marker_covariance")
  for (field in scientific_fields) {
    expect_identical(serial$draws[[field]], repeat_serial$draws[[field]])
    expect_identical(serial$draws[[field]], parallel$draws[[field]])
  }
  expect_identical(serial$final, repeat_serial$final)
  expect_identical(serial$final, parallel$final)
  expect_identical(serial$input$mcmc$task_seeds,
                   parallel$input$mcmc$task_seeds)
  expect_identical(parallel$diagnostics$workers$logical_task_order,
                   c("chain:0", "chain:1"))
  expect_identical(parallel$diagnostics$workers$diagnostics_rng_draws, 0L)
  capability <- sblr:::sparseLD_thread_info(2L)
  if (isTRUE(capability$openmp) &&
      capability$actual_threads_requested_region > 1L) {
    expect_gt(parallel$diagnostics$workers$actual_team_size, 1L)
    expect_gt(length(unique(
      parallel$diagnostics$workers$task_worker_ids)), 1L)
  }
})

test_that("Phase 4a authoritative covariance statistic excludes null markers", {
  fixture <- phase4a_fixture()
  on.exit(unlink(fixture$bed), add = TRUE)
  fit <- phase4a_fit(fixture)
  update <- fit$diagnostics$qualification$covariance_updates[[1L]]
  state <- fit$final$joint_states[1L, ]
  latent <- fit$final$latent_effects[1L, , ]
  included <- which(state != 0L)
  expected <- if (length(included)) {
    crossprod(latent[included, , drop = FALSE])
  } else matrix(0, 2L, 2L)
  expect_equal(update$statistic, expected, tolerance = 1e-13,
               ignore_attr = TRUE)
  expect_equal(update$active_marker_count, length(included))
  expect_equal(update$degrees_of_freedom, 4.5 + length(included))
  expect_equal(update$posterior_scale, fixture$Psi + expected,
               tolerance = 1e-13, ignore_attr = TRUE)
  expect_identical(
    drop(fit$final$marker_covariance[1L, , ]),
    drop(fit$diagnostics$convergence$marker_covariance[4L, 1L, , ]))
})

test_that("Phase 4a no-active covariance update is the proper prior reduction", {
  fixture <- phase4a_fixture()
  on.exit(unlink(fixture$bed), add = TRUE)
  fixture$Glist$rsids <- fixture$Glist$rsidsLD <- list("m1")
  fixture$Glist$chr <- list(1L)
  fixture$Glist$pos <- list(100)
  fixture$Glist$af <- fixture$Glist$af[1L]
  fixture$Glist$maf <- fixture$Glist$maf[1L]
  fit <- sblr:::.blr_phase4a_cheng_mt_bed(
    fixture$phenotype, fixture$Glist, fixture$Ve, fixture$Vb,
    1.5, fixture$Psi,
    initial_activity_pattern_probability =
      c(1 - 3e-12, 1e-12, 1e-12, 1e-12),
    update_marker_covariance = TRUE,
    update_activity_pattern_probability = FALSE,
    burn_in_iterations = 0L, sampling_iterations = 1L,
    thin_interval = 1L, chains = 1L, cores = 1L, seed = 91)
  update <- fit$diagnostics$qualification$covariance_updates[[1L]]
  expect_identical(update$active_marker_count, 0L)
  expect_equal(update$degrees_of_freedom, 1.5)
  expect_equal(update$statistic, matrix(0, 2L, 2L), tolerance = 0)
  expect_equal(update$posterior_scale, fixture$Psi, tolerance = 0,
               ignore_attr = TRUE)
  expect_silent(chol(drop(fit$draws$marker_covariance[1L, 1L, , ])))
})

test_that("Phase 4a validates its narrow fixed-residual boundary", {
  fixture <- phase4a_fixture()
  on.exit(unlink(fixture$bed), add = TRUE)
  invalid_ve <- list(
    y = fixture$phenotype, Glist = fixture$Glist,
    fixed_residual_covariance = matrix(c(1, 2, 2, 1), 2L),
    initial_marker_covariance = fixture$Vb,
    marker_covariance_prior_df = 4.5,
    marker_covariance_prior_scale = fixture$Psi,
    burn_in_iterations = 2L, sampling_iterations = 4L,
    thin_interval = 2L, chains = 1L, cores = 1L, seed = 17)
  expect_error(
    do.call(sblr:::.blr_phase4a_cheng_mt_bed, invalid_ve),
    "strictly positive definite")
  expect_error(
    sblr:::.blr_phase4a_cheng_mt_bed(
      fixture$phenotype[, 1L, drop = FALSE], fixture$Glist,
      fixture$Ve, fixture$Vb, 4.5, fixture$Psi),
    "exactly two traits|N x 2")
  invalid_df <- invalid_ve
  invalid_df$fixed_residual_covariance <- fixture$Ve
  invalid_df$marker_covariance_prior_df <- 1
  expect_error(
    do.call(sblr:::.blr_phase4a_cheng_mt_bed, invalid_df),
    "exceed T - 1")
  invalid_probability <- invalid_df
  invalid_probability$marker_covariance_prior_df <- 4.5
  invalid_probability$initial_activity_pattern_probability <-
    c(.5, .2, .2, .2)
  expect_error(
    do.call(sblr:::.blr_phase4a_cheng_mt_bed, invalid_probability),
    "sum to one")
  wrong_order <- c(`1_0` = .2, `0_0` = .5, `0_1` = .1, `1_1` = .2)
  invalid_probability$initial_activity_pattern_probability <- wrong_order
  expect_error(
    do.call(sblr:::.blr_phase4a_cheng_mt_bed, invalid_probability),
    "canonical activity-pattern order")
  misordered_ve <- fixture$Ve
  dimnames(misordered_ve) <- list(c("trait2", "trait1"),
                                  c("trait2", "trait1"))
  invalid_probability$initial_activity_pattern_probability <-
    c(.5, .2, .1, .2)
  invalid_probability$fixed_residual_covariance <- misordered_ve
  expect_error(
    do.call(sblr:::.blr_phase4a_cheng_mt_bed, invalid_probability),
    "declared trait order")
})

test_that("Phase 4b residual inverse-Wishart conditional uses N and E crossproduct", {
  residual <- cbind(
    c(-.8, .2, .5, -.1, .4, -.3),
    c(.3, -.4, .6, .2, -.2, .5))
  prior_df <- 1.5
  prior_scale <- matrix(c(.7, .14, .14, .6), 2L)
  expected_statistic <- crossprod(residual)
  expected_df <- prior_df + nrow(residual)
  expected_scale <- prior_scale + expected_statistic
  observed <- sblr:::mtblr_phase4b_residual_covariance_contract_internal(
    residual, prior_df, prior_scale, 6000L, 4294967295)

  expect_equal(observed$statistic, expected_statistic, tolerance = 1e-15)
  expect_equal(observed$degrees_of_freedom, expected_df, tolerance = 0)
  expect_equal(observed$scale, expected_scale, tolerance = 1e-15)
  expect_identical(dim(observed$draws), c(6000L, 2L, 2L))
  expected_mean <- expected_scale / (expected_df - 3)
  expect_equal(apply(observed$draws, c(2L, 3L), mean), expected_mean,
               tolerance = 0.04)
  for (draw in c(1L, 100L, 6000L)) {
    value <- observed$draws[draw, , ]
    expect_equal(value, t(value), tolerance = 0)
    expect_silent(chol(value))
  }
})

test_that("Phase 4b sampled residual covariance is one authoritative raw-v2 state", {
  fixture <- phase4a_fixture()
  on.exit(unlink(fixture$bed), add = TRUE)
  fit <- phase4b_fit(fixture, chains = 2L)

  expect_true(sblr:::validate_blr_raw_v2(fit))
  expect_identical(fit$model$residual_policy, "sampled_full")
  expect_identical(fit$model$update_order_version, 2L)
  expect_true(fit$input$prior$residual_covariance$sampled)
  expect_null(fit$input$prior$residual_covariance$fixed_value)
  expect_equal(fit$input$prior$residual_covariance$initial_value,
               fixture$Ve, ignore_attr = TRUE)
  expect_identical(
    fit$input$prior$residual_covariance$parameterization,
    "degrees_of_freedom_scale")
  expect_true(fit$input$prior$residual_covariance$proper)
  expect_false(fit$input$prior$residual_covariance$finite_mean)
  expect_true(fit$input$mcmc$update_flags$residual_covariance)
  expect_identical(dim(fit$draws$residual_covariance), c(2L, 2L, 2L, 2L))
  expect_identical(dim(fit$posterior$residual_covariance_mean), c(2L, 2L))
  expect_identical(dim(fit$final$residual_covariance), c(2L, 2L, 2L))
  expect_identical(
    dim(fit$draws$convergence$residual_covariance), c(4L, 2L, 2L, 2L))
  expect_equal(
    fit$posterior$residual_covariance_mean,
    apply(fit$draws$residual_covariance, c(3L, 4L), mean),
    tolerance = 1e-15, ignore_attr = TRUE)

  frequency <- rowMeans(fixture$dosage) / 2
  genotype <- t((fixture$dosage - 2 * frequency) /
                  sqrt(2 * frequency * (1 - frequency)))
  for (chain in 1:2) {
    realised <- fit$final$realised_effects[chain, , ]
    residual <- fixture$phenotype - genotype %*% realised
    update <- fit$diagnostics$qualification$
      residual_covariance_updates[[chain]]
    expect_equal(update$statistic, crossprod(residual), tolerance = 1e-12,
                 ignore_attr = TRUE)
    expect_equal(update$degrees_of_freedom, 1.5 + nrow(residual),
                 tolerance = 0)
    expect_equal(update$posterior_scale,
                 matrix(c(.60, .12, .12, .55), 2L) + crossprod(residual),
                 tolerance = 1e-12, ignore_attr = TRUE)
    expect_identical(update$update_count, 6L)
    expect_equal(
      fit$final$residual_covariance[chain, , ],
      fit$draws$convergence$residual_covariance[4L, chain, , ],
      tolerance = 0, ignore_attr = TRUE)
  }
  for (draw in 1:2) for (chain in 1:2) {
    expect_silent(chol(fit$draws$residual_covariance[draw, chain, , ]))
  }
})

test_that("Phase 4b scheduling and observational capture consume no RNG", {
  fixture <- phase4a_fixture()
  on.exit(unlink(fixture$bed), add = TRUE)
  serial <- phase4b_fit(fixture, chains = 2L, cores = 1L)
  repeated <- phase4b_fit(fixture, chains = 2L, cores = 1L)
  parallel <- phase4b_fit(fixture, chains = 2L, cores = 2L)
  without_trace <- phase4b_fit(
    fixture, chains = 2L, cores = 1L, keep_traces = FALSE)
  scientific <- c(
    "realised_effects", "latent_effects", "joint_states",
    "activity_pattern_parameters", "marker_covariance",
    "residual_covariance")
  for (field in scientific) {
    expect_identical(serial$draws[[field]], repeated$draws[[field]])
    expect_identical(serial$draws[[field]], parallel$draws[[field]])
    expect_identical(serial$draws[[field]], without_trace$draws[[field]])
  }
  expect_identical(serial$final, parallel$final)
  expect_identical(serial$final, without_trace$final)
  expect_null(without_trace$draws$convergence)
  expect_identical(parallel$diagnostics$workers$diagnostics_rng_draws, 0L)
  capability <- sblr:::sparseLD_thread_info(2L)
  if (isTRUE(capability$openmp) &&
      capability$actual_threads_requested_region > 1L) {
    expect_gt(parallel$diagnostics$workers$actual_team_size, 1L)
    expect_gt(length(unique(
      parallel$diagnostics$workers$task_worker_ids)), 1L)
  }
})

test_that("Phase 4b thinning selects draws without changing the chain", {
  fixture <- phase4a_fixture()
  on.exit(unlink(fixture$bed), add = TRUE)
  common <- list(
    y = fixture$phenotype, Glist = fixture$Glist,
    fixed_residual_covariance = NULL,
    initial_marker_covariance = fixture$Vb,
    marker_covariance_prior_df = 4.5,
    marker_covariance_prior_scale = fixture$Psi,
    burn_in_iterations = 2L, sampling_iterations = 6L,
    chains = 1L, cores = 1L, seed = 23, keep_traces = TRUE,
    residual_covariance_policy = "sampled_full",
    initial_residual_covariance = fixture$Ve,
    residual_covariance_prior_df = 1.5,
    residual_covariance_prior_scale = diag(c(.6, .55)))
  all_draws <- do.call(sblr:::.blr_phase4a_cheng_mt_bed,
                       c(common, list(thin_interval = 1L)))
  thinned <- do.call(sblr:::.blr_phase4a_cheng_mt_bed,
                     c(common, list(thin_interval = 2L)))
  expect_identical(all_draws$final, thinned$final)
  expect_identical(all_draws$draws$convergence,
                   thinned$draws$convergence)
  for (field in c("realised_effects", "latent_effects", "joint_states",
                  "activity_pattern_parameters", "marker_covariance",
                  "residual_covariance")) {
    value <- all_draws$draws[[field]]
    subscripts <- rep(list(TRUE), length(dim(value)))
    subscripts[[1L]] <- c(2L, 4L, 6L)
    selected <- do.call(`[`, c(list(value), subscripts, list(drop = FALSE)))
    expect_equal(selected, thinned$draws[[field]], tolerance = 0,
                 ignore_attr = TRUE)
  }
  expect_identical(thinned$input$mcmc$retained_transition_indices,
                   c(2L, 4L, 6L))
})

test_that("Phase 4b validates residual policy, prior, and raw presence rules", {
  fixture <- phase4a_fixture()
  on.exit(unlink(fixture$bed), add = TRUE)
  common <- list(
    y = fixture$phenotype, Glist = fixture$Glist,
    fixed_residual_covariance = NULL,
    initial_marker_covariance = fixture$Vb,
    marker_covariance_prior_df = 4.5,
    marker_covariance_prior_scale = fixture$Psi,
    burn_in_iterations = 0L, sampling_iterations = 2L,
    thin_interval = 1L, chains = 1L, cores = 1L, seed = 4,
    residual_covariance_policy = "sampled_full",
    initial_residual_covariance = fixture$Ve,
    residual_covariance_prior_df = 1.5,
    residual_covariance_prior_scale = diag(2))
  mutate_argument <- function(name, value) {
    candidate <- common
    candidate[name] <- list(value)
    candidate
  }
  expect_error(do.call(sblr:::.blr_phase4a_cheng_mt_bed,
                       mutate_argument("residual_covariance_prior_df", NULL)),
               "residual_covariance_prior_df")
  expect_error(do.call(sblr:::.blr_phase4a_cheng_mt_bed,
                       mutate_argument("residual_covariance_prior_df", 1)),
               "exceed T - 1")
  expect_error(do.call(sblr:::.blr_phase4a_cheng_mt_bed,
                       mutate_argument("residual_covariance_prior_scale",
                                       matrix(1, 2L, 3L))),
               "2 x 2")
  expect_error(do.call(sblr:::.blr_phase4a_cheng_mt_bed,
                       mutate_argument("residual_covariance_prior_scale",
                                       matrix(c(1, 2, 2, 1), 2L))),
               "strictly positive definite")
  expect_error(do.call(sblr:::.blr_phase4a_cheng_mt_bed,
                       mutate_argument("initial_residual_covariance",
                                       matrix(c(1, 2, 2, 1), 2L))),
               "strictly positive definite")
  fixed_with_prior <- common
  fixed_with_prior$residual_covariance_policy <- "fixed_full"
  fixed_with_prior$fixed_residual_covariance <- fixture$Ve
  expect_error(do.call(sblr:::.blr_phase4a_cheng_mt_bed, fixed_with_prior),
               "cannot carry sampled")

  fit <- do.call(sblr:::.blr_phase4a_cheng_mt_bed, common)
  invalid <- fit
  invalid$draws["residual_covariance"] <- list(NULL)
  expect_error(sblr:::validate_blr_raw_v2(invalid),
               "requires an initial state, retained draws")
  invalid <- fit
  invalid$posterior$residual_covariance_mean[1L, 1L] <- Inf
  expect_error(sblr:::validate_blr_raw_v2(invalid), "finite numeric")
  invalid <- fit
  invalid$final$residual_covariance[1L, 1L, 2L] <- 9
  expect_error(sblr:::validate_blr_raw_v2(invalid),
               "asymmetric|non-positive-definite")
  invalid <- fit
  dimnames(invalid$draws$residual_covariance)[[3L]] <- rev(
    dimnames(invalid$draws$residual_covariance)[[3L]])
  expect_error(sblr:::validate_blr_raw_v2(invalid), "Raw axis IDs")
  invalid <- fit
  invalid$schema$compatibility_id <- "phase0-v2"
  invalid$draws["residual_covariance"] <- list(NULL)
  expect_error(sblr:::validate_blr_raw_v2(invalid),
               "requires an initial state, retained draws")
})

test_that("Phase 4b leaves explicit fixed mode on the Phase 4a trajectory", {
  fixture <- phase4a_fixture()
  on.exit(unlink(fixture$bed), add = TRUE)
  default <- phase4a_fit(fixture)
  explicit <- phase4a_fit(fixture, residual_covariance_policy = "fixed_full")
  expect_identical(default, explicit)
  expect_identical(default$model$residual_policy, "fixed_full")
  expect_false(default$input$mcmc$update_flags$residual_covariance)
  expect_null(default$draws$residual_covariance)
  expect_null(default$posterior$residual_covariance_mean)
})

test_that("Phase 4b sampled chain agrees with the independent completed-active reference", {
  fixture <- phase4a_fixture()
  on.exit(unlink(fixture$bed), add = TRUE)
  research <- test_path("..", "research", "mtblr_covariance")
  reference <- new.env(parent = globalenv())
  for (file in c("mtblr_reference_model.R", "mtblr_pattern_reference.R",
                 "mtblr_pattern_samplers.R", "mtblr_sampled_residual.R")) {
    sys.source(file.path(research, file), envir = reference)
  }
  frequency <- rowMeans(fixture$dosage) / 2
  genotype <- t((fixture$dosage - 2 * frequency) /
                  sqrt(2 * frequency * (1 - frequency)))
  residual_scale <- matrix(c(.60, .12, .12, .55), 2L)
  native <- sblr:::.blr_phase4a_cheng_mt_bed(
    fixture$phenotype, fixture$Glist, NULL, fixture$Vb, 4.5, fixture$Psi,
    burn_in_iterations = 500L, sampling_iterations = 2500L,
    thin_interval = 1L, chains = 1L, cores = 1L, seed = 909,
    residual_covariance_policy = "sampled_full",
    initial_residual_covariance = fixture$Ve,
    residual_covariance_prior_df = 1.5,
    residual_covariance_prior_scale = residual_scale)
  completed <- reference$mtblr_pattern_sampler_sampled_residual(
    genotype, fixture$phenotype, reference$mt_pattern_space(), rep(1, 4L),
    4.5, fixture$Psi, 1.5, residual_scale, fixture$Vb, fixture$Ve,
    n_iter = 3000L, burn = 500L, seed = 1909L)
  summary <- reference$mt_pattern_summary(completed)

  # These tolerances were fixed from the scale of independent-chain Monte
  # Carlo uncertainty before comparing the two implementations.
  expect_equal(unname(native$posterior$pips),
               unname(summary$trait_pip), tolerance = .15,
               ignore_attr = TRUE)
  expect_equal(unname(native$posterior$realised_effect_mean),
               unname(summary$alpha_mean), tolerance = .14,
               ignore_attr = TRUE)
  expect_equal(unname(native$posterior$marker_covariance_mean),
               unname(summary$Vb_mean), tolerance = .20,
               ignore_attr = TRUE)
  expect_equal(
    unname(native$posterior$residual_covariance_mean),
    unname(apply(completed$Ve, c(1L, 2L), mean)), tolerance = .20,
    ignore_attr = TRUE)
  expect_equal(
    unname(colMeans(native$draws$activity_pattern_parameters[, 1L, ])),
    unname(summary$Pi_mean), tolerance = .12)
  expect_equal(
    unname(apply(native$derived$predictions[, 1L, , ], c(2L, 3L), mean)),
    unname(summary$fitted_mean), tolerance = .16)
})
