phase1_fixture_environment <- function() {
  environment <- new.env(parent = globalenv())
  sys.source(
    testthat::test_path("..", "research", "blr_framework_contract",
                        "blr_contract_fixtures.R"),
    envir = environment)
  environment
}

phase1_legacy_bayesc_raw <- function() {
  marker <- matrix(c(.1, 0, -.2), ncol = 1L)
  raw <- list(
    schema = list(class = "stblr_raw", version = 1L),
    meta = list(
      model = "bayesc", backend = "csr_bayesc", data_level = "summary",
      prior_type = "global", m = 3L, nt = 1L, n_trace = 5L,
      nit = 3L, nburn = 2L, nthin = 2L, nchains = 1L,
      keep_chains = FALSE, n_components = 2L, n_annotations = 0L,
      n_groups = 0L),
    marker = list(
      bm = marker, dm = matrix(c(.6, .1, .8), ncol = 1L), wy = marker,
      r = marker, b = marker, state = matrix(c(1, 0, 1), ncol = 1L)),
    trace = list(
      vbs = matrix(1:5 / 10, ncol = 1L),
      vgs = matrix(2:6 / 10, ncol = 1L),
      ves = matrix(6:10 / 10, ncol = 1L),
      vle = matrix(1:5 / 20, ncol = 1L),
      vld = matrix(1:5 / 25, ncol = 1L),
      pis = matrix(c(.2, .3, .4, .5, .6), ncol = 1L)),
    variance = list(
      covb = matrix(.5), covg = matrix(.6), cove = matrix(.7),
      vb = matrix(.5), vg = matrix(.6), ve = matrix(.7)),
    pi = list(final = matrix(c(.4, .6), nrow = 1L),
              mean = matrix(c(.5, .5), nrow = 1L),
              names = c("null", "active")),
    diagnostics = list(nsamples = 2L, n_used = 100L), chains = list(),
    prior = list(), group = list(), annotation = list(), component = list(),
    selection = list())
  class(raw) <- c("stblr_raw_v1", "stblr_raw", "list")
  raw
}

phase1_legacy_spec <- function(chains = 1L) {
  fixture <- phase1_fixture_environment()
  spec <- fixture$blr_spec_fixture("single_trait", "fixed_full")
  spec$data$operator_resources$bed_shared$operator_type <- "csr"
  spec$model$probability_policy <- "global"
  spec$model$marker_covariance_policy <- "traitwise_scalar"
  spec$model$residual_policy <- "scalar"
  spec$schema$compatibility_id <-
    "phase1-legacy-st-v1;seed=legacy_st_arithmetic_v1;retention=st_scalar_v1"
  spec$schema$seed_contract_version <- 0L
  spec$schema$retention_contract_version <- 0L
  spec$mcmc$burn_in_iterations <- 2L
  spec$mcmc$sampling_iterations <- 3L
  spec$mcmc$thin_interval <- 2L
  spec$mcmc$retained_transition_indices <- c(1L, 3L)
  spec$mcmc$retained_draws <- 2L
  spec$mcmc$chains <- as.integer(chains)
  spec$mcmc$task_seeds <- stats::setNames(
    seq_len(chains) + 100, paste0("chain", seq_len(chains)))
  spec
}

phase1_public_csr_prefix <- function(m) {
  prefix <- tempfile("phase1_csr_ld_")
  sblr:::.stblr_write_uint64_file(
    paste0(prefix, ".row_ptr.u64.bin"), rep(0, m + 1L))
  file.create(paste0(prefix, ".col_idx.u32.0based.bin"))
  file.create(paste0(prefix, ".values.f32.bin"))
  writeLines(c(
    "format=sparse_ld_csr", "storage=streamed_upper_triangle",
    "n_bed=NA", "n_used=NA", "n_samples_used=NA",
    paste0("n_variants=", m), "nnz=0", "triangle=upper",
    "diagonal=implicit_1",
    paste0("row_ptr_file=", prefix, ".row_ptr.u64.bin"),
    paste0("col_idx_file=", prefix, ".col_idx.u32.0based.bin"),
    paste0("values_file=", prefix, ".values.f32.bin"),
    "row_ptr_type=uint64", "col_idx_type=uint32", "values_type=float32",
    "index_base=0", "value=r"), paste0(prefix, ".meta.txt"))
  prefix
}

phase1_public_stats <- function(traits = c("T1")) {
  markers <- paste0("m", 1:3)
  wy <- lapply(seq_along(traits), function(index) {
    stats::setNames(c(2, -1, .5) * index, markers)
  })
  ww <- lapply(traits, function(.) stats::setNames(rep(50, 3), markers))
  names(wy) <- names(ww) <- traits
  list(wy = wy, ww = ww, yy = stats::setNames(rep(50, length(traits)), traits),
       n = 50L, m = 3L, marker_names = markers, trait_names = traits)
}

test_that("production validators accept the frozen Phase 0 fixtures", {
  fixture <- phase1_fixture_environment()
  for (mode in c("single_trait", "independent_traits", "joint_multitrait")) {
    spec <- fixture$blr_spec_fixture(
      mode, if (mode == "joint_multitrait") "sampled_full" else "fixed_full")
    expect_true(sblr:::validate_blr_resolved_spec(spec))
    expect_true(sblr:::validate_blr_raw_v2(fixture$blr_raw_fixture(mode)))
  }
  expect_true(sblr:::validate_blr_resolved_spec(
    fixture$blr_summary_spec_fixture()))
})

test_that("exact forwarded argument capture rejects unsafe names", {
  expect_error(sblr:::.blr_capture_forwarded_args(
    list(1), "pi_marker", what = "route(...)"), "unique, nonempty")
  duplicate <- structure(list(1, 2), names = c("pi_marker", "pi_marker"))
  expect_error(sblr:::.blr_capture_forwarded_args(
    duplicate, "pi_marker", what = "route(...)"), "unique, nonempty")
  expect_error(sblr:::.blr_capture_forwarded_args(
    list(pi_m = .2), "pi_marker", what = "route(...)"),
    "Unknown argument.*pi_m")
  expect_identical(sblr:::.blr_capture_forwarded_args(
    list(pi_marker = .2), "pi_marker", what = "route(...)"),
    list(pi_marker = .2))
  mapped <- sblr:::.blr_capture_forwarded_args(
    list(old_name = 2), "new_name", c(old_name = "new_name"), "route(...)")
  expect_identical(mapped[[1L]], 2)
  expect_identical(names(mapped), "new_name")
})

test_that("resolved spec rejects malformed policies and references", {
  fixture <- phase1_fixture_environment()
  spec <- fixture$blr_spec_fixture("independent_traits", "fixed_full")
  bad <- spec
  bad$compute$execution_mode <- "serial"
  bad$compute$parallelization <- "traits"
  expect_error(sblr:::validate_blr_resolved_spec(bad),
               "Invalid analysis/execution")
  bad <- spec
  bad$data$providers$p1$operator_resource_id <- "missing"
  expect_error(sblr:::validate_blr_resolved_spec(bad), "unknown operator")
  bad <- spec
  bad$mcmc$task_seeds <- c(1, 2)
  expect_error(sblr:::validate_blr_resolved_spec(bad), "task_seeds")
  bad <- spec
  bad$mcmc$retained_transition_indices <- 1L
  expect_error(sblr:::validate_blr_resolved_spec(bad), "retained indices")
  bad <- spec
  bad$prior$probability$pi <- .5
  expect_error(sblr:::validate_blr_resolved_spec(bad),
               "ambiguous probability name")
})

test_that("raw v2 retains size-one axes and probability meanings", {
  raw <- sblr:::convert_stblr_raw_v1_to_blr_raw_v2(
    phase1_legacy_bayesc_raw(), phase1_legacy_spec())
  expect_s3_class(raw, "blr_raw_v2")
  expect_true(sblr:::validate_blr_raw_v2(raw))
  expect_identical(dim(raw$posterior$realised_effect_mean), c(3L, 1L))
  expect_identical(attr(raw$posterior$realised_effect_mean, "dim_axis_names"),
                   c("marker", "trait"))
  expect_identical(dim(raw$draws$marker_variance), c(2L, 1L, 1L))
  expect_identical(dim(raw$final$realised_effects), c(1L, 3L, 1L))
  expect_identical(dim(raw$final$independent_trait_states), c(1L, 3L, 1L))
  expect_identical(raw$input$mcmc$retained_transition_indices, c(1L, 3L))
  expect_false(any(c("pi", "pis", "pim", "state_probabilities",
                     "pattern_probabilities") %in% names(raw$draws)))
  expect_true(all(c("latent_effects", "joint_states", "marker_covariance") %in%
                    names(raw$draws)))
  expect_null(raw$draws$latent_effects)
})

test_that("maintained public CSR fits attach validated raw v2 explicitly", {
  prefix <- phase1_public_csr_prefix(3L)
  on.exit(unlink(paste0(prefix, c(
    ".row_ptr.u64.bin", ".col_idx.u32.0based.bin", ".values.f32.bin",
    ".meta.txt"))), add = TRUE)
  fit <- stblr_csr(
    phase1_public_stats(c("T1", "T2")), ld_prefix = prefix,
    method = "sbayesc", nit = 3L, nburn = 1L, nthin = 2L,
    updateB = FALSE, updateE = FALSE, updatePi = FALSE,
    convergence = "none", seed = 31L)
  raw <- attr(fit, "blr_raw", exact = TRUE)
  spec <- attr(fit, "blr_resolved_spec", exact = TRUE)
  expect_s3_class(raw, "blr_raw_v2")
  expect_s3_class(spec, "blr_resolved_spec_v1")
  expect_true(sblr:::validate_blr_raw_v2(raw))
  expect_identical(spec$data$analysis_mode, "independent_traits")
  expect_identical(spec$data$providers$provider_1$sample_size, c(T1 = 50))
  expect_identical(spec$data$providers$provider_2$sample_size, c(T2 = 50))
  expect_true(is.numeric(spec$output$memory_estimate_bytes) &&
                length(spec$output$memory_estimate_bytes) == 1L &&
                is.finite(spec$output$memory_estimate_bytes))
  expect_true(all(c(
    "non_null_probability_initial",
    "non_null_probability_prior_shape_active",
    "non_null_probability_prior_shape_null") %in%
      names(spec$prior$probability$values)))
  expect_false(any(c("pi", "pis", "pim") %in%
                     names(spec$prior$probability$values)))
  expect_identical(dim(raw$posterior$pips), c(3L, 2L))
  expect_identical(dim(raw$draws$marker_variance), c(2L, 1L, 2L))
  expect_identical(dim(raw$final$realised_effects), c(1L, 3L, 2L))
  expect_identical(as.numeric(fit$bm),
                   as.numeric(raw$posterior$realised_effect_mean))
  expect_identical(as.numeric(fit$dm), as.numeric(raw$posterior$pips))
  expect_identical(fit$diagnostics$schema_v2_migration$status, "converted")
  consistency <- check_stblr_consistency(fit)
  expect_true(consistency$ok)
  expect_true(all(c("schema.blr_raw_v2", "schema.alias.bm",
                    "schema.alias.dm") %in% consistency$checks$check))
})

test_that("public BayesR conversion records component semantics without inventing draws", {
  prefix <- phase1_public_csr_prefix(3L)
  on.exit(unlink(paste0(prefix, c(
    ".row_ptr.u64.bin", ".col_idx.u32.0based.bin", ".values.f32.bin",
    ".meta.txt"))), add = TRUE)
  fit <- stblr_csr(
    phase1_public_stats("T1"), ld_prefix = prefix, method = "sbayesr",
    mixture_var = c(0, .2, 1), pi = c(.8, .1, .1),
    nit = 2L, nburn = 0L, updateB = FALSE, updateE = FALSE,
    updatePi = FALSE, convergence = "none", seed = 32L)
  raw <- attr(fit, "blr_raw", exact = TRUE)
  spec <- attr(fit, "blr_resolved_spec", exact = TRUE)
  expect_s3_class(raw, "blr_raw_v2")
  expect_identical(spec$model$state_space,
                   c("component_0", "component_1", "component_2"))
  expect_identical(unname(spec$prior$component_multipliers), c(0, .2, 1))
  expect_identical(
    attr(raw$posterior$traitwise_component_assignment_probabilities,
         "dim_axis_names"), c("marker", "trait", "component"))
  expect_identical(
    attr(raw$final$traitwise_component_probability_parameters,
         "dim_axis_names"), c("chain", "trait", "component"))
  expect_null(raw$draws$traitwise_component_probability_parameters)
  expect_identical(fit$diagnostics$schema_v2_migration$status, "converted")
})

test_that("supported fixed annotation fits use the same validated v2 bridge", {
  prefix <- phase1_public_csr_prefix(3L)
  on.exit(unlink(paste0(prefix, c(
    ".row_ptr.u64.bin", ".col_idx.u32.0based.bin", ".values.f32.bin",
    ".meta.txt"))), add = TRUE)
  annotation <- matrix(c(-1, 0, 1), ncol = 1L,
                       dimnames = list(paste0("m", 1:3), "A1"))
  fit <- stblr_csr_annot(
    phase1_public_stats("T1"),
    annotations = list(
      A = annotation, fixed_pi_marker = list(rep(.3, 3L)),
      fixed_vb_multiplier = list(rep(1, 3L))),
    annotation_model = "fixed_marker", ld_prefix = prefix,
    pi_init = .3, pi_prior_mean = .3, pi_prior_strength = 2,
    nit = 2L, nburn = 0L, updateB = FALSE, updateE = FALSE,
    updatePi = FALSE, convergence = "none", seed = 33L)
  raw <- attr(fit, "blr_raw", exact = TRUE)
  spec <- attr(fit, "blr_resolved_spec", exact = TRUE)
  expect_s3_class(raw, "blr_raw_v2")
  expect_true(sblr:::validate_blr_raw_v2(raw))
  expect_identical(spec$model$probability_policy, "fixed_marker")
  expect_identical(fit$diagnostics$schema_v2_migration$status, "converted")
})

test_that("public multi-chain ST remains explicit legacy when final states are unavailable", {
  prefix <- phase1_public_csr_prefix(3L)
  on.exit(unlink(paste0(prefix, c(
    ".row_ptr.u64.bin", ".col_idx.u32.0based.bin", ".values.f32.bin",
    ".meta.txt"))), add = TRUE)
  fit <- stblr_csr(
    phase1_public_stats("T1"), ld_prefix = prefix, method = "sbayesc",
    nit = 2L, nburn = 0L, nchains = 2L, ncores = 1L,
    updateB = FALSE, updateE = FALSE, updatePi = FALSE,
    convergence = "none", seed = 34L)
  expect_null(attr(fit, "blr_raw", exact = TRUE))
  expect_identical(fit$diagnostics$schema_v2_migration$status, "legacy")
  expect_match(fit$diagnostics$schema_v2_migration$reason,
               "does not expose every chain's final effect")
})

test_that("unsupported v1 semantics fail rather than being invented", {
  expect_error(sblr:::convert_stblr_raw_v1_to_blr_raw_v2(
    phase1_legacy_bayesc_raw(), phase1_legacy_spec(2L)),
    "does not expose every chain's final effect")
  mt <- list(schema = list(class = "mtblr_raw", version = 1L))
  expect_error(sblr:::convert_blr_raw_v1_to_v2(mt, phase1_legacy_spec()),
               "covariance hybrid")
  expect_error(sblr:::convert_blr_raw_v1_to_v2(
    list(schema = list(class = "unknown", version = 1L)),
    phase1_legacy_spec()), "No registered converter")

  wrong <- phase1_legacy_spec()
  wrong$model$family <- "bayesr"
  wrong$model$state_space <- c("component_0", "component_1")
  wrong$prior$component_multipliers <- c(component_0 = 0, component_1 = 1)
  expect_error(sblr:::convert_stblr_raw_v1_to_blr_raw_v2(
    phase1_legacy_bayesc_raw(), wrong), "model family")
  wrong <- phase1_legacy_spec()
  wrong$data$operator_resources$bed_shared$operator_type <- "bed"
  expect_error(sblr:::convert_stblr_raw_v1_to_blr_raw_v2(
    phase1_legacy_bayesc_raw(), wrong), "raw backend")
  wrong <- phase1_legacy_spec()
  wrong$model$probability_policy <- "fixed_marker"
  expect_error(sblr:::convert_stblr_raw_v1_to_blr_raw_v2(
    phase1_legacy_bayesc_raw(), wrong), "probability policy")
})

test_that("raw v2 rejects ambiguous names and unknown versions", {
  fixture <- phase1_fixture_environment()
  raw <- fixture$blr_raw_fixture("single_trait")
  raw$posterior$pi <- 0.5
  expect_error(sblr:::validate_blr_raw_v2(raw), "ambiguous probability")
  raw <- fixture$blr_raw_fixture("single_trait")
  raw$schema$version <- 99L
  expect_error(sblr:::validate_blr_raw_v2(raw), "schema\\$version")
  raw <- fixture$blr_raw_fixture("single_trait")
  raw$schema$compatibility_id <- "unknown-v2"
  expect_error(sblr:::validate_blr_raw_v2(raw), "compatibility_id")
})

test_that("cached provenance never invents Git metadata", {
  provenance <- sblr:::.blr_cached_provenance(refresh = TRUE)
  expect_true(is.character(provenance$package_version))
  expect_true(is.null(provenance$git_sha) ||
                grepl("^[0-9a-fA-F]{7,40}$", provenance$git_sha))
  expect_true("git_sha" %in% names(provenance))
  expect_true("dirty_build" %in% names(provenance))
  expect_true(is.null(provenance$dirty_build) ||
                is.logical(provenance$dirty_build))
  expect_identical(sblr:::.blr_cached_provenance(), provenance)
})

test_that("migrated public wrappers reject abbreviated explicit formals", {
  expect_error(stblr_csr_annot(
    NULL, annotations = matrix(0, 1L, 1L), ni = 2L),
    "exact argument names.*ni")
  expect_error(stblr_block_eigen(NULL, NULL, NULL, nb = 2L),
               "exact argument names.*nb")
  expect_error(stblr_csr(NULL, nch = 2L), "exact argument names.*nch")
  expect_error(stblr_bed(NULL, NULL, nth = 2L),
               "exact argument names.*nth")
  expect_error(make_summary_stats(NULL, NULL, nt = 2L),
               "exact argument names.*nt")
  expect_error(make_sparse_ld(NULL, nth = 2L),
               "exact argument names.*nth")
  expect_error(sparseLD_stream_CSR(NULL, 1L, list(), "x", nthread = 2L),
               "exact argument names.*nthread")
  expect_silent(sblr:::.blr_validate_exact_public_call(
    quote(stblr_csr(stats = NULL, nit = 2L)), stblr_csr, "stblr_csr()"))
})

test_that("resolved spec validates alleles regimes policies and priors", {
  fixture <- phase1_fixture_environment()
  spec <- fixture$blr_spec_fixture("independent_traits", "fixed_full")

  bad <- spec
  bad$data$global_alleles <- bad$data$global_alleles[c(2, 1, 3), ]
  expect_error(sblr:::validate_blr_resolved_spec(bad), "declared marker order")
  bad <- spec
  bad$data$likelihood_regime <- "mystery"
  expect_error(sblr:::validate_blr_resolved_spec(bad), "likelihood_regime")
  bad <- spec
  bad$data$providers$p1$likelihood_regime <- "mystery"
  expect_error(sblr:::validate_blr_resolved_spec(bad), "likelihood_regime")
  bad <- spec
  bad$data$providers$p1$likelihood_regime <- "common_sample"
  expect_error(sblr:::validate_blr_resolved_spec(bad), "must match")
  bad <- spec
  bad$model$probability_policy <- "future_policy"
  expect_error(sblr:::validate_blr_resolved_spec(bad), "probability_policy")
  bad <- spec
  bad$model$effect_storage_convention <- 1
  expect_error(sblr:::validate_blr_resolved_spec(bad),
               "effect_storage_convention")
  bad <- spec
  bad$prior$probability$alpha <- -1
  expect_error(sblr:::validate_blr_resolved_spec(bad),
               "probability.*alpha.*positive")
  bad <- spec
  bad$prior$scalar_variance$shape <- 0
  expect_error(sblr:::validate_blr_resolved_spec(bad),
               "scalar_variance.*positive")
  bad <- spec
  bad$prior$residual_covariance$fixed_value[1, 2] <- 2
  bad$prior$residual_covariance$fixed_value[2, 1] <- 2
  expect_error(sblr:::validate_blr_resolved_spec(bad),
               "positive-definite")
})

test_that("memory limit follows the frozen NULL nonnegative or Inf rule", {
  fixture <- phase1_fixture_environment()
  spec <- fixture$blr_spec_fixture("single_trait", "fixed_full")
  for (value in list(NULL, 0, 1024, Inf)) {
    candidate <- spec
    candidate$compute["memory_limit_bytes"] <- list(value)
    expect_true(sblr:::validate_blr_resolved_spec(candidate))
  }
  for (value in list(-1, -Inf, NA_real_, NaN, c(1, 2), "1")) {
    candidate <- spec
    candidate$compute$memory_limit_bytes <- value
    expect_error(sblr:::validate_blr_resolved_spec(candidate),
                 "memory_limit_bytes")
  }
})

test_that("raw v2 axes are bound to declared scientific IDs", {
  fixture <- phase1_fixture_environment()
  raw <- fixture$blr_raw_fixture("single_trait")
  mutate_axis <- function(x, axis, ids) {
    dimnames(x)[[match(axis, attr(x, "dim_axis_names"))]] <- ids
    x
  }
  bad <- raw
  bad$posterior$traitwise_state_probabilities <- mutate_axis(
    bad$posterior$traitwise_state_probabilities, "state",
    c("active", "null"))
  expect_error(sblr:::validate_blr_raw_v2(bad), "state")
  bad <- raw
  bad$posterior$traitwise_state_probabilities <- mutate_axis(
    bad$posterior$traitwise_state_probabilities, "state",
    c("null", "null"))
  expect_error(sblr:::validate_blr_raw_v2(bad), "state")
  bad <- raw
  bad$posterior$traitwise_state_probabilities <- fixture$.blr_probability_fixture(
    list(marker = raw$input$data$global_markers,
         trait = raw$input$data$trait_ids,
         state = c("null", "active", "other")))
  expect_error(sblr:::validate_blr_raw_v2(bad), "state")
  bad <- raw
  dimnames(bad$posterior$traitwise_state_probabilities)[[3L]] <- NULL
  expect_error(sblr:::validate_blr_raw_v2(bad), "state|stable dimnames")

  joint <- fixture$blr_raw_fixture("joint_multitrait")
  bad <- joint
  bad$posterior$joint_state_probabilities <- mutate_axis(
    bad$posterior$joint_state_probabilities, "joint_state",
    rev(joint$input$model$state_space))
  expect_error(sblr:::validate_blr_raw_v2(bad), "joint_state")
  bad <- joint
  bad$posterior$activity_pattern_probabilities <- mutate_axis(
    bad$posterior$activity_pattern_probabilities, "activity_pattern",
    c("00", "10", "01", "01"))
  expect_error(sblr:::validate_blr_raw_v2(bad), "activity_pattern")

  component <- raw
  component$input$model$family <- "bayesr"
  component$input$model$state_space <- c("component_0", "component_1")
  component$input$prior$component_multipliers <-
    c(component_0 = 0, component_1 = 1)
  component$model[names(component$input$model)] <- component$input$model
  dimnames(component$posterior$traitwise_state_probabilities)[[3L]] <-
    component$input$model$state_space
  dimnames(component$draws$traitwise_probability_parameters)[[4L]] <-
    component$input$model$state_space
  component$posterior$traitwise_component_assignment_probabilities <-
    fixture$.blr_probability_fixture(list(
      marker = component$input$data$global_markers,
      trait = component$input$data$trait_ids,
      component = component$input$model$state_space))
  expect_true(sblr:::validate_blr_raw_v2(component))
  bad <- component
  bad$posterior$traitwise_component_assignment_probabilities <- mutate_axis(
    bad$posterior$traitwise_component_assignment_probabilities, "component",
    rev(component$input$model$state_space))
  expect_error(sblr:::validate_blr_raw_v2(bad), "component")

  regional <- raw
  regional$input$data$statistical_regions <- stats::setNames(
    c("r1", "r1", "r2"), regional$input$data$global_markers)
  regional$input$model$marker_covariance_policy <- "regional"
  regional$model[names(regional$input$model)] <- regional$input$model
  regional$draws$regional_marker_covariance <- fixture$.blr_array(
    c(regional$input$mcmc$retained_draws, 1L, 2L, 1L, 1L),
    list(draw = paste0("draw", seq_len(regional$input$mcmc$retained_draws)),
         chain = "chain1", region = c("r1", "r2"),
         trait_row = "trait1", trait_col = "trait1"), value = 1)
  expect_true(sblr:::validate_blr_raw_v2(regional))
  bad <- regional
  bad$draws$regional_marker_covariance <- mutate_axis(
    bad$draws$regional_marker_covariance, "region", c("r2", "r1"))
  expect_error(sblr:::validate_blr_raw_v2(bad), "region")
})

test_that("compatibility IDs never bypass scientific validation", {
  fixture <- phase1_fixture_environment()
  for (compatibility_id in c("phase0-v2", "phase1-r-v2")) {
    raw <- fixture$blr_raw_fixture("joint_multitrait")
    raw$schema$compatibility_id <- compatibility_id
    raw$draws$marker_covariance[] <- 0
    expect_error(sblr:::validate_blr_raw_v2(raw),
                 "positive[- ]definite")
  }
})

test_that("raw provenance uses required-present NULL and exact types", {
  fixture <- phase1_fixture_environment()
  raw <- fixture$blr_raw_fixture("single_trait")
  expect_true(all(c("git_sha", "dirty_build", "compiler", "timestamp") %in%
                    names(raw$provenance)))
  expect_null(raw$provenance$git_sha)
  expect_null(raw$provenance$dirty_build)
  expect_true(sblr:::validate_blr_raw_v2(raw))
  bad <- raw
  bad$provenance["git_sha"] <- list(list("abcdef0"))
  expect_error(sblr:::validate_blr_raw_v2(bad), "git_sha")
  bad <- raw
  bad$provenance$dirty_build <- "FALSE"
  expect_error(sblr:::validate_blr_raw_v2(bad), "dirty_build")
  bad <- raw
  bad$provenance$git_sha <- "not-a-sha"
  expect_error(sblr:::validate_blr_raw_v2(bad), "git_sha")
  bad <- raw
  bad$provenance["compiler"] <- NULL
  expect_error(sblr:::validate_blr_raw_v2(bad), "missing required")
})

test_that("raw v2 preserves a one-trait one-chain one-state axis", {
  fixture <- phase1_fixture_environment()
  raw <- fixture$blr_raw_fixture("single_trait")
  raw$input$model$state_space <- "null"
  raw$input$model$null_state_index <- 1L
  raw$model[names(raw$input$model)] <- raw$input$model
  raw$posterior$traitwise_state_probabilities <- fixture$.blr_probability_fixture(
    list(marker = raw$input$data$global_markers,
         trait = raw$input$data$trait_ids, state = "null"))
  raw$draws$traitwise_probability_parameters <- fixture$.blr_probability_fixture(
    list(draw = paste0("draw", seq_len(raw$input$mcmc$retained_draws)),
         chain = "chain1", trait = raw$input$data$trait_ids, state = "null"))
  expect_true(sblr:::validate_blr_raw_v2(raw))
  expect_identical(dim(raw$posterior$traitwise_state_probabilities),
                   c(3L, 1L, 1L))
})
