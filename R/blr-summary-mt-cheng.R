# Internal Phase 6A independent-summary Cheng MT-BayesC-Pi qualification.

.blr_phase6a_native_resource <- function(resource) {
  .blr_validate_operator_resource(resource)
  type <- resource$operator_type
  if (identical(type, "csr")) {
    csr <- .blr_operator_storage_payload(resource$storage)
    return(list(
      operator_type = type,
      row_ptr = as.numeric(csr$row_ptr),
      column_index = as.integer(csr$column_index),
      values = as.numeric(csr$values),
      diagonal = as.numeric(csr$diagonal)))
  }
  if (type %in% c("full_rank_block_eigen",
                  "retained_rank_block_eigen")) {
    blocks <- lapply(resource$block_eigen$blocks, function(block) list(
      marker_indices = as.integer(block$marker_indices),
      eigenvectors = unname(block$eigenvectors),
      eigenvalues = as.numeric(block$eigenvalues)))
    diagonal <- numeric(length(resource$marker_ids))
    for (block in resource$block_eigen$blocks) {
      diagonal[block$marker_indices] <- rowSums(
        sweep(block$eigenvectors^2, 2L, block$eigenvalues, `*`))
    }
    return(list(
      operator_type = type, blocks = blocks, diagonal = diagonal))
  }
  stop("Phase 6A supports only CSR and block-eigen operator resources.",
       call. = FALSE)
}

.blr_phase6a_validate_collection <- function(collection, trait_ids) {
  if (!inherits(collection, "blr_provider_collection_v1")) {
    stop("Phase 6A requires a Phase 2 provider collection.", call. = FALSE)
  }
  trait_ids <- .blr_ids(trait_ids, "Phase 6A trait IDs")
  if (length(trait_ids) < 2L || length(trait_ids) > 12L) {
    stop("Phase 6A complete activity patterns require T in [2, 12].",
         call. = FALSE)
  }
  if (!identical(collection$analysis_mode, "joint_multitrait") ||
      !identical(collection$likelihood_regime, "independent_summary")) {
    stop(paste0(
      "Phase 6A requires joint_multitrait analysis with declared ",
      "independent_summary providers."), call. = FALSE)
  }
  provider_ids <- sort(names(collection$providers), method = "radix")
  resource_ids <- sort(names(collection$operator_resources), method = "radix")
  collection$providers <- collection$providers[provider_ids]
  collection$operator_resources <- collection$operator_resources[resource_ids]
  .blr_validate_global_marker_map(collection$global_marker_map)
  if (any(vapply(collection$providers, function(provider) {
      !identical(provider$likelihood_regime, "independent_summary") ||
        length(provider$trait_ids) != 1L
    }, logical(1)))) {
    stop("Every Phase 6A provider must own exactly one trait.", call. = FALSE)
  }
  represented <- unique(vapply(
    collection$providers, `[[`, character(1), "trait_ids"))
  if (!setequal(represented, trait_ids) ||
      any(!represented %in% trait_ids)) {
    stop("Every and only declared Phase 6A trait must have a provider.",
         call. = FALSE)
  }
  effect_scales <- vapply(
    collection$providers, `[[`, character(1), "effect_scale")
  if (length(unique(effect_scales)) != 1L) {
    stop("Phase 6A providers must use one compatible effect scale.",
         call. = FALSE)
  }
  for (provider in collection$providers) {
    if (!is.null(provider$overlap_group)) {
      stop(paste0(
        "Phase 6A does not interpret declared sample overlap as independent; ",
        "overlap-aware summary likelihoods are deferred."), call. = FALSE)
    }
    if (!identical(provider$residual_contract,
                   "fixed_provider_residual_scale")) {
      stop(paste0(
        "Phase 6A providers require residual_contract = ",
        "'fixed_provider_residual_scale'."), call. = FALSE)
    }
    resource <- collection$operator_resources[[
      provider$operator_resource_id]]
    .blr_validate_phase2_provider_alignment(
      provider, resource, collection$global_marker_map)
    .blr_validate_provider_statistics(provider, resource)
    if (!resource$operator_type %in% c(
        "csr", "full_rank_block_eigen", "retained_rank_block_eigen")) {
      stop("Phase 6A providers require CSR or block-eigen resources.",
           call. = FALSE)
    }
    if (!identical(resource$operator_scale, "cross_product")) {
      stop("Phase 6A resources must declare cross_product operator scale.",
           call. = FALSE)
    }
    phi <- provider$sufficient_statistics$residual_scale %||% NULL
    if (!is.numeric(phi) || length(phi) != 1L || is.na(phi) ||
        !is.finite(phi) || phi <= 0) {
      stop("Every Phase 6A provider requires one finite positive residual_scale.",
           call. = FALSE)
    }
  }
  collection
}

.blr_phase6a_native_inputs <- function(collection, trait_ids) {
  resource_ids <- names(collection$operator_resources)
  resources <- lapply(collection$operator_resources,
                      .blr_phase6a_native_resource)
  providers <- lapply(collection$providers, function(provider) {
    score <- provider$sufficient_statistics$score
    list(
      provider_id = provider$provider_id,
      trait_index = as.integer(match(provider$trait_ids, trait_ids) - 1L),
      resource_index = as.integer(match(
        provider$operator_resource_id, resource_ids) - 1L),
      local_to_global = as.integer(unname(provider$local_to_global) - 1L),
      score = as.numeric(score[, 1L]),
      residual_scale = as.numeric(
        provider$sufficient_statistics$residual_scale),
      sample_size = as.numeric(provider$sample_size[[1L]]))
  })
  list(resources = unname(resources), providers = unname(providers))
}

.blr_phase6a_memory_estimate <- function(
    collection, trait_count, chains, retained_draws, convergence_count,
    keep_traces, memory_limit_bytes = 256 * 1024^2, enforce = TRUE) {
  marker_count <- length(collection$global_marker_map$marker_ids)
  pattern_count <- bitwShiftL(1L, as.integer(trait_count))
  base <- .blr_phase5a_memory_estimate(
    marker_count = marker_count, trait_count = trait_count,
    observation_count = 0, chains = chains,
    retained_draws = retained_draws,
    convergence_count = convergence_count,
    sampled_residual = FALSE, keep_traces = keep_traces,
    source_sample_count = 0, selected_rows_used = FALSE,
    memory_limit_bytes = Inf, enforce = FALSE)
  components <- base$components[
    names(base$components) != "allocator_and_validation_headroom"]
  components[c("packed_bed_owner", "packed_bed_source_row_buffer")] <- 0
  product <- function(..., component) {
    .blr_phase5a_checked_product(
      c(...), component, trait_count, pattern_count)
  }
  sum_checked <- function(..., component) {
    .blr_phase5a_checked_sum(
      c(...), component, trait_count, pattern_count)
  }
  bytes <- function(count, width, component) {
    product(count, width, component = component)
  }
  resource_storage <- vapply(
    collection$operator_resources, function(resource) {
      if (identical(resource$operator_type, "csr")) {
        csr <- .blr_operator_storage_payload(resource$storage)
        return(sum_checked(
          bytes(length(csr$row_ptr), 8, "summary_csr_row_ptr"),
          bytes(length(csr$column_index), 12,
                "summary_csr_columns_values"),
          bytes(length(csr$diagonal), 8, "summary_csr_diagonal"),
          component = "summary_csr_resource"))
      }
      block_bytes <- vapply(resource$block_eigen$blocks, function(block) {
        sum_checked(
          bytes(length(block$marker_indices), 4,
                "summary_eigen_marker_indices"),
          bytes(length(block$eigenvectors), 8,
                "summary_eigen_factor"),
          bytes(length(block$eigenvalues), 8,
                "summary_eigen_values"),
          component = "summary_eigen_block")
      }, numeric(1))
      sum_checked(
        sum_checked(block_bytes, component = "summary_eigen_blocks"),
        bytes(length(resource$marker_ids), 16,
              "summary_eigen_diagonal_mapping"),
        component = "summary_eigen_resource")
    }, numeric(1))
  local_counts <- vapply(
    collection$providers, function(provider) length(provider$local_to_global),
    numeric(1))
  provider_count <- length(local_counts)
  components <- c(components,
    summary_operator_resources = sum_checked(
      resource_storage, component = "summary_operator_resources"),
    summary_provider_scores_maps_metadata = sum_checked(
      bytes(sum(local_counts), 20,
            "summary_provider_score_map_reference"),
      bytes(provider_count, 128, "summary_provider_metadata"),
      component = "summary_provider_scores_maps_metadata"),
    summary_chain_provider_residual_scores = bytes(product(
      chains, sum(local_counts), component = "summary_chain_provider_residuals"),
      8, "summary_chain_provider_residual_bytes"),
    summary_marker_provider_references = bytes(
      sum(local_counts), 16, "summary_marker_provider_reference_bytes"),
    summary_native_input_conversion_peak = sum_checked(
      sum_checked(resource_storage, component = "summary_resource_copy_peak"),
      bytes(sum(local_counts), 20, "summary_provider_copy_peak"),
      component = "summary_native_input_conversion_peak"))
  subtotal <- .blr_phase5a_checked_sum(
    unname(components), "phase6a_memory_subtotal", trait_count, pattern_count)
  headroom <- sum_checked(
    ceiling(subtotal / 10), 1024^2,
    component = "allocator_and_validation_headroom")
  components <- c(components, allocator_and_validation_headroom = headroom)
  estimate <- .blr_phase5a_checked_sum(
    unname(components), "phase6a_memory_estimate", trait_count, pattern_count)
  valid_limit <- is.null(memory_limit_bytes) || (
    is.numeric(memory_limit_bytes) && length(memory_limit_bytes) == 1L &&
    !is.na(memory_limit_bytes) && !is.nan(memory_limit_bytes) &&
    ((is.finite(memory_limit_bytes) && memory_limit_bytes >= 0) ||
       (is.infinite(memory_limit_bytes) && memory_limit_bytes > 0)))
  if (!valid_limit) {
    stop("memory_limit_bytes must be NULL, nonnegative, or positive Inf.",
         call. = FALSE)
  }
  exceeded <- !is.null(memory_limit_bytes) && is.finite(memory_limit_bytes) &&
    estimate > memory_limit_bytes
  if (isTRUE(enforce) && exceeded) {
    largest <- head(sort(components, decreasing = TRUE), 5L)
    stop(
      "Phase 6A memory preflight failed before native allocation: estimated ",
      format(estimate, scientific = FALSE), " bytes exceed limit ",
      format(memory_limit_bytes, scientific = FALSE), "; M = ", marker_count,
      ", T = ", trait_count, ", K = ", pattern_count,
      ", providers = ", provider_count, ", chains = ", chains,
      ", retained = ", retained_draws, ", convergence = ",
      convergence_count, ". Largest components: ",
      paste0(names(largest), "=", format(largest, scientific = FALSE),
             collapse = ", "), ".", call. = FALSE)
  }
  list(
    contract = "phase6a_summary_peak_incremental_allocation_v1",
    scope = "estimated_peak_incremental_bytes_for_this_fit",
    estimated_peak_incremental_bytes = estimate,
    limit_bytes = memory_limit_bytes, limit_exceeded = exceeded,
    dimensions = list(
      markers = marker_count, traits = trait_count, patterns = pattern_count,
      providers = provider_count, resources = length(resource_storage),
      chains = chains, retained_draws = retained_draws,
      convergence_iterations = convergence_count),
    output_policy = list(
      effect_draws = "full_qualification_draws",
      states = "joint_and_traitwise_activity",
      predictions = "unavailable_independent_summary",
      marker_pattern_probabilities = "mandatory_marker_by_pattern",
      convergence_traces_retained = keep_traces),
    components = components)
}

.blr_phase6a_resolved_spec <- function(
    collection, trait_ids, patterns, initial_marker_covariance,
    initial_probability, dirichlet_prior, prior_df, prior_scale,
    update_marker_covariance, update_probability,
    burn_in_iterations, sampling_iterations, thin_interval,
    chains, cores, seed, chain_seeds, keep_traces, memory_estimate) {
  retained <- .blr_retention_plan(
    burn_in_iterations, sampling_iterations, thin_interval,
    contract_version = 1L, retained_requested = TRUE)
  task_seeds <- .blr_task_seeds_v1(
    seed, "joint_multitrait", trait_ids, chains, chain_seeds)
  markers <- collection$global_marker_map$marker_ids
  new_blr_resolved_spec(
    schema = list(
      name = "blr_resolved_spec", version = 1L,
      compatibility_id = paste0(
        "phase6a-summary-cheng-qualification;seed=unified_fnv_splitmix_v1;",
        "retention=postburn_divisible_v1"),
      seed_contract_version = 1L, retention_contract_version = 1L,
      dimension_contract_version = 1L),
    data = list(
      analysis_mode = "joint_multitrait", trait_ids = trait_ids,
      global_markers = markers,
      global_alleles = collection$global_marker_map$alleles,
      operator_resources = collection$operator_resources,
      providers = collection$providers,
      provider_maps = lapply(collection$providers, `[[`, "local_to_global"),
      likelihood_regime = "independent_summary",
      statistical_regions = NULL),
    model = list(
      family = "bayesc", state_space = rownames(patterns),
      null_state_index = 1L, effect_storage_convention = "base_latent",
      probability_policy = "joint_activity_dirichlet",
      marker_scale_policy = "unit",
      marker_covariance_policy = "global_matrix",
      residual_policy = "diagonal", update_order_version = 1L),
    prior = list(
      probability = list(
        activity_patterns = patterns,
        activity_pattern_dirichlet = dirichlet_prior,
        initial_activity_pattern_probability = initial_probability,
        sampled = update_probability),
      component_multipliers = NULL,
      marker_multipliers = stats::setNames(rep(1, length(markers)), markers),
      scalar_variance = NULL,
      marker_covariance = list(
        degrees_of_freedom = if (update_marker_covariance) prior_df else NULL,
        scale = if (update_marker_covariance) prior_scale else NULL,
        fixed_value = if (update_marker_covariance) NULL else
          initial_marker_covariance,
        sampled = update_marker_covariance),
      residual_covariance = NULL, annotation = NULL),
    mcmc = list(
      burn_in_iterations = as.integer(burn_in_iterations),
      sampling_iterations = as.integer(sampling_iterations),
      thin_interval = as.integer(thin_interval),
      retained_draws = retained$retained_draws,
      retained_transition_indices = retained$post_burn,
      chains = as.integer(chains), seed = as.numeric(seed),
      task_seeds = task_seeds,
      update_flags = list(
        marker_covariance = update_marker_covariance,
        activity_pattern_probability = update_probability,
        residual_covariance = FALSE)),
    compute = list(
      execution_mode = if (cores == 1L) "serial" else "parallel",
      parallelization = if (cores == 1L) "none" else "chains",
      cores = as.integer(cores), scheduler_version = 1L,
      memory_limit_bytes = memory_estimate$limit_bytes,
      operator_numerical_controls = list(
        symmetry_tolerance = 1e-12,
        activity_pattern_count = nrow(patterns),
        transition_diagnostic_policy = "compact_occupancy_v1",
        dense_transition_diagnostics = FALSE,
        summary_operator_materialization = FALSE,
        fit_memory_estimate = memory_estimate)),
    output = list(
      posterior_summaries = TRUE,
      retained_parameters = c(
        "realised_effects", "latent_effects", "joint_states",
        "activity_pattern_parameters", "marker_covariance"),
      effect_draw_policy = "full_qualification_draws",
      state_draw_policy = "joint_activity_pattern_draws",
      convergence_policy = list(
        mode = "core", quantities = c(
          "marker_covariance", "activity_pattern_parameters",
          "active_marker_count")),
      derived_quantities = "none_independent_summary",
      preserve_chains = TRUE,
      memory_estimate_bytes =
        memory_estimate$estimated_peak_incremental_bytes))
}

.blr_phase6a_raw <- function(native, spec, initial_marker_covariance,
                             prior_df, prior_scale, dirichlet_prior,
                             keep_traces) {
  markers <- spec$data$global_markers
  traits <- spec$data$trait_ids
  patterns <- spec$prior$probability$activity_patterns
  state_ids <- rownames(patterns)
  draws <- paste0("draw", seq_len(spec$mcmc$retained_draws))
  chains <- paste0("chain", seq_len(spec$mcmc$chains))
  realised <- .blr_phase4a_bind_draws(native$chains, "realised_effects", list(
    draw = draws, chain = chains, marker = markers, trait = traits))
  latent <- .blr_phase4a_bind_draws(native$chains, "latent_effects", list(
    draw = draws, chain = chains, marker = markers, trait = traits))
  joint_state <- .blr_phase4a_bind_draws(native$chains, "joint_states", list(
    draw = draws, chain = chains, marker = markers))
  storage.mode(joint_state) <- "integer"
  marker_covariance <- .blr_phase4a_bind_draws(
    native$chains, "marker_covariance", list(
      draw = draws, chain = chains, trait_row = traits, trait_col = traits))
  pattern_parameter <- .blr_phase4a_bind_draws(
    native$chains, "activity_pattern_parameters", list(
      draw = draws, chain = chains, activity_pattern = state_ids))
  activity <- .blr_make_array(integer(length(realised)), list(
    draw = draws, chain = chains, marker = markers, trait = traits))
  for (trait in seq_along(traits)) {
    activity[, , , trait] <- patterns[joint_state + 1L, trait]
  }
  pattern_probability <- matrix(
    0, length(markers), length(state_ids),
    dimnames = list(markers, state_ids))
  for (marker in seq_along(markers)) {
    pattern_probability[marker, ] <- tabulate(
      as.integer(joint_state[, , marker]) + 1L,
      nbins = length(state_ids)) / (length(draws) * length(chains))
  }
  attr(pattern_probability, "dim_axis_names") <-
    c("marker", "activity_pattern")
  joint_probability <- pattern_probability
  attr(joint_probability, "dim_axis_names") <- c("marker", "joint_state")
  pips <- pattern_probability %*% patterns
  dimnames(pips) <- list(markers, traits)
  attr(pips, "dim_axis_names") <- c("marker", "trait")
  all_active <- which(rowSums(patterns) == length(traits))
  pleiotropic <- .blr_make_array(
    pattern_probability[, all_active], list(marker = markers))
  effect_mean <- apply(realised, c(3L, 4L), mean)
  dimnames(effect_mean) <- list(markers, traits)
  attr(effect_mean, "dim_axis_names") <- c("marker", "trait")
  covariance_mean <- apply(marker_covariance, c(3L, 4L), mean)
  dimnames(covariance_mean) <- list(traits, traits)
  attr(covariance_mean, "dim_axis_names") <- c("trait_row", "trait_col")

  final_realised <- .blr_phase4a_bind_final(
    native$chains, "final_realised_effects",
    list(chain = chains, marker = markers, trait = traits))
  final_latent <- .blr_phase4a_bind_final(
    native$chains, "final_latent_effects",
    list(chain = chains, marker = markers, trait = traits))
  final_state <- .blr_phase4a_bind_final(
    native$chains, "final_joint_states",
    list(chain = chains, marker = markers))
  storage.mode(final_state) <- "integer"
  final_covariance <- .blr_phase4a_bind_final(
    native$chains, "final_marker_covariance",
    list(chain = chains, trait_row = traits, trait_col = traits))
  final_probability <- .blr_phase4a_bind_final(
    native$chains, "final_activity_pattern_parameters",
    list(chain = chains, activity_pattern = state_ids))

  iterations <- paste0(
    "iteration", seq_len(spec$mcmc$sampling_iterations))
  convergence <- if (isTRUE(keep_traces)) list(
    transition_indices = seq_len(spec$mcmc$sampling_iterations),
    marker_covariance = .blr_phase4a_bind_draws(
      native$chains, "convergence_marker_covariance", list(
        iteration = iterations, chain = chains,
        trait_row = traits, trait_col = traits)),
    activity_pattern_parameters = .blr_phase4a_bind_draws(
      native$chains, "convergence_activity_pattern_parameters", list(
        iteration = iterations, chain = chains,
        activity_pattern = state_ids)),
    active_marker_count = .blr_phase4a_bind_draws(
      native$chains, "convergence_active_marker_count", list(
        iteration = iterations, chain = chains)),
    checkpoint = "completed_post_burn_iteration",
    thinning = "unthinned", rng_draws = 0L) else NULL
  covariance_updates <- lapply(native$chains, function(chain) list(
    degrees_of_freedom = chain$last_covariance_degrees_of_freedom,
    active_marker_count = chain$last_active_marker_count,
    statistic = chain$last_covariance_statistic,
    posterior_scale = chain$last_covariance_scale))
  provider_residuals <- stats::setNames(lapply(
    seq_along(native$chains), function(chain) stats::setNames(
      native$chains[[chain]]$final_provider_residual_score,
      names(spec$data$providers))), chains)
  provenance <- c(.blr_cached_provenance(), list(
    operator_resources = spec$data$operator_resources,
    marker_alignment = spec$data$provider_maps,
    seed_contract_version = spec$schema$seed_contract_version,
    task_seeds = spec$mcmc$task_seeds))
  approximations <- vapply(
    spec$data$operator_resources, `[[`, character(1), "approximation")
  warnings <- approximations[!approximations %in%
    c("exact_declared_operator", "exact_declared_block_diagonal_operator")]
  new_blr_raw_v2(
    model = c(list(analysis_mode = "joint_multitrait"), spec$model),
    input = spec,
    posterior = list(
      realised_effect_mean = effect_mean, latent_effect_mean = NULL,
      scaled_effect_mean = NULL, pips = pips,
      traitwise_state_probabilities = NULL,
      joint_state_probabilities = joint_probability,
      activity_pattern_probabilities = pattern_probability,
      traitwise_component_assignment_probabilities = NULL,
      joint_component_assignment_probabilities = NULL,
      marker_covariance_mean = covariance_mean,
      marker_variance_mean = NULL,
      residual_covariance_mean = NULL,
      residual_variance_mean = NULL, uncertainty = NULL,
      pleiotropic_probabilities = pleiotropic),
    draws = list(
      realised_effects = realised, latent_effects = latent,
      scaled_effects = NULL, independent_trait_states = NULL,
      joint_states = joint_state, traitwise_activity = activity,
      traitwise_probability_parameters = NULL,
      joint_probability_parameters = NULL,
      activity_pattern_parameters = pattern_parameter,
      traitwise_component_probability_parameters = NULL,
      joint_component_probability_parameters = NULL,
      marker_covariance = marker_covariance,
      residual_covariance = NULL, marker_variance = NULL,
      residual_variance = NULL, regional_marker_covariance = NULL,
      convergence = convergence),
    final = list(
      realised_effects = final_realised, latent_effects = final_latent,
      scaled_effects = NULL, independent_trait_states = NULL,
      joint_states = final_state,
      traitwise_probability_parameters = NULL,
      joint_probability_parameters = NULL,
      activity_pattern_parameters = final_probability,
      traitwise_component_probability_parameters = NULL,
      joint_component_probability_parameters = NULL,
      marker_covariance = final_covariance,
      residual_covariance = NULL, marker_variance = NULL,
      residual_variance = NULL, rng_continuation = NULL),
    derived = list(
      predictions = NULL, genetic_variance = NULL,
      genomic_covariance = NULL, operator_relative_quadratics = NULL,
      descriptive_bilinear_forms = NULL),
    diagnostics = list(
      convergence = convergence, acceptance = NULL, runtime = NULL,
      memory = spec$compute$operator_numerical_controls$fit_memory_estimate,
      workers = native$workers,
      numerical_safeguards = list(
        symmetry_tolerance = 1e-12, covariance_jitter = FALSE,
        log_weight_normalization = "log_sum_exp"),
      approximation_warnings = if (length(warnings)) warnings else NULL,
      qualification = list(
        status = "qualification_only",
        implementation = "phase6a_general_t_summary_cheng_mt_bayesc",
        likelihood_regime = "independent_summary",
        provider_residual_scales = stats::setNames(vapply(
          spec$data$providers, function(provider) {
            provider$sufficient_statistics$residual_scale
          }, numeric(1)), names(spec$data$providers)),
        operator_approximations = approximations,
        predictions_available = FALSE,
        update_order = c(
          "marker_sweep", "dirichlet", "marker_inverse_wishart",
          "convergence_capture", "retained_capture"),
        initial_marker_covariance = initial_marker_covariance,
        marker_covariance_prior = list(
          degrees_of_freedom = prior_df, scale = prior_scale),
        activity_pattern_dirichlet_prior = dirichlet_prior,
        covariance_updates = covariance_updates,
        final_provider_residual_scores = provider_residuals,
        transition_counts = NULL,
        pattern_occupancy_counts = stats::setNames(lapply(
          native$chains, function(chain) stats::setNames(
            chain$pattern_occupancy_counts, state_ids)), chains),
        pattern_change_counts = stats::setNames(vapply(
          native$chains, `[[`, numeric(1), "pattern_change_count"), chains),
        current_legacy_mt_route_used = FALSE)),
    provenance = provenance,
    compatibility_id = "phase1-r-v2",
    source_schema = list(
      name = "phase6a_summary_cheng_mt_bayesc_qualification",
      version = 1L),
    migration = list(
      status = "qualification_only", legacy_mt_conversion = FALSE))
}

.blr_cheng_mt_bayesc_summary_qualification <- function(
    collection, trait_ids, initial_marker_covariance,
    marker_covariance_prior_df, marker_covariance_prior_scale,
    initial_activity_pattern_probability = NULL,
    activity_pattern_dirichlet_prior = NULL,
    update_marker_covariance = TRUE,
    update_activity_pattern_probability = TRUE,
    burn_in_iterations = 100L, sampling_iterations = 200L,
    thin_interval = 1L, chains = 1L, cores = 1L, seed = 1,
    chain_seeds = NULL, keep_traces = TRUE,
    memory_limit_bytes = 256 * 1024^2,
    method = "bayesc", component_scales = NULL,
    initial_scale_probability = NULL, scale_dirichlet_prior = NULL,
    marker_multipliers = NULL) {
  trait_ids <- .blr_ids(trait_ids, "Phase 6A trait IDs")
  collection <- .blr_phase6a_validate_collection(collection, trait_ids)
  patterns <- .blr_phase4a_patterns(trait_ids)
  pattern_ids <- rownames(patterns)
  initial_marker_covariance <- .blr_phase4a_symmetric_spd(
    initial_marker_covariance, trait_ids, "initial_marker_covariance")
  marker_covariance_prior_scale <- .blr_phase4a_symmetric_spd(
    marker_covariance_prior_scale, trait_ids,
    "marker_covariance_prior_scale")
  if (!is.numeric(marker_covariance_prior_df) ||
      length(marker_covariance_prior_df) != 1L ||
      is.na(marker_covariance_prior_df) ||
      !is.finite(marker_covariance_prior_df) ||
      marker_covariance_prior_df <= length(trait_ids) - 1L) {
    stop("marker_covariance_prior_df must exceed T - 1.", call. = FALSE)
  }
  if (is.null(initial_activity_pattern_probability)) {
    initial_activity_pattern_probability <- rep(1 / nrow(patterns),
                                                nrow(patterns))
  }
  if (is.null(activity_pattern_dirichlet_prior)) {
    activity_pattern_dirichlet_prior <- rep(1, nrow(patterns))
  }
  initial_probability <- .blr_phase4a_probability(
    initial_activity_pattern_probability, pattern_ids,
    "initial_activity_pattern_probability", simplex = TRUE)
  dirichlet_prior <- .blr_phase4a_probability(
    activity_pattern_dirichlet_prior, pattern_ids,
    "activity_pattern_dirichlet_prior", simplex = FALSE)
  update_marker_covariance <- .blr_logical_scalar(
    update_marker_covariance, "update_marker_covariance")
  update_activity_pattern_probability <- .blr_logical_scalar(
    update_activity_pattern_probability,
    "update_activity_pattern_probability")
  keep_traces <- .blr_logical_scalar(keep_traces, "keep_traces")
  burn_in_iterations <- as.integer(.blr_scalar_whole(
    burn_in_iterations, "burn_in_iterations", 0, .Machine$integer.max))
  sampling_iterations <- as.integer(.blr_scalar_whole(
    sampling_iterations, "sampling_iterations", 1, .Machine$integer.max))
  thin_interval <- as.integer(.blr_scalar_whole(
    thin_interval, "thin_interval", 1, .Machine$integer.max))
  chains <- as.integer(.blr_scalar_whole(
    chains, "chains", 1, .Machine$integer.max))
  cores <- as.integer(.blr_scalar_whole(
    cores, "cores", 1, .Machine$integer.max))
  seed <- .blr_scalar_whole(seed, "seed", 0, 4294967295)
  .blr_chain_seed_base_uint32(chain_seeds, chains)
  retained <- .blr_retention_plan(
    burn_in_iterations, sampling_iterations, thin_interval,
    contract_version = 1L, retained_requested = TRUE)
  scale <- .blr_phase7_scale_spec(
    method, collection$global_marker_map$marker_ids, component_scales,
    initial_scale_probability, scale_dirichlet_prior, marker_multipliers)
  memory_estimate <- .blr_phase6a_memory_estimate(
    collection, length(trait_ids), chains, retained$retained_draws,
    sampling_iterations, keep_traces, memory_limit_bytes,
    enforce = !scale$enabled)
  memory_estimate <- .blr_phase7_memory_adjust(
    memory_estimate, length(collection$global_marker_map$marker_ids),
    length(trait_ids), length(scale$scales), chains,
    retained$retained_draws, sampling_iterations, keep_traces,
    memory_limit_bytes, enforce = TRUE, activity_patterns = patterns,
    concurrent_chains = min(chains, cores),
    workspace_enabled = scale$enabled && !(
      length(scale$scales) == 1L && scale$scales[[1L]] == 1 &&
      all(scale$marker_multipliers == 1)))
  spec <- .blr_phase6a_resolved_spec(
    collection, trait_ids, patterns, initial_marker_covariance,
    initial_probability, dirichlet_prior, marker_covariance_prior_df,
    marker_covariance_prior_scale, update_marker_covariance,
    update_activity_pattern_probability, burn_in_iterations,
    sampling_iterations, thin_interval, chains, cores, seed, chain_seeds,
    keep_traces, memory_estimate)
  native_input <- .blr_phase6a_native_inputs(collection, trait_ids)
  native <- mtblr_phase6a_summary_cheng_internal(
    operator_resources = native_input$resources,
    likelihood_providers = native_input$providers,
    global_marker_count = length(collection$global_marker_map$marker_ids),
    trait_count = length(trait_ids), activity_patterns = patterns,
    initial_marker_covariance = initial_marker_covariance,
    initial_activity_pattern_probability = as.numeric(initial_probability),
    activity_pattern_dirichlet_prior = as.numeric(dirichlet_prior),
    marker_covariance_prior_df = marker_covariance_prior_df,
    marker_covariance_prior_scale = marker_covariance_prior_scale,
    update_marker_covariance = update_marker_covariance,
    update_activity_pattern_probability =
      update_activity_pattern_probability,
    burn_in_iterations = burn_in_iterations,
    sampling_iterations = sampling_iterations,
    chains = chains, cores = cores,
    execution_contract = .blr_native_execution_contract(spec),
    native_memory_limit_bytes = if (is.null(memory_estimate$limit_bytes)) {
      Inf
    } else memory_estimate$limit_bytes,
    component_scales = unname(scale$scales),
    initial_scale_probability = unname(scale$initial),
    scale_dirichlet_prior = unname(scale$prior),
    marker_multipliers = unname(scale$marker_multipliers))
  raw <- .blr_phase6a_raw(
    native, spec, initial_marker_covariance,
    marker_covariance_prior_df, marker_covariance_prior_scale,
    dirichlet_prior, keep_traces)
  if (scale$enabled) {
    raw$input <- .blr_phase7_enrich_spec(raw$input, scale)
    raw <- .blr_phase7_enrich_raw(raw, native, scale, keep_traces)
  }
  validate_blr_raw_v2(raw)
  raw
}

.blr_pattern_scale_mt_bayesr_summary_qualification <- function(...) {
  .blr_cheng_mt_bayesc_summary_qualification(..., method = "bayesr")
}
