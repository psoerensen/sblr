.blr_phase4a_patterns <- function(trait_ids = c("trait1", "trait2"),
                                  maximum_traits = 12L) {
  trait_ids <- .blr_ids(trait_ids, "Cheng MT trait IDs")
  trait_count <- length(trait_ids)
  if (trait_count < 2L || trait_count > maximum_traits) {
    requested <- if (trait_count < 31L) 2^trait_count else Inf
    stop("Complete Cheng activity-pattern enumeration requires T in [2, ",
         maximum_traits, "]; requested T = ", trait_count,
         " (", format(requested, scientific = FALSE), " patterns).",
         call. = FALSE)
  }
  pattern_count <- bitwShiftL(1L, trait_count)
  state <- seq.int(0L, pattern_count - 1L)
  out <- vapply(seq_len(trait_count), function(trait) {
    bitwAnd(bitwShiftR(state, trait - 1L), 1L)
  }, integer(pattern_count))
  dimnames(out) <- list(
    apply(out, 1L, paste, collapse = "_"), trait_ids)
  storage.mode(out) <- "integer"
  out
}

.blr_phase5a_exact_integer_limit <- 2^53

.blr_phase5a_checked_product <- function(values, component, trait_count,
                                          pattern_count) {
  if (!is.numeric(values) || anyNA(values) || any(!is.finite(values)) ||
      any(values < 0) || any(values != floor(values)) ||
      any(values > .blr_phase5a_exact_integer_limit)) {
    stop("Phase 5A memory preflight cannot represent component '", component,
         "' exactly before sampling for T = ", trait_count, ", K = ",
         pattern_count, ".", call. = FALSE)
  }
  value <- 1
  for (factor in values) {
    if (value != 0 && factor > .blr_phase5a_exact_integer_limit / value) {
      stop("Phase 5A memory preflight overflowed component '", component,
           "' before sampling for T = ", trait_count, ", K = ",
           pattern_count, ".", call. = FALSE)
    }
    value <- value * factor
  }
  value
}

.blr_phase5a_checked_sum <- function(values, component, trait_count,
                                      pattern_count) {
  total <- 0
  for (value in values) {
    if (!is.numeric(value) || length(value) != 1L || is.na(value) ||
        !is.finite(value) || value < 0 || value != floor(value) ||
        value > .blr_phase5a_exact_integer_limit - total) {
      stop("Phase 5A memory preflight overflowed component '", component,
           "' before sampling for T = ", trait_count, ", K = ",
           pattern_count, ".", call. = FALSE)
    }
    total <- total + value
  }
  total
}

.blr_phase5a_packed_bed_allocation <- function(
    marker_count, selected_sample_count, source_sample_count,
    selected_rows_used, trait_count = 2L, pattern_count = 2^trait_count) {
  counts <- c(
    marker_count = marker_count,
    selected_sample_count = selected_sample_count,
    source_sample_count = source_sample_count)
  if (!is.numeric(counts) || anyNA(counts) || any(!is.finite(counts)) ||
      any(counts < 0) || any(counts != floor(counts)) ||
      any(counts > .blr_phase5a_exact_integer_limit) ||
      !is.logical(selected_rows_used) || length(selected_rows_used) != 1L ||
      is.na(selected_rows_used)) {
    stop("Phase 5A packed-BED dimensions and selection policy are invalid before sampling.",
         call. = FALSE)
  }
  if (isTRUE(selected_rows_used)) {
    if (source_sample_count < selected_sample_count) {
      stop("Phase 5A packed-BED source sample count cannot be smaller than the selected sample count.",
           call. = FALSE)
    }
  } else if (source_sample_count != selected_sample_count) {
    stop("Phase 5A all-sample packed-BED preparation requires equal source and selected sample counts.",
         call. = FALSE)
  }
  product <- function(..., component) {
    .blr_phase5a_checked_product(
      c(...), component, trait_count, pattern_count)
  }
  row_bytes <- floor(selected_sample_count / 4) +
    as.numeric(selected_sample_count %% 4 != 0)
  stride_blocks <- floor(row_bytes / 64) + as.numeric(row_bytes %% 64 != 0)
  aligned_stride <- product(
    stride_blocks, 64, component = "packed_bed_aligned_stride")
  owner_bytes <- product(
    marker_count, aligned_stride, component = "packed_bed_owner")
  source_row_bytes <- if (isTRUE(selected_rows_used)) {
    floor(source_sample_count / 4) +
      as.numeric(source_sample_count %% 4 != 0)
  } else 0
  list(
    selected_row_bytes = row_bytes,
    aligned_stride_bytes = aligned_stride,
    owner_bytes = owner_bytes,
    source_row_buffer_bytes = source_row_bytes,
    selected_rows_used = selected_rows_used)
}

.blr_phase5a_memory_estimate <- function(
    marker_count, trait_count, observation_count, chains, retained_draws,
    convergence_count, sampled_residual, keep_traces,
    source_sample_count = observation_count, selected_rows_used = FALSE,
    memory_limit_bytes = 256 * 1024^2, enforce = TRUE) {
  counts <- c(
    marker_count = marker_count, trait_count = trait_count,
    observation_count = observation_count,
    source_sample_count = source_sample_count, chains = chains,
    retained_draws = retained_draws, convergence_count = convergence_count)
  if (!is.numeric(counts) || anyNA(counts) || any(!is.finite(counts)) ||
      any(counts < 0) || any(counts != floor(counts)) ||
      trait_count < 2L || trait_count > 12L || chains < 1L) {
    stop("Phase 5A memory dimensions must be finite nonnegative whole values with T in [2, 12] and at least one chain.",
         call. = FALSE)
  }
  if (!is.logical(sampled_residual) || length(sampled_residual) != 1L ||
      is.na(sampled_residual) || !is.logical(keep_traces) ||
      length(keep_traces) != 1L || is.na(keep_traces) ||
      !is.logical(selected_rows_used) || length(selected_rows_used) != 1L ||
      is.na(selected_rows_used) ||
      !is.logical(enforce) || length(enforce) != 1L || is.na(enforce)) {
    stop("Phase 5A memory output-policy flags must be TRUE or FALSE.",
         call. = FALSE)
  }
  valid_limit <- is.null(memory_limit_bytes) || (
    is.numeric(memory_limit_bytes) && length(memory_limit_bytes) == 1L &&
    !is.na(memory_limit_bytes) && !is.nan(memory_limit_bytes) &&
    ((is.finite(memory_limit_bytes) && memory_limit_bytes >= 0) ||
       (is.infinite(memory_limit_bytes) && memory_limit_bytes > 0)))
  if (!valid_limit) {
    stop("memory_limit_bytes must be NULL, one finite nonnegative value, or positive Inf.",
         call. = FALSE)
  }

  pattern_count <- bitwShiftL(1L, as.integer(trait_count))
  packed_bed <- .blr_phase5a_packed_bed_allocation(
    marker_count = marker_count,
    selected_sample_count = observation_count,
    source_sample_count = source_sample_count,
    selected_rows_used = selected_rows_used,
    trait_count = trait_count, pattern_count = pattern_count)
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

  mt <- product(marker_count, trait_count, component = "marker_trait_count")
  nt <- product(observation_count, trait_count,
                component = "observation_trait_count")
  tt <- product(trait_count, trait_count,
                component = "trait_covariance_count")
  mk <- product(marker_count, pattern_count,
                component = "marker_pattern_count")
  cd <- product(chains, retained_draws,
                component = "chain_retained_count")
  cs <- product(chains, convergence_count,
                component = "chain_convergence_count")

  retained_effect_one <- bytes(product(cd, mt,
    component = "retained_effect_count"), 8,
    "retained_effect_bytes")
  retained_state <- bytes(product(cd, marker_count,
    component = "retained_state_count"), 4,
    "retained_state_bytes")
  retained_covariance <- bytes(product(cd, tt,
    1 + as.integer(sampled_residual),
    component = "retained_covariance_count"), 8,
    "retained_covariance_bytes")
  retained_probability <- bytes(product(cd, pattern_count,
    component = "retained_probability_count"), 8,
    "retained_probability_bytes")
  retained_prediction <- bytes(product(cd, nt,
    component = "retained_prediction_count"), 8,
    "retained_prediction_bytes")
  native_retained <- sum_checked(
    product(2, retained_effect_one,
      component = "two_native_retained_effect_arrays"),
    retained_state, retained_covariance,
    retained_probability, retained_prediction,
    component = "native_retained_output")

  convergence_covariance <- bytes(product(cs, tt,
    1 + as.integer(sampled_residual),
    component = "convergence_covariance_count"), 8,
    "convergence_covariance_bytes")
  convergence_probability <- bytes(product(cs, pattern_count,
    component = "convergence_probability_count"), 8,
    "convergence_probability_bytes")
  convergence_active <- bytes(cs, 4, "convergence_active_count_bytes")
  native_convergence <- sum_checked(
    convergence_covariance, convergence_probability, convergence_active,
    component = "native_convergence_output")

  final_effect_one <- bytes(product(chains, mt,
    component = "final_effect_count"), 8, "final_effect_bytes")
  final_state <- bytes(product(chains, marker_count,
    component = "final_state_count"), 4, "final_state_bytes")
  final_covariance <- bytes(product(chains, tt, 2,
    component = "final_covariance_count"), 8,
    "final_covariance_bytes")
  final_probability <- bytes(product(chains, pattern_count,
    component = "final_probability_count"), 8,
    "final_probability_bytes")
  final_prediction <- bytes(product(chains, nt,
    component = "final_prediction_count"), 8,
    "final_prediction_bytes")
  native_final <- sum_checked(
    product(2, final_effect_one,
      component = "two_native_final_effect_arrays"),
    final_state, final_covariance, final_probability,
    final_prediction, component = "native_final_output")

  components <- c(
    packed_bed_owner = packed_bed$owner_bytes,
    packed_bed_source_row_buffer = packed_bed$source_row_buffer_bytes,
    pattern_activity_metadata = sum_checked(
      bytes(product(pattern_count, trait_count,
        component = "pattern_activity_count"), 4,
        "pattern_activity_bytes"),
      bytes(pattern_count, 64, "pattern_id_bytes"),
      component = "pattern_activity_metadata"),
    pattern_parameters_and_counts = bytes(product(
      chains, pattern_count, 28,
      component = "pattern_parameters_and_counts"), 1,
      "pattern_parameters_and_counts_bytes"),
    pattern_conditional_workspace = bytes(product(
      chains, pattern_count,
      16 + 8 * trait_count + 8 * trait_count^2,
      component = "pattern_conditional_workspace"), 1,
      "pattern_conditional_workspace_bytes"),
    compact_transition_diagnostics = sum_checked(
      bytes(product(chains, pattern_count,
        component = "compact_occupancy_count"), 8,
        "compact_occupancy_bytes"),
      bytes(chains, 8, "compact_change_count_bytes"),
      component = "compact_transition_diagnostics"),
    dense_transition_diagnostics = 0,
    chain_mutable_effects = bytes(product(chains, mt, 2,
      component = "chain_mutable_effect_count"), 8,
      "chain_mutable_effect_bytes"),
    chain_mutable_states = bytes(product(chains, marker_count,
      component = "chain_mutable_state_count"), 4,
      "chain_mutable_state_bytes"),
    chain_residuals_and_marker_workspace = sum_checked(
      bytes(product(chains, nt,
        component = "chain_residual_count"), 8,
        "chain_residual_bytes"),
      bytes(product(chains, observation_count,
        component = "chain_marker_workspace_count"), 8,
        "chain_marker_workspace_bytes"),
      component = "chain_residuals_and_marker_workspace"),
    chain_iteration_indexing = bytes(product(
      chains, sum_checked(convergence_count, 1,
        component = "convergence_index_count"),
      component = "chain_iteration_index_count"), 4,
      "chain_iteration_index_bytes"),
    native_retained_output = native_retained,
    native_convergence_output = native_convergence,
    native_final_output = native_final,
    native_covariance_diagnostics = bytes(product(
      chains, tt, 4 + 2 * as.integer(sampled_residual),
      component = "native_covariance_diagnostic_count"), 8,
      "native_covariance_diagnostic_bytes"),
    native_marker_pattern_accumulator = 0,
    native_marker_pattern_return = 0,
    native_input_and_provider_copy = sum_checked(
      bytes(nt, 8, "native_phenotype_copy_bytes"),
      bytes(marker_count, 72, "native_marker_map_bytes"),
      component = "native_input_and_provider_copy"))

  rcpp_copy <- sum_checked(
    native_retained, native_convergence, native_final,
    components[["native_covariance_diagnostics"]],
    components[["compact_transition_diagnostics"]],
    component = "rcpp_native_return_copies")
  raw_retained <- sum_checked(
    native_retained,
    bytes(product(cd, mt,
      component = "retained_activity_count"), 4,
      "retained_activity_bytes"),
    component = "raw_retained_arrays")
  raw_convergence <- if (isTRUE(keep_traces)) native_convergence else 0
  raw_final <- native_final
  components <- c(components,
    rcpp_native_return_copies = rcpp_copy,
    raw_retained_arrays = raw_retained,
    raw_retained_binding_temporary = retained_effect_one,
    raw_convergence_arrays = raw_convergence,
    raw_convergence_binding_temporary = if (isTRUE(keep_traces)) {
      max(convergence_covariance, convergence_probability)
    } else 0,
    raw_final_arrays = raw_final,
    raw_final_binding_temporary = final_effect_one,
    marker_pattern_activity_probability_matrix = bytes(
      mk, 8, "marker_pattern_activity_probability_bytes"),
    marker_pattern_joint_state_probability_matrix = bytes(
      mk, 8, "marker_pattern_joint_state_probability_bytes"),
    marker_pattern_summary_outputs = sum_checked(
      bytes(mt, 8, "marker_pip_bytes"),
      bytes(marker_count, 8, "marker_pleiotropic_bytes"),
      component = "marker_pattern_summary_outputs"),
    posterior_effect_and_covariance_summaries = sum_checked(
      bytes(mt, 8, "posterior_effect_mean_bytes"),
      bytes(product(tt, 1 + as.integer(sampled_residual),
        component = "posterior_covariance_summary_count"), 8,
        "posterior_covariance_summary_bytes"),
      component = "posterior_effect_and_covariance_summaries"),
    resolved_spec_and_provider_metadata = sum_checked(
      bytes(marker_count, 96, "resolved_marker_metadata_bytes"),
      bytes(observation_count, 32, "resolved_observation_metadata_bytes"),
      bytes(nt, 8, "resolved_phenotype_bytes"),
      component = "resolved_spec_and_provider_metadata"))

  subtotal <- .blr_phase5a_checked_sum(
    unname(components), "estimated_peak_incremental_subtotal",
    trait_count, pattern_count)
  components <- c(components,
    allocator_and_validation_headroom = sum_checked(
      ceiling(subtotal / 10), 1024^2,
      component = "allocator_and_validation_headroom"))
  estimate <- .blr_phase5a_checked_sum(
    unname(components), "estimated_peak_incremental_bytes",
    trait_count, pattern_count)

  limit_exceeded <- !is.null(memory_limit_bytes) &&
    is.finite(memory_limit_bytes) && estimate > memory_limit_bytes
  if (isTRUE(enforce) && limit_exceeded) {
    largest <- head(sort(components, decreasing = TRUE), 5L)
    detail <- paste0(names(largest), "=", format(
      largest, scientific = FALSE, trim = TRUE), collapse = ", ")
    stop(
      "Phase 5A memory preflight failed before sampling: estimated peak ",
      "incremental bytes ", format(estimate, scientific = FALSE),
      " exceed configured limit ",
      format(memory_limit_bytes, scientific = FALSE), "; M = ",
      format(marker_count, scientific = FALSE, trim = TRUE), ", T = ",
      trait_count, ", K = ", pattern_count, ", chains = ",
      format(chains, scientific = FALSE, trim = TRUE), ", retained = ",
      format(retained_draws, scientific = FALSE, trim = TRUE),
      ", convergence = ",
      format(convergence_count, scientific = FALSE, trim = TRUE),
      ", selected N = ",
      format(observation_count, scientific = FALSE, trim = TRUE),
      ", source N = ",
      format(source_sample_count, scientific = FALSE, trim = TRUE),
      ", selected rows used = ", selected_rows_used,
      ", dense transition diagnostics requested = FALSE. Largest ",
      "components: ", detail, ". Complete enumeration and mandatory ",
      "markerwise pattern output scale as M * 2^T.", call. = FALSE)
  }

  list(
    contract = "phase5a_peak_incremental_fit_allocation_v1",
    scope = "estimated_peak_incremental_bytes_for_this_fit",
    estimated_peak_incremental_bytes = estimate,
    limit_bytes = memory_limit_bytes,
    limit_exceeded = limit_exceeded,
    dimensions = list(
      markers = as.numeric(marker_count), traits = as.integer(trait_count),
      patterns = as.integer(pattern_count),
      observations = as.numeric(observation_count),
      source_observations = as.numeric(source_sample_count),
      chains = as.numeric(chains), retained_draws = as.numeric(retained_draws),
      convergence_iterations = as.numeric(convergence_count)),
    packed_bed = packed_bed,
    output_policy = list(
      effect_draws = "full_qualification_draws",
      states = "joint_and_traitwise_activity",
      predictions = "retained_predictions",
      marker_pattern_probabilities = "mandatory_marker_by_pattern",
      convergence_traces_retained = keep_traces,
      sampled_residual_covariance = sampled_residual),
    dense_transition_diagnostics = FALSE,
    components = components)
}

.blr_phase4a_symmetric_spd <- function(x, trait_ids, what,
                                       tolerance = 1e-12) {
  if (!is.matrix(x) || !is.numeric(x) ||
      !identical(dim(x), rep.int(length(trait_ids), 2L)) ||
      any(!is.finite(x)) ||
      max(abs(x - t(x))) > tolerance) {
    stop(what, " must be a finite symmetric T x T matrix in declared trait order.",
         call. = FALSE)
  }
  if (!is.null(dimnames(x)) &&
      !identical(dimnames(x), list(trait_ids, trait_ids))) {
    stop(what, " dimnames must follow the declared trait order.",
         call. = FALSE)
  }
  value <- 0.5 * (x + t(x))
  if (inherits(try(chol(value), silent = TRUE), "try-error")) {
    stop(what, " must be strictly positive definite.", call. = FALSE)
  }
  dimnames(value) <- list(trait_ids, trait_ids)
  value
}

.blr_phase4a_probability <- function(x, pattern_ids, what,
                                     simplex = FALSE) {
  if (!is.numeric(x) || length(x) != length(pattern_ids) || anyNA(x) ||
      any(!is.finite(x)) || any(x <= 0)) {
    stop(what, " must contain one finite positive value per activity pattern.",
         call. = FALSE)
  }
  if (!is.null(names(x)) && !identical(names(x), pattern_ids)) {
    stop(what, " names must follow the canonical activity-pattern order.",
         call. = FALSE)
  }
  if (simplex && abs(sum(x) - 1) > 1e-10) {
    stop(what, " must sum to one.", call. = FALSE)
  }
  stats::setNames(as.numeric(x) / if (simplex) sum(x) else 1, pattern_ids)
}

.blr_phase4a_bed_collection <- function(dat, Glist, phenotype, trait_ids,
                                        residual_contract) {
  marker_ids <- as.character(dat$variable_names)
  marker_metadata <- .mtblr_bed_marker_metadata(dat, Glist)
  effect <- marker_metadata$effect_allele %||% rep(NA_character_, dat$m)
  other <- marker_metadata$other_allele %||% rep(NA_character_, dat$m)
  alleles <- data.frame(
    marker_id = marker_ids, effect = as.character(effect),
    other = as.character(other),
    coding = rep("effect_allele_count", dat$m),
    stringsAsFactors = FALSE)
  global_map <- .blr_new_global_marker_map(marker_ids, alleles)
  sample_ids <- rownames(phenotype)
  if (is.null(sample_ids)) {
    selected <- dat$rows %||% seq_len(dat$n)
    sample_ids <- if (!is.null(Glist$ids) && length(Glist$ids) == dat$n) {
      as.character(Glist$ids[selected])
    } else paste0("sample", seq_len(nrow(phenotype)))
  }
  rownames(phenotype) <- sample_ids
  selected_rows <- dat$rows %||% seq_len(dat$n)
  resource <- .blr_new_operator_resource(
    resource_id = "phase4a_common_bed", operator_type = "bed",
    marker_ids = marker_ids, alleles = alleles,
    genotype_coding = "effect_allele_count",
    centering = "allele_frequency_centered",
    standardization = "variance_standardized",
    operator_scale = "individual_standardized_genotypes",
    storage = .blr_new_operator_storage_ref(
      "packed_bed", list(
        bed_files = as.character(dat$bed_files),
        source_sample_count = as.integer(dat$n), sample_ids = sample_ids,
        selected_rows = as.integer(selected_rows),
        selected_columns = lapply(dat$cls, as.integer))),
    block_eigen = NULL, approximation = "exact_selected_genotypes",
    provenance = list(
      phase = "5a", ownership = "shared_immutable_packed_bed_reference"))
  provider <- .blr_new_likelihood_provider(
    provider_id = "phase4a_common_sample", trait_ids = trait_ids,
    operator_resource_id = resource$resource_id,
    local_to_global = stats::setNames(seq_along(marker_ids), marker_ids),
    sufficient_statistics = list(phenotype = phenotype),
    sample_size = stats::setNames(rep(nrow(phenotype), length(trait_ids)),
                                 trait_ids),
    likelihood_regime = "common_sample",
    residual_contract = residual_contract,
    population = NULL, effect_scale = "phenotype_native",
    overlap_group = NULL,
    provenance = list(phase = "5a", status = "qualification_only"))
  collection <- .blr_new_provider_collection(
    global_map, stats::setNames(list(resource), resource$resource_id),
    stats::setNames(list(provider), provider$provider_id),
    likelihood_regime = "common_sample",
    analysis_mode = "joint_multitrait")
  list(collection = collection, phenotype = phenotype,
       observation_ids = sample_ids)
}

.blr_phase4a_resolved_spec <- function(
    collection, observation_ids, trait_ids, marker_ids,
    residual_covariance_policy, fixed_residual_covariance,
    initial_residual_covariance, residual_covariance_prior_df,
    residual_covariance_prior_scale, initial_marker_covariance,
    initial_probability, dirichlet_prior, prior_df, prior_scale,
    update_marker_covariance, update_probability, activity_patterns,
    burn_in_iterations, sampling_iterations, thin_interval,
    chains, cores, seed, keep_traces, memory_estimate) {
  retained <- .blr_retention_plan(
    burn_in_iterations, sampling_iterations, thin_interval,
    contract_version = 1L, retained_requested = TRUE)
  task_seeds <- .blr_task_seeds_v1(
    seed, "joint_multitrait", trait_ids, chains)
  state_ids <- rownames(activity_patterns)
  data <- list(
    analysis_mode = "joint_multitrait", trait_ids = trait_ids,
    global_markers = marker_ids,
    global_alleles = collection$global_marker_map$alleles,
    operator_resources = collection$operator_resources,
    providers = collection$providers,
    provider_maps = lapply(collection$providers, `[[`, "local_to_global"),
    likelihood_regime = "common_sample", statistical_regions = NULL,
    observation_ids = observation_ids)
  new_blr_resolved_spec(
    schema = list(
      name = "blr_resolved_spec", version = 1L,
      compatibility_id = paste0(
        "phase5a-general-t-qualification;seed=unified_fnv_splitmix_v1;",
        "retention=postburn_divisible_v1"),
      seed_contract_version = 1L, retention_contract_version = 1L,
      dimension_contract_version = 1L),
    data = data,
    model = list(
      family = "bayesc", state_space = state_ids, null_state_index = 1L,
      effect_storage_convention = "base_latent",
      probability_policy = "joint_activity_dirichlet",
      marker_scale_policy = "unit",
      marker_covariance_policy = "global_matrix",
      residual_policy = residual_covariance_policy,
      update_order_version = if (identical(
        residual_covariance_policy, "sampled_full")) 2L else 1L),
    prior = list(
      probability = list(
        activity_patterns = activity_patterns,
        activity_pattern_dirichlet = dirichlet_prior,
        initial_activity_pattern_probability = initial_probability,
        sampled = update_probability),
      component_multipliers = NULL,
      marker_multipliers = stats::setNames(rep(1, length(marker_ids)),
                                            marker_ids),
      scalar_variance = NULL,
      marker_covariance = list(
        degrees_of_freedom = if (update_marker_covariance) prior_df else NULL,
        scale = if (update_marker_covariance) prior_scale else NULL,
        fixed_value = if (update_marker_covariance) NULL else
          initial_marker_covariance,
        sampled = update_marker_covariance),
      residual_covariance = if (identical(
        residual_covariance_policy, "sampled_full")) list(
          degrees_of_freedom = residual_covariance_prior_df,
          scale = residual_covariance_prior_scale,
          fixed_value = NULL, sampled = TRUE,
          initial_value = initial_residual_covariance,
          parameterization = "degrees_of_freedom_scale",
          proper = TRUE,
          finite_mean = residual_covariance_prior_df >
            length(trait_ids) + 1L) else list(
          degrees_of_freedom = NULL, scale = NULL,
          fixed_value = fixed_residual_covariance, sampled = FALSE),
      annotation = NULL),
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
        residual_covariance = identical(
          residual_covariance_policy, "sampled_full"))),
    compute = list(
      execution_mode = if (cores == 1L) "serial" else "parallel",
      parallelization = if (cores == 1L) "none" else "chains",
      cores = as.integer(cores), scheduler_version = 1L,
      memory_limit_bytes = memory_estimate$limit_bytes,
      operator_numerical_controls = list(
        symmetry_tolerance = 1e-12, genotype_materialization = FALSE,
        activity_pattern_count = memory_estimate$dimensions$patterns,
        dense_transition_diagnostics = FALSE,
        transition_diagnostic_policy = "compact_occupancy_v1",
        fit_memory_estimate = memory_estimate)),
    output = list(
      posterior_summaries = TRUE,
      retained_parameters = c(
        "realised_effects", "latent_effects", "joint_states",
        "activity_pattern_parameters", "marker_covariance", "predictions",
        if (identical(residual_covariance_policy, "sampled_full")) {
          "residual_covariance"
        } else character()),
      effect_draw_policy = "full_qualification_draws",
      state_draw_policy = "joint_activity_pattern_draws",
      convergence_policy = list(
        mode = "core", quantities = c(
          "marker_covariance", "activity_pattern_parameters",
          "active_marker_count",
          if (identical(residual_covariance_policy, "sampled_full")) {
            "residual_covariance"
          } else character())),
      derived_quantities = "retained_predictions",
      preserve_chains = TRUE,
      memory_estimate_bytes =
        memory_estimate$estimated_peak_incremental_bytes))
}

.blr_phase4a_bind_draws <- function(chains, field, axes) {
  out <- .blr_make_array(numeric(prod(lengths(axes))), axes)
  chain_axis <- match("chain", names(axes))
  for (chain in seq_along(chains)) {
    subscripts <- rep(list(TRUE), length(axes))
    subscripts[[chain_axis]] <- chain
    value <- chains[[chain]][[field]]
    target_dim <- lengths(axes)
    target_dim[[chain_axis]] <- 1L
    value <- array(value, dim = target_dim)
    out <- do.call(`[<-`, c(list(out), subscripts, list(value = value)))
  }
  attr(out, "dim_axis_names") <- names(axes)
  out
}

.blr_phase4a_bind_final <- function(chains, field, axes) {
  out <- .blr_make_array(numeric(prod(lengths(axes))), axes)
  chain_axis <- match("chain", names(axes))
  for (chain in seq_along(chains)) {
    subscripts <- rep(list(TRUE), length(axes))
    subscripts[[chain_axis]] <- chain
    target_dim <- lengths(axes)
    target_dim[[chain_axis]] <- 1L
    value <- array(chains[[chain]][[field]], dim = target_dim)
    out <- do.call(`[<-`, c(list(out), subscripts, list(value = value)))
  }
  attr(out, "dim_axis_names") <- names(axes)
  out
}

.blr_phase4a_raw <- function(native, spec, residual_covariance_policy,
                             fixed_residual_covariance,
                             initial_residual_covariance,
                             residual_covariance_prior_df,
                             residual_covariance_prior_scale,
                             initial_marker_covariance,
                             prior_df, prior_scale, dirichlet_prior,
                             keep_traces) {
  markers <- spec$data$global_markers
  traits <- spec$data$trait_ids
  observations <- spec$data$observation_ids
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
  sampled_residual <- identical(residual_covariance_policy, "sampled_full")
  residual_covariance <- if (sampled_residual) {
    .blr_phase4a_bind_draws(
      native$chains, "residual_covariance", list(
        draw = draws, chain = chains,
        trait_row = traits, trait_col = traits))
  } else NULL
  pattern_parameter <- .blr_phase4a_bind_draws(
    native$chains, "activity_pattern_parameters", list(
      draw = draws, chain = chains, activity_pattern = state_ids))
  predictions <- .blr_phase4a_bind_draws(native$chains, "predictions", list(
    draw = draws, chain = chains, observation = observations, trait = traits))
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
  pleiotropic <- .blr_make_array(pattern_probability[, all_active],
                                 list(marker = markers))
  effect_mean <- apply(realised, c(3L, 4L), mean)
  dimnames(effect_mean) <- list(markers, traits)
  attr(effect_mean, "dim_axis_names") <- c("marker", "trait")
  covariance_mean <- apply(marker_covariance, c(3L, 4L), mean)
  dimnames(covariance_mean) <- list(traits, traits)
  attr(covariance_mean, "dim_axis_names") <- c("trait_row", "trait_col")
  residual_covariance_mean <- if (sampled_residual) {
    value <- apply(residual_covariance, c(3L, 4L), mean)
    dimnames(value) <- list(traits, traits)
    attr(value, "dim_axis_names") <- c("trait_row", "trait_col")
    value
  } else NULL

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
  residual_chains <- if (sampled_residual) native$chains else
    lapply(seq_along(chains), function(index) list(
      final_residual_covariance = fixed_residual_covariance))
  final_residual <- .blr_phase4a_bind_final(
    residual_chains, "final_residual_covariance",
    list(chain = chains, trait_row = traits, trait_col = traits))

  iterations <- paste0("iteration", seq_len(spec$mcmc$sampling_iterations))
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
  if (isTRUE(keep_traces) && sampled_residual) {
    convergence <- append(convergence, list(
      residual_covariance = .blr_phase4a_bind_draws(
        native$chains, "convergence_residual_covariance", list(
          iteration = iterations, chain = chains,
          trait_row = traits, trait_col = traits))), after = 2L)
  }
  covariance_updates <- lapply(native$chains, function(chain) list(
    degrees_of_freedom = chain$last_covariance_degrees_of_freedom,
    active_marker_count = chain$last_active_marker_count,
    statistic = chain$last_covariance_statistic,
    posterior_scale = chain$last_covariance_scale))
  residual_covariance_updates <- if (sampled_residual) {
    lapply(native$chains, function(chain) list(
      degrees_of_freedom =
        chain$last_residual_covariance_degrees_of_freedom,
      observation_count = length(observations),
      statistic = chain$last_residual_covariance_statistic,
      posterior_scale = chain$last_residual_covariance_scale,
      update_count = chain$residual_covariance_update_count))
  } else NULL

  provenance <- c(.blr_cached_provenance(), list(
    operator_resources = spec$data$operator_resources,
    marker_alignment = spec$data$provider_maps,
    seed_contract_version = spec$schema$seed_contract_version,
    task_seeds = spec$mcmc$task_seeds))
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
      residual_covariance_mean = residual_covariance_mean,
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
      residual_covariance = residual_covariance, marker_variance = NULL,
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
      residual_covariance = final_residual,
      marker_variance = NULL, residual_variance = NULL,
      rng_continuation = NULL),
    derived = list(
      predictions = predictions, genetic_variance = NULL,
      genomic_covariance = NULL, operator_relative_quadratics = NULL,
      descriptive_bilinear_forms = NULL),
    diagnostics = list(
      convergence = convergence, acceptance = NULL, runtime = NULL,
      memory = spec$compute$operator_numerical_controls$fit_memory_estimate,
      workers = native$workers,
      numerical_safeguards = list(
        symmetry_tolerance = 1e-12, covariance_jitter = FALSE,
        log_weight_normalization = "log_sum_exp"),
      approximation_warnings = NULL,
      qualification = list(
        status = "qualification_only",
        implementation = "phase5a_general_t_cheng_mt_bayesc",
        trait_count = length(traits), pattern_count = length(state_ids),
        pattern_order = state_ids,
        update_order = if (sampled_residual) c(
          "marker_sweep", "dirichlet", "marker_inverse_wishart",
          "residual_inverse_wishart", "convergence_capture",
          "retained_capture") else
          c("marker_sweep", "dirichlet", "inverse_wishart"),
        residual_covariance_policy = residual_covariance_policy,
        fixed_residual_covariance = fixed_residual_covariance,
        initial_residual_covariance = initial_residual_covariance,
        residual_covariance_prior = if (sampled_residual) list(
          degrees_of_freedom = residual_covariance_prior_df,
          scale = residual_covariance_prior_scale,
          parameterization = "degrees_of_freedom_scale",
          proper = TRUE,
          finite_mean = residual_covariance_prior_df > length(traits) + 1L
        ) else NULL,
        initial_marker_covariance = initial_marker_covariance,
        marker_covariance_prior = list(
          degrees_of_freedom = prior_df, scale = prior_scale),
        activity_pattern_dirichlet_prior = dirichlet_prior,
        covariance_updates = covariance_updates,
        residual_covariance_updates = residual_covariance_updates,
        transition_counts = NULL,
        pattern_occupancy_counts = stats::setNames(
          lapply(native$chains, function(chain) {
            stats::setNames(chain$pattern_occupancy_counts, state_ids)
          }), chains),
        pattern_change_counts = stats::setNames(vapply(
          native$chains, `[[`, numeric(1), "pattern_change_count"), chains),
        genotype_contract = native$genotype_contract,
        current_legacy_mt_route_used = FALSE)),
    provenance = provenance,
    compatibility_id = "phase1-r-v2",
    source_schema = list(
      name = "phase5a_general_t_cheng_mt_bayesc_qualification",
      version = 1L),
    migration = list(
      status = "qualification_only", legacy_mt_conversion = FALSE))
}

.blr_phase4a_cheng_mt_bed <- function(
    y, Glist, fixed_residual_covariance, initial_marker_covariance,
    marker_covariance_prior_df, marker_covariance_prior_scale,
    initial_activity_pattern_probability = NULL,
    activity_pattern_dirichlet_prior = NULL,
    update_marker_covariance = TRUE,
    update_activity_pattern_probability = TRUE,
    burn_in_iterations = 100L, sampling_iterations = 200L,
    thin_interval = 1L, chains = 1L, cores = 1L, seed = 1,
    keep_traces = TRUE, chr = NULL, cls = NULL, rows = NULL,
    block_size = 1000L,
    residual_covariance_policy = "fixed_full",
    initial_residual_covariance = NULL,
    residual_covariance_prior_df = NULL,
    residual_covariance_prior_scale = NULL,
    memory_limit_bytes = 256 * 1024^2) {
  if (!is.list(Glist) || is.null(Glist$bedfiles)) {
    stop("Phase 4a requires one packed-BED Glist.", call. = FALSE)
  }
  dat <- .make_bed_marker_data(
    Glist = Glist, y = y, chr = chr, cls = cls,
    block_size = block_size, rows = rows)
  phenotype <- as.matrix(dat$y)
  if (!is.numeric(phenotype) || ncol(phenotype) < 2L ||
      nrow(phenotype) <= 1L || any(!is.finite(phenotype))) {
    stop("The Cheng MT phenotype must be a complete finite N x T matrix with T >= 2.",
         call. = FALSE)
  }
  trait_ids <- colnames(phenotype)
  if (is.null(trait_ids)) trait_ids <- paste0("trait", seq_len(ncol(phenotype)))
  trait_ids <- .blr_ids(trait_ids, "Phase 4a trait IDs")
  colnames(phenotype) <- trait_ids
  patterns <- .blr_phase4a_patterns(trait_ids)
  pattern_count <- nrow(patterns)
  if (!is.character(residual_covariance_policy) ||
      length(residual_covariance_policy) != 1L ||
      is.na(residual_covariance_policy) ||
      !residual_covariance_policy %in% c("fixed_full", "sampled_full")) {
    stop("residual_covariance_policy must be exactly 'fixed_full' or 'sampled_full'.",
         call. = FALSE)
  }
  sampled_residual <- identical(residual_covariance_policy, "sampled_full")
  if (sampled_residual) {
    if (!is.null(fixed_residual_covariance)) {
      stop("sampled_full residual covariance requires fixed_residual_covariance = NULL.",
           call. = FALSE)
    }
    initial_residual_covariance <- .blr_phase4a_symmetric_spd(
      initial_residual_covariance, trait_ids,
      "initial_residual_covariance")
    residual_covariance_prior_scale <- .blr_phase4a_symmetric_spd(
      residual_covariance_prior_scale, trait_ids,
      "residual_covariance_prior_scale")
    if (!is.numeric(residual_covariance_prior_df) ||
        length(residual_covariance_prior_df) != 1L ||
        is.na(residual_covariance_prior_df) ||
        !is.finite(residual_covariance_prior_df) ||
        residual_covariance_prior_df <= length(trait_ids) - 1L) {
      stop("residual_covariance_prior_df must exceed T - 1 for a proper inverse-Wishart prior.",
           call. = FALSE)
    }
  } else {
    if (!is.null(initial_residual_covariance) ||
        !is.null(residual_covariance_prior_df) ||
        !is.null(residual_covariance_prior_scale)) {
      stop("fixed_full residual covariance cannot carry sampled residual-covariance initial or prior fields.",
           call. = FALSE)
    }
    fixed_residual_covariance <- .blr_phase4a_symmetric_spd(
      fixed_residual_covariance, trait_ids, "fixed_residual_covariance")
    initial_residual_covariance <- fixed_residual_covariance
  }
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
    stop("marker_covariance_prior_df must exceed T - 1 for a proper inverse-Wishart prior.",
         call. = FALSE)
  }
  pattern_ids <- rownames(patterns)
  if (is.null(initial_activity_pattern_probability)) {
    initial_activity_pattern_probability <- rep(1 / pattern_count,
                                                pattern_count)
  }
  if (is.null(activity_pattern_dirichlet_prior)) {
    activity_pattern_dirichlet_prior <- rep(1, pattern_count)
  }
  initial_probability <- .blr_phase4a_probability(
    initial_activity_pattern_probability, pattern_ids,
    "initial_activity_pattern_probability", simplex = TRUE)
  dirichlet_prior <- .blr_phase4a_probability(
    activity_pattern_dirichlet_prior, pattern_ids,
    "activity_pattern_dirichlet_prior", simplex = FALSE)
  for (name in c("update_marker_covariance",
                 "update_activity_pattern_probability", "keep_traces")) {
    value <- get(name)
    if (!is.logical(value) || length(value) != 1L || is.na(value)) {
      stop(name, " must be TRUE or FALSE.", call. = FALSE)
    }
  }
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
  frequency <- unlist(dat$af, use.names = FALSE)
  if (length(frequency) != dat$m || any(!is.finite(frequency)) ||
      any(frequency <= 0 | frequency >= 1)) {
    stop("Phase 4a selected allele frequencies must lie inside (0, 1).",
         call. = FALSE)
  }

  retained <- .blr_retention_plan(
    burn_in_iterations, sampling_iterations, thin_interval,
    contract_version = 1L, retained_requested = TRUE)
  memory_estimate <- .blr_phase5a_memory_estimate(
    marker_count = dat$m, trait_count = length(trait_ids),
    observation_count = nrow(phenotype), chains = chains,
    retained_draws = retained$retained_draws,
    convergence_count = sampling_iterations,
    sampled_residual = sampled_residual, keep_traces = keep_traces,
    source_sample_count = dat$n,
    selected_rows_used = !is.null(dat$rows),
    memory_limit_bytes = memory_limit_bytes, enforce = TRUE)

  prepared <- .blr_phase4a_bed_collection(
    dat, Glist, phenotype, trait_ids,
    if (sampled_residual) {
      "sampled_full_inverse_wishart_residual_covariance"
    } else "fixed_full_residual_covariance")
  spec <- .blr_phase4a_resolved_spec(
    prepared$collection, prepared$observation_ids, trait_ids,
    dat$variable_names, residual_covariance_policy,
    fixed_residual_covariance, initial_residual_covariance,
    residual_covariance_prior_df, residual_covariance_prior_scale,
    initial_marker_covariance, initial_probability, dirichlet_prior,
    marker_covariance_prior_df, marker_covariance_prior_scale,
    update_marker_covariance, update_activity_pattern_probability, patterns,
    burn_in_iterations, sampling_iterations, thin_interval,
    chains, cores, seed, keep_traces, memory_estimate)
  execution_contract <- .blr_native_execution_contract(spec)
  native <- mtblr_phase4a_cheng_bed_internal(
    bed_files = as.character(dat$bed_files),
    source_sample_count = as.integer(dat$n),
    selected_columns = lapply(dat$cls, as.integer),
    selected_rows = if (is.null(dat$rows)) NULL else as.integer(dat$rows),
    allele_frequency = as.numeric(frequency), phenotype = phenotype,
    initial_residual_covariance = initial_residual_covariance,
    initial_marker_covariance = initial_marker_covariance,
    activity_patterns = patterns,
    initial_activity_pattern_probability = as.numeric(initial_probability),
    activity_pattern_dirichlet_prior = as.numeric(dirichlet_prior),
    marker_covariance_prior_df = marker_covariance_prior_df,
    marker_covariance_prior_scale = marker_covariance_prior_scale,
    update_marker_covariance = update_marker_covariance,
    update_activity_pattern_probability =
      update_activity_pattern_probability,
    update_residual_covariance = sampled_residual,
    residual_covariance_prior_df = if (sampled_residual) {
      residual_covariance_prior_df
    } else 0,
    residual_covariance_prior_scale = if (sampled_residual) {
      residual_covariance_prior_scale
    } else matrix(numeric(), 0L, 0L),
    burn_in_iterations = burn_in_iterations,
    sampling_iterations = sampling_iterations,
    chains = chains, cores = cores,
    execution_contract = execution_contract,
    native_memory_limit_bytes = if (is.null(memory_estimate$limit_bytes)) {
      Inf
    } else memory_estimate$limit_bytes)
  .blr_phase4a_raw(
    native, spec, residual_covariance_policy, fixed_residual_covariance,
    initial_residual_covariance, residual_covariance_prior_df,
    residual_covariance_prior_scale, initial_marker_covariance,
    marker_covariance_prior_df, marker_covariance_prior_scale,
    dirichlet_prior, keep_traces)
}

# Maintained qualification spelling for the generalized Cheng implementation.
# The historical Phase 4 name remains as the implementation boundary so the
# approved two-trait checkpoint can be compared without a duplicate sampler.
.blr_cheng_mt_bayesc_bed_qualification <- function(...) {
  .blr_phase4a_cheng_mt_bed(...)
}
