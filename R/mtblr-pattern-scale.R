.blr_phase7_scale_spec <- function(
    method, marker_ids, component_scales = NULL,
    initial_scale_probability = NULL, scale_dirichlet_prior = NULL,
    marker_multipliers = NULL) {
  method <- .blr_character_scalar(method, "method", c("bayesc", "bayesr"))
  marker_ids <- .blr_ids(marker_ids, "marker_ids")
  supplied <- list(component_scales, initial_scale_probability,
                   scale_dirichlet_prior, marker_multipliers)
  if (identical(method, "bayesc")) {
    if (any(!vapply(supplied, is.null, logical(1)))) {
      stop("component_scales, scale probability controls, and marker_multipliers require method = 'bayesr'.",
           call. = FALSE)
    }
    return(list(enabled = FALSE, scales = numeric(), scale_ids = character(),
                initial = numeric(), prior = numeric(),
                marker_multipliers = numeric()))
  }
  if (is.null(component_scales)) component_scales <- c(0.01, 0.1, 1)
  if (!is.numeric(component_scales) || !length(component_scales) ||
      anyNA(component_scales) || any(!is.finite(component_scales)) ||
      any(component_scales <= 0) || anyDuplicated(component_scales) ||
      is.unsorted(component_scales, strictly = TRUE)) {
    stop("component_scales must be finite, positive, unique, and strictly increasing; the null state is not a component scale.",
         call. = FALSE)
  }
  scale_ids <- paste0("scale_", seq_along(component_scales))
  validate_probability <- function(x, what, simplex) {
    if (is.null(x)) x <- rep(1, length(component_scales))
    if (!is.numeric(x) || length(x) != length(component_scales) ||
        anyNA(x) || any(!is.finite(x)) || any(x <= 0)) {
      stop(what, " must contain one finite positive value per positive scale.",
           call. = FALSE)
    }
    if (!is.null(names(x)) && !identical(names(x), scale_ids)) {
      stop(what, " names must follow the declared positive-scale order.",
           call. = FALSE)
    }
    x <- as.numeric(x)
    if (simplex) x <- x / sum(x)
    stats::setNames(x, scale_ids)
  }
  initial <- validate_probability(
    initial_scale_probability, "initial_scale_probability", TRUE)
  prior <- validate_probability(
    scale_dirichlet_prior, "scale_dirichlet_prior", FALSE)
  if (is.null(marker_multipliers)) marker_multipliers <- rep(1, length(marker_ids))
  if (!is.numeric(marker_multipliers) ||
      length(marker_multipliers) != length(marker_ids) ||
      anyNA(marker_multipliers) || any(!is.finite(marker_multipliers)) ||
      any(marker_multipliers <= 0) ||
      (!is.null(names(marker_multipliers)) &&
       !identical(names(marker_multipliers), marker_ids))) {
    stop("marker_multipliers must contain one finite positive q_j value in global marker order.",
         call. = FALSE)
  }
  list(enabled = TRUE,
       scales = stats::setNames(as.numeric(component_scales), scale_ids),
       scale_ids = scale_ids, initial = initial, prior = prior,
       marker_multipliers = stats::setNames(
         as.numeric(marker_multipliers), marker_ids))
}

.blr_phase7_enrich_spec <- function(spec, scale) {
  if (!isTRUE(scale$enabled)) return(spec)
  spec$model$family <- "bayesr"
  spec$model$effect_storage_convention <- "scaled_latent"
  spec$model$marker_scale_policy <- "component_marker"
  spec$prior$probability$scale_ids <- scale$scale_ids
  spec$prior$probability$scale_dirichlet <- scale$prior
  spec$prior$probability$initial_scale_probability <- scale$initial
  spec$prior$probability$scale_sampled <- length(scale$scales) > 1L
  spec$prior$component_multipliers <- scale$scales
  spec$prior$marker_multipliers <- scale$marker_multipliers
  spec$mcmc$update_flags$scale_probability <- length(scale$scales) > 1L
  spec$output$retained_parameters <- c(
    spec$output$retained_parameters, "component_assignments",
    "scale_probability_parameters")
  spec$output$state_draw_policy <- "activity_pattern_and_positive_scale_draws"
  spec
}

.blr_phase7_memory_adjust <- function(memory, marker_count, trait_count,
                                      scale_count, chains, retained_draws,
                                      convergence_count, keep_traces,
                                      memory_limit_bytes, enforce = TRUE,
                                      activity_patterns = NULL,
                                      concurrent_chains = chains,
                                      workspace_enabled = scale_count > 0L) {
  if (!scale_count) return(memory)
  valid_limit <- is.null(memory_limit_bytes) || (
    is.numeric(memory_limit_bytes) && length(memory_limit_bytes) == 1L &&
    !is.na(memory_limit_bytes) && !is.nan(memory_limit_bytes) &&
    ((is.finite(memory_limit_bytes) && memory_limit_bytes >= 0) ||
     (is.infinite(memory_limit_bytes) && memory_limit_bytes > 0)))
  if (!valid_limit) {
    stop("memory_limit_bytes must be NULL, one finite nonnegative value, or positive Inf.",
         call. = FALSE)
  }
  memory$limit_bytes <- memory_limit_bytes
  dimensions <- c(
    marker_count = marker_count, trait_count = trait_count,
    scale_count = scale_count, chains = chains,
    concurrent_chains = concurrent_chains, retained_draws = retained_draws,
    convergence_count = convergence_count)
  if (!is.numeric(dimensions) || anyNA(dimensions) ||
      any(!is.finite(dimensions)) || any(dimensions < 0) ||
      any(dimensions != floor(dimensions)) || trait_count < 2L ||
      scale_count < 1L || chains < 1L || concurrent_chains < 1L ||
      concurrent_chains > chains) {
    stop("Phase 7 memory dimensions must be finite whole values with T >= 2, K >= 1, and valid chain concurrency.",
         call. = FALSE)
  }
  if (is.null(activity_patterns)) {
    activity_patterns <- .blr_phase4a_patterns(
      paste0("trait", seq_len(trait_count)))
  }
  if (!is.matrix(activity_patterns) || ncol(activity_patterns) != trait_count ||
      anyNA(activity_patterns) ||
      any(!activity_patterns %in% c(0L, 1L))) {
    stop("Phase 7 memory preflight requires the declared binary activity-pattern matrix.",
         call. = FALSE)
  }
  active_sizes <- rowSums(activity_patterns)
  if (sum(active_sizes == 0L) != 1L) {
    stop("Phase 7 memory preflight requires exactly one null activity pattern.",
         call. = FALSE)
  }
  pattern_count <- nrow(activity_patterns)
  context <- function(component, state_count = "unavailable") {
    paste0("component '", component, "' before sampling for T = ",
           trait_count, ", K = ", scale_count, ", pattern count = ",
           pattern_count, ", state count = ", state_count)
  }
  product <- function(..., component, state_count = "unavailable") {
    tryCatch(
      .blr_phase5a_checked_product(
        c(...), component, trait_count, pattern_count),
      error = function(e) stop(
        "Phase 7 memory preflight overflowed ",
        context(component, state_count), ".", call. = FALSE))
  }
  sum_checked <- function(..., component, state_count = "unavailable") {
    tryCatch(
      .blr_phase5a_checked_sum(
        c(...), component, trait_count, pattern_count),
      error = function(e) stop(
        "Phase 7 memory preflight overflowed ",
        context(component, state_count), ".", call. = FALSE))
  }
  bytes <- function(count, width, component, state_count = "unavailable") {
    product(count, width, component = component, state_count = state_count)
  }

  nonnull_state_count <- product(
    sum(active_sizes > 0L), scale_count,
    component = "pattern_scale_nonnull_state_count")
  state_count <- sum_checked(
    nonnull_state_count, 1,
    component = "pattern_scale_state_count",
    state_count = nonnull_state_count)
  active_moment_elements <- sum_checked(
    active_sizes[active_sizes > 0L] + active_sizes[active_sizes > 0L]^2,
    component = "pattern_scale_active_moment_elements",
    state_count = state_count)

  # The native kernel retains two scalar probability tables, pattern and
  # scale indices, and one arma::vec plus one arma::mat object per state.
  # A fixed 256-byte conservative slot for each Armadillo object is shared
  # with the native guard; numeric payloads are counted separately.
  workspace_chains <- if (isTRUE(workspace_enabled)) concurrent_chains else 0
  concurrent_states <- product(
    workspace_chains, state_count,
    component = "pattern_scale_concurrent_states",
    state_count = state_count)
  workspace <- c(
    phase7_pattern_scale_state_tables = bytes(
      concurrent_states, 28, "pattern_scale_state_tables", state_count),
    phase7_pattern_scale_active_containers = bytes(
      concurrent_states, 512, "pattern_scale_active_containers", state_count),
    phase7_pattern_scale_active_numeric = bytes(product(
      workspace_chains, scale_count, active_moment_elements,
      component = "pattern_scale_active_numeric_elements",
      state_count = state_count), 8,
      "pattern_scale_active_numeric", state_count),
    phase7_pattern_scale_kernel_containers = bytes(
      workspace_chains, 256, "pattern_scale_kernel_containers", state_count))

  component <- c(
    phase7_base_effect_state = bytes(product(
      chains, marker_count, trait_count,
      component = "phase7_base_effect_count", state_count = state_count),
      8, "phase7_base_effect_state", state_count),
    phase7_component_state = bytes(product(
      chains, marker_count, component = "phase7_component_state_count",
      state_count = state_count), 4, "phase7_component_state", state_count),
    phase7_final_component_state = bytes(product(
      chains, marker_count, component = "phase7_final_component_state_count",
      state_count = state_count), 4, "phase7_final_component_state", state_count),
    phase7_scale_probability_state = bytes(product(
      chains, scale_count, component = "phase7_scale_probability_count",
      state_count = state_count), 8, "phase7_scale_probability_state", state_count),
    phase7_final_scale_probability = bytes(product(
      chains, scale_count, component = "phase7_final_scale_probability_count",
      state_count = state_count), 8, "phase7_final_scale_probability", state_count),
    phase7_retained_component_assignments = bytes(product(
      chains, retained_draws, marker_count,
      component = "phase7_retained_component_count", state_count = state_count),
      4, "phase7_retained_component_assignments", state_count),
    phase7_retained_scale_parameters = bytes(product(
      chains, retained_draws, scale_count,
      component = "phase7_retained_scale_count", state_count = state_count),
      8, "phase7_retained_scale_parameters", state_count),
    phase7_convergence_scale_parameters = if (keep_traces) bytes(product(
      chains, convergence_count, scale_count,
      component = "phase7_convergence_scale_count", state_count = state_count),
      8, "phase7_convergence_scale_parameters", state_count) else 0,
    phase7_marker_component_summary = bytes(product(
      marker_count, scale_count,
      component = "phase7_marker_component_summary_count",
      state_count = state_count), 8,
      "phase7_marker_component_summary", state_count),
    phase7_compact_scale_diagnostics = bytes(product(
      chains, sum_checked(scale_count, 1,
        component = "phase7_compact_scale_diagnostic_width",
        state_count = state_count),
      component = "phase7_compact_scale_diagnostic_count",
      state_count = state_count), 8,
      "phase7_compact_scale_diagnostics", state_count),
    phase7_native_raw_component_copy = sum_checked(
      bytes(product(chains, retained_draws, marker_count,
        component = "phase7_native_component_copy_count",
        state_count = state_count), 4,
        "phase7_native_component_copy", state_count),
      bytes(product(marker_count, scale_count,
        component = "phase7_raw_component_copy_count",
        state_count = state_count), 8,
        "phase7_raw_component_copy", state_count),
      component = "phase7_native_raw_component_copy",
      state_count = state_count))

  base_components <- memory$components[
    !names(memory$components) %in% c(
      "pattern_conditional_workspace", "allocator_and_validation_headroom")]
  components <- c(base_components, workspace, component)
  subtotal <- sum_checked(
    unname(components), component = "phase7_memory_subtotal",
    state_count = state_count)
  headroom <- sum_checked(
    ceiling(subtotal / 10), 1024^2,
    component = "allocator_and_validation_headroom",
    state_count = state_count)
  memory$components <- c(
    components, allocator_and_validation_headroom = headroom)
  memory$estimated_peak_incremental_bytes <- sum_checked(
    unname(memory$components), component = "phase7_memory_estimate",
    state_count = state_count)
  memory$contract <- "phase7_peak_incremental_fit_allocation_v1"
  memory$dimensions$positive_scales <- scale_count
  memory$dimensions$pattern_scale_nonnull_states <- nonnull_state_count
  memory$dimensions$pattern_scale_states <- state_count
  memory$dimensions$pattern_scale_workspace_chains <- workspace_chains
  memory$output_policy$marker_component_probabilities <-
    "marginal_marker_by_positive_scale"
  memory$limit_exceeded <- !is.null(memory$limit_bytes) &&
    is.finite(memory$limit_bytes) &&
    memory$estimated_peak_incremental_bytes > memory$limit_bytes
  if (enforce && memory$limit_exceeded) {
    largest <- head(sort(memory$components, decreasing = TRUE), 5L)
    detail <- paste0(names(largest), "=", format(
      largest, scientific = FALSE, trim = TRUE), collapse = ", ")
    stop("Phase 7 memory preflight failed before provider construction and sampling: estimated ",
         format(memory$estimated_peak_incremental_bytes, scientific = FALSE),
         " incremental bytes, exceeding limit ",
         format(memory$limit_bytes, scientific = FALSE),
         "; M = ", marker_count, ", T = ", trait_count,
         ", K = ", scale_count, ", non-null states = ",
         nonnull_state_count, ", state count = ", state_count,
         ", chains = ", chains, ", simultaneous workspaces = ",
         workspace_chains, ". Largest components: ", detail,
         ". Marker-by-pattern-by-scale output is not allocated.",
         call. = FALSE)
  }
  memory
}

.blr_phase7_enrich_raw <- function(raw, native, scale, keep_traces) {
  if (!isTRUE(scale$enabled)) return(raw)
  spec <- raw$input
  markers <- spec$data$global_markers
  chains <- paste0("chain", seq_len(spec$mcmc$chains))
  draws <- paste0("draw", seq_len(spec$mcmc$retained_draws))
  components <- scale$scale_ids
  assignment <- .blr_phase4a_bind_draws(
    native$chains, "component_assignments",
    list(draw = draws, chain = chains, marker = markers))
  storage.mode(assignment) <- "integer"
  omega <- .blr_phase4a_bind_draws(
    native$chains, "scale_probability_parameters",
    list(draw = draws, chain = chains, component = components))
  final_assignment <- .blr_phase4a_bind_final(
    native$chains, "final_component_assignments",
    list(chain = chains, marker = markers))
  storage.mode(final_assignment) <- "integer"
  final_omega <- .blr_phase4a_bind_final(
    native$chains, "final_scale_probability_parameters",
    list(chain = chains, component = components))
  marker_component <- matrix(
    0, length(markers), length(components),
    dimnames = list(markers, components))
  for (marker in seq_along(markers)) {
    selected <- as.integer(assignment[, , marker])
    marker_component[marker, ] <- tabulate(
      selected[selected >= 0L] + 1L, nbins = length(components)) /
      length(selected)
  }
  attr(marker_component, "dim_axis_names") <- c("marker", "component")
  omega_mean <- apply(omega, 3L, mean)
  omega_mean <- .blr_make_array(omega_mean, list(component = components))
  raw$input <- spec
  raw$model <- c(list(analysis_mode = "joint_multitrait"), spec$model)
  raw$posterior$joint_component_assignment_probabilities <- marker_component
  raw$posterior$joint_component_probability_parameter_mean <- omega_mean
  raw$draws$component_assignments <- assignment
  raw$draws$joint_component_probability_parameters <- omega
  raw$final$component_assignments <- final_assignment
  raw$final$joint_component_probability_parameters <- final_omega
  if (isTRUE(keep_traces) && !is.null(raw$draws$convergence)) {
    iterations <- paste0("iteration", seq_len(spec$mcmc$sampling_iterations))
    raw$draws$convergence$scale_probability_parameters <-
      .blr_phase4a_bind_draws(
        native$chains, "convergence_scale_probability_parameters",
        list(iteration = iterations, chain = chains,
             component = components))
  }
  raw$diagnostics$qualification$implementation <-
    sub("bayesc", "bayesr_pattern_scale",
        raw$diagnostics$qualification$implementation, fixed = TRUE)
  raw$diagnostics$qualification$component_scales <- scale$scales
  raw$diagnostics$qualification$scale_dirichlet_prior <- scale$prior
  raw$diagnostics$qualification$scale_occupancy_counts <-
    stats::setNames(lapply(native$chains, function(chain) stats::setNames(
      chain$scale_occupancy_counts, components)), chains)
  raw$diagnostics$qualification$scale_change_counts <- stats::setNames(
    vapply(native$chains, `[[`, numeric(1), "scale_change_count"), chains)
  raw$schema$source_schema$name <- sub(
    "bayesc", "bayesr_pattern_scale", raw$schema$source_schema$name,
    fixed = TRUE)
  validate_blr_raw_v2(raw)
  raw
}
