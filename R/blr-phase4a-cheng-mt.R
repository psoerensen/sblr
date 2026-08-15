.blr_phase4a_patterns <- function() {
  out <- rbind(`0_0` = c(0L, 0L), `1_0` = c(1L, 0L),
               `0_1` = c(0L, 1L), `1_1` = c(1L, 1L))
  colnames(out) <- c("trait1", "trait2")
  out
}

.blr_phase4a_symmetric_spd <- function(x, trait_ids, what,
                                       tolerance = 1e-12) {
  if (!is.matrix(x) || !is.numeric(x) ||
      !identical(dim(x), c(2L, 2L)) || any(!is.finite(x)) ||
      max(abs(x - t(x))) > tolerance) {
    stop(what, " must be a finite symmetric 2 x 2 matrix.", call. = FALSE)
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
    stop(what, " must contain four finite positive values.", call. = FALSE)
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

.blr_phase4a_bed_collection <- function(dat, Glist, phenotype, trait_ids) {
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
      phase = "4a", ownership = "shared_immutable_packed_bed_reference"))
  provider <- .blr_new_likelihood_provider(
    provider_id = "phase4a_common_sample", trait_ids = trait_ids,
    operator_resource_id = resource$resource_id,
    local_to_global = stats::setNames(seq_along(marker_ids), marker_ids),
    sufficient_statistics = list(phenotype = phenotype),
    sample_size = stats::setNames(rep(nrow(phenotype), 2L), trait_ids),
    likelihood_regime = "common_sample",
    residual_contract = "fixed_full_residual_covariance",
    population = NULL, effect_scale = "phenotype_native",
    overlap_group = NULL,
    provenance = list(phase = "4a", status = "qualification_only"))
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
    fixed_residual_covariance, initial_marker_covariance,
    initial_probability, dirichlet_prior, prior_df, prior_scale,
    update_marker_covariance, update_probability,
    burn_in_iterations, sampling_iterations, thin_interval,
    chains, cores, seed, keep_traces) {
  retained <- .blr_retention_plan(
    burn_in_iterations, sampling_iterations, thin_interval,
    contract_version = 1L, retained_requested = TRUE)
  task_seeds <- .blr_task_seeds_v1(
    seed, "joint_multitrait", trait_ids, chains)
  state_ids <- rownames(.blr_phase4a_patterns())
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
        "phase4a-qualification;seed=unified_fnv_splitmix_v1;",
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
      residual_policy = "fixed_full", update_order_version = 1L),
    prior = list(
      probability = list(
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
      residual_covariance = list(
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
        residual_covariance = FALSE)),
    compute = list(
      execution_mode = if (cores == 1L) "serial" else "parallel",
      parallelization = if (cores == 1L) "none" else "chains",
      cores = as.integer(cores), scheduler_version = 1L,
      memory_limit_bytes = NULL,
      operator_numerical_controls = list(
        symmetry_tolerance = 1e-12, genotype_materialization = FALSE)),
    output = list(
      posterior_summaries = TRUE,
      retained_parameters = c(
        "realised_effects", "latent_effects", "joint_states",
        "activity_pattern_parameters", "marker_covariance", "predictions"),
      effect_draw_policy = "full_qualification_draws",
      state_draw_policy = "joint_activity_pattern_draws",
      convergence_policy = list(
        mode = "core", quantities = c(
          "marker_covariance", "activity_pattern_parameters",
          "active_marker_count")),
      derived_quantities = "retained_predictions",
      preserve_chains = TRUE, memory_estimate_bytes = NA_real_))
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

.blr_phase4a_raw <- function(native, spec, fixed_residual_covariance,
                             initial_marker_covariance,
                             prior_df, prior_scale, dirichlet_prior,
                             keep_traces) {
  markers <- spec$data$global_markers
  traits <- spec$data$trait_ids
  observations <- spec$data$observation_ids
  patterns <- .blr_phase4a_patterns()
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
  predictions <- .blr_phase4a_bind_draws(native$chains, "predictions", list(
    draw = draws, chain = chains, observation = observations, trait = traits))
  activity <- .blr_make_array(integer(length(realised)), list(
    draw = draws, chain = chains, marker = markers, trait = traits))
  for (trait in seq_len(2L)) {
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
  pleiotropic <- .blr_make_array(pattern_probability[, "1_1"],
                                 list(marker = markers))
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
  residual_chains <- lapply(seq_along(chains), function(index) list(
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
  covariance_updates <- lapply(native$chains, function(chain) list(
    degrees_of_freedom = chain$last_covariance_degrees_of_freedom,
    active_marker_count = chain$last_active_marker_count,
    statistic = chain$last_covariance_statistic,
    posterior_scale = chain$last_covariance_scale))

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
      residual_covariance = final_residual,
      marker_variance = NULL, residual_variance = NULL,
      rng_continuation = NULL),
    derived = list(
      predictions = predictions, genetic_variance = NULL,
      genomic_covariance = NULL, operator_relative_quadratics = NULL,
      descriptive_bilinear_forms = NULL),
    diagnostics = list(
      convergence = convergence, acceptance = NULL, runtime = NULL,
      memory = NULL, workers = native$workers,
      numerical_safeguards = list(
        symmetry_tolerance = 1e-12, covariance_jitter = FALSE,
        log_weight_normalization = "log_sum_exp"),
      approximation_warnings = NULL,
      qualification = list(
        status = "qualification_only", pattern_order = state_ids,
        update_order = c("marker_sweep", "dirichlet", "inverse_wishart"),
        fixed_residual_covariance = fixed_residual_covariance,
        initial_marker_covariance = initial_marker_covariance,
        marker_covariance_prior = list(
          degrees_of_freedom = prior_df, scale = prior_scale),
        activity_pattern_dirichlet_prior = dirichlet_prior,
        covariance_updates = covariance_updates,
        transition_counts = lapply(native$chains, `[[`, "transition_counts"),
        genotype_contract = native$genotype_contract,
        current_legacy_mt_route_used = FALSE)),
    provenance = provenance,
    compatibility_id = "phase1-r-v2",
    source_schema = list(
      name = "phase4a_cheng_mt_bayesc_qualification", version = 1L),
    migration = list(
      status = "qualification_only", legacy_mt_conversion = FALSE))
}

.blr_phase4a_cheng_mt_bed <- function(
    y, Glist, fixed_residual_covariance, initial_marker_covariance,
    marker_covariance_prior_df, marker_covariance_prior_scale,
    initial_activity_pattern_probability = rep(0.25, 4L),
    activity_pattern_dirichlet_prior = rep(1, 4L),
    update_marker_covariance = TRUE,
    update_activity_pattern_probability = TRUE,
    burn_in_iterations = 100L, sampling_iterations = 200L,
    thin_interval = 1L, chains = 1L, cores = 1L, seed = 1,
    keep_traces = TRUE, chr = NULL, cls = NULL, rows = NULL,
    block_size = 1000L) {
  if (!is.list(Glist) || is.null(Glist$bedfiles)) {
    stop("Phase 4a requires one packed-BED Glist.", call. = FALSE)
  }
  dat <- .make_bed_marker_data(
    Glist = Glist, y = y, chr = chr, cls = cls,
    block_size = block_size, rows = rows)
  phenotype <- as.matrix(dat$y)
  if (!is.numeric(phenotype) || ncol(phenotype) != 2L ||
      nrow(phenotype) <= 1L || any(!is.finite(phenotype))) {
    stop("Phase 4a phenotype must be a complete finite N x 2 matrix.",
         call. = FALSE)
  }
  trait_ids <- colnames(phenotype)
  if (is.null(trait_ids)) trait_ids <- c("trait1", "trait2")
  trait_ids <- .blr_ids(trait_ids, "Phase 4a trait IDs")
  colnames(phenotype) <- trait_ids
  fixed_residual_covariance <- .blr_phase4a_symmetric_spd(
    fixed_residual_covariance, trait_ids, "fixed_residual_covariance")
  initial_marker_covariance <- .blr_phase4a_symmetric_spd(
    initial_marker_covariance, trait_ids, "initial_marker_covariance")
  marker_covariance_prior_scale <- .blr_phase4a_symmetric_spd(
    marker_covariance_prior_scale, trait_ids,
    "marker_covariance_prior_scale")
  if (!is.numeric(marker_covariance_prior_df) ||
      length(marker_covariance_prior_df) != 1L ||
      is.na(marker_covariance_prior_df) ||
      !is.finite(marker_covariance_prior_df) ||
      marker_covariance_prior_df <= 1) {
    stop("marker_covariance_prior_df must exceed T - 1 for a proper inverse-Wishart prior.",
         call. = FALSE)
  }
  pattern_ids <- rownames(.blr_phase4a_patterns())
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

  prepared <- .blr_phase4a_bed_collection(
    dat, Glist, phenotype, trait_ids)
  spec <- .blr_phase4a_resolved_spec(
    prepared$collection, prepared$observation_ids, trait_ids,
    dat$variable_names, fixed_residual_covariance,
    initial_marker_covariance, initial_probability, dirichlet_prior,
    marker_covariance_prior_df, marker_covariance_prior_scale,
    update_marker_covariance, update_activity_pattern_probability,
    burn_in_iterations, sampling_iterations, thin_interval,
    chains, cores, seed, keep_traces)
  execution_contract <- .blr_native_execution_contract(spec)
  native <- mtblr_phase4a_cheng_bed_internal(
    bed_files = as.character(dat$bed_files),
    source_sample_count = as.integer(dat$n),
    selected_columns = lapply(dat$cls, as.integer),
    selected_rows = if (is.null(dat$rows)) NULL else as.integer(dat$rows),
    allele_frequency = as.numeric(frequency), phenotype = phenotype,
    fixed_residual_covariance = fixed_residual_covariance,
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
    execution_contract = execution_contract)
  .blr_phase4a_raw(
    native, spec, fixed_residual_covariance, initial_marker_covariance,
    marker_covariance_prior_df, marker_covariance_prior_scale,
    dirichlet_prior, keep_traces)
}
