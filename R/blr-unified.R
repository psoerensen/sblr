.blr_scalar_integer <- function(value, name, minimum = 1L) {
  if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
      value != floor(value) || value < minimum || value > .Machine$integer.max) {
    stop(name, " must be one integer-compatible value >= ", minimum, ".",
         call. = FALSE)
  }
  as.integer(value)
}

.blr_logical_scalar <- function(value, name) {
  if (!is.logical(value) || length(value) != 1L || is.na(value)) {
    stop(name, " must be TRUE or FALSE.", call. = FALSE)
  }
  value
}

.blr_chain_controls <- function(nit = 1000, nburn = 500, nthin = 1,
                                seed = 1, nchains = 1L, ncores = 1L,
                                chain_seeds = NULL, keep_chains = FALSE) {
  nit <- .blr_scalar_integer(nit, "nit")
  nburn <- .blr_scalar_integer(nburn, "nburn", 0L)
  nthin <- .blr_scalar_integer(nthin, "nthin")
  seed <- .blr_scalar_integer(seed, "seed")
  nchains <- .blr_scalar_integer(nchains, "nchains")
  ncores <- .blr_scalar_integer(ncores, "ncores")
  keep_chains <- .blr_logical_scalar(keep_chains, "keep_chains")
  requested <- chain_seeds
  if (is.null(chain_seeds)) {
    native <- integer()
  } else {
    if (!is.numeric(chain_seeds) || length(chain_seeds) != nchains ||
        anyNA(chain_seeds) || any(!is.finite(chain_seeds)) ||
        any(chain_seeds != floor(chain_seeds)) ||
        any(chain_seeds < -2147483648) ||
        any(chain_seeds > .Machine$integer.max)) {
      stop("chain_seeds must contain exactly one signed 32-bit integer per chain.",
           call. = FALSE)
    }
    native <- unname(chain_seeds)
  }
  list(nit = nit, nburn = nburn, nthin = nthin, seed = seed,
       nchains = nchains, ncores = ncores,
       chain_seeds_requested = requested, chain_seeds_native = native,
       keep_chains = keep_chains)
}

.blr_convergence_controls <- function(convergence = c("auto", "none", "core", "extended"),
                                      convergence_control = NULL,
                                      nchains = 1L) {
  convergence <- match.arg(convergence)
  defaults <- list(
    warn = TRUE, rhat_threshold = 1.01,
    ess_per_chain_threshold = 100,
    mcse_mean_over_sd_threshold = 0.05, keep_traces = FALSE,
    extended_groups = NULL, selected_markers = NULL,
    selected_marker_quantities = c("b", "d"),
    full_probability_states = FALSE, aggregate_component_states = FALSE,
    max_trace_gb = 1,
    allow_large_traces = FALSE)
  if (is.null(convergence_control)) {
    resolved <- defaults
  } else {
    if (!is.list(convergence_control) || is.data.frame(convergence_control) ||
        is.null(names(convergence_control)) || any(!nzchar(names(convergence_control))) ||
        anyDuplicated(names(convergence_control)) ||
        any(!names(convergence_control) %in% names(defaults))) {
      stop("convergence_control must be a uniquely named list of supported controls.",
           call. = FALSE)
    }
    resolved <- utils::modifyList(defaults, convergence_control)
  }
  for (name in c("warn", "keep_traces")) {
    resolved[[name]] <- .blr_logical_scalar(
      resolved[[name]], paste0("convergence_control$", name))
  }
  for (name in c("rhat_threshold", "ess_per_chain_threshold",
                 "mcse_mean_over_sd_threshold")) {
    value <- resolved[[name]]
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
        value <= 0) {
      stop("convergence_control$", name,
           " must be one finite positive number.", call. = FALSE)
    }
  }
  if (convergence == "none" && resolved$keep_traces) {
    stop("convergence = 'none' cannot retain convergence traces.",
         call. = FALSE)
  }
  resolved <- .blr_validate_extended_controls(resolved, convergence)
  resolved$mode <- convergence
  resolved$requested <- convergence != "none"
  resolved$compute <- resolved$requested && nchains >= 2L
  resolved$thresholds <- resolved[c(
    "rhat_threshold", "ess_per_chain_threshold",
    "mcse_mean_over_sd_threshold")]
  resolved
}

.blr_convergence_bundle <- function(values, quantities, family, model,
                                    operator, scope = "core") {
  if (!is.array(values) || length(dim(values)) != 3L ||
      !is.data.frame(quantities) || nrow(quantities) != dim(values)[3L]) {
    stop("values and quantities do not define a scalar convergence bundle.",
         call. = FALSE)
  }
  quantities$family <- family
  quantities$model <- model
  quantities$operator <- operator
  defaults <- list(
    tier = 1L, parameter_name = NA_character_, trait2_index = -1L,
    marker_index = -1L, marker_id = NA_character_, component_index = -1L,
    component_name = NA_character_, pattern_index = -1L,
    pattern_name = NA_character_, annotation_index = -1L,
    annotation_name = NA_character_, stick_index = -1L,
    stick_name = NA_character_, is_intercept = FALSE,
    model_index = -1L, derived = FALSE, structural = FALSE,
    captured = TRUE)
  for (name in names(defaults)) {
    if (is.null(quantities[[name]])) quantities[[name]] <- defaults[[name]]
  }
  if (is.null(quantities$diagnostic_key)) {
    quantities$diagnostic_key <- paste(
      family, model, operator, quantities$group, quantities$parameter_name,
      quantities$trait_index,
      quantities$trait2_index, quantities$marker_index,
      quantities$component_index, quantities$pattern_index,
      quantities$annotation_index, quantities$stick_index,
      sep = ":")
  }
  order <- c(
    "quantity_index", "family", "model", "operator", "tier", "group",
    "parameter_name", "trait_index", "trait2_index", "marker_index",
    "marker_id", "component_index", "component_name", "pattern_index",
    "pattern_name", "annotation_index", "annotation_name", "stick_index",
    "stick_name", "is_intercept", "model_index", "updated", "derived",
    "structural", "captured", "diagnostic_key")
  quantities <- quantities[order]
  list(
    schema = list(class = "blr_convergence_trace_bundle", version = 1L),
    scope = scope, family = family, model = model, operator = operator,
    nchains = as.integer(dim(values)[2L]),
    postburn_draws_per_chain = as.integer(dim(values)[1L]),
    quantities = quantities, values = values)
}

.blr_public_model_error <- function(method, operator, valid) {
  data_level <- if (operator == "packed_bed") "individual" else
    "summary_statistics"
  stop(sprintf(
    paste0("method = '%s' is invalid for operator = '%s' and data_level = ",
           "'%s'. Valid method values: %s. The 's' prefix denotes summary ",
           "statistics; it does not activate maf_effect_s scaling."),
    paste(method, collapse = ","), operator, data_level,
    paste(sprintf("'%s'", valid), collapse = ", ")), call. = FALSE)
}

.blr_model_semantics <- function(method, operator, maf_effect_s = NULL,
                                 estimate_maf_effect_s = FALSE,
                                 probability_policy = NULL) {
  stopifnot(length(method) == 1L, length(operator) == 1L)
  data_level <- if (operator == "packed_bed") "individual" else
    "summary_statistics"
  prior_kernel <- switch(method,
    bayesc = "bayesc", sbayesc = "bayesc",
    bayesr = "bayesr", sbayesr = "bayesr",
    bayesrc = "bayesrc", sbayesrc = "bayesrc")
  maf_effect_s_active <- !is.null(maf_effect_s) || isTRUE(estimate_maf_effect_s)
  effect_scale_policy <- if (prior_kernel == "bayesc") {
    if (maf_effect_s_active) "maf_s" else "unit"
  } else {
    if (maf_effect_s_active) "component_maf_s" else "component"
  }
  list(
    model = method, prior_kernel = prior_kernel, data_level = data_level,
    effect_scale_policy = effect_scale_policy,
    maf_effect_s_active = maf_effect_s_active,
    probability_policy = probability_policy %||%
      if (prior_kernel == "bayesrc") "annotation_probit_stick" else "global",
    model_semantics_version = 2L,
    model_semantics = "s_prefix_means_summary_statistics")
}

.blr_resolve_st_model <- function(method, dots, supported, operator) {
  if (length(method) > 1L) method <- method[[1L]]
  if (!is.character(method) || length(method) != 1L || is.na(method) ||
      !method %in% supported) {
    .blr_public_model_error(method, operator, supported)
  }
  if (!is.null(dots$maf_effect_s) &&
      (!is.numeric(dots$maf_effect_s) || length(dots$maf_effect_s) != 1L ||
       !is.finite(dots$maf_effect_s))) {
    stop("maf_effect_s must be NULL or a finite numeric scalar.", call. = FALSE)
  }
  estimate <- isTRUE(dots$estimate_maf_effect_s %||% FALSE)
  semantics <- .blr_model_semantics(
    method, operator, dots$maf_effect_s, estimate)
  internal_kernel <- switch(method,
    sbayesc = "bayesc", sbayesr = "bayesr", sbayesrc = "sbayesrc", method)
  c(semantics, list(kernel = internal_kernel, dots = dots,
                    effect_scale = semantics$effect_scale_policy))
}

.blr_resolve_st_effect_maf <- function(
    effect_maf, allow_reference_maf_for_maf_effect_s,
    maf_effect_s_active, stats, Glist) {
  if (!is.logical(allow_reference_maf_for_maf_effect_s) ||
      length(allow_reference_maf_for_maf_effect_s) != 1L ||
      is.na(allow_reference_maf_for_maf_effect_s)) {
    stop("allow_reference_maf_for_maf_effect_s must be TRUE or FALSE.",
         call. = FALSE)
  }
  ids <- unlist(Glist$rsidsLD %||% list(), use.names = FALSE)
  validate <- function(x, source) {
    x <- as.numeric(x)
    if (!length(ids) || length(x) != length(ids) || any(!is.finite(x)) ||
        any(x <= 0 | x > 0.5)) {
      stop(source, " effect_maf must be aligned to the final LD order and lie in (0, 0.5].",
           call. = FALSE)
    }
    x
  }
  apply_to_glist <- function(x) {
    offset <- 0L
    for (chr in seq_along(Glist$rsidsLD)) {
      count <- length(Glist$rsidsLD[[chr]])
      values <- x[seq.int(offset + 1L, offset + count)]
      idx <- match(Glist$rsidsLD[[chr]], Glist$rsids[[chr]])
      if (anyNA(idx)) stop("Could not align effect_maf to Glist$rsids.", call. = FALSE)
      target <- Glist$maf[[chr]] %||% rep(NA_real_, length(Glist$rsids[[chr]]))
      target[idx] <- values
      Glist$maf[[chr]] <- target
      offset <- offset + count
    }
    Glist
  }
  result <- list(
    Glist = Glist, effect_maf_source = "not_requested",
    effect_maf_population = "not_applicable",
    effect_maf_alignment_status = "not_requested",
    effect_maf_fallback_used = FALSE)
  if (!isTRUE(maf_effect_s_active) && is.null(effect_maf)) return(result)
  summary_maf <- NULL
  metadata <- stats$marker_metadata
  if (is.data.frame(metadata) && "allele_frequency" %in% names(metadata)) {
    summary_maf <- pmin(metadata$allele_frequency,
                        1 - metadata$allele_frequency)
  }
  if (!is.null(effect_maf)) {
    value <- validate(effect_maf, "Explicit")
    result$effect_maf_source <- "explicit_effect_maf"
    result$effect_maf_population <- "user_declared"
  } else if (!is.null(summary_maf) && all(is.finite(summary_maf))) {
    value <- validate(summary_maf, "GWAS-summary")
    result$effect_maf_source <- "gwas_summary_allele_frequency"
    result$effect_maf_population <- "gwas_summary_population"
  } else if (identical(stats$source %||% NULL, "make_summary_stats")) {
    result$effect_maf_source <- "analysis_genotype_frequency"
    result$effect_maf_population <- "analysis_sample"
    result$effect_maf_alignment_status <- "by_construction"
    return(result)
  } else if (isTRUE(allow_reference_maf_for_maf_effect_s)) {
    result$effect_maf_source <- "reference_panel_frequency"
    result$effect_maf_population <- "reference_panel"
    result$effect_maf_alignment_status <- "aligned_to_final_marker_order"
    result$effect_maf_fallback_used <- TRUE
    return(result)
  } else {
    stop(paste0(
      "maf_effect_s requires explicit effect_maf or aligned GWAS-summary ",
      "MAF. Reference-panel fallback requires ",
      "allow_reference_maf_for_maf_effect_s = TRUE."), call. = FALSE)
  }
  result$Glist <- apply_to_glist(value)
  result$effect_maf_alignment_status <- "aligned_to_final_marker_order"
  result
}

.blr_validate_fit_model_semantics <- function(fit) {
  if (!identical(fit$input$model_semantics_version, 2L) ||
      !identical(fit$input$model_semantics,
                 "s_prefix_means_summary_statistics")) {
    stop("The fit lacks model-semantics version 2 and cannot be silently reinterpreted.",
         call. = FALSE)
  }
  invisible(fit)
}

.blr_model_capability_matrix <- function() {
  models <- c("bayesc", "sbayesc", "bayesr", "sbayesr",
              "bayesrc", "sbayesrc")
  operators <- c("csr", "block_eigen", "packed_bed")
  base <- expand.grid(
    family = c("stblr", "mtblr"), model = models,
    operator = operators, stringsAsFactors = FALSE)
  base$annotation_policy <- ifelse(
    base$model %in% c("bayesrc", "sbayesrc"),
    "annotation_probit_stick", "global")
  base$status <- "unsupported"
  supported <- list(
    stblr = list(
      csr = c("sbayesc", "sbayesr"),
      block_eigen = c("sbayesc", "sbayesr", "sbayesrc"),
      packed_bed = c("bayesc", "bayesr", "bayesrc")),
    mtblr = list(
      csr = c("sbayesc", "sbayesr"),
      block_eigen = c("sbayesc", "sbayesr"),
      packed_bed = c("bayesc", "bayesr")))
  for (family in names(supported)) {
    for (operator in names(supported[[family]])) {
      hit <- base$family == family & base$operator == operator &
        base$model %in% supported[[family]][[operator]]
      base$status[hit] <- "public_canonical"
    }
  }
  annotation <- expand.grid(
    family = "stblr", model = "sbayesc", operator = operators,
    annotation_policy = c("fixed_marker", "group", "learned_logistic"),
    stringsAsFactors = FALSE)
  annotation$status <- ifelse(
    annotation$operator == "csr", "public_supported", "unsupported")
  annotation_rc <- data.frame(
    family = "stblr", model = "sbayesrc", operator = "csr",
    annotation_policy = "annotation_probit_stick",
    status = "public_canonical", stringsAsFactors = FALSE)
  rbind(base, annotation, annotation_rc)
}

.blr_memory_estimate <- function(family, operator, m, nt, nchains, ncores,
                                 nit, ntrace, keep_chains,
                                 convergence_quantities = 0L,
                                 keep_traces = FALSE,
                                 operator_bytes = 0,
                                 sample_count = 0L) {
  workers <- min(as.integer(nchains), as.integer(ncores))
  shared_operator <- as.numeric(operator_bytes)
  private_per_worker <- 8 * (4 * m * nt + 3 * nt * nt + sample_count * nt)
  result_per_chain <- 8 * (4 * m * nt + 5 * ntrace * nt + 6 * nt * nt)
  convergence_capture <- 8 * nchains * nit * convergence_quantities
  convergence_workspace <- if (convergence_quantities > 0L) {
    8 * nchains * nit * 12
  } else 0
  compact <- if (keep_chains) result_per_chain * nchains else 0
  retained_convergence <- if (keep_traces) convergence_capture else 0
  formatted <- 8 * (4 * m * nt + 5 * ntrace * nt + 6 * nt * nt)
  total <- shared_operator + workers * private_per_worker +
    nchains * result_per_chain + convergence_capture +
    convergence_workspace + compact + retained_convergence + formatted
  list(
    estimate_kind = "analytical upper-bound estimate",
    measured_rss = FALSE, measured_peak_rss = FALSE,
    family = family, operator = operator,
    requested_workers = as.integer(ncores), used_workers = workers,
    shared_immutable_operator_data_bytes = shared_operator,
    private_sampler_state_per_worker_bytes = private_per_worker,
    result_state_per_logical_chain_bytes = result_per_chain,
    posterior_trace_capture_bytes = 8 * nchains * ntrace * 5 * nt,
    convergence_trace_capture_bytes = convergence_capture,
    convergence_workspace_bytes = convergence_workspace,
    retained_compact_chains_bytes = compact,
    retained_convergence_traces_bytes = retained_convergence,
    formatted_output_bytes = formatted,
    estimated_total_bytes = total,
    estimated_total_gib = total / 1024^3,
    execution_estimated_total_bytes = total,
    execution_estimated_total_gib = total / 1024^3)
}

.blr_memory_warning <- function(memory, memory_warning_gb,
                                convergence_mode = "none",
                                trace_capture = FALSE,
                                trace_retention = FALSE) {
  threshold <- as.numeric(memory_warning_gb)
  if (length(threshold) != 1L || !is.finite(threshold) || threshold <= 0) {
    stop("memory_warning_gb must be one finite positive number.",
         call. = FALSE)
  }
  if (is.finite(memory$estimated_total_gib) &&
      memory$estimated_total_gib > threshold) {
    warning(sprintf(
      paste0("BLR analytical upper-bound memory estimate is %.3f GiB ",
             "(threshold %.3f GiB; family=%s; operator=%s; ",
             "convergence=%s; trace capture=%s; trace retention=%s). ",
             "This is not measured RSS or measured peak RSS."),
      memory$estimated_total_gib, threshold, memory$family, memory$operator,
      convergence_mode, if (trace_capture) "yes" else "no",
      if (trace_retention) "yes" else "no"), call. = FALSE)
  }
  invisible(memory)
}

.blr_st_preflight_memory <- function(stats = NULL, y = NULL, Glist = NULL,
                                     operator, chain, conv,
                                     memory_warning_gb,
                                     trace_spec = NULL) {
  if (isTRUE(conv$aggregate_component_states) && !isTRUE(chain$keep_chains)) {
    stop(
      paste0(
        "aggregate_component_states requires keep_chains = TRUE so retained ",
        "draws keep an explicit chain dimension."
      ),
      call. = FALSE
    )
  }
  if (!is.null(stats)) {
    m <- as.integer(stats$m %||% length(stats$wy[[1L]]))
    nt <- as.integer(length(stats$wy %||% stats$yy))
  } else {
    marker_lists <- Glist$rsids %||% Glist$rsidsLD
    m <- as.integer(sum(lengths(marker_lists)))
    nt <- if (is.null(dim(y))) 1L else as.integer(ncol(y))
  }
  quantity_count <- if (conv$compute || conv$keep_traces) 5L * nt else 0L
  memory <- .blr_memory_estimate(
    "stblr", operator, m, nt, chain$nchains, chain$ncores, chain$nit,
    chain$nit + chain$nburn, chain$keep_chains,
    convergence_quantities = quantity_count,
    keep_traces = conv$keep_traces)
  if (!is.null(trace_spec) && identical(conv$mode, "extended")) {
    component_count <- as.integer(trace_spec$component_count %||% 0L)
    nselected <- length(trace_spec$markers %||% integer())
    probability_count <- if (isTRUE(trace_spec$probability))
      as.integer(trace_spec$probability_quantity_count %||% 0L) * nt else 0L
    annotation_count <- if (isTRUE(trace_spec$annotations))
      as.integer(trace_spec$annotation_quantity_count %||% 0L) * nt else 0L
    selected_b <- if (isTRUE(trace_spec$b)) nselected * nt else 0L
    selected_d <- if (isTRUE(trace_spec$d)) nselected * nt else 0L
    selected_component <- if (isTRUE(trace_spec$component))
      nselected * nt else 0L
    aggregate_component <- if (isTRUE(conv$aggregate_component_states))
      nt * (component_count + 1L + 3L * max(component_count - 1L, 0L)) else 0L
    extended <- .blr_extended_trace_memory(
      chain$nchains, chain$nit,
      probability_count + annotation_count + selected_b,
      selected_d + selected_component + aggregate_component, conv$keep_traces)
    extended$st_probability_trace_bytes <- 8 * chain$nchains * chain$nit *
      probability_count
    extended$st_annotation_group_trace_bytes <- 8 * chain$nchains * chain$nit *
      annotation_count
    extended$st_selected_b_trace_bytes <- 8 * chain$nchains * chain$nit * selected_b
    extended$st_selected_d_trace_bytes <- 4 * chain$nchains * chain$nit * selected_d
    extended$st_selected_component_trace_bytes <- 4 * chain$nchains * chain$nit *
      selected_component
    extended$st_aggregate_component_trace_bytes <- 4 * chain$nchains * chain$nit *
      aggregate_component
    counts <- c(
      chains = chain$nchains, draws = chain$nit,
      probability = probability_count, annotations = annotation_count,
      selected_b = selected_b, selected_d = selected_d,
      selected_component = selected_component,
      aggregate_component = aggregate_component)
    .blr_enforce_trace_guard(extended, conv, counts, nselected)
    memory <- .blr_add_extended_memory(
      memory, list(enabled = TRUE, memory = extended))
  }
  .blr_memory_warning(memory, memory_warning_gb, conv$mode,
                      conv$compute || conv$keep_traces, conv$keep_traces)
  memory
}

.blr_st_native_trace_spec <- function(conv, marker_ids, prior_kernel,
                                      annotations = FALSE,
                                      component_count = 0L,
                                      annotation_quantity_count = 0L) {
  selected <- .blr_resolve_selected_markers(conv$selected_markers, marker_ids)
  quantities <- conv$selected_marker_quantities %||% character()
  if ("component" %in% quantities &&
      !prior_kernel %in% c("bayesr", "bayesrc")) {
    stop("Selected component diagnostics are not applicable to BayesC models.",
         call. = FALSE)
  }
  probability_quantity_count <- if (prior_kernel == "bayesr") {
    if (component_count == 2L) 1L else as.integer(component_count)
  } else 1L
  list(
    markers = as.integer(selected$marker_index - 1L),
    selected = selected,
    probability = identical(conv$mode, "extended") &&
      "probability" %in% conv$extended_groups_resolved,
    probability_quantity_count = probability_quantity_count,
    annotations = identical(conv$mode, "extended") && annotations &&
      "annotations" %in% conv$extended_groups_resolved,
    annotation_quantity_count = as.integer(annotation_quantity_count),
    component_count = as.integer(component_count),
    b = "b" %in% quantities,
    d = "d" %in% quantities,
    component = isTRUE(conv$aggregate_component_states) ||
      "component" %in% quantities)
}

.blr_st_convergence_bundle <- function(chains, trait_names, model, operator,
                                       nit, nburn, updateB = TRUE,
                                       updateE = TRUE) {
  if (!is.list(chains) || length(chains) != length(trait_names) ||
      any(!vapply(chains, is.list, logical(1)))) {
    stop("STBLR convergence requires deterministic trait-by-chain records.",
         call. = FALSE)
  }
  nchains <- unique(lengths(chains))
  if (length(nchains) != 1L || nchains < 1L) {
    stop("Every trait must contain the same positive chain count.",
         call. = FALSE)
  }
  candidates <- list(
    vbs = list(group = "vbs", updated = updateB, derived = FALSE),
    vgs = list(group = "vgs", updated = TRUE, derived = TRUE),
    ves = list(group = "ves", updated = updateE, derived = FALSE),
    vle = list(group = "vle", updated = TRUE, derived = TRUE),
    vld = list(group = "vld", updated = TRUE, derived = TRUE))
  chain_trace <- function(chain, name) {
    chain[[name]] %||% (chain$trace %||% list())[[name]]
  }
  available <- names(candidates)[vapply(names(candidates), function(name) {
    any(vapply(chains, function(by_chain) {
      any(vapply(by_chain, function(chain) !is.null(chain_trace(chain, name)),
                 logical(1)))
    }, logical(1)))
  }, logical(1))]
  descriptors <- list()
  arrays <- list()
  index <- 0L
  for (name in available) for (trait in seq_along(trait_names)) {
    traces <- lapply(chains[[trait]], chain_trace, name = name)
    if (any(vapply(traces, is.null, logical(1)))) next
    matrix_value <- vapply(traces, function(trace) {
      trace <- as.numeric(trace)
      if (length(trace) == nit) return(trace)
      if (length(trace) < nburn + nit) {
        stop("STBLR chain trace is shorter than nburn + nit.", call. = FALSE)
      }
      trace[seq.int(nburn + 1L, nburn + nit)]
    }, numeric(nit))
    index <- index + 1L
    arrays[[index]] <- matrix_value
    descriptors[[index]] <- data.frame(
      quantity_index = index, group = candidates[[name]]$group,
      trait_index = trait, updated = candidates[[name]]$updated,
      derived = candidates[[name]]$derived, stringsAsFactors = FALSE)
  }
  if (!length(arrays)) stop("No eligible STBLR core traces were captured.",
                            call. = FALSE)
  values <- array(unlist(arrays, use.names = FALSE),
                  dim = c(nit, nchains, length(arrays)))
  .blr_convergence_bundle(
    values, do.call(rbind, descriptors), "stblr", model, operator)
}

.blr_flatten_st_chains <- function(chains, trait_names, model, operator,
                                   seeds = NULL) {
  if (is.null(chains)) return(NULL)
  result <- list()
  index <- 0L
  for (trait in seq_along(chains)) for (chain in seq_along(chains[[trait]])) {
    index <- index + 1L
    record <- chains[[trait]][[chain]]
    if (is.list(record$marker)) {
      for (name in c("bm", "dm", "b", "d", "state",
                     "component_probabilities", "comp_prob",
                     "dm_component_mean")) {
        if (is.null(record[[name]]) && !is.null(record$marker[[name]])) {
          record[[name]] <- record$marker[[name]]
        }
      }
      record$marker <- NULL
    }
    if (is.list(record$pi)) {
      record$final_pi <- record$final_pi %||% record$pi$final
      record$mean_pi <- record$mean_pi %||% record$pi$mean
      record$pi <- NULL
    }
    if (!is.null(record$comp_prob)) {
      record$component_probabilities <- record$comp_prob
      record$comp_prob <- NULL
    }
    if (!is.null(record$final_pi)) {
      record$pi_final <- record$final_pi
      record$final_pi <- NULL
    }
    if (!is.null(record$mean_pi)) {
      record$pi_mean <- record$mean_pi
      record$mean_pi <- NULL
    }
    record <- c(list(
      family = "stblr", model = model, operator = operator,
      trait_index = as.integer(trait), trait_name = trait_names[trait],
      chain_index = as.integer(chain),
      seed = if (length(seeds) >= index) seeds[index] else NA_real_,
      retained_draw_count = if (!is.null(record$vbs)) length(record$vbs) else
        if (!is.null(record$trace$vbs)) length(record$trace$vbs) else NA_integer_),
      record)
    result[[index]] <- record
  }
  names(result) <- paste0("task", seq_along(result))
  result
}

.blr_st_task_seeds <- function(chain, ntraits) {
  blr_scalar_seeds_cpp(
    as.integer(chain$seed), as.integer(ntraits), as.integer(chain$nchains),
    if (is.null(chain$chain_seeds_requested)) integer() else
      as.integer(chain$chain_seeds_requested))
}

.blr_operator_reduction_policy <- function(
    family, model, operator_a, operator_b, filtered = FALSE,
    residual = "diagonal") {
  stopifnot(family %in% c("stblr", "mtblr"),
            model %in% c("bayesc", "sbayesc", "bayesr", "sbayesr"))
  operators <- sort(c(operator_a, operator_b))
  if (identical(operators, sort(c("csr", "block_eigen")))) {
    return(if (isTRUE(filtered)) "approximate_filtered_operator" else
             "floating_point_equivalent")
  }
  if (identical(operators, sort(c("csr", "packed_bed")))) {
    if (family == "mtblr" && identical(residual, "full")) {
      return("not_comparable_residual_covariance")
    }
    return("comparable_when_likelihood_contracts_match")
  }
  "not_comparable_operator_contract"
}

.blr_attach_aggregate_component_traces <- function(fit) {
  bundle <- fit$convergence_traces
  if (is.null(bundle) || is.null(bundle$values) || is.null(bundle$quantities))
    return(fit)
  groups <- c(
    component_count = "component_count_trace",
    realized_active_count = "realized_active_count_trace",
    stick_eligible_count = "stick_eligible_count_trace",
    stick_continue_count = "stick_continue_count_trace",
    stick_stop_count = "stick_stop_count_trace")
  for (group in names(groups)) {
    index <- which(bundle$quantities$group == group)
    if (!length(index)) next
    value <- bundle$values[, , index, drop = FALSE]
    quantity_names <- bundle$quantities$component_name[index]
    stick <- bundle$quantities$stick_name[index]
    missing <- is.na(quantity_names) | !nzchar(quantity_names)
    quantity_names[missing] <- stick[missing]
    missing <- is.na(quantity_names) | !nzchar(quantity_names)
    quantity_names[missing] <- bundle$quantities$parameter_name[index][missing]
    dimnames(value) <- list(
      draw = paste0("Iter", seq_len(dim(value)[1L])),
      chain = paste0("chain", seq_len(dim(value)[2L])),
      quantity = make.unique(quantity_names))
    storage.mode(value) <- "integer"
    attr(value, "quantity_metadata") <- bundle$quantities[index, , drop = FALSE]
    fit[[groups[[group]]]] <- value
  }
  fit
}

.blr_finalize_fit <- function(fit, family, model, operator,
                              data = NULL, diagnostics = NULL,
                              memory_estimate = NULL) {
  stopifnot(family %in% c("stblr", "mtblr"),
            model %in% c("bayesc", "sbayesc", "bayesr", "sbayesr",
                         "bayesrc", "sbayesrc"),
            operator %in% c("csr", "block_eigen", "packed_bed",
                            "dense_reference"))
  fit <- .blr_attach_aggregate_component_traces(fit)
  fit$family <- family
  fit$model <- model
  fit$operator <- operator
  if (is.null(fit$input)) fit$input <- list()
  fit$input$family <- family
  fit$input$method <- model
  fit$input$model <- model
  fit$input$operator <- operator
  semantics <- .blr_model_semantics(
    model, operator, fit$input$maf_effect_s %||% NULL,
    fit$input$estimate_maf_effect_s %||% FALSE,
    fit$input$probability_policy %||% NULL)
  fit$input$prior_kernel <- semantics$prior_kernel
  fit$input$data_level <- semantics$data_level
  fit$input$effect_scale_policy <- semantics$effect_scale_policy
  fit$input$effect_scale <- semantics$effect_scale_policy
  fit$input$model_semantics_version <- semantics$model_semantics_version
  fit$input$model_semantics <- semantics$model_semantics
  fit$data <- data %||% fit$data %||% list()
  if (is.null(fit$data$alignment) && !is.null(fit$alignment)) {
    fit$data$alignment <- fit$alignment
  }
  if (is.null(fit$data$phenotype_preprocessing) &&
      !is.null(fit$phenotype_preprocessing)) {
    fit$data$phenotype_preprocessing <- fit$phenotype_preprocessing
  }
  nt <- if (!is.null(fit$bm) && length(dim(fit$bm)) == 2L) {
    ncol(fit$bm)
  } else as.integer(fit$input$nt %||% 1L)
  n_by_trait <- fit$data$n_by_trait %||% fit$data$sample_size %||%
    fit$input$n
  if (is.null(n_by_trait)) n_by_trait <- rep(NA_integer_, nt)
  n_by_trait <- rep(as.integer(n_by_trait), length.out = nt)
  fit$data$n_by_trait <- n_by_trait
  fit$data$n_total <- as.integer(fit$data$n_total %||%
    fit$input$n_total %||% NA_integer_)
  fit$data$n_used <- as.integer(fit$data$n_used %||%
    fit$input$n_used %||% if (length(unique(n_by_trait)) == 1L)
      n_by_trait[1L] else NA_integer_)
  fit$data$genotype_scale <- fit$data$genotype_scale %||%
    fit$input$genotype_scale %||% fit$data$scale %||%
    "standardized_genotype"
  fit$data$data_level <- semantics$data_level
  fit$data$effect_scale <- semantics$effect_scale_policy
  fit$data$model_semantics_version <- semantics$model_semantics_version
  fit$data$model_semantics <- semantics$model_semantics
  fit$data$effect_maf_source <- fit$data$effect_maf_source %||%
    fit$input$effect_maf_source %||% "not_requested"
  fit$data$effect_maf_alignment_status <-
    fit$data$effect_maf_alignment_status %||%
    fit$input$effect_maf_alignment_status %||% "not_requested"
  fit$data$effect_maf_fallback_used <- isTRUE(
    fit$data$effect_maf_fallback_used %||%
      fit$input$effect_maf_fallback_used %||% FALSE)
  fit$data$phenotype_scale <- fit$data$phenotype_scale %||%
    if (operator == "packed_bed") "centered_unscaled" else
      "summary_crossproduct"
  fit$data$ld_scale <- fit$data$ld_scale %||%
    if (operator %in% c("csr", "block_eigen")) "xtx_crossproduct" else
      "not_applicable_individual_data"
  fit$diagnostics <- diagnostics %||% fit$diagnostics %||% list()
  diagnostic_fields <- c(
    "log_cpo", "mean_log_cpo", "nsamples", "n_used", "pitrait",
    "pimarker", "ld_swap", "ld_swap_chains")
  for (name in diagnostic_fields) if (name %in% names(fit)) {
    target <- switch(name, pitrait = "trait_diagnostics",
                     pimarker = "marker_diagnostics", name)
    fit$diagnostics[[target]] <- fit[[name]]
    fit[[name]] <- NULL
  }
  fit$memory_estimate <- memory_estimate %||% fit$memory_estimate %||% list(
    estimate_kind = "analytical upper-bound estimate",
    measured_rss = FALSE, measured_peak_rss = FALSE,
    estimated_total_bytes = NA_real_, estimated_total_gib = NA_real_,
    execution_estimated_total_bytes = NA_real_,
    execution_estimated_total_gib = NA_real_)
  if (!"convergence" %in% names(fit)) fit["convergence"] <- list(NULL)
  if (!"convergence_traces" %in% names(fit)) {
    fit["convergence_traces"] <- list(NULL)
  }
  if (!"chains" %in% names(fit)) fit["chains"] <- list(NULL)

  rename <- c(
    pi = "pi_final", pim = "pi_mean", pis = "pi_trace",
    comp_prob = "component_probabilities",
    covb = "cov_b_mean", covg = "cov_g_mean", cove = "cov_e_mean",
    vb = "cov_b_final", vg = "cov_g_final", ve = "cov_e_final",
    bm_sd = "bm_chain_mean_sd", bm_min = "bm_chain_mean_min",
    bm_max = "bm_chain_mean_max", dm_sd = "dm_chain_mean_sd",
    dm_min = "dm_chain_mean_min", dm_max = "dm_chain_mean_max",
    maf_effect_s_sd = "maf_effect_s_chain_mean_sd",
    maf_effect_s_min = "maf_effect_s_chain_mean_min",
    maf_effect_s_max = "maf_effect_s_chain_mean_max")
  for (old in names(rename)) if (old %in% names(fit)) {
    fit[[rename[[old]]]] <- fit[[old]]
    fit[[old]] <- NULL
  }
  for (old in c("final_pi", "mean_pi")) fit[[old]] <- NULL
  for (old in c("alignment", "bed_diagnostics", "phenotype_preprocessing")) {
    fit[[old]] <- NULL
  }
  fit$input$memory_estimate <- NULL
  class(fit) <- c(paste0(family, "_fit"), "blr_fit", "list")
  fit
}

.blr_finalize_st_public <- function(fit, model, operator, chain, conv,
                                    memory_warning_gb, verbose = FALSE,
                                    memory = NULL) {
  trait_names <- colnames(fit$bm) %||% paste0("T", seq_len(ncol(fit$bm)))
  original_chains <- fit$chains
  if (isTRUE(conv$full_probability_states)) {
    stop("full_probability_states is available only for MT BayesR/SBayesR.",
         call. = FALSE)
  }
  unavailable_result <- function() {
    groups <- c("vbs", "vgs", "ves", "vle", "vld")
    present <- c(!is.null(fit$vbs), !is.null(fit$vgs), !is.null(fit$ves),
                 !is.null(fit$vle), !is.null(fit$vld))
    groups <- groups[present]
    updated <- c(fit$input$updateB %||% TRUE, TRUE,
                 fit$input$updateE %||% TRUE, TRUE, TRUE)[present]
    .blr_convergence_unavailable(
      trait_names, TRUE, TRUE, chain$nchains, chain$nit, conv$thresholds,
      groups = groups, updated = updated)
  }
  if (conv$compute) {
    bundle <- tryCatch(
      .blr_st_convergence_bundle(
        original_chains, trait_names, model, operator,
        chain$nit, chain$nburn,
        updateB = fit$input$updateB %||% TRUE,
        updateE = fit$input$updateE %||% TRUE),
      error = function(error) {
        if (identical(conditionMessage(error),
                      "No eligible STBLR core traces were captured.")) NULL
        else stop(error)
      })
    extended_bundle <- .blr_st_extended_bundle(
      original_chains, trait_names, model, operator, chain$nit,
      chain$nburn, fit, conv)
    if (is.null(bundle) && is.null(extended_bundle)) {
      fit$convergence <- unavailable_result()
    } else {
      bundle <- if (is.null(bundle)) extended_bundle else if (is.null(extended_bundle))
        bundle else .blr_merge_convergence_bundles(bundle, extended_bundle)
      fit$convergence <- .blr_convergence_tier1(
        bundle, trait_names, conv$thresholds, conv$keep_traces)
      if (conv$keep_traces) {
        fit$convergence_traces <- bundle
        fit$convergence_traces$quantities$quantity <-
          fit$convergence$summary$quantity
        dimnames(fit$convergence_traces$values) <- list(
          paste0("Iter", seq_len(chain$nit)),
          paste0("chain", seq_len(chain$nchains)),
          fit$convergence$summary$quantity)
      }
    }
  } else if (conv$mode == "none") {
    fit$convergence <- .blr_convergence_not_requested(
      trait_names, TRUE, TRUE, chain$nchains, chain$nit, conv$thresholds)
  } else {
    fit$convergence <- unavailable_result()
  }
  task_seeds <- .blr_st_task_seeds(chain, length(trait_names))
  fit$chains <- if (chain$keep_chains) .blr_flatten_st_chains(
    original_chains, trait_names, model, operator, task_seeds) else NULL
  fit$input$chain_seeds_requested <- chain$chain_seeds_requested
  fit$input$task_seeds_resolved <- as.numeric(task_seeds)
  fit$input$nchains <- chain$nchains
  fit$input$ncores <- chain$ncores
  fit$input$keep_chains <- chain$keep_chains
  fit$input$convergence <- conv$mode
  fit$input$convergence_control <- conv[c(
    "warn", "rhat_threshold", "ess_per_chain_threshold",
    "mcse_mean_over_sd_threshold", "keep_traces",
    "extended_groups_requested", "extended_groups_resolved",
    "selected_markers", "selected_marker_quantities",
    "full_probability_states", "aggregate_component_states",
    "max_trace_gb", "allow_large_traces")]
  fit$input$memory_warning_gb <- memory_warning_gb
  if (isTRUE(conv$warn) && conv$mode != "none" &&
      !(conv$mode == "auto" && chain$nchains == 1L)) {
    messages <- .blr_convergence_warning_messages(
      fit$convergence, conv$mode,
      family = "stblr", operator = operator)
    if (length(messages)) warning(messages[1L], call. = FALSE)
  }
  if (isTRUE(verbose)) {
    print(fit$input[c("model", "nchains", "ncores", "convergence")])
  }
  if (is.null(memory)) memory <- .blr_memory_estimate(
    "stblr", operator, nrow(fit$bm), ncol(fit$bm), chain$nchains,
    chain$ncores, chain$nit, chain$nit + chain$nburn,
    chain$keep_chains,
    convergence_quantities = if (conv$compute || conv$keep_traces)
      nrow(fit$convergence$summary) else 0L,
    keep_traces = conv$keep_traces)
  native_diagnostics <- fit$diagnostics %||% list()
  if (!isTRUE(fit$input$updateLDswap)) {
    fit$ld_swap <- NULL
    fit$ld_swap_chains <- NULL
  }
  if (is.null(fit$ld_swap)) native_diagnostics$ld_swap <- NULL
  if (is.null(fit$ld_swap_chains)) native_diagnostics$ld_swap_chains <- NULL
  out <- .blr_finalize_fit(
    fit, "stblr", model, operator,
    data = list(
      marker_ids = rownames(fit$bm), trait_names = trait_names,
      data_level = fit$input$data_level %||% NA_character_,
      alignment = fit$alignment %||% fit$input$alignment %||% NULL,
      phenotype_preprocessing = fit$phenotype_preprocessing %||% NULL,
      provenance = fit$input[c(
        "bed_files", "ld_prefix", "marker_policy", "scale", "rows",
        "chr", "cls")]),
    diagnostics = list(
      backend = fit$input$backend %||% NA_character_,
      native = native_diagnostics,
      packed_bed = fit$bed_diagnostics %||% NULL),
    memory_estimate = memory)
  .blr_phase1_finalize_st(out, chain, model, operator)
}
