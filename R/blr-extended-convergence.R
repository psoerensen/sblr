.blr_extended_group_values <- c(
  "covariance", "probability", "maf_effect_s", "annotations")
.blr_selected_marker_values <- c("b", "d", "component")

.blr_validate_extended_controls <- function(resolved, mode) {
  groups <- resolved$extended_groups
  if (!is.null(groups)) {
    if (!is.character(groups) || anyNA(groups) || any(!nzchar(groups)) ||
        anyDuplicated(groups) ||
        any(!groups %in% .blr_extended_group_values)) {
      stop("convergence_control$extended_groups must be NULL or a unique subset of covariance, probability, maf_effect_s, and annotations.",
           call. = FALSE)
    }
  }
  markers <- resolved$selected_markers
  if (!is.null(markers)) {
    invalid_shortcut <- is.logical(markers) ||
      (is.character(markers) && length(markers) == 1L &&
       markers %in% c("all", "*"))
    if (invalid_shortcut || !length(markers) || anyNA(markers) ||
        anyDuplicated(markers) ||
        !(is.character(markers) ||
          (is.numeric(markers) && all(is.finite(markers)) &&
           all(markers == floor(markers))))) {
      stop("convergence_control$selected_markers must be unique marker IDs or one-based integer indices; all-marker shortcuts are not supported.",
           call. = FALSE)
    }
  }
  quantities <- resolved$selected_marker_quantities
  if (!is.character(quantities) || !length(quantities) || anyNA(quantities) ||
      any(!nzchar(quantities)) || anyDuplicated(quantities) ||
      any(!quantities %in% .blr_selected_marker_values)) {
    stop("convergence_control$selected_marker_quantities must be a unique nonempty subset of b, d, and component.",
         call. = FALSE)
  }
  for (name in c("full_probability_states", "aggregate_component_states",
                 "allow_large_traces")) {
    resolved[[name]] <- .blr_logical_scalar(
      resolved[[name]], paste0("convergence_control$", name))
  }
  if (!is.numeric(resolved$max_trace_gb) ||
      length(resolved$max_trace_gb) != 1L ||
      !is.finite(resolved$max_trace_gb) || resolved$max_trace_gb <= 0) {
    stop("convergence_control$max_trace_gb must be one finite positive number.",
         call. = FALSE)
  }
  used_extended <- !is.null(groups) || !is.null(markers) ||
    isTRUE(resolved$full_probability_states) ||
    isTRUE(resolved$aggregate_component_states) ||
    !identical(quantities, c("b", "d")) ||
    !identical(resolved$max_trace_gb, 1) ||
    isTRUE(resolved$allow_large_traces)
  if (mode == "none" && used_extended) {
    stop("Extended convergence controls cannot be used with convergence = 'none'.",
         call. = FALSE)
  }
  if (!is.null(markers) && mode != "extended") {
    stop("Selected-marker diagnostics require convergence = 'extended'.",
         call. = FALSE)
  }
  if (!is.null(groups) && mode != "extended") {
    stop("extended_groups requires convergence = 'extended'.",
         call. = FALSE)
  }
  if (is.null(markers) && !identical(quantities, c("b", "d"))) {
    stop("selected_marker_quantities requires selected_markers.",
         call. = FALSE)
  }
  if (isTRUE(resolved$full_probability_states) && mode != "extended") {
    stop("full_probability_states requires convergence = 'extended'.",
         call. = FALSE)
  }
  if (isTRUE(resolved$aggregate_component_states) && mode != "extended") {
    stop("aggregate_component_states requires convergence = 'extended'.",
         call. = FALSE)
  }
  resolved$extended_groups_requested <- groups
  resolved$extended_groups_resolved <- if (mode == "extended")
    groups %||% .blr_extended_group_values else character()
  resolved
}

.blr_resolve_selected_markers <- function(selected_markers, marker_ids) {
  if (is.null(selected_markers)) return(data.frame(
    requested_order = integer(), requested_marker = character(),
    marker_index = integer(), marker_id = character(),
    stringsAsFactors = FALSE))
  if (!is.character(marker_ids) || anyNA(marker_ids) ||
      any(!nzchar(marker_ids)) || anyDuplicated(marker_ids)) {
    stop("Final marker metadata must contain unique nonempty marker IDs.",
         call. = FALSE)
  }
  if (is.character(selected_markers)) {
    index <- match(selected_markers, marker_ids)
    if (anyNA(index)) stop("Unknown selected marker ID(s): ",
                           paste(selected_markers[is.na(index)], collapse = ", "),
                           call. = FALSE)
    requested <- selected_markers
  } else {
    index <- as.integer(selected_markers)
    if (any(index < 1L | index > length(marker_ids))) {
      stop("Selected marker indices must refer to the final one-based marker order.",
           call. = FALSE)
    }
    requested <- as.character(selected_markers)
  }
  data.frame(
    requested_order = seq_along(index), requested_marker = requested,
    marker_index = index, marker_id = marker_ids[index],
    stringsAsFactors = FALSE)
}

.blr_extended_trace_memory <- function(nchains, nit, numeric_quantities,
                                       state_quantities, keep_traces,
                                       descriptor_count = NULL) {
  values <- c(nchains, nit, numeric_quantities, state_quantities)
  if (any(!is.finite(values)) || any(values < 0) ||
      any(values != floor(values)) || nchains < 1 || nit < 1) {
    stop("Extended trace memory dimensions must be nonnegative integers.",
         call. = FALSE)
  }
  descriptors <- descriptor_count %||% (numeric_quantities + state_quantities)
  numeric_bytes <- 8 * nchains * nit * numeric_quantities
  state_bytes <- 4 * nchains * nit * state_quantities
  capture <- numeric_bytes + state_bytes
  workspace <- if (numeric_quantities + state_quantities) 8 * nchains * nit * 6 else 0
  retained <- if (isTRUE(keep_traces))
    8 * nchains * nit * (numeric_quantities + state_quantities) else 0
  descriptor_bytes <- 256 * descriptors
  summary_bytes <- 8 * 48 * descriptors
  total <- capture + workspace + retained + descriptor_bytes + summary_bytes
  list(
    numeric_trace_capture_bytes = numeric_bytes,
    state_trace_capture_bytes = state_bytes,
    trace_capture_bytes = capture,
    workspace_bytes = workspace,
    retained_trace_bytes = retained,
    descriptor_metadata_bytes = descriptor_bytes,
    summary_output_bytes = summary_bytes,
    estimated_total_bytes = total,
    estimated_total_gib = total / 1024^3,
    measured_rss = FALSE, measured_peak_rss = FALSE)
}

.blr_enforce_trace_guard <- function(memory, controls, counts,
                                     selected_marker_count) {
  guarded_bytes <- memory$trace_capture_bytes + memory$retained_trace_bytes
  if (guarded_bytes / 1024^3 <= controls$max_trace_gb) {
    return(invisible(FALSE))
  }
  detail <- sprintf(
    paste0("requested diagnostic trace capture: groups=%s; chains=%d; ",
           "draws/chain=%d; selected markers=%d; estimated=%.6f GiB; ",
           "threshold=%.6f GiB; retained=%s"),
    paste(names(counts), counts, sep = "=", collapse = ","),
    counts[["chains"]], counts[["draws"]], selected_marker_count,
    guarded_bytes / 1024^3, controls$max_trace_gb,
    if (controls$keep_traces) "yes" else "no")
  if (!isTRUE(controls$allow_large_traces)) {
    stop(detail, "; set allow_large_traces = TRUE to proceed without truncation.",
         call. = FALSE)
  }
  warning(detail, "; proceeding because allow_large_traces = TRUE.",
          call. = FALSE)
  invisible(TRUE)
}

.blr_add_extended_memory <- function(memory, extended) {
  if (is.null(extended) || !isTRUE(extended$enabled)) return(memory)
  bytes <- extended$memory$estimated_total_bytes
  memory$extended_diagnostics <- extended$memory
  memory$estimated_total_bytes <- memory$estimated_total_bytes + bytes
  memory$estimated_total_gib <- memory$estimated_total_bytes / 1024^3
  memory$execution_estimated_total_bytes <-
    (memory$execution_estimated_total_bytes %||% memory$estimated_total_bytes - bytes) + bytes
  memory$execution_estimated_total_gib <-
    memory$execution_estimated_total_bytes / 1024^3
  memory
}

.blr_convergence_group_overview <- function(summary) {
  groups <- c("core", "covariance", "probability", "maf_effect_s",
              "annotations", "selected_markers")
  category <- ifelse(summary$tier == 1L, "core",
    ifelse(summary$tier == 3L, "selected_markers",
      ifelse(summary$group %in% c("cov_b", "cov_g", "cov_e"),
             "covariance",
      ifelse(grepl("pi", summary$group, fixed = TRUE), "probability",
      ifelse(summary$group == "maf_effect_s", "maf_effect_s", "annotations")))))
  setNames(lapply(groups, function(group) {
    x <- summary[category == group, , drop = FALSE]
    if (!nrow(x)) return(list(
      requested = FALSE, applicable = FALSE, captured = FALSE,
      n_computed = 0L, n_unavailable = 0L, n_flagged = 0L,
      overall_status = "not_requested", max_rhat = NA_real_,
      min_ess_bulk_per_chain = NA_real_, min_ess_tail_per_chain = NA_real_,
      max_mcse_mean_over_sd = NA_real_))
    overview <- .blr_convergence_overview(x)
    list(
      requested = TRUE,
      applicable = any(!x$status %in% c("not_applicable", "structural_zero")),
      captured = any(x$captured), n_computed = overview$n_computed,
      n_unavailable = overview$n_unavailable,
      n_flagged = overview$n_flagged,
      overall_status = overview$overall_status,
      max_rhat = overview$max_rhat,
      min_ess_bulk_per_chain = overview$min_ess_bulk_per_chain,
      min_ess_tail_per_chain = overview$min_ess_tail_per_chain,
      max_mcse_mean_over_sd = overview$max_mcse_mean_over_sd)
  }), groups)
}

.blr_merge_convergence_bundles <- function(...) {
  bundles <- Filter(Negate(is.null), list(...))
  if (!length(bundles)) stop("At least one convergence bundle is required.",
                             call. = FALSE)
  bundles <- lapply(bundles, .blr_validate_convergence_trace_bundle)
  reference <- bundles[[1L]]
  if (any(vapply(bundles, function(x)
      x$nchains != reference$nchains ||
      x$postburn_draws_per_chain != reference$postburn_draws_per_chain ||
      x$family != reference$family || x$model != reference$model ||
      x$operator != reference$operator, logical(1)))) {
    stop("Convergence bundles cannot be merged across executions.", call. = FALSE)
  }
  quantities <- do.call(rbind, lapply(bundles, `[[`, "quantities"))
  if (anyDuplicated(quantities$diagnostic_key)) {
    stop("Duplicate convergence diagnostic_key values are not allowed.",
         call. = FALSE)
  }
  values <- array(unlist(lapply(bundles, `[[`, "values"), use.names = FALSE),
                  dim = c(reference$postburn_draws_per_chain,
                          reference$nchains, nrow(quantities)))
  quantities$quantity_index <- seq_len(nrow(quantities))
  .blr_convergence_bundle(values, quantities, reference$family,
                          reference$model, reference$operator,
                          scope = if (any(quantities$tier > 1L)) "extended" else "core")
}

.blr_mtblr_extended_plan <- function(controls, marker_ids, trait_names,
                                     model, mixture = NULL, bayesrc = NULL,
                                     updateB = TRUE, updateE = TRUE,
                                     updatePi = TRUE,
                                     residual_covariance = "diagonal",
                                     nchains, nit) {
  selected <- .blr_resolve_selected_markers(
    controls$selected_markers, marker_ids)
  enabled <- identical(controls$mode, "extended")
  groups <- if (enabled) controls$extended_groups_resolved else character()
  mixture_parameters <- mixture$model_parameters$mixture %||% list()
  annotation_parameters <- bayesrc$model_parameters$annotations %||% list()
  prior_kernel <- .blr_model_semantics(model,
    if (grepl("^s", model)) "csr" else "packed_bed")$prior_kernel
  component_names <- annotation_parameters$component_names %||%
    mixture_parameters$component_names %||% character()
  joint_names <- mixture_parameters$joint_state_names %||% character()
  pattern_names <- annotation_parameters$pattern_names
  if (is.null(pattern_names) && length(joint_names)) {
    patterns <- mixture_parameters$trait_patterns
    active_count <- if (is.matrix(patterns)) sum(rowSums(patterns) > 0L) else 0L
    pattern_names <- paste0("pattern_", seq_len(active_count))
  }
  if (is.null(pattern_names)) pattern_names <- character()
  if (!length(pattern_names) && prior_kernel == "bayesc")
    pattern_names <- mixture$patterns$names %||% character()
  pattern_group <- "pattern_pi"
  if (prior_kernel == "bayesc" && is.matrix(mixture$patterns$matrix) &&
      nrow(mixture$patterns$matrix) == 2L) {
    null <- which(rowSums(mixture$patterns$matrix != 0L) == 0L)
    if (length(null) == 1L) {
      pattern_names <- pattern_names[-null]
      pattern_group <- "pi_active"
    }
  }
  annotation_names <- annotation_parameters$processed_annotation_names %||%
    annotation_parameters$annotation_column_names %||% character()
  stick_names <- annotation_parameters$stick_names %||% character()
  if (!length(stick_names) && length(component_names) > 1L)
    stick_names <- paste0(component_names[-length(component_names)], "_stick")
  selected_quantities <- if (nrow(selected))
    controls$selected_marker_quantities else character()
  if (isTRUE(controls$full_probability_states) && prior_kernel != "bayesr") {
    stop("full_probability_states is available only for MT BayesR/SBayesR.",
         call. = FALSE)
  }
  if (isTRUE(controls$aggregate_component_states) &&
      !prior_kernel %in% c("bayesr", "bayesrc")) {
    stop("aggregate_component_states is available only for MT BayesR/SBayesR and BayesRC/SBayesRC.",
         call. = FALSE)
  }
  if ("component" %in% selected_quantities && prior_kernel == "bayesc") {
    stop("Selected component diagnostics are not applicable to BayesC models.",
         call. = FALSE)
  }
  nt <- length(trait_names)
  covariance_count <- if ("covariance" %in% groups) 3L * nt * (nt - 1L) / 2L else 0L
  probability_count <- 0L
  if ("probability" %in% groups) {
    probability_count <- length(pattern_names)
    if (prior_kernel == "bayesr") probability_count <-
      probability_count + length(component_names)
    if (isTRUE(controls$full_probability_states) && prior_kernel == "bayesr")
      probability_count <- probability_count + length(joint_names)
  }
  annotation_count <- if ("annotations" %in% groups && prior_kernel == "bayesrc")
    length(annotation_names) * length(stick_names) + length(stick_names) else 0L
  selected_b <- if ("b" %in% selected_quantities) nrow(selected) * nt else 0L
  selected_d <- if ("d" %in% selected_quantities) nrow(selected) * nt else 0L
  selected_component <- if ("component" %in% selected_quantities)
    nrow(selected) else 0L
  aggregate_component <- if (isTRUE(controls$aggregate_component_states)) {
    length(component_names) + 1L + 3L * max(length(component_names) - 1L, 0L)
  } else 0L
  numeric_quantities <- covariance_count + probability_count + annotation_count + selected_b
  state_quantities <- selected_d + selected_component + aggregate_component
  memory <- .blr_extended_trace_memory(
    nchains, nit, numeric_quantities, state_quantities,
    controls$keep_traces)
  counts <- c(chains = nchains, draws = nit, covariance = covariance_count,
              probability = probability_count, annotations = annotation_count,
              selected_b = selected_b, selected_d = selected_d,
              selected_component = selected_component,
              aggregate_component = aggregate_component)
  .blr_enforce_trace_guard(memory, controls, counts, nrow(selected))
  list(
    enabled = enabled, groups = groups, selected = selected,
    selected_quantities = selected_quantities,
    component_names = component_names, pattern_names = pattern_names,
    pattern_group = pattern_group,
    joint_names = joint_names, annotation_names = annotation_names,
    stick_names = stick_names, trait_names = trait_names,
    prior_kernel = prior_kernel, updateB = updateB, updateE = updateE,
    updatePi = updatePi,
    updateAlpha = isTRUE(bayesrc$updateAlpha),
    residual_covariance = residual_covariance,
    aggregate_component_states = isTRUE(controls$aggregate_component_states),
    memory = memory, counts = counts,
    native = list(
      convergence_covariance = enabled && "covariance" %in% groups,
      convergence_probability = enabled && "probability" %in% groups,
      convergence_annotations = enabled && "annotations" %in% groups,
      convergence_full_probability = enabled &&
        isTRUE(controls$full_probability_states),
      convergence_markers = as.integer(selected$marker_index - 1L),
      convergence_b = "b" %in% selected_quantities,
      convergence_d = "d" %in% selected_quantities,
      convergence_component = isTRUE(controls$aggregate_component_states) ||
        "component" %in% selected_quantities))
}

.blr_mtblr_extended_bundle <- function(raws, plan, model, operator,
                                       nit, nburn) {
  if (is.null(plan) || !isTRUE(plan$enabled)) return(NULL)
  captures <- lapply(raws, `[[`, "convergence_capture")
  if (any(vapply(captures, is.null, logical(1)))) {
    stop("MTBLR extended diagnostics were requested but native capture is absent.",
         call. = FALSE)
  }
  nt <- length(plan$trait_names); nchains <- length(raws)
  arrays <- list(); descriptors <- list(); index <- 0L
  append_matrix <- function(field, descriptor_rows, state = FALSE) {
    capture <- lapply(captures, `[[`, field)
    if (!length(descriptor_rows)) return(invisible(NULL))
    if (any(vapply(capture, is.null, logical(1))))
      stop("Missing native extended trace field: ", field, call. = FALSE)
    matrix_value <- array(NA_real_, c(nit, nchains, length(descriptor_rows)))
    for (chain in seq_len(nchains)) {
      x <- as.matrix(capture[[chain]])
      if (nrow(x) < nburn + nit || ncol(x) != length(descriptor_rows))
        stop("Invalid native extended trace dimensions for ", field,
             ": observed ", nrow(x), " x ", ncol(x), "; expected at least ",
             nburn + nit, " x ", length(descriptor_rows), ".",
             call. = FALSE)
      matrix_value[, chain, ] <- x[seq.int(nburn + 1L, nburn + nit), , drop = FALSE]
    }
    for (q in seq_along(descriptor_rows)) {
      index <<- index + 1L
      arrays[[index]] <<- matrix_value[, , q, drop = FALSE]
      row <- utils::modifyList(list(
        tier = 2L, group = NA_character_, parameter_name = NA_character_,
        trait_index = -1L, trait2_index = -1L, marker_index = -1L,
        marker_id = NA_character_, component_index = -1L,
        component_name = NA_character_, pattern_index = -1L,
        pattern_name = NA_character_, annotation_index = -1L,
        annotation_name = NA_character_, stick_index = -1L,
        stick_name = NA_character_, is_intercept = FALSE,
        model_index = -1L, updated = TRUE, derived = FALSE,
        structural = FALSE, captured = TRUE), descriptor_rows[[q]])
      row$quantity_index <- index
      descriptors[[index]] <<- as.data.frame(row, stringsAsFactors = FALSE)
    }
    invisible(NULL)
  }
  pair_rows <- function(group, updated, structural = FALSE) {
    result <- list()
    for (col in seq_len(nt)) if (col < nt) for (row in seq.int(col + 1L, nt))
      result[[length(result) + 1L]] <- list(
        tier = 2L, group = group, parameter_name = group,
        trait_index = col, trait2_index = row, updated = updated,
        derived = group == "cov_g", structural = structural,
        captured = !structural)
    result
  }
  if (isTRUE(plan$aggregate_component_states)) {
    component_rows <- lapply(seq_along(plan$component_names), function(component)
      list(tier = 2L, group = "component_count",
           parameter_name = "component_count", trait_index = -1L,
           component_index = component,
           component_name = plan$component_names[[component]], updated = TRUE,
           derived = TRUE))
    append_matrix("component_count", component_rows, state = TRUE)
    append_matrix("realized_active_count", list(list(
      tier = 2L, group = "realized_active_count",
      parameter_name = "realized_active_count", trait_index = -1L,
      updated = TRUE, derived = TRUE)), state = TRUE)
    stick_rows <- function(group) lapply(seq_along(plan$stick_names), function(stick)
      list(tier = 2L, group = group, parameter_name = group,
           trait_index = -1L, stick_index = stick,
           stick_name = plan$stick_names[[stick]], updated = TRUE,
           derived = TRUE))
    append_matrix("stick_eligible_count",
                  stick_rows("stick_eligible_count"), state = TRUE)
    append_matrix("stick_continue_count",
                  stick_rows("stick_continue_count"), state = TRUE)
    append_matrix("stick_stop_count", stick_rows("stick_stop_count"),
                  state = TRUE)
  }
  if ("covariance" %in% plan$groups) {
    append_matrix("cov_b", pair_rows("cov_b", plan$updateB))
    append_matrix("cov_g", pair_rows("cov_g", TRUE))
    structural_e <- identical(plan$residual_covariance, "diagonal")
    append_matrix("cov_e", pair_rows("cov_e", plan$updateE,
                                      structural_e))
  }
  named_rows <- function(group, names, updated, identity) lapply(
    seq_along(names), function(i) c(list(
      tier = 2L, group = group, parameter_name = group,
      trait_index = -1L, updated = updated),
      setNames(list(i, names[[i]]), c(paste0(identity, "_index"),
                                      paste0(identity, "_name")))))
  if ("probability" %in% plan$groups) {
    if (plan$prior_kernel == "bayesr") append_matrix(
      "component_pi", named_rows("component_pi", plan$component_names,
                                  plan$updatePi, "component"))
    append_matrix("pattern_pi", named_rows(
      plan$pattern_group, plan$pattern_names, plan$updatePi, "pattern"))
    if (length(plan$joint_names) && isTRUE(
        length(captures[[1L]]$joint_pi) > 0L)) append_matrix(
      "joint_pi", named_rows("joint_pi", plan$joint_names,
                              plan$updatePi, "pattern"))
  }
  if ("annotations" %in% plan$groups && plan$prior_kernel == "bayesrc") {
    alpha_rows <- list()
    for (stick in seq_along(plan$stick_names))
      for (annotation in seq_along(plan$annotation_names))
        alpha_rows[[length(alpha_rows) + 1L]] <- list(
          tier = 2L, group = "annotations", parameter_name = "alpha",
          trait_index = -1L, annotation_index = annotation,
          annotation_name = plan$annotation_names[[annotation]],
          stick_index = stick, stick_name = plan$stick_names[[stick]],
          is_intercept = grepl("^\\(?intercept\\)?$",
                               plan$annotation_names[[annotation]],
                               ignore.case = TRUE),
          updated = plan$updateAlpha)
    append_matrix("alpha", alpha_rows)
    append_matrix("sigmaSqAlpha", lapply(seq_along(plan$stick_names),
      function(stick) list(tier = 2L, group = "annotations",
        parameter_name = "sigmaSqAlpha", trait_index = -1L,
        stick_index = stick, stick_name = plan$stick_names[[stick]],
        updated = plan$updateAlpha)))
  } else if ("annotations" %in% plan$groups) {
    index <- index + 1L
    arrays[[index]] <- array(NA_real_, c(nit, nchains, 1L))
    descriptors[[index]] <- data.frame(
      quantity_index = index, tier = 2L, group = "annotations",
      parameter_name = "annotations", trait_index = -1L,
      updated = FALSE, derived = FALSE, captured = FALSE,
      stringsAsFactors = FALSE)
  }
  if ("maf_effect_s" %in% plan$groups) {
    index <- index + 1L
    arrays[[index]] <- array(NA_real_, c(nit, nchains, 1L))
    descriptors[[index]] <- data.frame(
      quantity_index = index, tier = 2L, group = "maf_effect_s",
      parameter_name = "maf_effect_s", trait_index = -1L,
      updated = FALSE, derived = FALSE, captured = FALSE,
      stringsAsFactors = FALSE)
  }
  selected_rows <- function(parameter, by_trait = TRUE) {
    result <- list()
    for (marker in seq_len(nrow(plan$selected))) {
      traits <- if (by_trait) seq_len(nt) else -1L
      for (trait in traits) result[[length(result) + 1L]] <- list(
        tier = 3L, group = parameter, parameter_name = parameter,
        trait_index = trait, marker_index = plan$selected$marker_index[[marker]],
        marker_id = plan$selected$marker_id[[marker]], updated = TRUE)
    }
    result
  }
  if ("b" %in% plan$selected_quantities)
    append_matrix("b", selected_rows("b"))
  if ("d" %in% plan$selected_quantities)
    append_matrix("d", selected_rows("d"))
  if ("component" %in% plan$selected_quantities)
    append_matrix("component", selected_rows("component", FALSE))
  if (!length(arrays)) return(NULL)
  descriptor_defaults <- list(
    tier = 2L, group = NA_character_, parameter_name = NA_character_,
    trait_index = -1L, trait2_index = -1L, marker_index = -1L,
    marker_id = NA_character_, component_index = -1L,
    component_name = NA_character_, pattern_index = -1L,
    pattern_name = NA_character_, annotation_index = -1L,
    annotation_name = NA_character_, stick_index = -1L,
    stick_name = NA_character_, is_intercept = FALSE, model_index = -1L,
    updated = TRUE, derived = FALSE, structural = FALSE, captured = TRUE)
  descriptors <- lapply(descriptors, function(row) {
    for (name in names(descriptor_defaults)) if (is.null(row[[name]]))
      row[[name]] <- descriptor_defaults[[name]]
    row[c("quantity_index", names(descriptor_defaults))]
  })
  values <- array(unlist(arrays, use.names = FALSE), c(nit, nchains, length(arrays)))
  .blr_convergence_bundle(values, do.call(rbind, descriptors), "mtblr",
                          model, operator, scope = "extended")
}

.blr_mtblr_convert_native_bundle <- function(bundle, plan, model,
                                              operator = "packed_bed") {
  q <- bundle$quantities
  if (!is.data.frame(q)) stop("Invalid native convergence descriptors.",
                              call. = FALSE)
  descriptor <- data.frame(
    quantity_index = seq_len(nrow(q)), tier = q$tier %||% 1L,
    group = as.character(q$group), parameter_name = as.character(q$group),
    trait_index = as.integer(q$trait_index),
    trait2_index = as.integer(q$trait2_index %||% -1L),
    marker_index = as.integer(q$marker_index %||% -1L),
    marker_id = NA_character_, component_index = as.integer(q$component_index %||% -1L),
    component_name = NA_character_, pattern_index = as.integer(q$pattern_index %||% -1L),
    pattern_name = NA_character_, annotation_index = -1L,
    annotation_name = NA_character_, stick_index = as.integer(q$stick_index %||% -1L),
    stick_name = NA_character_, is_intercept = FALSE, model_index = -1L,
    updated = as.logical(q$updated), derived = as.logical(q$derived),
    structural = as.logical(q$structural %||% FALSE),
    captured = !as.logical(q$structural %||% FALSE),
    stringsAsFactors = FALSE)
  binary_pattern <- descriptor$group == "pattern_pi" &
    identical(plan$pattern_group, "pi_active")
  descriptor$group[binary_pattern] <- "pi_active"
  descriptor$parameter_name[binary_pattern] <- "pi_active"
  valid_marker <- descriptor$marker_index > 0L
  if (any(valid_marker)) descriptor$marker_id[valid_marker] <-
    plan$selected$marker_id[match(descriptor$marker_index[valid_marker],
                                  plan$selected$marker_index)]
  valid_component <- descriptor$component_index > 0L
  if (any(valid_component)) descriptor$component_name[valid_component] <-
    plan$component_names[descriptor$component_index[valid_component]]
  valid_pattern <- descriptor$pattern_index > 0L
  for (i in which(valid_pattern)) {
    source <- if (descriptor$group[[i]] == "joint_pi") plan$joint_names else
      plan$pattern_names
    descriptor$pattern_name[[i]] <- source[[descriptor$pattern_index[[i]]]]
  }
  alpha <- descriptor$group == "alpha"
  annotation_count <- length(plan$annotation_names)
  if (any(alpha) && annotation_count) for (i in which(alpha)) {
    flat <- as.integer(q$annotation_index[[i]]) - 1L
    annotation <- flat %% annotation_count + 1L
    stick <- flat %/% annotation_count + 1L
    descriptor$annotation_index[[i]] <- annotation
    descriptor$annotation_name[[i]] <- plan$annotation_names[[annotation]]
    descriptor$stick_index[[i]] <- stick
    descriptor$stick_name[[i]] <- plan$stick_names[[stick]]
    descriptor$is_intercept[[i]] <- grepl(
      "^\\(?intercept\\)?$", plan$annotation_names[[annotation]],
      ignore.case = TRUE)
  }
  sigma <- descriptor$group == "sigmaSqAlpha"
  if (any(sigma)) descriptor$stick_name[sigma] <-
    plan$stick_names[descriptor$stick_index[sigma]]
  .blr_convergence_bundle(bundle$values, descriptor, "mtblr", model,
                          operator, scope = bundle$scope)
}

.blr_st_extended_bundle <- function(chains, trait_names, model, operator,
                                    nit, nburn, fit, controls) {
  if (!identical(controls$mode, "extended")) return(NULL)
  groups <- controls$extended_groups_resolved
  prior_kernel <- fit$input$prior_kernel %||%
    .blr_model_semantics(model, operator)$prior_kernel
  nchains <- unique(lengths(chains))
  if (length(nchains) != 1L || nchains < 1L) return(NULL)
  trace_value <- function(chain, name) chain[[name]] %||%
    (chain$trace %||% list())[[name]]
  arrays <- list(); descriptors <- list(); index <- 0L
  append_trait_trace <- function(name, group, updated, derived = FALSE) {
    for (trait in seq_along(trait_names)) {
      traces <- lapply(chains[[trait]], trace_value, name = name)
      if (any(vapply(traces, is.null, logical(1)))) next
      x <- vapply(traces, function(value) {
        value <- as.numeric(value)
        if (length(value) == nit) return(value)
        if (length(value) < nburn + nit) return(rep(NA_real_, nit))
        value[seq.int(nburn + 1L, nburn + nit)]
      }, numeric(nit))
      index <<- index + 1L; arrays[[index]] <<- x
      descriptors[[index]] <<- data.frame(
        quantity_index = index, tier = 2L, group = group,
        parameter_name = group, trait_index = trait, updated = updated,
        derived = derived, captured = TRUE, stringsAsFactors = FALSE)
    }
  }
  append_native_matrix <- function(field, group, updated, tier = 2L,
                                   state = FALSE, component_names = NULL,
                                   stick_names = NULL, selected = NULL) {
    for (trait in seq_along(trait_names)) {
      matrices <- lapply(chains[[trait]], function(chain) {
        value <- (chain$convergence_trace %||% list())[[field]]
        if (is.null(value)) return(NULL)
        value <- as.matrix(value)
        if (nrow(value) != nit) return(NULL)
        value
      })
      if (any(vapply(matrices, is.null, logical(1)))) next
      quantity_count <- ncol(matrices[[1L]])
      if (!quantity_count || any(vapply(matrices, ncol, integer(1)) !=
                                 quantity_count)) next
      for (quantity in seq_len(quantity_count)) {
        x <- vapply(matrices, function(value) as.numeric(value[, quantity]),
                    numeric(nit))
        index <<- index + 1L; arrays[[index]] <<- x
        descriptor <- data.frame(
          quantity_index = index, tier = tier, group = group,
          parameter_name = group, trait_index = trait, updated = updated,
          derived = FALSE, captured = TRUE, stringsAsFactors = FALSE)
        if (!is.null(component_names)) {
          descriptor$component_index <- quantity
          descriptor$component_name <- component_names[[quantity]]
        }
        if (!is.null(stick_names)) {
          descriptor$stick_index <- quantity
          descriptor$stick_name <- stick_names[[quantity]]
        }
        if (!is.null(selected)) {
          descriptor$marker_index <- selected$marker_index[[quantity]]
          descriptor$marker_id <- selected$marker_id[[quantity]]
        }
        descriptors[[index]] <<- descriptor
      }
    }
  }
  if (isTRUE(controls$aggregate_component_states)) {
    component_matrix <- if (is.list(fit$component_probabilities))
      fit$component_probabilities[[1L]] else fit$component_probabilities
    mixture_var <- fit$input$mixture_var %||% fit$input$gamma %||%
      fit$mixture_var
    component_count <- length(mixture_var)
    if (!component_count && is.matrix(component_matrix))
      component_count <- ncol(component_matrix)
    component_names <- colnames(component_matrix) %||%
      names(mixture_var) %||%
      paste0("component_", seq_len(component_count) - 1L)
    stick_names <- paste0(component_names[-length(component_names)], "_stick")
    append_native_matrix("component_count", "component_count", TRUE,
                         state = TRUE, component_names = component_names)
    append_native_matrix("realized_active_count", "realized_active_count",
                         TRUE, state = TRUE)
    append_native_matrix("stick_eligible_count", "stick_eligible_count",
                         TRUE, state = TRUE, stick_names = stick_names)
    append_native_matrix("stick_continue_count", "stick_continue_count",
                         TRUE, state = TRUE, stick_names = stick_names)
    append_native_matrix("stick_stop_count", "stick_stop_count",
                         TRUE, state = TRUE, stick_names = stick_names)
  }
  if ("probability" %in% groups) {
    if (prior_kernel == "bayesc") append_trait_trace(
      "pis", "pi_active", fit$input$updatePi %||% TRUE)
    if (prior_kernel == "bayesr") {
      component_matrix <- if (is.list(fit$component_probabilities))
        fit$component_probabilities[[1L]] else fit$component_probabilities
      mixture_var <- fit$input$mixture_var %||% fit$mixture_var %||%
        fit$input$gamma
      component_count <- length(mixture_var)
      if (!component_count && is.matrix(component_matrix))
        component_count <- ncol(component_matrix)
      component_names <- colnames(component_matrix) %||%
        names(mixture_var) %||%
        paste0("component_", seq_len(component_count) - 1L)
      physical_names <- if (length(component_names) == 2L)
        component_names[2L] else component_names
      append_native_matrix(
        "component_pi", "component_pi", fit$input$updatePi %||% TRUE,
        component_names = physical_names)
    }
  }
  selected <- .blr_resolve_selected_markers(
    controls$selected_markers,
    fit$data$marker_metadata$marker_id %||% rownames(fit$bm) %||%
      paste0("m", seq_len(nrow(fit$bm))))
  if (nrow(selected)) {
    quantities <- controls$selected_marker_quantities
    if ("b" %in% quantities) append_native_matrix(
      "b", "selected_b", TRUE, tier = 3L, selected = selected)
    if ("d" %in% quantities) append_native_matrix(
      "d", "selected_d", TRUE, tier = 3L, state = TRUE,
      selected = selected)
    if ("component" %in% quantities) append_native_matrix(
      "component", "selected_component", TRUE, tier = 3L, state = TRUE,
      selected = selected)
  }
  if ("annotations" %in% groups &&
      identical(prior_kernel, "bayesrc")) {
    annotation_names <- fit$input$annotation_names %||%
      colnames(fit$input$A) %||% character()
    component_matrix <- if (is.list(fit$component_probabilities))
      fit$component_probabilities[[1L]] else fit$component_probabilities
    mixture_var <- fit$input$mixture_var %||% fit$input$gamma %||%
      fit$mixture_var
    component_count <- length(mixture_var)
    if (!component_count && is.matrix(component_matrix))
      component_count <- ncol(component_matrix)
    component_names <- colnames(component_matrix) %||%
      names(mixture_var) %||%
      paste0("component_", seq_len(component_count) - 1L)
    stick_names <- paste0(component_names[-length(component_names)], "_stick")
    for (trait in seq_along(trait_names)) {
      chain_values <- lapply(chains[[trait]], function(chain)
        (chain$convergence_trace %||% list())$alpha)
      if (length(annotation_names) && length(stick_names) &&
          !any(vapply(chain_values, is.null, logical(1)))) {
        for (stick in seq_along(stick_names))
          for (annotation in seq_along(annotation_names)) {
            q <- (stick - 1L) * length(annotation_names) + annotation
            x <- vapply(chain_values, function(value)
              as.numeric(as.matrix(value)[, q]), numeric(nit))
            index <- index + 1L; arrays[[index]] <- x
            descriptors[[index]] <- data.frame(
              quantity_index = index, tier = 2L, group = "annotations",
              parameter_name = "alpha", trait_index = trait,
              annotation_index = annotation,
              annotation_name = annotation_names[[annotation]],
              stick_index = stick, stick_name = stick_names[[stick]],
              is_intercept = grepl("^\\(?intercept\\)?$",
                annotation_names[[annotation]], ignore.case = TRUE),
              updated = fit$input$updateAlpha %||% TRUE,
              derived = FALSE, captured = TRUE, stringsAsFactors = FALSE)
          }
      }
      sigma_values <- lapply(chains[[trait]], function(chain)
        (chain$convergence_trace %||% list())$sigmaSqAlpha)
      if (length(stick_names) &&
          !any(vapply(sigma_values, is.null, logical(1)))) {
        for (stick in seq_along(stick_names)) {
          x <- vapply(sigma_values, function(value)
            as.numeric(as.matrix(value)[, stick]), numeric(nit))
          index <- index + 1L; arrays[[index]] <- x
          descriptors[[index]] <- data.frame(
            quantity_index = index, tier = 2L, group = "annotations",
            parameter_name = "sigmaSqAlpha", trait_index = trait,
            stick_index = stick, stick_name = stick_names[[stick]],
            updated = fit$input$updateAlpha %||% TRUE,
            derived = FALSE, captured = TRUE, stringsAsFactors = FALSE)
        }
      }
    }
  } else if ("annotations" %in% groups &&
             (!is.null(fit$input$learn_pi_annot) ||
              !is.null(fit$input$learn_vb_annot))) {
    annotation_names <- fit$input$annotation_names %||%
      colnames(fit$input$A) %||% character()
    fields <- list(
      inclusion_coefficient = list(updated = isTRUE(fit$input$learn_pi_annot)),
      variance_coefficient = list(updated = isTRUE(fit$input$learn_vb_annot)))
    for (field in names(fields)) for (trait in seq_along(trait_names)) {
      matrices <- lapply(chains[[trait]], function(chain)
        (chain$convergence_trace %||% list())[[field]])
      if (any(vapply(matrices, is.null, logical(1))) || !length(annotation_names)) next
      for (annotation in seq_along(annotation_names)) {
        x <- vapply(matrices, function(value)
          as.numeric(as.matrix(value)[, annotation]), numeric(nit))
        index <- index + 1L; arrays[[index]] <- x
        descriptors[[index]] <- data.frame(
          quantity_index = index, tier = 2L, group = "annotations",
          parameter_name = field, trait_index = trait,
          annotation_index = annotation,
          annotation_name = annotation_names[[annotation]],
          is_intercept = grepl("^\\(?intercept\\)?$",
            annotation_names[[annotation]], ignore.case = TRUE),
          updated = fields[[field]]$updated, derived = FALSE,
          captured = TRUE, stringsAsFactors = FALSE)
      }
    }
  } else if ("annotations" %in% groups &&
             length(fit$input$group_names %||% character())) {
    group_names <- fit$input$group_names %||%
      paste0("group_", seq_len(fit$input$ngroup %||% 0L))
    for (field in c("group_pi", "group_vb")) {
      parameter <- if (field == "group_pi") "group_pi" else
        "group_variance_multiplier"
      updated <- if (field == "group_pi") fit$input$updatePi %||% TRUE else
        fit$input$updateGroupVb %||% FALSE
      for (trait in seq_along(trait_names)) {
        matrices <- lapply(chains[[trait]], function(chain)
          (chain$convergence_trace %||% list())[[field]])
        if (any(vapply(matrices, is.null, logical(1)))) next
        for (group in seq_along(group_names)) {
          x <- vapply(matrices, function(value)
            as.numeric(as.matrix(value)[, group]), numeric(nit))
          index <- index + 1L; arrays[[index]] <- x
          descriptors[[index]] <- data.frame(
            quantity_index = index, tier = 2L, group = "annotations",
            parameter_name = parameter, trait_index = trait,
            pattern_index = group, pattern_name = group_names[[group]],
            updated = updated, derived = FALSE, captured = TRUE,
            stringsAsFactors = FALSE)
        }
      }
    }
  }
  if ("maf_effect_s" %in% groups) {
    sampled <- isTRUE(fit$input$estimate_maf_effect_s)
    if (sampled) append_trait_trace("maf_effect_s", "maf_effect_s", TRUE)
    else if (!is.null(fit$input$maf_effect_s)) {
      for (trait in seq_along(trait_names)) {
        index <- index + 1L
        arrays[[index]] <- matrix(rep(fit$input$maf_effect_s[[min(
          trait, length(fit$input$maf_effect_s))]], nit * nchains), nit, nchains)
        descriptors[[index]] <- data.frame(
          quantity_index = index, tier = 2L, group = "maf_effect_s",
          parameter_name = "maf_effect_s", trait_index = trait,
          updated = FALSE, derived = FALSE, captured = TRUE,
          stringsAsFactors = FALSE)
      }
    } else for (trait in seq_along(trait_names)) {
      index <- index + 1L
      arrays[[index]] <- matrix(NA_real_, nit, nchains)
      descriptors[[index]] <- data.frame(
        quantity_index = index, tier = 2L, group = "maf_effect_s",
        parameter_name = "maf_effect_s", trait_index = trait,
        updated = FALSE, derived = FALSE, captured = FALSE,
        stringsAsFactors = FALSE)
    }
  }
  if (!length(arrays)) return(NULL)
  values <- array(unlist(arrays, use.names = FALSE), c(nit, nchains, length(arrays)))
  descriptor_names <- unique(unlist(lapply(descriptors, names), use.names = FALSE))
  descriptor_default <- list(
    marker_index = -1L, marker_id = NA_character_,
    component_index = -1L, component_name = NA_character_,
    annotation_index = -1L, annotation_name = NA_character_,
    stick_index = -1L, stick_name = NA_character_, is_intercept = FALSE,
    pattern_index = -1L, pattern_name = NA_character_)
  descriptors <- lapply(descriptors, function(descriptor) {
    for (name in setdiff(descriptor_names, names(descriptor)))
      descriptor[[name]] <- descriptor_default[[name]] %||% NA
    descriptor[descriptor_names]
  })
  .blr_convergence_bundle(values, do.call(rbind, descriptors), "stblr",
                          model, operator, scope = "extended")
}
