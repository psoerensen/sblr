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

.blr_convergence_controls <- function(convergence = c("auto", "none", "core"),
                                      convergence_control = NULL,
                                      nchains = 1L) {
  convergence <- match.arg(convergence)
  defaults <- list(
    warn = TRUE, rhat_threshold = 1.01,
    ess_per_chain_threshold = 100,
    mcse_mean_over_sd_threshold = 0.05, keep_traces = FALSE)
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
  resolved$mode <- convergence
  resolved$requested <- convergence != "none"
  resolved$compute <- resolved$requested && nchains >= 2L
  resolved$thresholds <- resolved[c(
    "rhat_threshold", "ess_per_chain_threshold",
    "mcse_mean_over_sd_threshold")]
  resolved
}

.blr_convergence_bundle <- function(values, quantities, family, model,
                                    operator) {
  if (!is.array(values) || length(dim(values)) != 3L ||
      !is.data.frame(quantities) || nrow(quantities) != dim(values)[3L]) {
    stop("values and quantities do not define a scalar convergence bundle.",
         call. = FALSE)
  }
  quantities$family <- family
  quantities$model <- model
  quantities$operator <- operator
  defaults <- list(
    tier = 1L, trait2_index = -1L, marker_index = -1L,
    model_index = -1L, derived = FALSE, structural = FALSE,
    captured = TRUE)
  for (name in names(defaults)) {
    if (is.null(quantities[[name]])) quantities[[name]] <- defaults[[name]]
  }
  if (is.null(quantities$diagnostic_key)) {
    quantities$diagnostic_key <- paste(
      family, model, operator, quantities$group, quantities$trait_index,
      sep = ":")
  }
  order <- c(
    "quantity_index", "family", "model", "operator", "tier", "group",
    "trait_index", "trait2_index", "marker_index", "model_index",
    "updated", "derived", "structural", "captured", "diagnostic_key")
  quantities <- quantities[order]
  list(
    schema = list(class = "blr_convergence_trace_bundle", version = 1L),
    scope = "core", family = family, model = model, operator = operator,
    nchains = as.integer(dim(values)[2L]),
    postburn_draws_per_chain = as.integer(dim(values)[1L]),
    quantities = quantities, values = values)
}

.blr_resolve_st_model <- function(method, dots, supported) {
  method <- match.arg(method, supported)
  scaled <- method %in% c("sbayesc", "sbayesr", "sbayesrc")
  kernel <- switch(method, sbayesc = "bayesc", sbayesr = "bayesr",
                   sbayesrc = "sbayesrc", method)
  selection_controls <- c(
    "selection_s", "estimate_selection_s", "selection_s_init",
    "selection_s_prior", "selection_s_proposal_sd")
  supplied <- intersect(names(dots), selection_controls)
  estimate <- isTRUE(dots$estimate_selection_s %||% FALSE)
  if (!scaled && length(supplied)) {
    stop("maf_s controls require an S model: sbayesc, sbayesr, or sbayesrc.",
         call. = FALSE)
  }
  if (scaled && is.null(dots$selection_s) && !estimate) {
    dots$selection_s <- 0
  }
  list(
    model = method, kernel = kernel, dots = dots,
    effect_scale = switch(
      method, bayesc = "unit", sbayesc = "maf_s",
      bayesr = "component", sbayesr = "component_maf_s",
      bayesrc = "component", sbayesrc = "component_maf_s"),
    probability_policy = switch(
      method, bayesc = "global", sbayesc = "global",
      bayesr = "global", sbayesr = "global",
      bayesrc = "annotation_probit_stick",
      sbayesrc = "annotation_probit_stick"))
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
      csr = c("bayesc", "sbayesc", "bayesr", "sbayesr", "sbayesrc"),
      block_eigen = c("bayesc", "sbayesc", "bayesr", "sbayesr",
                      "sbayesrc"),
      packed_bed = c("bayesc", "bayesr", "bayesrc")),
    mtblr = list(
      csr = "bayesc", block_eigen = "bayesc", packed_bed = "bayesc"))
  for (family in names(supported)) {
    for (operator in names(supported[[family]])) {
      hit <- base$family == family & base$operator == operator &
        base$model %in% supported[[family]][[operator]]
      base$status[hit] <- ifelse(
        base$model[hit] %in% c("sbayesc", "sbayesr", "sbayesrc"),
        "public_supported", "public_canonical")
    }
  }
  annotation <- expand.grid(
    family = "stblr", model = "bayesc", operator = operators,
    annotation_policy = c("fixed_marker", "group", "learned_logistic"),
    stringsAsFactors = FALSE)
  annotation$status <- ifelse(
    annotation$operator == "csr", "public_supported", "unsupported")
  rbind(base, annotation)
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
                                     memory_warning_gb) {
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
  .blr_memory_warning(memory, memory_warning_gb, conv$mode,
                      conv$compute || conv$keep_traces, conv$keep_traces)
  memory
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
  blr_phase4_scalar_seeds_cpp(
    as.integer(chain$seed), as.integer(ntraits), as.integer(chain$nchains),
    if (is.null(chain$chain_seeds_requested)) integer() else
      as.integer(chain$chain_seeds_requested))
}

.blr_operator_reduction_policy <- function(
    family, model, operator_a, operator_b, filtered = FALSE,
    residual = "diagonal") {
  stopifnot(family %in% c("stblr", "mtblr"), model == "bayesc")
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

.blr_finalize_fit <- function(fit, family, model, operator,
                              data = NULL, diagnostics = NULL,
                              memory_estimate = NULL) {
  stopifnot(family %in% c("stblr", "mtblr"),
            model %in% c("bayesc", "sbayesc", "bayesr", "sbayesr",
                         "bayesrc", "sbayesrc"),
            operator %in% c("csr", "block_eigen", "packed_bed",
                            "dense_reference"))
  fit$family <- family
  fit$model <- model
  fit$operator <- operator
  if (is.null(fit$input)) fit$input <- list()
  fit$input$family <- family
  fit$input$model <- model
  fit$input$operator <- operator
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
  fit$data$effect_scale <- fit$input$effect_scale %||%
    switch(model, sbayesc = "maf_s", sbayesr = "component_maf_s",
           sbayesrc = "component_maf_s", bayesr = "component",
           bayesrc = "component", "unit")
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
    selection_s_sd = "selection_s_chain_mean_sd",
    selection_s_min = "selection_s_chain_mean_min",
    selection_s_max = "selection_s_chain_mean_max")
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
    if (is.null(bundle)) {
      fit$convergence <- unavailable_result()
    } else {
      fit$convergence <- .blr_convergence_tier1(
        bundle, trait_names, conv$thresholds, conv$keep_traces)
      if (conv$keep_traces) {
        fit$convergence_traces <- bundle
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
    "mcse_mean_over_sd_threshold", "keep_traces")]
  fit$input$memory_warning_gb <- memory_warning_gb
  if (isTRUE(conv$warn) && conv$mode != "none" &&
      !(conv$mode == "auto" && chain$nchains == 1L)) {
    messages <- .blr_convergence_warning_messages(
      fit$convergence, if (conv$mode == "core") "core" else "auto",
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
  .blr_finalize_fit(
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
}
