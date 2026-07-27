.mtblr_resolve_public_method <- function(method, operator) {
  valid <- switch(operator,
    packed_bed = c("bayesc", "bayesr"),
    csr = c("sbayesc", "sbayesr"),
    block_eigen = c("sbayesc", "sbayesr"))
  if (!is.character(method) || length(method) != 1L || is.na(method) ||
      !method %in% valid) {
    .blr_public_model_error(method, operator, valid)
  }
  .blr_model_semantics(method, operator)
}

.mtblr_bayesr_spec <- function(prior_kernel, pattern_spec, marker_frequency,
                               marker_count, mixture_var,
                               joint_pi = NULL, joint_pi_prior = NULL,
                               component = NULL, selection_s = NULL,
                               estimate_selection_s = FALSE,
                               selection_s_init = NULL,
                               selection_s_prior = NULL,
                               selection_s_proposal_sd = NULL) {
  if (!is.character(prior_kernel) || length(prior_kernel) != 1L ||
      is.na(prior_kernel) || !prior_kernel %in% c("bayesc", "bayesr")) {
    stop("prior_kernel must be exactly 'bayesc' or 'bayesr'.",
         call. = FALSE)
  }
  if (prior_kernel == "bayesc") {
    supplied <- list(mixture_var = mixture_var, joint_pi = joint_pi,
      joint_pi_prior = joint_pi_prior, component = component,
      selection_s = selection_s, selection_s_init = selection_s_init,
      selection_s_prior = selection_s_prior,
      selection_s_proposal_sd = selection_s_proposal_sd)
    if (any(!vapply(supplied, is.null, logical(1))) ||
        !identical(estimate_selection_s, FALSE)) {
      stop("BayesR component and selection-S controls require the BayesR prior kernel.",
           call. = FALSE)
    }
    return(list(patterns = pattern_spec, method_code = 4L,
      joint_component = integer(), joint_multiplier = numeric(),
      joint_names = character(), component_count = 0L,
      marker_scale = numeric(), component_init = integer(),
      pi_prior = numeric(), selection_s = NULL,
      mixture_var = NULL, model_parameters = NULL))
  }
  if (is.null(mixture_var)) mixture_var <- c(0, 0.01, 0.1, 1)
  if (!is.numeric(mixture_var) || length(mixture_var) < 2L ||
      any(!is.finite(mixture_var)) || sum(mixture_var == 0) != 1L ||
      mixture_var[[1L]] != 0 || any(mixture_var[-1L] <= 0) ||
      anyDuplicated(mixture_var[-1L]) || is.unsorted(mixture_var[-1L], strictly = TRUE)) {
    stop("mixture_var must contain one leading zero and unique ascending positive multipliers.",
         call. = FALSE)
  }
  if (!is.logical(estimate_selection_s) || length(estimate_selection_s) != 1L ||
      is.na(estimate_selection_s)) stop("estimate_selection_s must be TRUE or FALSE.", call. = FALSE)
  if (estimate_selection_s) {
    stop("Sampled selection_s is not implemented for the joint MT prior; supply a fixed scalar selection_s.",
         call. = FALSE)
  }
  if (!is.null(selection_s_init) || !is.null(selection_s_prior) ||
      !is.null(selection_s_proposal_sd)) {
    stop("selection_s_init, selection_s_prior, and selection_s_proposal_sd require estimate_selection_s = TRUE, which is unsupported for MTBLR in Phase 19.",
         call. = FALSE)
  }
  selection_s_active <- !is.null(selection_s)
  if (selection_s_active) {
    if (!is.numeric(selection_s) || length(selection_s) != 1L ||
        !is.finite(selection_s)) stop("selection_s must be one finite scalar.", call. = FALSE)
    if (is.null(marker_frequency) || length(marker_frequency) != marker_count ||
        any(!is.finite(marker_frequency)) || any(marker_frequency <= 0 | marker_frequency >= 1)) {
      stop("Fixed selection_s requires aligned allele frequencies strictly inside (0, 1).",
           call. = FALSE)
    }
    marker_scale <- (2 * marker_frequency * (1 - marker_frequency)) ^
      (selection_s + 1)
  } else {
    selection_s <- NULL
    marker_scale <- rep(1, marker_count)
  }
  patterns <- pattern_spec$matrix
  null <- which(rowSums(patterns) == 0L)
  if (length(null) != 1L) stop("BayesR requires exactly one null trait pattern.", call. = FALSE)
  active <- which(rowSums(patterns) > 0L)
  positive <- mixture_var[-1L]
  state_pattern <- rbind(patterns[null, , drop = FALSE],
    patterns[rep(active, each = length(positive)), , drop = FALSE])
  state_pattern_index <- c(null, rep(active, each = length(positive)))
  state_component <- c(0L, rep(seq_along(positive), times = length(active)))
  state_multiplier <- c(0, rep(positive, times = length(active)))
  state_names <- c("null", unlist(lapply(active, function(index)
    paste0(pattern_spec$names[index], "__component_", seq_along(positive)))))
  if (nrow(state_pattern) > 4096L) {
    stop("The requested pattern-by-component state space exceeds 4096 joint states.",
         call. = FALSE)
  }
  if (is.null(joint_pi)) {
    joint_pi <- c(pattern_spec$probabilities[null], unlist(lapply(active,
      function(index) rep(pattern_spec$probabilities[index] / length(positive),
                          length(positive)))))
  }
  if (!is.numeric(joint_pi) || length(joint_pi) != nrow(state_pattern) ||
      any(!is.finite(joint_pi)) || any(joint_pi < 0) || sum(joint_pi) <= 0) {
    stop("joint_pi must be finite, nonnegative, and match the joint-state count.",
         call. = FALSE)
  }
  joint_pi <- as.numeric(joint_pi / sum(joint_pi))
  if (is.null(joint_pi_prior)) joint_pi_prior <- rep(1, length(joint_pi))
  if (!is.numeric(joint_pi_prior) || length(joint_pi_prior) != length(joint_pi) ||
      any(!is.finite(joint_pi_prior)) || any(joint_pi_prior <= 0)) {
    stop("joint_pi_prior must be finite, positive, and match the joint-state count.",
         call. = FALSE)
  }
  if (is.null(component)) component <- integer(marker_count)
  if (!is.numeric(component) || length(component) != marker_count ||
      any(!is.finite(component)) || any(component != as.integer(component)) ||
      any(component < 0L | component >= length(mixture_var))) {
    stop("component must contain one valid integer component index per marker.",
         call. = FALSE)
  }
  component <- as.integer(component)
  list(patterns = list(matrix = state_pattern,
      native = lapply(seq_len(nrow(state_pattern)), function(i)
        as.integer(state_pattern[i, ])), probabilities = joint_pi,
      names = state_names), method_code = 5L,
    joint_component = state_component,
    joint_multiplier = state_multiplier, joint_names = state_names,
    component_count = as.integer(length(mixture_var)),
    marker_scale = marker_scale, component_init = component,
    pi_prior = as.numeric(joint_pi_prior), selection_s = selection_s,
    mixture_var = mixture_var,
    model_parameters = list(mixture = list(
      component_multipliers = mixture_var,
      component_names = paste0("component_", seq_along(mixture_var) - 1L),
      trait_patterns = patterns, joint_state_names = state_names,
      joint_pattern_index = as.integer(state_pattern_index),
      joint_component_index = as.integer(state_component),
      selection_s = selection_s, selection_s_active = selection_s_active,
      selection_s_estimated = FALSE)))
}

.mtblr_bayesr_format_fit <- function(fit, model_parameters) {
  if (is.null(model_parameters)) return(fit)
  mixture <- model_parameters$mixture
  component_names <- mixture$component_names
  colnames(fit$component_probabilities) <- component_names
  if (!is.null(fit$chains)) {
    fit$chains <- lapply(fit$chains, function(chain) {
      if (!is.null(chain$marker$component_probabilities))
        colnames(chain$marker$component_probabilities) <- component_names
      chain
    })
  }
  pattern_count <- nrow(mixture$trait_patterns)
  component_count <- length(component_names)
  pi_final <- fit$pi_final %||% fit$pi
  pi_mean <- fit$pi_mean %||% fit$pim
  marginal <- function(probability, index, count) {
    out <- numeric(count)
    for (i in seq_along(probability)) out[index[[i]]] <-
      out[index[[i]]] + probability[[i]]
    out
  }
  mixture$pattern_pi_final <- marginal(
    pi_final, mixture$joint_pattern_index, pattern_count)
  mixture$pattern_pi_mean <- marginal(
    pi_mean, mixture$joint_pattern_index, pattern_count)
  mixture$component_pi_final <- marginal(
    pi_final, mixture$joint_component_index + 1L, component_count)
  mixture$component_pi_mean <- marginal(
    pi_mean, mixture$joint_component_index + 1L, component_count)
  names(mixture$component_pi_final) <- names(mixture$component_pi_mean) <-
    component_names
  model_parameters$mixture <- mixture
  fit$model_parameters <- model_parameters
  fit
}

.mtblr_bayesr_enrich_raw <- function(raw, method, model_parameters) {
  raw$meta$model <- method
  raw$meta$prior_kernel <- switch(
    method, bayesc = "bayesc", sbayesc = "bayesc",
    bayesr = "bayesr", sbayesr = "bayesr")
  raw$meta$model_semantics_version <- 2L
  raw$meta$model_semantics <- "s_prefix_means_summary_statistics"
  if (!is.null(model_parameters)) raw$model$mixture <- model_parameters$mixture
  raw
}

.mtblr_bayesr_frequency <- function(marker_metadata) {
  if (is.null(marker_metadata) || is.null(marker_metadata$allele_frequency)) NULL
  else as.numeric(marker_metadata$allele_frequency)
}

.mtblr_resolve_selection_maf <- function(
    selection_maf, selection_s_active, marker_count,
    summary_marker_metadata = NULL, analysis_frequency = NULL,
    reference_marker_metadata = NULL,
    allow_reference_maf_for_selection_s = FALSE) {
  if (!is.logical(allow_reference_maf_for_selection_s) ||
      length(allow_reference_maf_for_selection_s) != 1L ||
      is.na(allow_reference_maf_for_selection_s)) {
    stop("allow_reference_maf_for_selection_s must be TRUE or FALSE.",
         call. = FALSE)
  }
  validate <- function(x, source) {
    x <- as.numeric(x)
    if (length(x) != marker_count || any(!is.finite(x)) ||
        any(x <= 0 | x >= 1)) {
      stop(source, " selection MAF must be aligned, finite, length m, and strictly inside (0, 1).",
           call. = FALSE)
    }
    x
  }
  source <- "not_requested"
  population <- "not_applicable"
  fallback <- FALSE
  value <- NULL
  if (!isTRUE(selection_s_active) && is.null(selection_maf)) {
    return(list(
      values = NULL, selection_maf_source = source,
      selection_maf_population = population,
      selection_maf_alignment_status = "not_requested",
      selection_maf_fallback_used = FALSE))
  }
  if (!is.null(selection_maf)) {
    value <- validate(selection_maf, "Explicit")
    source <- "explicit_selection_maf"
    population <- "user_declared"
  } else if (!is.null(summary_marker_metadata) &&
             "allele_frequency" %in% names(summary_marker_metadata) &&
             all(is.finite(summary_marker_metadata$allele_frequency))) {
    value <- validate(summary_marker_metadata$allele_frequency, "GWAS-summary")
    source <- "gwas_summary_allele_frequency"
    population <- "gwas_summary_population"
  } else if (!is.null(analysis_frequency)) {
    value <- validate(analysis_frequency, "Analysis-genotype")
    source <- "analysis_genotype_frequency"
    population <- "analysis_sample"
  } else if (isTRUE(selection_s_active)) {
    reference <- .mtblr_bayesr_frequency(reference_marker_metadata)
    if (is.null(reference) || !isTRUE(allow_reference_maf_for_selection_s)) {
      stop(paste0(
        "selection_s requires explicit selection_maf or aligned GWAS-summary ",
        "allele frequency. Reference-panel MAF fallback is disabled unless ",
        "allow_reference_maf_for_selection_s = TRUE."), call. = FALSE)
    }
    value <- validate(reference, "Reference-panel")
    source <- "reference_panel_frequency"
    population <- "reference_panel"
    fallback <- TRUE
  }
  list(
    values = value, selection_maf_source = source,
    selection_maf_population = population,
    selection_maf_alignment_status = if (is.null(value)) "not_requested" else
      "aligned_to_final_marker_order",
    selection_maf_fallback_used = fallback)
}

.mtblr_bayesr_initialization <- function(beta, b, state, component,
                                          patterns, m, nt, method) {
  matrix_value <- function(value, name, mode = "double") {
    if (is.null(value)) return(NULL)
    if (is.list(value)) value <- do.call(cbind, value)
    value <- as.matrix(value)
    if (!identical(dim(value), c(m, nt)) || any(!is.finite(value)))
      stop(name, " must be a finite m by nt matrix or trait list.", call. = FALSE)
    if (mode == "integer" && any(!value %in% 0:1))
      stop(name, " must be binary.", call. = FALSE)
    if (mode == "integer") storage.mode(value) <- "integer"
    value
  }
  b <- matrix_value(b, "b") %||% matrix(0, m, nt)
  if (method %in% c("bayesc", "sbayesc")) {
    if (!is.null(beta) || !is.null(state) || length(component))
      stop("beta, state, and component initialization require BayesR or SBayesR.",
           call. = FALSE)
    return(list(beta = matrix(0, m, nt), b = b,
                state = matrix(0L, m, nt), component = integer(),
                policy = "phase18_bayesc_initialization"))
  }
  state <- matrix_value(state, "state", "integer")
  if (is.null(state)) {
    if (any(b != 0)) stop("Nonzero BayesR b initialization requires explicit state and component.",
                          call. = FALSE)
    state <- matrix(0L, m, nt)
  }
  beta <- matrix_value(beta, "beta")
  if (is.null(beta)) beta <- b
  component <- as.integer(component)
  null <- rowSums(state) == 0L
  if (any((component == 0L) != null))
    stop("component must be zero exactly for null marker patterns.", call. = FALSE)
  allowed <- apply(patterns, 1L, paste, collapse = "_")
  if (any(!apply(state, 1L, paste, collapse = "_") %in% allowed))
    stop("Every initialized state row must equal a supplied trait pattern.",
         call. = FALSE)
  if (any(b[state == 0L] != 0) || any(b != state * beta))
    stop("Initialized effective b must equal state * beta exactly.", call. = FALSE)
  list(beta = beta, b = b, state = state, component = component,
       policy = if (all(component == 0L)) "all_null" else "explicit_joint_state")
}

.mtblr_bayesr_memory <- function(memory, method, m, nchains, ncores,
                                 ntrace, state_count, component_count) {
  if (method %in% c("bayesc", "sbayesc")) {
    memory$bayesr_components_bytes <- 0
    return(memory)
  }
  workers <- min(as.integer(nchains), as.integer(ncores))
  shared <- 16 * state_count
  private <- workers * (4 * m + 16 * state_count)
  chain <- nchains * (4 * m + 8 * m * component_count +
                       8 * ntrace * state_count)
  formatted <- 4 * m + 8 * m * component_count
  additional <- shared + private + chain + formatted
  memory$bayesr_shared_state_descriptors_bytes <- shared
  memory$bayesr_private_worker_state_bytes <- private
  memory$bayesr_chain_result_bytes <- chain
  memory$bayesr_component_output_bytes <- formatted
  memory$bayesr_components_bytes <- additional
  memory$estimated_total_bytes <- memory$estimated_total_bytes + additional
  memory$estimated_total_gib <- memory$estimated_total_bytes / 1024^3
  if (!is.null(memory$execution_estimated_total_bytes)) {
    memory$execution_estimated_total_bytes <-
      memory$execution_estimated_total_bytes + additional
    memory$execution_estimated_total_gib <-
      memory$execution_estimated_total_bytes / 1024^3
  }
  memory
}
