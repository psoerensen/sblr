.blr_fit_raw_v2 <- function(fit) {
  raw <- attr(fit, "blr_raw", exact = TRUE)
  if (is.null(raw)) return(NULL)
  if (!inherits(raw, "blr_raw") || !identical(raw$schema$version, 2L)) {
    stop("The fit carries an invalid raw-v2 object.", call. = FALSE)
  }
  raw
}

.blr_fit_quantity <- function(fit, quantity, state) {
  raw <- .blr_fit_raw_v2(fit)
  formatted <- list(
    realised_effects = "bm", pips = "dm",
    component_probabilities = "component_probabilities",
    activity_pattern_probabilities = "activity_pattern_probabilities",
    pleiotropic_probabilities = "pleiotropic_probabilities",
    effect_covariance = "cov_b_mean",
    residual_covariance = "cov_e_mean",
    predictions = "predictions",
    pip_chain_mean_sd = "dm_chain_mean_sd",
    effect_chain_mean_sd = "bm_chain_mean_sd",
    mean_component_assignment = "dm_component_mean")
  raw_name <- list(
    retained = c(
      realised_effects = "realised_effects",
      latent_effects = "latent_effects",
      states = "independent_trait_states",
      joint_states = "joint_states",
      component_assignments = "component_assignments",
      effect_covariance = "marker_covariance",
      residual_covariance = "residual_covariance",
      predictions = "predictions",
      activity_probability_parameters = "activity_pattern_parameters",
      scale_probability_parameters = "joint_component_probability_parameters"),
    final = c(
      realised_effects = "realised_effects",
      latent_effects = "latent_effects",
      states = "independent_trait_states",
      joint_states = "joint_states",
      component_assignments = "component_assignments",
      effect_covariance = "marker_covariance",
      residual_covariance = "residual_covariance",
      activity_probability_parameters = "activity_pattern_parameters",
      scale_probability_parameters = "joint_component_probability_parameters"))
  if (identical(state, "posterior")) {
    field <- formatted[[quantity]]
    if (is.null(field)) {
      stop("quantity = '", quantity,
           "' is not a stored posterior fit quantity.", call. = FALSE)
    }
    return(fit[[field]])
  }
  if (is.null(raw)) {
    stop("Retained and final extraction requires a fit carrying raw-v2.",
         call. = FALSE)
  }
  field <- unname(raw_name[[state]][quantity])
  if (is.na(field)) {
    stop("quantity = '", quantity, "' is unavailable for state = '", state,
         "'.", call. = FALSE)
  }
  raw[[if (identical(state, "retained")) "draws" else "final"]][[field]]
}

.blr_quantity_axes <- function(quantity, state, value) {
  if (is.null(value) || is.list(value) && !is.array(value)) return(NULL)
  if (identical(state, "posterior")) {
    return(switch(quantity,
      realised_effects = c("marker", "trait"),
      pips = c("marker", "trait"),
      component_probabilities = c("marker", "component"),
      activity_pattern_probabilities = c("marker", "activity_pattern"),
      pleiotropic_probabilities = "marker",
      effect_covariance = c("trait_row", "trait_col"),
      residual_covariance = c("trait_row", "trait_col"),
      predictions = c("observation", "trait"),
      pip_chain_mean_sd = c("marker", "trait"),
      effect_chain_mean_sd = c("marker", "trait"),
      mean_component_assignment = c("marker", "trait"), NULL))
  }
  prefix <- if (identical(state, "retained")) c("draw", "chain") else
    "chain"
  switch(quantity,
    realised_effects = c(prefix, "marker", "trait"),
    latent_effects = c(prefix, "marker", "trait"),
    states = c(prefix, "marker", "trait"),
    joint_states = c(prefix, "marker"),
    component_assignments = c(prefix, "marker"),
    effect_covariance = c(prefix, "trait_row", "trait_col"),
    residual_covariance = c(prefix, "trait_row", "trait_col"),
    predictions = c(prefix, "observation", "trait"),
    activity_probability_parameters = c(prefix, "activity_pattern"),
    scale_probability_parameters = c(prefix, "component"), NULL)
}

.blr_select_axis <- function(value, axes, axis, selection) {
  if (is.null(selection) || is.null(value)) return(value)
  position <- match(axis, axes)
  if (is.na(position)) {
    stop("The extracted quantity has no ", axis, " axis.", call. = FALSE)
  }
  ids <- dimnames(value)[[position]]
  index <- if (is.character(selection)) {
    if (is.null(ids) || anyNA(match(selection, ids))) {
      stop("Unknown ", axis, " identifier.", call. = FALSE)
    }
    match(selection, ids)
  } else {
    index <- as.integer(selection)
    if (!length(index) || anyNA(index) || any(index < 1L) ||
        any(index > dim(value)[[position]])) {
      stop(axis, " selection is outside the declared axis.", call. = FALSE)
    }
    index
  }
  subscripts <- rep(list(quote(expr = )), length(dim(value)))
  subscripts[[position]] <- index
  do.call(`[`, c(list(value), subscripts, list(drop = FALSE)))
}

#' Extract a canonical posterior quantity
#'
#' Retrieves a declared scientific quantity from a formatted STBLR or MTBLR
#' fit. Extraction never computes a posterior summary and never drops array
#' dimensions. A quantity that is scientifically unavailable is returned as
#' `NULL` when its canonical field is present-but-unavailable.
#'
#' @param fit A formatted `stblr_fit` or `mtblr_fit`.
#' @param quantity Scientific quantity name. Common posterior quantities are
#'   `"realised_effects"`, `"pips"`, `"component_probabilities"`,
#'   `"effect_covariance"`, `"residual_covariance"`, and `"predictions"`.
#' @param state One of `"posterior"`, `"retained"`, or `"final"`.
#' @param markers,traits,chains,draws,components,activity_patterns Optional
#'   integer positions or exact IDs on the corresponding axis.
#' @return The requested object with all declared axes retained, or `NULL`.
#' @export
extract_posterior <- function(
    fit, quantity, state = c("posterior", "retained", "final"),
    markers = NULL, traits = NULL, chains = NULL, draws = NULL,
    components = NULL, activity_patterns = NULL) {
  if (!inherits(fit, c("stblr_fit", "mtblr_fit", "blr_fit"))) {
    stop("fit must inherit from stblr_fit, mtblr_fit, or blr_fit.",
         call. = FALSE)
  }
  quantity <- .blr_character_scalar(quantity, "quantity")
  state <- match.arg(state)
  value <- .blr_fit_quantity(fit, quantity, state)
  if (is.null(value)) return(NULL)
  if (is.list(value) && !is.array(value)) {
    if (!is.null(traits)) {
      ids <- names(value)
      index <- if (is.character(traits)) match(traits, ids) else
        as.integer(traits)
      if (!length(index) || anyNA(index) || any(index < 1L) ||
          any(index > length(value))) {
        stop("Unknown trait selection.", call. = FALSE)
      }
      value <- value[index]
    }
    return(value)
  }
  if (is.null(dim(value))) {
    value <- array(value, dim = length(value), dimnames = list(names(value)))
  }
  axes <- .blr_quantity_axes(quantity, state, value)
  if (!is.null(traits) && !"trait" %in% axes &&
      all(c("trait_row", "trait_col") %in% axes)) {
    value <- .blr_select_axis(value, axes, "trait_row", traits)
    value <- .blr_select_axis(value, axes, "trait_col", traits)
    traits <- NULL
  }
  selections <- list(
    marker = markers, trait = traits, chain = chains, draw = draws,
    component = components, activity_pattern = activity_patterns)
  for (axis in names(selections)) {
    value <- .blr_select_axis(value, axes, axis, selections[[axis]])
  }
  value
}

#' Extract fit diagnostics
#'
#' Returns one canonical read-only view of sampler, execution, convergence,
#' provider, and provenance diagnostics. The view assembles fields already
#' stored in the validated raw-v2 object or formatted fit; it does not
#' recompute scientific quantities.
#'
#' @param fit A formatted `stblr_fit` or `mtblr_fit`.
#' @param quantity Optional exact diagnostic namespace: `"sampler"`,
#'   `"execution"`, `"convergence"`, `"providers"`, or `"provenance"`.
#' @return The complete diagnostic view, one named namespace, or `NULL` for an
#'   unavailable value within a returned namespace.
#' @export
extract_diagnostics <- function(fit, quantity = NULL) {
  if (!inherits(fit, c("stblr_fit", "mtblr_fit", "blr_fit"))) {
    stop("fit must inherit from stblr_fit, mtblr_fit, or blr_fit.",
         call. = FALSE)
  }
  raw <- .blr_fit_raw_v2(fit)
  raw_input <- raw$input %||% list()
  raw_mcmc <- raw_input$mcmc %||% list()
  raw_compute <- raw_input$compute %||% list()
  raw_data <- raw_input$data %||% list()
  sampler <- raw$diagnostics %||% fit$diagnostics %||% list()
  workers <- sampler$workers %||% fit$diagnostics$workers %||%
    fit$diagnostics$worker %||% NULL
  provenance <- raw$provenance %||% list(
    operator_resources = fit$data$operator_resources %||% NULL,
    marker_alignment = fit$data$alignment %||% NULL,
    seed_contract_version = fit$input$seed_contract_version %||% NULL,
    task_seeds = fit$input$task_seeds_resolved %||% NULL)
  diagnostics <- list(
    sampler = sampler,
    execution = list(
      task_seeds = raw_mcmc$task_seeds %||%
        fit$input$task_seeds_resolved %||% provenance$task_seeds %||% NULL,
      logical_task_ids = fit$input$logical_task_ids %||%
        workers$logical_task_order %||% NULL,
      retained_transition_indices =
        raw_mcmc$retained_transition_indices %||%
          fit$input$retained_transition_indices %||% NULL,
      convergence_iteration_indices =
        fit$input$convergence_iteration_indices %||%
          raw$draws$convergence$transition_indices %||% NULL,
      seed_contract_version = provenance$seed_contract_version %||%
        fit$input$seed_contract_version %||% NULL,
      retention_contract_version =
        raw_input$schema$retention_contract_version %||%
          fit$input$retention_contract_version %||% NULL,
      execution_contract_version =
        raw_input$schema$execution_contract_version %||% NULL,
      scheduler_version = raw_compute$scheduler_version %||%
        workers$scheduler_version %||% fit$input$scheduler_version %||% NULL,
      workers = workers),
    convergence = list(
      result = fit$convergence %||% sampler$convergence %||% NULL,
      traces = fit$convergence_traces %||%
        raw$draws$convergence %||% NULL),
    providers = list(
      providers = raw_data$providers %||% fit$data$providers %||% NULL,
      operator_resources = raw_data$operator_resources %||%
        fit$data$operator_resources %||% NULL,
      provider_maps = raw_data$provider_maps %||%
        provenance$marker_alignment %||% NULL,
      likelihood_regime = raw_data$likelihood_regime %||% NULL),
    provenance = provenance)
  if (is.null(quantity)) return(diagnostics)
  quantity <- .blr_character_scalar(quantity, "quantity")
  if (is.null(names(diagnostics)) || !quantity %in% names(diagnostics)) {
    stop("Unknown diagnostic namespace: ", quantity, ".", call. = FALSE)
  }
  diagnostics[[quantity]]
}

.blr_summarise_retained <- function(fit, quantity, prob, ...) {
  value <- extract_posterior(fit, quantity, state = "retained", ...)
  if (is.null(value)) return(NULL)
  if (!is.array(value) || length(dim(value)) < 2L) {
    stop("The retained quantity must have draw and chain axes.", call. = FALSE)
  }
  if (!is.numeric(prob) || length(prob) != 1L || !is.finite(prob) ||
      prob <= 0 || prob >= 1) {
    stop("prob must be a finite scalar in (0, 1).", call. = FALSE)
  }
  sample_count <- prod(dim(value)[1:2])
  cell_dims <- dim(value)[-(1:2)]
  cell_names <- dimnames(value)[-(1:2)]
  axes <- .blr_quantity_axes(quantity, "retained", value)[-(1:2)]
  matrix_value <- matrix(value, nrow = sample_count)
  alpha <- (1 - prob) / 2
  statistics <- t(vapply(seq_len(ncol(matrix_value)), function(column) {
    x <- matrix_value[, column]
    x <- x[is.finite(x)]
    if (!length(x)) return(rep(NA_real_, 6L))
    interval <- stats::quantile(
      x, c(alpha, 1 - alpha), names = FALSE, type = 8)
    c(length(x), mean(x), stats::sd(x), stats::median(x), interval)
  }, numeric(6L)))
  colnames(statistics) <- c(
    "n", "mean", "sd", "median", "q_lower", "q_upper")
  if (!length(cell_dims)) return(as.data.frame(statistics))
  index <- do.call(expand.grid, c(lapply(cell_dims, seq_len),
                                 KEEP.OUT.ATTRS = FALSE,
                                 stringsAsFactors = FALSE))
  names(index) <- axes
  for (axis in seq_along(axes)) {
    ids <- cell_names[[axis]]
    if (!is.null(ids)) index[[axis]] <- ids[index[[axis]]]
  }
  cbind(index, as.data.frame(statistics), row.names = NULL)
}
