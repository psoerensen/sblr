.mtblr_bayesrc_empty <- function() list(
  annotations = matrix(numeric(), 0L, 0L),
  alpha_init = matrix(numeric(), 0L, 0L), sigma_alpha_init = numeric(),
  pattern_pi_init = numeric(), pattern_pi_prior = numeric(),
  updateAlpha = FALSE, intercept_flat = TRUE, sigma_alpha_a = 2,
  sigma_alpha_b = 2, pi_floor = 1e-12, alpha_update_every = 1L,
  metadata = NULL, model_parameters = NULL, maf_annotation_overlap = FALSE)

.mtblr_bayesrc_prior_probabilities <- function(annotations, alpha,
                                                pi_floor = 1e-12) {
  eta <- annotations %*% alpha
  stick <- stats::pnorm(eta)
  out <- matrix(0, nrow(annotations), ncol(alpha) + 1L)
  remaining <- rep(1, nrow(annotations))
  for (k in seq_len(ncol(alpha))) {
    out[, k] <- remaining * (1 - stick[, k])
    remaining <- remaining * stick[, k]
  }
  out[, ncol(out)] <- remaining
  out <- pmax(out, pi_floor)
  out / rowSums(out)
}

.mtblr_calibration_inputs <- function(mixture, bayesrc, marker_count) {
 patterns <- mixture$patterns$matrix
 state_count <- nrow(patterns)
 gamma <- if (length(mixture$joint_multiplier) == state_count) {
  as.numeric(mixture$joint_multiplier)
 } else rep(1, state_count)
 marker_scale <- if (length(mixture$marker_scale) == marker_count) {
  as.numeric(mixture$marker_scale)
 } else rep(1, marker_count)
 if (is.null(bayesrc$model_parameters)) {
  initial <- as.numeric(mixture$patterns$probabilities)
  prior <- if (length(mixture$pi_prior) == state_count) {
   as.numeric(mixture$pi_prior / sum(mixture$pi_prior))
  } else initial
  return(list(
   patterns = patterns, initial = initial, prior = prior, gamma = gamma,
   marker_scale = marker_scale,
   component_probability_source = if (length(mixture$pi_prior))
    "joint_initial_and_dirichlet_prior_mean" else "joint_pattern_initial",
   annotation_probability_policy = "not_applicable"))
 }
 component_probability <- .mtblr_bayesrc_prior_probabilities(
  bayesrc$annotations, bayesrc$alpha_init, bayesrc$pi_floor)
 joint <- matrix(0, marker_count, state_count)
 pattern_index <- mixture$model_parameters$mixture$joint_pattern_index
 component_index <- mixture$model_parameters$mixture$joint_component_index
 active_patterns <- which(rowSums(
  mixture$model_parameters$mixture$trait_patterns) > 0L)
 for (s in seq_len(state_count)) {
  component <- component_index[s]
  if (component == 0L) {
   joint[, s] <- component_probability[, 1L]
  } else {
   conditional_index <- match(pattern_index[s], active_patterns)
   joint[, s] <- component_probability[, component + 1L] *
    bayesrc$pattern_pi_init[conditional_index]
  }
 }
 joint <- joint / rowSums(joint)
 list(
  patterns = patterns, initial = joint, prior = joint, gamma = gamma,
  marker_scale = marker_scale,
  component_probability_source = "resolved_marker_joint_state_probabilities",
  annotation_probability_policy = "resolved_initial_annotation_probabilities")
}

.mtblr_bayesrc_controls <- function(
    prior_kernel, annotations, marker_ids, pattern_spec, mixture,
    add_intercept = TRUE, standardize_annotations = TRUE,
    center_binary_annotations = FALSE, alpha_init = NULL,
    sigmaSqAlpha_init = NULL, intercept_flat = TRUE,
    sigmaSqAlpha_a = 2, sigmaSqAlpha_b = 2, pi_floor = 1e-12,
    alpha_update_every = 1L, updateAlpha = TRUE) {
  supplied <- list(annotations, alpha_init, sigmaSqAlpha_init)
  if (!identical(prior_kernel, "bayesrc")) {
    if (any(!vapply(supplied, is.null, logical(1))) ||
        !identical(add_intercept, TRUE) ||
        !identical(standardize_annotations, TRUE) ||
        !identical(center_binary_annotations, FALSE) ||
        !identical(intercept_flat, TRUE) || !identical(sigmaSqAlpha_a, 2) ||
        !identical(sigmaSqAlpha_b, 2) || !identical(pi_floor, 1e-12) ||
        !identical(as.integer(alpha_update_every), 1L) ||
        !identical(updateAlpha, TRUE)) {
      stop("Annotation controls require method = 'bayesrc' or 'sbayesrc'.",
           call. = FALSE)
    }
    return(.mtblr_bayesrc_empty())
  }
  if (is.null(annotations))
    stop("annotations is required for MT BayesRC/SBayesRC.", call. = FALSE)
  for (name in c("add_intercept", "standardize_annotations",
                 "center_binary_annotations", "intercept_flat", "updateAlpha")) {
    value <- get(name)
    if (!is.logical(value) || length(value) != 1L || is.na(value))
      stop(name, " must be TRUE or FALSE.", call. = FALSE)
  }
  for (name in c("sigmaSqAlpha_a", "sigmaSqAlpha_b", "pi_floor")) {
    value <- get(name)
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
        value <= 0 || (name == "pi_floor" && value >= 0.5))
      stop(name, " must be a valid finite positive scalar.", call. = FALSE)
  }
  alpha_update_every <- .blr_scalar_integer(
    alpha_update_every, "alpha_update_every")
  annotation_ids <- NULL
  if (is.data.frame(annotations)) {
    id_columns <- names(annotations)[tolower(names(annotations)) %in%
      c("marker_id", "marker", "rsid", "rsids")]
    if (length(id_columns) == 1L) annotation_ids <- annotations[[id_columns]]
  }
  if (is.null(annotation_ids)) {
    annotation_rows <- rownames(annotations)
    default_rows <- !is.null(annotation_rows) &&
      identical(annotation_rows, as.character(seq_len(nrow(annotations))))
    if (!is.null(annotation_rows) && !default_rows)
      annotation_ids <- annotation_rows
  }
  if (is.null(annotation_ids)) {
    stop("MT BayesRC annotations require explicit marker IDs in row names or one marker-ID column; external row order is never trusted.",
         call. = FALSE)
  }
  info <- .stblr_align_bed_bayesrc_annotations(
    annotations, marker_ids, add_intercept, standardize_annotations,
    center_binary_annotations)
  if (info$alignment$unused_annotation_rows != 0L) {
    stop("annotations contains marker IDs outside the final selected marker set; extra rows are not silently dropped.",
         call. = FALSE)
  }
  A <- info$A
  patterns <- pattern_spec$matrix
  active <- which(rowSums(patterns) > 0L)
  omega <- as.numeric(pattern_spec$probabilities[active])
  omega <- omega / sum(omega)
  component_probability <- numeric(mixture$component_count)
  for (state in seq_along(mixture$patterns$probabilities)) {
    index <- mixture$joint_component[[state]] + 1L
    component_probability[[index]] <- component_probability[[index]] +
      mixture$patterns$probabilities[[state]]
  }
  initialized <- .stblr_initialize_bed_bayesrc_prior(
    A, mixture$mixture_var, component_probability,
    annot_alpha_init = alpha_init,
    annot_sigma_sq_alpha_init = sigmaSqAlpha_init,
    pi_floor = pi_floor)
  maf_names <- c("maf", "allele_frequency", "heterozygosity",
                 "log_heterozygosity")
  overlap <- any(tolower(colnames(A)) %in% maf_names)
  component_names <- paste0("component_", seq_len(mixture$component_count) - 1L)
  metadata <- list(
    annotation_policy = "annotation_probit_stick",
    annotation_source = "public_annotations",
    annotation_marker_alignment_status = if (isTRUE(info$alignment$matched_by_id))
      "exact_marker_id_match" else "by_construction_order",
    annotation_column_names = colnames(A),
    processed_annotation_names = colnames(A),
    annotation_column_types = ifelse(vapply(seq_len(ncol(A)), function(j)
      all(A[, j] %in% c(0, 1)), logical(1)), "binary", "continuous"),
    annotation_processed_column_means = stats::setNames(
      colMeans(A), colnames(A)),
    annotation_processed_column_sds = stats::setNames(vapply(
      seq_len(ncol(A)), function(j) stats::sd(A[, j]), numeric(1)),
      colnames(A)),
    annotation_standardization_policy = if (standardize_annotations)
      "continuous_z_score" else "none",
    annotation_intercept_policy = if (info$preprocessing$intercept_added)
      "added_once" else "preserved_or_not_requested",
    annotation_preprocessing = c(info$preprocessing, info$alignment),
    maf_annotation_overlap = overlap)
  list(
    annotations = A, alpha_init = initialized$annot_alpha_init,
    sigma_alpha_init = initialized$annot_sigma_sq_alpha_init,
    pattern_pi_init = omega, pattern_pi_prior = rep(1, length(omega)),
    updateAlpha = updateAlpha, intercept_flat = intercept_flat,
    sigma_alpha_a = sigmaSqAlpha_a, sigma_alpha_b = sigmaSqAlpha_b,
    pi_floor = pi_floor, alpha_update_every = alpha_update_every,
    metadata = metadata, maf_annotation_overlap = overlap,
    model_parameters = list(annotations = c(metadata, list(
      component_names = component_names,
      stick_names = paste0("step_", seq_len(mixture$component_count - 1L)),
      pattern_names = pattern_spec$names[active]))))
}

.mtblr_bayesrc_enrich_raw <- function(raw, bayesrc, method, updatePi) {
  if (is.null(bayesrc$model_parameters)) return(raw)
  if (is.null(raw$annotations))
    stop("Native MT BayesRC result omitted annotations.", call. = FALSE)
  raw$annotations$policy <- "annotation_probit_stick"
  raw$annotations$metadata <- bayesrc$metadata
  raw$annotations$coefficient_jitter_attempts <- 0L
  raw$annotations$coefficient_max_jitter <- 0
  raw$annotations$coefficient_status <- if (bayesrc$updateAlpha)
    "computed" else "not_updated"
  raw$annotations$pattern_probability_status <- if (isTRUE(updatePi))
    "computed" else "not_updated"
  raw$meta$model <- method
  raw$meta$prior_kernel <- "bayesrc"
  raw$model$annotations <- bayesrc$model_parameters$annotations
  raw$pi$final <- raw$pi$mean <- raw$pi$trace <- NULL
  raw
}

.mtblr_bayesrc_format_fit <- function(fit, raw_annotations, bayesrc) {
  if (is.null(bayesrc$model_parameters)) return(fit)
  parameters <- bayesrc$model_parameters
  parameters$annotations <- c(parameters$annotations, raw_annotations[c(
    "annotation_coefficients_final", "annotation_coefficients_mean",
    "annotation_variances_final", "annotation_variances_mean",
    "pattern_pi_final", "pattern_pi_mean", "pattern_pi_trace",
    "prior_component_probabilities", "annotation_updates_attempted",
    "annotation_updates_completed", "coefficient_jitter_attempts",
    "coefficient_max_jitter", "coefficient_status",
    "pattern_probability_status")])
  parameters$annotations$prior_active_probabilities <-
    1 - parameters$annotations$prior_component_probabilities[, 1L]
  fit$model_parameters <- utils::modifyList(fit$model_parameters %||% list(),
                                             parameters)
  fit["pi_final"] <- list(NULL)
  fit["pi_mean"] <- list(NULL)
  fit["pi_trace"] <- list(NULL)
  fit
}

.mtblr_bayesrc_memory <- function(memory, bayesrc, marker_count, nchains,
                                  worker_count, total_iterations) {
  if (is.null(bayesrc$model_parameters)) {
    additions <- list(
      bayesrc_annotations_requested = FALSE,
      bayesrc_shared_annotation_bytes = 0,
      bayesrc_private_annotation_bytes_per_worker = 0,
      bayesrc_chain_annotation_result_bytes = 0,
      bayesrc_formatted_annotation_output_bytes = 0,
      bayesrc_estimated_total_bytes = 0,
      bayesrc_estimated_total_gib = 0)
  } else {
    q <- ncol(bayesrc$annotations)
    components <- ncol(bayesrc$alpha_init) + 1L
    shared <- 8 * marker_count * q
    private <- 8 * (q * (components - 1L) + (components - 1L) +
      marker_count * (components - 1L) + marker_count * components)
    chain_result <- nchains * 8 * (2 * q * (components - 1L) +
      2 * (components - 1L) + 2 * marker_count * components +
      2 * total_iterations)
    formatted <- 8 * marker_count * components
    total <- shared + min(nchains, worker_count) * private + chain_result +
      formatted
    additions <- list(
      bayesrc_annotations_requested = TRUE,
      bayesrc_shared_annotation_bytes = shared,
      bayesrc_private_annotation_bytes_per_worker = private,
      bayesrc_chain_annotation_result_bytes = chain_result,
      bayesrc_formatted_annotation_output_bytes = formatted,
      bayesrc_estimated_total_bytes = total,
      bayesrc_estimated_total_gib = total / 1024^3)
    memory$estimated_total_bytes <- memory$estimated_total_bytes + total
    memory$estimated_total_gib <- memory$estimated_total_bytes / 1024^3
    memory$execution_estimated_total_bytes <-
      memory$execution_estimated_total_bytes + total
    memory$execution_estimated_total_gib <-
      memory$execution_estimated_total_bytes / 1024^3
  }
  c(memory, additions)
}
