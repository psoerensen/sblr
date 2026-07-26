phase17w_strict_lower <- function(nt) {
  stopifnot(length(nt) == 1L, is.finite(nt), nt == as.integer(nt), nt >= 1L)
  out <- vector("list", nt * (nt - 1L) / 2L)
  at <- 0L
  for (column in seq_len(nt)) {
    if (column == nt) next
    for (row in seq.int(column + 1L, nt)) {
      at <- at + 1L
      out[[at]] <- c(row = row, column = column)
    }
  }
  if (!length(out)) matrix(integer(), nrow = 0L, ncol = 2L,
                            dimnames = list(NULL, c("row", "column"))) else
    do.call(rbind, out)
}

phase17w_covariance_descriptors <- function(traits, residual = "full",
                                            updateB = TRUE, updateE = TRUE) {
  pairs <- phase17w_strict_lower(length(traits))
  groups <- c("B_cov", "G_cov", "E_cov")
  out <- do.call(rbind, lapply(groups, function(group) {
    if (!nrow(pairs)) return(NULL)
    structural <- group == "E_cov" && residual == "diagonal"
    updated <- switch(group, B_cov = updateB, G_cov = TRUE,
                      E_cov = updateE && !structural)
    status <- if (structural) "structural_zero" else
      if (!updated) "not_updated" else "eligible"
    data.frame(
      tier = "2A", group = group,
      trait = traits[pairs[, "row"]], trait2 = traits[pairs[, "column"]],
      row = pairs[, "row"], column = pairs[, "column"],
      updated = updated, captured = updated, derived = group == "G_cov",
      structural = structural, status = status,
      quantity = sprintf("%s[%s,%s]", substr(group, 1L, 1L),
                         traits[pairs[, "row"]], traits[pairs[, "column"]]),
      stringsAsFactors = FALSE)
  }))
  if (is.null(out)) data.frame() else out
}

phase17w_null_index <- function(models) {
  models <- as.matrix(models)
  hit <- which(rowSums(models) == 0L)
  if (length(hit) != 1L) stop("exactly one null pattern is required")
  hit
}

phase17w_resolve_patterns <- function(models, model_names, selection) {
  models <- as.matrix(models)
  stopifnot(length(model_names) == nrow(models), !anyDuplicated(model_names))
  if (is.null(selection)) return(integer())
  if (!length(selection) || anyNA(selection) || anyDuplicated(selection))
    stop("pattern selection must be nonempty, nonmissing, and unique")
  if (is.character(selection)) {
    hit <- match(selection, model_names)
    if (anyNA(hit)) stop("unknown model name")
    return(as.integer(hit))
  }
  if (!is.numeric(selection) || any(!is.finite(selection)) ||
      any(selection != as.integer(selection)) ||
      any(selection < 1L | selection > nrow(models)))
    stop("model indices must be valid one-based indices")
  as.integer(selection)
}

phase17w_probability_plan <- function(models, model_names, probability = "mass",
                                      selection = NULL, updatePi = TRUE) {
  probability <- match.arg(probability, c("none", "mass", "selected", "all"))
  null <- phase17w_null_index(models)
  selected <- phase17w_resolve_patterns(models, model_names, selection)
  if (probability == "selected" && !length(selected))
    stop("selected probability scope requires models")
  if (probability != "selected" && length(selected))
    stop("model selection is valid only for selected probability scope")
  physical <- switch(probability, none = integer(), mass = null,
                     selected = selected, all = seq_len(nrow(models)))
  list(
    probability = probability, null_index = null,
    primary_mass_quantity = if (probability == "mass") "pi_active" else NULL,
    complement = if (probability == "mass") "pi_null" else NULL,
    diagnostic_key = if (probability == "mass") "pi_mass:null_active" else
      paste0("pi_pattern:", physical),
    physical_model_indices = unique(physical),
    updated = updatePi, captured = updatePi && length(physical) > 0L)
}

phase17w_resolve_markers <- function(marker_ids, marker_ids_request = NULL,
                                     marker_indices = NULL) {
  if (!is.null(marker_ids_request) && !is.null(marker_indices))
    stop("choose marker IDs or marker indices, not both")
  if (is.null(marker_ids_request) && is.null(marker_indices)) return(integer())
  if (!is.null(marker_ids_request)) {
    if (!is.character(marker_ids_request) || !length(marker_ids_request) ||
        anyNA(marker_ids_request) || any(!nzchar(marker_ids_request)) ||
        anyDuplicated(marker_ids_request)) stop("invalid marker IDs")
    hits <- lapply(marker_ids_request, function(id) which(marker_ids == id))
    if (any(lengths(hits) == 0L)) stop("unknown marker ID")
    if (any(lengths(hits) > 1L)) stop("ambiguous marker ID")
    return(as.integer(unlist(hits, use.names = FALSE)))
  }
  if (!is.numeric(marker_indices) || !length(marker_indices) ||
      anyNA(marker_indices) || any(!is.finite(marker_indices)) ||
      any(marker_indices != as.integer(marker_indices)) ||
      any(marker_indices < 1L | marker_indices > length(marker_ids)) ||
      anyDuplicated(marker_indices)) stop("invalid marker indices")
  as.integer(marker_indices)
}

phase17w_extended_memory <- function(C, N, T, updateB = TRUE, updateE = TRUE,
                                     updatePi = TRUE, residual = "full",
                                     probability = "mass", P = 2L,
                                     selected_patterns = integer(),
                                     null_index = 1L, K = 0L,
                                     marker_quantities = c("b", "d"),
                                     keep_traces = FALSE) {
  qoff <- T * (T - 1) / 2
  b_cov <- as.numeric(updateB) * 8 * C * N * qoff
  g_cov <- 8 * C * N * qoff
  e_cov <- as.numeric(updateE && residual == "full") * 8 * C * N * qoff
  models <- switch(probability,
    none = integer(), mass = null_index,
    selected = unique(as.integer(selected_patterns)), all = seq_len(P))
  probability_bytes <- as.numeric(updatePi) * 8 * C * N * length(unique(models))
  mass <- as.numeric(updatePi && probability == "mass") * 8 * C * N
  selected <- as.numeric(updatePi && probability == "selected") *
    8 * C * N * length(unique(selected_patterns))
  full <- as.numeric(updatePi && probability == "all") * 8 * C * N * P
  marker_b <- as.numeric("b" %in% marker_quantities) * 8 * C * N * K * T
  marker_d <- as.numeric("d" %in% marker_quantities) * 4 * C * N * K * T
  capture <- b_cov + g_cov + e_cov + probability_bytes + marker_b + marker_d
  c(B_cov_bytes = b_cov, G_cov_bytes = g_cov, E_cov_bytes = e_cov,
    pi_mass_bytes = mass, selected_pi_bytes = selected,
    full_pi_bytes = full, probability_unique_bytes = probability_bytes,
    selected_b_bytes = marker_b, selected_d_bytes = marker_d,
    captured_trace_bytes = capture,
    retained_trace_bytes = if (keep_traces) capture else 0,
    workspace_bytes = 8 * C * N * 8)
}

phase17w_quantity_counts <- function(T, probability_rows = 0L, K = 0L,
                                     marker_quantities = c("b", "d")) {
  qoff <- T * (T - 1L) / 2L
  as.integer(c(tier1 = 3L * T, covariance = 3L * qoff,
    probability = probability_rows,
    marker_b = as.integer("b" %in% marker_quantities) * K * T,
    marker_d = as.integer("d" %in% marker_quantities) * K * T)) |>
    stats::setNames(c("tier1", "covariance", "probability", "marker_b", "marker_d"))
}

phase17w_resolved_limit_gb <- function(memory_warning_gb = 8,
                                       trace_limit_gb = NULL) {
  if (!is.null(trace_limit_gb)) return(trace_limit_gb)
  if (!is.finite(memory_warning_gb)) return(2)
  min(2, max(.25, .25 * memory_warning_gb))
}

phase17w_large_request <- function(capture_bytes, memory_warning_gb = 8,
                                   trace_limit_gb = NULL,
                                   allow_large_traces = FALSE) {
  limit <- phase17w_resolved_limit_gb(memory_warning_gb, trace_limit_gb)
  exceeds <- capture_bytes > limit * 1024^3
  list(limit_gb = limit, exceeds = exceeds,
       decision = if (!exceeds || allow_large_traces) "allow" else
         "require_explicit_override")
}

phase17w_validate_future_control <- function(convergence = "extended",
                                             extended = list()) {
  if (!identical(convergence, "extended") && length(extended))
    stop("extended controls require convergence='extended'")
  defaults <- list(covariance = "off_diagonal", probability = "mass",
    probability_models = NULL, marker_ids = NULL, marker_indices = NULL,
    marker_quantities = c("b", "d"), allow_large_traces = FALSE,
    trace_limit_gb = NULL)
  if (!is.list(extended) || is.data.frame(extended) ||
      (length(extended) && (is.null(names(extended)) || any(!nzchar(names(extended))) ||
                            anyDuplicated(names(extended))))) stop("invalid extended control")
  if (any(!names(extended) %in% names(defaults))) stop("unknown extended control")
  supplied <- names(extended)
  out <- utils::modifyList(defaults, extended)
  out$covariance <- match.arg(out$covariance, c("off_diagonal", "none"))
  out$probability <- match.arg(out$probability, c("mass", "none", "selected", "all"))
  if (out$probability == "selected" && is.null(out$probability_models))
    stop("selected probability scope requires probability_models")
  if (out$probability != "selected" && !is.null(out$probability_models))
    stop("probability_models conflicts with probability scope")
  if (!is.null(out$marker_ids) && !is.null(out$marker_indices))
    stop("marker selection representations conflict")
  selected <- !is.null(out$marker_ids) || !is.null(out$marker_indices)
  if (!selected && "marker_quantities" %in% supplied)
    stop("marker quantities require marker selection")
  if (selected && (!length(out$marker_quantities) ||
      any(!out$marker_quantities %in% c("b", "d")) ||
      anyDuplicated(out$marker_quantities))) stop("invalid marker quantities")
  if (!is.logical(out$allow_large_traces) || length(out$allow_large_traces) != 1L ||
      is.na(out$allow_large_traces)) stop("allow_large_traces must be logical")
  if (!is.null(out$trace_limit_gb) &&
      (length(out$trace_limit_gb) != 1L || !is.numeric(out$trace_limit_gb) ||
       !is.finite(out$trace_limit_gb) || out$trace_limit_gb <= 0))
    stop("trace_limit_gb must be finite and positive")
  out
}

phase17w_warning_groups <- function(summary) {
  flagged <- with(summary, rhat_flag | ess_bulk_flag | ess_tail_flag | mcse_flag |
                    status == "constant_chain_mismatch")
  table(factor(summary$group[flagged], levels = unique(summary$group)))
}
