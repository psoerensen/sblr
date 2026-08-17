# ---- consolidated from tests/testthat/helper-mtblr-bed-convergence-contract.R ----
blr_convergence_contract_postburn <- function(chains, nburn) {
  stopifnot(is.matrix(chains), nburn >= 0, nburn < ncol(chains))
  chains[, seq.int(nburn + 1L, ncol(chains)), drop = FALSE]
}

blr_convergence_contract_split_chains <- function(chains) {
  stopifnot(is.matrix(chains))
  half <- ncol(chains) %/% 2L
  if (half < 1L) return(matrix(numeric(), 2L * nrow(chains), 0L))
  rbind(chains[, seq_len(half), drop = FALSE],
        chains[, ncol(chains) - half + seq_len(half), drop = FALSE])
}

blr_convergence_contract_rank_normalize <- function(x) {
  ranks <- rank(as.vector(x), ties.method = "average")
  stats::qnorm((ranks - 3 / 8) / (length(ranks) + 1 / 4))
}

blr_convergence_contract_rhat_basic <- function(split) {
  m <- nrow(split); n <- ncol(split)
  within <- apply(split, 1L, stats::var)
  w <- mean(within)
  if (!is.finite(w) || w <= 0) return(NA_real_)
  b <- n * stats::var(rowMeans(split))
  sqrt((((n - 1) / n) * w + b / n) / w)
}

blr_convergence_contract_rank_rhat <- function(chains) {
  split <- blr_convergence_contract_split_chains(chains)
  z <- matrix(blr_convergence_contract_rank_normalize(split), nrow(split), ncol(split))
  blr_convergence_contract_rhat_basic(z)
}

blr_convergence_contract_folded_rhat <- function(chains) {
  folded <- abs(chains - stats::median(as.vector(chains)))
  blr_convergence_contract_rank_rhat(folded)
}

blr_convergence_contract_rhat <- function(chains) {
  max(blr_convergence_contract_rank_rhat(chains), blr_convergence_contract_folded_rhat(chains), na.rm = TRUE)
}

blr_convergence_contract_autocovariance <- function(x) {
  n <- length(x)
  variance <- stats::var(x)
  if (variance == 0) return(rep(0, n))
  size <- 2L * stats::nextn(n)
  centered <- c(x - mean(x), rep(0, size - n))
  transformed <- fft(centered)
  autocov <- Re(fft(Mod(transformed)^2, inverse = TRUE))[seq_len(n)]
  autocov / autocov[1L] * variance * (n - 1) / n
}

blr_convergence_contract_ess_split <- function(split) {
  n <- ncol(split); m <- nrow(split)
  if (m < 2L || n < 3L) return(NA_real_)
  x <- t(split)
  acov <- vapply(seq_len(m), function(i) blr_convergence_contract_autocovariance(x[, i]),
                  numeric(n))
  acov_means <- rowMeans(acov)
  w <- acov_means[1L] * n / (n - 1)
  var_plus <- w * (n - 1) / n + stats::var(colMeans(x))
  if (!is.finite(var_plus) || var_plus <= 0) return(NA_real_)
  rho <- numeric(n); at <- 0L; even <- 1
  rho[1L] <- even
  odd <- 1 - (w - acov_means[2L]) / var_plus
  rho[2L] <- odd
  while (at < n - 5L && is.finite(even + odd) && even + odd > 0) {
    at <- at + 2L
    even <- 1 - (w - acov_means[at + 1L]) / var_plus
    odd <- 1 - (w - acov_means[at + 2L]) / var_plus
    if (is.finite(even + odd) && even + odd >= 0) {
      rho[at + 1L] <- even
      rho[at + 2L] <- odd
    }
  }
  max_at <- at
  if (even > 0) rho[max_at + 1L] <- even
  at <- 0L
  while (at <= max_at - 4L) {
    at <- at + 2L
    if (rho[at + 1L] + rho[at + 2L] > rho[at - 1L] + rho[at]) {
      rho[at + 1L] <- (rho[at - 1L] + rho[at]) / 2
      rho[at + 2L] <- rho[at + 1L]
    }
  }
  retained <- if (max_at == 0L) rho[1L] else sum(rho[seq_len(max_at)])
  tau <- -1 + 2 * retained + rho[max_at + 1L]
  bound <- 1 / log10(m * n)
  if (tau < bound) tau <- bound
  if (!is.finite(tau) || tau <= 0) return(NA_real_)
  m * n / tau
}

blr_convergence_contract_ess_bulk <- function(chains) {
  split <- blr_convergence_contract_split_chains(chains)
  z <- matrix(blr_convergence_contract_rank_normalize(split), nrow(split), ncol(split))
  blr_convergence_contract_ess_split(z)
}

blr_convergence_contract_ess_tail <- function(chains) {
  pooled <- as.vector(chains)
  limits <- stats::quantile(pooled, c(0.05, 0.95), names = FALSE, type = 7)
  values <- vapply(limits, function(q) {
    blr_convergence_contract_ess_split(blr_convergence_contract_split_chains(1L * (chains <= q)))
  }, numeric(1))
  if (anyNA(values)) NA_real_ else min(values)
}

blr_convergence_contract_ess_mean <- function(chains) {
  blr_convergence_contract_ess_split(blr_convergence_contract_split_chains(chains))
}

blr_convergence_contract_mcse_mean <- function(chains) {
  posterior_sd <- stats::sd(as.vector(chains))
  ess <- blr_convergence_contract_ess_mean(chains)
  if (!is.finite(posterior_sd) || posterior_sd <= 0 || !is.finite(ess)) {
    return(c(ess_mean = NA_real_, posterior_sd = posterior_sd,
             mcse_mean = NA_real_, mcse_mean_over_sd = NA_real_))
  }
  mcse <- posterior_sd / sqrt(ess)
  c(ess_mean = ess, posterior_sd = posterior_sd, mcse_mean = mcse,
    mcse_mean_over_sd = mcse / posterior_sd)
}

blr_convergence_contract_status <- function(chains, updated = TRUE) {
  if (!updated) return("not_updated")
  if (any(!is.finite(chains))) return("nonfinite")
  if (nrow(chains) < 2L) return("unavailable_single_chain")
  if (ncol(chains) < 4L) return("insufficient_draws")
  if (length(unique(as.vector(chains))) == 1L) return("constant")
  chain_constant <- apply(chains, 1L, function(x) length(unique(x)) == 1L)
  if (any(chain_constant)) return("constant_chain_mismatch")
  if (nrow(chains) < 4L) "computed_fewer_than_four_chains" else "computed"
}

blr_convergence_contract_flags <- function(rhat, ess_bulk, ess_tail, mcse_relative, nchains,
                           rhat_threshold = 1.01,
                           ess_per_chain_threshold = 100,
                           mcse_threshold = 0.05) {
  c(rhat_flag = is.finite(rhat) && rhat > rhat_threshold,
    ess_bulk_flag = is.finite(ess_bulk) &&
      ess_bulk < ess_per_chain_threshold * nchains,
    ess_tail_flag = is.finite(ess_tail) &&
      ess_tail < ess_per_chain_threshold * nchains,
    mcse_flag = is.finite(mcse_relative) && mcse_relative > mcse_threshold)
}

blr_convergence_contract_overview <- function(summary) {
  computed <- summary$status %in% c("computed", "computed_fewer_than_four_chains")
  flagged <- rowSums(summary[c("rhat_flag", "ess_bulk_flag",
                               "ess_tail_flag", "mcse_flag")]) > 0L
  list(overall_status = if (!any(computed)) "unavailable" else if (any(flagged))
         "warning" else if (all(computed)) "ok" else "partial",
       n_computed = sum(computed), n_unavailable = sum(!computed),
       n_flagged = sum(flagged), max_rhat = max(summary$rhat, na.rm = TRUE),
       min_ess_bulk = min(summary$ess_bulk, na.rm = TRUE),
       min_ess_tail = min(summary$ess_tail, na.rm = TRUE),
       max_mcse_mean_over_sd = max(summary$mcse_mean_over_sd, na.rm = TRUE),
       fewer_than_four_chains = summary$nchains[1L] < 4L)
}

blr_convergence_contract_memory <- function(nchains, nit, nt, nmodels, selected_markers) {
  q <- nt * (nt + 1) / 2
  c(tier1_trace_bytes = 8 * nchains * nit * 5 * nt,
    covariance_trace_bytes = 8 * nchains * nit * 3 * q,
    probability_trace_bytes = 8 * nchains * nit * 2,
    full_pi_trace_bytes = 8 * nchains * nit * nmodels,
    selected_b_trace_bytes = 8 * nchains * nit * selected_markers * nt,
    selected_d_trace_bytes = 4 * nchains * nit * selected_markers * nt,
    per_quantity_workspace_bytes = 8 * nchains * nit)
}

blr_convergence_contract_scope_contract <- function() {
  list(tier1 = c("vbs", "vgs", "ves", "vle", "vld"),
       tier2_requires_new_traces = c("B_lower", "G_lower", "E_lower",
                                    "pi_null", "pi_active"),
       tier3_opt_in = c("selected_marker_b", "selected_marker_d"),
       all_markers_default = FALSE, full_pi_default = FALSE,
       keep_chains_required = FALSE, diagnostic_thinning = FALSE,
       trace_retention = "independent_bundle")
}

blr_convergence_contract_fixtures <- function() {
  half <- seq(-2, 2, length.out = 20)
  base <- t(vapply(0:3, function(i) {
    perm <- c(half[(seq_along(half) + 5L * i - 1L) %% 20L + 1L])
    c(perm, rev(perm))
  }, numeric(40)))
  list(well_mixed = base,
       shifted = base + c(0, 0, 1.5, 1.5),
       scales = base * c(1, 1, 2.5, 2.5),
       drift = base + rep(seq(0, 2, length.out = 40), each = 4),
       positive_ar = t(vapply(1:4, function(i) cumsum(sin((1:40 + i) / 7)),
                              numeric(40))),
       negative_ar = t(vapply(1:4, function(i) (-1)^(1:40) + i / 100,
                              numeric(40))),
       poor_tail = rbind(base[1:3, ], c(base[4, 1:35], rep(8, 5))),
       tied = matrix(rep(c(0, 0, 1, 1), 40), 4),
       binary = matrix(rep(c(0, 1), 80), 4),
       constant = matrix(2, 4, 40),
       one_constant = rbind(rep(0, 40), base[2:4, ]),
       nonfinite = { x <- base; x[1, 1] <- NA_real_; x },
       one_chain = base[1, , drop = FALSE],
       two_chains = base[1:2, , drop = FALSE],
       four_chains = base,
       odd = base[, 1:39], short = base[, 1:3])
}

# ---- consolidated from tests/testthat/helper-mtblr-bed-convergence-engine.R ----
blr_convergence_engine_args <- function(case, nchains = 2L, ncores = 1L,
                          chain_seeds = integer(),
                          keep_chains = FALSE) {
  blr_bed_chains_args(case, nchains = nchains, ncores = ncores,
                chain_seeds = chain_seeds, keep_chains = keep_chains)
}

blr_convergence_engine_native_call <- function(case, nchains = 2L, ncores = 1L,
                                 chain_seeds = integer(),
                                 keep_chains = FALSE) {
  do.call(sblr:::mtblr_bed_convergence_trace_internal,
          blr_convergence_engine_args(case, nchains, ncores, chain_seeds, keep_chains))
}

blr_convergence_engine_without_timing <- function(x) {
  if (!is.null(x$raw)) x$raw <- blr_bed_chains_without_timing(x$raw)
  x
}

blr_convergence_engine_bundle_from_chains <- function(chains, updated = rep(TRUE, 5L),
                                        group = c("vbs", "vgs", "ves",
                                                  "vle", "vld")) {
  stopifnot(is.matrix(chains), length(updated) == length(group))
  nit <- ncol(chains)
  nchains <- nrow(chains)
  values <- array(NA_real_, c(nit, nchains, length(group)))
  for (quantity in seq_along(group)) values[, , quantity] <- t(chains)
  list(
    schema = list(class = "mtblr_convergence_trace_bundle", version = 1L),
    scope = "core", nchains = as.integer(nchains),
    postburn_draws_per_chain = as.integer(nit),
    quantities = data.frame(
      quantity_index = seq_along(group), group = group,
      trait_index = rep.int(1L, length(group)),
      updated = updated, stringsAsFactors = FALSE),
    values = values)
}

blr_convergence_engine_fixture_metrics <- function(chains) {
  sblr:::.blr_convergence_scalar(t(chains))
}

# ---- consolidated from tests/testthat/helper-mtblr-bed-public-convergence.R ----
blr_public_convergence_public_args <- function(case, nchains = 1L, ncores = 1L,
                                 convergence = "auto", warn = FALSE,
                                 keep_traces = FALSE, keep_chains = FALSE,
                                 nit = 8L, nburn = 3L, nthin = 1L, ...) {
  blr_bed_public_chains_public_args(
    case, nchains = nchains, ncores = ncores, keep_chains = keep_chains,
    convergence = convergence, nit = nit, nburn = nburn, nthin = nthin,
    convergence_control = list(warn = warn, keep_traces = keep_traces), ...)
}

blr_public_convergence_without_additions <- function(fit) {
  fit$convergence <- NULL
  fit$convergence_traces <- NULL
  fit$input[c(
    "convergence", "convergence_requested", "convergence_scope",
    "convergence_trace_route", "convergence_warning_enabled",
    "convergence_warning_emitted", "convergence_keep_traces",
    "convergence_thresholds", "convergence_status",
    "convergence_computed", "convergence_memory_estimate")] <- NULL
  fit
}

blr_public_convergence_warning_fixture <- function(status = "computed", flagged = FALSE,
                                     nchains = 4L, nit = 10L) {
  result <- sblr:::.blr_convergence_unavailable(
    "T1", TRUE, TRUE, nchains, nit, sblr:::.blr_convergence_control())
  row <- result$summary[1L, ]
  row$status <- status
  if (status %in% c("computed", "computed_fewer_than_four_chains",
                    "computed_partial")) {
    row$rhat <- if (flagged) 1.2 else 1
    row$rhat_rank <- row$rhat_folded <- row$rhat
    row$rhat_available <- TRUE
    row$ess_bulk <- row$ess_tail <- row$ess_mean <- 1000
    row$ess_q05 <- row$ess_q95 <- 1000
    row$ess_bulk_per_chain <- row$ess_tail_per_chain <- 1000 / nchains
    row$ess_bulk_available <- row$ess_tail_available <- TRUE
    row$ess_mean_available <- TRUE
    row$posterior_sd <- 1
    row$mcse_mean <- 1 / sqrt(1000)
    row$mcse_mean_over_sd <- row$mcse_mean
    row$mcse_mean_available <- TRUE
    row$rhat_flag <- flagged
  }
  result$summary <- row
  result$availability <- list(
    rhat = sum(row$rhat_available),
    ess_bulk = sum(row$ess_bulk_available),
    ess_tail = sum(row$ess_tail_available),
    ess_mean = sum(row$ess_mean_available),
    mcse_mean = sum(row$mcse_mean_available))
  result$overview <- sblr:::.blr_convergence_overview(row)
  result$overall_status <- result$overview$overall_status
  result$computed <- any(unlist(result$availability) > 0L)
  sblr:::.blr_validate_convergence_result(result)
}

# ---- consolidated from tests/testthat/helper-mtblr-bed-extended-convergence-contract.R ----
blr_extended_contract_strict_lower <- function(nt) {
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

blr_extended_contract_covariance_descriptors <- function(traits, residual = "full",
                                            updateB = TRUE, updateE = TRUE) {
  pairs <- blr_extended_contract_strict_lower(length(traits))
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

blr_extended_contract_null_index <- function(models) {
  models <- as.matrix(models)
  hit <- which(rowSums(models) == 0L)
  if (length(hit) != 1L) stop("exactly one null pattern is required")
  hit
}

blr_extended_contract_resolve_patterns <- function(models, model_names, selection) {
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

blr_extended_contract_probability_plan <- function(models, model_names, probability = "mass",
                                      selection = NULL, updatePi = TRUE) {
  probability <- match.arg(probability, c("none", "mass", "selected", "all"))
  null <- blr_extended_contract_null_index(models)
  selected <- blr_extended_contract_resolve_patterns(models, model_names, selection)
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

blr_extended_contract_resolve_markers <- function(marker_ids, marker_ids_request = NULL,
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

blr_extended_contract_extended_memory <- function(C, N, T, updateB = TRUE, updateE = TRUE,
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

blr_extended_contract_quantity_counts <- function(T, probability_rows = 0L, K = 0L,
                                     marker_quantities = c("b", "d")) {
  qoff <- T * (T - 1L) / 2L
  as.integer(c(tier1 = 3L * T, covariance = 3L * qoff,
    probability = probability_rows,
    marker_b = as.integer("b" %in% marker_quantities) * K * T,
    marker_d = as.integer("d" %in% marker_quantities) * K * T)) |>
    stats::setNames(c("tier1", "covariance", "probability", "marker_b", "marker_d"))
}

blr_extended_contract_resolved_limit_gb <- function(memory_warning_gb = 8,
                                       trace_limit_gb = NULL) {
  if (!is.null(trace_limit_gb)) return(trace_limit_gb)
  if (!is.finite(memory_warning_gb)) return(2)
  min(2, max(.25, .25 * memory_warning_gb))
}

blr_extended_contract_large_request <- function(capture_bytes, memory_warning_gb = 8,
                                   trace_limit_gb = NULL,
                                   allow_large_traces = FALSE) {
  limit <- blr_extended_contract_resolved_limit_gb(memory_warning_gb, trace_limit_gb)
  exceeds <- capture_bytes > limit * 1024^3
  list(limit_gb = limit, exceeds = exceeds,
       decision = if (!exceeds || allow_large_traces) "allow" else
         "require_explicit_override")
}

blr_extended_contract_validate_future_control <- function(convergence = "extended",
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

blr_extended_contract_warning_groups <- function(summary) {
  flagged <- with(summary, rhat_flag | ess_bulk_flag | ess_tail_flag | mcse_flag |
                    status == "constant_chain_mismatch")
  table(factor(summary$group[flagged], levels = unique(summary$group)))
}

# ---- consolidated from tests/testthat/helper-mtblr-bed-multichain-contract.R ----
blr_bed_multichain_uint32 <- function(x) {
  modulus <- 2^32
  ((as.numeric(x) %% modulus) + modulus) %% modulus
}

blr_bed_multichain_tasks <- function(nchains, nt = 1L) {
  stopifnot(length(nchains) == 1L, nchains >= 1L, nt >= 1L)
  data.frame(chain = seq.int(0L, as.integer(nchains) - 1L),
             result_slot = seq.int(0L, as.integer(nchains) - 1L))
}

blr_bed_multichain_seeds <- function(seed, nchains, chain_seeds = NULL) {
  if (!is.null(chain_seeds)) {
    stopifnot(length(chain_seeds) == nchains, all(is.finite(chain_seeds)))
    return(vapply(chain_seeds, blr_bed_multichain_uint32, numeric(1)))
  }
  vapply(seq.int(0, nchains - 1), function(chain) {
    blr_bed_multichain_uint32(as.numeric(seed) + 9176 * chain)
  }, numeric(1))
}

blr_bed_multichain_chain <- function(index, retained, bm, dm, final_state,
                           trace_offset = 0, failed = FALSE, error = "") {
  bm <- as.matrix(bm); dm <- as.matrix(dm)
  nt <- ncol(bm)
  list(
    chain = as.integer(index), failed = failed, error = error,
    retained = as.integer(retained),
    bm_acc = bm * retained, dm_acc = dm * retained,
    covb_acc = diag(nt) * retained * index,
    covg_acc = diag(nt) * retained * (index + 1),
    cove_acc = diag(nt) * retained * (index + 2),
    pi_acc = c(.25, .75) * retained,
    vbs = matrix(seq_len(3 * nt) + trace_offset, nt),
    vgs = matrix(seq_len(3 * nt) + trace_offset + 10, nt),
    ves = matrix(seq_len(3 * nt) + trace_offset + 20, nt),
    b = bm + index, state = as.matrix(final_state), r = bm - index,
    B = diag(nt) * index, G = diag(nt) * (index + 1),
    E = diag(nt) * (index + 2), pi_final = c(.1, .9),
    seed = blr_bed_multichain_uint32(index), seconds = index / 10
  )
}

blr_bed_multichain_stability <- function(values) {
  arr <- simplify2array(values)
  if (length(values) == 1L) {
    return(list(sd = array(0, dim(values[[1]])),
                min = values[[1]], max = values[[1]]))
  }
  margin <- seq_len(length(dim(arr)) - 1L)
  list(sd = apply(arr, margin, stats::sd),
       min = apply(arr, margin, min),
       max = apply(arr, margin, max))
}

blr_bed_multichain_aggregate <- function(chains, keep_chains = FALSE) {
  failures <- chains[vapply(chains, `[[`, logical(1), "failed")]
  if (length(failures)) {
    failures <- failures[order(vapply(failures, `[[`, integer(1), "chain"))]
    stop(paste(vapply(failures, function(x) {
      sprintf("chain %d: %s", x$chain, x$error)
    }, character(1)), collapse = "; "), call. = FALSE)
  }
  chains <- chains[order(vapply(chains, `[[`, integer(1), "chain"))]
  retained <- sum(vapply(chains, `[[`, integer(1), "retained"))
  sum_field <- function(field) Reduce(`+`, lapply(chains, `[[`, field))
  mean_trace <- function(field) sum_field(field) / length(chains)
  bm_chain <- lapply(chains, function(x) x$bm_acc / x$retained)
  dm_chain <- lapply(chains, function(x) x$dm_acc / x$retained)
  compact <- function(x) x[c("chain", "seed", "bm_acc", "dm_acc", "b",
                              "state", "vbs", "vgs", "ves", "B", "G", "E",
                              "pi_final", "seconds")]
  list(
    bm = sum_field("bm_acc") / retained,
    dm = sum_field("dm_acc") / retained,
    covb = sum_field("covb_acc") / retained,
    covg = sum_field("covg_acc") / retained,
    cove = sum_field("cove_acc") / retained,
    pi_mean = sum_field("pi_acc") / retained,
    vbs = mean_trace("vbs"), vgs = mean_trace("vgs"),
    ves = mean_trace("ves"),
    b = chains[[1]]$b, state = chains[[1]]$state,
    r = chains[[1]]$r, B = chains[[1]]$B, G = chains[[1]]$G,
    E = chains[[1]]$E, pi_final = chains[[1]]$pi_final,
    bm_stability = blr_bed_multichain_stability(bm_chain),
    dm_stability = blr_bed_multichain_stability(dm_chain),
    chains = if (keep_chains) setNames(lapply(chains, compact),
                                      paste0("chain", seq_along(chains))) else NULL,
    retained = retained
  )
}

blr_bed_multichain_memory <- function(n, m, nt, nmodels, trace_length,
                            nchains, ncores, keep_chains = FALSE) {
  packed <- m * ceiling(n / 4)
  shared <- packed + 8 * (n * nt + 5 * m + m * nt) + 4 * m
  private <- 8 * (n * nt + 2 * m * nt + nt^2 * 4 + n) + 4 * m * nt
  result <- 8 * (4 * m * nt + 3 * nt * trace_length + 4 * nt^2 + nmodels)
  retained <- 8 * (3 * m * nt + 3 * nt * trace_length + 3 * nt^2 + nmodels)
  workers <- min(as.integer(nchains), as.integer(ncores))
  retained_total <- if (keep_chains) nchains * retained else 0
  total <- shared + workers * private + nchains * result + retained_total
  list(shared_bytes = shared, private_state_bytes_per_chain = private,
       result_bytes_per_chain = result,
       retained_chain_bytes_per_chain = retained,
       worker_count = workers, nchains = nchains,
       estimated_concurrent_bytes = shared + workers * private,
       estimated_retained_output_bytes = retained_total,
       estimated_total_bytes = total, estimated_total_gib = total / 2^30,
       estimate_kind = "analytical upper-bound estimate; not measured RSS or peak RSS")
}

blr_bed_multichain_future_controls <- function() {
  list(nchains = 1L, ncores = 1L, chain_seeds = NULL, keep_chains = FALSE,
       task_topology = "one complete joint MT chain",
       openmp_unavailable = "warn once and run serial",
       raw_schema = "mtblr_raw version 1",
       final_state_policy = "primary_chain",
       primary_chain = 1L,
       posterior_summary_policy = "pooled_retained_samples",
       trace_policy = "iterationwise_chain_mean")
}

blr_bed_multichain_future_internal_signature <- function() {
  "mtblr_bed_chains_internal(current arguments, int nchains, int ncores, std::vector<int> chain_seeds, bool keep_chains)"
}

