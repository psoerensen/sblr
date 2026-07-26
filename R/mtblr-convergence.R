.mtblr_convergence_control <- function(control = NULL) {
  defaults <- list(
    rhat_threshold = 1.01,
    ess_per_chain_threshold = 100,
    mcse_mean_over_sd_threshold = 0.05
  )
  if (is.null(control)) return(defaults)
  if (!is.list(control) || is.null(names(control)) ||
      any(!names(control) %in% names(defaults)) ||
      anyDuplicated(names(control))) {
    stop("convergence control must be a named list of supported thresholds.",
         call. = FALSE)
  }
  out <- utils::modifyList(defaults, control)
  for (name in names(out)) {
    value <- out[[name]]
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
        value <= 0) {
      stop(name, " must be one finite positive number.", call. = FALSE)
    }
  }
  out
}

.mtblr_bed_convergence_controls <- function(convergence,
                                             convergence_control,
                                             nchains) {
  defaults <- list(
    warn = TRUE,
    rhat_threshold = 1.01,
    ess_per_chain_threshold = 100,
    mcse_mean_over_sd_threshold = 0.05,
    keep_traces = FALSE
  )
  if (is.null(convergence_control)) {
    out <- defaults
  } else {
    if (!is.list(convergence_control) || is.data.frame(convergence_control) ||
        is.null(names(convergence_control)) ||
        any(!nzchar(names(convergence_control))) ||
        anyDuplicated(names(convergence_control)) ||
        any(!names(convergence_control) %in% names(defaults))) {
      stop("convergence_control must be a uniquely named list containing only supported controls.",
           call. = FALSE)
    }
    out <- utils::modifyList(defaults, convergence_control)
  }
  for (name in c("warn", "keep_traces")) {
    if (!is.logical(out[[name]]) || length(out[[name]]) != 1L ||
        is.na(out[[name]])) {
      stop("convergence_control$", name,
           " must be TRUE or FALSE.", call. = FALSE)
    }
  }
  for (name in c("rhat_threshold", "ess_per_chain_threshold",
                 "mcse_mean_over_sd_threshold")) {
    value <- out[[name]]
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
        value <= 0) {
      stop("convergence_control$", name,
           " must be one finite positive number.", call. = FALSE)
    }
  }
  if (identical(convergence, "none") && isTRUE(out$keep_traces)) {
    stop("convergence = 'none' cannot retain convergence traces.",
         call. = FALSE)
  }
  out$diagnostic_requested <- !identical(convergence, "none")
  out$trace_route_required <- isTRUE(out$keep_traces) ||
    (out$diagnostic_requested && nchains >= 2L)
  out$thresholds <- out[c(
    "rhat_threshold", "ess_per_chain_threshold",
    "mcse_mean_over_sd_threshold")]
  out
}

.mtblr_convergence_summary_columns <- function() c(
  "quantity", "group", "trait", "trait2", "marker_id", "model_name",
  "updated", "status", "rhat", "rhat_rank", "rhat_folded",
  "ess_bulk", "ess_bulk_per_chain", "ess_q05", "ess_q95",
  "ess_tail", "ess_tail_per_chain", "ess_mean", "posterior_sd",
  "mcse_mean", "mcse_mean_over_sd", "nchains", "draws_per_chain",
  "split_draws_per_chain", "rhat_available", "ess_bulk_available",
  "ess_tail_available", "ess_mean_available", "mcse_mean_available",
  "ess_stability_bound_applied", "rhat_flag", "ess_bulk_flag",
  "ess_tail_flag", "mcse_flag")

.mtblr_convergence_summary_row <- function(group, trait, updated, status,
                                            nchains, nit) {
  data.frame(
    quantity = paste0(group, "[", trait, "]"), group = group,
    trait = trait, trait2 = NA_character_, marker_id = NA_character_,
    model_name = NA_character_, updated = updated, status = status,
    rhat = NA_real_, rhat_rank = NA_real_, rhat_folded = NA_real_,
    ess_bulk = NA_real_, ess_bulk_per_chain = NA_real_,
    ess_q05 = NA_real_, ess_q95 = NA_real_, ess_tail = NA_real_,
    ess_tail_per_chain = NA_real_, ess_mean = NA_real_,
    posterior_sd = NA_real_, mcse_mean = NA_real_,
    mcse_mean_over_sd = NA_real_, nchains = as.integer(nchains),
    draws_per_chain = as.integer(nit),
    split_draws_per_chain = as.integer(nit %/% 2L),
    rhat_available = FALSE, ess_bulk_available = FALSE,
    ess_tail_available = FALSE, ess_mean_available = FALSE,
    mcse_mean_available = FALSE, ess_stability_bound_applied = FALSE,
    rhat_flag = FALSE, ess_bulk_flag = FALSE, ess_tail_flag = FALSE,
    mcse_flag = FALSE, stringsAsFactors = FALSE, check.names = FALSE)
}

.mtblr_convergence_empty_overview <- function(nchains,
                                               status = "not_requested") {
  list(
    overall_status = status, n_quantities = 0L, n_computed = 0L,
    n_partial = 0L, n_unavailable = 0L, n_not_updated = 0L,
    n_flagged = 0L, max_rhat = NA_real_,
    max_rhat_quantity = NA_character_, min_ess_bulk = NA_real_,
    min_ess_bulk_per_chain = NA_real_,
    min_ess_bulk_quantity = NA_character_, min_ess_tail = NA_real_,
    min_ess_tail_per_chain = NA_real_,
    min_ess_tail_quantity = NA_character_,
    max_mcse_mean_over_sd = NA_real_, max_mcse_quantity = NA_character_,
    n_constant_chain_mismatch = 0L,
    n_ess_stability_bound_applied = 0L,
    fewer_than_four_chains = nchains < 4L)
}

.mtblr_convergence_algorithm <- function() list(
  authority = "Vehtari et al. 2021",
  implementation_target = "posterior 1.6.1",
  reference_package = "posterior", reference_version = "1.6.1",
  contract_amendments = c(
    "Blom rank denominator is S + 1/4",
    "ESS stability bound is nominal times log10(nominal)",
    "R-hat needs four draws; ESS and MCSE need six"))

.mtblr_convergence_not_requested <- function(trait_names, updateB, updateE,
                                               nchains, nit, control) {
  template <- .mtblr_convergence_summary_row(
    "G_diag", trait_names[1L], TRUE, "unavailable_single_chain",
    nchains, nit)[0L, ]
  result <- list(
    version = 1L, requested = FALSE, computed = FALSE, scope = "none",
    overall_status = "not_requested", nchains = as.integer(nchains),
    postburn_draws_per_chain = as.integer(nit),
    split_draws_per_chain = as.integer(nit %/% 2L),
    total_postburn_draws = as.integer(nchains * nit), thresholds = control,
    availability = list(rhat = 0L, ess_bulk = 0L, ess_tail = 0L,
                        ess_mean = 0L, mcse_mean = 0L),
    summary = template,
    overview = .mtblr_convergence_empty_overview(nchains),
    warning_messages = character(),
    trace_retention = list(retained = FALSE, burnin_included = FALSE,
                           additional_thinning = FALSE),
    algorithm = .mtblr_convergence_algorithm())
  .mtblr_validate_convergence_result(result)
}

.mtblr_convergence_unavailable <- function(trait_names, updateB, updateE,
                                            nchains, nit, control,
                                            reason = "unavailable_single_chain",
                                            keep_traces = FALSE) {
  groups <- rep(c("B_diag", "G_diag", "E_diag"), each = length(trait_names))
  traits <- rep(trait_names, 3L)
  updated <- c(rep(isTRUE(updateB), length(trait_names)),
               rep(TRUE, length(trait_names)),
               rep(isTRUE(updateE), length(trait_names)))
  rows <- lapply(seq_along(groups), function(i) {
    .mtblr_convergence_summary_row(
      groups[i], traits[i], updated[i],
      if (updated[i]) reason else "not_updated", nchains, nit)
  })
  summary <- do.call(rbind, rows)
  rownames(summary) <- NULL
  overview <- .mtblr_convergence_overview(summary)
  result <- list(
    version = 1L, requested = TRUE, computed = FALSE, scope = "core",
    overall_status = overview$overall_status, nchains = as.integer(nchains),
    postburn_draws_per_chain = as.integer(nit),
    split_draws_per_chain = as.integer(nit %/% 2L),
    total_postburn_draws = as.integer(nchains * nit), thresholds = control,
    availability = list(rhat = 0L, ess_bulk = 0L, ess_tail = 0L,
                        ess_mean = 0L, mcse_mean = 0L),
    summary = summary, overview = overview, warning_messages = character(),
    trace_retention = list(retained = isTRUE(keep_traces),
                           burnin_included = FALSE,
                           additional_thinning = FALSE),
    algorithm = .mtblr_convergence_algorithm())
  .mtblr_validate_convergence_result(result)
}

.mtblr_convergence_split_chains <- function(x) {
  if (!is.matrix(x) || !is.numeric(x)) {
    stop("convergence draws must be a numeric iteration-by-chain matrix.",
         call. = FALSE)
  }
  half <- nrow(x) %/% 2L
  if (!half) return(matrix(numeric(), 0L, 2L * ncol(x)))
  cbind(x[seq_len(half), , drop = FALSE],
        x[seq.int(nrow(x) - half + 1L, nrow(x)), , drop = FALSE])
}

.mtblr_convergence_rank_normalize <- function(x) {
  ranks <- rank(as.vector(x), ties.method = "average")
  # posterior 1.6.1 backtransform_ranks(): S - 2c + 1 = S + 1/4.
  values <- stats::qnorm((ranks - 3 / 8) / (length(ranks) + 1 / 4))
  array(values, dim = dim(x), dimnames = dimnames(x))
}

.mtblr_convergence_fold <- function(x) {
  abs(x - stats::median(as.vector(x)))
}

.mtblr_convergence_is_constant <- function(x) {
  diff(range(x)) < .Machine$double.eps
}

.mtblr_convergence_rhat_basic <- function(x) {
  if (!is.matrix(x) || nrow(x) < 2L || ncol(x) < 2L ||
      any(!is.finite(x)) || .mtblr_convergence_is_constant(x)) {
    return(NA_real_)
  }
  chain_means <- colMeans(x)
  chain_vars <- apply(x, 2L, stats::var)
  within <- mean(chain_vars)
  if (!is.finite(within) || within <= 0) return(NA_real_)
  between <- nrow(x) * stats::var(chain_means)
  sqrt((between / within + nrow(x) - 1) / nrow(x))
}

.mtblr_convergence_autocovariance <- function(x) {
  n <- length(x)
  variance <- stats::var(x)
  if (!is.finite(variance) || variance == 0) return(rep(0, n))
  padded_length <- 2L * stats::nextn(n)
  centered <- c(x - mean(x), rep.int(0, padded_length - n))
  transformed <- stats::fft(centered)
  autocov <- Re(stats::fft(Mod(transformed)^2, inverse = TRUE)[seq_len(n)])
  autocov / autocov[1L] * variance * (n - 1) / n
}

.mtblr_convergence_ess <- function(x) {
  unavailable <- list(value = NA_real_, stability_bound_applied = FALSE)
  if (!is.matrix(x) || nrow(x) < 3L || any(!is.finite(x)) ||
      .mtblr_convergence_is_constant(x)) {
    return(unavailable)
  }
  niterations <- nrow(x)
  nchains <- ncol(x)
  autocov <- vapply(seq_len(nchains), function(chain) {
    .mtblr_convergence_autocovariance(x[, chain])
  }, numeric(niterations))
  if (is.null(dim(autocov))) autocov <- matrix(autocov, ncol = 1L)
  autocov_means <- rowMeans(autocov)
  mean_variance <- autocov_means[1L] * niterations / (niterations - 1)
  variance_plus <- mean_variance * (niterations - 1) / niterations
  if (nchains > 1L) variance_plus <- variance_plus + stats::var(colMeans(x))
  if (!is.finite(variance_plus) || variance_plus <= 0) return(unavailable)

  rho <- numeric(niterations)
  lag <- 0L
  rho_even <- 1
  rho[1L] <- rho_even
  rho_odd <- 1 - (mean_variance - autocov_means[2L]) / variance_plus
  rho[2L] <- rho_odd
  while (lag < nrow(autocov) - 5L && is.finite(rho_even + rho_odd) &&
         rho_even + rho_odd > 0) {
    lag <- lag + 2L
    rho_even <- 1 -
      (mean_variance - autocov_means[lag + 1L]) / variance_plus
    rho_odd <- 1 -
      (mean_variance - autocov_means[lag + 2L]) / variance_plus
    if (is.finite(rho_even + rho_odd) && rho_even + rho_odd >= 0) {
      rho[lag + 1L] <- rho_even
      rho[lag + 2L] <- rho_odd
    }
  }
  max_lag <- lag
  if (is.finite(rho_even) && rho_even > 0) rho[max_lag + 1L] <- rho_even
  lag <- 0L
  while (lag <= max_lag - 4L) {
    lag <- lag + 2L
    if (rho[lag + 1L] + rho[lag + 2L] >
        rho[lag - 1L] + rho[lag]) {
      rho[lag + 1L] <- (rho[lag - 1L] + rho[lag]) / 2
      rho[lag + 2L] <- rho[lag + 1L]
    }
  }
  # This explicit form reproduces posterior's 1:max_lag behavior at zero.
  retained <- if (max_lag == 0L) rho[1L] else sum(rho[seq_len(max_lag)])
  tau <- -1 + 2 * retained + rho[max_lag + 1L]
  nominal <- nchains * niterations
  bound <- 1 / log10(nominal)
  bounded <- is.finite(tau) && tau < bound
  if (bounded) tau <- bound
  if (!is.finite(tau) || tau <= 0) return(unavailable)
  list(value = nominal / tau, stability_bound_applied = bounded)
}

.mtblr_convergence_scalar <- function(draws, updated = TRUE,
                                      control = NULL) {
  control <- .mtblr_convergence_control(control)
  if (!is.matrix(draws) || !is.numeric(draws) || !nrow(draws) ||
      !ncol(draws)) {
    stop("scalar convergence draws must be a nonempty numeric matrix.",
         call. = FALSE)
  }
  nchains <- ncol(draws)
  ndraws <- nrow(draws)
  finite <- all(is.finite(draws))
  constant <- finite && .mtblr_convergence_is_constant(draws)
  chain_constant <- finite && any(vapply(seq_len(nchains), function(chain) {
    .mtblr_convergence_is_constant(draws[, chain])
  }, logical(1)))

  base_status <- if (!updated) "not_updated" else if (!finite) {
    "nonfinite"
  } else if (constant) {
    "constant"
  } else if (nchains < 2L) {
    "unavailable_single_chain"
  } else if (ndraws < 4L) {
    "insufficient_draws"
  } else if (chain_constant) {
    "constant_chain_mismatch"
  } else {
    NA_character_
  }
  rhat_available <- updated && finite && !constant && nchains >= 2L &&
    ndraws >= 4L
  ess_available <- updated && finite && !constant && nchains >= 2L &&
    ndraws >= 6L

  rhat_rank <- rhat_folded <- rhat <- NA_real_
  if (rhat_available) {
    split <- .mtblr_convergence_split_chains(draws)
    rhat_rank <- .mtblr_convergence_rhat_basic(
      .mtblr_convergence_rank_normalize(split))
    folded_split <- .mtblr_convergence_split_chains(
      .mtblr_convergence_fold(draws))
    rhat_folded <- .mtblr_convergence_rhat_basic(
      .mtblr_convergence_rank_normalize(folded_split))
    rhat <- max(rhat_rank, rhat_folded)
    rhat_available <- is.finite(rhat)
  }

  ess_bulk <- ess_q05 <- ess_q95 <- ess_tail <- ess_mean <- NA_real_
  bound_applied <- FALSE
  if (ess_available) {
    split <- .mtblr_convergence_split_chains(draws)
    bulk <- .mtblr_convergence_ess(
      .mtblr_convergence_rank_normalize(split))
    thresholds <- stats::quantile(
      as.vector(draws), c(0.05, 0.95), names = FALSE, type = 7)
    q05 <- .mtblr_convergence_ess(.mtblr_convergence_split_chains(
      1 * (draws <= thresholds[1L])))
    q95 <- .mtblr_convergence_ess(.mtblr_convergence_split_chains(
      1 * (draws <= thresholds[2L])))
    mean_ess <- .mtblr_convergence_ess(split)
    ess_bulk <- bulk$value
    ess_q05 <- q05$value
    ess_q95 <- q95$value
    ess_tail <- if (all(is.finite(c(ess_q05, ess_q95)))) {
      min(ess_q05, ess_q95)
    } else NA_real_
    ess_mean <- mean_ess$value
    bound_applied <- any(c(
      bulk$stability_bound_applied, q05$stability_bound_applied,
      q95$stability_bound_applied, mean_ess$stability_bound_applied))
  }
  bulk_available <- is.finite(ess_bulk)
  tail_available <- is.finite(ess_tail)
  mean_available <- is.finite(ess_mean)
  posterior_sd <- mcse_mean <- mcse_relative <- NA_real_
  if (mean_available) {
    posterior_sd <- stats::sd(as.vector(draws))
    if (is.finite(posterior_sd) && posterior_sd > 0) {
      mcse_mean <- posterior_sd / sqrt(ess_mean)
      mcse_relative <- mcse_mean / posterior_sd
    }
  }
  mcse_available <- is.finite(mcse_mean)

  if (is.na(base_status)) {
    complete <- rhat_available && bulk_available && tail_available &&
      mean_available && mcse_available
    base_status <- if (!complete) "computed_partial" else if (nchains < 4L) {
      "computed_fewer_than_four_chains"
    } else "computed"
  }
  mismatch <- identical(base_status, "constant_chain_mismatch")
  rhat_flag <- mismatch ||
    (rhat_available && rhat > control$rhat_threshold)
  ess_bulk_flag <- bulk_available &&
    ess_bulk / nchains < control$ess_per_chain_threshold
  ess_tail_flag <- tail_available &&
    ess_tail / nchains < control$ess_per_chain_threshold
  mcse_flag <- mcse_available &&
    mcse_relative > control$mcse_mean_over_sd_threshold

  list(
    status = base_status,
    rhat = rhat, rhat_rank = rhat_rank, rhat_folded = rhat_folded,
    ess_bulk = ess_bulk,
    ess_bulk_per_chain = if (bulk_available) ess_bulk / nchains else NA_real_,
    ess_q05 = ess_q05, ess_q95 = ess_q95,
    ess_tail = ess_tail,
    ess_tail_per_chain = if (tail_available) ess_tail / nchains else NA_real_,
    ess_mean = ess_mean, posterior_sd = posterior_sd,
    mcse_mean = mcse_mean, mcse_mean_over_sd = mcse_relative,
    nchains = nchains, draws_per_chain = ndraws,
    split_draws_per_chain = ndraws %/% 2L,
    rhat_available = rhat_available,
    ess_bulk_available = bulk_available,
    ess_tail_available = tail_available,
    ess_mean_available = mean_available,
    mcse_mean_available = mcse_available,
    ess_stability_bound_applied = bound_applied,
    rhat_flag = rhat_flag, ess_bulk_flag = ess_bulk_flag,
    ess_tail_flag = ess_tail_flag, mcse_flag = mcse_flag
  )
}

.mtblr_validate_convergence_trace_bundle <- function(bundle,
                                                      nt = NULL,
                                                      updateB = NULL,
                                                      updateE = NULL) {
  if (!is.list(bundle) ||
      !identical(bundle$schema$class, "mtblr_convergence_trace_bundle") ||
      !identical(as.integer(bundle$schema$version), 1L) ||
      !identical(bundle$scope, "core")) {
    stop("Invalid MT BLR convergence trace-bundle schema.", call. = FALSE)
  }
  nchains <- as.integer(bundle$nchains)
  nit <- as.integer(bundle$postburn_draws_per_chain)
  quantities <- bundle$quantities
  if (length(nchains) != 1L || is.na(nchains) || nchains < 1L ||
      length(nit) != 1L || is.na(nit) || nit < 1L ||
      !is.data.frame(quantities) ||
      !identical(names(quantities),
                 c("quantity_index", "group", "trait_index", "updated")) ||
      !is.array(bundle$values) ||
      !identical(dim(bundle$values),
                 c(nit, nchains, nrow(quantities)))) {
    stop("Invalid MT BLR convergence trace-bundle dimensions.", call. = FALSE)
  }
  inferred_nt <- nrow(quantities) %/% 3L
  if (!is.null(nt) && !identical(as.integer(nt), as.integer(inferred_nt))) {
    stop("Convergence trait count does not match the trace bundle.",
         call. = FALSE)
  }
  groups <- rep(c("B_diag", "G_diag", "E_diag"), each = inferred_nt)
  traits <- rep(seq_len(inferred_nt), 3L)
  expected_updated <- c(
    rep(if (is.null(updateB)) quantities$updated[1L] else isTRUE(updateB),
        inferred_nt),
    rep(TRUE, inferred_nt),
    rep(if (is.null(updateE)) utils::tail(quantities$updated, 1L) else isTRUE(updateE),
        inferred_nt)
  )
  if (inferred_nt < 1L || nrow(quantities) != 3L * inferred_nt ||
      !identical(as.integer(quantities$quantity_index),
                 seq_len(nrow(quantities))) ||
      !identical(as.character(quantities$group), groups) ||
      !identical(as.integer(quantities$trait_index), traits) ||
      !identical(as.logical(quantities$updated), expected_updated) ||
      any(!is.finite(bundle$values[, , quantities$updated, drop = FALSE]))) {
    stop("Invalid MT BLR convergence trace-bundle quantities.",
         call. = FALSE)
  }
  bundle
}

.mtblr_convergence_extreme <- function(values, quantities, fn) {
  eligible <- which(is.finite(values))
  if (!length(eligible)) return(list(value = NA_real_, quantity = NA_character_))
  local <- fn(values[eligible])
  index <- eligible[which(values[eligible] == local)[1L]]
  list(value = unname(local), quantity = quantities[index])
}

.mtblr_convergence_overview <- function(summary) {
  computed <- summary$status %in% c(
    "computed", "computed_fewer_than_four_chains",
    "computed_partial", "constant_chain_mismatch")
  partial <- summary$status == "computed_partial"
  not_updated <- summary$status == "not_updated"
  flagged <- rowSums(summary[c(
    "rhat_flag", "ess_bulk_flag", "ess_tail_flag", "mcse_flag")]) > 0L
  rhat <- .mtblr_convergence_extreme(summary$rhat, summary$quantity, max)
  bulk <- .mtblr_convergence_extreme(summary$ess_bulk, summary$quantity, min)
  bulk_pc <- .mtblr_convergence_extreme(
    summary$ess_bulk_per_chain, summary$quantity, min)
  tail_ess <- .mtblr_convergence_extreme(summary$ess_tail, summary$quantity, min)
  tail_pc <- .mtblr_convergence_extreme(
    summary$ess_tail_per_chain, summary$quantity, min)
  mcse <- .mtblr_convergence_extreme(
    summary$mcse_mean_over_sd, summary$quantity, max)
  status <- if (!any(computed)) "unavailable" else if (any(flagged)) {
    "warning"
  } else if (any(partial) || any(!computed & !not_updated)) {
    "partial"
  } else "ok"
  list(
    overall_status = status,
    n_quantities = nrow(summary),
    n_computed = sum(computed),
    n_partial = sum(partial),
    n_unavailable = sum(!computed & !not_updated),
    n_not_updated = sum(not_updated),
    n_flagged = sum(flagged),
    max_rhat = rhat$value, max_rhat_quantity = rhat$quantity,
    min_ess_bulk = bulk$value,
    min_ess_bulk_per_chain = bulk_pc$value,
    min_ess_bulk_quantity = bulk$quantity,
    min_ess_tail = tail_ess$value,
    min_ess_tail_per_chain = tail_pc$value,
    min_ess_tail_quantity = tail_ess$quantity,
    max_mcse_mean_over_sd = mcse$value,
    max_mcse_quantity = mcse$quantity,
    n_constant_chain_mismatch =
      sum(summary$status == "constant_chain_mismatch"),
    n_ess_stability_bound_applied =
      sum(summary$ess_stability_bound_applied),
    fewer_than_four_chains = summary$nchains[1L] < 4L
  )
}

.mtblr_convergence_tier1 <- function(bundle, trait_names,
                                      control = NULL,
                                      keep_traces = FALSE) {
  bundle <- .mtblr_validate_convergence_trace_bundle(bundle)
  control <- .mtblr_convergence_control(control)
  nt <- nrow(bundle$quantities) %/% 3L
  if (!is.character(trait_names) || length(trait_names) != nt ||
      anyNA(trait_names) || any(!nzchar(trait_names)) ||
      anyDuplicated(trait_names)) {
    stop("trait_names must uniquely name every convergence trait.",
         call. = FALSE)
  }
  rows <- lapply(seq_len(nrow(bundle$quantities)), function(quantity) {
    descriptor <- bundle$quantities[quantity, ]
    draws <- matrix(
      bundle$values[, , quantity, drop = FALSE],
      nrow = bundle$postburn_draws_per_chain,
      ncol = bundle$nchains)
    metric <- .mtblr_convergence_scalar(
      draws,
      updated = descriptor$updated, control = control)
    identity <- list(
      quantity = paste0(descriptor$group, "[",
                        trait_names[descriptor$trait_index], "]"),
      group = as.character(descriptor$group),
      trait = trait_names[descriptor$trait_index],
      trait2 = NA_character_, marker_id = NA_character_,
      model_name = NA_character_, updated = as.logical(descriptor$updated)
    )
    as.data.frame(c(identity, metric), stringsAsFactors = FALSE,
                  optional = TRUE)
  })
  summary <- do.call(rbind, rows)
  rownames(summary) <- NULL
  overview <- .mtblr_convergence_overview(summary)
  availability <- list(
    rhat = sum(summary$rhat_available),
    ess_bulk = sum(summary$ess_bulk_available),
    ess_tail = sum(summary$ess_tail_available),
    ess_mean = sum(summary$ess_mean_available),
    mcse_mean = sum(summary$mcse_mean_available)
  )
  result <- list(
    version = 1L, requested = TRUE,
    computed = any(unlist(availability, use.names = FALSE) > 0L),
    scope = "core", overall_status = overview$overall_status,
    nchains = as.integer(bundle$nchains),
    postburn_draws_per_chain =
      as.integer(bundle$postburn_draws_per_chain),
    split_draws_per_chain =
      as.integer(bundle$postburn_draws_per_chain %/% 2L),
    total_postburn_draws =
      as.integer(bundle$nchains * bundle$postburn_draws_per_chain),
    thresholds = control, availability = availability,
    summary = summary, overview = overview,
    warning_messages = if (overview$n_flagged) {
      paste0(overview$n_flagged,
             " Tier 1 convergence quantities exceed advisory thresholds.")
    } else character(),
    trace_retention = list(
      retained = isTRUE(keep_traces), burnin_included = FALSE,
      additional_thinning = FALSE),
    algorithm = .mtblr_convergence_algorithm()
  )
  .mtblr_validate_convergence_result(result)
}

.mtblr_validate_convergence_result <- function(result) {
  required <- c(
    "version", "requested", "computed", "scope", "overall_status",
    "nchains", "postburn_draws_per_chain", "split_draws_per_chain",
    "total_postburn_draws", "thresholds", "availability", "summary",
    "overview", "warning_messages", "trace_retention", "algorithm")
  summary_columns <- .mtblr_convergence_summary_columns()
  statuses <- c(
    "not_updated", "nonfinite", "constant", "unavailable_single_chain",
    "insufficient_draws", "constant_chain_mismatch", "computed_partial",
    "computed_fewer_than_four_chains", "computed")
  summary <- result$summary
  valid_flags <- if (is.data.frame(summary) &&
                     all(c("status", "rhat", "ess_bulk_per_chain",
                           "ess_tail_per_chain", "mcse_mean_over_sd",
                           "rhat_flag", "ess_bulk_flag", "ess_tail_flag",
                           "mcse_flag") %in% names(summary))) {
    identical(as.logical(summary$rhat_flag),
      summary$status == "constant_chain_mismatch" |
        (is.finite(summary$rhat) &
           summary$rhat > result$thresholds$rhat_threshold)) &&
      identical(as.logical(summary$ess_bulk_flag),
        is.finite(summary$ess_bulk_per_chain) &
          summary$ess_bulk_per_chain <
            result$thresholds$ess_per_chain_threshold) &&
      identical(as.logical(summary$ess_tail_flag),
        is.finite(summary$ess_tail_per_chain) &
          summary$ess_tail_per_chain <
            result$thresholds$ess_per_chain_threshold) &&
      identical(as.logical(summary$mcse_flag),
        is.finite(summary$mcse_mean_over_sd) &
          summary$mcse_mean_over_sd >
            result$thresholds$mcse_mean_over_sd_threshold)
  } else FALSE
  expected_availability <- if (is.data.frame(summary)) list(
    rhat = sum(summary$rhat_available),
    ess_bulk = sum(summary$ess_bulk_available),
    ess_tail = sum(summary$ess_tail_available),
    ess_mean = sum(summary$ess_mean_available),
    mcse_mean = sum(summary$mcse_mean_available)) else NULL
  not_requested <- is.list(result) && identical(result$scope, "none") &&
    identical(result$overall_status, "not_requested")
  overview_consistent <- if (not_requested && is.data.frame(summary) &&
                             nrow(summary) == 0L) {
    identical(.mtblr_convergence_empty_overview(result$nchains),
              result$overview)
  } else if (is.data.frame(summary) && nrow(summary)) {
    identical(.mtblr_convergence_overview(summary), result$overview)
  } else FALSE
  if (!is.list(result) || !all(required %in% names(result)) ||
      !identical(as.integer(result$version), 1L) ||
      !result$scope %in% c("none", "core") ||
      !result$overall_status %in% c(
        "not_requested", "unavailable", "ok", "warning", "partial") ||
      !is.data.frame(result$summary) ||
      !identical(names(result$summary), summary_columns) ||
      any(!result$summary$status %in% statuses) ||
      any(result$summary$rhat_available & !is.finite(result$summary$rhat)) ||
      any(result$summary$ess_bulk_available &
          (!is.finite(result$summary$ess_bulk) |
           result$summary$ess_bulk < 0)) ||
      any(result$summary$ess_tail_available &
          (!is.finite(result$summary$ess_tail) |
           result$summary$ess_tail < 0)) ||
      any(result$summary$mcse_mean_available &
          (!is.finite(result$summary$mcse_mean) |
           result$summary$mcse_mean < 0)) ||
      !valid_flags ||
      !identical(result$availability, expected_availability) ||
      !overview_consistent ||
      !identical(as.integer(result$total_postburn_draws),
                 as.integer(result$nchains *
                              result$postburn_draws_per_chain)) ||
      !identical(as.integer(result$split_draws_per_chain),
                 as.integer(result$postburn_draws_per_chain %/% 2L)) ||
      !identical(result$overall_status, result$overview$overall_status) ||
      !is.logical(result$trace_retention$retained) ||
      length(result$trace_retention$retained) != 1L ||
      is.na(result$trace_retention$retained) ||
      !identical(result$trace_retention$burnin_included, FALSE) ||
      !identical(result$trace_retention$additional_thinning, FALSE) ||
      (identical(result$scope, "none") &&
       (!identical(result$requested, FALSE) ||
        !identical(result$computed, FALSE) || nrow(summary) != 0L ||
        !identical(result$trace_retention$retained, FALSE))) ||
      (!isTRUE(result$requested) && !identical(result$scope, "none")) ||
      (identical(result$overall_status, "not_requested") &&
       !identical(result$scope, "none"))) {
    stop("Invalid MT BLR convergence result.", call. = FALSE)
  }
  result
}

.mtblr_convergence_format_value <- function(value) {
  if (length(value) != 1L || !is.finite(value)) "NA" else
    formatC(value, digits = 6L, format = "fg", flag = "#")
}

.mtblr_convergence_format_quantity <- function(value) {
  if (length(value) != 1L || is.na(value) || !nzchar(value)) "NA" else value
}

.mtblr_convergence_warning_messages <- function(convergence,
                                                mode = c("auto", "core")) {
  mode <- match.arg(mode)
  convergence <- .mtblr_validate_convergence_result(convergence)
  overview <- convergence$overview
  advisory <- overview$n_flagged > 0L ||
    identical(convergence$overall_status, "partial") ||
    (identical(mode, "core") &&
     identical(convergence$overall_status, "unavailable"))
  if (!advisory) return(character())
  if (identical(convergence$overall_status, "unavailable")) {
    reasons <- unique(convergence$summary$status[
      convergence$summary$status != "not_updated"])
    return(sprintf(
      paste0("MT BED Tier 1 convergence diagnostics are unavailable: ",
             "nchains=%d, nit=%d, reason=%s; review fit$convergence."),
      convergence$nchains, convergence$postburn_draws_per_chain,
      if (length(reasons)) paste(reasons, collapse = ",") else "unavailable"))
  }
  sprintf(
    paste0("MT BED Tier 1 convergence advisory requires review: status=%s; ",
           "flagged quantities=%d; max R-hat=%s (%s); ",
           "min bulk ESS/chain=%s (%s); min tail ESS/chain=%s (%s); ",
           "max relative MCSE=%s (%s); constant-chain mismatches=%d; ",
           "fewer than four chains=%s; review fit$convergence."),
    convergence$overall_status, overview$n_flagged,
    .mtblr_convergence_format_value(overview$max_rhat),
    .mtblr_convergence_format_quantity(overview$max_rhat_quantity),
    .mtblr_convergence_format_value(overview$min_ess_bulk_per_chain),
    .mtblr_convergence_format_quantity(overview$min_ess_bulk_quantity),
    .mtblr_convergence_format_value(overview$min_ess_tail_per_chain),
    .mtblr_convergence_format_quantity(overview$min_ess_tail_quantity),
    .mtblr_convergence_format_value(overview$max_mcse_mean_over_sd),
    .mtblr_convergence_format_quantity(overview$max_mcse_quantity),
    overview$n_constant_chain_mismatch,
    if (isTRUE(overview$fewer_than_four_chains)) "yes" else "no")
}

.mtblr_bed_convergence_internal <- function(native_result, trait_names,
                                             updateB, updateE,
                                             control = NULL,
                                             keep_traces = FALSE) {
  if (!is.list(native_result) ||
      !identical(names(native_result), c("raw", "trace_bundle"))) {
    stop("Expected an MT BED native convergence-trace result.",
         call. = FALSE)
  }
  raw <- .validate_mtblr_raw(native_result$raw)
  bundle <- .mtblr_validate_convergence_trace_bundle(
    native_result$trace_bundle, nt = length(trait_names),
    updateB = updateB, updateE = updateE)
  convergence <- .mtblr_convergence_tier1(
    bundle, trait_names, control = control, keep_traces = keep_traces)
  raw$diagnostics$convergence <- convergence
  traces <- NULL
  if (isTRUE(keep_traces)) {
    values <- bundle$values
    quantity_names <- convergence$summary$quantity
    dimnames(values) <- list(
      paste0("Iter", seq_len(dim(values)[1L])),
      paste0("chain", seq_len(dim(values)[2L])),
      quantity_names)
    traces <- list(
      scope = "core",
      postburn_draws_per_chain =
        as.integer(bundle$postburn_draws_per_chain),
      quantities = transform(
        bundle$quantities, quantity = quantity_names),
      values = values)
  }
  list(raw = raw, convergence_traces = traces)
}

.mtblr_convergence_memory_estimate <- function(nchains, nit, nt,
                                                keep_traces = FALSE) {
  values <- c(nchains = nchains, nit = nit, nt = nt)
  if (any(lengths(as.list(values)) != 1L) ||
      any(!is.finite(values)) || any(values < 1) ||
      any(values != floor(values)) ||
      !is.logical(keep_traces) || length(keep_traces) != 1L ||
      is.na(keep_traces)) {
    stop("nchains, nit, and nt must be positive integers and keep_traces logical.",
         call. = FALSE)
  }
  trace_bytes <- 8 * nchains * nit * 3 * nt
  # One quantity: raw/split/rank/autocovariance/rho working arrays.
  workspace <- 8 * nchains * nit * 6
  summary <- 8 * 35 * 3 * nt
  retained <- if (keep_traces) trace_bytes else 0
  total <- trace_bytes + workspace + summary + retained
  list(
    trace_capture_bytes = trace_bytes,
    workspace_bytes_per_quantity = workspace,
    maximum_workspace_bytes = workspace,
    summary_output_bytes = summary,
    retained_trace_bytes = retained,
    estimated_total_bytes = total,
    estimated_total_gib = total / 1024^3,
    measured_rss = FALSE,
    measured_peak_rss = FALSE
  )
}
