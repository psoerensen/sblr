phase17v_public_args <- function(case, nchains = 1L, ncores = 1L,
                                 convergence = "auto", warn = FALSE,
                                 keep_traces = FALSE, keep_chains = FALSE,
                                 nit = 8L, nburn = 3L, nthin = 1L, ...) {
  phase17s_public_args(
    case, nchains = nchains, ncores = ncores, keep_chains = keep_chains,
    convergence = convergence, nit = nit, nburn = nburn, nthin = nthin,
    convergence_control = list(warn = warn, keep_traces = keep_traces), ...)
}

phase17v_without_additions <- function(fit) {
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

phase17v_internal_diagnostic <- function(args, trait_names,
                                         updateB = FALSE,
                                         updateE = FALSE,
                                         keep_traces = FALSE) {
  native <- do.call(
    sblr:::mtblr_bed_convergence_trace_internal,
    phase17s_internal_args(args))
  sblr:::.mtblr_bed_convergence_internal(
    native, trait_names, updateB, updateE,
    control = sblr:::.mtblr_convergence_control(),
    keep_traces = keep_traces)
}

phase17v_warning_fixture <- function(status = "computed", flagged = FALSE,
                                     nchains = 4L, nit = 10L) {
  result <- sblr:::.mtblr_convergence_unavailable(
    "T1", TRUE, TRUE, nchains, nit, sblr:::.mtblr_convergence_control())
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
  result$overview <- sblr:::.mtblr_convergence_overview(row)
  result$overall_status <- result$overview$overall_status
  result$computed <- any(unlist(result$availability) > 0L)
  sblr:::.mtblr_validate_convergence_result(result)
}
