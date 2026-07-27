.mtblr_summary_chain_seeds <- function(chain) {
  if (!is.null(chain$chain_seeds_requested)) {
    return(as.numeric(chain$chain_seeds_requested))
  }
  seeds <- (as.numeric(chain$seed) + 9176 * (seq_len(chain$nchains) - 1)) %%
    4294967296
  ifelse(seeds > .Machine$integer.max, seeds - 4294967296, seeds)
}

.mtblr_summary_chain_record <- function(raw, chain, seed, model, operator) {
  list(
    family = "mtblr", model = model, operator = operator,
    trait_index = NA_integer_, trait_name = NA_character_,
    chain_index = as.integer(chain), seed = seed,
    retained_draw_count = as.integer(raw$diagnostics$marker),
    marker = raw$marker[c("bm", "dm", "b", "state")],
    trace = raw$trace[c("vbs", "vgs", "ves", "vle", "vld")],
    variance = list(
      cov_b_mean = raw$variance$covb,
      cov_g_mean = raw$variance$covg,
      cov_e_mean = raw$variance$cove,
      cov_b_final = raw$variance$vb,
      cov_g_final = raw$variance$vg,
      cov_e_final = raw$variance$ve),
    pi = list(pi_final = raw$pi$final, pi_mean = raw$pi$mean))
}

.mtblr_summary_pool_raw <- function(raws, keep_chains, seeds, operator) {
  if (!length(raws)) stop("At least one MTBLR chain is required.", call. = FALSE)
  raws <- lapply(raws, .validate_mtblr_raw)
  out <- raws[[1L]]
  nchains <- length(raws)
  weights <- vapply(raws, function(x) as.numeric(x$diagnostics$marker),
                    numeric(1))
  if (any(!is.finite(weights) | weights <= 0)) {
    stop("Every MTBLR chain must retain posterior marker draws.",
         call. = FALSE)
  }
  weighted_matrix <- function(path, weight = weights) {
    values <- lapply(raws, function(x) x[[path[[1L]]]][[path[[2L]]]])
    if (sum(weight) <= 0) return(values[[1L]])
    Reduce(`+`, Map(`*`, values, weight)) / sum(weight)
  }
  simple_mean <- function(path) {
    values <- lapply(raws, function(x) x[[path[[1L]]]][[path[[2L]]]])
    Reduce(`+`, values) / nchains
  }
  out$marker$bm <- weighted_matrix(c("marker", "bm"))
  out$marker$dm <- weighted_matrix(c("marker", "dm"))
  for (name in c("vbs", "vgs", "ves", "vle", "vld")) {
    out$trace[[name]] <- simple_mean(c("trace", name))
  }
  for (name in c("covb", "covg", "cove")) {
    count <- vapply(raws, function(x) as.numeric(x$diagnostics[[name]]),
                    numeric(1))
    out$variance[[name]] <- weighted_matrix(c("variance", name), count)
  }
  pi_count <- vapply(raws, function(x) as.numeric(x$diagnostics$pi), numeric(1))
  out$pi$mean <- weighted_matrix(c("pi", "mean"), pi_count)

  bm <- simplify2array(lapply(raws, function(x) x$marker$bm))
  dm <- simplify2array(lapply(raws, function(x) x$marker$dm))
  if (nchains == 1L) {
    out$marker$bm_sd <- out$marker$bm * 0
    out$marker$dm_sd <- out$marker$dm * 0
  } else {
    out$marker$bm_sd <- apply(bm, c(1L, 2L), stats::sd)
    out$marker$dm_sd <- apply(dm, c(1L, 2L), stats::sd)
  }
  out$marker$bm_min <- apply(bm, c(1L, 2L), min)
  out$marker$bm_max <- apply(bm, c(1L, 2L), max)
  out$marker$dm_min <- apply(dm, c(1L, 2L), min)
  out$marker$dm_max <- apply(dm, c(1L, 2L), max)
  out$diagnostics$multichain <- list(
    topology = "one_complete_joint_mt_model_per_chain",
    aggregation = "pooled_retained_draws",
    trace_policy = "iterationwise_chain_mean",
    final_state_policy = "primary_chain_1",
    primary_chain = 1L, nchains = as.integer(nchains),
    requested_cores = NA_integer_, used_workers = 1L,
    chain_seeds = seeds,
    implementation = "serial_R_orchestration_over_validated_native_chain")
  out$chains <- NULL
  records <- if (keep_chains) lapply(seq_along(raws), function(i) {
    .mtblr_summary_chain_record(raws[[i]], i, seeds[[i]], "bayesc", operator)
  }) else NULL
  if (!is.null(records)) names(records) <- paste0("chain", seq_along(records))
  list(raw = out, chains = records)
}

.mtblr_summary_convergence_bundle <- function(raws, trait_names, model,
                                               operator, nit, nburn,
                                               updateB, updateE) {
  nt <- length(trait_names)
  groups <- rep(c("vbs", "vgs", "ves", "vle", "vld"), each = nt)
  fields <- rep(c("vbs", "vgs", "ves", "vle", "vld"), each = nt)
  traits <- rep(seq_len(nt), 5L)
  values <- array(NA_real_, c(nit, length(raws), 5L * nt))
  for (q in seq_len(5L * nt)) for (chain in seq_along(raws)) {
    trace <- raws[[chain]]$trace[[fields[[q]]]][, traits[[q]]]
    values[, chain, q] <- trace[seq.int(nburn + 1L, nburn + nit)]
  }
  quantities <- data.frame(
    quantity_index = seq_len(5L * nt), group = groups,
    trait_index = traits,
    updated = c(rep(updateB, nt), rep(TRUE, nt), rep(updateE, nt),
                rep(TRUE, 2L * nt)),
    derived = c(rep(FALSE, nt), rep(TRUE, nt), rep(FALSE, nt),
                rep(TRUE, 2L * nt)))
  .blr_convergence_bundle(values, quantities, "mtblr", model, operator)
}

.mtblr_summary_multichain <- function(native_execution, chain, conv,
                                      trait_names, model, operator,
                                      updateB, updateE) {
  if (!is.list(native_execution) ||
      length(native_execution$raws) != chain$nchains ||
      !identical(as.integer(native_execution$operator_preparations), 1L)) {
    stop("MT summary native execution must return one prepared chain set.",
         call. = FALSE)
  }
  seeds <- as.numeric(native_execution$chain_seeds)
  raws <- native_execution$raws
  pooled <- .mtblr_summary_pool_raw(raws, chain$keep_chains, seeds, operator)
  pooled$raw$diagnostics$multichain$requested_cores <- chain$ncores
  pooled$raw$diagnostics$multichain$used_workers <-
    as.integer(native_execution$used_workers)
  pooled$raw$diagnostics$multichain$implementation <-
    "native_static_chain_dispatch_shared_operator_preparation"
  pooled$raw$diagnostics$multichain$operator_preparations <- 1L
  traces <- NULL
  if (conv$compute || conv$keep_traces) {
    bundle <- .mtblr_summary_convergence_bundle(
      raws, trait_names, model, operator, chain$nit, chain$nburn,
      updateB, updateE)
    convergence <- .blr_convergence_tier1(
      bundle, trait_names, conv$thresholds, conv$keep_traces)
    if (conv$keep_traces) {
      traces <- bundle
      dimnames(traces$values) <- list(
        paste0("Iter", seq_len(chain$nit)),
        paste0("chain", seq_len(chain$nchains)),
        convergence$summary$quantity)
    }
  } else if (conv$mode == "none") {
    convergence <- .blr_convergence_not_requested(
      trait_names, updateB, updateE, chain$nchains, chain$nit,
      conv$thresholds)
  } else {
    convergence <- .blr_convergence_unavailable(
      trait_names, updateB, updateE, chain$nchains, chain$nit,
      conv$thresholds)
  }
  pooled$raw$diagnostics$convergence <- convergence
  c(pooled, list(convergence = convergence, convergence_traces = traces,
                 seeds = seeds,
                 used_workers = as.integer(native_execution$used_workers)))
}
