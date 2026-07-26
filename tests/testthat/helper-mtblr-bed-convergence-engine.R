phase17u_args <- function(case, nchains = 2L, ncores = 1L,
                          chain_seeds = integer(),
                          keep_chains = FALSE) {
  phase17r_args(case, nchains = nchains, ncores = ncores,
                chain_seeds = chain_seeds, keep_chains = keep_chains)
}

phase17u_native_call <- function(case, nchains = 2L, ncores = 1L,
                                 chain_seeds = integer(),
                                 keep_chains = FALSE) {
  do.call(sblr:::mtblr_bed_convergence_trace_internal,
          phase17u_args(case, nchains, ncores, chain_seeds, keep_chains))
}

phase17u_diagnose <- function(native, trait_names, updateB = TRUE,
                              updateE = TRUE, control = NULL,
                              keep_traces = FALSE) {
  sblr:::.mtblr_bed_convergence_internal(
    native, trait_names, updateB, updateE, control, keep_traces)
}

phase17u_without_timing <- function(x) {
  if (!is.null(x$raw)) x$raw <- phase17r_without_timing(x$raw)
  x
}

phase17u_bundle_from_chains <- function(chains, updated = rep(TRUE, 3L),
                                        group = c("B_diag", "G_diag",
                                                  "E_diag")) {
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

phase17u_fixture_metrics <- function(chains) {
  sblr:::.mtblr_convergence_scalar(t(chains))
}
