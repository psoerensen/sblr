phase17r_seed_oracle <- function(seed, nchains) {
  base <- as.double(seed) %% 2^32
  (base + 9176 * (seq_len(nchains) - 1)) %% 2^32
}

phase17r_args <- function(case, nchains = 1L, ncores = 1L,
                          chain_seeds = integer(), keep_chains = FALSE) {
  c(phase17o_args(case), list(
    nchains = as.integer(nchains), ncores = as.integer(ncores),
    chain_seeds = as.integer(chain_seeds), keep_chains = keep_chains,
    joint_component = integer(), joint_multiplier = numeric(),
    joint_names = character(), component_count = 0L,
    marker_scale = numeric(), pi_prior = numeric(),
    component_init = integer()
  ))
}

phase17r_call <- function(case, nchains = 1L, ncores = 1L,
                          chain_seeds = integer(), keep_chains = FALSE) {
  do.call(sblr:::mtblr_bed_chains_internal,
          phase17r_args(case, nchains, ncores, chain_seeds, keep_chains))
}

phase17r_numerical_raw <- function(raw) {
  list(
    marker = raw$marker[c("bm", "dm", "wy", "r", "b", "state", "order")],
    trace = raw$trace,
    variance = raw$variance,
    pi = raw$pi,
    diagnostics = raw$diagnostics[c("marker", "covb", "covg", "cove", "pi")]
  )
}

phase17r_without_timing <- function(raw) {
  bed <- raw$diagnostics$mt_bed
  bed[c("requested_cores", "used_workers", "chain_seconds", "seconds_mean",
        "seconds_max", "dispatch_seconds")] <- NULL
  raw$diagnostics$mt_bed <- bed
  if (!is.null(raw$chains)) {
    raw$chains <- lapply(raw$chains, function(chain) {
      chain$diagnostics$seconds <- NULL
      chain
    })
  }
  raw
}
