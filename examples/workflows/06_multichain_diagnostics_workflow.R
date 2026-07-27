## 06_multichain_diagnostics_workflow.R
##
## Multi-chain and chain-summary diagnostics workflow.
##
## Demonstrates:
##   - nchains, chain_seeds, keep_chains
##   - stable single-chain and multi-chain marker summaries
##   - fit$chains and fit$ld_swap_chains
##   - check_stblr_consistency(require_chain_summaries = TRUE)

library(sblr)
source("./examples/workflows/00_workflow_helpers.R")

chr <- 1L
nthreads <- 4L
seed <- 10L

data_dir <- workflow_data_dir()
Glist <- read_example_glist()

sim <- mtsim(
  Glist = Glist,
  chr = chr,
  rsids = Glist$rsidsLD[[chr]],
  nt = 3,
  n_shared = 30,
  n_specific = 10,
  h2 = c(0.4, 0.5, 0.3),
  rg = matrix(c(1.0, 0.7, 0.3, 0.7, 1.0, 0.5, 0.3, 0.5, 1.0), nrow = 3),
  re = 0,
  seed = 1
)
y <- as.matrix(scale(sim$y))

stats <- make_summary_stats(Glist = Glist, y = y, chr = chr, nthreads = nthreads)
Glist <- make_sparse_ld(
  Glist = Glist,
  out_prefix = file.path(data_dir, "ld_multichain"),
  chr = chr,
  pos_bp = NULL,
  max_distance_bp = 0,
  max_distance_variants = 1000,
  r2_threshold = 0.001,
  block_size = 1024,
  nthreads = nthreads
)

## -------------------------------------------------------------------------
## 1. Single-chain fit
## -------------------------------------------------------------------------
## For nchains == 1, marker chain-summary fields are represented as degenerate
## summaries: *_sd = 0 and *_min == *_max == corresponding mean. This is a
## formatter convention for stable fit structure, not a sampler change.

fit_single <- stblr_csr(
  stats = stats,
  Glist = Glist,
  method = "bayesc",
  nit = 500,
  nburn = 100,
  ncores = nthreads,
  seed = seed,
  nchains = 1
)

check_fit(fit_single, require_chain_summaries = TRUE)

stopifnot(all(fit_single$dm_chain_mean_sd == 0, na.rm = TRUE))
stopifnot(all(fit_single$bm_chain_mean_sd == 0, na.rm = TRUE))
stopifnot(isTRUE(all.equal(fit_single$dm_chain_mean_min, fit_single$dm)))
stopifnot(isTRUE(all.equal(fit_single$dm_chain_mean_max, fit_single$dm)))
stopifnot(isTRUE(all.equal(fit_single$bm_chain_mean_min, fit_single$bm)))
stopifnot(isTRUE(all.equal(fit_single$bm_chain_mean_max, fit_single$bm)))

## -------------------------------------------------------------------------
## 2. Multi-chain fit
## -------------------------------------------------------------------------

fit_multi <- stblr_csr(
  stats = stats,
  Glist = Glist,
  method = "bayesc",
  nit = 500,
  nburn = 100,
  ncores = nthreads,
  seed = seed,
  nchains = 2,
  chain_seeds = c(101, 202),
  keep_chains = TRUE
)

check_fit(fit_multi, require_chain_summaries = TRUE)

length(fit_multi$chains)
summary(as.vector(fit_multi$dm_chain_mean_sd))
summary(as.vector(fit_multi$bm_chain_mean_sd))

## -------------------------------------------------------------------------
## 3. Multi-chain LD-swap diagnostics
## -------------------------------------------------------------------------

fit_multi_ldswap <- stblr_csr(
  stats = stats,
  Glist = Glist,
  method = "bayesc",
  nit = 500,
  nburn = 100,
  ncores = nthreads,
  seed = seed,
  nchains = 2,
  chain_seeds = c(101, 202),
  keep_chains = TRUE,
  updateLDswap = TRUE,
  ld_swap_prob = 0.10,
  ld_swap_r2 = 0.05,
  ld_swap_moves = 5
)

check_fit(
  fit_multi_ldswap,
  require_chain_summaries = TRUE,
  require_ld_swap = TRUE
)

fit_multi_ldswap$ld_swap
fit_multi_ldswap$ld_swap_chains

fits <- list(
  single_chain = fit_single,
  multi_chain = fit_multi,
  multi_chain_ldswap = fit_multi_ldswap
)

fit_metadata_table(fits)
fit_field_inventory(fits)
ld_swap_summary_table(fits)

## -------------------------------------------------------------------------
## 4. Optional multi-chain BayesR example
## -------------------------------------------------------------------------

if (workflow_should_run_heavy()) {
  fit_bayesr_multi <- stblr_csr(
    stats = stats,
    Glist = Glist,
    method = "bayesr",
    nit = 1000,
    nburn = 250,
    ncores = nthreads,
    seed = seed,
    nchains = 2,
    chain_seeds = c(303, 404),
    keep_chains = TRUE
  )

  check_fit(fit_bayesr_multi, require_chain_summaries = TRUE)
  summarise_components(fit_bayesr_multi)
}

invisible(fits)
