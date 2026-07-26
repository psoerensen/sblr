# Public joint multivariate packed-BED workflow.
# Supply an existing small BED-backed Glist named `Glist` and phenotype matrix
# `Y` whose row names match Glist individual IDs.

Y <- as.matrix(Y)

fit_full <- mtblr_bed(
  y = Y,
  Glist = Glist,
  center = TRUE,
  residual_covariance = "full",
  models = "restrictive",
  nit = 50,
  nburn = 20,
  seed = 17
)

fit_diagonal <- mtblr_bed(
  y = Y,
  Glist = Glist,
  center = TRUE,
  residual_covariance = "diagonal",
  models = NULL,
  nit = 50,
  nburn = 20,
  seed = 17
)

fit_full$re
fit_full$bed_diagnostics
fit_full$phenotype_preprocessing
fit_full$memory_estimate

# Default single-chain execution preserves all pre-17S numerical fields.
fit1 <- mtblr_bed(
  y = Y, Glist = Glist, nit = 50, nburn = 20, seed = 101
)

fit_serial <- mtblr_bed(
  y = Y, Glist = Glist, nit = 50, nburn = 20,
  nchains = 4, ncores = 1, seed = 101
)
fit_parallel <- mtblr_bed(
  y = Y, Glist = Glist, nit = 50, nburn = 20,
  nchains = 4, ncores = 2, seed = 101
)
fit_retained <- mtblr_bed(
  y = Y, Glist = Glist, nit = 50, nburn = 20,
  nchains = 4, ncores = 2,
  chain_seeds = c(101L, 9277L, 18453L, 27629L),
  keep_chains = TRUE
)
fit_retained$nchains
fit_retained$chain_seeds
fit_retained$chain_diagnostics
fit_retained$bm_sd
fit_retained$dm_sd
fit_retained$chains
fit_retained$memory_estimate

# Tier 1 convergence diagnostics use only post-burn B/G/E diagonal traces.
# These short runs illustrate the interface; they are not convergence advice.
fit_auto <- mtblr_bed(
  y = Y, Glist = Glist, nit = 20, nburn = 10, seed = 101
)
fit_none <- mtblr_bed(
  y = Y, Glist = Glist, nit = 20, nburn = 10, seed = 101,
  convergence = "none"
)
fit_core <- mtblr_bed(
  y = Y, Glist = Glist, nit = 20, nburn = 10, seed = 101,
  nchains = 4, ncores = 2, convergence = "core",
  convergence_control = list(warn = FALSE)
)
fit_custom <- mtblr_bed(
  y = Y, Glist = Glist, nit = 20, nburn = 10, seed = 101,
  nchains = 4, convergence = "core",
  convergence_control = list(
    warn = FALSE, rhat_threshold = 1.01,
    ess_per_chain_threshold = 100,
    mcse_mean_over_sd_threshold = 0.05,
    keep_traces = TRUE)
)
fit_auto$convergence
fit_core$convergence$summary
fit_core$convergence$overview
fit_custom$convergence_traces
fit_custom$memory_estimate
