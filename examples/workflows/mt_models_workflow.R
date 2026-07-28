# Small iterations below illustrate the interface; they are not convergence
# recommendations.
fit_csr <- mtblr_csr(stats, ld_prefix = ld_prefix,
  ld_metadata = ld_metadata, method = "sbayesr",
  mixture_var = c(0, 0.01, 0.1, 1), nchains = 2, ncores = 2,
  nit = 100, nburn = 50,
  convergence_control = list(warn = FALSE))

fit_s <- mtblr_csr(stats, ld_prefix = ld_prefix,
  ld_metadata = ld_metadata, method = "sbayesr", maf_effect_s = -1,
  mixture_var = c(0, 0.01, 0.1, 1), nit = 100, nburn = 50)

fit_block <- mtblr_block_eigen(stats, Glist, block_start = 1L,
  method = "sbayesr", mixture_var = c(0, 0.01, 0.1, 1),
  nit = 100, nburn = 50)

fit_bed <- mtblr_bed(y, Glist, method = "bayesr",
  mixture_var = c(0, 0.01, 0.1, 1), nit = 100, nburn = 50)

fit_csr$component_probabilities
fit_csr$pi_mean
fit_csr$model_parameters$mixture
fit_csr[c("vbs", "vgs", "ves", "vle", "vld")]
fit_csr$convergence
