# canonical MT BayesRC/SBayesRC workflow (small inputs assumed).
# `annotations` must carry explicit marker IDs in row names or a marker_id
# column. Mixed binary and continuous columns share one preprocessing contract.

annotation <- data.frame(
  marker_id = stats$marker_names,
  coding = as.integer(seq_along(stats$marker_names) %% 2L),
  score = seq(-1, 1, length.out = length(stats$marker_names))
)

fit_csr <- mtblr_csr(
  stats, ld_prefix = ld_prefix, ld_metadata = ld_metadata,
  method = "sbayesrc", annotations = annotation,
  mixture_var = c(0, 0.01, 0.1, 1),
  annotation_intercept_prior = list(
    distribution = "normal", mean = "initial_mixture", sd = 1),
  nchains = 2L, ncores = 2L, convergence = "auto"
)

fit_block <- mtblr_block_eigen(
  stats, Glist, block_start = 1L, method = "sbayesrc",
  annotations = annotation, mixture_var = c(0, 0.01, 0.1, 1),
  annotation_intercept_prior = list(mean = "initial_mixture", sd = 1)
)

fit_bed <- mtblr_bed(
  y, Glist, method = "bayesrc", annotations = annotation,
  mixture_var = c(0, 0.01, 0.1, 1),
  annotation_intercept_prior = list(mean = "initial_mixture", sd = 1),
  maf_effect_s = -1
)

# Posterior state probabilities and annotation-driven prior probabilities are
# deliberately different objects.
fit_csr$component_probabilities
fit_csr$model_parameters$annotations$prior_component_probabilities
fit_csr$model_parameters$annotations$pattern_pi_mean
fit_csr[c("vbs", "vgs", "ves", "vle", "vld")]
fit_csr$convergence
