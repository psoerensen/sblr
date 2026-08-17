# Corrected Cheng MT-BayesC-Pi workflow for complete common-sample packed-BED
# data. The small iteration counts illustrate the public interface; they are
# not convergence recommendations. Phenotypes must already be centred or
# residualized and must have stable row and column names.
stopifnot(
  is.matrix(y), is.numeric(y), ncol(y) >= 2L,
  !is.null(rownames(y)), !is.null(colnames(y)), all(is.finite(y))
)

trait_ids <- colnames(y)
trait_count <- ncol(y)
residual_variances <- pmax(apply(y, 2L, stats::var), .Machine$double.eps)
fixed_residual_covariance <- diag(residual_variances, trait_count)
dimnames(fixed_residual_covariance) <- list(trait_ids, trait_ids)

initial_marker_covariance <- 0.05 * fixed_residual_covariance
marker_covariance_prior_scale <- initial_marker_covariance

fit_bed <- mtblr_bed(
  y = y,
  Glist = Glist,
  trait_ids = trait_ids,
  sample_ids = rownames(y),
  method = "bayesc",
  residual_covariance_policy = "fixed_full",
  fixed_residual_covariance = fixed_residual_covariance,
  initial_marker_covariance = initial_marker_covariance,
  marker_covariance_prior_degrees_of_freedom = trait_count + 2,
  marker_covariance_prior_scale = marker_covariance_prior_scale,
  nit = 20L,
  nburn = 10L,
  nthin = 2L,
  nchains = 1L,
  ncores = 1L,
  seed = 17
)

extract_posterior(fit_bed, "pips")
extract_posterior(fit_bed, "realised_effects")
extract_posterior(fit_bed, "activity_pattern_probabilities")
extract_posterior(fit_bed, "pleiotropic_probabilities")
extract_posterior(fit_bed, "effect_covariance")
extract_posterior(fit_bed, "residual_covariance", state = "final")
