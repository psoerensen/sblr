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
