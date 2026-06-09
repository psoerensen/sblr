# Basic sblr workflow using summary statistics
#
# This is a template because the package does not yet include a small example
# dataset. Replace the placeholder inputs below with summary statistics and LD
# matrices from your study before running the fitting call.

library(sblr)

# 1. Prepare or load summary statistics for each trait:
#    Xy: X'y vectors, one vector per trait.
#    yy: y'y values, one value per trait.
#    n: sample sizes, one value per trait.
#
# Xy <- list(trait1 = xy_trait1, trait2 = xy_trait2)
# yy <- c(trait1 = yy_trait1, trait2 = yy_trait2)
# n <- c(trait1 = n_trait1, trait2 = n_trait2)

# 2. Prepare LD-derived X'X matrices.
#    Each matrix should correspond to the markers and ordering used in Xy.
#    Its diagonal contains marker X'X values; off-diagonal entries represent LD.
#
# XX <- list(trait1 = xtx_trait1, trait2 = xtx_trait2)

# 3. Fit the model. Adjust prior settings and MCMC lengths for the analysis.
#
# fit <- sblr(
#   yy = yy,
#   Xy = Xy,
#   XX = XX,
#   n = n,
#   h2 = 0.5,
#   pi = 0.001,
#   method = "bayesC",
#   nit = 1000,
#   nburn = 500,
#   verbose = FALSE
# )

# 4. Inspect posterior summaries.
#
# fit$bm    # posterior mean marker effects
# fit$dm    # posterior inclusion probabilities
# fit$vg    # genetic variance estimates
# fit$ve    # residual variance estimates
# fit$covg  # genetic covariance estimates
