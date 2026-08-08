# Tiny exact-state reference for the Phase-4C global selection-model block.
# This is a conditional fixed-z Gibbs kernel, not a genomic production option.

.sbs4c_global_model_transition <- function(model_probability) {
  probability <- as.numeric(model_probability)
  stopifnot(all(is.finite(probability)), all(probability >= 0),
            abs(sum(probability) - 1) < 1e-12)
  matrix(rep(probability, times = length(probability)),
         nrow = length(probability), byrow = TRUE)
}

.sbs4c_detailed_balance_error <- function(model_probability, transition) {
  flow <- diag(as.numeric(model_probability)) %*% transition
  max(abs(flow - t(flow)))
}
