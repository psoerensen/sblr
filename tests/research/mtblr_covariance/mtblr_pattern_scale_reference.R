# Independent Phase 7 pattern-by-scale reference.
# This file deliberately uses only base R linear algebra and does not call
# package or candidate-production helpers.

mtblr_pattern_scale_states <- function(patterns, scales) {
  patterns <- as.matrix(patterns)
  null <- which(rowSums(patterns) == 0L)
  if (length(null) != 1L || null != 1L || any(scales <= 0) ||
      any(!is.finite(scales))) stop("invalid pattern-by-scale declaration")
  active <- which(rowSums(patterns) > 0L)
  data.frame(
    pattern = c(null, rep(active, each = length(scales))),
    scale = c(NA_integer_, rep(seq_along(scales), times = length(active))))
}

mtblr_pattern_scale_conditional <- function(
    h, likelihood_precision, Vb, patterns, pattern_probability,
    scales, scale_probability, q = 1) {
  patterns <- as.matrix(patterns)
  states <- mtblr_pattern_scale_states(patterns, scales)
  if (length(h) != nrow(Vb) || any(dim(likelihood_precision) != dim(Vb)) ||
      length(pattern_probability) != nrow(patterns) ||
      length(scale_probability) != length(scales) || q <= 0) {
    stop("incompatible conditional inputs")
  }
  pattern_probability <- pattern_probability / sum(pattern_probability)
  scale_probability <- scale_probability / sum(scale_probability)
  means <- covariances <- vector("list", nrow(states))
  log_weight <- numeric(nrow(states))
  log_weight[1L] <- log(pattern_probability[1L])
  for (state in 2:nrow(states)) {
    pattern <- states$pattern[state]
    scale <- states$scale[state]
    active <- which(patterns[pattern, ] == 1L)
    prior <- scales[scale] * q * Vb[active, active, drop = FALSE]
    precision <- solve(prior) +
      likelihood_precision[active, active, drop = FALSE]
    covariance <- solve(precision)
    mean <- as.numeric(covariance %*% h[active])
    means[[state]] <- mean
    covariances[[state]] <- covariance
    log_weight[state] <- log(pattern_probability[pattern]) +
      log(scale_probability[scale]) - determinant(prior, TRUE)$modulus / 2 -
      determinant(precision, TRUE)$modulus / 2 +
      drop(crossprod(h[active], mean)) / 2
  }
  probability <- exp(log_weight - max(log_weight))
  probability <- probability / sum(probability)
  list(states = states, log_weight = log_weight, probability = probability,
       active_mean = means, active_covariance = covariances)
}

mtblr_pattern_scale_completion <- function(active_value, active, Vb,
                                            scale, q = 1) {
  inactive <- setdiff(seq_len(nrow(Vb)), active)
  if (!length(inactive)) {
    return(list(mean = numeric(), covariance = matrix(numeric(), 0L, 0L)))
  }
  aa <- Vb[active, active, drop = FALSE]
  ia <- Vb[inactive, active, drop = FALSE]
  ai <- Vb[active, inactive, drop = FALSE]
  ii <- Vb[inactive, inactive, drop = FALSE]
  list(
    mean = as.numeric(ia %*% solve(aa, active_value)),
    covariance = scale * q * (ii - ia %*% solve(aa, ai)))
}

mtblr_pattern_scale_vb_statistic <- function(theta, pattern, scale, q) {
  theta <- as.matrix(theta)
  pattern <- as.integer(pattern)
  active <- which(pattern != 1L)
  if (!length(active)) return(matrix(0, ncol(theta), ncol(theta)))
  out <- matrix(0, ncol(theta), ncol(theta))
  for (marker in active) {
    out <- out + tcrossprod(theta[marker, ]) / (scale[marker] * q[marker])
  }
  out
}
