# Independent small-matrix reference for complete Cheng activity patterns.
# This file intentionally does not call package production helpers.

mtblr_general_t_patterns <- function(trait_ids) {
  stopifnot(is.character(trait_ids), length(trait_ids) >= 2L,
            !anyNA(trait_ids), !anyDuplicated(trait_ids))
  grid <- expand.grid(rep(list(0:1), length(trait_ids)),
                      KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  out <- as.matrix(grid)
  storage.mode(out) <- "integer"
  dimnames(out) <- list(apply(out, 1L, paste, collapse = "_"), trait_ids)
  out
}

mtblr_general_t_marker_reference <- function(score, marker_sum_squares,
                                              Vb, Ve, probability,
                                              patterns) {
  trait_count <- length(score)
  stopifnot(identical(dim(Vb), c(trait_count, trait_count)),
            identical(dim(Ve), c(trait_count, trait_count)),
            ncol(patterns) == trait_count,
            nrow(patterns) == length(probability))
  omega <- solve(Ve)
  h <- as.numeric(omega %*% score)
  log_weight <- log(probability)
  means <- covariances <- vector("list", nrow(patterns))
  for (state in seq_len(nrow(patterns))[-1L]) {
    active <- which(patterns[state, ] == 1L)
    prior <- Vb[active, active, drop = FALSE]
    precision <- solve(prior) + marker_sum_squares *
      omega[active, active, drop = FALSE]
    covariance <- solve(precision)
    mean <- as.numeric(covariance %*% h[active])
    log_weight[state] <- log_weight[state] -
      0.5 * as.numeric(determinant(prior, logarithm = TRUE)$modulus) -
      0.5 * as.numeric(determinant(precision, logarithm = TRUE)$modulus) +
      0.5 * sum(h[active] * mean)
    means[[state]] <- mean
    covariances[[state]] <- covariance
  }
  weight <- exp(log_weight - max(log_weight))
  list(log_weight = log_weight, probability = weight / sum(weight),
       active_mean = means, active_covariance = covariances)
}

mtblr_general_t_completion_reference <- function(Vb, pattern, active_value) {
  active <- which(pattern == 1L)
  inactive <- which(pattern == 0L)
  stopifnot(length(active) > 0L, length(active) == length(active_value))
  if (!length(inactive)) {
    return(list(coefficient = matrix(numeric(), 0L, length(active)),
                mean = numeric(), covariance = matrix(numeric(), 0L, 0L)))
  }
  Vaa <- Vb[active, active, drop = FALSE]
  Via <- Vb[inactive, active, drop = FALSE]
  coefficient <- Via %*% solve(Vaa)
  covariance <- Vb[inactive, inactive, drop = FALSE] -
    coefficient %*% Vb[active, inactive, drop = FALSE]
  list(coefficient = coefficient,
       mean = as.numeric(coefficient %*% active_value),
       covariance = 0.5 * (covariance + t(covariance)))
}

mtblr_general_t_covariance_conditionals <- function(
    Y, X, realised_effects, latent_effects, states,
    marker_prior_df, marker_prior_scale,
    residual_prior_df, residual_prior_scale) {
  residual <- Y - X %*% realised_effects
  active <- states != 0L
  marker_statistic <- if (any(active)) {
    crossprod(latent_effects[active, , drop = FALSE])
  } else matrix(0, ncol(Y), ncol(Y))
  list(
    residual = residual,
    marker_degrees_of_freedom = marker_prior_df + sum(active),
    marker_scale = marker_prior_scale + marker_statistic,
    residual_degrees_of_freedom = residual_prior_df + nrow(Y),
    residual_scale = residual_prior_scale + crossprod(residual))
}
