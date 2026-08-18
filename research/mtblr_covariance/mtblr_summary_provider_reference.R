# Independent dense reference for the Phase 6A provider likelihood.

mt_summary_pattern_reference <- function(score, diagonal, marker_covariance,
                                         pattern_probability, patterns) {
  stopifnot(
    is.numeric(score), is.numeric(diagonal), length(score) == ncol(patterns),
    length(diagonal) == length(score), all(diagonal >= 0),
    length(pattern_probability) == nrow(patterns))
  log_weight <- log(pattern_probability)
  active_mean <- vector("list", nrow(patterns))
  active_covariance <- vector("list", nrow(patterns))
  for (state in seq_len(nrow(patterns))) {
    active <- which(patterns[state, ] == 1L)
    if (!length(active)) next
    prior <- marker_covariance[active, active, drop = FALSE]
    precision <- solve(prior) + diag(diagonal[active], length(active))
    covariance <- solve(precision)
    mean <- drop(covariance %*% score[active])
    log_weight[state] <- log_weight[state] -
      determinant(prior, logarithm = TRUE)$modulus / 2 -
      determinant(precision, logarithm = TRUE)$modulus / 2 +
      sum(score[active] * mean) / 2
    active_mean[[state]] <- mean
    active_covariance[[state]] <- covariance
  }
  probability <- exp(log_weight - max(log_weight))
  probability <- probability / sum(probability)
  list(
    log_weight = as.numeric(log_weight), probability = probability,
    active_mean = active_mean, active_covariance = active_covariance)
}

mt_summary_provider_residual <- function(score, cross_product,
                                         local_effect) {
  as.numeric(score - cross_product %*% local_effect)
}

mt_summary_aggregate_marker <- function(providers, global_effect, marker,
                                        trait_count) {
  score <- diagonal <- numeric(trait_count)
  residuals <- vector("list", length(providers))
  for (index in seq_along(providers)) {
    provider <- providers[[index]]
    local_effect <- global_effect[provider$map, provider$trait]
    residual <- mt_summary_provider_residual(
      provider$score, provider$cross_product, local_effect)
    residuals[[index]] <- residual
    local <- match(marker, provider$map)
    if (is.na(local)) next
    restored <- residual[local] +
      provider$cross_product[local, local] * global_effect[marker, provider$trait]
    score[provider$trait] <- score[provider$trait] +
      restored / provider$residual_scale
    diagonal[provider$trait] <- diagonal[provider$trait] +
      provider$cross_product[local, local] / provider$residual_scale
  }
  list(score = score, diagonal = diagonal, residuals = residuals)
}

mt_summary_log_likelihood <- function(providers, global_effect) {
  sum(vapply(providers, function(provider) {
    effect <- global_effect[provider$map, provider$trait]
    drop((-crossprod(effect, provider$cross_product %*% effect) +
      2 * crossprod(effect, provider$score)) /
      (2 * provider$residual_scale))
  }, numeric(1)))
}
