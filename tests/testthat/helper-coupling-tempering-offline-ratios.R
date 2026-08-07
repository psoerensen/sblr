# Independent test-only oracle for the historical coupling-tempering ratios.
# Package tests must not source development scripts from tools/ because tools/
# is intentionally excluded from installed source packages.

.offline_tempered_component_probability <- function(
    annotation, alpha, baseline, lambda, floor = 1e-12) {
  eta <- sweep(lambda * annotation %*% alpha, 2L,
               (1 - lambda) * baseline, `+`)
  continuation <- pnorm(eta)
  probability <- matrix(0, nrow(annotation), ncol(alpha) + 1L)
  remaining <- rep(1, nrow(annotation))
  for (stick in seq_len(ncol(alpha))) {
    probability[, stick] <- remaining * (1 - continuation[, stick])
    remaining <- remaining * continuation[, stick]
  }
  probability[, ncol(alpha) + 1L] <- remaining
  probability <- pmax(probability, floor)
  probability / rowSums(probability)
}

.offline_log_allocation_prior <- function(
    annotation, alpha, baseline, lambda, component, floor = 1e-12) {
  probability <- .offline_tempered_component_probability(
    annotation, alpha, baseline, lambda, floor)
  sum(log(probability[cbind(seq_len(nrow(probability)), component + 1L)]))
}

.offline_log_nonintercept_prior <- function(alpha, sigma_sq) {
  stopifnot(nrow(alpha) >= 1L, ncol(alpha) == length(sigma_sq))
  if (nrow(alpha) == 1L) return(0)
  sum(vapply(seq_len(ncol(alpha)), function(stick) {
    sum(dnorm(alpha[-1L, stick], 0, sqrt(sigma_sq[[stick]]), log = TRUE))
  }, numeric(1L)))
}

.offline_complete_exchange_ratio <- function(
    annotation, baseline, lambda_a, lambda_b,
    component_a, alpha_a, component_b, alpha_b, floor = 1e-12) {
  .offline_log_allocation_prior(
    annotation, alpha_b, baseline, lambda_a, component_b, floor) +
    .offline_log_allocation_prior(
      annotation, alpha_a, baseline, lambda_b, component_a, floor) -
    .offline_log_allocation_prior(
      annotation, alpha_a, baseline, lambda_a, component_a, floor) -
    .offline_log_allocation_prior(
      annotation, alpha_b, baseline, lambda_b, component_b, floor)
}

.offline_alpha_sigma_exchange_ratio <- function(
    annotation, baseline, lambda_a, lambda_b,
    component_a, alpha_a, sigma_a,
    component_b, alpha_b, sigma_b,
    prior_df, prior_scale, floor = 1e-12) {
  # Exchanging alpha and sigma together merely permutes complete hierarchy
  # prior factors, so those terms cancel exactly.
  .offline_log_allocation_prior(
    annotation, alpha_b, baseline, lambda_a, component_a, floor) +
    .offline_log_allocation_prior(
      annotation, alpha_a, baseline, lambda_b, component_b, floor) -
    .offline_log_allocation_prior(
      annotation, alpha_a, baseline, lambda_a, component_a, floor) -
    .offline_log_allocation_prior(
      annotation, alpha_b, baseline, lambda_b, component_b, floor)
}

.offline_alpha_only_exchange_ratio <- function(
    annotation, baseline, lambda_a, lambda_b,
    component_a, alpha_a, sigma_a,
    component_b, alpha_b, sigma_b, floor = 1e-12) {
  allocation <- .offline_log_allocation_prior(
    annotation, alpha_b, baseline, lambda_a, component_a, floor) +
    .offline_log_allocation_prior(
      annotation, alpha_a, baseline, lambda_b, component_b, floor) -
    .offline_log_allocation_prior(
      annotation, alpha_a, baseline, lambda_a, component_a, floor) -
    .offline_log_allocation_prior(
      annotation, alpha_b, baseline, lambda_b, component_b, floor)
  hierarchy <- .offline_log_nonintercept_prior(alpha_b, sigma_a) +
    .offline_log_nonintercept_prior(alpha_a, sigma_b) -
    .offline_log_nonintercept_prior(alpha_a, sigma_a) -
    .offline_log_nonintercept_prior(alpha_b, sigma_b)
  allocation + hierarchy
}

