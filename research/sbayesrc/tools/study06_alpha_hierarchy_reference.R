# Independent R reference calculations for the Study 06 alpha-hierarchy audit.
# This file is development tooling, not package API.

alpha_hierarchy_stick_state <- function(component, nstep) {
  component <- as.integer(component)
  stopifnot(all(component >= 0L), nstep >= 1L)
  indicator <- outer(component, seq_len(nstep) - 1L, `>`)
  eligible <- lapply(seq_len(nstep), function(j) {
    if (j == 1L) seq_along(component) else which(indicator[, j - 1L])
  })
  list(
    indicator = indicator,
    eligible = eligible,
    continuation = vapply(seq_len(nstep), function(j)
      sum(indicator[eligible[[j]], j]), integer(1)),
    stopping = vapply(seq_len(nstep), function(j)
      sum(!indicator[eligible[[j]], j]), integer(1)))
}

alpha_hierarchy_scalar_conditionals <- function(X, z, alpha, sigma_sq,
                                                intercept_mean,
                                                intercept_sd) {
  X <- as.matrix(X)
  z <- as.numeric(z)
  alpha <- as.numeric(alpha)
  stopifnot(nrow(X) == length(z), ncol(X) == length(alpha),
            length(sigma_sq) == 1L, sigma_sq > 0,
            length(intercept_mean) == 1L,
            length(intercept_sd) == 1L, intercept_sd > 0)
  residual <- z - drop(X %*% alpha)
  out <- matrix(NA_real_, ncol(X), 2L,
                dimnames = list(colnames(X), c("mean", "variance")))
  for (k in seq_len(ncol(X))) {
    x <- X[, k]
    diagonal <- sum(x * x)
    likelihood_rhs <- sum(x * residual) + diagonal * alpha[k]
    prior_precision <- if (k == 1L) intercept_sd^-2 else sigma_sq^-1
    prior_mean <- if (k == 1L) intercept_mean else 0
    variance <- 1 / (diagonal + prior_precision)
    out[k, ] <- c(
      variance * (likelihood_rhs + prior_precision * prior_mean), variance)
  }
  out
}

alpha_hierarchy_variance_conditional <- function(alpha_non_intercept, a, b) {
  q <- length(alpha_non_intercept)
  list(
    df = q + a,
    numerator = sum(alpha_non_intercept^2) + b,
    inverse_gamma_shape = (q + a) / 2,
    inverse_gamma_scale = (sum(alpha_non_intercept^2) + b) / 2)
}

alpha_hierarchy_prior_mapping <- function(a, b) {
  stopifnot(length(a) == 1L, length(b) == 1L, a > 0, b > 0)
  data.frame(
    a = a,
    b = b,
    nu0 = a,
    scale0 = b / a,
    inverse_gamma_shape = a / 2,
    inverse_gamma_scale = b / 2,
    prior_mean = if (a > 2) b / (a - 2) else Inf,
    prior_variance = if (a > 4)
      2 * b^2 / ((a - 2)^2 * (a - 4)) else Inf,
    survival_tail_power = a / 2,
    stringsAsFactors = FALSE)
}

alpha_hierarchy_component_probability <- function(stick_probability) {
  stick_probability <- as.matrix(stick_probability)
  n <- nrow(stick_probability)
  nstep <- ncol(stick_probability)
  out <- matrix(0, n, nstep + 1L)
  remaining <- rep(1, n)
  for (j in seq_len(nstep)) {
    out[, j] <- remaining * (1 - stick_probability[, j])
    remaining <- remaining * stick_probability[, j]
  }
  out[, nstep + 1L] <- remaining
  out
}

alpha_hierarchy_convergence <- function(draws) {
  if (!requireNamespace("posterior", quietly = TRUE))
    stop("posterior is required for convergence summaries.")
  x <- posterior::as_draws_array(draws)
  summary <- posterior::summarise_draws(
    x, "rhat", "ess_bulk", "ess_tail", "mcse_mean", "sd")
  summary$relative_mcse <- summary$mcse_mean / summary$sd
  as.data.frame(summary)
}
