# Independent test-only reference for the standard SBayesRC probit-alpha
# conditional. These helpers do not call package implementation code.

.sbayesrc_alpha_reference_design <- function(annotation) {
  annotation <- as.matrix(annotation)
  stopifnot(
    is.numeric(annotation), nrow(annotation) > 0L,
    all(is.finite(annotation))
  )
  cbind(Intercept = 1, annotation)
}

.sbayesrc_alpha_exact_posterior <- function(
    latent, annotation, tau2, intercept_mean = 0,
    intercept_precision = 0) {
  design <- .sbayesrc_alpha_reference_design(annotation)
  latent <- as.numeric(latent)
  stopifnot(
    length(latent) == nrow(design), all(is.finite(latent)),
    length(tau2) == 1L, is.finite(tau2), tau2 > 0,
    length(intercept_mean) == 1L, is.finite(intercept_mean),
    length(intercept_precision) == 1L,
    is.finite(intercept_precision), intercept_precision >= 0
  )
  prior_precision <- c(intercept_precision, rep(tau2^-1, ncol(annotation)))
  prior_mean <- c(intercept_mean, rep(0, ncol(annotation)))
  precision <- crossprod(design) + diag(prior_precision)
  rhs <- crossprod(design, latent) + prior_precision * prior_mean
  covariance <- solve(precision)
  mean <- drop(covariance %*% rhs)
  names(mean) <- colnames(design)
  dimnames(covariance) <- list(colnames(design), colnames(design))
  list(mean = mean, covariance = covariance, precision = precision)
}

.sbayesrc_alpha_scalar_moments <- function(
    latent, annotation, alpha, tau2, intercept_mean = 0,
    intercept_precision = 0) {
  design <- .sbayesrc_alpha_reference_design(annotation)
  latent <- as.numeric(latent)
  alpha <- as.numeric(alpha)
  stopifnot(length(latent) == nrow(design), length(alpha) == ncol(design))
  residual <- latent - drop(design %*% alpha)
  out <- matrix(NA_real_, ncol(design), 2L,
                dimnames = list(colnames(design), c("mean", "variance")))
  for (column in seq_len(ncol(design))) {
    x <- design[, column]
    diagonal <- sum(x * x)
    likelihood_rhs <- sum(x * residual) + diagonal * alpha[column]
    prior_precision <- if (column == 1L) intercept_precision else tau2^-1
    prior_mean <- if (column == 1L) intercept_mean else 0
    variance <- 1 / (diagonal + prior_precision)
    out[column, ] <- c(
      variance * (likelihood_rhs + prior_precision * prior_mean),
      variance
    )
  }
  out
}

.sbayesrc_alpha_scalar_sweep <- function(
    latent, annotation, alpha, tau2, intercept_mean = 0,
    intercept_precision = 0) {
  design <- .sbayesrc_alpha_reference_design(annotation)
  latent <- as.numeric(latent)
  alpha <- as.numeric(alpha)
  stopifnot(length(latent) == nrow(design), length(alpha) == ncol(design))
  residual <- latent - drop(design %*% alpha)
  for (column in seq_len(ncol(design))) {
    x <- design[, column]
    old <- alpha[column]
    diagonal <- sum(x * x)
    likelihood_rhs <- sum(x * residual) + diagonal * old
    prior_precision <- if (column == 1L) intercept_precision else tau2^-1
    prior_mean <- if (column == 1L) intercept_mean else 0
    variance <- 1 / (diagonal + prior_precision)
    mean <- variance *
      (likelihood_rhs + prior_precision * prior_mean)
    alpha[column] <- stats::rnorm(1L, mean, sqrt(variance))
    residual <- residual + x * (old - alpha[column])
  }
  alpha
}

.sbayesrc_alpha_scalar_chain <- function(
    latent, annotation, tau2, iterations, burn, initial = NULL,
    intercept_mean = 0, intercept_precision = 0) {
  design <- .sbayesrc_alpha_reference_design(annotation)
  stopifnot(iterations > burn, burn >= 0)
  alpha <- if (is.null(initial)) numeric(ncol(design)) else as.numeric(initial)
  stopifnot(length(alpha) == ncol(design))
  draws <- matrix(NA_real_, iterations - burn, ncol(design),
                  dimnames = list(NULL, colnames(design)))
  retained <- 0L
  for (iteration in seq_len(iterations)) {
    alpha <- .sbayesrc_alpha_scalar_sweep(
      latent, annotation, alpha, tau2, intercept_mean, intercept_precision
    )
    if (iteration > burn) {
      retained <- retained + 1L
      draws[retained, ] <- alpha
    }
  }
  draws
}

.sbayesrc_alpha_blocked_draws <- function(
    latent, annotation, tau2, draws, intercept_mean = 0,
    intercept_precision = 0) {
  posterior <- .sbayesrc_alpha_exact_posterior(
    latent, annotation, tau2, intercept_mean, intercept_precision
  )
  chol_precision <- chol(posterior$precision)
  standard_normal <- matrix(
    stats::rnorm(length(posterior$mean) * draws),
    nrow = length(posterior$mean)
  )
  sampled <- backsolve(chol_precision, standard_normal)
  sampled <- sweep(sampled, 1L, posterior$mean, `+`)
  result <- t(sampled)
  colnames(result) <- names(posterior$mean)
  result
}

.sbayesrc_alpha_reference_component_probability <- function(stick) {
  stick <- as.matrix(stick)
  stopifnot(all(is.finite(stick)), all(stick >= 0), all(stick <= 1))
  probability <- matrix(0, nrow(stick), ncol(stick) + 1L)
  remaining <- rep(1, nrow(stick))
  for (column in seq_len(ncol(stick))) {
    probability[, column] <- remaining * (1 - stick[, column])
    remaining <- remaining * stick[, column]
  }
  probability[, ncol(probability)] <- remaining
  probability
}

