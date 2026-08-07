# Development references for the exact probit scale-sandwich derivation.
# Production samplers do not call these R helpers.

.sbayesrc_px_factor <- function(latent, design, prior_mean, prior_precision) {
  latent <- as.numeric(latent)
  design <- as.matrix(design)
  prior_mean <- as.numeric(prior_mean)
  prior_precision <- as.numeric(prior_precision)
  if (!length(latent) || nrow(design) != length(latent) ||
      ncol(design) != length(prior_mean) ||
      length(prior_precision) != length(prior_mean) ||
      anyNA(c(latent, design, prior_mean, prior_precision)) ||
      any(!is.finite(c(latent, design, prior_mean, prior_precision))) ||
      any(prior_precision <= 0)) {
    stop("Invalid SBayesRC PX factor inputs.", call. = FALSE)
  }
  precision <- crossprod(design) + diag(prior_precision)
  chol_precision <- chol(precision)
  xtz <- drop(crossprod(design, latent))
  prior_rhs <- prior_precision * prior_mean
  solve_precision <- function(rhs) {
    backsolve(chol_precision,
      forwardsolve(t(chol_precision), rhs))
  }
  solved_xtz <- solve_precision(xtz)
  list(
    precision = precision,
    chol_precision = chol_precision,
    xtz = xtz,
    prior_rhs = prior_rhs,
    a = sum(latent^2) - drop(crossprod(xtz, solved_xtz)),
    b = drop(crossprod(xtz, solve_precision(prior_rhs)))
  )
}

.sbayesrc_px_log_scale_ratio <- function(log_scale, latent, design,
                                          prior_mean, prior_precision) {
  if (length(log_scale) != 1L || !is.finite(log_scale)) {
    stop("log_scale must be one finite value.", call. = FALSE)
  }
  factor <- .sbayesrc_px_factor(
    latent, design, prior_mean, prior_precision)
  scale <- exp(log_scale)
  if (!is.finite(scale) || scale <= 0) {
    return(-Inf)
  }
  length(latent) * log_scale - 0.5 * (
    factor$a * (scale^2 - 1) - 2 * factor$b * (scale - 1))
}

.sbayesrc_px_log_latent_marginal <- function(latent, design, prior_mean,
                                              prior_precision) {
  factor <- .sbayesrc_px_factor(
    latent, design, prior_mean, prior_precision)
  rhs <- factor$xtz + factor$prior_rhs
  solved <- backsolve(factor$chol_precision,
    forwardsolve(t(factor$chol_precision), rhs))
  -0.5 * (sum(latent^2) + sum(prior_precision * prior_mean^2) -
    drop(crossprod(rhs, solved))) - sum(log(diag(factor$chol_precision)))
}

.sbayesrc_px_alpha_conditional <- function(latent, design, prior_mean,
                                            prior_precision) {
  factor <- .sbayesrc_px_factor(
    latent, design, prior_mean, prior_precision)
  rhs <- factor$xtz + factor$prior_rhs
  mean <- backsolve(factor$chol_precision,
    forwardsolve(t(factor$chol_precision), rhs))
  list(mean = drop(mean), covariance = chol2inv(factor$chol_precision))
}
