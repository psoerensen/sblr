# Development/reference implementation used for posterior-correctness and
# sampler-feasibility audits. Not a production sampler option and not part of
# the supported public API. These helpers provide an independently testable
# implementation of the exact finite-subset formulas documented in docs/dev/.

.bayesrc_coordinated_component_probabilities <- function(annotation, alpha,
                                                           probability_floor = 1e-12) {
  annotation <- as.matrix(annotation)
  alpha <- as.matrix(alpha)
  if (nrow(annotation) < 1L || ncol(annotation) != nrow(alpha) ||
      ncol(alpha) < 1L || anyNA(annotation) || anyNA(alpha) ||
      any(!is.finite(annotation)) || any(!is.finite(alpha)) ||
      length(probability_floor) != 1L || !is.finite(probability_floor) ||
      probability_floor <= 0 || probability_floor >= 1) {
    stop("Invalid coordinated BayesRC probability inputs.", call. = FALSE)
  }
  continuation <- stats::pnorm(annotation %*% alpha)
  marker_count <- nrow(annotation)
  stick_count <- ncol(alpha)
  probability <- matrix(0, marker_count, stick_count + 1L)
  remaining <- rep(1, marker_count)
  for (stick in seq_len(stick_count)) {
    probability[, stick] <- remaining * (1 - continuation[, stick])
    remaining <- remaining * continuation[, stick]
  }
  probability[, stick_count + 1L] <- remaining
  probability <- pmax(probability, probability_floor)
  probability / rowSums(probability)
}

.bayesrc_coordinated_log_alpha_prior <- function(alpha, intercept_mean,
                                                  intercept_sd,
                                                  sigma_sq_alpha) {
  alpha <- as.numeric(alpha)
  if (!length(alpha) || length(intercept_mean) != 1L ||
      length(intercept_sd) != 1L || length(sigma_sq_alpha) != 1L ||
      any(!is.finite(c(alpha, intercept_mean, intercept_sd, sigma_sq_alpha))) ||
      intercept_sd <= 0 || sigma_sq_alpha <= 0) {
    stop("Invalid coordinated BayesRC alpha-prior inputs.", call. = FALSE)
  }
  stats::dnorm(alpha[[1L]], intercept_mean, intercept_sd, log = TRUE) +
    if (length(alpha) > 1L) {
      sum(stats::dnorm(alpha[-1L], 0, sqrt(sigma_sq_alpha), log = TRUE))
    } else 0
}

.bayesrc_coordinated_subset_states <- function(score, operator, residual_variance,
                                                marker_variance, gamma,
                                                marker_probability,
                                                marker_scale = NULL) {
  score <- as.numeric(score)
  operator <- as.matrix(operator)
  marker_probability <- as.matrix(marker_probability)
  gamma <- as.numeric(gamma)
  subset_size <- length(score)
  component_count <- length(gamma)
  if (is.null(marker_scale)) marker_scale <- rep(1, subset_size)
  marker_scale <- as.numeric(marker_scale)
  if (!subset_size || component_count < 2L || gamma[[1L]] != 0 ||
      any(gamma[-1L] <= 0) || nrow(operator) != subset_size ||
      ncol(operator) != subset_size || nrow(marker_probability) != subset_size ||
      ncol(marker_probability) != component_count ||
      length(marker_scale) != subset_size || any(marker_probability <= 0) ||
      any(!is.finite(c(score, operator, marker_probability, marker_scale,
                       residual_variance, marker_variance, gamma))) ||
      residual_variance <= 0 || marker_variance <= 0 || any(marker_scale <= 0)) {
    stop("Invalid coordinated BayesRC subset inputs.", call. = FALSE)
  }
  configurations <- as.matrix(expand.grid(rep(list(0:(component_count - 1L)),
                                               subset_size)))
  storage.mode(configurations) <- "integer"
  log_weight <- numeric(nrow(configurations))
  means <- vector("list", nrow(configurations))
  covariance <- vector("list", nrow(configurations))
  for (state in seq_len(nrow(configurations))) {
    component <- configurations[state, ]
    active <- which(component > 0L)
    log_weight[[state]] <- sum(log(marker_probability[
      cbind(seq_len(subset_size), component + 1L)]))
    means[[state]] <- numeric(subset_size)
    covariance[[state]] <- matrix(0, subset_size, subset_size)
    if (length(active)) {
      prior_variance <- marker_variance * marker_scale[active] *
        gamma[component[active] + 1L]
      precision <- operator[active, active, drop = FALSE] / residual_variance +
        diag(1 / prior_variance, length(active))
      chol_precision <- chol(precision)
      log_determinant_precision <- 2 * sum(log(diag(chol_precision)))
      h <- score[active] / residual_variance
      mean_active <- backsolve(chol_precision,
        forwardsolve(t(chol_precision), h))
      covariance_active <- chol2inv(chol_precision)
      log_weight[[state]] <- log_weight[[state]] - 0.5 *
        (sum(log(prior_variance)) + log_determinant_precision) +
        0.5 * drop(crossprod(h, mean_active))
      means[[state]][active] <- mean_active
      covariance[[state]][active, active] <- covariance_active
    }
  }
  maximum <- max(log_weight)
  probability <- exp(log_weight - maximum)
  probability <- probability / sum(probability)
  list(component = configurations, log_weight = log_weight,
       probability = probability,
       log_normalizer = maximum + log(sum(exp(log_weight - maximum))),
       mean = means, covariance = covariance)
}

.bayesrc_coordinated_log_mh <- function(alpha_old, alpha_new, component_outside,
                                        annotation_outside, alpha_all_old,
                                        alpha_all_new, subset_old, subset_new,
                                        intercept_mean, intercept_sd,
                                        sigma_sq_alpha, log_q_reverse = 0,
                                        log_q_forward = 0,
                                        probability_floor = 1e-12) {
  component_outside <- as.integer(component_outside)
  old_probability <- .bayesrc_coordinated_component_probabilities(
    annotation_outside, alpha_all_old, probability_floor)
  new_probability <- .bayesrc_coordinated_component_probabilities(
    annotation_outside, alpha_all_new, probability_floor)
  if (length(component_outside) != nrow(old_probability) ||
      any(component_outside < 0L) ||
      any(component_outside >= ncol(old_probability))) {
    stop("Invalid outside allocation state.", call. = FALSE)
  }
  index <- cbind(seq_along(component_outside), component_outside + 1L)
  .bayesrc_coordinated_log_alpha_prior(
    alpha_new, intercept_mean, intercept_sd, sigma_sq_alpha) -
    .bayesrc_coordinated_log_alpha_prior(
      alpha_old, intercept_mean, intercept_sd, sigma_sq_alpha) +
    sum(log(new_probability[index]) - log(old_probability[index])) +
    subset_new$log_normalizer - subset_old$log_normalizer +
    log_q_reverse - log_q_forward
}

.bayesrc_coordinated_beta_fixed_support <- function(component, beta,
                                                     tolerance = 0) {
  component <- as.integer(component)
  beta <- as.numeric(beta)
  if (length(component) != length(beta) || any(component < 0L) ||
      length(tolerance) != 1L || !is.finite(tolerance) || tolerance < 0) {
    stop("Invalid beta-fixed support inputs.", call. = FALSE)
  }
  all((component == 0L & abs(beta) <= tolerance) |
      (component > 0L & abs(beta) > tolerance))
}
