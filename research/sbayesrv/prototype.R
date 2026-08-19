# Independent SBayesRV Gate 1 reference calculations.
# This file is research code. It neither calls package internals nor writes output.

.sbayesrv_stop <- function(...) stop(..., call. = FALSE)

.sbayesrv_named_vector <- function(x, ids, name) {
  supplied_names <- names(x)
  x <- as.numeric(x)
  if (length(x) != length(ids) || any(!is.finite(x))) {
    .sbayesrv_stop(name, " must contain one finite value per declared ID.")
  }
  if (!is.null(supplied_names)) {
    if (anyDuplicated(supplied_names) || !setequal(supplied_names, ids)) {
      .sbayesrv_stop(name, " names must match the declared IDs exactly.")
    }
    x <- x[match(ids, supplied_names)]
  }
  names(x) <- ids
  x
}

.sbayesrv_align_theta <- function(theta, annotation_ids) {
  supplied_names <- names(theta)
  theta <- as.numeric(theta)
  if (length(theta) != length(annotation_ids) || any(!is.finite(theta))) {
    .sbayesrv_stop("theta must contain one finite value per annotation.")
  }
  if (!is.null(supplied_names)) {
    if (anyDuplicated(supplied_names) || !setequal(supplied_names, annotation_ids)) {
      .sbayesrv_stop("theta names must match the annotation IDs exactly.")
    }
    theta <- theta[match(annotation_ids, supplied_names)]
  }
  names(theta) <- annotation_ids
  theta
}

sbayesrv_preprocess_annotations <- function(
    annotations, marker_ids, annotation_types = NULL) {
  if (!is.matrix(annotations) && !is.data.frame(annotations)) {
    .sbayesrv_stop("annotations must be a numeric matrix or data frame.")
  }
  marker_ids <- as.character(marker_ids)
  if (!length(marker_ids) || anyNA(marker_ids) || any(!nzchar(marker_ids)) ||
      anyDuplicated(marker_ids)) {
    .sbayesrv_stop("marker_ids must be explicit, nonempty, and unique.")
  }
  raw <- as.matrix(annotations)
  storage.mode(raw) <- "double"
  if (nrow(raw) != length(marker_ids) || ncol(raw) < 1L ||
      any(!is.finite(raw))) {
    .sbayesrv_stop(
      "annotations must be finite with one row per marker and at least one column.")
  }
  annotation_ids <- colnames(raw)
  if (is.null(annotation_ids) || anyNA(annotation_ids) ||
      any(!nzchar(annotation_ids)) || anyDuplicated(annotation_ids)) {
    .sbayesrv_stop("annotation column IDs must be explicit, nonempty, and unique.")
  }
  raw_marker_ids <- rownames(raw)
  if (!is.null(raw_marker_ids) &&
      !identical(as.character(raw_marker_ids), marker_ids)) {
    .sbayesrv_stop("annotation row IDs must match marker_ids exactly.")
  }
  rownames(raw) <- marker_ids

  inferred_types <- vapply(seq_len(ncol(raw)), function(column) {
    values <- unique(raw[, column])
    if (length(values) < 2L) {
      .sbayesrv_stop("annotation columns must not be constant.")
    }
    if (all(values %in% c(0, 1))) "binary" else "continuous"
  }, character(1))
  if (is.null(annotation_types)) {
    annotation_types <- inferred_types
  } else {
    supplied_names <- names(annotation_types)
    annotation_types <- as.character(annotation_types)
    if (!is.null(supplied_names)) {
      if (anyDuplicated(supplied_names) ||
          !setequal(supplied_names, annotation_ids)) {
        .sbayesrv_stop("annotation_types names must match annotation IDs exactly.")
      }
      annotation_types <- annotation_types[match(annotation_ids, supplied_names)]
    }
    if (length(annotation_types) != ncol(raw) ||
        any(!annotation_types %in% c("binary", "continuous"))) {
      .sbayesrv_stop("annotation_types must declare binary or continuous per column.")
    }
    if (!identical(unname(annotation_types), unname(inferred_types))) {
      .sbayesrv_stop(
        "annotation_types must follow the binary 0/1 versus continuous convention.")
    }
  }
  names(annotation_types) <- annotation_ids

  processed <- matrix(
    0, nrow(raw), ncol(raw), dimnames = list(marker_ids, annotation_ids))
  centers <- scales <- numeric(ncol(raw))
  transform <- character(ncol(raw))
  for (column in seq_len(ncol(raw))) {
    values <- raw[, column]
    centers[column] <- mean(values)
    if (annotation_types[column] == "binary") {
      scales[column] <- 1
      transform[column] <- "center_only"
    } else {
      scales[column] <- stats::sd(values)
      if (!is.finite(scales[column]) || scales[column] <= 0) {
        .sbayesrv_stop("continuous annotation columns require positive finite SD.")
      }
      transform[column] <- "center_and_sample_sd"
    }
    processed[, column] <- (values - centers[column]) / scales[column]
  }
  if (qr(processed, tol = sqrt(.Machine$double.eps))$rank < ncol(processed)) {
    .sbayesrv_stop(
      "processed annotations are duplicate or rank deficient.")
  }
  if (any(!is.finite(processed))) {
    .sbayesrv_stop("processed annotations contain a non-finite value.")
  }

  list(
    X = processed,
    transform = data.frame(
      annotation_id = annotation_ids,
      type = unname(annotation_types),
      center = centers,
      scale = scales,
      transform = transform,
      stringsAsFactors = FALSE
    )
  )
}

sbayesrv_eta_q <- function(theta, X) {
  X <- as.matrix(X)
  if (!is.numeric(X) || nrow(X) < 1L || ncol(X) < 1L ||
      is.null(rownames(X)) || is.null(colnames(X)) || any(!is.finite(X))) {
    .sbayesrv_stop("X must be a finite named marker-by-annotation matrix.")
  }
  theta <- .sbayesrv_align_theta(theta, colnames(X))
  eta <- as.numeric(X %*% theta)
  names(eta) <- rownames(X)
  if (any(!is.finite(eta))) .sbayesrv_stop("eta contains a non-finite value.")
  q <- exp(eta)
  names(q) <- rownames(X)
  if (any(!is.finite(q)) || any(q <= 0)) {
    .sbayesrv_stop("exp(eta) is outside the positive finite range.")
  }
  list(eta = eta, q = q, geometric_mean_q = exp(mean(eta)))
}

.sbayesrv_model_inputs <- function(X, beta, component, v_b, gamma) {
  X <- as.matrix(X)
  marker_ids <- rownames(X)
  if (is.null(marker_ids) || is.null(colnames(X)) || any(!is.finite(X))) {
    .sbayesrv_stop("X must retain finite marker and annotation dimensions.")
  }
  beta <- .sbayesrv_named_vector(beta, marker_ids, "beta")
  component_names <- names(component)
  component <- as.numeric(component)
  if (length(component) != nrow(X) || any(!is.finite(component)) ||
      any(component != floor(component))) {
    .sbayesrv_stop("component must contain one integer state per marker.")
  }
  if (!is.null(component_names)) {
    if (anyDuplicated(component_names) || !setequal(component_names, marker_ids)) {
      .sbayesrv_stop("component names must match marker IDs exactly.")
    }
    component <- component[match(marker_ids, component_names)]
  }
  component <- as.integer(component)
  names(component) <- marker_ids
  gamma <- as.numeric(gamma)
  if (length(gamma) < 2L || gamma[1L] != 0 || any(!is.finite(gamma)) ||
      any(gamma[-1L] <= 0) || any(component < 0L) ||
      any(component >= length(gamma))) {
    .sbayesrv_stop(
      "gamma must start at zero and component states must index its positive entries.")
  }
  if (length(v_b) != 1L || !is.finite(v_b) || v_b <= 0) {
    .sbayesrv_stop("v_b must be a positive finite scalar.")
  }
  list(X = X, beta = beta, component = component, v_b = v_b, gamma = gamma)
}

.sbayesrv_scaled_square <- function(effect, log_variance) {
  out <- numeric(length(effect))
  nonzero <- effect != 0
  out[nonzero] <- exp(2 * log(abs(effect[nonzero])) - log_variance[nonzero])
  if (any(!is.finite(out))) {
    .sbayesrv_stop("conditional theta quadratic term is non-finite.")
  }
  out
}

sbayesrv_conditional_theta_log_posterior <- function(
    theta, X, beta, component, v_b, gamma, sigma_theta = 0.7) {
  input <- .sbayesrv_model_inputs(X, beta, component, v_b, gamma)
  if (length(sigma_theta) != 1L || !is.finite(sigma_theta) || sigma_theta <= 0) {
    .sbayesrv_stop("sigma_theta must be a positive finite scalar.")
  }
  theta <- .sbayesrv_align_theta(theta, colnames(input$X))
  eta <- sbayesrv_eta_q(theta, input$X)$eta
  active <- which(input$component > 0L)
  log_likelihood <- 0
  if (length(active)) {
    log_variance <- log(input$v_b) +
      log(input$gamma[input$component[active] + 1L]) + eta[active]
    ratio <- .sbayesrv_scaled_square(input$beta[active], log_variance)
    log_likelihood <- -0.5 * sum(eta[active] + ratio)
  }
  log_likelihood - 0.5 * sum(theta^2) / sigma_theta^2
}

sbayesrv_conditional_theta_gradient <- function(
    theta, X, beta, component, v_b, gamma, sigma_theta = 0.7) {
  input <- .sbayesrv_model_inputs(X, beta, component, v_b, gamma)
  if (length(sigma_theta) != 1L || !is.finite(sigma_theta) || sigma_theta <= 0) {
    .sbayesrv_stop("sigma_theta must be a positive finite scalar.")
  }
  theta <- .sbayesrv_align_theta(theta, colnames(input$X))
  eta <- sbayesrv_eta_q(theta, input$X)$eta
  active <- which(input$component > 0L)
  gradient <- -theta / sigma_theta^2
  if (length(active)) {
    log_variance <- log(input$v_b) +
      log(input$gamma[input$component[active] + 1L]) + eta[active]
    ratio <- .sbayesrv_scaled_square(input$beta[active], log_variance)
    gradient <- gradient + as.numeric(crossprod(
      input$X[active, , drop = FALSE], 0.5 * (ratio - 1)))
  }
  names(gradient) <- colnames(input$X)
  gradient
}

.sbayesrv_row_log_sum_exp <- function(value) {
  value <- as.matrix(value)
  maximum <- apply(value, 1L, max)
  if (any(!is.finite(maximum))) {
    .sbayesrv_stop("each marker must have at least one finite component weight.")
  }
  maximum + log(rowSums(exp(value - maximum)))
}

.sbayesrv_summary_inputs <- function(score, precision, q, v_b, gamma, pi) {
  marker_ids <- names(score)
  score <- as.numeric(score)
  if (is.null(marker_ids)) marker_ids <- paste0("marker_", seq_along(score))
  precision <- .sbayesrv_named_vector(precision, marker_ids, "precision")
  q <- .sbayesrv_named_vector(q, marker_ids, "q")
  names(score) <- marker_ids
  if (any(!is.finite(score)) || any(precision < 0) || any(q <= 0)) {
    .sbayesrv_stop("score must be finite, precision non-negative, and q positive.")
  }
  gamma_ids <- names(gamma)
  gamma <- as.numeric(gamma)
  pi <- as.numeric(pi)
  if (length(gamma) < 2L || length(pi) != length(gamma) || gamma[1L] != 0 ||
      any(!is.finite(gamma)) || any(gamma[-1L] <= 0) ||
      any(!is.finite(pi)) || any(pi < 0) || sum(pi) <= 0) {
    .sbayesrv_stop("gamma and pi must define matching null-first BayesR components.")
  }
  pi <- pi / sum(pi)
  if (length(v_b) != 1L || !is.finite(v_b) || v_b <= 0) {
    .sbayesrv_stop("v_b must be a positive finite scalar.")
  }
  component_ids <- gamma_ids
  if (is.null(component_ids)) component_ids <- paste0("component_", seq_along(gamma) - 1L)
  list(
    marker_ids = marker_ids, component_ids = component_ids, score = score,
    precision = precision, q = q, v_b = v_b, gamma = gamma, pi = pi
  )
}

sbayesrv_independent_summary_terms_q <- function(
    score, precision, q, v_b, gamma, pi) {
  input <- .sbayesrv_summary_inputs(score, precision, q, v_b, gamma, pi)
  M <- length(input$score)
  K <- length(input$gamma)
  dims <- list(input$marker_ids, input$component_ids)
  tau2 <- matrix(0, M, K, dimnames = dims)
  for (component in seq.int(2L, K)) {
    tau2[, component] <- input$v_b * input$gamma[component] * input$q
  }
  posterior_variance <- posterior_mean <- matrix(0, M, K, dimnames = dims)
  component_log_marginal <- matrix(0, M, K, dimnames = dims)
  eta_derivative <- matrix(0, M, K, dimnames = dims)
  for (component in seq.int(2L, K)) {
    active_variance <- tau2[, component]
    denominator <- 1 + input$precision * active_variance
    posterior_variance[, component] <- active_variance / denominator
    posterior_mean[, component] <- input$score * active_variance / denominator
    component_log_marginal[, component] <-
      -0.5 * log1p(input$precision * active_variance) +
      0.5 * input$score^2 * active_variance / denominator
    eta_derivative[, component] <- 0.5 * (
      -input$precision * active_variance / denominator +
        input$score^2 * active_variance / denominator^2)
  }
  component_probability <- matrix(
    rep(input$pi, each = M), M, K, dimnames = dims)
  component_log_weight <- component_log_marginal + log(component_probability)
  marker_log_likelihood <- .sbayesrv_row_log_sum_exp(component_log_weight)
  names(marker_log_likelihood) <- input$marker_ids
  responsibilities <- exp(component_log_weight - marker_log_likelihood)
  dimnames(responsibilities) <- dims
  if (any(!is.finite(component_log_marginal)) ||
      any(!is.finite(marker_log_likelihood)) ||
      any(!is.finite(responsibilities))) {
    .sbayesrv_stop("independent-summary calculations produced a non-finite value.")
  }
  list(
    q = input$q,
    tau2 = tau2,
    component_probability = component_probability,
    posterior_mean = posterior_mean,
    posterior_variance = posterior_variance,
    component_log_marginal = component_log_marginal,
    component_log_weight = component_log_weight,
    marker_log_likelihood = marker_log_likelihood,
    responsibilities = responsibilities,
    eta_derivative = eta_derivative,
    log_likelihood = sum(marker_log_likelihood)
  )
}

sbayesrv_independent_summary_terms <- function(
    theta, X, score, precision, v_b, gamma, pi) {
  X <- as.matrix(X)
  theta <- .sbayesrv_align_theta(theta, colnames(X))
  scale <- sbayesrv_eta_q(theta, X)
  terms <- sbayesrv_independent_summary_terms_q(
    score, precision, scale$q, v_b, gamma, pi)
  terms$eta <- scale$eta
  terms$geometric_mean_q <- scale$geometric_mean_q
  terms
}

sbayesrv_collapsed_theta_log_posterior <- function(
    theta, X, score, precision, v_b, gamma, pi, sigma_theta = 0.7) {
  if (length(sigma_theta) != 1L || !is.finite(sigma_theta) || sigma_theta <= 0) {
    .sbayesrv_stop("sigma_theta must be a positive finite scalar.")
  }
  theta <- .sbayesrv_align_theta(theta, colnames(X))
  terms <- sbayesrv_independent_summary_terms(
    theta, X, score, precision, v_b, gamma, pi)
  terms$log_likelihood - 0.5 * sum(theta^2) / sigma_theta^2
}

sbayesrv_collapsed_theta_gradient <- function(
    theta, X, score, precision, v_b, gamma, pi, sigma_theta = 0.7) {
  if (length(sigma_theta) != 1L || !is.finite(sigma_theta) || sigma_theta <= 0) {
    .sbayesrv_stop("sigma_theta must be a positive finite scalar.")
  }
  X <- as.matrix(X)
  theta <- .sbayesrv_align_theta(theta, colnames(X))
  terms <- sbayesrv_independent_summary_terms(
    theta, X, score, precision, v_b, gamma, pi)
  marker_gradient <- rowSums(terms$responsibilities * terms$eta_derivative)
  gradient <- as.numeric(crossprod(X, marker_gradient)) - theta / sigma_theta^2
  names(gradient) <- colnames(X)
  gradient
}

sbayesrv_central_gradient <- function(fn, theta, relative_step = 1e-6) {
  theta <- as.numeric(theta)
  if (!length(theta) || any(!is.finite(theta)) ||
      length(relative_step) != 1L || !is.finite(relative_step) ||
      relative_step <= 0) {
    .sbayesrv_stop("finite-difference inputs are invalid.")
  }
  out <- numeric(length(theta))
  for (column in seq_along(theta)) {
    step <- relative_step * max(1, abs(theta[column]))
    plus <- minus <- theta
    plus[column] <- plus[column] + step
    minus[column] <- minus[column] - step
    out[column] <- (fn(plus) - fn(minus)) / (2 * step)
  }
  out
}

sbayesrv_elliptical_slice_theta <- function(
    theta, log_likelihood, sigma_theta = 0.7, max_evaluations = 10000L) {
  theta_names <- names(theta)
  theta <- as.numeric(theta)
  if (!length(theta) || any(!is.finite(theta)) ||
      length(sigma_theta) != 1L || !is.finite(sigma_theta) || sigma_theta <= 0) {
    .sbayesrv_stop("elliptical-slice inputs require finite theta and positive prior SD.")
  }
  current <- log_likelihood(theta)
  if (!is.finite(current)) .sbayesrv_stop("current theta likelihood is non-finite.")
  direction <- stats::rnorm(length(theta), 0, sigma_theta)
  threshold <- current + log(max(stats::runif(1L), .Machine$double.xmin))
  angle <- stats::runif(1L, 0, 2 * pi)
  lower <- angle - 2 * pi
  upper <- angle
  for (evaluation in seq_len(as.integer(max_evaluations))) {
    proposal <- theta * cos(angle) + direction * sin(angle)
    value <- tryCatch(log_likelihood(proposal), error = function(error) -Inf)
    if (is.finite(value) && value > threshold) {
      names(proposal) <- theta_names
      return(proposal)
    }
    if (angle < 0) lower <- angle else upper <- angle
    angle <- stats::runif(1L, lower, upper)
  }
  .sbayesrv_stop("elliptical slice update exceeded its evaluation guard.")
}

sbayesrv_conditional_theta_update <- function(
    theta, X, beta, component, v_b, gamma, sigma_theta = 0.7,
    max_evaluations = 10000L) {
  input <- .sbayesrv_model_inputs(X, beta, component, v_b, gamma)
  theta <- .sbayesrv_align_theta(theta, colnames(input$X))
  if (length(sigma_theta) != 1L || !is.finite(sigma_theta) || sigma_theta <= 0) {
    .sbayesrv_stop("sigma_theta must be a positive finite scalar.")
  }
  if (!any(input$component > 0L)) {
    draw <- stats::rnorm(length(theta), 0, sigma_theta)
    names(draw) <- names(theta)
    return(draw)
  }
  log_likelihood <- function(value) {
    sbayesrv_conditional_theta_log_posterior(
      value, input$X, input$beta, input$component, input$v_b, input$gamma,
      sigma_theta
    ) + 0.5 * sum(value^2) / sigma_theta^2
  }
  sbayesrv_elliptical_slice_theta(
    theta, log_likelihood, sigma_theta, max_evaluations)
}
