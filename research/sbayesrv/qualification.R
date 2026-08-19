# Deterministic SBayesRV Gate 1 qualification.
# Source prototype.R first, then call run_sbayesrv_qualification().

.sbayesrv_gate_row <- function(fixture, gate, value, tolerance, pass) {
  data.frame(
    fixture = fixture,
    gate = gate,
    value = as.numeric(value),
    tolerance = as.numeric(tolerance),
    pass = isTRUE(pass),
    stringsAsFactors = FALSE
  )
}

.sbayesrv_max_abs <- function(x, y) max(abs(as.numeric(x) - as.numeric(y)))

run_sbayesrv_qualification <- function(stop_on_failure = TRUE) {
  marker_ids <- paste0("m", seq_len(10L))
  binary <- rep(c(0, 1), each = 5L)
  continuous <- c(-4, -2.5, -1.5, -0.5, 0, 0.7, 1.4, 2.2, 3.1, 5)
  proxy <- continuous + c(0.2, -0.1, 0.15, -0.2, 0.05,
                           0.1, -0.15, 0.2, -0.05, -0.2)
  raw <- list(
    binary = matrix(binary, ncol = 1L,
                    dimnames = list(marker_ids, "binary")),
    continuous = matrix(continuous, ncol = 1L,
                        dimnames = list(marker_ids, "continuous")),
    correlated = cbind(binary = binary, continuous = continuous, proxy = proxy)
  )
  rownames(raw$correlated) <- marker_ids
  processed <- list(
    binary = sbayesrv_preprocess_annotations(
      raw$binary, marker_ids, c(binary = "binary")),
    continuous = sbayesrv_preprocess_annotations(
      raw$continuous, marker_ids, c(continuous = "continuous")),
    correlated = sbayesrv_preprocess_annotations(
      raw$correlated, marker_ids,
      c(binary = "binary", continuous = "continuous", proxy = "continuous"))
  )

  score <- setNames(
    c(0.15, -0.35, 0.8, 1.15, -0.7, 0.25, 1.4, -1.1, 0.55, -0.2),
    marker_ids)
  precision <- setNames(c(35, 42, 55, 63, 48, 39, 70, 61, 46, 52), marker_ids)
  gamma <- c(null = 0, small = 0.01, medium = 0.1, large = 1)
  pi_component <- c(0.72, 0.15, 0.09, 0.04)
  v_b <- 0.06
  sigma_theta <- 0.7
  beta <- setNames(
    c(0, 0.025, -0.08, 0.21, 0, -0.045, 0.16, 0, -0.12, 0.035),
    marker_ids)
  component <- setNames(c(0L, 1L, 2L, 3L, 0L, 1L, 3L, 0L, 2L, 1L), marker_ids)

  fixtures <- list(
    null = list(X = processed$correlated$X,
                theta = c(binary = 0, continuous = 0, proxy = 0)),
    one_binary = list(X = processed$binary$X, theta = c(binary = 0.65)),
    one_continuous = list(
      X = processed$continuous$X, theta = c(continuous = -0.45)),
    correlated = list(
      X = processed$correlated$X,
      theta = c(binary = 0.35, continuous = -0.25, proxy = 0.18))
  )

  rows <- list()
  add <- function(fixture, gate, value, tolerance, pass = value <= tolerance) {
    rows[[length(rows) + 1L]] <<-
      .sbayesrv_gate_row(fixture, gate, value, tolerance, pass)
  }

  for (fixture_name in names(processed)) {
    column_error <- max(abs(colMeans(processed[[fixture_name]]$X)))
    add(fixture_name, "processed_column_mean", column_error, 2e-15)
  }

  for (fixture_name in names(fixtures)) {
    fixture <- fixtures[[fixture_name]]
    scale <- sbayesrv_eta_q(fixture$theta, fixture$X)
    add(fixture_name, "eta_dimension", length(scale$eta) - nrow(fixture$X), 0,
        identical(names(scale$eta), rownames(fixture$X)))
    add(fixture_name, "q_dimension", length(scale$q) - nrow(fixture$X), 0,
        identical(names(scale$q), rownames(fixture$X)))
    add(fixture_name, "geometric_mean_q", abs(mean(log(scale$q))), 2e-14)

    conditional_fn <- function(value) {
      sbayesrv_conditional_theta_log_posterior(
        value, fixture$X, beta, component, v_b, gamma, sigma_theta)
    }
    conditional_analytic <- sbayesrv_conditional_theta_gradient(
      fixture$theta, fixture$X, beta, component, v_b, gamma, sigma_theta)
    # At these O(1)-to-O(10) log-posterior scales, a relative 1e-6 central
    # step balances second-order truncation against double-precision
    # cancellation. The 2e-7 gate is over 100 times the observed interior
    # discrepancy while remaining materially stricter than the model signal.
    conditional_numeric <- sbayesrv_central_gradient(
      conditional_fn, fixture$theta, relative_step = 1e-6)
    add(
      fixture_name, "conditional_gradient",
      .sbayesrv_max_abs(conditional_analytic, conditional_numeric), 2e-7)

    collapsed_fn <- function(value) {
      sbayesrv_collapsed_theta_log_posterior(
        value, fixture$X, score, precision, v_b, gamma, pi_component,
        sigma_theta)
    }
    collapsed_analytic <- sbayesrv_collapsed_theta_gradient(
      fixture$theta, fixture$X, score, precision, v_b, gamma, pi_component,
      sigma_theta)
    collapsed_numeric <- sbayesrv_central_gradient(
      collapsed_fn, fixture$theta, relative_step = 1e-6)
    add(
      fixture_name, "collapsed_gradient",
      .sbayesrv_max_abs(collapsed_analytic, collapsed_numeric), 2e-7)

    terms <- sbayesrv_independent_summary_terms(
      fixture$theta, fixture$X, score, precision, v_b, gamma, pi_component)
    add(fixture_name, "component_probability_rows",
        max(abs(rowSums(terms$component_probability) - 1)), 2e-15)
    add(fixture_name, "responsibility_rows",
        max(abs(rowSums(terms$responsibilities) - 1)), 2e-15)
  }

  multi <- sbayesrv_independent_summary_terms(
    fixtures$correlated$theta, fixtures$correlated$X,
    score, precision, v_b, gamma, pi_component)
  marker <- 7L
  active_component <- 4L
  tau2_direct <- v_b * gamma[active_component] * multi$q[marker]
  variance_direct <- 1 / (precision[marker] + 1 / tau2_direct)
  mean_direct <- variance_direct * score[marker]
  log_marginal_direct <- -0.5 * log(1 + precision[marker] * tau2_direct) +
    0.5 * score[marker]^2 * tau2_direct /
      (1 + precision[marker] * tau2_direct)
  add("multiple_components", "active_conditional_variance",
      abs(multi$posterior_variance[marker, active_component] - variance_direct),
      2e-15)
  add("multiple_components", "active_conditional_mean",
      abs(multi$posterior_mean[marker, active_component] - mean_direct), 2e-15)
  add("multiple_components", "scalar_marginal_weight",
      abs(multi$component_log_marginal[marker, active_component] -
            log_marginal_direct), 2e-15)

  zero_terms <- sbayesrv_independent_summary_terms(
    fixtures$null$theta, fixtures$null$X,
    score, precision, v_b, gamma, pi_component)
  ordinary_terms <- sbayesrv_independent_summary_terms_q(
    score, precision, setNames(rep(1, length(score)), marker_ids),
    v_b, gamma, pi_component)
  add("theta_zero", "ordinary_sbayesr_log_weights",
      .sbayesrv_max_abs(zero_terms$component_log_weight,
                        ordinary_terms$component_log_weight), 0)
  add("theta_zero", "ordinary_sbayesr_responsibilities",
      .sbayesrv_max_abs(zero_terms$responsibilities,
                        ordinary_terms$responsibilities), 0)

  fixed_q <- setNames(exp(seq(-0.8, 0.8, length.out = length(marker_ids))),
                      marker_ids)
  fixed_terms <- sbayesrv_independent_summary_terms_q(
    score, precision, fixed_q, v_b, gamma, pi_component)
  expected_tau <- outer(fixed_q, gamma, function(q, multiplier) {
    v_b * q * multiplier
  })
  dimnames(expected_tau) <- dimnames(fixed_terms$tau2)
  add("fixed_q", "fixed_marker_multiplier",
      .sbayesrv_max_abs(fixed_terms$tau2, expected_tau), 2e-17)

  one_gamma <- c(null = 0, active = 1)
  one_pi <- c(0.8, 0.2)
  one_terms <- sbayesrv_independent_summary_terms_q(
    score, precision, fixed_q, v_b, one_gamma, one_pi)
  one_direct <- -0.5 * log1p(precision * v_b * fixed_q) +
    0.5 * score^2 * v_b * fixed_q / (1 + precision * v_b * fixed_q)
  add("one_positive_component", "bayesc_style_log_marginal",
      .sbayesrv_max_abs(one_terms$component_log_marginal[, 2L], one_direct),
      2e-15)

  no_active <- integer(length(marker_ids))
  no_active_theta <- fixtures$correlated$theta
  no_active_log <- sbayesrv_conditional_theta_log_posterior(
    no_active_theta, fixtures$correlated$X, numeric(length(marker_ids)),
    no_active, v_b, gamma, sigma_theta)
  no_active_gradient <- sbayesrv_conditional_theta_gradient(
    no_active_theta, fixtures$correlated$X, numeric(length(marker_ids)),
    no_active, v_b, gamma, sigma_theta)
  add("no_active", "conditional_prior_log_posterior",
      abs(no_active_log + 0.5 * sum(no_active_theta^2) / sigma_theta^2),
      2e-15)
  add("no_active", "conditional_prior_gradient",
      .sbayesrv_max_abs(no_active_gradient,
                        -no_active_theta / sigma_theta^2), 2e-15)

  marker_order <- c(10, 2, 7, 1, 5, 8, 3, 9, 4, 6)
  permuted <- sbayesrv_independent_summary_terms(
    fixtures$correlated$theta,
    fixtures$correlated$X[marker_order, , drop = FALSE],
    score[marker_order], precision[marker_order], v_b, gamma, pi_component)
  add("alignment", "marker_permutation_log_likelihood",
      abs(multi$log_likelihood - permuted$log_likelihood), 2e-14)
  add("alignment", "marker_permutation_responsibilities",
      .sbayesrv_max_abs(
        multi$responsibilities[marker_order, , drop = FALSE],
        permuted$responsibilities), 2e-15)

  annotation_order <- c("proxy", "binary", "continuous")
  X_permuted <- fixtures$correlated$X[, annotation_order, drop = FALSE]
  theta_permuted <- fixtures$correlated$theta[annotation_order]
  annotation_gradient <- sbayesrv_collapsed_theta_gradient(
    theta_permuted, X_permuted, score, precision, v_b, gamma,
    pi_component, sigma_theta)
  reference_gradient <- sbayesrv_collapsed_theta_gradient(
    fixtures$correlated$theta, fixtures$correlated$X, score, precision,
    v_b, gamma, pi_component, sigma_theta)
  add("alignment", "annotation_permutation_posterior",
      abs(sbayesrv_collapsed_theta_log_posterior(
        theta_permuted, X_permuted, score, precision, v_b, gamma,
        pi_component, sigma_theta) -
        sbayesrv_collapsed_theta_log_posterior(
          fixtures$correlated$theta, fixtures$correlated$X, score, precision,
          v_b, gamma, pi_component, sigma_theta)), 2e-14)
  add("alignment", "annotation_permutation_gradient",
      .sbayesrv_max_abs(annotation_gradient[colnames(fixtures$correlated$X)],
                        reference_gradient), 2e-14)

  strong_fixtures <- list(
    strong_positive = c(continuous = 6),
    strong_negative = c(continuous = -6)
  )
  for (fixture_name in names(strong_fixtures)) {
    theta <- strong_fixtures[[fixture_name]]
    terms <- sbayesrv_independent_summary_terms(
      theta, processed$continuous$X, score, precision,
      v_b, gamma, pi_component)
    analytic <- sbayesrv_collapsed_theta_gradient(
      theta, processed$continuous$X, score, precision,
      v_b, gamma, pi_component, sigma_theta)
    numeric <- sbayesrv_central_gradient(function(value) {
      sbayesrv_collapsed_theta_log_posterior(
        value, processed$continuous$X, score, precision,
        v_b, gamma, pi_component, sigma_theta)
    }, theta, relative_step = 1e-6)
    add(fixture_name, "finite_strong_eta", 0, 0,
        all(is.finite(c(terms$eta, terms$q, terms$component_log_weight,
                        terms$responsibilities))))
    add(fixture_name, "collapsed_gradient",
        .sbayesrv_max_abs(analytic, numeric), 2e-7)
  }

  conditional_likelihood <- function(value) {
    sbayesrv_conditional_theta_log_posterior(
      value, fixtures$correlated$X, beta, component, v_b, gamma,
      sigma_theta) + 0.5 * sum(value^2) / sigma_theta^2
  }
  set.seed(20260819)
  draw_one <- sbayesrv_elliptical_slice_theta(
    fixtures$correlated$theta, conditional_likelihood, sigma_theta)
  set.seed(20260819)
  draw_two <- sbayesrv_elliptical_slice_theta(
    fixtures$correlated$theta, conditional_likelihood, sigma_theta)
  add("elliptical_slice", "same_seed", as.numeric(!identical(draw_one, draw_two)),
      0, identical(draw_one, draw_two))

  set.seed(81)
  empty_one <- sbayesrv_conditional_theta_update(
    fixtures$correlated$theta, fixtures$correlated$X,
    numeric(length(marker_ids)), no_active, v_b, gamma, sigma_theta)
  set.seed(81)
  empty_two <- sbayesrv_conditional_theta_update(
    fixtures$correlated$theta, fixtures$correlated$X,
    numeric(length(marker_ids)), no_active, v_b, gamma, sigma_theta)
  add("no_active", "same_seed_prior_draw",
      as.numeric(!identical(empty_one, empty_two)), 0,
      identical(empty_one, empty_two))

  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  if (isTRUE(stop_on_failure) && any(!result$pass)) {
    failures <- result[!result$pass, c("fixture", "gate", "value", "tolerance")]
    stop(
      paste0("SBayesRV deterministic qualification failed:\n",
             paste(capture.output(print(failures, row.names = FALSE)),
                   collapse = "\n")),
      call. = FALSE
    )
  }
  result
}
