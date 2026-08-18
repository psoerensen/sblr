# Independent R reference for the continuous-alpha collapsed posterior.
# Research code only: no package API, variance learning, selection, or LD model.

normal_log_density <- function(x, mean = 0, variance) {
  if (any(!is.finite(x)) || any(!is.finite(variance)) || any(variance <= 0)) {
    stop("Normal inputs and variances must be finite, with positive variances.")
  }
  -0.5 * (log(2 * pi * variance) + (x - mean)^2 / variance)
}

row_log_sum_exp <- function(x) {
  x <- as.matrix(x)
  largest <- apply(x, 1L, max)
  largest + log(rowSums(exp(x - largest)))
}

make_collapsed_alpha_model <- function(
    summary_beta, summary_se, V_beta, mixture_scale, annotation,
    continuation_center, prior_mean, prior_sd,
    marker_id = names(summary_beta), annotation_id = colnames(annotation),
    stick_id = NULL, component_id = NULL) {
  annotation <- as.matrix(annotation)
  storage.mode(annotation) <- "double"
  summary_beta <- as.numeric(summary_beta)
  summary_se <- as.numeric(summary_se)
  mixture_scale <- as.numeric(mixture_scale)
  n_marker <- length(summary_beta)
  n_component <- length(mixture_scale)
  n_stick <- n_component - 1L
  n_annotation <- ncol(annotation)
  continuation_center <- as.matrix(continuation_center)
  prior_mean <- as.matrix(prior_mean)
  prior_sd <- as.matrix(prior_sd)

  if (n_marker < 1L || length(summary_se) != n_marker ||
      nrow(annotation) != n_marker || n_component < 2L) {
    stop("Summary, annotation, and mixture dimensions are inconsistent.")
  }
  if (any(!is.finite(c(summary_beta, summary_se, annotation))) ||
      any(summary_se <= 0) || length(V_beta) != 1L ||
      !is.finite(V_beta) || V_beta <= 0) {
    stop("Summary inputs must be finite; SE and V_beta must be positive.")
  }
  if (mixture_scale[1L] != 0 || any(!is.finite(mixture_scale)) ||
      any(mixture_scale[-1L] <= 0)) {
    stop("The null mixture scale must be zero and active scales positive.")
  }
  target_dimension <- c(n_stick, n_annotation + 1L)
  if (!identical(dim(continuation_center), c(n_stick, n_annotation)) ||
      !identical(dim(prior_mean), target_dimension) ||
      !identical(dim(prior_sd), target_dimension) ||
      any(!is.finite(c(continuation_center, prior_mean, prior_sd))) ||
      any(prior_sd <= 0)) {
    stop("Center and prior dimensions must be stick by annotation/coefficient.")
  }

  marker_id <- marker_id %||% paste0("marker", seq_len(n_marker))
  annotation_id <- annotation_id %||% paste0("annotation", seq_len(n_annotation))
  stick_id <- stick_id %||% paste0("stick", seq_len(n_stick))
  component_id <- component_id %||% paste0("component", 0:(n_component - 1L))
  if (length(marker_id) != n_marker || length(annotation_id) != n_annotation ||
      length(stick_id) != n_stick || length(component_id) != n_component ||
      anyDuplicated(marker_id) || anyDuplicated(annotation_id) ||
      anyDuplicated(stick_id) || anyDuplicated(component_id)) {
    stop("Scientific IDs must be unique and match their dimensions.")
  }

  design <- lapply(seq_len(n_stick), function(stick) {
    centered <- if (n_annotation) {
      sweep(annotation, 2L, continuation_center[stick, ], "-")
    } else {
      matrix(numeric(), n_marker, 0L)
    }
    cbind(intercept = 1, centered)
  })
  coefficient_id <- c("intercept", annotation_id)
  design_variance <- outer(summary_se^2, rep(1, n_component)) +
    outer(rep(V_beta, n_marker), mixture_scale)
  component_log_density <- -0.5 * (
    log(2 * pi * design_variance) + summary_beta^2 / design_variance
  )
  dimnames(component_log_density) <- list(marker_id, component_id)
  dimnames(prior_mean) <- dimnames(prior_sd) <- list(stick_id, coefficient_id)
  dimnames(continuation_center) <- list(stick_id, annotation_id)

  structure(list(
    summary_beta = summary_beta,
    summary_se = summary_se,
    V_beta = V_beta,
    mixture_scale = mixture_scale,
    annotation = annotation,
    continuation_center = continuation_center,
    prior_mean = prior_mean,
    prior_sd = prior_sd,
    design = design,
    component_log_density = component_log_density,
    marker_id = marker_id,
    annotation_id = annotation_id,
    coefficient_id = coefficient_id,
    stick_id = stick_id,
    component_id = component_id,
    n_marker = n_marker,
    n_annotation = n_annotation,
    n_stick = n_stick,
    n_component = n_component,
    parameter_dimension = prod(target_dimension)
  ), class = "collapsed_alpha_model")
}

`%||%` <- function(x, y) if (is.null(x)) y else x

as_theta_matrix <- function(theta, model) {
  theta <- if (is.matrix(theta)) as.numeric(t(theta)) else as.numeric(theta)
  if (length(theta) != model$parameter_dimension || any(!is.finite(theta))) {
    stop("theta must be a finite stick-by-coefficient parameter vector.")
  }
  out <- matrix(theta, model$n_stick, model$n_annotation + 1L, byrow = TRUE)
  dimnames(out) <- dimnames(model$prior_mean)
  out
}

flatten_theta <- function(theta) as.numeric(t(as.matrix(theta)))

theta_to_alpha <- function(theta, model) {
  theta <- as_theta_matrix(theta, model)
  alpha <- t(theta)
  if (model$n_annotation) {
    for (stick in seq_len(model$n_stick)) {
      alpha[1L, stick] <- theta[stick, 1L] -
        sum(model$continuation_center[stick, ] * theta[stick, -1L])
    }
  }
  dimnames(alpha) <- list(model$coefficient_id, model$stick_id)
  alpha
}

alpha_to_theta <- function(alpha, model) {
  alpha <- as.matrix(alpha)
  if (!identical(dim(alpha), c(model$n_annotation + 1L, model$n_stick)) ||
      any(!is.finite(alpha))) {
    stop("alpha must be a finite coefficient-by-stick matrix.")
  }
  theta <- t(alpha)
  if (model$n_annotation) {
    for (stick in seq_len(model$n_stick)) {
      theta[stick, 1L] <- alpha[1L, stick] +
        sum(model$continuation_center[stick, ] * alpha[-1L, stick])
    }
  }
  dimnames(theta) <- dimnames(model$prior_mean)
  theta
}

continuation_probabilities <- function(theta, model) {
  theta <- as_theta_matrix(theta, model)
  eta <- matrix(
    NA_real_, model$n_marker, model$n_stick,
    dimnames = list(model$marker_id, model$stick_id)
  )
  for (stick in seq_len(model$n_stick)) {
    eta[, stick] <- as.numeric(model$design[[stick]] %*% theta[stick, ])
  }
  if (any(!is.finite(eta))) {
    stop("The continuation linear predictors must be finite.")
  }
  probability <- matrix(
    pnorm(eta), model$n_marker, model$n_stick,
    dimnames = dimnames(eta)
  )
  log_probability <- matrix(
    pnorm(eta, log.p = TRUE), model$n_marker, model$n_stick,
    dimnames = dimnames(eta)
  )
  log_survival_probability <- matrix(
    pnorm(eta, lower.tail = FALSE, log.p = TRUE),
    model$n_marker, model$n_stick, dimnames = dimnames(eta)
  )
  list(
    probability = probability,
    log_probability = log_probability,
    log_survival_probability = log_survival_probability,
    linear_predictor = eta
  )
}

stick_breaking_log_probabilities <- function(log_q, log_1mq,
                                             marker_id = NULL,
                                             component_id = NULL) {
  log_q <- as.matrix(log_q)
  log_1mq <- as.matrix(log_1mq)
  if (!identical(dim(log_q), dim(log_1mq)) || nrow(log_q) < 1L ||
      ncol(log_q) < 1L || any(is.na(log_q)) || any(is.na(log_1mq)) ||
      any(log_q > 0) || any(log_1mq > 0)) {
    stop("Log continuation probabilities must have valid matching dimensions.")
  }
  n_marker <- nrow(log_q)
  n_stick <- ncol(log_q)
  marker_id <- marker_id %||% rownames(log_q)
  component_id <- component_id %||% paste0("component", 0:n_stick)
  log_probability <- matrix(
    NA_real_, n_marker, n_stick + 1L,
    dimnames = list(marker_id, component_id)
  )
  log_remaining <- rep(0, n_marker)
  for (stick in seq_len(n_stick)) {
    log_probability[, stick] <- log_remaining + log_1mq[, stick]
    log_remaining <- log_remaining + log_q[, stick]
  }
  log_probability[, n_stick + 1L] <- log_remaining
  probability <- matrix(
    exp(log_probability), n_marker, n_stick + 1L,
    dimnames = dimnames(log_probability)
  )
  probability <- probability / rowSums(probability)
  list(probability = probability, log_probability = log_probability)
}

stick_breaking_probabilities <- function(continuation) {
  continuation <- as.matrix(continuation)
  if (nrow(continuation) < 1L || ncol(continuation) < 1L ||
      any(!is.finite(continuation)) ||
      any(continuation < 0 | continuation > 1)) {
    stop("Continuation probabilities must be a marker-by-stick matrix in [0, 1].")
  }
  stick_breaking_log_probabilities(
    log(continuation), log1p(-continuation), rownames(continuation)
  )
}

collapsed_marker_terms <- function(theta, model) {
  continuation <- continuation_probabilities(theta, model)
  mixture <- stick_breaking_log_probabilities(
    continuation$log_probability,
    continuation$log_survival_probability,
    model$marker_id,
    model$component_id
  )
  log_weight <- mixture$log_probability + model$component_log_density
  marker_log_likelihood <- row_log_sum_exp(log_weight)
  responsibility <- exp(log_weight - marker_log_likelihood)
  dimnames(log_weight) <- dimnames(responsibility) <-
    list(model$marker_id, model$component_id)
  names(marker_log_likelihood) <- model$marker_id
  list(
    continuation_probability = continuation$probability,
    log_continuation_probability = continuation$log_probability,
    log_continuation_survival_probability =
      continuation$log_survival_probability,
    linear_predictor = continuation$linear_predictor,
    component_probability = mixture$probability,
    log_component_probability = mixture$log_probability,
    component_log_density = model$component_log_density,
    component_log_weight = log_weight,
    marker_log_likelihood = marker_log_likelihood,
    responsibility = responsibility,
    log_likelihood = sum(marker_log_likelihood)
  )
}

alpha_log_prior <- function(theta, model, include_constant = TRUE) {
  theta <- as_theta_matrix(theta, model)
  standardized <- (theta - model$prior_mean) / model$prior_sd
  kernel <- -0.5 * sum(standardized^2)
  if (!include_constant) return(kernel)
  kernel - sum(log(model$prior_sd)) -
    0.5 * model$parameter_dimension * log(2 * pi)
}

collapsed_log_posterior_details <- function(theta, model) {
  marker <- collapsed_marker_terms(theta, model)
  prior_kernel <- alpha_log_prior(theta, model, include_constant = FALSE)
  prior <- alpha_log_prior(theta, model, include_constant = TRUE)
  list(
    marker = marker,
    log_likelihood = marker$log_likelihood,
    log_prior = prior,
    log_prior_kernel = prior_kernel,
    log_posterior = marker$log_likelihood + prior,
    log_posterior_kernel = marker$log_likelihood + prior_kernel
  )
}

collapsed_log_posterior_and_gradient <- function(theta, model,
                                                 include_constant = TRUE) {
  theta <- as_theta_matrix(theta, model)
  details <- collapsed_log_posterior_details(theta, model)
  r <- details$marker$responsibility
  log_q <- details$marker$log_continuation_probability
  log_1mq <- details$marker$log_continuation_survival_probability
  eta <- details$marker$linear_predictor
  gradient <- matrix(0, model$n_stick, model$n_annotation + 1L,
                     dimnames = dimnames(model$prior_mean))
  for (stick in seq_len(model$n_stick)) {
    success <- rowSums(r[, (stick + 1L):model$n_component, drop = FALSE])
    failure <- r[, stick]
    log_phi <- dnorm(eta[, stick], log = TRUE)
    score <- success * exp(log_phi - log_q[, stick]) -
      failure * exp(log_phi - log_1mq[, stick])
    gradient[stick, ] <- crossprod(model$design[[stick]], score)
  }
  gradient <- gradient - (theta - model$prior_mean) / model$prior_sd^2
  value <- if (include_constant) details$log_posterior else
    details$log_posterior_kernel
  list(value = value, gradient = flatten_theta(gradient), details = details)
}

collapsed_log_posterior <- function(theta, model, include_constant = TRUE) {
  collapsed_log_posterior_and_gradient(theta, model, include_constant)$value
}

check_collapsed_gradient <- function(theta, model, step = 5e-5) {
  theta <- if (is.matrix(theta)) flatten_theta(theta) else as.numeric(theta)
  analytic <- collapsed_log_posterior_and_gradient(theta, model)$gradient
  finite_difference <- vapply(seq_along(theta), function(index) {
    plus_two <- plus_one <- minus_one <- minus_two <- theta
    plus_two[index] <- plus_two[index] + 2 * step
    plus_one[index] <- plus_one[index] + step
    minus_one[index] <- minus_one[index] - step
    minus_two[index] <- minus_two[index] - 2 * step
    (-collapsed_log_posterior(plus_two, model) +
       8 * collapsed_log_posterior(plus_one, model) -
       8 * collapsed_log_posterior(minus_one, model) +
       collapsed_log_posterior(minus_two, model)) / (12 * step)
  }, numeric(1))
  list(
    analytic = analytic,
    finite_difference = finite_difference,
    discrepancy = analytic - finite_difference,
    maximum_absolute_discrepancy = max(abs(analytic - finite_difference))
  )
}

regularize_laplace_metric <- function(metric, minimum_eigenvalue = 1e-6,
                                      maximum_condition = 1e8) {
  metric <- 0.5 * (as.matrix(metric) + t(as.matrix(metric)))
  if (any(!is.finite(metric)) || nrow(metric) != ncol(metric)) {
    stop("Laplace precision metric must be finite and square.")
  }
  decomposition <- eigen(metric, symmetric = TRUE)
  maximum <- max(decomposition$values)
  if (!is.finite(maximum) || maximum <= 0) {
    stop("The observed posterior precision has no positive eigenvalue.")
  }
  floor_value <- max(minimum_eigenvalue, maximum / maximum_condition)
  regularized_value <- pmax(decomposition$values, floor_value)
  precision <- decomposition$vectors %*%
    (regularized_value * t(decomposition$vectors))
  covariance <- decomposition$vectors %*%
    ((1 / regularized_value) * t(decomposition$vectors))
  precision <- 0.5 * (precision + t(precision))
  covariance <- 0.5 * (covariance + t(covariance))
  root <- tryCatch(t(chol(covariance)), error = function(error) NULL)
  if (is.null(root)) stop("The regularized Laplace covariance is not SPD.")
  list(
    precision = precision,
    covariance = covariance,
    root = root,
    original_eigenvalue = decomposition$values,
    regularized_eigenvalue = regularized_value,
    eigenvalue_floor = floor_value,
    n_regularized = sum(decomposition$values < floor_value),
    maximum_added_precision = max(regularized_value - decomposition$values)
  )
}

find_collapsed_mode <- function(model, start = model$prior_mean,
                                minimum_eigenvalue = 1e-6,
                                maximum_condition = 1e8) {
  start <- flatten_theta(as_theta_matrix(start, model))
  objective <- function(x) -collapsed_log_posterior(x, model)
  gradient <- function(x) -collapsed_log_posterior_and_gradient(x, model)$gradient
  fit <- optim(start, objective, gradient, method = "BFGS",
               control = list(maxit = 2000L, reltol = 1e-12))
  if (fit$convergence != 0L || !is.finite(fit$value) ||
      any(!is.finite(fit$par))) {
    stop("Collapsed posterior mode failed: ", fit$message %||% fit$convergence)
  }
  observed_precision <- optimHess(fit$par, objective, gradient)
  metric <- regularize_laplace_metric(
    observed_precision, minimum_eigenvalue, maximum_condition
  )
  c(list(
    mode = as_theta_matrix(fit$par, model),
    log_posterior = -fit$value,
    observed_precision = observed_precision,
    optimizer = fit
  ), metric)
}

whiten_theta <- function(theta, geometry) {
  theta <- if (is.matrix(theta)) flatten_theta(theta) else as.numeric(theta)
  as.numeric(forwardsolve(geometry$root,
                          theta - flatten_theta(geometry$mode)))
}

inverse_whiten_theta <- function(position, geometry) {
  flatten_theta(geometry$mode) + as.numeric(geometry$root %*% position)
}

whitened_target <- function(position, model, geometry) {
  theta <- inverse_whiten_theta(position, geometry)
  evaluated <- tryCatch(
    collapsed_log_posterior_and_gradient(theta, model),
    error = function(error) NULL
  )
  if (is.null(evaluated) || !is.finite(evaluated$value) ||
      any(!is.finite(evaluated$gradient))) {
    return(list(value = -Inf, gradient = rep(NA_real_, length(position)),
                theta = theta))
  }
  list(
    value = evaluated$value,
    gradient = as.numeric(t(geometry$root) %*% evaluated$gradient),
    theta = theta
  )
}

hamiltonian <- function(log_posterior, momentum) {
  -log_posterior + 0.5 * sum(momentum^2)
}

leapfrog_integrate <- function(position, momentum, step_size, n_step,
                               model, geometry) {
  if (length(step_size) != 1L || !is.finite(step_size) || step_size <= 0 ||
      length(n_step) != 1L || n_step < 1L || n_step != as.integer(n_step)) {
    stop("Leapfrog step size and count must be positive and finite/integer.")
  }
  current <- whitened_target(position, model, geometry)
  if (!is.finite(current$value)) {
    return(list(valid = FALSE, position = position, momentum = momentum,
                target = current))
  }
  p <- momentum + 0.5 * step_size * current$gradient
  z <- position
  proposed <- current
  for (step in seq_len(as.integer(n_step))) {
    z <- z + step_size * p
    proposed <- whitened_target(z, model, geometry)
    if (!is.finite(proposed$value) || any(!is.finite(proposed$gradient))) {
      return(list(valid = FALSE, position = z, momentum = p,
                  target = proposed))
    }
    if (step < n_step) p <- p + step_size * proposed$gradient
  }
  p <- p + 0.5 * step_size * proposed$gradient
  list(valid = TRUE, position = z, momentum = p, target = proposed)
}

hmc_transition <- function(position, step_size, n_step, model, geometry) {
  current <- whitened_target(position, model, geometry)
  if (!is.finite(current$value)) stop("The current HMC state is not finite.")
  initial_momentum <- rnorm(length(position))
  proposal <- leapfrog_integrate(
    position, initial_momentum, step_size, n_step, model, geometry
  )
  current_hamiltonian <- hamiltonian(current$value, initial_momentum)
  if (proposal$valid) {
    proposed_hamiltonian <- hamiltonian(
      proposal$target$value, proposal$momentum
    )
    log_acceptance_ratio <- current_hamiltonian - proposed_hamiltonian
  } else {
    proposed_hamiltonian <- Inf
    log_acceptance_ratio <- -Inf
  }
  acceptance_probability <- exp(min(0, log_acceptance_ratio))
  accepted <- log(runif(1L)) < log_acceptance_ratio
  list(
    position = if (accepted) proposal$position else position,
    accepted = accepted,
    acceptance_probability = acceptance_probability,
    log_acceptance_ratio = log_acceptance_ratio,
    current_hamiltonian = current_hamiltonian,
    proposed_hamiltonian = proposed_hamiltonian,
    energy_error = proposed_hamiltonian - current_hamiltonian,
    proposal_valid = proposal$valid
  )
}

run_collapsed_hmc_chain <- function(
    model, geometry, seed, n_iteration = 800L, warmup = 300L, thin = 1L,
    step_size = 0.12, n_step = 10L, target_acceptance = 0.75,
    initial_theta = NULL) {
  n_iteration <- as.integer(n_iteration)
  warmup <- as.integer(warmup)
  thin <- as.integer(thin)
  if (n_iteration <= warmup || warmup < 0L || thin < 1L) {
    stop("Require n_iteration > warmup >= 0 and thin >= 1.")
  }
  set.seed(seed)
  position <- if (is.null(initial_theta)) {
    rnorm(model$parameter_dimension, sd = 1.25)
  } else {
    whiten_theta(initial_theta, geometry)
  }
  keep_iteration <- which(seq_len(n_iteration) > warmup &
                            ((seq_len(n_iteration) - warmup - 1L) %% thin == 0L))
  draw <- matrix(NA_real_, length(keep_iteration), model$parameter_dimension)
  accepted <- logical(n_iteration)
  acceptance_probability <- energy_error <- rep(NA_real_, n_iteration)
  proposal_valid <- logical(n_iteration)
  step_trace <- rep(NA_real_, n_iteration)
  keep_index <- 0L
  current_step <- step_size
  for (iteration in seq_len(n_iteration)) {
    transition <- hmc_transition(
      position, current_step, n_step, model, geometry
    )
    position <- transition$position
    accepted[iteration] <- transition$accepted
    acceptance_probability[iteration] <- transition$acceptance_probability
    energy_error[iteration] <- transition$energy_error
    proposal_valid[iteration] <- transition$proposal_valid
    step_trace[iteration] <- current_step
    if (iteration <= warmup) {
      learning_rate <- 0.05 / sqrt(iteration)
      current_step <- exp(log(current_step) + learning_rate *
                            (transition$acceptance_probability - target_acceptance))
      current_step <- min(max(current_step, 0.005), 0.5)
    }
    if (iteration %in% keep_iteration) {
      keep_index <- keep_index + 1L
      draw[keep_index, ] <- inverse_whiten_theta(position, geometry)
    }
  }
  colnames(draw) <- as.vector(t(outer(
    model$stick_id, model$coefficient_id, paste, sep = ":"
  )))
  post <- seq.int(warmup + 1L, n_iteration)
  list(
    theta = draw,
    alpha = lapply(seq_len(nrow(draw)), function(index) {
      theta_to_alpha(draw[index, ], model)
    }),
    seed = seed,
    settings = list(n_iteration = n_iteration, warmup = warmup, thin = thin,
                    n_step = n_step, initial_step_size = step_size,
                    final_step_size = current_step,
                    target_acceptance = target_acceptance),
    diagnostics = list(
      acceptance_rate = mean(accepted[post]),
      mean_acceptance_probability = mean(acceptance_probability[post]),
      mean_absolute_energy_error = mean(abs(energy_error[post]), na.rm = TRUE),
      maximum_absolute_energy_error = max(abs(energy_error[post]), na.rm = TRUE),
      invalid_proposals = sum(!proposal_valid[post]),
      accepted = accepted,
      acceptance_probability = acceptance_probability,
      energy_error = energy_error,
      step_size = step_trace
    )
  )
}

run_collapsed_hmc_chains <- function(model, geometry, seeds, ...) {
  lapply(as.integer(seeds), function(seed) {
    run_collapsed_hmc_chain(model, geometry, seed = seed, ...)
  })
}

split_rhat <- function(chains) {
  chains <- as.matrix(chains)
  n <- nrow(chains)
  if (ncol(chains) < 2L || n < 4L) return(NA_real_)
  half <- floor(n / 2L)
  split <- cbind(chains[seq_len(half), , drop = FALSE],
                 chains[(n - half + 1L):n, , drop = FALSE])
  within <- mean(apply(split, 2L, var))
  between <- half * var(colMeans(split))
  if (!is.finite(within) || within <= 0) return(NA_real_)
  sqrt(((half - 1) / half * within + between / half) / within)
}

effective_sample_size <- function(chains) {
  chains <- as.matrix(chains)
  n <- nrow(chains)
  m <- ncol(chains)
  if (n < 4L) return(NA_real_)
  autocorrelation <- vapply(seq_len(m), function(chain) {
    acf(chains[, chain], lag.max = n - 1L, plot = FALSE,
        demean = TRUE)$acf[-1L]
  }, numeric(n - 1L))
  rho <- rowMeans(autocorrelation)
  pair <- rho[seq.int(1L, length(rho) - 1L, by = 2L)] +
    rho[seq.int(2L, length(rho), by = 2L)]
  positive <- which(pair <= 0 | !is.finite(pair))
  if (length(positive)) pair <- pair[seq_len(positive[1L] - 1L)]
  tau <- 1 + 2 * sum(pair)
  min(m * n, m * n / max(tau, 1))
}

summarize_hmc_chains <- function(chains) {
  if (length(chains) < 2L) stop("At least two chains are required.")
  n_draw <- unique(vapply(chains, function(chain) nrow(chain$theta), integer(1)))
  if (length(n_draw) != 1L) stop("Chains must retain the same number of draws.")
  parameter <- colnames(chains[[1L]]$theta)
  parameter_summary <- do.call(rbind, lapply(seq_along(parameter), function(j) {
    x <- vapply(chains, function(chain) chain$theta[, j], numeric(n_draw))
    data.frame(
      parameter = parameter[j],
      mean = mean(x),
      sd = sd(as.numeric(x)),
      rhat = split_rhat(x),
      ess = effective_sample_size(x),
      stringsAsFactors = FALSE
    )
  }))
  chain_summary <- data.frame(
    chain = seq_along(chains),
    seed = vapply(chains, `[[`, integer(1), "seed"),
    acceptance_rate = vapply(chains, function(x) x$diagnostics$acceptance_rate,
                              numeric(1)),
    mean_acceptance_probability = vapply(
      chains, function(x) x$diagnostics$mean_acceptance_probability, numeric(1)
    ),
    final_step_size = vapply(
      chains, function(x) x$settings$final_step_size, numeric(1)
    ),
    mean_absolute_energy_error = vapply(
      chains, function(x) x$diagnostics$mean_absolute_energy_error, numeric(1)
    ),
    maximum_absolute_energy_error = vapply(
      chains, function(x) x$diagnostics$maximum_absolute_energy_error, numeric(1)
    ),
    invalid_proposals = vapply(
      chains, function(x) x$diagnostics$invalid_proposals, integer(1)
    )
  )
  list(parameter = parameter_summary, chain = chain_summary)
}

estimate_continuation_centers <- function(
    annotation, component_log_density, baseline_continuation,
    intercept_prior_sd = 2) {
  annotation <- as.matrix(annotation)
  component_log_density <- as.matrix(component_log_density)
  n_stick <- ncol(component_log_density) - 1L
  if (length(baseline_continuation) != n_stick ||
      any(baseline_continuation <= 0 | baseline_continuation >= 1)) {
    stop("Baseline continuation probabilities must match the sticks.")
  }
  prior_mean <- matrix(qnorm(baseline_continuation), n_stick, 1L)
  objective <- function(intercept) {
    eta <- matrix(intercept, nrow(annotation), n_stick, byrow = TRUE)
    log_q <- matrix(
      pnorm(eta, log.p = TRUE), nrow(annotation), n_stick
    )
    log_1mq <- matrix(
      pnorm(eta, lower.tail = FALSE, log.p = TRUE),
      nrow(annotation), n_stick
    )
    log_pi <- stick_breaking_log_probabilities(log_q, log_1mq)$log_probability
    log_weight <- log_pi + component_log_density
    marker_ll <- row_log_sum_exp(log_weight)
    responsibility <- exp(log_weight - marker_ll)
    log_phi <- dnorm(eta, log = TRUE)
    gradient <- numeric(n_stick)
    for (stick in seq_len(n_stick)) {
      success <- rowSums(responsibility[, (stick + 1L):(n_stick + 1L),
                                         drop = FALSE])
      failure <- responsibility[, stick]
      gradient[stick] <- sum(
        success * exp(log_phi[, stick] - log_q[, stick]) -
          failure * exp(log_phi[, stick] - log_1mq[, stick])
      )
    }
    prior_kernel <- -0.5 * sum(((intercept - prior_mean) /
                                  intercept_prior_sd)^2)
    gradient <- gradient - (intercept - prior_mean) / intercept_prior_sd^2
    list(value = sum(marker_ll) + prior_kernel, gradient = gradient,
         responsibility = responsibility)
  }
  fit <- optim(as.numeric(prior_mean), function(x) -objective(x)$value,
               function(x) -objective(x)$gradient, method = "BFGS",
               control = list(maxit = 1000L, reltol = 1e-12))
  if (fit$convergence != 0L || !is.finite(fit$value)) {
    stop("Intercept-only continuation-center fit failed.")
  }
  responsibility <- objective(fit$par)$responsibility
  reach <- matrix(1, nrow(annotation), n_stick)
  if (n_stick > 1L) {
    for (stick in 2:n_stick) {
      reach[, stick] <- rowSums(
        responsibility[, stick:(n_stick + 1L), drop = FALSE]
      )
    }
  }
  center <- t(vapply(seq_len(n_stick), function(stick) {
    colSums(annotation * reach[, stick]) / sum(reach[, stick])
  }, numeric(ncol(annotation))))
  list(center = center, intercept_mode = fit$par, reach_weight = reach)
}
