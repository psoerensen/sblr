# Single research entry point for deterministic R/Python qualification and a
# bounded gsim smoke analysis. All generated files stay under output/.

standardize_annotation <- function(x) {
  x <- as.matrix(x)
  center <- colMeans(x)
  scale <- apply(x, 2L, sd)
  if (any(!is.finite(scale)) || any(scale <= 1e-10)) {
    stop("Annotation columns must have positive finite sample SD.")
  }
  sweep(sweep(x, 2L, center, "-"), 2L, scale, "/")
}

make_annotation_alpha_fixture <- function(scenario, n_marker = 24L) {
  scenario <- match.arg(scenario, c(
    "null", "sparse_independent", "correlated_proxy", "later_stick"
  ))
  index <- seq_len(n_marker)
  signal <- qnorm((index - 0.5) / n_marker)
  proxy_noise <- scale(sin(index * 1.7))[, 1L]
  proxy <- 0.96 * scale(signal)[, 1L] +
    sqrt(1 - 0.96^2) * proxy_noise
  rare <- as.numeric(index <= max(3L, round(0.15 * n_marker)))
  annotation <- standardize_annotation(cbind(
    signal = signal, proxy = proxy, rare_binary = rare
  ))
  rownames(annotation) <- paste0("marker", index)
  slope <- matrix(0, 3L, 3L,
                  dimnames = list(colnames(annotation), paste0("stick", 1:3)))
  if (scenario == "sparse_independent") {
    slope["signal", ] <- c(0.70, 0.42, 0.25)
    slope["rare_binary", ] <- c(0.45, 0.25, 0.15)
  } else if (scenario == "correlated_proxy") {
    slope["signal", ] <- c(0.75, 0.45, 0.25)
  } else if (scenario == "later_stick") {
    slope["signal", ] <- c(0, 0.65, 0.42)
    slope["rare_binary", ] <- c(0, 0.35, 0.25)
  }
  target_continuation <- c(0.12, 0.50, 0.45)
  alpha <- rbind(intercept = numeric(3L), slope)
  continuation <- matrix(NA_real_, n_marker, 3L)
  reach <- rep(1, n_marker)
  for (stick in 1:3) {
    offset <- annotation %*% slope[, stick]
    objective <- function(intercept) {
      sum(reach * pnorm(intercept + offset)) / sum(reach) -
        target_continuation[stick]
    }
    alpha[1L, stick] <- uniroot(objective, c(-12, 12), tol = 1e-12)$root
    continuation[, stick] <- pnorm(alpha[1L, stick] + offset)
    reach <- reach * continuation[, stick]
  }
  colnames(continuation) <- paste0("stick", 1:3)
  rownames(continuation) <- rownames(annotation)
  component_probability <- stick_breaking_probabilities(continuation)$probability
  colnames(component_probability) <- paste0("component", 0:3)
  rownames(component_probability) <- rownames(annotation)
  mixture_scale <- c(0, 0.01, 0.10, 1.0)
  V_beta <- 0.08
  expected_scale <- as.numeric(component_probability %*% mixture_scale)
  summary_se <- 0.035 + 0.003 * (index %% 3L)
  summary_beta <- sqrt(V_beta * expected_scale + summary_se^2) *
    (0.55 * sin(index * sqrt(2)) + 0.20 * cos(index / 3))
  list(
    scenario = scenario,
    marker_id = rownames(annotation),
    annotation_id = colnames(annotation),
    annotation = annotation,
    alpha = alpha,
    continuation = continuation,
    component_probability = component_probability,
    target_continuation = target_continuation,
    mixture_scale = mixture_scale,
    V_beta = V_beta,
    summary_beta = summary_beta,
    summary_se = summary_se
  )
}

make_tail_qualification_fixture <- function(tail = c("negative", "positive")) {
  tail <- match.arg(tail)
  intercept <- if (tail == "negative") -8 else 8
  annotation <- matrix(c(-1, 1), 2L, 1L,
                       dimnames = list(c("tail_marker1", "tail_marker2"),
                                       "tail_annotation"))
  alpha <- rbind(intercept = rep(intercept, 3L), tail_annotation = 0)
  colnames(alpha) <- paste0("stick", 1:3)
  continuation <- matrix(
    pnorm(cbind(1, annotation) %*% alpha), 2L, 3L,
    dimnames = list(rownames(annotation), colnames(alpha))
  )
  mixture <- stick_breaking_probabilities(continuation)$probability
  dimnames(mixture) <- list(rownames(annotation), paste0("component", 0:3))
  list(
    scenario = paste0("tail_", tail),
    marker_id = rownames(annotation),
    annotation_id = colnames(annotation),
    annotation = annotation,
    alpha = alpha,
    continuation = continuation,
    component_probability = mixture,
    target_continuation = c(0.12, 0.50, 0.45),
    mixture_scale = c(0, 0.01, 0.10, 1.0),
    V_beta = 0.08,
    summary_beta = if (tail == "negative") rep(0.3, 2L) else rep(0, 2L),
    summary_se = if (tail == "negative") rep(0.01, 2L) else rep(1e-8, 2L),
    estimate_center = FALSE,
    qualification_theta = alpha
  )
}

make_size_one_qualification_fixture <- function(n_stick = c(3L, 1L)) {
  n_stick <- as.integer(match.arg(as.character(n_stick), c("3", "1")))
  target <- if (n_stick == 1L) 0.30 else c(0.12, 0.50, 0.45)
  mixture_scale <- if (n_stick == 1L) c(0, 1) else c(0, 0.01, 0.10, 1.0)
  annotation <- matrix(
    0.25, 1L, 1L,
    dimnames = list("singleton_marker", "singleton_annotation")
  )
  alpha <- rbind(
    intercept = qnorm(target),
    singleton_annotation = seq_len(n_stick) / 10
  )
  colnames(alpha) <- paste0("stick", seq_len(n_stick))
  continuation <- matrix(
    pnorm(cbind(1, annotation) %*% alpha), 1L, n_stick,
    dimnames = list(rownames(annotation), colnames(alpha))
  )
  mixture <- stick_breaking_probabilities(continuation)$probability
  dimnames(mixture) <- list(
    rownames(annotation), paste0("component", 0:n_stick)
  )
  list(
    scenario = paste0("size_one_", n_stick, "_stick"),
    marker_id = rownames(annotation),
    annotation_id = colnames(annotation),
    annotation = annotation,
    alpha = alpha,
    continuation = continuation,
    component_probability = mixture,
    target_continuation = target,
    mixture_scale = mixture_scale,
    V_beta = 0.08,
    summary_beta = 0.04,
    summary_se = 0.03,
    estimate_center = FALSE
  )
}

fixture_model <- function(fixture, estimate_center = TRUE) {
  n_component <- length(fixture$mixture_scale)
  n_stick <- n_component - 1L
  variance <- outer(fixture$summary_se^2, rep(1, n_component)) +
    outer(rep(fixture$V_beta, length(fixture$summary_beta)),
          fixture$mixture_scale)
  log_density <- -0.5 * (log(2 * pi * variance) +
                           fixture$summary_beta^2 / variance)
  center <- if (estimate_center) {
    estimate_continuation_centers(
      fixture$annotation, log_density, fixture$target_continuation
    )$center
  } else {
    matrix(0, n_stick, ncol(fixture$annotation))
  }
  prior_mean <- matrix(0, n_stick, ncol(fixture$annotation) + 1L)
  prior_mean[, 1L] <- qnorm(fixture$target_continuation)
  prior_sd <- matrix(0.75, n_stick, ncol(prior_mean))
  prior_sd[, 1L] <- 1.5
  make_collapsed_alpha_model(
    fixture$summary_beta, fixture$summary_se, fixture$V_beta,
    fixture$mixture_scale, fixture$annotation, center,
    prior_mean, prior_sd, fixture$marker_id, fixture$annotation_id,
    paste0("stick", seq_len(n_stick)),
    paste0("component", 0:(n_component - 1L))
  )
}

write_numeric_matrix <- function(value, path) {
  write.table(as.matrix(value), path, sep = ",", row.names = FALSE,
              col.names = FALSE, quote = FALSE, na = "")
}

read_numeric_matrix <- function(path) {
  as.matrix(read.csv(path, header = FALSE, check.names = FALSE))
}

find_working_python <- function() {
  candidate <- unique(unname(Sys.which(c("python", "python3"))))
  candidate <- candidate[nzchar(candidate)]
  for (command in candidate) {
    probe <- tryCatch(
      suppressWarnings(system2(command, "--version", stdout = TRUE,
                               stderr = TRUE)),
      error = function(error) structure(character(), status = 1L)
    )
    status <- attr(probe, "status") %||% 0L
    if (identical(status, 0L) && any(grepl("Python", probe, fixed = TRUE))) {
      return(command)
    }
  }
  ""
}

maximum_discrepancy <- function(x, y) {
  if (!identical(dim(as.matrix(x)), dim(as.matrix(y)))) return(Inf)
  max(abs(as.numeric(x) - as.numeric(y)))
}

parity_row <- function(scenario, quantity, discrepancy, tolerance,
                       comparable = TRUE) {
  data.frame(
    scenario = scenario,
    quantity = quantity,
    maximum_absolute_discrepancy = discrepancy,
    tolerance = tolerance,
    comparable = comparable,
    status = if (!comparable) "INFORMATIONAL" else
      if (is.finite(discrepancy) && discrepancy <= tolerance) "PASS" else "FAIL",
    stringsAsFactors = FALSE
  )
}

run_python_parity <- function(fixture, model, geometry, theta, python,
                              python_dir, output_dir) {
  scenario_dir <- file.path(output_dir, "python_parity", fixture$scenario)
  input_dir <- file.path(scenario_dir, "input")
  oracle_dir <- file.path(scenario_dir, "oracle")
  dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(oracle_dir, recursive = TRUE, showWarnings = FALSE)
  position <- seq_len(model$parameter_dimension) /
    (5 * model$parameter_dimension)
  momentum <- cos(seq_len(model$parameter_dimension)) / 3
  control <- c(0.025, 5)
  write_numeric_matrix(fixture$annotation, file.path(input_dir, "annotation.csv"))
  write_numeric_matrix(model$continuation_center, file.path(input_dir, "centers.csv"))
  write_numeric_matrix(as_theta_matrix(theta, model), file.path(input_dir, "theta.csv"))
  write_numeric_matrix(model$prior_mean, file.path(input_dir, "prior_mean.csv"))
  write_numeric_matrix(model$prior_sd, file.path(input_dir, "prior_sd.csv"))
  write_numeric_matrix(cbind(fixture$summary_beta, fixture$summary_se),
                       file.path(input_dir, "summary.csv"))
  write_numeric_matrix(c(fixture$V_beta, fixture$mixture_scale),
                       file.path(input_dir, "mixture.csv"))
  write_numeric_matrix(geometry$mode, file.path(input_dir, "mode_r.csv"))
  write_numeric_matrix(geometry$covariance,
                       file.path(input_dir, "covariance_r.csv"))
  write_numeric_matrix(position, file.path(input_dir, "position.csv"))
  write_numeric_matrix(momentum, file.path(input_dir, "momentum.csv"))
  write_numeric_matrix(control, file.path(input_dir, "hmc_control.csv"))
  script <- file.path(python_dir, "parity_oracle.py")
  status <- system2(
    python,
    c(shQuote(normalizePath(script)), shQuote(normalizePath(input_dir)),
      shQuote(normalizePath(oracle_dir))),
    stdout = TRUE, stderr = TRUE
  )
  if (!is.null(attr(status, "status")) && attr(status, "status") != 0L) {
    stop("Python parity adapter failed: ", paste(status, collapse = "\n"))
  }

  terms <- collapsed_marker_terms(theta, model)
  details <- collapsed_log_posterior_details(theta, model)
  evaluated <- collapsed_log_posterior_and_gradient(theta, model)
  leapfrog <- leapfrog_integrate(
    position, momentum, control[1L], control[2L], model, geometry
  )
  current <- whitened_target(position, model, geometry)
  current_h <- hamiltonian(current$value, momentum)
  proposed_h <- hamiltonian(leapfrog$target$value, leapfrog$momentum)
  log_accept <- current_h - proposed_h
  scalar_python <- as.numeric(read_numeric_matrix(
    file.path(oracle_dir, "scalars.csv")
  ))
  comparison <- list(
    parity_row(fixture$scenario, "linear_predictor",
      maximum_discrepancy(terms$linear_predictor,
                          read_numeric_matrix(file.path(oracle_dir, "linear_predictor.csv"))),
      2e-11),
    parity_row(fixture$scenario, "continuation_probability",
      maximum_discrepancy(terms$continuation_probability,
                          read_numeric_matrix(file.path(oracle_dir, "continuation_probability.csv"))),
      2e-11),
    parity_row(fixture$scenario, "log_continuation_probability",
      maximum_discrepancy(terms$log_continuation_probability,
                          read_numeric_matrix(file.path(oracle_dir, "log_continuation_probability.csv"))),
      3e-11),
    parity_row(fixture$scenario, "log_continuation_survival_probability",
      maximum_discrepancy(terms$log_continuation_survival_probability,
                          read_numeric_matrix(file.path(oracle_dir, "log_continuation_survival_probability.csv"))),
      3e-11),
    parity_row(fixture$scenario, "component_probability",
      maximum_discrepancy(terms$component_probability,
                          read_numeric_matrix(file.path(oracle_dir, "component_probability.csv"))),
      2e-11),
    parity_row(fixture$scenario, "log_component_probability",
      maximum_discrepancy(terms$log_component_probability,
                          read_numeric_matrix(file.path(oracle_dir, "log_component_probability.csv"))),
      3e-11),
    parity_row(fixture$scenario, "component_log_density",
      maximum_discrepancy(terms$component_log_density,
                          read_numeric_matrix(file.path(oracle_dir, "component_log_density.csv"))),
      2e-11),
    parity_row(fixture$scenario, "component_log_weight",
      maximum_discrepancy(terms$component_log_weight,
                          read_numeric_matrix(file.path(oracle_dir, "component_log_weight.csv"))),
      3e-11),
    parity_row(fixture$scenario, "marker_log_likelihood",
      maximum_discrepancy(terms$marker_log_likelihood,
                          read_numeric_matrix(file.path(oracle_dir, "marker_log_likelihood.csv"))),
      3e-11),
    parity_row(fixture$scenario, "responsibility",
      maximum_discrepancy(terms$responsibility,
                          read_numeric_matrix(file.path(oracle_dir, "responsibility.csv"))),
      3e-11),
    parity_row(fixture$scenario, "total_log_likelihood",
      abs(details$log_likelihood - scalar_python[1L]), 5e-10),
    parity_row(fixture$scenario, "log_prior_kernel",
      abs(details$log_prior_kernel - scalar_python[2L]), 2e-11),
    parity_row(fixture$scenario, "log_prior",
      abs(details$log_prior - scalar_python[3L]), 2e-11),
    parity_row(fixture$scenario, "log_posterior_kernel",
      abs(details$log_posterior_kernel - scalar_python[4L]), 5e-10),
    parity_row(fixture$scenario, "log_posterior",
      abs(details$log_posterior - scalar_python[5L]), 5e-10),
    parity_row(fixture$scenario, "analytic_gradient",
      maximum_discrepancy(as_theta_matrix(evaluated$gradient, model),
                          read_numeric_matrix(file.path(oracle_dir, "gradient.csv"))),
      2e-8),
    parity_row(fixture$scenario, "laplace_mode",
      maximum_discrepancy(geometry$mode,
                          read_numeric_matrix(file.path(oracle_dir, "mode_python.csv"))),
      2e-4),
    parity_row(fixture$scenario, "laplace_covariance_algorithm_difference",
      maximum_discrepancy(geometry$covariance,
                          read_numeric_matrix(file.path(oracle_dir, "covariance_python.csv"))),
      NA_real_, comparable = FALSE),
    parity_row(fixture$scenario, "whitening",
      maximum_discrepancy(whiten_theta(theta, geometry),
                          read_numeric_matrix(file.path(oracle_dir, "whitened.csv"))),
      2e-10),
    parity_row(fixture$scenario, "inverse_whitening",
      maximum_discrepancy(theta,
                          read_numeric_matrix(file.path(oracle_dir, "inverse_whitened.csv"))),
      2e-10),
    parity_row(fixture$scenario, "leapfrog_position",
      maximum_discrepancy(leapfrog$position,
                          read_numeric_matrix(file.path(oracle_dir, "leapfrog_position.csv"))),
      2e-9),
    parity_row(fixture$scenario, "leapfrog_momentum",
      maximum_discrepancy(leapfrog$momentum,
                          read_numeric_matrix(file.path(oracle_dir, "leapfrog_momentum.csv"))),
      2e-9),
    parity_row(fixture$scenario, "current_hamiltonian",
      abs(current_h - scalar_python[6L]), 2e-9),
    parity_row(fixture$scenario, "proposed_hamiltonian",
      abs(proposed_h - scalar_python[7L]), 2e-9),
    parity_row(fixture$scenario, "log_acceptance_ratio",
      abs(log_accept - scalar_python[8L]), 3e-9),
    parity_row(fixture$scenario, "acceptance_probability",
      abs(exp(min(0, log_accept)) - scalar_python[9L]), 3e-9)
  )
  do.call(rbind, comparison)
}

run_deterministic_qualification <- function(output_dir, python, python_dir) {
  fixtures <- c(
    lapply(c("null", "sparse_independent", "correlated_proxy", "later_stick"),
           make_annotation_alpha_fixture),
    list(
      make_tail_qualification_fixture("negative"),
      make_tail_qualification_fixture("positive"),
      make_size_one_qualification_fixture(3L),
      make_size_one_qualification_fixture(1L)
    )
  )
  checks <- list()
  parity <- list()
  first <- NULL
  for (fixture in fixtures) {
    scenario <- fixture$scenario
    model <- fixture_model(
      fixture, estimate_center = fixture$estimate_center %||% TRUE
    )
    theta_truth <- alpha_to_theta(fixture$alpha, model)
    theta <- if (!is.null(fixture$qualification_theta)) {
      flatten_theta(alpha_to_theta(fixture$qualification_theta, model))
    } else {
      flatten_theta(0.35 * theta_truth + 0.65 * model$prior_mean) +
        seq_len(model$parameter_dimension) / 1000
    }
    geometry <- find_collapsed_mode(model)
    gradient_points <- list(
      theta,
      theta + sin(seq_along(theta)) / 20,
      flatten_theta(model$prior_mean) + cos(seq_along(theta)) / 25
    )
    gradient_error <- max(vapply(gradient_points, function(point) {
      check_collapsed_gradient(point, model)$maximum_absolute_discrepancy
    }, numeric(1)))
    terms <- collapsed_marker_terms(theta, model)
    plateau_change <- NA_real_
    if (startsWith(scenario, "tail_")) {
      plus <- minus <- theta
      plus[1L] <- plus[1L] + 5e-5
      minus[1L] <- minus[1L] - 5e-5
      plateau_change <- abs(
        collapsed_marker_terms(plus, model)$log_likelihood -
          collapsed_marker_terms(minus, model)$log_likelihood
      )
    }
    alpha_roundtrip <- theta_to_alpha(alpha_to_theta(fixture$alpha, model), model)
    position <- seq_len(model$parameter_dimension) /
      (4 * model$parameter_dimension)
    momentum <- sin(seq_len(model$parameter_dimension)) / 4
    leapfrog <- leapfrog_integrate(position, momentum, 0.02, 6L, model, geometry)
    reverse <- leapfrog_integrate(
      leapfrog$position, -leapfrog$momentum, 0.02, 6L, model, geometry
    )
    current <- whitened_target(position, model, geometry)
    log_accept <- hamiltonian(current$value, momentum) -
      hamiltonian(leapfrog$target$value, leapfrog$momentum)
    full_difference <- (
      leapfrog$target$value - 0.5 * sum(leapfrog$momentum^2)
    ) - (current$value - 0.5 * sum(momentum^2))
    checks[[scenario]] <- data.frame(
      scenario = scenario,
      maximum_gradient_error = gradient_error,
      maximum_component_sum_error = max(abs(rowSums(
        terms$component_probability) - 1)),
      maximum_responsibility_sum_error = max(abs(rowSums(
        terms$responsibility) - 1)),
      maximum_probability_bound_error = max(
        0, -min(terms$component_probability),
        max(terms$component_probability) - 1
      ),
      finite_log_posterior = is.finite(collapsed_log_posterior(theta, model)),
      finite_log_component_probability = all(is.finite(
        terms$log_component_probability
      )),
      tail_likelihood_change = plateau_change,
      declared_dimensions_preserved =
        identical(dim(terms$continuation_probability),
                  c(model$n_marker, model$n_stick)) &&
        identical(dim(terms$component_probability),
                  c(model$n_marker, model$n_component)) &&
        identical(dim(terms$component_log_density),
                  c(model$n_marker, model$n_component)) &&
        identical(dim(terms$responsibility),
                  c(model$n_marker, model$n_component)),
      declared_ids_preserved =
        identical(dimnames(terms$continuation_probability),
                  list(model$marker_id, model$stick_id)) &&
        identical(dimnames(terms$component_probability),
                  list(model$marker_id, model$component_id)) &&
        identical(dimnames(terms$responsibility),
                  list(model$marker_id, model$component_id)),
      alpha_theta_roundtrip_error = max(abs(alpha_roundtrip - fixture$alpha)),
      whitening_roundtrip_error = max(abs(
        inverse_whiten_theta(whiten_theta(theta, geometry), geometry) - theta
      )),
      leapfrog_position_reversal_error = max(abs(reverse$position - position)),
      leapfrog_momentum_reversal_error = max(abs(reverse$momentum + momentum)),
      full_hamiltonian_acceptance_error = abs(log_accept - full_difference),
      metric_regularized_eigenvalues = geometry$n_regularized,
      metric_eigenvalue_floor = geometry$eigenvalue_floor,
      stringsAsFactors = FALSE
    )
    if (nzchar(python)) {
      parity[[scenario]] <- run_python_parity(
        fixture, model, geometry, theta, python, python_dir, output_dir
      )
    } else {
      parity[[scenario]] <- data.frame(
        scenario = scenario,
        quantity = "python_oracle_execution",
        maximum_absolute_discrepancy = NA_real_, tolerance = NA_real_,
        comparable = FALSE, status = "SKIP_PYTHON_UNAVAILABLE",
        stringsAsFactors = FALSE
      )
    }
    if (is.null(first)) first <- list(model = model, geometry = geometry)
  }
  reproducible_a <- run_collapsed_hmc_chain(
    first$model, first$geometry, seed = 20260831L,
    n_iteration = 100L, warmup = 40L, thin = 2L,
    step_size = 0.10, n_step = 8L
  )
  reproducible_b <- run_collapsed_hmc_chain(
    first$model, first$geometry, seed = 20260831L,
    n_iteration = 100L, warmup = 40L, thin = 2L,
    step_size = 0.10, n_step = 8L
  )
  reproducibility <- data.frame(
    scenario = "null",
    maximum_same_seed_draw_error = max(abs(
      reproducible_a$theta - reproducible_b$theta
    )),
    identical_acceptance_trace = identical(
      reproducible_a$diagnostics$accepted,
      reproducible_b$diagnostics$accepted
    )
  )
  list(checks = do.call(rbind, checks), parity = do.call(rbind, parity),
       reproducibility = reproducibility)
}

make_gsim_annotation_adapter <- function(component_probability) {
  ordered_index <- c(4L, 3L, 2L, 1L)
  ordered <- component_probability[, ordered_index, drop = FALSE]
  continuation <- matrix(NA_real_, nrow(ordered), 3L)
  for (stick in 1:3) {
    continuation[, stick] <- ordered[, stick] /
      rowSums(ordered[, stick:4, drop = FALSE])
  }
  eta <- qlogis(continuation)
  baseline <- plogis(colMeans(eta))
  ordered_global <- numeric(4L)
  remaining <- 1
  for (stick in 1:3) {
    ordered_global[stick] <- remaining * baseline[stick]
    remaining <- remaining * (1 - baseline[stick])
  }
  ordered_global[4L] <- remaining
  global_probability <- numeric(4L)
  global_probability[ordered_index] <- ordered_global
  annotation <- eta
  colnames(annotation) <- paste0("adapter_stick", 1:3)
  alpha <- diag(3L)
  rownames(alpha) <- colnames(annotation)
  colnames(alpha) <- paste0("stick", 1:3)
  list(annotation = annotation, alpha = alpha, pi = global_probability,
       target_probability = component_probability)
}

git_value <- function(repository, arguments) {
  value <- system2("git", c("-C", shQuote(normalizePath(repository)), arguments),
                   stdout = TRUE, stderr = TRUE)
  if (!is.null(attr(value, "status")) && attr(value, "status") != 0L) {
    stop("Git provenance query failed for ", repository)
  }
  paste(value, collapse = "\n")
}

run_gsim_smoke <- function(sblrbench_root, output_dir) {
  gsim_environment <- new.env(parent = globalenv())
  source(file.path(sblrbench_root, "R", "gsim_internal.R"),
         local = gsim_environment)
  source(file.path(sblrbench_root, "R", "gsim.R"),
         local = gsim_environment)
  gsim_head <- git_value(sblrbench_root, c("rev-parse", "HEAD"))
  scenarios <- c("null", "sparse_independent", "correlated_proxy",
                 "later_stick")
  gsim_seed <- setNames(20260920L + seq_along(scenarios), scenarios)
  scenario_summary <- parameter_summary <- chain_summary <- list()
  truth_surface <- truth_alpha <- manifest <- list()

  for (scenario in scenarios) {
    fixture <- make_annotation_alpha_fixture(scenario, n_marker = 60L)
    adapter <- make_gsim_annotation_adapter(fixture$component_probability)
    marker_id <- paste0("m", seq_len(60L))
    rownames(adapter$annotation) <- marker_id
    simulation <- gsim_environment$gsim(
      A = adapter$annotation,
      architecture = "bayesr",
      annotation_model = "sbayesrc",
      alpha = adapter$alpha,
      pi = adapter$pi,
      mixture_variances = fixture$mixture_scale,
      n = 400L,
      m = 60L,
      h2 = 0.40,
      seed = gsim_seed[[scenario]],
      scale_effects = FALSE,
      compute_sumstats = TRUE,
      return_marker_probabilities = TRUE
    )
    probability_error <- max(abs(
      simulation$marker_probabilities - fixture$component_probability
    ))
    summary <- simulation$sumstats
    summary <- summary[match(marker_id, summary$rsid), ]
    smoke_fixture <- fixture
    smoke_fixture$marker_id <- marker_id
    rownames(smoke_fixture$annotation) <- marker_id
    rownames(smoke_fixture$continuation) <- marker_id
    rownames(smoke_fixture$component_probability) <- marker_id
    smoke_fixture$summary_beta <- summary$beta
    smoke_fixture$summary_se <- summary$se
    smoke_fixture$V_beta <- 1
    model <- fixture_model(smoke_fixture)
    geometry <- find_collapsed_mode(model)
    chain_seed <- gsim_seed[[scenario]] + 1000L + 1:3
    chains <- run_collapsed_hmc_chains(
      model, geometry, seeds = chain_seed,
      n_iteration = 500L, warmup = 200L, thin = 2L,
      step_size = 0.25, n_step = 8L, target_acceptance = 0.75
    )
    diagnostics <- summarize_hmc_chains(chains)
    draws <- do.call(rbind, lapply(chains, `[[`, "theta"))
    theta_truth <- alpha_to_theta(smoke_fixture$alpha, model)
    posterior_mean <- colMeans(draws)
    truth_flat <- flatten_theta(theta_truth)
    parameter_summary[[scenario]] <- data.frame(
      scenario = scenario,
      parameter = colnames(draws),
      truth = truth_flat,
      posterior_mean = posterior_mean,
      posterior_sd = apply(draws, 2L, sd),
      stringsAsFactors = FALSE
    )
    posterior_q <- matrix(0, model$n_marker, model$n_stick)
    posterior_pi <- matrix(0, model$n_marker, model$n_component)
    for (draw_index in seq_len(nrow(draws))) {
      term <- collapsed_marker_terms(draws[draw_index, ], model)
      posterior_q <- posterior_q + term$continuation_probability
      posterior_pi <- posterior_pi + term$component_probability
    }
    posterior_q <- posterior_q / nrow(draws)
    posterior_pi <- posterior_pi / nrow(draws)
    alpha_mean <- theta_to_alpha(posterior_mean, model)
    scenario_summary[[scenario]] <- data.frame(
      scenario = scenario,
      gsim_seed = gsim_seed[[scenario]],
      n_marker = model$n_marker,
      n_annotation = model$n_annotation,
      realized_causal_markers = simulation$settings$n_causal,
      n_chain = length(chains),
      retained_draws_per_chain = nrow(chains[[1L]]$theta),
      gsim_probability_adapter_error = probability_error,
      alpha_rmse = sqrt(mean((alpha_mean - smoke_fixture$alpha)^2)),
      continuation_probability_rmse = sqrt(mean(
        (posterior_q - smoke_fixture$continuation)^2
      )),
      component_probability_rmse = sqrt(mean(
        (posterior_pi - smoke_fixture$component_probability)^2
      )),
      maximum_rhat = max(diagnostics$parameter$rhat, na.rm = TRUE),
      minimum_ess = min(diagnostics$parameter$ess, na.rm = TRUE),
      mean_acceptance_rate = mean(diagnostics$chain$acceptance_rate),
      maximum_absolute_energy_error = max(
        diagnostics$chain$maximum_absolute_energy_error
      ),
      invalid_proposals = sum(diagnostics$chain$invalid_proposals),
      stringsAsFactors = FALSE
    )
    chain_summary[[scenario]] <- transform(
      diagnostics$chain, scenario = scenario
    )
    truth_surface[[scenario]] <- data.frame(
      scenario = scenario,
      marker_id = marker_id,
      annotation_signal = smoke_fixture$annotation[, "signal"],
      annotation_proxy = smoke_fixture$annotation[, "proxy"],
      annotation_rare_binary = smoke_fixture$annotation[, "rare_binary"],
      summary_beta = smoke_fixture$summary_beta,
      summary_se = smoke_fixture$summary_se,
      q_stick1 = smoke_fixture$continuation[, 1L],
      q_stick2 = smoke_fixture$continuation[, 2L],
      q_stick3 = smoke_fixture$continuation[, 3L],
      pi_component0 = smoke_fixture$component_probability[, 1L],
      pi_component1 = smoke_fixture$component_probability[, 2L],
      pi_component2 = smoke_fixture$component_probability[, 3L],
      pi_component3 = smoke_fixture$component_probability[, 4L]
    )
    truth_alpha[[scenario]] <- data.frame(
      scenario = scenario,
      coefficient = rep(rownames(smoke_fixture$alpha), 3L),
      stick = rep(colnames(smoke_fixture$alpha), each = nrow(smoke_fixture$alpha)),
      alpha = as.numeric(smoke_fixture$alpha)
    )
    manifest[[scenario]] <- data.frame(
      scenario = scenario,
      gsim_seed = gsim_seed[[scenario]],
      sblrbench_head = gsim_head,
      realized_causal_markers = simulation$settings$n_causal,
      mixture_scale = paste(smoke_fixture$mixture_scale, collapse = ";"),
      V_beta = smoke_fixture$V_beta,
      marker_ids = paste(marker_id, collapse = ";"),
      annotation_ids = paste(smoke_fixture$annotation_id, collapse = ";"),
      stringsAsFactors = FALSE
    )
  }
  list(
    scenario = do.call(rbind, scenario_summary),
    parameter = do.call(rbind, parameter_summary),
    chain = do.call(rbind, chain_summary),
    truth_surface = do.call(rbind, truth_surface),
    truth_alpha = do.call(rbind, truth_alpha),
    manifest = do.call(rbind, manifest),
    sblrbench_head = gsim_head
  )
}

write_analysis_summary <- function(path, decision, python, deterministic, smoke,
                                   sblrbench_head) {
  parity_status <- table(deterministic$parity$status)
  lines <- c(
    "# Continuous-alpha prototype qualification",
    "",
    paste0("Decision: **", decision, "**"),
    "",
    paste0("Python: ", if (nzchar(python)) python else "unavailable"),
    paste0("sblrbench HEAD: `", sblrbench_head, "`"),
    "",
    "## Deterministic checks",
    "",
    paste0("Maximum analytic/finite-difference gradient error: ",
           format(max(deterministic$checks$maximum_gradient_error), digits = 6)),
    paste0("Maximum probability normalization error: ",
           format(max(deterministic$checks$maximum_component_sum_error), digits = 6)),
    paste0("Same-seed draw error: ",
           deterministic$reproducibility$maximum_same_seed_draw_error),
    paste0("Parity statuses: ", paste(names(parity_status), parity_status,
                                      sep = "=", collapse = ", ")),
    "",
    "## Bounded gsim smoke",
    "",
    paste(capture.output(print(smoke$scenario, row.names = FALSE)),
          collapse = "\n"),
    "",
    "These smoke results are research diagnostics, not accepted benchmark evidence."
  )
  writeLines(lines, path, useBytes = TRUE)
}

run_annotation_alpha_analysis <- function(project_root = ".") {
  project_root <- normalizePath(project_root, mustWork = TRUE)
  research_root <- file.path(project_root, "research", "sbayesrc",
    "continuous_alpha_hmc")
  output_dir <- file.path(research_root, "output")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  output_normalized <- normalizePath(output_dir, mustWork = TRUE)
  if (!startsWith(output_normalized,
                  normalizePath(research_root, mustWork = TRUE))) {
    stop("Resolved output path escaped the annotation-alpha workspace.")
  }
  source(file.path(research_root, "prototype.R"),
         local = environment(make_annotation_alpha_fixture))
  sblrbench_root <- normalizePath(file.path(project_root, "..", "sblrbench"),
                                  mustWork = TRUE)
  python_dir <- file.path(research_root, "python")
  python <- find_working_python()
  deterministic <- run_deterministic_qualification(
    output_dir, python, python_dir
  )
  smoke <- run_gsim_smoke(sblrbench_root, output_dir)
  write.csv(deterministic$checks,
            file.path(output_dir, "deterministic_r_checks.csv"), row.names = FALSE)
  write.csv(deterministic$parity,
            file.path(output_dir, "deterministic_python_parity.csv"), row.names = FALSE)
  write.csv(deterministic$reproducibility,
            file.path(output_dir, "r_seed_reproducibility.csv"), row.names = FALSE)
  write.csv(smoke$scenario, file.path(output_dir, "gsim_smoke_scenarios.csv"),
            row.names = FALSE)
  write.csv(smoke$parameter, file.path(output_dir, "gsim_smoke_parameters.csv"),
            row.names = FALSE)
  write.csv(smoke$chain, file.path(output_dir, "gsim_smoke_chains.csv"),
            row.names = FALSE)
  write.csv(smoke$truth_surface, file.path(output_dir, "gsim_smoke_truth.csv"),
            row.names = FALSE)
  write.csv(smoke$truth_alpha, file.path(output_dir, "gsim_smoke_alpha_truth.csv"),
            row.names = FALSE)
  write.csv(smoke$manifest, file.path(output_dir, "gsim_smoke_manifest.csv"),
            row.names = FALSE)

  deterministic_pass <-
    max(deterministic$checks$maximum_gradient_error) <= 2e-5 &&
    max(deterministic$checks$maximum_component_sum_error) <= 1e-12 &&
    max(deterministic$checks$maximum_responsibility_sum_error) <= 1e-12 &&
    max(deterministic$checks$maximum_probability_bound_error) == 0 &&
    all(deterministic$checks$finite_log_posterior) &&
    all(deterministic$checks$finite_log_component_probability) &&
    all(deterministic$checks$declared_dimensions_preserved) &&
    all(deterministic$checks$declared_ids_preserved) &&
    all(deterministic$checks$tail_likelihood_change[
      !is.na(deterministic$checks$tail_likelihood_change)
    ] > 0) &&
    max(deterministic$checks$whitening_roundtrip_error) <= 1e-10 &&
    max(deterministic$checks$leapfrog_position_reversal_error) <= 1e-9 &&
    max(deterministic$checks$full_hamiltonian_acceptance_error) <= 1e-12 &&
    deterministic$reproducibility$maximum_same_seed_draw_error == 0 &&
    deterministic$reproducibility$identical_acceptance_trace
  parity_pass <- nzchar(python) &&
    !any(deterministic$parity$status == "FAIL")
  decision <- if (deterministic_pass && parity_pass) {
    "READY FOR FINAL INDEPENDENT CONTINUOUS-ALPHA PROTOTYPE VERIFICATION"
  } else if (!nzchar(python)) {
    "NOT READY — Python interpreter unavailable; deterministic oracle parity was not executed"
  } else {
    "NOT READY — deterministic mathematical qualification failed"
  }
  write_analysis_summary(
    file.path(output_dir, "summary.md"), decision, python,
    deterministic, smoke, smoke$sblrbench_head
  )
  list(decision = decision, python = python,
       deterministic = deterministic, smoke = smoke,
       output_dir = output_normalized)
}

if (sys.nframe() == 0L) {
  result <- run_annotation_alpha_analysis(".")
  cat(result$decision, "\n")
  cat("Output:", result$output_dir, "\n")
}
