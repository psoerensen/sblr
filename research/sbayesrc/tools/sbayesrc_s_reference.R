# SBayesRC-S Phase-1 fixed-z reference
#
# Development/reference implementation.
# Not a production SBayesRC-S sampler.
# Not part of the supported public API.
#
# This script is intentionally self-contained and uses base R only. It
# independently enumerates the fixed-z posterior and validates a collapsed
# selection sampler against that exact finite-state oracle.

log_sum_exp <- function(x) {
  maximum <- max(x)
  maximum + log(sum(exp(x - maximum)))
}

model_states <- function(annotation_count) {
  states <- as.matrix(expand.grid(rep(list(0:1), annotation_count)))
  storage.mode(states) <- "integer"
  rownames(states) <- apply(states, 1L, paste0, collapse = "")
  states
}

component_probability <- function(q) {
  q <- as.matrix(q)
  probability <- matrix(0, nrow(q), ncol(q) + 1L)
  remaining <- rep(1, nrow(q))
  for (stick in seq_len(ncol(q))) {
    probability[, stick] <- remaining * (1 - q[, stick])
    remaining <- remaining * q[, stick]
  }
  probability[, ncol(probability)] <- remaining
  probability
}

stick_model <- function(z, annotation, selected, tau2) {
  selected <- as.logical(selected)
  design <- cbind(Intercept = 1, annotation[, selected, drop = FALSE])
  precision <- crossprod(design) +
    diag(c(0, rep(tau2^-1, sum(selected))), nrow = ncol(design))
  rhs <- drop(crossprod(design, z))
  chol_precision <- chol(precision)
  mean <- backsolve(chol_precision, forwardsolve(t(chol_precision), rhs))
  covariance <- chol2inv(chol_precision)
  log_marginal <- -0.5 * sum(selected) * log(tau2) -
    sum(log(diag(chol_precision))) + 0.5 * sum(rhs * mean)
  list(
    log_marginal = log_marginal,
    mean = mean,
    covariance = covariance,
    precision = precision,
    chol_precision = chol_precision,
    selected = which(selected)
  )
}

exact_reference <- function(z, annotation, eligible, pi_a, tau2) {
  annotation <- as.matrix(annotation)
  annotation_count <- ncol(annotation)
  stick_count <- length(z)
  states <- model_states(annotation_count)
  log_weight <- numeric(nrow(states))
  posterior <- vector("list", nrow(states))
  for (model in seq_len(nrow(states))) {
    selected <- states[model, ] == 1L
    log_weight[model] <- sum(dbinom(states[model, ], 1L, pi_a, log = TRUE))
    posterior[[model]] <- vector("list", stick_count)
    for (stick in seq_len(stick_count)) {
      rows <- eligible[[stick]]
      posterior[[model]][[stick]] <- stick_model(
        z[[stick]], annotation[rows, , drop = FALSE], selected, tau2[stick]
      )
      log_weight[model] <- log_weight[model] +
        posterior[[model]][[stick]]$log_marginal
    }
  }
  model_probability <- exp(log_weight - log_sum_exp(log_weight))
  names(model_probability) <- rownames(states)
  annotation_pip <- drop(crossprod(model_probability, states))
  names(annotation_pip) <- colnames(annotation)

  alpha_mean <- matrix(0, annotation_count, stick_count,
                       dimnames = list(colnames(annotation), paste0("stick", seq_len(stick_count))))
  alpha_included_numerator <- alpha_mean
  intercept_mean <- numeric(stick_count)
  q_by_model <- array(NA_real_, c(nrow(annotation), stick_count, nrow(states)))
  pi_by_model <- array(NA_real_, c(nrow(annotation), stick_count + 1L, nrow(states)))
  for (model in seq_len(nrow(states))) {
    selected <- states[model, ] == 1L
    for (stick in seq_len(stick_count)) {
      post <- posterior[[model]][[stick]]
      intercept_mean[stick] <- intercept_mean[stick] +
        model_probability[model] * post$mean[1L]
      if (any(selected)) {
        alpha_mean[selected, stick] <- alpha_mean[selected, stick] +
          model_probability[model] * post$mean[-1L]
        alpha_included_numerator[selected, stick] <-
          alpha_included_numerator[selected, stick] +
          model_probability[model] * post$mean[-1L]
      }
      prediction_design <- cbind(1, annotation[, selected, drop = FALSE])
      eta_mean <- drop(prediction_design %*% post$mean)
      eta_variance <- rowSums(
        (prediction_design %*% post$covariance) * prediction_design
      )
      q_by_model[, stick, model] <- pnorm(
        eta_mean / sqrt(1 + eta_variance)
      )
    }
    pi_by_model[, , model] <- component_probability(q_by_model[, , model])
  }
  alpha_mean_given_inclusion <- alpha_included_numerator /
    matrix(annotation_pip, annotation_count, stick_count)
  q_mean <- apply(
    q_by_model * rep(model_probability, each = nrow(annotation) * stick_count),
    c(1L, 2L), sum
  )
  pi_mean <- apply(
    pi_by_model * rep(
      model_probability, each = nrow(annotation) * (stick_count + 1L)
    ),
    c(1L, 2L), sum
  )
  list(
    states = states,
    model_probability = model_probability,
    annotation_pip = annotation_pip,
    intercept_mean = intercept_mean,
    alpha_mean = alpha_mean,
    alpha_mean_given_inclusion = alpha_mean_given_inclusion,
    q_mean = q_mean,
    pi_mean = pi_mean,
    posterior = posterior,
    log_weight = log_weight
  )
}

log_bf_formula <- function(x, residual, tau2) {
  s <- sum(x * x)
  t <- sum(x * residual)
  -0.5 * log1p(tau2 * s) +
    0.5 * tau2 * t^2 / (1 + tau2 * s)
}

log_bf_direct <- function(x, residual, tau2) {
  covariance <- diag(length(x)) + tau2 * tcrossprod(x)
  chol_covariance <- chol(covariance)
  solved <- backsolve(
    chol_covariance,
    forwardsolve(t(chol_covariance), residual)
  )
  -sum(log(diag(chol_covariance))) -
    0.5 * sum(residual * solved) + 0.5 * sum(residual^2)
}

draw_from_stick_model <- function(posterior) {
  draw <- posterior$mean + backsolve(
    posterior$chol_precision, rnorm(length(posterior$mean))
  )
  list(intercept = draw[1L], slopes = draw[-1L])
}

run_chain <- function(z, annotation, eligible, pi_a, tau2,
                      iterations, burn, initial_delta) {
  annotation_count <- ncol(annotation)
  stick_count <- length(z)
  delta <- as.integer(initial_delta)
  alpha <- matrix(0, annotation_count, stick_count)
  intercept <- numeric(stick_count)
  retained <- iterations - burn
  delta_draws <- matrix(NA_integer_, retained, annotation_count)
  alpha_draws <- array(NA_real_, c(retained, annotation_count, stick_count))
  intercept_draws <- matrix(NA_real_, retained, stick_count)
  q_sum <- matrix(0, nrow(annotation), stick_count)
  pi_sum <- matrix(0, nrow(annotation), stick_count + 1L)
  states <- model_states(annotation_count)
  posterior_cache <- lapply(seq_len(nrow(states)), function(model) {
    selected <- states[model, ] == 1L
    lapply(seq_len(stick_count), function(stick) {
      rows <- eligible[[stick]]
      stick_model(
        z[[stick]], annotation[rows, , drop = FALSE], selected, tau2[stick]
      )
    })
  })

  kept <- 0L
  for (iteration in seq_len(iterations)) {
    # Each update samples the joint conditional of delta_j and all of its
    # stick-specific slopes, with that slope block analytically integrated in
    # the Bernoulli probability and then regenerated conditional on inclusion.
    for (j in seq_len(annotation_count)) {
      log_odds <- qlogis(pi_a)
      moments <- vector("list", stick_count)
      for (stick in seq_len(stick_count)) {
        rows <- eligible[[stick]]
        other <- setdiff(seq_len(annotation_count), j)
        residual <- z[[stick]] - intercept[stick]
        if (length(other)) {
          residual <- residual - drop(
            annotation[rows, other, drop = FALSE] %*% alpha[other, stick]
          )
        }
        x <- annotation[rows, j]
        s <- sum(x * x)
        t <- sum(x * residual)
        variance <- 1 / (s + tau2[stick]^-1)
        moments[[stick]] <- c(mean = variance * t, variance = variance)
        log_odds <- log_odds + log_bf_formula(x, residual, tau2[stick])
      }
      delta[j] <- rbinom(1L, 1L, plogis(log_odds))
      if (delta[j] == 1L) {
        for (stick in seq_len(stick_count)) {
          alpha[j, stick] <- rnorm(
            1L, moments[[stick]]["mean"], sqrt(moments[[stick]]["variance"])
          )
        }
      } else {
        alpha[j, ] <- 0
      }
    }

    # A separate exact Gibbs step redraws the intercept and every included
    # slope jointly, conditional on the completed selection sweep.
    selected <- delta == 1L
    model <- match(paste0(delta, collapse = ""), rownames(states))
    for (stick in seq_len(stick_count)) {
      draw <- draw_from_stick_model(posterior_cache[[model]][[stick]])
      intercept[stick] <- draw$intercept
      alpha[, stick] <- 0
      if (any(selected)) alpha[selected, stick] <- draw$slopes
    }

    if (iteration > burn) {
      kept <- kept + 1L
      delta_draws[kept, ] <- delta
      alpha_draws[kept, , ] <- alpha
      intercept_draws[kept, ] <- intercept
      q <- pnorm(sweep(annotation %*% alpha, 2L, intercept, `+`))
      q_sum <- q_sum + q
      pi_sum <- pi_sum + component_probability(q)
    }
  }
  list(
    delta_draws = delta_draws,
    alpha_draws = alpha_draws,
    intercept_draws = intercept_draws,
    q_mean = q_sum / retained,
    pi_mean = pi_sum / retained
  )
}

pool_chains <- function(chains, annotation_names) {
  delta <- do.call(rbind, lapply(chains, `[[`, "delta_draws"))
  chain_sizes <- vapply(chains, function(x) dim(x$alpha_draws)[1L], integer(1L))
  alpha <- array(NA_real_, c(sum(chain_sizes), dim(chains[[1L]]$alpha_draws)[2:3]))
  offset <- 0L
  for (chain in seq_along(chains)) {
    rows <- offset + seq_len(chain_sizes[chain])
    alpha[rows, , ] <- chains[[chain]]$alpha_draws
    offset <- offset + chain_sizes[chain]
  }
  intercept <- do.call(rbind, lapply(chains, `[[`, "intercept_draws"))
  colnames(delta) <- annotation_names
  keys <- apply(delta, 1L, paste0, collapse = "")
  levels <- rownames(model_states(ncol(delta)))
  model_probability <- as.numeric(table(factor(keys, levels = levels))) / nrow(delta)
  names(model_probability) <- levels
  alpha_mean <- apply(alpha, c(2L, 3L), mean)
  alpha_mean_given_inclusion <- alpha_mean
  for (j in seq_len(ncol(delta))) {
    included <- delta[, j] == 1L
    alpha_mean_given_inclusion[j, ] <- apply(
      alpha[included, j, , drop = FALSE], 3L, mean
    )
  }
  list(
    delta = delta,
    alpha = alpha,
    intercept = intercept,
    model_probability = model_probability,
    annotation_pip = colMeans(delta),
    alpha_mean = alpha_mean,
    alpha_mean_given_inclusion = alpha_mean_given_inclusion,
    intercept_mean = colMeans(intercept),
    q_mean = Reduce(`+`, lapply(chains, `[[`, "q_mean")) / length(chains),
    pi_mean = Reduce(`+`, lapply(chains, `[[`, "pi_mean")) / length(chains)
  )
}

make_fixture <- function(seed = 20260811L) {
  set.seed(seed)
  observations <- 180L
  annotation <- cbind(
    enriched_binary = as.numeric(seq_len(observations) %% 5L == 0L),
    continuous_signal = as.numeric(scale(sin(seq(0, 4 * pi, length.out = observations)))),
    null_annotation = as.numeric(scale(cos(seq(0, 7 * pi, length.out = observations))))
  )
  eligible <- list(
    seq_len(observations),
    seq.int(1L, observations, by = 2L),
    seq.int(1L, observations, by = 4L)
  )
  true_intercept <- c(-0.8, -0.25, 0.1)
  true_alpha <- rbind(
    enriched_binary = c(0.4, 0.3, 0.2),
    continuous_signal = c(0.22, 0.16, 0.1),
    null_annotation = c(0, 0, 0)
  )
  z <- lapply(seq_along(eligible), function(stick) {
    rows <- eligible[[stick]]
    true_intercept[stick] + drop(annotation[rows, ] %*% true_alpha[, stick]) +
      rnorm(length(rows), 0, 1)
  })
  list(
    annotation = annotation, eligible = eligible, z = z,
    pi_a = 0.35, tau2 = c(0.8, 0.8, 0.8),
    true_delta = c(1L, 1L, 0L), true_alpha = true_alpha,
    true_intercept = true_intercept
  )
}

switch_counts <- function(delta_draws) {
  if (nrow(delta_draws) < 2L) return(rep(0L, ncol(delta_draws)))
  colSums(delta_draws[-1L, , drop = FALSE] !=
            delta_draws[-nrow(delta_draws), , drop = FALSE])
}

run_qualification <- function() {
  settings <- list(
    fixture_seed = 20260811L,
    chain_seeds = 20260901:20260904,
    iterations = 12000L,
    burn = 2000L,
    pi_a = 0.35,
    tau2 = c(0.8, 0.8, 0.8),
    pip_tolerance = 0.02,
    total_variation_tolerance = 0.03,
    model_probability_tolerance = 0.02,
    alpha_mean_tolerance = 0.035,
    q_pi_tolerance = 0.012
  )
  fixture <- make_fixture(settings$fixture_seed)
  exact <- exact_reference(
    fixture$z, fixture$annotation, fixture$eligible,
    fixture$pi_a, fixture$tau2
  )
  initial <- list(c(0L, 0L, 0L), c(1L, 1L, 1L),
                  c(1L, 0L, 1L), c(0L, 1L, 0L))
  started <- proc.time()[["elapsed"]]
  chains <- lapply(seq_along(settings$chain_seeds), function(chain) {
    set.seed(settings$chain_seeds[chain])
    run_chain(
      fixture$z, fixture$annotation, fixture$eligible,
      fixture$pi_a, fixture$tau2,
      settings$iterations, settings$burn, initial[[chain]]
    )
  })
  runtime <- proc.time()[["elapsed"]] - started
  pooled <- pool_chains(chains, colnames(fixture$annotation))

  metrics <- list(
    max_annotation_pip_error = max(abs(
      pooled$annotation_pip - exact$annotation_pip
    )),
    model_total_variation = 0.5 * sum(abs(
      pooled$model_probability - exact$model_probability
    )),
    max_model_probability_error = max(abs(
      pooled$model_probability - exact$model_probability
    )),
    max_alpha_mean_error = max(abs(pooled$alpha_mean - exact$alpha_mean)),
    max_alpha_given_inclusion_mean_error = max(abs(
      pooled$alpha_mean_given_inclusion - exact$alpha_mean_given_inclusion
    )),
    max_intercept_mean_error = max(abs(
      pooled$intercept_mean - exact$intercept_mean
    )),
    max_q_mean_error = max(abs(pooled$q_mean - exact$q_mean)),
    max_pi_mean_error = max(abs(pooled$pi_mean - exact$pi_mean)),
    chain_annotation_pip = do.call(rbind, lapply(
      chains, function(x) colMeans(x$delta_draws)
    )),
    switching_counts = do.call(rbind, lapply(
      chains, function(x) switch_counts(x$delta_draws)
    )),
    runtime_seconds = runtime
  )
  colnames(metrics$chain_annotation_pip) <- colnames(fixture$annotation)
  colnames(metrics$switching_counts) <- colnames(fixture$annotation)

  # Independent direct check of the rank-one Gaussian Bayes factor.
  rows <- fixture$eligible[[2L]]
  x <- fixture$annotation[rows, "continuous_signal"]
  residual <- fixture$z[[2L]] - fixture$true_intercept[2L] -
    fixture$annotation[rows, "enriched_binary"] * 0.17
  metrics$log_bf_crosscheck_error <- abs(
    log_bf_formula(x, residual, fixture$tau2[2L]) -
      log_bf_direct(x, residual, fixture$tau2[2L])
  )

  # Zero-information annotation.
  zero_fixture <- fixture
  zero_fixture$annotation[, "null_annotation"] <- 0
  zero_exact <- exact_reference(
    zero_fixture$z, zero_fixture$annotation, zero_fixture$eligible,
    zero_fixture$pi_a, zero_fixture$tau2
  )
  set.seed(20260920L)
  zero_chain <- run_chain(
    zero_fixture$z, zero_fixture$annotation, zero_fixture$eligible,
    zero_fixture$pi_a, zero_fixture$tau2,
    12000L, 2000L, c(0L, 0L, 0L)
  )
  metrics$zero_exact_pip <- zero_exact$annotation_pip["null_annotation"]
  metrics$zero_mcmc_pip <- mean(zero_chain$delta_draws[, 3L])

  # Column permutation.
  permutation <- c(3L, 1L, 2L)
  permutation_exact <- exact_reference(
    fixture$z, fixture$annotation[, permutation, drop = FALSE],
    fixture$eligible, fixture$pi_a, fixture$tau2
  )
  set.seed(20260930L)
  permutation_chain <- run_chain(
    fixture$z, fixture$annotation[, permutation, drop = FALSE],
    fixture$eligible, fixture$pi_a, fixture$tau2,
    12000L, 2000L, c(1L, 0L, 1L)
  )
  metrics$permutation_exact_error <- max(abs(
    permutation_exact$annotation_pip - exact$annotation_pip[permutation]
  ))
  metrics$permutation_mcmc_error <- max(abs(
    colMeans(permutation_chain$delta_draws) -
      permutation_exact$annotation_pip
  ))

  # Duplicate-column exchange symmetry.
  duplicate <- cbind(
    copy_a = fixture$annotation[, "continuous_signal"],
    copy_b = fixture$annotation[, "continuous_signal"],
    null_annotation = fixture$annotation[, "null_annotation"]
  )
  duplicate_exact <- exact_reference(
    fixture$z, duplicate, fixture$eligible, fixture$pi_a, fixture$tau2
  )
  set.seed(20260940L)
  duplicate_chain <- run_chain(
    fixture$z, duplicate, fixture$eligible, fixture$pi_a, fixture$tau2,
    12000L, 2000L, c(0L, 1L, 0L)
  )
  duplicate_mcmc_pip <- colMeans(duplicate_chain$delta_draws)
  metrics$duplicate_exact_symmetry_error <- abs(
    duplicate_exact$annotation_pip[1L] - duplicate_exact$annotation_pip[2L]
  )
  metrics$duplicate_mcmc_symmetry_error <- abs(
    duplicate_mcmc_pip[1L] - duplicate_mcmc_pip[2L]
  )

  qualification <- c(
    model_probability_normalized = abs(sum(exact$model_probability) - 1) < 1e-14,
    bayes_factor_crosscheck = metrics$log_bf_crosscheck_error <= 1e-10,
    annotation_pip = metrics$max_annotation_pip_error <= settings$pip_tolerance,
    model_total_variation = metrics$model_total_variation <=
      settings$total_variation_tolerance,
    model_probability = metrics$max_model_probability_error <=
      settings$model_probability_tolerance,
    alpha_mean = metrics$max_alpha_mean_error <= settings$alpha_mean_tolerance,
    alpha_given_inclusion = metrics$max_alpha_given_inclusion_mean_error <=
      settings$alpha_mean_tolerance,
    intercept_mean = metrics$max_intercept_mean_error <=
      settings$alpha_mean_tolerance,
    q_mean = metrics$max_q_mean_error <= settings$q_pi_tolerance,
    pi_mean = metrics$max_pi_mean_error <= settings$q_pi_tolerance,
    zero_exact = abs(metrics$zero_exact_pip - fixture$pi_a) <= 1e-12,
    zero_mcmc = abs(metrics$zero_mcmc_pip - fixture$pi_a) <= 0.02,
    permutation_exact = metrics$permutation_exact_error <= 1e-12,
    permutation_mcmc = metrics$permutation_mcmc_error <= 0.02,
    duplicate_exact = metrics$duplicate_exact_symmetry_error <= 1e-12,
    duplicate_mcmc = metrics$duplicate_mcmc_symmetry_error <= 0.035
  )

  list(
    exact_model_prob = exact$model_probability,
    exact_annotation_pip = exact$annotation_pip,
    mcmc_model_prob = pooled$model_probability,
    mcmc_annotation_pip = pooled$annotation_pip,
    exact_alpha_mean = exact$alpha_mean,
    exact_alpha_mean_given_inclusion = exact$alpha_mean_given_inclusion,
    exact_intercept_mean = exact$intercept_mean,
    mcmc_alpha_mean = pooled$alpha_mean,
    mcmc_alpha_mean_given_inclusion = pooled$alpha_mean_given_inclusion,
    mcmc_intercept_mean = pooled$intercept_mean,
    exact_q_mean = exact$q_mean,
    mcmc_q_mean = pooled$q_mean,
    exact_pi_mean = exact$pi_mean,
    mcmc_pi_mean = pooled$pi_mean,
    delta_draws = pooled$delta,
    alpha_draws = pooled$alpha,
    settings = settings,
    seed = settings$fixture_seed,
    fixture = list(
      observations = nrow(fixture$annotation),
      sticks = length(fixture$z),
      annotations = ncol(fixture$annotation),
      eligible_counts = lengths(fixture$eligible),
      annotation_names = colnames(fixture$annotation),
      true_delta = fixture$true_delta,
      true_alpha = fixture$true_alpha,
      true_intercept = fixture$true_intercept
    ),
    metrics = metrics,
    qualification = qualification,
    passed = all(qualification)
  )
}

print_summary <- function(result) {
  cat("SBayesRC-S Phase-1 exact R reference\n")
  cat("=====================================\n")
  cat("Fixed-z target; shared delta; fixed pi_A and tau2; flat intercept.\n\n")
  cat("Exact model probabilities:\n")
  print(round(result$exact_model_prob, 6))
  cat("\nExact and MCMC annotation PIPs:\n")
  print(round(rbind(
    exact = result$exact_annotation_pip,
    mcmc = result$mcmc_annotation_pip
  ), 6))
  cat("\nQualification metrics:\n")
  printable <- result$metrics[vapply(result$metrics, is.numeric, logical(1L)) &
                                lengths(result$metrics) == 1L]
  print(unlist(printable))
  cat("\nChain-specific annotation PIPs:\n")
  print(round(result$metrics$chain_annotation_pip, 4))
  cat("\nAnnotation switching counts:\n")
  print(result$metrics$switching_counts)
  cat("\nGates:\n")
  print(result$qualification)
  cat("\nOverall:", if (result$passed) "PASS" else "FAIL", "\n")
}

result <- run_qualification()
print_summary(result)

output_directory <- file.path("results", "local", "sbayesrc_s_reference")
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
compact <- result
compact$delta_draws <- NULL
compact$alpha_draws <- NULL
saveRDS(compact, file.path(output_directory, "phase1_qualification_summary.rds"))
write.csv(
  data.frame(
    annotation = names(result$exact_annotation_pip),
    exact_pip = result$exact_annotation_pip,
    mcmc_pip = result$mcmc_annotation_pip,
    absolute_error = abs(result$mcmc_annotation_pip - result$exact_annotation_pip)
  ),
  file.path(output_directory, "annotation_pip_summary.csv"),
  row.names = FALSE
)

if (!result$passed) stop("SBayesRC-S Phase-1 reference qualification failed")
