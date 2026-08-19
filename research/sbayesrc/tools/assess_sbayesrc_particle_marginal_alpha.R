#!/usr/bin/env Rscript

# Development-only estimator and correlated pseudo-marginal feasibility audit.
# This script does not run an MCMC chain or modify a production sampler.

sblr_root <- normalizePath(".", winslash = "/", mustWork = TRUE)
output_root <- file.path(sblr_root, "results", "local",
                         "sbayesrc_particle_marginal_alpha")
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
pkgload::load_all(sblr_root, quiet = TRUE)

log_mean_exp <- function(x) {
  maximum <- max(x)
  maximum + log(mean(exp(x - maximum)))
}

summarise_estimates <- function(log_estimate, exact_log = NA_real_) {
  shifted <- exp(log_estimate - max(log_estimate))
  cv <- stats::sd(shifted) / mean(shifted)
  ratio <- if (is.finite(exact_log)) exp(log_mean_exp(log_estimate) - exact_log) else NA_real_
  data.frame(mean_lhat_over_exact = ratio,
             variance_log_lhat = stats::var(log_estimate),
             sd_log_lhat = stats::sd(log_estimate),
             cv_lhat = cv)
}

# Tiny exact likelihood-scale gate.
set.seed(912001)
Q_tiny <- matrix(stats::rnorm(15), 3L, 5L)
w_tiny <- stats::rnorm(3L)
A_tiny <- cbind(1, c(0, 1, 0, 1, 1))
alpha_tiny <- matrix(c(-0.8, 0.7), 2L, 1L)
gamma_tiny <- c(0, 0.15)
pi_tiny <- sbayesrc_marker_pi(A_tiny, alpha_tiny, gamma_tiny)
exact_tiny <- sblr:::.sbayesrc_exact_block_allocation(
  Q_tiny, w_tiny, pi_tiny, gamma_tiny, vb = 0.25, ve = 1.2)
tiny_rows <- list()
for (particles in c(4L, 8L, 16L, 32L, 64L)) {
  repetitions <- 400L
  values <- vector("list", repetitions)
  started <- proc.time()[["elapsed"]]
  for (repetition in seq_len(repetitions)) {
    set.seed(912000L + particles * 1000L + repetition)
    values[[repetition]] <- sblr:::.sbayesrc_sis_block_likelihood(
      Q_tiny, w_tiny, pi_tiny, gamma_tiny, 0.25, 1.2, particles,
      retain_paths = FALSE)
  }
  seconds <- proc.time()[["elapsed"]] - started
  log_estimate <- vapply(values, `[[`, numeric(1), "log_likelihood")
  summary <- summarise_estimates(log_estimate, exact_tiny$log_normalizer)
  tiny_rows[[length(tiny_rows) + 1L]] <- cbind(
    particles = particles, repetitions = repetitions,
    seconds = seconds, seconds_per_estimate = seconds / repetitions,
    summary,
    mean_final_ess = mean(vapply(values, `[[`, numeric(1), "final_ess")),
    mean_minimum_prefix_ess = mean(vapply(
      values, `[[`, numeric(1), "minimum_prefix_ess")))
}
tiny_summary <- do.call(rbind, tiny_rows)
utils::write.csv(tiny_summary, file.path(output_root, "tiny_unbiasedness.csv"),
                 row.names = FALSE)

# Retained Study 06 large-block state. The factor was created and validated by
# the preceding BLOCK-MIX-R4 work; the sibling repository is only read here.
block_path <- file.path(sblr_root, "results", "local",
                        "sbayesrc_block_particle", "frozen_block_1_contract.rds")
bundle_path <- file.path(sblr_root, "..", "sblrbench", "results", "local",
  "06_annotation_models", "large_feasibility", "prepared_bundle.rds")
captured_path <- file.path(sblr_root, "..", "sblrbench", "results", "local",
  "06_annotation_models", "large_feasibility", "continuation",
  "b0_iteration0_no_updateE.rds")
stopifnot(file.exists(block_path), file.exists(bundle_path), file.exists(captured_path))
block <- readRDS(block_path)
bundle <- readRDS(bundle_path)
captured <- readRDS(captured_path)
stopifnot(ncol(block$Q) == 500L, nrow(bundle$annotations) == 37991L)

gamma <- as.numeric(captured$mixture_var)
vb <- as.numeric(captured$vbs[1L, 1L])
ve <- block$yy / (block$n - 1)
alpha_truth <- bundle$calibrated$alpha
alpha_neutral <- make_sbayesrc_alpha_init(
  A = bundle$annotations, gamma = gamma,
  pi_init = 1 - bundle$calibrated$expected_pi[[1L]],
  active_comp_weights = bundle$calibrated$expected_pi[-1L] /
    sum(bundle$calibrated$expected_pi[-1L]))$alpha_init

prepare_block <- function(marker_count, alpha) {
  selected <- seq_len(marker_count)
  excluded <- setdiff(seq_len(500L), selected)
  beta <- as.numeric(captured$b[seq_len(500L), 1L])
  adjusted_w <- block$w
  if (length(excluded))
    adjusted_w <- adjusted_w - as.numeric(
      block$Q[, excluded, drop = FALSE] %*% beta[excluded])
  list(Q = block$Q[, selected, drop = FALSE], w = adjusted_w,
       probability = sbayesrc_marker_pi(
         bundle$annotations[selected, , drop = FALSE], alpha, gamma))
}

# Independent-estimator variance at neutral and truth-like annotation states.
independent_rows <- list()
for (marker_count in c(100L, 500L)) {
  for (state in c("neutral", "truth")) {
    alpha <- if (state == "neutral") alpha_neutral else alpha_truth
    fixture <- prepare_block(marker_count, alpha)
    for (particles in c(8L, 16L, 32L, 64L, 128L)) {
      repetitions <- if (marker_count == 100L) 20L else 10L
      estimates <- vector("list", repetitions)
      started <- proc.time()[["elapsed"]]
      for (repetition in seq_len(repetitions)) {
        set.seed(913000L + marker_count * 100L + particles * 10L +
                   match(state, c("neutral", "truth")) * 100000L + repetition)
        estimates[[repetition]] <- sblr:::.sbayesrc_sis_block_likelihood(
          fixture$Q, fixture$w, fixture$probability, gamma, vb, ve, particles,
          retain_paths = FALSE)
      }
      seconds <- proc.time()[["elapsed"]] - started
      log_estimate <- vapply(estimates, `[[`, numeric(1), "log_likelihood")
      independent_rows[[length(independent_rows) + 1L]] <- cbind(
        markers = marker_count, state = state, particles = particles,
        repetitions = repetitions, seconds = seconds,
        seconds_per_estimate = seconds / repetitions,
        summarise_estimates(log_estimate),
        mean_final_ess = mean(vapply(estimates, `[[`, numeric(1), "final_ess")),
        mean_minimum_prefix_ess = mean(vapply(
          estimates, `[[`, numeric(1), "minimum_prefix_ess")))
    }
  }
}
independent_summary <- do.call(rbind, independent_rows)
utils::write.csv(independent_summary,
                 file.path(output_root, "independent_block_variance.csv"),
                 row.names = FALSE)

# Directly measure the correlated log-ratio noise. Proposals are fixed before
# inspection: one complete stick, one non-intercept stick, then all sticks.
proposal_registry <- expand.grid(
  proposal = c("stick_all", "stick_nonintercept", "all_sticks"),
  scale = c(0.02, 0.05, 0.10, 0.25),
  rho = c(0, 0.9, 0.99, 0.999), stringsAsFactors = FALSE)
correlated_rows <- list()
marker_count <- 500L
particles <- 16L
A <- bundle$annotations[seq_len(marker_count), , drop = FALSE]
base_fixture <- prepare_block(marker_count, alpha_truth)
for (row in seq_len(nrow(proposal_registry))) {
  entry <- proposal_registry[row, ]
  repetitions <- 4L
  log_current <- log_proposed <- log_difference <- numeric(repetitions)
  alpha_jump <- numeric(repetitions)
  set.seed(914000L + row * 100L)
  direction <- stats::rnorm(length(alpha_truth))
  direction <- direction / sqrt(sum(direction^2))
  mask <- matrix(FALSE, nrow(alpha_truth), ncol(alpha_truth))
  if (entry$proposal == "stick_all") mask[, 1L] <- TRUE
  if (entry$proposal == "stick_nonintercept") mask[-1L, 1L] <- TRUE
  if (entry$proposal == "all_sticks") mask[,] <- TRUE
  proposed <- alpha_truth
  proposed[mask] <- proposed[mask] + entry$scale * direction[mask]
  started <- proc.time()[["elapsed"]]
  for (repetition in seq_len(repetitions)) {
    set.seed(914000L + row * 100L + repetition)
    auxiliary <- sblr:::.sbayesrc_particle_auxiliary(particles, marker_count)
    correlated <- if (entry$rho == 0) {
      sblr:::.sbayesrc_particle_auxiliary(particles, marker_count)
    } else {
      sblr:::.sbayesrc_correlate_auxiliary(auxiliary, entry$rho)
    }
    current <- sblr:::.sbayesrc_sis_block_likelihood(
      base_fixture$Q, base_fixture$w, base_fixture$probability,
      gamma, vb, ve, particles, auxiliary, retain_paths = FALSE)
    proposed_probability <- sbayesrc_marker_pi(A, proposed, gamma)
    candidate <- sblr:::.sbayesrc_sis_block_likelihood(
      base_fixture$Q, base_fixture$w, proposed_probability,
      gamma, vb, ve, particles, correlated, retain_paths = FALSE)
    log_current[[repetition]] <- current$log_likelihood
    log_proposed[[repetition]] <- candidate$log_likelihood
    log_difference[[repetition]] <- candidate$log_likelihood - current$log_likelihood
    alpha_jump[[repetition]] <- sqrt(sum((proposed - alpha_truth)^2))
  }
  seconds <- proc.time()[["elapsed"]] - started
  correlated_rows[[row]] <- data.frame(
    proposal = entry$proposal, scale = entry$scale, rho = entry$rho,
    markers = marker_count, particles = particles, repetitions = repetitions,
    seconds = seconds, mean_alpha_jump = mean(alpha_jump),
    correlation_log_lhat = stats::cor(log_current, log_proposed),
    sd_log_likelihood_difference = stats::sd(log_difference),
    median_abs_log_likelihood_difference = stats::median(abs(log_difference)),
    median_log_likelihood_difference = stats::median(log_difference),
    mean_noise_only_acceptance_proxy = mean(pmin(1, exp(log_difference))))
}
correlated_summary <- do.call(rbind, correlated_rows)
utils::write.csv(correlated_summary,
  file.path(output_root, "correlated_log_ratio_screen.csv"), row.names = FALSE)

saveRDS(list(tiny = tiny_summary, independent = independent_summary,
             correlated = correlated_summary),
        file.path(output_root, "particle_marginal_feasibility.rds"))
print(tiny_summary)
print(independent_summary)
print(correlated_summary)
