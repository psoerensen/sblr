#!/usr/bin/env Rscript

# Estimator-only feasibility on five evenly spaced frozen blocks from the
# 76-block large Study 06 design. No MCMC is run.

sblr_root <- normalizePath(".", winslash = "/", mustWork = TRUE)
output_root <- file.path(sblr_root, "results", "local",
                         "sbayesrc_particle_marginal_alpha")
cache <- readRDS(file.path(output_root,
  "study06_large_representative_blocks.rds"))
bundle <- readRDS(file.path(sblr_root, "..", "sblrbench", "results", "local",
  "06_annotation_models", "large_feasibility", "prepared_bundle.rds"))
pkgload::load_all(sblr_root, quiet = TRUE)

evaluate_blocks <- function(alpha, particles, auxiliary = NULL) {
  if (is.null(auxiliary)) auxiliary <- lapply(cache$factor, function(Q)
    sblr:::.sbayesrc_particle_auxiliary(particles, ncol(Q)))
  values <- lapply(seq_along(cache$factor), function(index) {
    probability <- sbayesrc_marker_pi(cache$annotations[[index]], alpha,
                                      cache$gamma)
    sblr:::.sbayesrc_sis_block_likelihood(
      cache$factor[[index]], as.numeric(cache$transformed_score[[index]]),
      probability, cache$gamma, cache$vb, cache$block_ve[[index]], particles,
      auxiliary[[index]], retain_paths = FALSE)
  })
  list(log = vapply(values, `[[`, numeric(1), "log_likelihood"),
       auxiliary = auxiliary)
}

alpha <- cache$alpha_learned
variance_rows <- list()
for (particles in c(8L, 16L, 32L, 64L)) {
  repetitions <- 10L
  values <- matrix(NA_real_, repetitions, length(cache$factor))
  started <- proc.time()[["elapsed"]]
  for (repetition in seq_len(repetitions)) {
    set.seed(941000L + particles * 100L + repetition)
    values[repetition, ] <- evaluate_blocks(alpha, particles)$log
  }
  seconds <- proc.time()[["elapsed"]] - started
  block_variance <- apply(values, 2L, stats::var)
  variance_rows[[length(variance_rows) + 1L]] <- data.frame(
    particles = particles, repetitions = repetitions, seconds = seconds,
    seconds_per_five_blocks = seconds / repetitions,
    projected_seconds_per_76_blocks = seconds / repetitions * 76 / 5,
    mean_block_variance = mean(block_variance),
    maximum_block_variance = max(block_variance),
    projected_variance_log_lhat_76 = mean(block_variance) * 76,
    projected_sd_log_lhat_76 = sqrt(mean(block_variance) * 76))
}
variance_summary <- do.call(rbind, variance_rows)
utils::write.csv(variance_summary, file.path(output_root,
  "large_representative_estimator_variance.csv"), row.names = FALSE)

registry <- expand.grid(
  proposal = c("stick_all", "stick_nonintercept", "all_sticks"),
  scale = c(0.05, 0.10, 0.25), rho = c(0, 0.99, 0.999),
  stringsAsFactors = FALSE)
proposal_rows <- list()
particles <- 16L
neutral <- make_sbayesrc_alpha_init(
  bundle$annotations, cache$gamma, pi_init = 0.03,
  active_comp_weights = c(0.5, 1 / 3, 1 / 6))$alpha_init
for (row in seq_len(nrow(registry))) {
  entry <- registry[row, ]
  set.seed(942000L + row * 100L)
  direction <- matrix(stats::rnorm(length(alpha)), nrow(alpha), ncol(alpha))
  direction <- direction / sqrt(sum(direction^2))
  mask <- matrix(FALSE, nrow(alpha), ncol(alpha))
  if (entry$proposal == "stick_all") mask[, 1L] <- TRUE
  if (entry$proposal == "stick_nonintercept") mask[-1L, 1L] <- TRUE
  if (entry$proposal == "all_sticks") mask[,] <- TRUE
  proposed <- alpha
  proposed[mask] <- proposed[mask] + entry$scale * direction[mask]
  repetitions <- 4L
  difference <- numeric(repetitions)
  started <- proc.time()[["elapsed"]]
  for (repetition in seq_len(repetitions)) {
    set.seed(942000L + row * 100L + repetition)
    auxiliary <- lapply(cache$factor, function(Q)
      sblr:::.sbayesrc_particle_auxiliary(particles, ncol(Q)))
    proposed_auxiliary <- if (entry$rho == 0) {
      lapply(cache$factor, function(Q)
        sblr:::.sbayesrc_particle_auxiliary(particles, ncol(Q)))
    } else lapply(auxiliary, sblr:::.sbayesrc_correlate_auxiliary,
                  rho = entry$rho)
    current <- evaluate_blocks(alpha, particles, auxiliary)
    candidate <- evaluate_blocks(proposed, particles, proposed_auxiliary)
    difference[[repetition]] <- sum(candidate$log - current$log)
  }
  seconds <- proc.time()[["elapsed"]] - started
  prior_ratio <- sblr:::.sbayesrc_alpha_log_prior(
    proposed, neutral[1L, ], rep(1, ncol(alpha)), cache$sigma_sq_alpha) -
    sblr:::.sbayesrc_alpha_log_prior(
      alpha, neutral[1L, ], rep(1, ncol(alpha)), cache$sigma_sq_alpha)
  projected_ratio <- difference * 76 / length(cache$factor) + prior_ratio
  probability_current <- sbayesrc_marker_pi(bundle$annotations, alpha, cache$gamma)
  probability_proposed <- sbayesrc_marker_pi(bundle$annotations, proposed, cache$gamma)
  proposal_rows[[row]] <- data.frame(
    proposal = entry$proposal, scale = entry$scale, rho = entry$rho,
    particles = particles, repetitions = repetitions, seconds = seconds,
    projected_seconds_per_76_block_pair = seconds / repetitions * 76 / 5,
    alpha_jump = sqrt(sum((proposed - alpha)^2)),
    expected_active_change = sum(1 - probability_proposed[, 1L]) -
      sum(1 - probability_current[, 1L]),
    sd_five_block_log_difference = stats::sd(difference),
    projected_sd_log_difference = stats::sd(difference) * sqrt(76 / 5),
    projected_median_log_mh = stats::median(projected_ratio),
    projected_mean_acceptance = mean(pmin(1, exp(projected_ratio))))
}
proposal_summary <- do.call(rbind, proposal_rows)
utils::write.csv(proposal_summary, file.path(output_root,
  "large_representative_alpha_proposals.csv"), row.names = FALSE)
saveRDS(list(variance = variance_summary, proposals = proposal_summary),
        file.path(output_root, "large_representative_feasibility.rds"))
print(variance_summary)
print(proposal_summary)
