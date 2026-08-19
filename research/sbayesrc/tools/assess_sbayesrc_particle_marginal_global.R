#!/usr/bin/env Rscript

# Estimator-only global alpha feasibility screen for the immutable 1,500-marker
# Study 06 block model. No MCMC transition is run.

sblr_root <- normalizePath(".", winslash = "/", mustWork = TRUE)
output_root <- file.path(sblr_root, "results", "local",
                         "sbayesrc_particle_marginal_alpha")
cache <- readRDS(file.path(output_root, "study06_small_block_contract.rds"))
stopifnot(identical(cache$specification_hash,
  "241c15afab8fefc571e38e625130de6e4ab58b958c30cff369e26998ee30fa56"),
  identical(cache$truth_hash,
  "169d52bef390022a9106d7e61b200493869b40cbb76f1c8d911eebbc80fea1eb"))
pkgload::load_all(sblr_root, quiet = TRUE)

block_probability <- function(alpha) {
  probability <- sbayesrc_marker_pi(cache$annotations, alpha, cache$gamma)
  lapply(seq_along(cache$factor), function(block) {
    index <- ((block - 1L) * 100L + 1L):(block * 100L)
    probability[index, , drop = FALSE]
  })
}

evaluate_global <- function(alpha, particles, auxiliary = NULL,
                            retain_paths = FALSE) {
  probability <- block_probability(alpha)
  if (is.null(auxiliary))
    auxiliary <- lapply(cache$factor, function(Q)
      sblr:::.sbayesrc_particle_auxiliary(particles, ncol(Q)))
  values <- lapply(seq_along(cache$factor), function(block) {
    sblr:::.sbayesrc_sis_block_likelihood(
      cache$factor[[block]], as.numeric(cache$transformed_score[[block]]),
      probability[[block]], cache$gamma, cache$vb, cache$block_ve[[block]],
      particles, auxiliary[[block]], retain_paths = retain_paths)
  })
  list(log_likelihood = sum(vapply(values, `[[`, numeric(1), "log_likelihood")),
       block_log_likelihood = vapply(values, `[[`, numeric(1), "log_likelihood"),
       auxiliary = auxiliary, values = values)
}

alpha <- cache$alpha_pooled
neutral <- make_sbayesrc_alpha_init(
  cache$annotations, cache$gamma, pi_init = 0.12,
  active_comp_weights = c(0.5, 1 / 3, 1 / 6))$alpha_init
intercept_mean <- neutral[1L, ]
intercept_sd <- rep(1, ncol(alpha))

# Independent genome-wide estimator variance and its blockwise accumulation.
independent_rows <- list()
block_rows <- list()
for (particles in c(8L, 16L, 32L, 64L)) {
  repetitions <- 12L
  global <- numeric(repetitions)
  per_block <- matrix(NA_real_, repetitions, length(cache$factor))
  started <- proc.time()[["elapsed"]]
  for (repetition in seq_len(repetitions)) {
    set.seed(921000L + particles * 100L + repetition)
    value <- evaluate_global(alpha, particles)
    global[[repetition]] <- value$log_likelihood
    per_block[repetition, ] <- value$block_log_likelihood
  }
  seconds <- proc.time()[["elapsed"]] - started
  independent_rows[[length(independent_rows) + 1L]] <- data.frame(
    particles = particles, repetitions = repetitions, seconds = seconds,
    seconds_per_global_estimate = seconds / repetitions,
    variance_log_lhat = stats::var(global), sd_log_lhat = stats::sd(global),
    sum_block_variances = sum(apply(per_block, 2L, stats::var)),
    maximum_block_variance = max(apply(per_block, 2L, stats::var)))
  block_rows[[length(block_rows) + 1L]] <- data.frame(
    particles = particles, block = seq_len(ncol(per_block)),
    variance_log_lhat = apply(per_block, 2L, stats::var),
    sd_log_lhat = apply(per_block, 2L, stats::sd))
}
independent <- do.call(rbind, independent_rows)
block_variance <- do.call(rbind, block_rows)
utils::write.csv(independent, file.path(output_root,
  "small_global_independent_variance.csv"), row.names = FALSE)
utils::write.csv(block_variance, file.path(output_root,
  "small_block_variance_contributions.csv"), row.names = FALSE)

# Fixed proposal grid. Each row uses a single predeclared direction across
# repeated auxiliary draws so variance estimates isolate particle randomness.
registry <- expand.grid(
  proposal = c("stick_all", "stick_nonintercept", "all_sticks"),
  scale = c(0.02, 0.05, 0.10, 0.25), rho = c(0, 0.99, 0.999),
  stringsAsFactors = FALSE)
proposal_rows <- list()
particles <- 16L
for (row in seq_len(nrow(registry))) {
  entry <- registry[row, ]
  set.seed(922000L + row * 100L)
  direction <- matrix(stats::rnorm(length(alpha)), nrow(alpha), ncol(alpha))
  direction <- direction / sqrt(sum(direction^2))
  mask <- matrix(FALSE, nrow(alpha), ncol(alpha))
  if (entry$proposal == "stick_all") mask[, 1L] <- TRUE
  if (entry$proposal == "stick_nonintercept") mask[-1L, 1L] <- TRUE
  if (entry$proposal == "all_sticks") mask[,] <- TRUE
  proposed <- alpha
  proposed[mask] <- proposed[mask] + entry$scale * direction[mask]
  repetitions <- 4L
  current_log <- proposed_log <- log_ratio <- numeric(repetitions)
  started <- proc.time()[["elapsed"]]
  for (repetition in seq_len(repetitions)) {
    set.seed(922000L + row * 100L + repetition)
    auxiliary <- lapply(cache$factor, function(Q)
      sblr:::.sbayesrc_particle_auxiliary(particles, ncol(Q)))
    proposed_auxiliary <- if (entry$rho == 0) {
      lapply(cache$factor, function(Q)
        sblr:::.sbayesrc_particle_auxiliary(particles, ncol(Q)))
    } else {
      lapply(auxiliary, sblr:::.sbayesrc_correlate_auxiliary, rho = entry$rho)
    }
    current <- evaluate_global(alpha, particles, auxiliary)
    candidate <- evaluate_global(proposed, particles, proposed_auxiliary)
    current_log[[repetition]] <- current$log_likelihood
    proposed_log[[repetition]] <- candidate$log_likelihood
    log_ratio[[repetition]] <- sblr:::.sbayesrc_particle_marginal_log_ratio(
      alpha, proposed, current$log_likelihood, candidate$log_likelihood,
      intercept_mean, intercept_sd, cache$sigma_sq_alpha)
  }
  seconds <- proc.time()[["elapsed"]] - started
  proposal_rows[[row]] <- data.frame(
    proposal = entry$proposal, scale = entry$scale, rho = entry$rho,
    particles = particles, repetitions = repetitions, seconds = seconds,
    alpha_jump = sqrt(sum((proposed - alpha)^2)),
    expected_active_change = sum(vapply(block_probability(proposed), function(p)
      sum(1 - p[, 1L]), numeric(1L))) -
      sum(vapply(block_probability(alpha), function(p)
        sum(1 - p[, 1L]), numeric(1L))),
    correlation_log_lhat = stats::cor(current_log, proposed_log),
    sd_log_likelihood_difference = stats::sd(proposed_log - current_log),
    median_log_mh_ratio = stats::median(log_ratio),
    minimum_log_mh_ratio = min(log_ratio), maximum_log_mh_ratio = max(log_ratio),
    mean_acceptance = mean(pmin(1, exp(log_ratio))))
}
proposal_summary <- do.call(rbind, proposal_rows)
utils::write.csv(proposal_summary, file.path(output_root,
  "small_global_alpha_proposal_screen.csv"), row.names = FALSE)

# Ordinary-chain alpha step sizes provide a scale reference only; no fit is run.
fit <- readRDS(file.path(sblr_root, "results", "local",
  "study06_kernel_composition_audit", "block_eigen_PX_screen_fit_result.rds"))
chain_steps <- lapply(seq_along(fit$result$native_fit$chains), function(chain) {
  trace <- fit$result$native_fit$chains[[chain]]$alpha
  if (length(dim(trace)) != 3L) return(NULL)
  steps <- vapply(2:dim(trace)[1L], function(iteration) sqrt(sum(
    (trace[iteration, , ] - trace[iteration - 1L, , ])^2)), numeric(1L))
  data.frame(chain = chain, median_step = stats::median(steps),
             q95_step = stats::quantile(steps, 0.95), maximum_step = max(steps))
})
chain_steps <- do.call(rbind, chain_steps)
if (!is.null(chain_steps)) utils::write.csv(chain_steps, file.path(output_root,
  "ordinary_alpha_step_scale.csv"), row.names = FALSE)

saveRDS(list(independent = independent, block_variance = block_variance,
             proposals = proposal_summary, ordinary_steps = chain_steps),
        file.path(output_root, "small_global_feasibility.rds"))
print(independent)
print(proposal_summary)
print(chain_steps)
