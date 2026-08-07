#!/usr/bin/env Rscript

# Short, development-only selected-path correlated PMMH screen. The non-alpha
# global state is fixed at the retained Study 06 reference values. This is not a
# production sampler or a qualification fit.

sblr_root <- normalizePath(".", winslash = "/", mustWork = TRUE)
output_root <- file.path(sblr_root, "results", "local",
                         "sbayesrc_particle_marginal_alpha")
cache <- readRDS(file.path(output_root, "study06_small_block_contract.rds"))
pkgload::load_all(sblr_root, quiet = TRUE)

block_probability <- function(alpha) {
  probability <- sbayesrc_marker_pi(cache$annotations, alpha, cache$gamma)
  lapply(seq_along(cache$factor), function(block) {
    index <- ((block - 1L) * 100L + 1L):(block * 100L)
    probability[index, , drop = FALSE]
  })
}

evaluate_global <- function(alpha, particles, auxiliary) {
  probability <- block_probability(alpha)
  values <- lapply(seq_along(cache$factor), function(block) {
    sblr:::.sbayesrc_sis_block_likelihood(
      cache$factor[[block]], as.numeric(cache$transformed_score[[block]]),
      probability[[block]], cache$gamma, cache$vb, cache$block_ve[[block]],
      particles, auxiliary[[block]], retain_paths = TRUE)
  })
  list(log_likelihood = sum(vapply(values, `[[`, numeric(1), "log_likelihood")),
       values = values, auxiliary = auxiliary)
}

neutral <- make_sbayesrc_alpha_init(
  cache$annotations, cache$gamma, pi_init = 0.12,
  active_comp_weights = c(0.5, 1 / 3, 1 / 6))$alpha_init
intercept_mean <- neutral[1L, ]
intercept_sd <- rep(1, ncol(neutral))
particles <- 8L
rho <- 0.99
proposal_scale <- 0.10
iterations <- 300L
burn <- 100L
seeds <- c(931121L, 931222L, 931323L, 931424L)

run_chain <- function(seed, chain) {
  set.seed(seed)
  alpha <- cache$alpha_pooled + matrix(stats::rnorm(length(cache$alpha_pooled),
                                                     0, 0.03),
                                       nrow(cache$alpha_pooled))
  auxiliary <- lapply(cache$factor, function(Q)
    sblr:::.sbayesrc_particle_auxiliary(particles, ncol(Q)))
  current <- evaluate_global(alpha, particles, auxiliary)
  alpha_trace <- array(NA_real_, c(iterations, nrow(alpha), ncol(alpha)))
  occupancy <- matrix(NA_integer_, iterations, length(cache$gamma))
  expected_active <- numeric(iterations)
  accepted <- logical(iterations)
  log_ratio <- numeric(iterations)
  alpha_jump <- active_jump <- numeric(iterations)
  started <- proc.time()[["elapsed"]]
  current_component <- unlist(lapply(current$values, `[[`, "component"),
                              use.names = FALSE)
  for (iteration in seq_len(iterations)) {
    stick <- (iteration - 1L) %% ncol(alpha) + 1L
    proposed <- alpha
    proposed[, stick] <- proposed[, stick] +
      stats::rnorm(nrow(alpha), 0, proposal_scale)
    proposed_auxiliary <- lapply(auxiliary,
      sblr:::.sbayesrc_correlate_auxiliary, rho = rho)
    candidate <- evaluate_global(proposed, particles, proposed_auxiliary)
    ratio <- sblr:::.sbayesrc_particle_marginal_log_ratio(
      alpha, proposed, current$log_likelihood, candidate$log_likelihood,
      intercept_mean, intercept_sd, cache$sigma_sq_alpha)
    proposed_component <- unlist(lapply(candidate$values, `[[`, "component"),
                                 use.names = FALSE)
    if (log(stats::runif(1L)) < min(0, ratio)) {
      alpha_jump[[iteration]] <- sqrt(sum((proposed - alpha)^2))
      active_jump[[iteration]] <- sum(proposed_component > 0L) -
        sum(current_component > 0L)
      alpha <- proposed
      auxiliary <- proposed_auxiliary
      current <- candidate
      current_component <- proposed_component
      accepted[[iteration]] <- TRUE
    }
    log_ratio[[iteration]] <- ratio
    alpha_trace[iteration, , ] <- alpha
    occupancy[iteration, ] <- tabulate(current_component + 1L,
                                        nbins = length(cache$gamma))
    probability <- block_probability(alpha)
    expected_active[[iteration]] <- sum(vapply(probability, function(p)
      sum(1 - p[, 1L]), numeric(1L)))
  }
  seconds <- proc.time()[["elapsed"]] - started
  list(alpha = alpha_trace, occupancy = occupancy,
       expected_active = expected_active, accepted = accepted,
       log_ratio = log_ratio, alpha_jump = alpha_jump,
       active_jump = active_jump, seconds = seconds, chain = chain,
       seed = seed, final_component = current_component,
       final_beta = unlist(lapply(current$values, `[[`, "beta"), use.names = FALSE))
}

chains <- lapply(seq_along(seeds), function(chain) run_chain(seeds[[chain]], chain))
saveRDS(chains, file.path(output_root, "small_reference_pmmh_chains.rds"),
        compress = FALSE)

summary <- do.call(rbind, lapply(chains, function(x) data.frame(
  chain = x$chain, seed = x$seed, iterations = iterations, burn = burn,
  retained = iterations - burn, seconds = x$seconds,
  acceptance = mean(x$accepted), accepted_moves = sum(x$accepted),
  median_accepted_alpha_jump = stats::median(x$alpha_jump[x$accepted]),
  maximum_accepted_alpha_jump = max(x$alpha_jump),
  median_absolute_active_jump = stats::median(abs(x$active_jump[x$accepted])),
  maximum_absolute_active_jump = max(abs(x$active_jump)),
  mean_active = mean(rowSums(x$occupancy[-seq_len(burn), -1L, drop = FALSE])),
  sd_active = stats::sd(rowSums(x$occupancy[-seq_len(burn), -1L, drop = FALSE])),
  mean_expected_active = mean(x$expected_active[-seq_len(burn)]),
  stringsAsFactors = FALSE)))
utils::write.csv(summary, file.path(output_root, "small_reference_pmmh_summary.csv"),
                 row.names = FALSE)

alpha_draws <- do.call(rbind, lapply(chains, function(x) {
  index <- (burn + 1L):iterations
  dimensions <- dim(x$alpha)
  data.frame(chain = x$chain, iteration = index,
    value = as.vector(x$alpha[index, , ]),
    coefficient = rep(rep(seq_len(dimensions[[2L]]), each = length(index)),
                      dimensions[[3L]]),
    stick = rep(seq_len(dimensions[[3L]]),
                each = length(index) * dimensions[[2L]]))
}))
utils::write.csv(alpha_draws, file.path(output_root,
  "small_reference_pmmh_alpha.csv"), row.names = FALSE)
print(summary)
