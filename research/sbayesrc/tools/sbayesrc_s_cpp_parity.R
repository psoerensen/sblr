# Development/reference implementation.
# SBayesRC-S Phase 4A C++/R hierarchy qualification.
# Not a production sampler. Not part of the supported public API.

load_reference_environment <- function() {
  devtools::load_all(quiet = TRUE)
  source(file.path("tests", "testthat", "helper-sbayesrc-s-reference.R"),
         local = .GlobalEnv)
  source(file.path("tests", "testthat", "helper-sbayesrc-alpha-reference.R"),
         local = .GlobalEnv)
  source(file.path("tests", "testthat", "helper-sbayesrc-s-probit-reference.R"),
         local = .GlobalEnv)
  source(file.path("tests", "testthat", "helper-sbayesrc-s-hierarchical-reference.R"),
         local = .GlobalEnv)
  source(file.path("tests", "testthat", "helper-sbayesrc-s-cpp-reference.R"),
         local = .GlobalEnv)
}

run_pair <- function(fixture, initials, seeds, iterations, burn,
                     a_pi = 1, b_pi = 9, fixed_delta = NULL) {
  initial_tau2 <- if (!is.null(fixture$tau2)) fixture$tau2 else fixture$tau_true
  cpp <- lapply(seq_along(initials), function(i) {
    .sbs4a_cpp_chain(
      fixture, seeds[i], iterations, burn, initials[[i]],
      if (is.null(fixed_delta)) integer() else fixed_delta,
      a_pi = a_pi, b_pi = b_pi
    )
  })
  r <- lapply(seq_along(initials), function(i) {
    set.seed(seeds[i] + 10000L)
    .sbs3_run_chain(
      fixture$annotation, fixture$eligible, iterations, burn, initials[[i]],
      0.35, initial_tau2, outcome = fixture$outcome,
      learn_pi = TRUE, learn_tau = TRUE,
      a_pi = a_pi, b_pi = b_pi, a_tau = 3, b_tau = 1.6,
      fixed_delta = fixed_delta,
      intercept_prior = fixture$intercept_prior
    )
  })
  list(cpp = cpp, r = r, cpp_summary = .sbs4a_cpp_summary(cpp),
       r_summary = .sbs3_pool(r, colnames(fixture$annotation)))
}

comparison <- function(pair) {
  cs <- pair$cpp_summary
  rs <- pair$r_summary
  alpha_sd <- apply(rs$alpha, c(2L, 3L), stats::sd)
  c(
    max_pip_error = max(abs(cs$annotation_pip - rs$annotation_pip)),
    max_q_error = max(abs(cs$q_mean - rs$q_mean)),
    max_pi_error = max(abs(cs$component_probability_mean - rs$pi_mean)),
    max_normalized_alpha_error = max(
      abs(cs$alpha_mean - rs$alpha_mean) / pmax(alpha_sd, 0.05)
    ),
    pi_A_mean_error = abs(cs$pi_a_mean - mean(rs$pi_a)),
    max_tau2_mean_error = max(abs(cs$tau2_mean - colMeans(rs$tau2))),
    max_cpp_chain_pip_range = max(apply(cs$chain_pip, 2L, diff_range))
  )
}

diff_range <- function(x) diff(range(x))

load_reference_environment()
started <- proc.time()[["elapsed"]]

small <- .sbs2_fixture(observations = 300L, seed = 20270850L)
small_initial <- list(c(0L, 0L, 0L), c(1L, 1L, 1L),
                      c(1L, 0L, 1L), c(0L, 1L, 0L))
small_pair <- run_pair(small, small_initial, 20270851:20270854,
                       8000L, 1500L, a_pi = 1, b_pi = 1)
small_metrics <- comparison(small_pair)

moderate <- .sbs3_moderate_fixture(
  observations = 260L, annotation_count = 12L, seed = 20270860L
)
moderate_pair <- run_pair(
  moderate,
  list(rep(0L, 12L), rep(1L, 12L), rep(c(1L, 0L), 6L),
       rep(c(0L, 1L), 6L)),
  20270861:20270864, 10000L, 1500L, a_pi = 1, b_pi = 9
)
moderate_metrics <- comparison(moderate_pair)

# Structural bridges and invariances use the same observed outcomes, priors,
# and hyperparameter schedule on both sides.
included_pair <- run_pair(
  small, list(rep(1L, 3L), rep(1L, 3L)), 20270871:20270872,
  6000L, 1000L, a_pi = 1, b_pi = 1, fixed_delta = rep(1L, 3L)
)
included_metrics <- comparison(included_pair)

permutation <- c(3L, 1L, 2L)
permuted <- small
permuted$annotation <- small$annotation[, permutation, drop = FALSE]
permutation_cpp <- .sbs4a_cpp_summary(lapply(seq_along(small_initial), function(i) {
  .sbs4a_cpp_chain(permuted, 20270880L + i, 7000L, 1250L,
                   small_initial[[i]][permutation], a_pi = 1, b_pi = 1)
}))
permutation_error <- max(abs(
  permutation_cpp$annotation_pip[order(permutation)] -
    small_pair$cpp_summary$annotation_pip
))

duplicate <- small
duplicate$annotation[, 3L] <- duplicate$annotation[, 1L]
duplicate_cpp <- .sbs4a_cpp_summary(list(
  .sbs4a_cpp_chain(duplicate, 20270891L, 9000L, 1500L, rep(0L, 3L),
                   a_pi = 1, b_pi = 1),
  .sbs4a_cpp_chain(duplicate, 20270892L, 9000L, 1500L, rep(1L, 3L),
                   a_pi = 1, b_pi = 1)
))
duplicate_pip_difference <- abs(
  duplicate_cpp$annotation_pip[1L] - duplicate_cpp$annotation_pip[3L]
)

qualification <- list(
  small = small_metrics,
  moderate = moderate_metrics,
  all_included = included_metrics,
  permutation_pip_error = permutation_error,
  duplicate_pip_difference = duplicate_pip_difference,
  small_cpp_pip = small_pair$cpp_summary$annotation_pip,
  small_r_pip = small_pair$r_summary$annotation_pip,
  moderate_cpp_pip = moderate_pair$cpp_summary$annotation_pip,
  moderate_r_pip = moderate_pair$r_summary$annotation_pip,
  runtime_seconds = proc.time()[["elapsed"]] - started
)
qualification$pass <- isTRUE(
  small_metrics[["max_pip_error"]] <= 0.03 &&
    small_metrics[["max_q_error"]] <= 0.02 &&
    small_metrics[["max_pi_error"]] <= 0.02 &&
    small_metrics[["max_normalized_alpha_error"]] <= 0.15 &&
    small_metrics[["pi_A_mean_error"]] <= 0.05 &&
    small_metrics[["max_tau2_mean_error"]] <= 0.10 &&
    moderate_metrics[["max_pip_error"]] <= 0.03 &&
    moderate_metrics[["max_q_error"]] <= 0.02 &&
    moderate_metrics[["max_pi_error"]] <= 0.02 &&
    moderate_metrics[["max_normalized_alpha_error"]] <= 0.15 &&
    moderate_metrics[["pi_A_mean_error"]] <= 0.05 &&
    moderate_metrics[["max_tau2_mean_error"]] <= 0.10 &&
    included_metrics[["max_pip_error"]] <= 1e-12 &&
    included_metrics[["max_q_error"]] <= 0.02 &&
    included_metrics[["max_pi_error"]] <= 0.02 &&
    permutation_error <= 0.03 && duplicate_pip_difference <= 0.05
)

output_dir <- file.path("results", "local", "sbayesrc_s_reference", "phase4A2")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(qualification, file.path(output_dir, "phase4a_cpp_r_parity.rds"))
write.csv(
  data.frame(metric = names(unlist(qualification[1:5])),
             value = unname(unlist(qualification[1:5]))),
  file.path(output_dir, "phase4a_cpp_r_parity.csv"), row.names = FALSE
)
print(qualification)
if (!qualification$pass) stop("Phase 4A qualification gates did not pass")
