# Development/reference implementation.
# SBayesRC-S Phase 3D proper-intercept qualification.
# Not a production sampler. Not part of the supported public API.

source(file.path("tests", "testthat", "helper-sbayesrc-s-reference.R"))
source(file.path("tests", "testthat", "helper-sbayesrc-s-probit-reference.R"))
source(file.path("tests", "testthat", "helper-sbayesrc-s-hierarchical-reference.R"))

pool_observed <- function(chains, names) .sbs3_pool(chains, names)
run_hierarchy <- function(fixture, seeds, iterations, burn, initials,
                          method = "primary") {
  lapply(seq_along(seeds), function(i) {
    set.seed(seeds[i])
    .sbs3_run_chain(
      fixture$annotation, fixture$eligible, iterations, burn, initials[[i]],
      fixture$pi_a %||% 0.35,
      fixture$tau2 %||% fixture$tau_true,
      outcome = fixture$outcome, learn_pi = TRUE, learn_tau = TRUE,
      a_pi = 1, b_pi = 1, a_tau = 3, b_tau = 1.6,
      method = method, intercept_prior = fixture$intercept_prior
    )
  })
}

`%||%` <- function(x, y) if (is.null(x)) y else x
started <- proc.time()[["elapsed"]]
chosen_prior <- .sbs_intercept_prior(3L, rep(0.25, 4L), sd = 1)
candidate_prior <- list(
  zero_sd1 = .sbs_intercept_prior(3L, rep(0.25, 4L), mean = 0, sd = 1),
  mixture_sd05 = .sbs_intercept_prior(3L, rep(0.25, 4L), sd = 0.5),
  mixture_sd1 = chosen_prior,
  mixture_sd2 = .sbs_intercept_prior(3L, rep(0.25, 4L), sd = 2)
)
prior_predictive <- do.call(rbind, lapply(names(candidate_prior), function(name) {
  cbind(candidate = name, stick = seq_len(3L),
        .sbs_intercept_prior_predictive(candidate_prior[[name]]))
}))

fixed <- .sbs_fixture()
fixed$intercept_prior <- chosen_prior
exact <- .sbs_exact_posterior(
  fixed$z, fixed$annotation, fixed$eligible, fixed$pi_a, fixed$tau2,
  chosen_prior
)
initials <- list(rep(0L, 3L), rep(1L, 3L), c(1L, 0L, 1L), c(0L, 1L, 0L))
fixed_chains <- lapply(seq_along(initials), function(i) {
  set.seed(20271020L + i)
  .sbs_mcmc_chain(
    fixed$z, fixed$annotation, fixed$eligible, fixed$pi_a, fixed$tau2,
    6000L, 1000L, initials[[i]], chosen_prior
  )
})
fixed_pool <- .sbs_pool_chains(fixed_chains, colnames(fixed$annotation))
fixed_metrics <- c(
  max_pip_error = max(abs(fixed_pool$annotation_pip - exact$annotation_pip)),
  model_tv = 0.5 * sum(abs(fixed_pool$model_probability -
                             exact$model_probability)),
  max_alpha_error = max(abs(fixed_pool$alpha_mean - exact$alpha_mean)),
  max_intercept_error = max(abs(fixed_pool$intercept_mean -
                                  exact$intercept_mean)),
  max_q_error = max(abs(fixed_pool$q_mean - exact$q_mean)),
  max_pi_error = max(abs(fixed_pool$pi_mean - exact$pi_mean))
)

wide_prior <- .sbs_intercept_prior(3L, rep(0.25, 4L), sd = 100)
wide_exact <- .sbs_exact_posterior(
  fixed$z, fixed$annotation, fixed$eligible, fixed$pi_a, fixed$tau2,
  wide_prior
)
nonempty_prior_comparison <- list(
  chosen_vs_wide_max_pip = max(abs(exact$annotation_pip -
                                      wide_exact$annotation_pip)),
  chosen_vs_wide_max_q = max(abs(exact$q_mean - wide_exact$q_mean)),
  chosen_annotation_pip = exact$annotation_pip,
  wide_annotation_pip = wide_exact$annotation_pip
)

observed <- .sbs2_fixture(observations = 300L, seed = 20271040L)
observed$intercept_prior <- chosen_prior
primary <- run_hierarchy(observed, 20271041:20271044, 4500L, 750L, initials)
direct <- run_hierarchy(observed, 20271051:20271054, 4500L, 750L, initials,
                        method = "direct")
primary_pool <- pool_observed(primary, colnames(observed$annotation))
direct_pool <- pool_observed(direct, colnames(observed$annotation))
observed_metrics <- c(
  max_pip_difference = max(abs(primary_pool$annotation_pip -
                                 direct_pool$annotation_pip)),
  max_q_difference = max(abs(primary_pool$q_mean - direct_pool$q_mean)),
  max_pi_difference = max(abs(primary_pool$pi_mean - direct_pool$pi_mean)),
  max_chain_pip_range = max(apply(primary_pool$chain_pip, 2L,
                                  function(x) diff(range(x))))
)

moderate <- .sbs3_moderate_fixture(
  observations = 240L, annotation_count = 12L, seed = 20271060L
)
moderate$pi_a <- 0.35
moderate$intercept_prior <- chosen_prior
moderate_chains <- run_hierarchy(
  moderate, 20271061:20271062, 3500L, 500L,
  list(rep(0L, 12L), rep(1L, 12L))
)
moderate_pool <- pool_observed(moderate_chains, colnames(moderate$annotation))

empty_annotation <- observed$annotation
empty_outcome <- list(rep(0L, nrow(empty_annotation)), integer(), integer())
empty_eligible <- list(seq_len(nrow(empty_annotation)), integer(), integer())
empty_fixture <- observed
empty_fixture$outcome <- empty_outcome
empty_fixture$eligible <- empty_eligible
empty_chains <- run_hierarchy(
  empty_fixture, 20271071:20271072, 2000L, 300L,
  list(rep(0L, 3L), rep(1L, 3L))
)
empty_pool <- pool_observed(empty_chains, colnames(empty_annotation))
empty_metrics <- c(
  all_finite = all(is.finite(c(empty_pool$alpha, empty_pool$intercept,
                               empty_pool$pi_a, empty_pool$tau2,
                               empty_pool$q_mean, empty_pool$pi_mean))),
  max_pi_row_error = max(abs(rowSums(empty_pool$pi_mean) - 1)),
  later_intercept_mean_error = max(abs(
    colMeans(empty_pool$intercept)[2:3] - chosen_prior$mean[2:3]
  )),
  later_intercept_variance_error = max(abs(
    apply(empty_pool$intercept[, 2:3, drop = FALSE], 2L, stats::var) -
      chosen_prior$variance[2:3]
  ))
)

repeated <- lapply(seq_len(5L), function(replicate) {
  fixture <- .sbs3_moderate_fixture(
    observations = 180L, annotation_count = 12L,
    seed = 20271080L + replicate, hierarchical = TRUE
  )
  fixture$pi_a <- 0.35
  fixture$intercept_prior <- chosen_prior
  chains <- run_hierarchy(
    fixture, 20271100L + 2L * replicate + 0:1, 1800L, 300L,
    list(rep(0L, 12L), rep(1L, 12L))
  )
  pooled <- pool_observed(chains, colnames(fixture$annotation))
  data.frame(
    replicate = replicate, truth = fixture$true_delta,
    pip = pooled$annotation_pip,
    pi_a = mean(pooled$pi_a), tau2 = mean(pooled$tau2),
    max_chain_pip_range = max(apply(pooled$chain_pip, 2L,
                                    function(x) diff(range(x))))
  )
})
repeated <- do.call(rbind, repeated)

qualification <- list(
  prior = chosen_prior,
  prior_predictive = prior_predictive,
  fixed_metrics = fixed_metrics,
  nonempty_prior_comparison = nonempty_prior_comparison,
  observed_metrics = observed_metrics,
  observed_annotation_pip = primary_pool$annotation_pip,
  moderate_annotation_pip = moderate_pool$annotation_pip,
  moderate_chain_pip = moderate_pool$chain_pip,
  empty_metrics = empty_metrics,
  repeated = repeated,
  runtime_seconds = proc.time()[["elapsed"]] - started
)
qualification$pass <- isTRUE(
  fixed_metrics[["max_pip_error"]] <= 0.02 &&
    fixed_metrics[["model_tv"]] <= 0.03 &&
    fixed_metrics[["max_alpha_error"]] <= 0.04 &&
    fixed_metrics[["max_intercept_error"]] <= 0.04 &&
    fixed_metrics[["max_q_error"]] <= 0.015 &&
    fixed_metrics[["max_pi_error"]] <= 0.015 &&
    observed_metrics[["max_pip_difference"]] <= 0.03 &&
    observed_metrics[["max_q_difference"]] <= 0.02 &&
    observed_metrics[["max_pi_difference"]] <= 0.02 &&
    as.logical(empty_metrics[["all_finite"]]) &&
    empty_metrics[["max_pi_row_error"]] <= 1e-12 &&
    empty_metrics[["later_intercept_mean_error"]] <= 0.08 &&
    empty_metrics[["later_intercept_variance_error"]] <= 0.12 &&
    all(is.finite(c(moderate_pool$annotation_pip, moderate_pool$pi_a,
                    moderate_pool$tau2, repeated$pip)))
)

output_dir <- file.path(
  "results", "local", "sbayesrc_s_reference", "phase3D"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(qualification, file.path(output_dir, "phase3d_qualification.rds"))
write.csv(prior_predictive, file.path(output_dir, "prior_predictive.csv"),
          row.names = FALSE)
write.csv(repeated, file.path(output_dir, "repeated_hierarchy.csv"),
          row.names = FALSE)
print(qualification)
if (!qualification$pass) stop("Phase 3D qualification gates did not pass")
