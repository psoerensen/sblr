# SBayesRC-S Phase-3 hierarchical hyperparameter reference
#
# Development/reference implementation.
# SBayesRC-S Phase 3 hierarchical hyperparameter sampler.
# Not a production sampler.
# Not part of the supported public API.

phase1_helper <- file.path("tests", "testthat", "helper-sbayesrc-s-reference.R")
phase2_helper <- file.path(
  "tests", "testthat", "helper-sbayesrc-s-probit-reference.R"
)
phase3_helper <- file.path(
  "tests", "testthat", "helper-sbayesrc-s-hierarchical-reference.R"
)
if (!all(file.exists(c(phase1_helper, phase2_helper, phase3_helper)))) {
  stop("Run this reference from the sblr repository root")
}
sys.source(phase1_helper, envir = environment())
sys.source(phase2_helper, envir = environment())
sys.source(phase3_helper, envir = environment())

initial_states <- function(annotation_count) {
  list(
    rep(0L, annotation_count),
    rep(1L, annotation_count),
    as.integer(seq_len(annotation_count) %% 2L == 1L),
    as.integer(seq_len(annotation_count) %% 3L == 0L)
  )
}

run_hierarchical_chains <- function(fixture, seeds, iterations, burn,
                                    learn_pi, collapsed_pi, learn_tau,
                                    a_pi, b_pi, a_tau, b_tau,
                                    fixed_z = FALSE) {
  initial <- initial_states(ncol(fixture$annotation))
  lapply(seq_along(seeds), function(chain) {
    set.seed(seeds[chain])
    .sbs3_run_chain(
      fixture$annotation, fixture$eligible, iterations, burn,
      initial[[chain]], pi_a = 0.35, tau2 = rep(0.8, 3L),
      outcome = if (fixed_z) NULL else fixture$outcome,
      fixed_z = if (fixed_z) fixture$latent_true else NULL,
      learn_pi = learn_pi, collapsed_pi = collapsed_pi,
      a_pi = a_pi, b_pi = b_pi,
      learn_tau = learn_tau, a_tau = a_tau, b_tau = b_tau,
      initial_alpha = matrix(
        stats::rnorm(ncol(fixture$annotation) * 3L, 0, 0.2),
        ncol(fixture$annotation), 3L
      ),
      initial_intercept = stats::rnorm(3L, 0, 0.3)
    )
  })
}

chain_rhat <- function(chains, field) {
  .sbs2_rhat(lapply(chains, function(x) as.matrix(x[[field]])))
}

chain_ess <- function(chains, field) {
  if (!requireNamespace("coda", quietly = TRUE)) return(NA_real_)
  values <- coda::mcmc.list(lapply(chains, function(x) {
    coda::mcmc(as.matrix(x[[field]]))
  }))
  as.numeric(coda::effectiveSize(values))
}

summarize_hierarchy <- function(chains, annotation_names) {
  pooled <- .sbs3_pool(chains, annotation_names)
  pooled$pi_rhat <- chain_rhat(chains, "pi_a_draws")
  pooled$tau_rhat <- chain_rhat(chains, "tau2_draws")
  pooled$included_rhat <- chain_rhat(chains, "included_draws")
  pooled$pi_ess <- chain_ess(chains, "pi_a_draws")
  pooled$tau_ess <- chain_ess(chains, "tau2_draws")
  pooled$included_ess <- chain_ess(chains, "included_draws")
  pooled$pip_range <- apply(pooled$chain_pip, 2L, function(x) diff(range(x)))
  pooled
}

settings <- list(
  phase2_head = "a51bf7c95ef79c8e3ed649f935329f6999776e09",
  a_tau = 3,
  b_tau = 1.6,
  tau_prior_mean = 0.8,
  phase3a_pip_gate = 0.025,
  phase3a_route_gate = 0.035,
  phase3b_pip_gate = 0.035,
  phase3b_tau_mean_gate = 0.10,
  phase3c_route_pip_gate = 0.07,
  phase3c_route_tau_gate = 0.18,
  phase3c_rhat_gate = 1.10,
  phase3c_pip_range_gate = 0.12
)

cat("SBayesRC-S Phase 3A: learned pi_A\n")
cat("===================================\n")
fixture <- .sbs2_fixture(observations = 240L, seed = 20270410L)
exact_a <- .sbs3_exact_pi_posterior(
  fixture$latent_true, fixture$annotation, fixture$eligible,
  fixture$tau2, 1, 1
)
started <- proc.time()[["elapsed"]]
explicit_a <- run_hierarchical_chains(
  fixture, 20270411:20270414, 5500L, 1000L,
  TRUE, FALSE, FALSE, 1, 1, settings$a_tau, settings$b_tau, TRUE
)
collapsed_a <- run_hierarchical_chains(
  fixture, 20270421:20270424, 5500L, 1000L,
  TRUE, TRUE, FALSE, 1, 1, settings$a_tau, settings$b_tau, TRUE
)
phase3a_runtime <- proc.time()[["elapsed"]] - started
pool_explicit_a <- summarize_hierarchy(explicit_a, colnames(fixture$annotation))
pool_collapsed_a <- summarize_hierarchy(collapsed_a, colnames(fixture$annotation))

observed_explicit_a <- run_hierarchical_chains(
  fixture, 20270431:20270434, 4500L, 750L,
  TRUE, FALSE, FALSE, 1, 1, settings$a_tau, settings$b_tau
)
observed_collapsed_a <- run_hierarchical_chains(
  fixture, 20270441:20270444, 4500L, 750L,
  TRUE, TRUE, FALSE, 1, 1, settings$a_tau, settings$b_tau
)
pool_observed_explicit_a <- summarize_hierarchy(
  observed_explicit_a, colnames(fixture$annotation)
)
pool_observed_collapsed_a <- summarize_hierarchy(
  observed_collapsed_a, colnames(fixture$annotation)
)

sensitivity_a <- lapply(list(uniform = c(1, 1), sparse = c(1, 9)), function(prior) {
  chains <- run_hierarchical_chains(
    fixture, 20270450L + seq_len(4L) + prior[2L], 3500L, 500L,
    TRUE, TRUE, FALSE, prior[1L], prior[2L],
    settings$a_tau, settings$b_tau
  )
  summary <- summarize_hierarchy(chains, colnames(fixture$annotation))
  list(
    annotation_pip = summary$annotation_pip,
    pi_mean = mean(summary$pi_a),
    expected_included = mean(summary$included)
  )
})

phase3a_metrics <- list(
  exact_pip = exact_a$annotation_pip,
  explicit_pip = pool_explicit_a$annotation_pip,
  collapsed_pip = pool_collapsed_a$annotation_pip,
  max_explicit_exact_pip_error = max(abs(
    pool_explicit_a$annotation_pip - exact_a$annotation_pip
  )),
  max_collapsed_exact_pip_error = max(abs(
    pool_collapsed_a$annotation_pip - exact_a$annotation_pip
  )),
  max_explicit_collapsed_pip_difference = max(abs(
    pool_explicit_a$annotation_pip - pool_collapsed_a$annotation_pip
  )),
  observed_max_pip_difference = max(abs(
    pool_observed_explicit_a$annotation_pip -
      pool_observed_collapsed_a$annotation_pip
  )),
  observed_max_q_difference = max(abs(
    pool_observed_explicit_a$q_mean - pool_observed_collapsed_a$q_mean
  )),
  expected_included_identity_error = max(
    abs(mean(pool_explicit_a$included) - sum(pool_explicit_a$annotation_pip)),
    abs(mean(pool_collapsed_a$included) - sum(pool_collapsed_a$annotation_pip))
  ),
  explicit_pi_mean = mean(pool_explicit_a$pi_a),
  collapsed_pi_mean = mean(pool_collapsed_a$pi_a),
  explicit_pi_rhat = pool_explicit_a$pi_rhat,
  collapsed_pi_rhat = pool_collapsed_a$pi_rhat,
  sensitivity = sensitivity_a,
  runtime_seconds = phase3a_runtime
)
phase3a_qualification <- c(
  exact_normalization = abs(sum(exact_a$model_probability) - 1) <= 1e-14,
  explicit_exact = phase3a_metrics$max_explicit_exact_pip_error <=
    settings$phase3a_pip_gate,
  collapsed_exact = phase3a_metrics$max_collapsed_exact_pip_error <=
    settings$phase3a_pip_gate,
  route_agreement = phase3a_metrics$max_explicit_collapsed_pip_difference <=
    settings$phase3a_route_gate,
  observed_route_agreement = phase3a_metrics$observed_max_pip_difference <= 0.05,
  expected_included_identity = phase3a_metrics$expected_included_identity_error <= 1e-13,
  pi_mixing = max(c(pool_explicit_a$pi_rhat, pool_collapsed_a$pi_rhat),
                  na.rm = TRUE) <= 1.05,
  prior_sensitivity = sensitivity_a$sparse$expected_included <
    sensitivity_a$uniform$expected_included
)
if (!all(phase3a_qualification)) {
  print(phase3a_metrics)
  print(phase3a_qualification)
  stop("SBS3A-R2: pi_A hierarchy blocked")
}
cat("SBS3A-R1 PASS\n")

cat("\nSBayesRC-S Phase 3B: learned tau2\n")
cat("==================================\n")
one_annotation <- fixture$annotation[, "continuous_signal", drop = FALSE]
oracle_b <- .sbs3_tau_quadrature(
  fixture$latent_true[1L], one_annotation, fixture$eligible[1L],
  pi_a = 0.35, settings$a_tau, settings$b_tau
)
set.seed(20270501L)
chain_b_oracle <- .sbs3_run_chain(
  one_annotation, fixture$eligible[1L], 22000L, 3000L,
  0L, 0.35, 0.8, fixed_z = fixture$latent_true[1L],
  learn_tau = TRUE, a_tau = settings$a_tau, b_tau = settings$b_tau
)
started <- proc.time()[["elapsed"]]
chains_b <- run_hierarchical_chains(
  fixture, 20270511:20270514, 6500L, 1000L,
  FALSE, FALSE, TRUE, 1, 1, settings$a_tau, settings$b_tau
)
phase3b_runtime <- proc.time()[["elapsed"]] - started
pool_b <- summarize_hierarchy(chains_b, colnames(fixture$annotation))
concentrated_b <- run_hierarchical_chains(
  fixture, 20270521:20270524, 4000L, 750L,
  FALSE, FALSE, TRUE, 1, 1, 10, 7.2
)
pool_concentrated_b <- summarize_hierarchy(
  concentrated_b, colnames(fixture$annotation)
)

# Fixed pi_A, learned tau2 repeated behavior. The generating slopes are fixed,
# not slab draws, so this is an empirical-regularization diagnostic rather than
# a claim of hyperparameter truth recovery.
repeat_b <- vector("list", 10L)
for (replicate in seq_len(10L)) {
  generated_b <- .sbs2_fixture(
    observations = 180L, seed = 20270540L + replicate
  )
  chains_b_rep <- run_hierarchical_chains(
    generated_b, 20270560L + replicate * 3L + seq_len(2L),
    1200L, 250L, FALSE, FALSE, TRUE, 1, 1,
    settings$a_tau, settings$b_tau
  )
  pool_b_rep <- summarize_hierarchy(
    chains_b_rep, colnames(generated_b$annotation)
  )
  informative <- generated_b$true_delta == 1L
  repeat_b[[replicate]] <- data.frame(
    replicate = replicate,
    tau_mean = mean(pool_b_rep$tau2),
    informative_pip = mean(pool_b_rep$annotation_pip[informative]),
    null_pip = mean(pool_b_rep$annotation_pip[!informative]),
    max_tau_rhat = max(pool_b_rep$tau_rhat, na.rm = TRUE)
  )
}
repeat_b <- do.call(rbind, repeat_b)

set.seed(20270531L)
empty_draws <- replicate(30000L, .sbs3_draw_tau2(
  numeric(3L), integer(3L), settings$a_tau, settings$b_tau
))
phase3b_metrics <- list(
  quadrature_pip = oracle_b$pip,
  mcmc_pip = mean(chain_b_oracle$delta_draws),
  quadrature_tau_mean = oracle_b$mean_tau,
  mcmc_tau_mean = mean(chain_b_oracle$tau2_draws),
  observed_tau_mean = colMeans(pool_b$tau2),
  observed_tau_rhat = pool_b$tau_rhat,
  observed_tau_ess = pool_b$tau_ess,
  concentrated_tau_mean = colMeans(pool_concentrated_b$tau2),
  repeated_tau_mean = mean(repeat_b$tau_mean),
  repeated_informative_pip = mean(repeat_b$informative_pip),
  repeated_null_pip = mean(repeat_b$null_pip),
  repeated_max_tau_rhat = max(repeat_b$max_tau_rhat),
  empty_model_prior_mean = mean(empty_draws),
  annotation_pip = pool_b$annotation_pip,
  runtime_seconds = phase3b_runtime
)
phase3b_qualification <- c(
  quadrature_pip = abs(phase3b_metrics$mcmc_pip - oracle_b$pip) <=
    settings$phase3b_pip_gate,
  quadrature_tau = abs(phase3b_metrics$mcmc_tau_mean -
                         oracle_b$mean_tau) <= settings$phase3b_tau_mean_gate,
  empty_prior = abs(phase3b_metrics$empty_model_prior_mean - 0.8) <= 0.025,
  finite_positive = all(is.finite(pool_b$tau2)) && all(pool_b$tau2 > 0),
  tau_mixing = max(pool_b$tau_rhat, na.rm = TRUE) <= 1.05,
  repeated_stability = max(repeat_b$max_tau_rhat) <= 1.15 &&
    mean(repeat_b$informative_pip) > mean(repeat_b$null_pip),
  probability_contract = max(abs(rowSums(pool_b$pi_mean) - 1)) <= 1e-12
)
if (!all(phase3b_qualification)) {
  print(phase3b_metrics)
  print(phase3b_qualification)
  stop("SBS3B-R2: slab-variance hierarchy blocked")
}
cat("SBS3B-R1 PASS\n")

cat("\nSBayesRC-S Phase 3C: joint hierarchy\n")
cat("=====================================\n")
started <- proc.time()[["elapsed"]]
chains_c <- run_hierarchical_chains(
  fixture, 20270601:20270604, 7500L, 1250L,
  TRUE, FALSE, TRUE, 1, 1, settings$a_tau, settings$b_tau
)
chains_c_collapsed <- run_hierarchical_chains(
  fixture, 20270611:20270614, 7500L, 1250L,
  TRUE, TRUE, TRUE, 1, 1, settings$a_tau, settings$b_tau
)
pool_c <- summarize_hierarchy(chains_c, colnames(fixture$annotation))
pool_c_collapsed <- summarize_hierarchy(
  chains_c_collapsed, colnames(fixture$annotation)
)

moderate <- .sbs3_moderate_fixture(
  observations = 260L, annotation_count = 12L, seed = 20270620L
)
moderate_chains <- run_hierarchical_chains(
  moderate, 20270621:20270624, 4200L, 700L,
  TRUE, TRUE, TRUE, 1, 9, settings$a_tau, settings$b_tau
)
pool_moderate <- summarize_hierarchy(
  moderate_chains, colnames(moderate$annotation)
)

delta_moderate <- pool_moderate$delta
joint_proxy <- c(
  signal_only = mean(delta_moderate[, 2L] == 1L & delta_moderate[, 4L] == 0L),
  proxy_only = mean(delta_moderate[, 2L] == 0L & delta_moderate[, 4L] == 1L),
  both = mean(delta_moderate[, 2L] == 1L & delta_moderate[, 4L] == 1L),
  neither = mean(delta_moderate[, 2L] == 0L & delta_moderate[, 4L] == 0L)
)

repeat_records <- vector("list", 20L)
repeat_failures <- 0L
for (replicate in seq_len(20L)) {
  generated <- .sbs3_moderate_fixture(
    observations = 140L, annotation_count = 10L,
    seed = 20270700L + replicate, hierarchical = TRUE,
    pi_true = 0.25, tau_true = rep(0.8, 3L)
  )
  repeat_chains <- run_hierarchical_chains(
    generated, 20271000L + replicate * 10L + seq_len(4L),
    900L, 200L, TRUE, TRUE, TRUE, 1, 3,
    settings$a_tau, settings$b_tau
  )
  repeat_pool <- summarize_hierarchy(
    repeat_chains, colnames(generated$annotation)
  )
  if (max(repeat_pool$pip_range) > 0.30 ||
      max(repeat_pool$tau_rhat, na.rm = TRUE) > 1.25) {
    repeat_failures <- repeat_failures + 1L
  }
  informative <- generated$true_delta == 1L
  alpha_bias <- numeric(0L)
  alpha_covered <- logical(0L)
  for (j in which(informative)) {
    included_draw <- repeat_pool$delta[, j] == 1L
    for (stick in seq_len(3L)) {
      values <- repeat_pool$alpha[included_draw, j, stick]
      if (length(values)) {
        interval <- stats::quantile(values, c(0.025, 0.975))
        alpha_bias <- c(alpha_bias, mean(values) - generated$true_alpha[j, stick])
        alpha_covered <- c(
          alpha_covered,
          interval[1L] <= generated$true_alpha[j, stick] &&
            interval[2L] >= generated$true_alpha[j, stick]
        )
      }
    }
  }
  pi_interval <- stats::quantile(repeat_pool$pi_a, c(0.025, 0.975))
  tau_interval <- apply(
    repeat_pool$tau2, 2L, stats::quantile, probs = c(0.025, 0.975)
  )
  repeat_records[[replicate]] <- data.frame(
    replicate = replicate,
    annotation = colnames(generated$annotation),
    truth_delta = generated$true_delta,
    pip = repeat_pool$annotation_pip,
    pi_true = generated$pi_true,
    pi_mean = mean(repeat_pool$pi_a),
    tau_true = mean(generated$tau_true),
    tau_mean = mean(repeat_pool$tau2),
    included_mean = mean(repeat_pool$included),
    informative_alpha_bias = if (length(alpha_bias)) mean(alpha_bias) else NA_real_,
    informative_alpha_coverage = if (length(alpha_covered)) {
      mean(alpha_covered)
    } else NA_real_,
    pi_covered = pi_interval[1L] <= generated$pi_true &&
      pi_interval[2L] >= generated$pi_true,
    tau_coverage = mean(
      tau_interval[1L, ] <= generated$tau_true &
        tau_interval[2L, ] >= generated$tau_true
    )
  )
}
repeat_records <- do.call(rbind, repeat_records)
phase3c_runtime <- proc.time()[["elapsed"]] - started

correlation_c <- c(
  pi_M = stats::cor(pool_c$pi_a, pool_c$included),
  pi_tau1 = stats::cor(pool_c$pi_a, pool_c$tau2[, 1L]),
  M_tau1 = stats::cor(pool_c$included, pool_c$tau2[, 1L]),
  pi_tau2 = stats::cor(pool_c$pi_a, pool_c$tau2[, 2L]),
  M_tau2 = stats::cor(pool_c$included, pool_c$tau2[, 2L]),
  pi_tau3 = stats::cor(pool_c$pi_a, pool_c$tau2[, 3L]),
  M_tau3 = stats::cor(pool_c$included, pool_c$tau2[, 3L])
)
bfdr <- data.frame(
  threshold = c(0.5, 0.8, 0.9),
  bfdr = vapply(c(0.5, 0.8, 0.9), function(x) {
    .sbs3_bfdr(pool_moderate$annotation_pip, x)
  }, numeric(1L)),
  selected = vapply(c(0.5, 0.8, 0.9), function(x) {
    sum(pool_moderate$annotation_pip >= x)
  }, integer(1L))
)

phase3c_metrics <- list(
  explicit_pip = pool_c$annotation_pip,
  collapsed_pip = pool_c_collapsed$annotation_pip,
  max_route_pip_difference = max(abs(
    pool_c$annotation_pip - pool_c_collapsed$annotation_pip
  )),
  max_route_tau_difference = max(abs(
    colMeans(pool_c$tau2) - colMeans(pool_c_collapsed$tau2)
  )),
  max_route_q_difference = max(abs(pool_c$q_mean - pool_c_collapsed$q_mean)),
  pi_mean = mean(pool_c$pi_a),
  pi_interval = stats::quantile(pool_c$pi_a, c(0.025, 0.5, 0.975)),
  tau_mean = colMeans(pool_c$tau2),
  tau_interval = apply(pool_c$tau2, 2L, stats::quantile,
                       probs = c(0.025, 0.5, 0.975)),
  included_mean = mean(pool_c$included),
  included_distribution = table(pool_c$included) / length(pool_c$included),
  pi_rhat = pool_c$pi_rhat,
  tau_rhat = pool_c$tau_rhat,
  included_rhat = pool_c$included_rhat,
  pi_ess = pool_c$pi_ess,
  tau_ess = pool_c$tau_ess,
  included_ess = pool_c$included_ess,
  pip_range = pool_c$pip_range,
  correlations = correlation_c,
  moderate_pip = pool_moderate$annotation_pip,
  moderate_pi_mean = mean(pool_moderate$pi_a),
  moderate_tau_mean = colMeans(pool_moderate$tau2),
  moderate_included_mean = mean(pool_moderate$included),
  moderate_pip_range = pool_moderate$pip_range,
  moderate_hyper_rhat = c(pool_moderate$pi_rhat, pool_moderate$tau_rhat),
  proxy_joint = joint_proxy,
  bfdr = bfdr,
  repeated_informative_pip_mean = mean(
    repeat_records$pip[repeat_records$truth_delta == 1L]
  ),
  repeated_null_pip_mean = mean(
    repeat_records$pip[repeat_records$truth_delta == 0L]
  ),
  repeated_pi_mean = mean(unique(
    repeat_records[c("replicate", "pi_mean")]$pi_mean
  )),
  repeated_tau_mean = mean(unique(
    repeat_records[c("replicate", "tau_mean")]$tau_mean
  )),
  repeated_alpha_bias = mean(unique(
    repeat_records[c("replicate", "informative_alpha_bias")]
  )$informative_alpha_bias, na.rm = TRUE),
  repeated_alpha_coverage = mean(unique(
    repeat_records[c("replicate", "informative_alpha_coverage")]
  )$informative_alpha_coverage, na.rm = TRUE),
  repeated_pi_coverage = mean(unique(
    repeat_records[c("replicate", "pi_covered")]
  )$pi_covered),
  repeated_tau_coverage = mean(unique(
    repeat_records[c("replicate", "tau_coverage")]
  )$tau_coverage),
  repeated_convergence_failures = repeat_failures,
  runtime_seconds = phase3c_runtime
)

phase3c_qualification <- c(
  route_pip = phase3c_metrics$max_route_pip_difference <=
    settings$phase3c_route_pip_gate,
  route_tau = phase3c_metrics$max_route_tau_difference <=
    settings$phase3c_route_tau_gate,
  route_q = phase3c_metrics$max_route_q_difference <= 0.04,
  primary_hyper_mixing = max(c(pool_c$pi_rhat, pool_c$tau_rhat),
                             na.rm = TRUE) <= settings$phase3c_rhat_gate,
  primary_pip_stability = max(pool_c$pip_range) <=
    settings$phase3c_pip_range_gate,
  moderate_finite = all(is.finite(c(
    pool_moderate$annotation_pip, pool_moderate$pi_a, pool_moderate$tau2
  ))),
  moderate_probability = max(abs(rowSums(pool_moderate$pi_mean) - 1)) <= 1e-12,
  moderate_discrimination = mean(pool_moderate$annotation_pip[1:3]) >
    mean(pool_moderate$annotation_pip[5:12]),
  repeat_discrimination = phase3c_metrics$repeated_informative_pip_mean >
    phase3c_metrics$repeated_null_pip_mean,
  repeat_convergence = repeat_failures <= 4L,
  phase3a = all(phase3a_qualification),
  phase3b = all(phase3b_qualification)
)

result <- list(
  settings = settings,
  phase3a = list(
    metrics = phase3a_metrics,
    qualification = phase3a_qualification,
    decision = "SBS3A-R1"
  ),
  phase3b = list(
    metrics = phase3b_metrics,
    qualification = phase3b_qualification,
    decision = "SBS3B-R1"
  ),
  phase3c = list(
    metrics = phase3c_metrics,
    qualification = phase3c_qualification,
    primary_chain_pip = pool_c$chain_pip,
    primary_switching = pool_c$switching,
    moderate_chain_pip = pool_moderate$chain_pip,
    moderate_switching = pool_moderate$switching,
    repeated = repeat_records,
    decision = if (all(phase3c_qualification)) "SBS3-R1" else "SBS3-R2"
  ),
  passed = all(phase3c_qualification)
)

cat("Phase 3A metrics:\n")
print(phase3a_metrics[!vapply(phase3a_metrics, is.list, logical(1L))])
cat("Phase 3A gates:\n")
print(phase3a_qualification)
cat("Phase 3B metrics:\n")
print(phase3b_metrics)
cat("Phase 3B gates:\n")
print(phase3b_qualification)
cat("Phase 3C primary chain PIPs:\n")
print(round(pool_c$chain_pip, 4))
cat("Phase 3C metrics:\n")
print(phase3c_metrics)
cat("Phase 3C gates:\n")
print(phase3c_qualification)
cat("Overall:", if (result$passed) "SBS3-R1 PASS" else "SBS3-R2 FAIL", "\n")

output_directory <- file.path(
  "results", "local", "sbayesrc_s_reference", "phase3"
)
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
compact <- result
compact$phase3c$repeated <- NULL
saveRDS(compact, file.path(output_directory, "phase3C_joint_summary.rds"))
saveRDS(result$phase3a, file.path(output_directory, "phase3A_pi_summary.rds"))
saveRDS(result$phase3b, file.path(output_directory, "phase3B_tau_summary.rds"))
write.csv(
  data.frame(
    annotation = colnames(fixture$annotation),
    pip = pool_c$annotation_pip,
    chain_range = pool_c$pip_range
  ),
  file.path(output_directory, "phase3_annotation_pips.csv"), row.names = FALSE
)
write.csv(
  repeat_records,
  file.path(output_directory, "phase3_repeated_recovery.csv"), row.names = FALSE
)
write.csv(
  data.frame(
    annotation = colnames(moderate$annotation),
    pip = pool_moderate$annotation_pip,
    truth_delta = moderate$true_delta,
    chain_range = pool_moderate$pip_range
  ),
  file.path(output_directory, "phase3_moderate_annotation_pips.csv"),
  row.names = FALSE
)
write.csv(
  data.frame(pair = names(joint_proxy), probability = joint_proxy),
  file.path(output_directory, "phase3_joint_inclusion.csv"), row.names = FALSE
)
write.csv(bfdr, file.path(output_directory, "phase3_bfdr.csv"), row.names = FALSE)

if (!result$passed) stop("SBS3-R2: joint hierarchy mixing unresolved")
