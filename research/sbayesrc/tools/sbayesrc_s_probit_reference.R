# SBayesRC-S Phase-2 observed-d / Albert-Chib reference
#
# Development/reference implementation.
# Not a production SBayesRC-S sampler.
# Not part of the supported public API.
#
# The script explicitly reuses the committed, test-only Phase-1 mathematics
# and the independent Phase-2 test helpers. Package tests never source tools/.

phase1_helper <- file.path("tests", "testthat", "helper-sbayesrc-s-reference.R")
phase2_helper <- file.path(
  "tests", "testthat", "helper-sbayesrc-s-probit-reference.R"
)
standard_helper <- file.path(
  "tests", "testthat", "helper-sbayesrc-alpha-reference.R"
)
if (!file.exists(phase1_helper) || !file.exists(phase2_helper) ||
    !file.exists(standard_helper)) {
  stop("Run this reference from the sblr repository root")
}
sys.source(phase1_helper, envir = environment())
sys.source(standard_helper, envir = environment())
sys.source(phase2_helper, envir = environment())

pool_phase2 <- function(chains, annotation_names) {
  pooled <- .sbs_pool_chains(chains, annotation_names)
  pooled$chain_pip <- do.call(rbind, lapply(
    chains, function(x) colMeans(x$delta_draws)
  ))
  colnames(pooled$chain_pip) <- annotation_names
  pooled$pip_range <- apply(pooled$chain_pip, 2L, function(x) diff(range(x)))
  pooled$switching <- lapply(chains, function(x) {
    .sbs2_switch_diagnostics(x$delta_draws)
  })
  pooled$alpha_rhat <- matrix(
    .sbs2_rhat(lapply(chains, function(x) {
      matrix(x$alpha_draws, nrow = dim(x$alpha_draws)[1L])
    })),
    nrow = dim(chains[[1L]]$alpha_draws)[2L],
    ncol = dim(chains[[1L]]$alpha_draws)[3L]
  )
  pooled$intercept_rhat <- .sbs2_rhat(lapply(chains, `[[`, "intercept_draws"))
  pooled$alpha_sd <- apply(pooled$alpha, c(2L, 3L), stats::sd)
  pooled$alpha_ess <- if (requireNamespace("coda", quietly = TRUE)) {
    chain_list <- coda::mcmc.list(lapply(chains, function(x) {
      coda::mcmc(matrix(x$alpha_draws, nrow = dim(x$alpha_draws)[1L]))
    }))
    matrix(
      as.numeric(coda::effectiveSize(chain_list)),
      nrow = dim(chains[[1L]]$alpha_draws)[2L],
      ncol = dim(chains[[1L]]$alpha_draws)[3L]
    )
  } else {
    matrix(
      NA_real_, dim(chains[[1L]]$alpha_draws)[2L],
      dim(chains[[1L]]$alpha_draws)[3L]
    )
  }
  pooled$alpha_interval <- apply(
    pooled$alpha, c(2L, 3L), stats::quantile, probs = c(0.025, 0.975)
  )
  pooled$alpha_given_inclusion_sd <- pooled$alpha_sd
  pooled$alpha_given_inclusion_interval <- array(
    NA_real_, c(2L, ncol(pooled$delta), dim(pooled$alpha)[3L])
  )
  for (j in seq_len(ncol(pooled$delta))) {
    included <- pooled$delta[, j] == 1L
    for (stick in seq_len(dim(pooled$alpha)[3L])) {
      values <- pooled$alpha[included, j, stick]
      pooled$alpha_given_inclusion_sd[j, stick] <- stats::sd(values)
      pooled$alpha_given_inclusion_interval[, j, stick] <-
        stats::quantile(values, c(0.025, 0.975))
    }
  }
  pooled
}

run_chains <- function(fixture, seeds, iterations, burn, method,
                       initial, fixed_delta = NULL) {
  lapply(seq_along(seeds), function(chain) {
    set.seed(seeds[chain])
    .sbs2_run_chain(
      fixture$outcome, fixture$annotation, fixture$eligible,
      fixture$pi_a, fixture$tau2, iterations, burn,
      initial_delta = initial[[chain]], method = method,
      fixed_delta = fixed_delta,
      initial_alpha = matrix(
        stats::rnorm(ncol(fixture$annotation) * length(fixture$outcome),
                     0, 0.25),
        ncol(fixture$annotation), length(fixture$outcome)
      ),
      initial_intercept = stats::rnorm(length(fixture$outcome), 0, 0.5)
    )
  })
}

compare_pools <- function(primary, direct) {
  scale <- pmax(0.05, 0.5 * (primary$alpha_sd + direct$alpha_sd))
  list(
    max_pip_difference = max(abs(
      primary$annotation_pip - direct$annotation_pip
    )),
    max_normalized_alpha_mean_difference = max(abs(
      primary$alpha_mean - direct$alpha_mean
    ) / scale),
    max_unconditional_alpha_difference = max(abs(
      primary$alpha_mean - direct$alpha_mean
    )),
    max_conditional_alpha_difference = max(abs(
      primary$alpha_mean_given_inclusion -
        direct$alpha_mean_given_inclusion
    )),
    max_intercept_difference = max(abs(
      primary$intercept_mean - direct$intercept_mean
    )),
    max_q_difference = max(abs(primary$q_mean - direct$q_mean)),
    max_pi_difference = max(abs(primary$pi_mean - direct$pi_mean))
  )
}

run_phase1_bridge <- function(fixture) {
  exact <- .sbs_exact_posterior(
    fixture$latent_true, fixture$annotation, fixture$eligible,
    fixture$pi_a, fixture$tau2
  )
  initial <- list(c(0L, 0L, 0L), c(1L, 1L, 1L),
                  c(1L, 0L, 1L), c(0L, 1L, 0L))
  chains <- lapply(seq_along(initial), function(chain) {
    set.seed(20261100L + chain)
    .sbs_mcmc_chain(
      fixture$latent_true, fixture$annotation, fixture$eligible,
      fixture$pi_a, fixture$tau2, 5000L, 1000L, initial[[chain]]
    )
  })
  sampled <- .sbs_pool_chains(chains, colnames(fixture$annotation))
  list(
    max_pip_error = max(abs(sampled$annotation_pip - exact$annotation_pip)),
    model_tv = 0.5 * sum(abs(
      sampled$model_probability - exact$model_probability
    )),
    max_alpha_error = max(abs(sampled$alpha_mean - exact$alpha_mean)),
    max_q_error = max(abs(sampled$q_mean - exact$q_mean)),
    max_pi_error = max(abs(sampled$pi_mean - exact$pi_mean)),
    passed = max(abs(sampled$annotation_pip - exact$annotation_pip)) <= 0.025 &&
      0.5 * sum(abs(sampled$model_probability - exact$model_probability)) <= 0.035
  )
}

run_all_included_bridge <- function(fixture) {
  set.seed(20261201L)
  selection <- .sbs2_run_chain(
    fixture$outcome, fixture$annotation, fixture$eligible,
    fixture$pi_a, fixture$tau2, 7000L, 1000L, rep(1L, 3L),
    fixed_delta = rep(1L, 3L)
  )
  set.seed(20261202L)
  continuous <- .sbs2_standard_continuous_chain(
    fixture$outcome, fixture$annotation, fixture$eligible,
    fixture$tau2, 7000L, 1000L
  )
  selection_coefficient <- array(NA_real_, c(6000L, 4L, 3L))
  selection_coefficient[, 1L, ] <- selection$intercept_draws
  selection_coefficient[, -1L, ] <- selection$alpha_draws
  selection_mean <- apply(selection_coefficient, c(2L, 3L), mean)
  continuous_mean <- apply(continuous$coefficient_draws, c(2L, 3L), mean)
  selection_sd <- apply(selection_coefficient, c(2L, 3L), stats::sd)
  continuous_sd <- apply(continuous$coefficient_draws, c(2L, 3L), stats::sd)
  list(
    max_alpha_mean_difference = max(abs(selection_mean - continuous_mean)),
    max_alpha_sd_difference = max(abs(selection_sd - continuous_sd)),
    max_q_difference = max(abs(selection$q_mean - continuous$q_mean)),
    max_pi_difference = max(abs(selection$pi_mean - continuous$pi_mean))
  )
}

run_all_excluded_bridge <- function(fixture) {
  set.seed(20261210L)
  result <- .sbs2_run_chain(
    fixture$outcome, fixture$annotation, fixture$eligible,
    fixture$pi_a, fixture$tau2, 4000L, 500L, rep(0L, 3L),
    fixed_delta = rep(0L, 3L)
  )
  list(
    max_absolute_slope = max(abs(result$alpha_draws)),
    max_pi_row_sum_error = max(abs(rowSums(result$pi_mean) - 1)),
    finite_intercepts = all(is.finite(result$intercept_draws))
  )
}

run_zero_column <- function(fixture) {
  fixture$annotation[, 3L] <- 0
  chains <- run_chains(
    fixture, 20261221:20261224, 5000L, 1000L, "primary",
    list(c(0L, 0L, 0L), c(1L, 1L, 1L),
         c(1L, 0L, 1L), c(0L, 1L, 0L))
  )
  pips <- vapply(chains, function(x) mean(x$delta_draws[, 3L]), numeric(1L))
  list(chain_pip = pips, pooled_pip = mean(pips), prior = fixture$pi_a)
}

run_permutation <- function(fixture, primary) {
  permutation <- c(3L, 1L, 2L)
  permuted <- fixture
  permuted$annotation <- fixture$annotation[, permutation, drop = FALSE]
  chains <- run_chains(
    permuted, 20261231:20261234, 7000L, 1000L, "primary",
    list(c(0L, 0L, 0L), c(1L, 1L, 1L),
         c(1L, 0L, 1L), c(0L, 1L, 0L))
  )
  pooled <- pool_phase2(chains, colnames(permuted$annotation))
  list(
    max_pip_difference = max(abs(
      pooled$annotation_pip - primary$annotation_pip[permutation]
    )),
    max_alpha_difference = max(abs(
      pooled$alpha_mean - primary$alpha_mean[permutation, , drop = FALSE]
    ))
  )
}

run_duplicate <- function(fixture) {
  duplicate <- fixture
  duplicate$annotation <- cbind(
    copy_a = fixture$annotation[, "continuous_signal"],
    copy_b = fixture$annotation[, "continuous_signal"],
    null_annotation = fixture$annotation[, "null_annotation"]
  )
  chains <- run_chains(
    duplicate, 20261301:20261304, 7000L, 1000L, "primary",
    list(c(0L, 0L, 0L), c(1L, 1L, 1L),
         c(1L, 0L, 1L), c(0L, 1L, 0L))
  )
  pooled <- pool_phase2(chains, colnames(duplicate$annotation))
  list(
    pooled_pip = pooled$annotation_pip,
    chain_pip = pooled$chain_pip,
    symmetry_difference = abs(
      pooled$annotation_pip[1L] - pooled$annotation_pip[2L]
    ),
    switching = pooled$switching
  )
}

run_signal_case <- function(signal_scale, seed, chain_seed) {
  fixture <- .sbs2_fixture(
    observations = 320L, seed = seed, signal_scale = signal_scale
  )
  set.seed(chain_seed)
  chain <- .sbs2_run_chain(
    fixture$outcome, fixture$annotation, fixture$eligible,
    fixture$pi_a, fixture$tau2, 6000L, 1000L, rep(0L, 3L)
  )
  list(
    pip = colMeans(chain$delta_draws),
    switching = .sbs2_switch_diagnostics(chain$delta_draws)
  )
}

run_prior_only <- function() {
  records <- vector("list", 10L)
  for (replicate in seq_len(10L)) {
    fixture <- .sbs2_fixture(
      observations = 260L, seed = 20261400L + replicate,
      signal = c(0, 0, 0)
    )
    # Use exchangeable, standardized null columns so this diagnostic tests
    # model/column symmetry rather than the scale dependence of a common slab.
    set.seed(20261450L + replicate)
    fixture$annotation <- apply(
      matrix(stats::rnorm(260L * 3L), 260L, 3L), 2L,
      function(x) as.numeric(scale(x))
    )
    colnames(fixture$annotation) <- c(
      "null_a", "null_b", "null_c"
    )
    set.seed(20261500L + replicate)
    chain <- .sbs2_run_chain(
      fixture$outcome, fixture$annotation, fixture$eligible,
      fixture$pi_a, fixture$tau2, 2200L, 400L, rep(0L, 3L)
    )
    records[[replicate]] <- colMeans(chain$delta_draws)
  }
  values <- do.call(rbind, records)
  list(mean_pip = colMeans(values), pip = values,
       annotation_range = diff(range(colMeans(values))))
}

run_repeated_simulation <- function(replicates = 20L) {
  records <- list()
  position <- 0L
  convergence_failures <- 0L
  for (replicate in seq_len(replicates)) {
    fixture <- .sbs2_fixture(
      observations = 300L, seed = 20262000L + replicate
    )
    initial <- list(c(0L, 0L, 0L), c(1L, 1L, 1L),
                    c(1L, 0L, 1L), c(0L, 1L, 0L))
    chains <- run_chains(
      fixture, 20263000L + replicate * 10L + seq_len(4L),
      1200L, 200L, "primary", initial
    )
    chain_pip <- do.call(rbind, lapply(
      chains, function(x) colMeans(x$delta_draws)
    ))
    if (max(apply(chain_pip, 2L, function(x) diff(range(x)))) > 0.15) {
      convergence_failures <- convergence_failures + 1L
    }
    pooled <- .sbs_pool_chains(chains, colnames(fixture$annotation))
    for (j in seq_len(3L)) {
      included <- pooled$delta[, j] == 1L
      for (stick in seq_len(3L)) {
        values <- pooled$alpha[included, j, stick]
        interval <- stats::quantile(values, c(0.025, 0.975))
        position <- position + 1L
        records[[position]] <- data.frame(
          replicate = replicate,
          annotation = colnames(fixture$annotation)[j],
          truth_delta = fixture$true_delta[j],
          stick = stick,
          pip = pooled$annotation_pip[j],
          conditional_alpha_mean = mean(values),
          alpha_truth = fixture$true_alpha[j, stick],
          conditional_bias = mean(values) - fixture$true_alpha[j, stick],
          covered = interval[1L] <= fixture$true_alpha[j, stick] &&
            interval[2L] >= fixture$true_alpha[j, stick]
        )
      }
    }
  }
  records <- do.call(rbind, records)
  annotation_level <- unique(records[c(
    "replicate", "annotation", "truth_delta", "pip"
  )])
  split_pip <- split(annotation_level$pip, annotation_level$truth_delta)
  summaries <- lapply(split_pip, function(x) c(
    mean = mean(x), median = stats::median(x),
    q05 = unname(stats::quantile(x, 0.05)),
    q95 = unname(stats::quantile(x, 0.95)),
    fraction_gt_0_5 = mean(x > 0.5), fraction_gt_0_9 = mean(x > 0.9)
  ))
  list(
    records = records,
    annotation_level = annotation_level,
    pip_summary = summaries,
    informative_conditional_bias = mean(
      records$conditional_bias[records$truth_delta == 1L]
    ),
    informative_coverage = mean(records$covered[records$truth_delta == 1L]),
    null_conditional_bias = mean(
      records$conditional_bias[records$truth_delta == 0L]
    ),
    convergence_failures = convergence_failures
  )
}

settings <- list(
  phase1_head = "d6daa179d28c7307abffe77b49fdbf15ee6b33bc",
  fixture_seed = 20261001L,
  chain_seeds = 20261071:20261074,
  iterations = 9000L,
  burn = 1500L,
  pi_a = 0.35,
  tau2 = rep(0.8, 3L),
  pip_comparison_gate = 0.03,
  normalized_alpha_gate = 0.15,
  q_pi_gate = 0.02,
  rhat_gate = 1.05,
  pip_range_gate = 0.05
)

fixture <- .sbs2_fixture(observations = 420L, seed = settings$fixture_seed)
initial <- list(c(0L, 0L, 0L), c(1L, 1L, 1L),
                c(1L, 0L, 1L), c(0L, 1L, 0L))

phase1_bridge <- run_phase1_bridge(fixture)

primary_started <- proc.time()[["elapsed"]]
primary_chains <- run_chains(
  fixture, settings$chain_seeds, settings$iterations, settings$burn,
  "primary", initial
)
primary_runtime <- proc.time()[["elapsed"]] - primary_started
primary <- pool_phase2(primary_chains, colnames(fixture$annotation))

direct_started <- proc.time()[["elapsed"]]
direct_chains <- run_chains(
  fixture, 20261081:20261084, settings$iterations, settings$burn,
  "direct", initial
)
direct_runtime <- proc.time()[["elapsed"]] - direct_started
direct <- pool_phase2(direct_chains, colnames(fixture$annotation))
comparison <- compare_pools(primary, direct)

all_included <- run_all_included_bridge(fixture)
all_excluded <- run_all_excluded_bridge(fixture)
zero_column <- run_zero_column(fixture)
permutation <- run_permutation(fixture, primary)
duplicate <- run_duplicate(fixture)
strong_signal <- run_signal_case(1.8, 20261320L, 20261321L)
weak_signal <- run_signal_case(0.35, 20261330L, 20261331L)
prior_only <- run_prior_only()

repeated_started <- proc.time()[["elapsed"]]
repeated <- run_repeated_simulation(20L)
repeated_runtime <- proc.time()[["elapsed"]] - repeated_started

intermediate <- which(primary$annotation_pip >= 0.1 &
                        primary$annotation_pip <= 0.9)
switch_gate <- if (length(intermediate)) {
  all(vapply(primary$switching, function(x) {
    all(x[intermediate, "zero_to_one"] > 0 &
          x[intermediate, "one_to_zero"] > 0)
  }, logical(1L)))
} else TRUE

qualification <- c(
  phase1_bridge = phase1_bridge$passed,
  primary_direct_pip = comparison$max_pip_difference <=
    settings$pip_comparison_gate,
  primary_direct_alpha = comparison$max_normalized_alpha_mean_difference <=
    settings$normalized_alpha_gate,
  primary_direct_q = comparison$max_q_difference <= settings$q_pi_gate,
  primary_direct_pi = comparison$max_pi_difference <= settings$q_pi_gate,
  alpha_rhat = max(primary$alpha_rhat, na.rm = TRUE) <= settings$rhat_gate,
  intercept_rhat = max(primary$intercept_rhat, na.rm = TRUE) <=
    settings$rhat_gate,
  pip_range = max(primary$pip_range) <= settings$pip_range_gate,
  intermediate_switching = switch_gate,
  all_included = all_included$max_alpha_mean_difference <= 0.06 &&
    all_included$max_q_difference <= 0.02 &&
    all_included$max_pi_difference <= 0.02,
  all_excluded = all_excluded$max_absolute_slope == 0 &&
    all_excluded$max_pi_row_sum_error <= 1e-12 &&
    all_excluded$finite_intercepts,
  zero_column = abs(zero_column$pooled_pip - fixture$pi_a) <= 0.025,
  permutation = permutation$max_pip_difference <= 0.035,
  duplicate = duplicate$symmetry_difference <= 0.06,
  strong_signal = all(strong_signal$pip[1:2] > 0.8),
  weak_signal = any(weak_signal$pip[1:2] > 0.1 &
                      weak_signal$pip[1:2] < 0.9),
  repeated_discrimination =
    repeated$pip_summary[["1"]]["mean"] >
      repeated$pip_summary[["0"]]["mean"],
  repeated_convergence = repeated$convergence_failures <= 2L
)

result <- list(
  settings = settings,
  fixture = list(
    observations = nrow(fixture$annotation),
    sticks = length(fixture$outcome),
    eligible_counts = lengths(fixture$eligible),
    continuation_counts = vapply(fixture$outcome, sum, integer(1L)),
    annotations = colnames(fixture$annotation),
    true_delta = fixture$true_delta,
    true_alpha = fixture$true_alpha,
    true_intercept = fixture$true_intercept
  ),
  phase1_bridge = phase1_bridge,
  primary = primary,
  direct = direct,
  comparison = comparison,
  primary_runtime_seconds = primary_runtime,
  direct_runtime_seconds = direct_runtime,
  all_included = all_included,
  all_excluded = all_excluded,
  zero_column = zero_column,
  permutation = permutation,
  duplicate = duplicate,
  strong_signal = strong_signal,
  weak_signal = weak_signal,
  prior_only = prior_only,
  repeated = repeated,
  repeated_runtime_seconds = repeated_runtime,
  qualification = qualification,
  passed = all(qualification)
)

cat("SBayesRC-S Phase-2 observed-d reference\n")
cat("========================================\n")
cat("Eligible counts:", paste(result$fixture$eligible_counts, collapse = "/"), "\n")
cat("Continuation counts:",
    paste(result$fixture$continuation_counts, collapse = "/"), "\n\n")
cat("Primary chain PIPs:\n")
print(round(primary$chain_pip, 4))
cat("Pooled primary/direct PIPs:\n")
print(round(rbind(primary = primary$annotation_pip,
                  direct = direct$annotation_pip), 4))
cat("PIP ranges:\n")
print(round(primary$pip_range, 4))
cat("Primary/direct comparison:\n")
print(unlist(comparison))
cat("Primary alpha R-hat:\n")
print(primary$alpha_rhat)
cat("Primary intercept R-hat:\n")
print(primary$intercept_rhat)
cat("All-included bridge:\n")
print(unlist(all_included))
cat("Zero-column result:\n")
print(zero_column)
cat("Strong/weak PIPs:\n")
print(rbind(strong = strong_signal$pip, weak = weak_signal$pip))
cat("Repeated-simulation PIP summaries:\n")
print(repeated$pip_summary)
cat("Qualification:\n")
print(qualification)
cat("Overall:", if (result$passed) "PASS" else "FAIL", "\n")

output_directory <- file.path(
  "results", "local", "sbayesrc_s_reference", "phase2"
)
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
compact <- result
compact$primary$delta <- compact$primary$alpha <- compact$primary$intercept <- NULL
compact$direct$delta <- compact$direct$alpha <- compact$direct$intercept <- NULL
compact$repeated$records <- NULL
saveRDS(compact, file.path(output_directory, "phase2_primary_summary.rds"))
write.csv(
  data.frame(
    annotation = colnames(fixture$annotation),
    pooled_pip = primary$annotation_pip,
    direct_pip = direct$annotation_pip,
    pip_range = primary$pip_range
  ),
  file.path(output_directory, "phase2_annotation_pips.csv"), row.names = FALSE
)
write.csv(
  repeated$annotation_level,
  file.path(output_directory, "phase2_repeated_recovery.csv"), row.names = FALSE
)
write.csv(
  data.frame(
    chain = rep(seq_along(primary$switching), each = 3L),
    annotation = rep(colnames(fixture$annotation), times = 4L),
    do.call(rbind, primary$switching)
  ),
  file.path(output_directory, "phase2_chain_diagnostics.csv"), row.names = FALSE
)

if (!result$passed) stop("SBayesRC-S Phase-2 qualification failed")
