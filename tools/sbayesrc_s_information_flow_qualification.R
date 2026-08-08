# Development/reference implementation.
# SBayesRC/SBayesRC-S Phase 4D hard-versus-soft information audit.
# Diagnostic only: no transition probability consumes these summaries.

devtools::load_all(quiet = TRUE)
source(file.path("tests", "testthat", "helper-sbayesrc-s-reference.R"))
source(file.path("tests", "testthat", "helper-sbayesrc-s-genomic-reference.R"))

output_dir <- file.path(
  "results", "local", "sbayesrc_s_reference", "phase4D"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

master_seed <- 20271200L
chain_seeds <- 20271101:20271104
initials <- list(
  excluded = c(0L, 0L, 0L), included = c(1L, 1L, 1L),
  mixed_1 = c(1L, 0L, 1L), mixed_2 = c(0L, 1L, 0L)
)
fixed_truth <- c(1L, 1L, 0L)
retained <- 1800L
burnin <- 400L

.rhat <- function(chains) {
  chains <- lapply(chains, as.matrix)
  n <- min(vapply(chains, nrow, integer(1L)))
  chains <- lapply(chains, function(x) x[seq_len(n), , drop = FALSE])
  vapply(seq_len(ncol(chains[[1L]])), function(column) {
    values <- lapply(chains, function(x) x[, column])
    if (all(vapply(values, function(x) all(x == values[[1L]][1L]),
                   logical(1L)))) return(1)
    means <- vapply(values, mean, numeric(1L))
    within <- mean(vapply(values, stats::var, numeric(1L)))
    if (!is.finite(within) || within <= 0) return(Inf)
    between <- n * stats::var(means)
    sqrt((((n - 1) / n) * within + between / n) / within)
  }, numeric(1L))
}

.ess <- function(x) {
  x <- as.numeric(x)
  if (length(x) < 4L || !is.finite(stats::var(x)) || stats::var(x) == 0)
    return(length(x))
  correlation <- as.numeric(stats::acf(
    x, plot = FALSE, lag.max = min(1000L, floor(length(x) / 2L)),
    demean = TRUE
  )$acf)[-1L]
  pair_count <- floor(length(correlation) / 2L)
  if (pair_count == 0L) return(length(x))
  pair <- correlation[2L * seq_len(pair_count) - 1L] +
    correlation[2L * seq_len(pair_count)]
  positive <- pair[cumprod(pair > 0) == 1]
  min(length(x), length(x) / max(1, 1 + 2 * sum(positive)))
}

.safe_cor <- function(x, y) {
  if (stats::sd(x) == 0 || stats::sd(y) == 0) return(NA_real_)
  stats::cor(x, y)
}
.chain <- function(fit) fit$chains[[1L]][[1L]]
.flow <- function(fit) .chain(fit)$information_flow
.trace_active <- function(flow, kind, marker_count) {
  marker_count - flow[[paste0(kind, "_component_count_trace")]][, 1L]
}

.summarize_fits <- function(name, fits, elapsed, model) {
  flow <- lapply(fits, .flow)
  chains <- lapply(fits, .chain)
  marker_count <- nrow(flow[[1L]]$rb_comp_prob)
  component_count <- ncol(flow[[1L]]$rb_comp_prob)
  stick_count <- component_count - 1L
  prior_active <- lapply(flow, .trace_active, "prior", marker_count)
  soft_active <- lapply(flow, .trace_active, "rb", marker_count)
  hard_active <- lapply(flow, .trace_active, "hard", marker_count)
  correlation <- do.call(rbind, Map(function(prior, soft, hard) c(
    prior_soft = .safe_cor(prior, soft), soft_hard = .safe_cor(soft, hard),
    prior_hard = .safe_cor(prior, hard)
  ), prior_active, soft_active, hard_active))
  hard_probability <- lapply(chains, function(x) x$component$prob)
  rb_probability <- lapply(flow, `[[`, "rb_comp_prob")
  hard_pip <- lapply(hard_probability, function(x) 1 - x[, 1L])
  rb_pip <- lapply(rb_probability, function(x) 1 - x[, 1L])
  hard_pip_range <- apply(do.call(cbind, hard_pip), 1L, function(x) diff(range(x)))
  rb_pip_range <- apply(do.call(cbind, rb_pip), 1L, function(x) diff(range(x)))
  hard_component_range <- apply(
    array(unlist(hard_probability), c(marker_count, component_count, length(fits))),
    c(1L, 2L), function(x) diff(range(x))
  )
  rb_component_range <- apply(
    array(unlist(rb_probability), c(marker_count, component_count, length(fits))),
    c(1L, 2L), function(x) diff(range(x))
  )
  alpha <- lapply(chains, function(x) as.matrix(x$convergence_trace$alpha))

  stick <- lapply(seq_len(stick_count), function(stick_index) {
    q_column <- 3L * stick_index
    hard_q <- lapply(flow, function(x) x$hard_stick_trace[, q_column])
    soft_q <- lapply(flow, function(x) x$soft_stick_trace[, q_column])
    data.frame(
      scenario = name, stick = stick_index,
      hard_eligible_mean = mean(vapply(flow, function(x)
        mean(x$hard_stick_trace[, q_column - 2L]), numeric(1L))),
      soft_eligible_mean = mean(vapply(flow, function(x)
        mean(x$soft_stick_trace[, q_column - 2L]), numeric(1L))),
      hard_success_mean = mean(vapply(flow, function(x)
        mean(x$hard_stick_trace[, q_column - 1L]), numeric(1L))),
      soft_success_mean = mean(vapply(flow, function(x)
        mean(x$soft_stick_trace[, q_column - 1L]), numeric(1L))),
      hard_q_mean = mean(vapply(hard_q, mean, numeric(1L))),
      soft_q_mean = mean(vapply(soft_q, mean, numeric(1L))),
      hard_q_range = diff(range(vapply(hard_q, mean, numeric(1L)))),
      soft_q_range = diff(range(vapply(soft_q, mean, numeric(1L)))),
      hard_q_rhat = max(.rhat(hard_q)), soft_q_rhat = max(.rhat(soft_q)),
      min_hard_q_ess = min(vapply(hard_q, .ess, numeric(1L))),
      min_soft_q_ess = min(vapply(soft_q, .ess, numeric(1L)))
    )
  })
  stick <- do.call(rbind, stick)
  hard_annotation_mean <- do.call(rbind, lapply(flow, function(x)
    colMeans(x$hard_annotation_information)))
  soft_annotation_mean <- do.call(rbind, lapply(flow, function(x)
    colMeans(x$soft_annotation_information)))
  annotation <- data.frame(
    scenario = name, index = seq_len(ncol(hard_annotation_mean)),
    hard_chain_range = apply(hard_annotation_mean, 2L, function(x) diff(range(x))),
    soft_chain_range = apply(soft_annotation_mean, 2L, function(x) diff(range(x)))
  )
  gain <- do.call(rbind, lapply(flow, function(x) colMeans(x$information_gain)))
  colnames(gain) <- flow[[1L]]$information_gain_layout
  active <- data.frame(
    scenario = name, model = model, marker_count = marker_count,
    prior_chain_means = paste(round(vapply(prior_active, mean, numeric(1L)), 4),
                              collapse = "/"),
    soft_chain_means = paste(round(vapply(soft_active, mean, numeric(1L)), 4),
                             collapse = "/"),
    hard_chain_means = paste(round(vapply(hard_active, mean, numeric(1L)), 4),
                             collapse = "/"),
    prior_active_range = diff(range(vapply(prior_active, mean, numeric(1L)))),
    soft_active_range = diff(range(vapply(soft_active, mean, numeric(1L)))),
    hard_active_range = diff(range(vapply(hard_active, mean, numeric(1L)))),
    prior_active_rhat = max(.rhat(prior_active)),
    soft_active_rhat = max(.rhat(soft_active)),
    hard_active_rhat = max(.rhat(hard_active)),
    min_prior_active_ess = min(vapply(prior_active, .ess, numeric(1L))),
    min_soft_active_ess = min(vapply(soft_active, .ess, numeric(1L))),
    min_hard_active_ess = min(vapply(hard_active, .ess, numeric(1L))),
    max_alpha_rhat = max(.rhat(alpha)),
    min_alpha_ess = min(vapply(alpha, function(x)
      min(apply(x, 2L, .ess)), numeric(1L))),
    mean_prior_soft_correlation = mean(correlation[, "prior_soft"], na.rm = TRUE),
    mean_soft_hard_correlation = mean(correlation[, "soft_hard"], na.rm = TRUE),
    mean_prior_hard_correlation = mean(correlation[, "prior_hard"], na.rm = TRUE),
    hard_snp_range_median = median(hard_pip_range),
    hard_snp_range_p90 = unname(stats::quantile(hard_pip_range, 0.9)),
    hard_snp_range_max = max(hard_pip_range),
    rb_snp_range_median = median(rb_pip_range),
    rb_snp_range_p90 = unname(stats::quantile(rb_pip_range, 0.9)),
    rb_snp_range_max = max(rb_pip_range),
    hard_component_range_median = median(hard_component_range),
    hard_component_range_p90 = unname(stats::quantile(hard_component_range, 0.9)),
    hard_component_range_max = max(hard_component_range),
    rb_component_range_median = median(rb_component_range),
    rb_component_range_p90 = unname(stats::quantile(rb_component_range, 0.9)),
    rb_component_range_max = max(rb_component_range),
    hard_rb_pip_mad = mean(abs(Reduce(`+`, hard_pip) / length(hard_pip) -
                                 Reduce(`+`, rb_pip) / length(rb_pip))),
    hard_rb_pip_max = max(abs(Reduce(`+`, hard_pip) / length(hard_pip) -
                                     Reduce(`+`, rb_pip) / length(rb_pip))),
    entropy_mean = mean(gain[, "entropy_mean"]),
    entropy_median = mean(gain[, "entropy_median"]),
    kl_mean = mean(gain[, "kl_mean"]),
    kl_median = mean(gain[, "kl_median"]),
    kl_p90 = mean(gain[, "kl_p90"]),
    tv_mean = mean(gain[, "tv_mean"]),
    runtime_seconds = elapsed
  )
  list(name = name, active = active, stick = stick, annotation = annotation,
       fits = fits, correlations = correlation,
       hard_probability = hard_probability, rb_probability = rb_probability)
}

.run_backend <- function(fixture, name, fixed_delta = integer(),
                         update_pi_A = TRUE, update_tau2 = TRUE,
                         selection_enabled = TRUE, initial = initials) {
  started <- proc.time()[["elapsed"]]
  fits <- lapply(seq_along(initial), function(chain) .sbs4b_run(
    fixture, chain_seeds[chain], retained, burnin,
    fixed_delta = fixed_delta, updateB = TRUE, updateE = FALSE,
    initial_delta = initial[[chain]], update_hierarchy = TRUE,
    update_pi_A = update_pi_A, update_tau2 = update_tau2,
    selection_enabled = selection_enabled, information_diagnostics = TRUE
  ))
  .summarize_fits(name, fits, proc.time()[["elapsed"]] - started,
                  if (selection_enabled) "SBayesRC-S" else "SBayesRC")
}

.frozen_information_chain <- function(fixture, seed, initial_delta) {
  component <- as.integer(fixture$comp_init[[1L]])
  eligible <- list(seq_along(component), which(component > 0L),
                   which(component > 1L))
  outcome <- list(as.integer(component > 0L),
                  as.integer(component[component > 0L] > 1L),
                  as.integer(component[component > 1L] > 2L))
  hierarchy <- .st_bayesrc_selection_hierarchy(
    fixture$annotation, eligible, outcome, fixture$alpha_init,
    as.integer(initial_delta), 0.35, fixture$tau2_init,
    1, 1, 3, 1.6, fixture$intercept_prior, 1e-12,
    retained + burnin, burnin, seed, integer()
  )
  marker_count <- nrow(fixture$A)
  component_count <- length(fixture$gamma)
  hard_prob <- matrix(0, marker_count, component_count)
  hard_prob[cbind(seq_len(marker_count), component + 1L)] <- 1
  prior_sum <- rb_sum <- matrix(0, marker_count, component_count)
  prior_trace <- rb_trace <- hard_trace <- matrix(0, retained, component_count)
  hard_trace[,] <- matrix(rep(colSums(hard_prob), each = retained), retained)
  hard_stick <- soft_stick <- matrix(0, retained, 3L * (component_count - 1L))
  hard_annotation <- soft_annotation <- matrix(
    0, retained, 3L * (component_count - 1L) * ncol(fixture$A)
  )
  information_gain <- matrix(0, retained, 9L)
  score <- fixture$wy[[1L]]
  wi <- fixture$ww[[1L]]
  vb <- fixture$B[1L, 1L]
  ve <- fixture$E[1L, 1L]
  log_bf <- matrix(0, marker_count, component_count)
  for (k in 2:component_count) {
    vbk <- vb * fixture$gamma[k]
    denominator <- ve + wi * vbk
    log_bf[, k] <- 0.5 * log(ve / denominator) +
      0.5 * score^2 * vbk / (ve * denominator)
  }
  for (draw in seq_len(retained)) {
    alpha <- matrix(hierarchy$alpha_draws[draw, , ], nrow = ncol(fixture$A))
    q <- stats::pnorm(fixture$A %*% alpha)
    prior <- t(apply(q, 1L, sblr:::.sbayesrc_stick_to_component_prob))
    logp <- log(pmax(prior, 1e-300)) + log_bf
    rb <- exp(logp - apply(logp, 1L, max))
    rb <- rb / rowSums(rb)
    summary <- .st_bayesrc_information_summary(
      log(rb), prior, component, fixture$A
    )
    prior_sum <- prior_sum + prior
    rb_sum <- rb_sum + rb
    prior_trace[draw, ] <- colSums(prior)
    rb_trace[draw, ] <- colSums(rb)
    hard_stick[draw, ] <- summary$hard_stick_trace
    soft_stick[draw, ] <- summary$soft_stick_trace
    hard_annotation[draw, ] <- summary$hard_annotation_information
    soft_annotation[draw, ] <- summary$soft_annotation_information
    information_gain[draw, ] <- summary$information_gain
  }
  flow <- list(
    prior_comp_prob = prior_sum / retained, rb_comp_prob = rb_sum / retained,
    rb_dm = 1 - rb_sum[, 1L] / retained,
    prior_component_count_trace = prior_trace,
    rb_component_count_trace = rb_trace,
    hard_component_count_trace = hard_trace,
    hard_stick_trace = hard_stick, soft_stick_trace = soft_stick,
    hard_annotation_information = hard_annotation,
    soft_annotation_information = soft_annotation,
    information_gain = information_gain,
    information_gain_layout = c(
      "entropy_mean", "entropy_median", "entropy_p90", "kl_mean",
      "kl_median", "kl_p90", "tv_mean", "tv_median", "tv_p90"
    )
  )
  list(marker = list(dm = rowSums(hard_prob[, -1L, drop = FALSE])),
       component = list(prob = hard_prob),
       annotation = list(),
       convergence_trace = list(alpha = matrix(
         hierarchy$alpha_draws, nrow = retained
       )), information_flow = flow)
}

.run_frozen <- function(fixture, name = "B0_frozen_allocations") {
  started <- proc.time()[["elapsed"]]
  raw <- lapply(seq_along(initials), function(chain)
    .frozen_information_chain(fixture, chain_seeds[chain], initials[[chain]]))
  fits <- lapply(raw, function(chain) list(chains = list(list(chain))))
  .summarize_fits(name, fits, proc.time()[["elapsed"]] - started,
                  "SBayesRC-S frozen allocations")
}

.material_stability <- function(hard_range, soft_range, hard_rhat, soft_rhat,
                                hard_ess, soft_ess) {
  criteria <- c(
    is.finite(hard_range) && hard_range > 0 && soft_range <= 0.5 * hard_range,
    is.finite(hard_rhat) && hard_rhat > 1.05 && soft_rhat <= 1.05,
    is.finite(hard_ess) && hard_ess > 0 && soft_ess >= 2 * hard_ess
  )
  sum(criteria) >= 2L
}

fixture_160 <- .sbs4b_fixture(160L, 20270930L)
B0 <- .run_frozen(fixture_160)
B1 <- .run_backend(
  fixture_160, "B1_fixed_delta_all", rep(1L, 3L), FALSE, FALSE,
  initial = rep(list(rep(1L, 3L)), 4L)
)
B2 <- .run_backend(
  fixture_160, "B2_fixed_delta_truth", fixed_truth, FALSE, FALSE,
  initial = rep(list(fixed_truth), 4L)
)
B3 <- .run_backend(fixture_160, "B3_full_hierarchy")
standard <- .run_backend(
  fixture_160, "standard_SBayesRC", integer(), FALSE, FALSE,
  selection_enabled = FALSE, initial = rep(list(fixed_truth), 4L)
)
audit <- list(B0 = B0, B1 = B1, B2 = B2, B3 = B3, standard = standard)

# Marker-count ladder: new deterministic simulations, never marker replication.
marker_count <- c(160L, 500L, 2000L)
marker_ladder <- lapply(marker_count, function(m) {
  fixture <- .sbs4b_fixture(m, master_seed + m, annotation_contrast = 1)
  .run_backend(
    fixture, paste0("marker_M", m), fixed_truth, FALSE, FALSE,
    initial = rep(list(fixed_truth), 4L)
  )
})

# Annotation-strength ladder is preregistered as 0x, 0.5x, 1x, and 2x. The
# same seed fixes annotations and noise; only the generating contrast changes.
contrast <- c(0, 0.5, 1, 2)
contrast_name <- c("null", "weak", "baseline", "strong")
strength_ladder <- Map(function(scale, label) {
  fixture <- .sbs4b_fixture(
    160L, master_seed + 160L, annotation_contrast = scale
  )
  .run_backend(
    fixture, paste0("strength_", label), fixed_truth, FALSE, FALSE,
    initial = rep(list(fixed_truth), 4L)
  )
}, contrast, contrast_name)

# Phase 4D-D selection rule frozen before its C5 result: use the largest marker
# count (M=2000), because it maximizes aggregate annotation information without
# changing the baseline contrast or inference prior. Baseline C5 is B3.
fixture_2000_c5 <- .sbs4b_fixture(
  2000L, master_seed + 2000L, annotation_contrast = 1
)
C5_2000 <- .run_backend(fixture_2000_c5, "C5_marker_M2000")

all_results <- c(audit, marker_ladder, strength_ladder,
                 list(C5_2000 = C5_2000))
active_table <- do.call(rbind, lapply(all_results, `[[`, "active"))
stick_table <- do.call(rbind, lapply(all_results, `[[`, "stick"))
annotation_table <- do.call(rbind, lapply(all_results, `[[`, "annotation"))

key <- list(B1 = B1, B2 = B2, B3 = B3)
stability <- unlist(lapply(key, function(result) {
  active <- result$active
  values <- .material_stability(
    active$hard_active_range, active$soft_active_range,
    active$hard_active_rhat, active$soft_active_rhat,
    active$min_hard_active_ess, active$min_soft_active_ess
  )
  stick <- apply(result$stick, 1L, function(x) .material_stability(
    as.numeric(x[["hard_q_range"]]), as.numeric(x[["soft_q_range"]]),
    as.numeric(x[["hard_q_rhat"]]), as.numeric(x[["soft_q_rhat"]]),
    as.numeric(x[["min_hard_q_ess"]]), as.numeric(x[["min_soft_q_ess"]])
  ))
  c(active = values, setNames(stick, paste0("stick", seq_along(stick))))
}))
information_decision <- if (sum(stability) >= ceiling(length(stability) / 2)) {
  "SBS4D-R1"
} else if (sum(stability) <= floor(length(stability) / 4)) {
  "SBS4D-R2"
} else {
  "SBS4D-R3"
}

largest <- marker_ladder[[3L]]$active
baseline_scale <- strength_ladder[[3L]]$active
strong_scale <- strength_ladder[[4L]]$active
scale_qualifier <- character()
if (largest$hard_active_rhat <= 1.05 && largest$soft_active_rhat <= 1.05 &&
    largest$max_alpha_rhat <= 1.05) scale_qualifier <- c(scale_qualifier, "SCALE-D1")
if (largest$hard_active_rhat > 1.05 && largest$soft_active_rhat <= 1.05 &&
    largest$soft_active_range <= 0.5 * largest$hard_active_range)
  scale_qualifier <- c(scale_qualifier, "SCALE-D2")
if (largest$hard_active_rhat > 1.05 && largest$soft_active_rhat > 1.05)
  scale_qualifier <- c(scale_qualifier, "SCALE-D3")
if (strong_scale$hard_active_range > 1.25 * baseline_scale$hard_active_range &&
    strong_scale$max_alpha_rhat > baseline_scale$max_alpha_rhat)
  scale_qualifier <- c(scale_qualifier, "SCALE-D4")
if (length(scale_qualifier) == 0L) scale_qualifier <- "SCALE-D5"

write.csv(active_table, file.path(output_dir, "hard_soft_160_summary.csv"),
          row.names = FALSE)
write.csv(stick_table, file.path(output_dir, "stick_information_summary.csv"),
          row.names = FALSE)
write.csv(annotation_table,
          file.path(output_dir, "annotation_information_summary.csv"),
          row.names = FALSE)
write.csv(do.call(rbind, lapply(marker_ladder, `[[`, "active")),
          file.path(output_dir, "marker_scale_ladder.csv"), row.names = FALSE)
write.csv(do.call(rbind, lapply(strength_ladder, `[[`, "active")),
          file.path(output_dir, "annotation_strength_ladder.csv"),
          row.names = FALSE)
write.csv(standard$active,
          file.path(output_dir, "standard_sbayesrc_bridge.csv"), row.names = FALSE)
write.csv(rbind(B3$active, C5_2000$active),
          file.path(output_dir, "full_hierarchy_confirmation.csv"),
          row.names = FALSE)
write.csv(active_table[, c("scenario", "runtime_seconds")],
          file.path(output_dir, "runtime_summary.csv"), row.names = FALSE)

result <- list(
  settings = list(master_seed = master_seed, chain_seeds = chain_seeds,
                  retained = retained, burnin = burnin,
                  fixed_delta = fixed_truth,
                  marker_count = marker_count, contrast = contrast,
                  full_hierarchy_confirmation_rule = "largest M=2000"),
  audit = audit, marker_ladder = marker_ladder,
  strength_ladder = strength_ladder, C5_2000 = C5_2000,
  stability = stability, information_decision = information_decision,
  scale_qualifier = scale_qualifier,
  active_table = active_table, stick_table = stick_table,
  annotation_table = annotation_table
)
saveRDS(result, file.path(output_dir, "phase4d_information_flow.rds"))
print(active_table)
print(stick_table)
cat("material-stability keys:", sum(stability), "/", length(stability), "\n")
cat("information decision:", information_decision, "\n")
cat("scale qualifier:", paste(scale_qualifier, collapse = ", "), "\n")
