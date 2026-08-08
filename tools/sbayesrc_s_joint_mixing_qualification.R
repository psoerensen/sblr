# Development/reference implementation.
# SBayesRC-S Phase 4C coupling and invariant-schedule qualification.
# Not a supported public sampler or API.

devtools::load_all(quiet = TRUE)
source(file.path("tests", "testthat", "helper-sbayesrc-s-genomic-reference.R"))

output_dir <- file.path(
  "results", "local", "sbayesrc_s_reference", "phase4C"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

initials <- list(
  excluded = c(0L, 0L, 0L), included = c(1L, 1L, 1L),
  mixed_1 = c(1L, 0L, 1L), mixed_2 = c(0L, 1L, 0L)
)

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
  if (length(correlation) < 2L) return(length(x))
  pair_count <- floor(length(correlation) / 2L)
  pair <- correlation[2L * seq_len(pair_count) - 1L] +
    correlation[2L * seq_len(pair_count)]
  positive <- pair[cumprod(pair > 0) == 1]
  min(length(x), length(x) / max(1, 1 + 2 * sum(positive)))
}

.safe_max <- function(x) if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)
.safe_min <- function(x) if (all(is.na(x))) NA_real_ else min(x, na.rm = TRUE)
.safe_cor <- function(x, y) {
  if (stats::sd(x) == 0 || stats::sd(y) == 0) return(NA_real_)
  stats::cor(x, y)
}
.chain <- function(fit) fit$chains[[1L]][[1L]]
.trace <- function(fit) .chain(fit)$convergence_trace

.summarize_genomic <- function(name, fits, elapsed) {
  traces <- lapply(fits, .trace)
  annotation <- lapply(fits, function(x) .chain(x)$annotation)
  delta <- lapply(traces, function(x) as.matrix(x$annotation_delta))
  alpha <- lapply(traces, function(x) as.matrix(x$alpha))
  pi_a <- lapply(traces, function(x) as.matrix(x$annotation_pi_A))
  tau2 <- lapply(traces, function(x) as.matrix(x$sigmaSqAlpha))
  included <- lapply(traces, function(x) as.matrix(x$annotation_included_count))
  active <- lapply(traces, function(x) as.matrix(x$realized_active_count))
  component <- lapply(traces, function(x) as.matrix(x$component_count))
  chain_pip <- do.call(rbind, lapply(delta, colMeans))
  pooled_pip <- colMeans(chain_pip)
  pip_range <- apply(chain_pip, 2L, function(x) diff(range(x)))
  correlations <- do.call(rbind, Map(function(d, p, tau, m, a, comp) {
    data.frame(
      M_active = .safe_cor(m[, 1L], a[, 1L]),
      piA_active = .safe_cor(p[, 1L], a[, 1L]),
      tau1_active = .safe_cor(tau[, 1L], a[, 1L]),
      delta1_active = .safe_cor(d[, 1L], a[, 1L]),
      delta1_component1 = .safe_cor(d[, 1L], comp[, 2L])
    )
  }, delta, pi_a, tau2, included, active, component))
  list(
    name = name, elapsed_seconds = elapsed,
    chain_pip = chain_pip, pooled_pip = pooled_pip, pip_range = pip_range,
    switches = do.call(rbind, lapply(annotation, function(x)
      as.numeric(x$annotation_switches))),
    alpha_rhat = .rhat(alpha), alpha_ess = vapply(alpha, function(x)
      min(apply(x, 2L, .ess)), numeric(1L)),
    pi_a_rhat = .rhat(pi_a), pi_a_ess = vapply(pi_a, function(x) .ess(x[, 1L]),
                                               numeric(1L)),
    tau2_rhat = .rhat(tau2), tau2_ess = vapply(tau2, function(x)
      min(apply(x, 2L, .ess)), numeric(1L)),
    M_rhat = .rhat(included), M_ess = vapply(included, function(x) .ess(x[, 1L]),
                                             numeric(1L)),
    active_rhat = .rhat(active), active_ess = vapply(active, function(x)
      .ess(x[, 1L]), numeric(1L)),
    component_rhat = .rhat(component), component_ess = vapply(component,
      function(x) min(apply(x, 2L, .ess)), numeric(1L)),
    correlations = correlations,
    beta = Reduce(`+`, lapply(fits, function(x) x$marker$bm[, 1L])) /
      length(fits),
    snp_pip = Reduce(`+`, lapply(fits, function(x) x$marker$dm[, 1L])) /
      length(fits),
    empty = lapply(annotation, `[[`, "empty_stick_diagnostics")
  )
}

.row <- function(x, stage) data.frame(
  stage = stage, configuration = x$name,
  max_pip_range = max(x$pip_range),
  max_alpha_rhat = .safe_max(x$alpha_rhat),
  min_alpha_ess = .safe_min(x$alpha_ess),
  pi_a_rhat = .safe_max(x$pi_a_rhat), min_pi_a_ess = .safe_min(x$pi_a_ess),
  max_tau2_rhat = .safe_max(x$tau2_rhat), min_tau2_ess = .safe_min(x$tau2_ess),
  M_rhat = .safe_max(x$M_rhat), min_M_ess = .safe_min(x$M_ess),
  active_rhat = .safe_max(x$active_rhat), min_active_ess = .safe_min(x$active_ess),
  max_component_rhat = .safe_max(x$component_rhat),
  min_component_ess = .safe_min(x$component_ess),
  runtime_seconds = x$elapsed_seconds
)

.run_genomic <- function(fixture, name, nit, nburn,
                         fixed_delta = integer(), update_pi_A = TRUE,
                         update_tau2 = TRUE, hierarchy_sweeps = 1L,
                         genomic_sweeps = 1L, initial = initials) {
  started <- proc.time()[["elapsed"]]
  fits <- lapply(seq_along(initial), function(i) .sbs4b_run(
    fixture, 20271000L + 100L * match(name, unique(name)) + i,
    nit, nburn, fixed_delta = fixed_delta, updateB = TRUE, updateE = FALSE,
    initial_delta = initial[[i]], update_hierarchy = TRUE,
    update_pi_A = update_pi_A, update_tau2 = update_tau2,
    hierarchy_sweeps = hierarchy_sweeps, genomic_sweeps = genomic_sweeps
  ))
  .summarize_genomic(name, fits, proc.time()[["elapsed"]] - started)
}

.build_sticks <- function(component) {
  list(
    eligible = list(seq_along(component), which(component > 0L),
                    which(component > 1L)),
    outcome = list(as.integer(component > 0L),
                   as.integer(component[component > 0L] > 1L),
                   as.integer(component[component > 1L] > 2L))
  )
}

.run_c0 <- function(fixture, nit = 1800L, nburn = 400L) {
  sticks <- .build_sticks(as.integer(fixture$comp_init[[1L]]))
  chains <- lapply(seq_along(initials), function(i) .st_bayesrc_selection_hierarchy(
    fixture$annotation, sticks$eligible, sticks$outcome,
    fixture$alpha_init, as.integer(initials[[i]]), 0.35, fixture$tau2_init,
    1, 1, 3, 1.6, fixture$intercept_prior, 1e-12,
    nit + nburn, nburn, 20271100L + i, integer()
  ))
  alpha <- lapply(chains, function(x) {
    value <- x$alpha_draws
    matrix(value, nrow = dim(value)[1L])
  })
  delta <- lapply(chains, function(x) as.matrix(x$delta_draws))
  pi_a <- lapply(chains, function(x) as.matrix(x$pi_A_draws))
  tau2 <- lapply(chains, function(x) as.matrix(x$tau2_draws))
  included <- lapply(chains, function(x) as.matrix(x$included_draws))
  chain_pip <- do.call(rbind, lapply(delta, colMeans))
  list(
    name = "C0_frozen_allocations", elapsed_seconds = NA_real_,
    chain_pip = chain_pip, pooled_pip = colMeans(chain_pip),
    pip_range = apply(chain_pip, 2L, function(x) diff(range(x))),
    switches = do.call(rbind, lapply(chains, function(x) colSums(x$switches))),
    alpha_rhat = .rhat(alpha), alpha_ess = vapply(alpha, function(x)
      min(apply(x, 2L, .ess)), numeric(1L)),
    pi_a_rhat = .rhat(pi_a), pi_a_ess = vapply(pi_a, function(x) .ess(x[, 1L]),
                                               numeric(1L)),
    tau2_rhat = .rhat(tau2), tau2_ess = vapply(tau2, function(x)
      min(apply(x, 2L, .ess)), numeric(1L)),
    M_rhat = .rhat(included), M_ess = vapply(included, function(x) .ess(x[, 1L]),
                                             numeric(1L)),
    active_rhat = NA_real_, active_ess = NA_real_, component_rhat = NA_real_,
    component_ess = NA_real_, correlations = NULL, empty = NULL
  )
}

fixture <- .sbs4b_fixture(160L, 20270930L)

# 4C-A: same scientific fixture, only conditional/dynamic state differs.
c0_started <- proc.time()[["elapsed"]]
C0 <- .run_c0(fixture)
C0$elapsed_seconds <- proc.time()[["elapsed"]] - c0_started
C1_all <- .run_genomic(
  fixture, "C1_fixed_delta_all", 1800L, 400L, rep(1L, 3L),
  initial = rep(list(rep(1L, 3L)), 4L)
)
C1_truth <- .run_genomic(
  fixture, "C1_fixed_delta_truth", 1800L, 400L, c(1L, 1L, 0L),
  initial = rep(list(c(1L, 1L, 0L)), 4L)
)
C2 <- .run_genomic(
  fixture, "C2_dynamic_delta_fixed_pi_tau", 1800L, 400L,
  update_pi_A = FALSE, update_tau2 = FALSE
)
C3 <- .run_genomic(
  fixture, "C3_dynamic_pi_fixed_tau", 1800L, 400L,
  update_pi_A = TRUE, update_tau2 = FALSE
)
C4 <- .run_genomic(
  fixture, "C4_dynamic_tau_fixed_pi", 1800L, 400L,
  update_pi_A = FALSE, update_tau2 = TRUE
)
C5 <- .run_genomic(fixture, "C5_full", 1800L, 400L)
ladder <- list(C0 = C0, C1_all = C1_all, C1_truth = C1_truth,
               C2 = C2, C3 = C3, C4 = C4, C5 = C5)
ladder_table <- do.call(rbind, lapply(ladder, .row, stage = "coupling"))
write.csv(ladder_table, file.path(output_dir, "coupling_ladder.csv"), row.names = FALSE)

# Chain-length screen. The 1x result is exactly C5; scientific inputs and seeds
# remain fixed while burn-in and retained length scale together.
length_2x <- .run_genomic(fixture, "C5_length_2x", 3600L, 800L)
length_4x <- .run_genomic(fixture, "C5_length_4x", 7200L, 1600L)
lengths <- list(C5, length_2x, length_4x)
length_table <- do.call(rbind, lapply(lengths, .row, stage = "chain_length"))
length_table$length_multiplier <- c(1, 2, 4)
write.csv(length_table, file.path(output_dir, "chain_length_screen.csv"),
          row.names = FALSE)

# 4C-B: invariant repeated Gibbs kernels, using the existing schedule contract.
S1_H2 <- .run_genomic(fixture, "S1_H2_G1", 1800L, 400L,
                      hierarchy_sweeps = 2L)
S1_H5 <- .run_genomic(fixture, "S1_H5_G1", 1800L, 400L,
                      hierarchy_sweeps = 5L)
S2_G2 <- .run_genomic(fixture, "S2_H1_G2", 1800L, 400L,
                      genomic_sweeps = 2L)
S2_G5 <- .run_genomic(fixture, "S2_H1_G5", 1800L, 400L,
                      genomic_sweeps = 5L)
schedules <- list(S0 = C5, S1_H2 = S1_H2, S1_H5 = S1_H5,
                  S2_G2 = S2_G2, S2_G5 = S2_G5)
schedule_table <- do.call(rbind, lapply(schedules, .row, stage = "schedule"))
schedule_table$min_alpha_ess_per_second <- schedule_table$min_alpha_ess /
  schedule_table$runtime_seconds
schedule_table$min_active_ess_per_second <- schedule_table$min_active_ess /
  schedule_table$runtime_seconds
schedule_table$min_M_ess_per_second <- schedule_table$min_M_ess /
  schedule_table$runtime_seconds
write.csv(schedule_table, file.path(output_dir, "schedule_screen.csv"),
          row.names = FALSE)

qualifies <- function(x) isTRUE(
  max(x$pip_range) <= 0.10 && .safe_max(x$alpha_rhat) <= 1.05 &&
    .safe_max(x$M_rhat) <= 1.05 && .safe_max(x$active_rhat) <= 1.05 &&
    .safe_max(x$component_rhat) <= 1.05 && .safe_min(x$alpha_ess) >= 100 &&
    .safe_min(x$M_ess) >= 100 && .safe_min(x$active_ess) >= 100
)
qualified_schedules <- names(Filter(qualifies, schedules))

result <- list(
  settings = list(marker_count = 160L, retained = 1800L, burnin = 400L,
                  initials = initials, pip_range_gate = 0.10,
                  rhat_gate = 1.05, minimum_ess_gate = 100),
  ladder = ladder, ladder_table = ladder_table,
  chain_length = lengths, chain_length_table = length_table,
  schedules = schedules, schedule_table = schedule_table,
  qualified_schedules = qualified_schedules
)
saveRDS(result, file.path(output_dir, "phase4c_diagnostic_screen.rds"))
print(ladder_table)
print(length_table)
print(schedule_table)
cat("qualified schedules:", paste(qualified_schedules, collapse = ", "), "\n")
