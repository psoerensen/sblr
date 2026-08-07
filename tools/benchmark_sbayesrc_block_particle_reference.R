#!/usr/bin/env Rscript

# Development-only fixed-alpha particle-Gibbs path-diversity screen. It reads
# the preserved large B0 state and never writes to sblrbench.

sblr_root <- normalizePath(".", winslash = "/", mustWork = TRUE)
bench_root <- normalizePath("../sblrbench", winslash = "/", mustWork = TRUE)
output_root <- file.path(sblr_root, "results", "local", "sbayesrc_block_particle")
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
pkgload::load_all(sblr_root, quiet = TRUE)

bundle_path <- file.path(bench_root, "results/local/06_annotation_models",
                         "large_feasibility/prepared_bundle.rds")
captured_path <- file.path(bench_root, "results/local/06_annotation_models",
                           "large_feasibility/continuation",
                           "b0_iteration0_no_updateE.rds")
stopifnot(file.exists(bundle_path), file.exists(captured_path))
bundle <- readRDS(bundle_path)
captured <- readRDS(captured_path)

cache_path <- file.path(output_root, "frozen_block_1_contract.rds")
if (file.exists(cache_path)) {
  block <- readRDS(cache_path)
} else {
  st <- sblr:::.mtblr_normalize_stats(bundle$gwas$stats)
  provenance <- st$genotype_provenance[[1L]]
  provenance$cls <- unname(provenance$cls)
  provenance$af <- unname(provenance$af)
  old_directory <- setwd(bench_root)
  on.exit(setwd(old_directory), add = TRUE)
  reference <- sblr:::.mtblr_block_eigen_reference(bundle$glist, provenance)
  contract <- sblr:::stblr_block_low_rank_contract_internal(
    reference$bed_files, reference$n_bed, reference$cls, reference$rows,
    reference$af, as.integer(bundle$blocks$block_start - 1L),
    matrix(st$wy[[1L]], nrow = 1L), as.numeric(captured$b[, 1L]),
    bundle$spec$block$eigen_prop, st$yy[[1L]], 0)
  setwd(old_directory)
  block <- list(Q = contract$factor[[1L]],
                w = as.numeric(contract$transformed_score[[1L]]),
                n = st$n, yy = st$yy[[1L]])
  saveRDS(block, cache_path, compress = FALSE)
}

gamma <- as.numeric(captured$mixture_var)
probability_all <- matrix(as.numeric(captured$pi_final[1L, ]),
                          nrow = 500L, ncol = length(gamma), byrow = TRUE)
component_all <- as.integer(captured$component[seq_len(500L), 1L])
beta_all <- as.numeric(captured$b[seq_len(500L), 1L])
vb <- as.numeric(captured$vbs[1L, 1L])
ve <- block$yy / (block$n - 1)
particle_counts <- c(8L, 16L, 32L, 64L)

run_size <- function(marker_count, repetitions) {
  selected <- seq_len(marker_count)
  excluded <- setdiff(seq_len(500L), selected)
  Q <- block$Q[, selected, drop = FALSE]
  w <- block$w
  if (length(excluded))
    w <- w - as.numeric(block$Q[, excluded, drop = FALSE] %*% beta_all[excluded])
  rows <- list()
  for (particles in particle_counts) {
    component <- component_all[selected]
    beta <- beta_all[selected]
    for (repetition in seq_len(repetitions)) {
      set.seed(860000L + marker_count * 100L + particles * 10L + repetition)
      started <- proc.time()[["elapsed"]]
      update <- sblr:::.sbayesrc_particle_block_step(
        Q, w, component, beta, probability_all[selected, , drop = FALSE],
        gamma, vb, ve, particles = particles, resampling_threshold = 0.5,
        retain_diagnostics = TRUE)
      seconds <- proc.time()[["elapsed"]] - started
      diagnostic <- update$diagnostics
      rows[[length(rows) + 1L]] <- data.frame(
        markers = marker_count, particles = particles, repetition = repetition,
        seconds = seconds, minimum_particle_ess = diagnostic$minimum_ess,
        median_particle_ess = diagnostic$median_ess,
        resampling_count = diagnostic$resampling_count,
        minimum_ancestor_diversity = min(diagnostic$ancestor_diversity),
        final_unique_paths = diagnostic$final_unique_paths,
        selected_reference_path = diagnostic$selected_reference_path,
        allocation_changes = diagnostic$allocation_changes,
        active_count_jump = diagnostic$active_count_jump)
      component <- update$component
      beta <- update$beta
    }
  }
  do.call(rbind, rows)
}

# This is a state-based path-diversity screen, not a posterior chain.  A small
# fixed number of refreshes keeps the pure-R reference tractable while covering
# every preregistered particle count.
result <- rbind(run_size(100L, 5L), run_size(500L, 2L))
utils::write.csv(result, file.path(output_root, "particle_path_screen.csv"),
                 row.names = FALSE)
summary <- aggregate(cbind(seconds, minimum_particle_ess, median_particle_ess,
  resampling_count, minimum_ancestor_diversity, final_unique_paths,
  selected_reference_path, allocation_changes, active_count_jump) ~
    markers + particles, result, mean)
utils::write.csv(summary, file.path(output_root, "particle_path_summary.csv"),
                 row.names = FALSE)
print(summary)
