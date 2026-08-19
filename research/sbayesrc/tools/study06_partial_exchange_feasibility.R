#!/usr/bin/env Rscript

# Offline-only audit of retained coupling-tempering evidence. This script never
# calls an sblr sampler and deliberately fails closed when replica state needed
# for an exact ratio is absent.

tempered_component_probability <- function(annotation, alpha, baseline, lambda,
                                           floor = 1e-12) {
  eta <- sweep(lambda * annotation %*% alpha, 2L,
               (1 - lambda) * baseline, `+`)
  continuation <- pnorm(eta)
  markers <- nrow(annotation)
  sticks <- ncol(alpha)
  probability <- matrix(0, markers, sticks + 1L)
  remaining <- rep(1, markers)
  for (stick in seq_len(sticks)) {
    probability[, stick] <- remaining * (1 - continuation[, stick])
    remaining <- remaining * continuation[, stick]
  }
  probability[, sticks + 1L] <- remaining
  probability <- pmax(probability, floor)
  probability / rowSums(probability)
}

log_allocation_prior <- function(annotation, alpha, baseline, lambda, component,
                                 floor = 1e-12) {
  probability <- tempered_component_probability(
    annotation, alpha, baseline, lambda, floor)
  sum(log(probability[cbind(seq_len(nrow(probability)), component + 1L)]))
}

log_intercept_prior <- function(intercept, mean, sd = 1) {
  sum(dnorm(intercept, mean, sd, log = TRUE))
}

log_nonintercept_prior <- function(alpha, sigma_sq) {
  stopifnot(nrow(alpha) >= 1L, ncol(alpha) == length(sigma_sq))
  if (nrow(alpha) == 1L) return(0)
  sum(vapply(seq_len(ncol(alpha)), function(stick) {
    sum(dnorm(alpha[-1L, stick], 0, sqrt(sigma_sq[stick]), log = TRUE))
  }, numeric(1L)))
}

log_sigma_prior <- function(sigma_sq, prior_df, prior_scale) {
  # p(sigma^2) proportional to (sigma^2)^-(nu/2+1)
  # exp[-scale/(2 sigma^2)] under the package's scaled-inverse-chi-square form.
  sum(-(prior_df / 2 + 1) * log(sigma_sq) -
        prior_scale / (2 * sigma_sq))
}

log_hierarchy_prior <- function(alpha, sigma_sq, intercept_mean,
                                prior_df, prior_scale) {
  log_intercept_prior(alpha[1L, ], intercept_mean) +
    log_nonintercept_prior(alpha, sigma_sq) +
    log_sigma_prior(sigma_sq, prior_df, prior_scale)
}

complete_exchange_ratio <- function(annotation, baseline, lambda_a, lambda_b,
                                    component_a, alpha_a,
                                    component_b, alpha_b, floor = 1e-12) {
  log_allocation_prior(annotation, alpha_b, baseline, lambda_a, component_b, floor) +
    log_allocation_prior(annotation, alpha_a, baseline, lambda_b, component_a, floor) -
    log_allocation_prior(annotation, alpha_a, baseline, lambda_a, component_a, floor) -
    log_allocation_prior(annotation, alpha_b, baseline, lambda_b, component_b, floor)
}

alpha_sigma_exchange_ratio <- function(annotation, baseline, lambda_a, lambda_b,
                                       component_a, alpha_a, sigma_a,
                                       component_b, alpha_b, sigma_b,
                                       prior_df, prior_scale, floor = 1e-12) {
  allocation <-
    log_allocation_prior(annotation, alpha_b, baseline, lambda_a, component_a, floor) +
    log_allocation_prior(annotation, alpha_a, baseline, lambda_b, component_b, floor) -
    log_allocation_prior(annotation, alpha_a, baseline, lambda_a, component_a, floor) -
    log_allocation_prior(annotation, alpha_b, baseline, lambda_b, component_b, floor)
  hierarchy <-
    log_hierarchy_prior(alpha_b, sigma_b, baseline, prior_df, prior_scale) +
    log_hierarchy_prior(alpha_a, sigma_a, baseline, prior_df, prior_scale) -
    log_hierarchy_prior(alpha_a, sigma_a, baseline, prior_df, prior_scale) -
    log_hierarchy_prior(alpha_b, sigma_b, baseline, prior_df, prior_scale)
  allocation + hierarchy
}

alpha_only_exchange_ratio <- function(annotation, baseline, lambda_a, lambda_b,
                                      component_a, alpha_a, sigma_a,
                                      component_b, alpha_b, sigma_b,
                                      floor = 1e-12) {
  allocation <-
    log_allocation_prior(annotation, alpha_b, baseline, lambda_a, component_a, floor) +
    log_allocation_prior(annotation, alpha_a, baseline, lambda_b, component_b, floor) -
    log_allocation_prior(annotation, alpha_a, baseline, lambda_a, component_a, floor) -
    log_allocation_prior(annotation, alpha_b, baseline, lambda_b, component_b, floor)
  # Intercept-prior products cancel; non-intercept coefficients are evaluated
  # under the destination replica's retained sigmaSqAlpha.
  hierarchy <- log_nonintercept_prior(alpha_b, sigma_a) +
    log_nonintercept_prior(alpha_a, sigma_b) -
    log_nonintercept_prior(alpha_a, sigma_a) -
    log_nonintercept_prior(alpha_b, sigma_b)
  allocation + hierarchy
}

sha256_file <- function(path) {
  if (!requireNamespace("digest", quietly = TRUE))
    stop("The digest package is required for SHA-256 provenance.")
  unname(digest::digest(path, algo = "sha256", file = TRUE, serialize = FALSE))
}

main <- function(root = ".") {
  root <- normalizePath(root, winslash = "/", mustWork = TRUE)
  evidence_root <- file.path(root, "results", "local",
    "study06_bed_coupling_tempering_screen")
  output_root <- file.path(evidence_root, "partial_exchange_feasibility")
  dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
  fit_path <- file.path(evidence_root, "study06_bed_tempered_fit.rds")
  if (!file.exists(fit_path)) stop("Retained tempered fit is unavailable.")
  fit <- readRDS(fit_path)
  native <- fit$result$native_fit
  chains <- native$chains
  if (length(chains) != 4L) stop("Expected four retained ensembles.")

  expected_spec <-
    "241c15afab8fefc571e38e625130de6e4ab58b958c30cff369e26998ee30fa56"
  expected_truth <-
    "169d52bef390022a9106d7e61b200493869b40cbb76f1c8d911eebbc80fea1eb"
  provenance <- read.csv(file.path(evidence_root, "provenance.csv"),
                         check.names = FALSE)
  if (!identical(provenance$spec_hash[[1L]], expected_spec) ||
      !identical(provenance$truth_hash[[1L]], expected_truth))
    stop("Retained Study 06 identity mismatch.")

  field_rows <- list()
  for (ensemble in seq_along(chains)) {
    chain <- chains[[ensemble]]
    coupling <- chain$coupling_tempering
    target <- chain$convergence_trace
    retained <- list(
      target_component = target$component,
      target_effect = target$b,
      target_indicator = target$d,
      target_alpha = target$alpha,
      target_sigmaSqAlpha = target$sigmaSqAlpha,
      all_slot_identity = coupling$replica_identity,
      all_slot_active_count = coupling$active_count,
      all_slot_expected_active = coupling$expected_active_count,
      exchange_record = coupling$swap)
    for (field in names(retained)) {
      value <- retained[[field]]
      field_rows[[length(field_rows) + 1L]] <- data.frame(
        ensemble = ensemble, field = field,
        dimension = paste(dim(value) %||% length(value), collapse = "x"),
        retained = TRUE, stringsAsFactors = FALSE)
    }
    for (field in c(
      "slot0_component", "slot05_component", "slot0_effect", "slot05_effect",
      "slot0_residual", "slot05_residual", "slot0_alpha", "slot05_alpha",
      "slot0_sigmaSqAlpha", "slot05_sigmaSqAlpha", "slot0_marker_probability",
      "slot05_marker_probability", "slot_rng_state")) {
      field_rows[[length(field_rows) + 1L]] <- data.frame(
        ensemble = ensemble, field = field, dimension = NA_character_,
        retained = FALSE, stringsAsFactors = FALSE)
    }
  }
  inventory <- do.call(rbind, field_rows)
  write.csv(inventory, file.path(output_root, "checkpoint_inventory.csv"),
            row.names = FALSE)

  # Exact scalar compatibility that can be evaluated without reconstructing
  # missing replica allocations or hierarchy states.
  compatibility <- list()
  for (ensemble in seq_along(chains)) {
    coupling <- chains[[ensemble]]$coupling_tempering
    swap <- coupling$swap
    for (attempt in seq_len(nrow(swap))) {
      iteration <- as.integer(swap[attempt, 1L])
      lower <- as.integer(swap[attempt, 2L]) + 1L
      upper <- lower + 1L
      compatibility[[length(compatibility) + 1L]] <- data.frame(
        ensemble = ensemble, attempt = attempt, iteration = iteration,
        lower_lambda = coupling$lambda[lower], upper_lambda = coupling$lambda[upper],
        lower_active = coupling$active_count[iteration, lower],
        upper_active = coupling$active_count[iteration, upper],
        active_difference = coupling$active_count[iteration, upper] -
          coupling$active_count[iteration, lower],
        lower_expected_active = coupling$expected_active_count[iteration, lower],
        upper_expected_active = coupling$expected_active_count[iteration, upper],
        expected_active_difference =
          coupling$expected_active_count[iteration, upper] -
          coupling$expected_active_count[iteration, lower],
        recorded_log_ratio = swap[attempt, 5L],
        recorded_probability = swap[attempt, 4L])
    }
  }
  compatibility <- do.call(rbind, compatibility)
  write.csv(compatibility, file.path(output_root, "retained_count_compatibility.csv"),
            row.names = FALSE)
  count_summary <- aggregate(
    cbind(active_difference, expected_active_difference, recorded_log_ratio) ~
      lower_lambda + upper_lambda,
    compatibility,
    function(x) c(mean = mean(x), sd = sd(x), median = median(x),
                  minimum = min(x), maximum = max(x)))
  write.csv(count_summary, file.path(output_root, "retained_count_summary.csv"),
            row.names = FALSE)
  count_correlation <- do.call(rbind, lapply(
    split(compatibility,
      paste(compatibility$lower_lambda, compatibility$upper_lambda, sep = "_")),
    function(x) data.frame(
      lower_lambda = x$lower_lambda[1L], upper_lambda = x$upper_lambda[1L],
      correlation_abs_active_difference_log_ratio =
        cor(abs(x$active_difference), x$recorded_log_ratio),
      correlation_abs_expected_difference_log_ratio =
        cor(abs(x$expected_active_difference), x$recorded_log_ratio))))
  write.csv(count_correlation,
    file.path(output_root, "retained_count_correlations.csv"), row.names = FALSE)

  feasibility <- data.frame(
    proposal = c("complete", paste0("P", 1:9),
      "mean_sparsity_matched", "mean_preserving_path", "marker_count_scaling"),
    exact_from_retained_state = c(FALSE, rep(FALSE, 9L), FALSE, FALSE, FALSE),
    classification = c(
      "unavailable", rep("unavailable", 8L), "invalid_or_unavailable",
      "unavailable", "unavailable", "unavailable"),
    missing_state = c(
      rep("slot-0 and slot-0.5 allocations and alpha at exchange times", 1L),
      rep("slot-0 and slot-0.5 allocations, alpha, and sigmaSqAlpha at exchange times", 7L),
      "slot-0 and slot-0.5 complete marker/effect/residual/variance state",
      "slot-0 and slot-0.5 allocations/effects; spike-support consistency",
      "source allocations and alpha for both adjacent slots",
      "source allocations and alpha for both adjacent slots",
      "markerwise exchange contributions for both adjacent slots"),
    stringsAsFactors = FALSE)
  write.csv(feasibility, file.path(output_root, "proposal_feasibility.csv"),
            row.names = FALSE)

  baseline_files <- c(
    "R/RcppExports.R", "R/sparse_ld_bed_helper.R", "src/RcppExports.cpp",
    "src/blr_bed_bayesrc_core_impl.h", "src/blr_bed_bayesrc_types.h",
    "src/blr_mt_bayesrc_types.h", "src/st_bayesrc_annotation_prior.h",
    "src/stblr_cpg_omp_bed_marker_scheduled_chains_bayesrc.cpp",
    "docs/dev/study06_bed_coupling_tempering_decision.json",
    "docs/dev/study06_bed_coupling_tempering_screen.md",
    "tests/testthat/test-bayesrc-coupling-tempering.R",
    "research/sbayesrc/tools/coupling_tempering_tiny_reference.R",
    "research/sbayesrc/tools/study06_bed_coupling_tempering_screen.R")
  evidence_files <- file.path("results", "local",
    "study06_bed_coupling_tempering_screen", c(
      "pre_partial_exchange_audit.patch", "provenance.csv",
      "study06_bed_coupling_tempering_decision.local.json",
      "study06_bed_tempered_component_trace.rds", "study06_bed_tempered_fit.rds",
      "swap_summary.csv", "active_count_summary.csv", "target_convergence.csv",
      "tiny/tiny_reference_result.rds"))
  manifest_paths <- c(baseline_files, evidence_files)
  absolute <- file.path(root, manifest_paths)
  if (any(!file.exists(absolute)))
    stop("A required provenance file is missing: ",
         paste(manifest_paths[!file.exists(absolute)], collapse = ", "))
  manifest <- data.frame(
    path = manifest_paths,
    bytes = file.info(absolute)$size,
    sha256 = vapply(absolute, sha256_file, character(1L)),
    stringsAsFactors = FALSE)
  write.csv(manifest, file.path(output_root, "preserved_file_hashes.csv"),
            row.names = FALSE)

  decision <- list(
    schema_version = 1L,
    source_head = "8908267a68a46267fcccb910850b4f6380bfa978",
    specification_hash = expected_spec,
    truth_hash = expected_truth,
    exact_complete_ratio_reproduced = FALSE,
    exact_complete_ratio_status = "unavailable_missing_non_target_replica_state",
    retained_exact_fields = c("slot active count", "slot expected active count",
      "replica identity", "recorded exchange log ratio", "target-level histories"),
    missing_minimum_state = c("component allocation by slot at every exchange attempt",
      "alpha and sigmaSqAlpha by slot at every exchange attempt"),
    required_for_P8 = c("effects", "residual", "effect variance",
      "residual variance", "other persistent variance state by slot"),
    primary_decision = "F6",
    interpretation = "retained state is insufficient for exact decomposition or partial-exchange audit",
    mcmc_ran = FALSE,
    native_sampler_changed = FALSE,
    sibling_modified = FALSE)
  jsonlite::write_json(decision,
    file.path(output_root, "partial_exchange_decision.local.json"),
    pretty = TRUE, auto_unbox = TRUE)
  invisible(list(inventory = inventory, compatibility = compatibility,
                 feasibility = feasibility, manifest = manifest,
                 decision = decision))
}

`%||%` <- function(x, y) if (is.null(x)) y else x

if (sys.nframe() == 0L) main()
