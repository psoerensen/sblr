#!/usr/bin/env Rscript

# Development-only Study 06 alpha-hierarchy audit. The sibling benchmark
# repository is read-only; compact outputs are written below this repository's
# ignored results/local directory.

parse_option <- function(name, default) {
  prefix <- paste0("--", name, "=")
  value <- commandArgs(trailingOnly = TRUE)
  value <- value[startsWith(value, prefix)]
  if (!length(value)) return(default)
  substring(value[[length(value)]], nchar(prefix) + 1L)
}

sblr_root <- normalizePath(parse_option("sblr-root", "."), winslash = "/",
                           mustWork = TRUE)
bench_root <- normalizePath(parse_option("sblrbench-root", "../sblrbench"),
                            winslash = "/", mustWork = TRUE)
phase <- parse_option("phase", "frozen")
condition <- parse_option("condition", "")
if (!phase %in% c("frozen", "trace-smoke", "dynamic"))
  stop("phase must be frozen, trace-smoke, or dynamic.")
output_root <- file.path(sblr_root, "results", "local",
                         "study06_alpha_hierarchy_audit")
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

isolated_library <- file.path(
  bench_root, "results", "local", "current_benchmark_refresh", "rlib")
if (!dir.exists(isolated_library))
  stop("The Study 06 isolated benchmark library is unavailable.")
.libPaths(c(isolated_library, .libPaths()))
pkgload::load_all(bench_root, quiet = TRUE)
source(file.path(sblr_root, "tools", "study06_alpha_hierarchy_reference.R"),
       local = FALSE)
source(file.path(bench_root, "studies", "06_annotation_models",
                 "power-isolation.R"), local = FALSE)

expected_spec_hash <-
  "241c15afab8fefc571e38e625130de6e4ab58b958c30cff369e26998ee30fa56"
expected_truth_hash <-
  "169d52bef390022a9106d7e61b200493869b40cbb76f1c8d911eebbc80fea1eb"
spec <- read_benchmark_spec(file.path(
  bench_root, "studies", "06_annotation_models", "spec.R"))
if (!identical(benchmark_annotation_spec_hash(spec), expected_spec_hash))
  stop("Study 06 specification hash mismatch.")
glist_cache <- file.path(bench_root, "results", "local",
  "06_annotation_models", "checkpoints", "data", "human_glist.rds")
if (!nzchar(Sys.getenv("SBLR_BENCH_GLIST", "")) && file.exists(glist_cache)) {
  glist <- readRDS(glist_cache)
  benchmark_bedfiles <- glist$bedfiles
  relative <- !grepl("^(?:[A-Za-z]:[/\\\\]|/)", glist$bedfiles)
  glist$bedfiles[relative] <- file.path(bench_root, glist$bedfiles[relative])
  local_glist <- file.path(output_root, "human_glist_absolute_paths.rds")
  saveRDS(glist, local_glist)
  Sys.setenv(SBLR_BENCH_GLIST = local_glist)
}

data <- prepare_prediction_data(spec, output_root)
logic <- .annotation_logic(spec)
annotations <- logic$construct_annotation_design(data$markers$marker_ids, spec)
annotation_truth <- logic$construct_annotation_truth(annotations, spec)
seed_grid <- benchmark_annotation_seeds(spec, "benchmark", "qualification")
seed_row <- seed_grid[
  seed_grid$scenario == "informative_annotations" &
    seed_grid$replicate == 1L & seed_grid$method == "st_bed_bayesrc", ,
  drop = FALSE]
coordinate <- as.list(seed_row[1L, c(
  "scenario", "replicate", "component_seed", "effect_seed", "residual_seed")])
simulation <- logic$simulate_annotation_architecture(
  coordinate, data$scaled$all, data$split$train_rows,
  annotations, annotation_truth, spec)
stats <- benchmark_summary_stats(simulation, data$ld_glist, data$split,
                                 spec$data)
if (exists("benchmark_bedfiles", inherits = FALSE) &&
    length(stats$bed_files) == length(benchmark_bedfiles)) {
  # Preserve the committed semantic identity while the in-memory Glist uses
  # absolute paths solely to read the sibling BED from a read-only checkout.
  stats$bed_files <- benchmark_bedfiles
}
bundle <- list(spec = spec, simulation = simulation, stats = stats,
  marker_truth = simulation$extras$marker_truth)
truth_hash <- benchmark_hash_object(study06_power_truth_identity(data, bundle))
if (!identical(truth_hash, expected_truth_hash)) {
  saveRDS(stats, file.path(output_root, "regenerated_stats_hash_mismatch.rds"))
  utils::write.csv(data.frame(
    truth_hash = truth_hash,
    marker_order_hash = benchmark_hash_object(data$markers$marker_ids),
    phenotype_hash = benchmark_hash_object(simulation$truth$phenotypes),
    effect_hash = benchmark_hash_object(simulation$truth$effects),
    summary_statistics_hash = benchmark_hash_object(stats),
    block_hash = benchmark_hash_object(list(data$block_start, data$marker_panel)),
    component_counts = paste(simulation$extras$component_counts,
                             collapse = "/"),
    realized_h2 = simulation$extras$realized_heritability),
    file.path(output_root, "truth_hash_mismatch.csv"), row.names = FALSE)
  stop("Study 06 shared truth hash mismatch: ", truth_hash)
}
if (!identical(as.numeric(annotations[, 1L]), rep(1, nrow(annotations))))
  stop("Study 06 annotation column one is not an exact all-ones intercept.")

# Truth regeneration intentionally uses the pinned optimized package build
# above. Load the working source only after the exact semantic hash passes.
pkgload::load_all(sblr_root, quiet = TRUE)

write_provenance <- function() {
  data.frame(
    phase = phase,
    sblr_head = system2("git", c("-C", sblr_root, "rev-parse", "HEAD"),
                        stdout = TRUE),
    sblrbench_head = system2("git", c("-C", bench_root, "rev-parse", "HEAD"),
                             stdout = TRUE),
    package_version = as.character(utils::packageVersion("sblr")),
    spec_hash = expected_spec_hash, truth_hash = truth_hash,
    individuals = length(data$sample_ids), training = length(data$split$train_ids),
    validation = length(data$split$test_ids), markers = nrow(annotations),
    stringsAsFactors = FALSE)
}
utils::write.csv(write_provenance(), file.path(output_root, "provenance.csv"),
                 row.names = FALSE)

if (phase == "trace-smoke") {
  methods <- resolve_benchmark_methods(spec)
  names(methods) <- vapply(methods, `[[`, character(1), "id")
  controls <- study06_power_controls(
    spec, data$markers$marker_ids[1L], seed_row$fit_seed,
    seed_row$chain_seeds[[1L]], FALSE, FALSE)
  controls$nit <- 20L
  result <- fit_annotation_method(
    methods[["st_bed_bayesr"]], controls, simulation, stats,
    data$ld_glist, data$split, annotations, annotation_truth,
    data$block_start)
  component <- study06_component_trace(result, data$markers$marker_ids[1L])
  if (!identical(dim(component), c(20L, 4L, 1L)))
    stop("Selected-component trace smoke dimensions are invalid.")
  message("Exact Study 06 selected-component trace smoke passed.")
  quit(save = "no", status = 0L)
}

if (phase == "dynamic") {
  valid_conditions <- c("D1_bed_fixed", "D1_block_fixed",
                        "D2_bed_production", "D2_block_production")
  if (!condition %in% valid_conditions)
    stop("--condition must be one of: ", paste(valid_conditions, collapse = ", "))
  methods <- resolve_benchmark_methods(spec)
  names(methods) <- vapply(methods, `[[`, character(1), "id")
  is_bed <- grepl("_bed_", condition)
  is_fixed <- grepl("_fixed$", condition)
  method_id <- if (is_bed) "st_bed_bayesrc" else
    "st_block_eigen_sbayesrc"
  method <- methods[[method_id]]
  controls <- study06_power_controls(
    spec, data$markers$marker_ids, seed_row$fit_seed,
    seed_row$chain_seeds[[1L]], TRUE, TRUE)
  controls$convergence_control$max_trace_gb <- 4
  controls$convergence_control$allow_large_traces <- TRUE
  controls$.diagnostic_updateSigmaSqAlpha <- !is_fixed
  if (is_fixed) {
    controls$sigmaSqAlpha_init <- rep(1, 3L)
  } else {
    controls$sigmaSqAlpha_a <- 4
    controls$sigmaSqAlpha_b <- 4
  }
  fit_id <- condition
  fit_row <- data.frame(fit_id = fit_id, update_alpha = TRUE,
                        stringsAsFactors = FALSE)
  fit_stats <- stats
  if (!is_bed && length(fit_stats$bed_files) == length(data$ld_glist$bedfiles))
    fit_stats$bed_files <- data$ld_glist$bedfiles
  started <- Sys.time()
  result <- fit_annotation_method(
    method, controls, simulation, fit_stats, data$ld_glist, data$split,
    annotations, annotation_truth, data$block_start)
  elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  component <- study06_component_trace(result, data$markers$marker_ids)
  expected_dimension <- c(spec$qualification$maximum_history,
                          spec$qualification$nchains, nrow(annotations))
  if (!identical(dim(component), expected_dimension))
    stop("Full component-trace dimensions are invalid for ", fit_id, ".")
  occupancy <- study06_occupancy_summary(component, bundle$marker_truth, fit_id)
  draws <- study06_power_required_traces(
    result, component, occupancy, annotations, bundle$marker_truth,
    fit_row, spec)
  selected <- study06_selected_diagnostics(draws, spec)

  probabilities <- extract_marker_probabilities(result)$posterior_inclusion
  if (is.matrix(probabilities)) probabilities <- probabilities[, 1L]
  pip <- as.numeric(probabilities[match(data$markers$marker_ids,
                                       names(probabilities))])
  effect <- extract_marker_effects(result)[, 1L]
  effect <- effect[match(data$markers$marker_ids, names(effect))]
  prediction <- as.numeric(data$scaled$test %*% effect)
  genetic_truth <- as.numeric(simulation$truth$genetic_values[
    data$split$test_ids, 1L])
  phenotype <- as.numeric(simulation$truth$phenotypes[
    data$split$test_ids, 1L])
  metrics <- study06_power_metrics(
    pip, effect, bundle$marker_truth, prediction, genetic_truth, phenotype,
    fit_id)

  alpha <- extract_annotation_coefficient_traces(result, 4L)
  alpha$quantity <- ifelse(alpha$parameter == "alpha",
    paste("alpha", alpha$annotation, alpha$stick, sep = ":"),
    paste("sigmaSqAlpha", alpha$stick, sep = ":"))
  prior <- summarise_drawwise_annotation_prior(
    alpha, annotations, spec$controls$simulation$mixture_var,
    marker_truth = bundle$marker_truth, retain_marker_summary = FALSE)
  if (!identical(prior$status, "ok")) stop(prior$reason, call. = FALSE)

  selected_draws <- selected$rows
  if (is_fixed)
    selected_draws <- selected_draws[
      !grepl("^sigmaSqAlpha:", selected_draws$quantity), , drop = FALSE]
  selected_draws$condition <- condition
  runtime <- data.frame(
    condition = condition, route = if (is_bed) "BED" else "block_eigen_0.995",
    variance_policy = if (is_fixed) "fixed_1" else "sampled_nu0_4_scale0_1",
    iterations = spec$qualification$maximum_history,
    chains = spec$qualification$nchains, seconds = elapsed,
    selected_burnin = selected$selected$burnin,
    selected_retained = selected$selected$retained,
    stringsAsFactors = FALSE)

  keep <- draws$iteration > selected$selected$burnin &
    draws$iteration <= selected$selected$burnin + selected$selected$retained
  scientific <- draws[keep & (grepl("^prior_", draws$quantity) |
    grepl("^occupancy_", draws$quantity) |
    draws$quantity %in% c("effect_variance", "genetic_variance",
      "residual_variance", "heritability", "expected_active_count")), ]
  scientific_summary <- aggregate(value ~ quantity, scientific,
    function(x) c(mean = mean(x), sd = stats::sd(x),
      minimum = min(x), maximum = max(x)))
  scientific_value <- scientific_summary$value
  scientific_summary <- data.frame(
    condition = condition, quantity = scientific_summary$quantity,
    mean = scientific_value[, "mean"], sd = scientific_value[, "sd"],
    minimum = scientific_value[, "minimum"],
    maximum = scientific_value[, "maximum"], row.names = NULL)

  # Geometry summaries use the selected qualification window and preserve
  # iteration/chain pairing without retaining a second raw history.
  alpha_keep <- alpha[alpha$iteration > selected$selected$burnin &
    alpha$iteration <= selected$selected$burnin +
      selected$selected$retained, ]
  alpha_wide <- reshape(alpha_keep[c("iteration", "chain", "quantity", "value")],
    idvar = c("iteration", "chain"), timevar = "quantity", direction = "wide")
  names(alpha_wide) <- sub("^value\\.", "", names(alpha_wide))
  occupancy_wide <- occupancy[
    occupancy$iteration > selected$selected$burnin &
      occupancy$iteration <= selected$selected$burnin +
        selected$selected$retained, ]
  geometry_data <- merge(alpha_wide, occupancy_wide,
    by = c("iteration", "chain"), sort = FALSE)
  expected_active <- prior$draws[
    prior$draws$iteration > selected$selected$burnin &
      prior$draws$iteration <= selected$selected$burnin +
        selected$selected$retained,
    c("iteration", "chain", "expected_active")]
  geometry_data <- merge(geometry_data, expected_active,
    by = c("iteration", "chain"), sort = FALSE)
  core_wide <- reshape(draws[keep & draws$quantity %in%
    c("effect_variance", "genetic_variance", "residual_variance",
      "heritability"), c("iteration", "chain", "quantity", "value")],
    idvar = c("iteration", "chain"), timevar = "quantity", direction = "wide")
  names(core_wide) <- sub("^value\\.", "", names(core_wide))
  geometry_data <- merge(geometry_data, core_wide,
    by = c("iteration", "chain"), sort = FALSE)
  alpha_columns <- grep("^alpha:", names(geometry_data), value = TRUE)
  sigma_columns <- grep("^sigmaSqAlpha:", names(geometry_data), value = TRUE)
  geometry_rows <- list()
  add_correlation <- function(x, y, label, chain) {
    ok <- is.finite(x) & is.finite(y)
    value <- if (sum(ok) > 2L && stats::sd(x[ok]) > 0 &&
      stats::sd(y[ok]) > 0) stats::cor(x[ok], y[ok]) else NA_real_
    data.frame(condition = condition, chain = chain, comparison = label,
               correlation = value, stringsAsFactors = FALSE)
  }
  for (chain in seq_len(spec$qualification$nchains)) {
    z <- geometry_data[geometry_data$chain == chain, , drop = FALSE]
    for (stick in paste0("component_", 0:2, "_stick")) {
      nonintercept <- alpha_columns[grepl(paste0(":", stick, "$"), alpha_columns) &
        !grepl(":Intercept:", alpha_columns)]
      norm <- sqrt(rowSums(z[, nonintercept, drop = FALSE]^2))
      sigma_name <- paste0("sigmaSqAlpha:", stick)
      if (sigma_name %in% sigma_columns) {
        geometry_rows[[length(geometry_rows) + 1L]] <- add_correlation(
          z[[sigma_name]], norm, paste0(stick, ":sigma_vs_alpha_norm"), chain)
        geometry_rows[[length(geometry_rows) + 1L]] <- add_correlation(
          z[[sigma_name]], z$expected_active,
          paste0(stick, ":sigma_vs_expected_active"), chain)
      }
      for (alpha_name in alpha_columns[grepl(paste0(":", stick, "$"),
                                             alpha_columns)]) {
        for (occupancy_name in paste0("traced_component_", 0:3))
          geometry_rows[[length(geometry_rows) + 1L]] <- add_correlation(
            z[[alpha_name]], z[[occupancy_name]],
            paste0(alpha_name, ":vs_", occupancy_name), chain)
        geometry_rows[[length(geometry_rows) + 1L]] <- add_correlation(
          z[[alpha_name]], z$traced_active_count,
          paste0(alpha_name, ":vs_traced_active_count"), chain)
      }
    }
    for (variance_name in c("effect_variance", "genetic_variance",
                            "residual_variance", "heritability"))
      geometry_rows[[length(geometry_rows) + 1L]] <- add_correlation(
        z[[variance_name]], z$traced_active_count,
        paste0(variance_name, ":vs_traced_active_count"), chain)
  }
  geometry <- do.call(rbind, geometry_rows)

  chain_fields <- c("expected_active", "traced_component_0",
    "traced_component_1", "traced_component_2", "traced_component_3",
    "traced_active_count", "causal_active_count", "effect_variance",
    "genetic_variance", "residual_variance", "heritability")
  chain_rows <- list()
  for (chain in seq_len(spec$qualification$nchains)) {
    z <- geometry_data[geometry_data$chain == chain, , drop = FALSE]
    for (name in chain_fields)
      chain_rows[[length(chain_rows) + 1L]] <- data.frame(
        condition = condition, chain = chain, quantity = name,
        mean = mean(z[[name]]), sd = stats::sd(z[[name]]),
        minimum = min(z[[name]]), maximum = max(z[[name]]),
        stringsAsFactors = FALSE)
  }
  chain_summary <- do.call(rbind, chain_rows)
  native_fit <- .benchmark_native_fit(result)
  drift <- native_fit$diagnostics$native$low_rank_residual
  drift_summary <- if (is.null(drift)) {
    data.frame(condition = condition, applicable = FALSE,
      maximum_absolute_drift = NA_real_, rebuild_count = NA_real_)
  } else {
    data.frame(condition = condition, applicable = TRUE,
      maximum_absolute_drift = max(
        drift$low_rank_residual_max_abs_drift, na.rm = TRUE),
      rebuild_count = sum(drift$low_rank_residual_rebuild_count, na.rm = TRUE))
  }

  prefix <- file.path(output_root, condition)
  utils::write.csv(runtime, paste0(prefix, "_runtime.csv"), row.names = FALSE)
  utils::write.csv(selected_draws, paste0(prefix, "_convergence.csv"),
                   row.names = FALSE)
  utils::write.csv(scientific_summary, paste0(prefix, "_summaries.csv"),
                   row.names = FALSE)
  utils::write.csv(metrics, paste0(prefix, "_power_metrics.csv"),
                   row.names = FALSE)
  utils::write.csv(study06_component_recovery(
    pip, bundle$marker_truth, fit_id), paste0(prefix, "_component_recovery.csv"),
    row.names = FALSE)
  utils::write.csv(geometry, paste0(prefix, "_geometry.csv"), row.names = FALSE)
  utils::write.csv(chain_summary, paste0(prefix, "_chain_summary.csv"),
                   row.names = FALSE)
  utils::write.csv(drift_summary, paste0(prefix, "_residual_drift.csv"),
                   row.names = FALSE)
  utils::write.csv(aggregate(value ~ parameter + annotation + stick, alpha_keep,
    function(x) c(mean = mean(x), sd = stats::sd(x), minimum = min(x),
      maximum = max(x))), paste0(prefix, "_alpha_summary.csv"), row.names = FALSE)
  message("Dynamic alpha-hierarchy audit completed: ", condition)
  quit(save = "no", status = 0L)
}

summarise_prior_draws <- function(alpha_array, A) {
  iteration_count <- dim(alpha_array)[1L]
  chain_count <- dim(alpha_array)[2L]
  result <- array(NA_real_, c(iteration_count, chain_count, 4L),
    dimnames = list(NULL, NULL, c("expected_active", "enriched_contrast",
      "continuous_contrast", "null_contrast")))
  enriched <- A[, "enriched_binary"] == 1
  for (chain in seq_len(chain_count)) {
    chunks <- split(seq_len(iteration_count),
                    ceiling(seq_len(iteration_count) / 128L))
    for (index in chunks) {
      coefficient <- t(matrix(
        alpha_array[index, chain, seq_len(ncol(A)), drop = TRUE],
        nrow = length(index), ncol = ncol(A)))
      nonnull <- stats::pnorm(A %*% coefficient)
      result[index, chain, "expected_active"] <- colSums(nonnull)
      result[index, chain, "enriched_contrast"] <-
        colMeans(nonnull[enriched, , drop = FALSE]) -
        colMeans(nonnull[!enriched, , drop = FALSE])
      intercept <- coefficient[match("Intercept", rownames(coefficient)), ]
      for (name in c("continuous_signal", "null_annotation")) {
        value <- coefficient[match(name, rownames(coefficient)), ]
        result[index, chain, sub("_signal|_annotation", "_contrast", name)] <-
          stats::pnorm(intercept + value) - stats::pnorm(intercept - value)
      }
    }
  }
  result
}

run_frozen <- function(id, a, b, update_sigma) {
  gamma <- spec$controls$simulation$mixture_var
  initial <- sblr::make_sbayesrc_alpha_init(
    annotations, gamma = gamma, pi_init = 0.001,
    sigmaSqAlpha_init = rep(1, length(gamma) - 1L))
  prior <- sblr:::.sbayesrc_resolve_intercept_prior(
    initial$component_prob_init)
  native <- prior$native
  if (!update_sigma)
    native <- rbind(native, update_sigmaSqAlpha = rep(0, ncol(native)))
  total <- 12000L
  burnin <- 3000L
  seeds <- c(860101L, 860202L, 860303L, 860404L)
  started <- proc.time()[[3L]]
  chains <- getFromNamespace(".st_bayesrc_frozen_hierarchy_chains", "sblr")(
    annotations, bundle$marker_truth$component_index,
    initial$alpha_init, initial$sigmaSqAlpha_init, native,
    a, b, 1e-12, total, seeds, 4L)
  elapsed <- proc.time()[[3L]] - started
  keep <- seq.int(burnin + 1L, total)
  alpha <- array(NA_real_, c(length(keep), length(chains),
    length(initial$alpha_init)))
  sigma <- array(NA_real_, c(length(keep), length(chains),
    length(initial$sigmaSqAlpha_init)))
  probability <- array(NA_real_, c(length(keep), length(chains),
    ncol(chains[[1L]]$probability_summary)))
  for (chain in seq_along(chains)) {
    alpha[, chain, ] <- chains[[chain]]$alpha[keep, , drop = FALSE]
    sigma[, chain, ] <- chains[[chain]]$sigmaSqAlpha[keep, , drop = FALSE]
    probability[, chain, ] <-
      chains[[chain]]$probability_summary[keep, , drop = FALSE]
  }
  alpha_names <- unlist(lapply(seq_len(ncol(initial$alpha_init)), function(j)
    paste0("alpha:stick", j, ":", rownames(initial$alpha_init))))
  dimnames(alpha)[[3L]] <- alpha_names
  dimnames(sigma)[[3L]] <- paste0("sigmaSqAlpha:stick",
                                  seq_len(dim(sigma)[3L]))
  probability_names <- c("expected_active",
    paste0("mean_component_", 0:3), paste0("min_component_", 0:3),
    paste0("max_component_", 0:3))
  dimnames(probability)[[3L]] <- probability_names
  prior_summary <- summarise_prior_draws(alpha, annotations)
  prior_contrasts <- prior_summary[, , -1L, drop = FALSE]
  monitored <- array(c(alpha, probability, prior_contrasts),
    c(length(keep), length(chains),
      dim(alpha)[3L] + dim(probability)[3L] + dim(prior_contrasts)[3L]))
  dimnames(monitored)[[3L]] <- c(
    dimnames(alpha)[[3L]], probability_names,
    dimnames(prior_contrasts)[[3L]])
  if (update_sigma) {
    combined <- array(c(monitored, sigma),
      c(dim(monitored)[1L], dim(monitored)[2L],
        dim(monitored)[3L] + dim(sigma)[3L]))
    dimnames(combined)[[3L]] <- c(
      dimnames(monitored)[[3L]], dimnames(sigma)[[3L]])
    monitored <- combined
  }
  convergence <- alpha_hierarchy_convergence(monitored)
  convergence$condition <- id
  convergence$stick <- ifelse(grepl("stick[123]", convergence$variable),
    sub(".*stick([123]).*", "\\1", convergence$variable), NA_character_)
  deterministic <- is.finite(convergence$sd) & convergence$sd == 0
  convergence$pass <- deterministic | (convergence$rhat <= 1.01 &
    convergence$ess_bulk >= 400 & convergence$ess_tail >= 400 &
    convergence$relative_mcse <= 0.05)
  list(convergence = convergence,
       runtime = data.frame(condition = id, a = a, b = b,
         update_sigmaSqAlpha = update_sigma, total = total, burnin = burnin,
         retained = length(keep), chains = length(chains), seconds = elapsed,
         pass = all(convergence$pass), stringsAsFactors = FALSE),
       alpha = alpha, sigma = sigma, probability = probability)
}

frozen <- list(
  F1_current = run_frozen("F1_current", 2, 2, TRUE),
  F2_production = run_frozen("F2_production", 4, 4, TRUE),
  F3_fixed = run_frozen("F3_fixed", 2, 2, FALSE))
convergence <- do.call(rbind, lapply(frozen, `[[`, "convergence"))
runtime <- do.call(rbind, lapply(frozen, `[[`, "runtime"))
utils::write.csv(convergence, file.path(output_root, "frozen_convergence.csv"),
                 row.names = FALSE)
utils::write.csv(runtime, file.path(output_root, "frozen_runtime.csv"),
                 row.names = FALSE)
utils::write.csv(simulation$extras$stick_counts,
  file.path(output_root, "frozen_stick_counts.csv"), row.names = FALSE)

true_alpha <- simulation$extras$true_alpha
true_probability <- sblr::sbayesrc_marker_pi(
  annotations, true_alpha, spec$controls$simulation$mixture_var)
f4_error <- max(abs(true_probability - simulation$extras$true_marker_prior))
if (!is.finite(f4_error) || f4_error > 1e-12)
  stop("F4 true-alpha marker probabilities do not reproduce truth.")
utils::write.csv(data.frame(
  condition = "F4_true_alpha", maximum_absolute_probability_error = f4_error,
  pass = f4_error <= 1e-12), file.path(output_root, "frozen_true_alpha.csv"),
  row.names = FALSE)
message("Frozen alpha-hierarchy audit completed: ", output_root)
