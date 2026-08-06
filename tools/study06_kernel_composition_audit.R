#!/usr/bin/env Rscript

# Development-only, resumable Study 06 kernel-composition audit. Raw marker
# traces are written only below results/local/ and are never qualification data.

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
phase <- parse_option("phase", "smoke")
condition <- parse_option("condition", "")
if (!phase %in% c("smoke", "fit", "aggregate"))
  stop("phase must be smoke, fit, or aggregate.")

output_root <- file.path(sblr_root, "results", "local",
                         "study06_kernel_composition_audit")
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

schedule_registry <- data.frame(
  schedule = c("S1", "H5", "H20", "A5", "A20"),
  allocation_updates = c(1L, 1L, 1L, 5L, 20L),
  annotation_updates = c(1L, 5L, 20L, 1L, 1L),
  stringsAsFactors = FALSE
)
fit_registry <- merge(
  data.frame(route = c("bed", "block_eigen"), stringsAsFactors = FALSE),
  schedule_registry, by = NULL, sort = FALSE)
fit_registry$fit_id <- paste(fit_registry$route, fit_registry$schedule, sep = "_")
fit_registry <- fit_registry[match(
  as.vector(outer(c("bed", "block_eigen"), schedule_registry$schedule,
                  paste, sep = "_")), fit_registry$fit_id), ]
rownames(fit_registry) <- NULL
utils::write.csv(fit_registry, file.path(output_root, "fit_registry.csv"),
                 row.names = FALSE)

isolated_library <- file.path(
  bench_root, "results", "local", "current_benchmark_refresh", "rlib")
if (!dir.exists(isolated_library))
  stop("The Study 06 isolated benchmark library is unavailable.")
.libPaths(c(isolated_library, .libPaths()))
pkgload::load_all(bench_root, quiet = TRUE)
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
    length(stats$bed_files) == length(benchmark_bedfiles))
  stats$bed_files <- benchmark_bedfiles
bundle <- list(spec = spec, simulation = simulation, stats = stats,
               marker_truth = simulation$extras$marker_truth)
truth_hash <- benchmark_hash_object(study06_power_truth_identity(data, bundle))
if (!identical(truth_hash, expected_truth_hash))
  stop("Study 06 shared truth hash mismatch: ", truth_hash)
if (!identical(as.numeric(annotations[, 1L]), rep(1, nrow(annotations))))
  stop("Study 06 annotation column one is not an exact all-ones intercept.")

# Truth is established using the pinned benchmark library before loading the
# working package source whose internal diagnostic schedule is under test.
pkgload::load_all(sblr_root, quiet = TRUE)

methods <- resolve_benchmark_methods(spec)
names(methods) <- vapply(methods, `[[`, character(1), "id")

write_provenance <- function() {
  utils::write.csv(data.frame(
    sblr_head = system2("git", c("-C", sblr_root, "rev-parse", "HEAD"),
                        stdout = TRUE),
    sblrbench_head = system2("git", c("-C", bench_root, "rev-parse", "HEAD"),
                             stdout = TRUE),
    package_version = as.character(utils::packageVersion("sblr")),
    spec_hash = expected_spec_hash, truth_hash = truth_hash,
    fit_seed = 701020L, chain_seeds = "701121/701222/701323/701424",
    individuals = length(data$sample_ids), training = length(data$split$train_ids),
    validation = length(data$split$test_ids), markers = nrow(annotations),
    stringsAsFactors = FALSE), file.path(output_root, "provenance.csv"),
    row.names = FALSE)
}
write_provenance()

make_controls <- function(row, smoke = FALSE) {
  markers <- if (smoke) data$markers$marker_ids[seq_len(12L)] else
    data$markers$marker_ids
  controls <- study06_power_controls(
    spec, markers, 701020L, c(701121L, 701222L, 701323L, 701424L), TRUE, TRUE)
  controls$.diagnostic_allocation_updates_per_cycle <- row$allocation_updates
  controls$.diagnostic_annotation_updates_per_cycle <- row$annotation_updates
  controls$convergence_control$selected_marker_quantities <-
    if (smoke) c("b", "d", "component") else c("b", "d", "component")
  controls$convergence_control$allow_large_traces <- TRUE
  controls$convergence_control$max_trace_gb <- if (smoke) 1 else 8
  if (smoke) {
    controls$nit <- 20L
    controls$nburn <- 5L
  }
  controls
}

fit_one <- function(row, smoke = FALSE) {
  method_id <- if (row$route == "bed") "st_bed_bayesrc" else
    "st_block_eigen_sbayesrc"
  fit_stats <- stats
  if (row$route == "block_eigen" &&
      length(fit_stats$bed_files) == length(data$ld_glist$bedfiles))
    fit_stats$bed_files <- data$ld_glist$bedfiles
  controls <- make_controls(row, smoke)
  started <- Sys.time()
  result <- fit_annotation_method(
    methods[[method_id]], controls, simulation, fit_stats, data$ld_glist,
    data$split, annotations, annotation_truth, data$block_start)
  list(result = result,
       seconds = as.numeric(difftime(Sys.time(), started, units = "secs")),
       controls = controls)
}

check_fit <- function(result, marker_ids) {
  native <- .benchmark_native_fit(result)
  alpha <- extract_annotation_coefficient_traces(result, ncol(annotations))
  if (any(!is.finite(alpha$value))) stop("Non-finite annotation trace.")
  prior <- native$component_probabilities
  if (is.list(prior)) prior <- prior[[1L]]
  prior <- as.matrix(prior)
  if (any(!is.finite(prior)) || any(prior < 0) ||
      max(abs(rowSums(prior) - 1)) > 1e-10)
    stop("Invalid marker-specific component probabilities.")
  component <- study06_component_trace(result, marker_ids)
  if (any(component < 0 | component > 3 | component != floor(component)))
    stop("Invalid component trace.")
  drift <- native$diagnostics$native$low_rank_residual
  if (!is.null(drift) && any(!is.finite(
      drift$low_rank_residual_max_abs_drift)))
    stop("Invalid low-rank residual drift diagnostic.")
  invisible(component)
}

if (phase == "smoke") {
  rows <- vector("list", nrow(fit_registry))
  for (index in seq_len(nrow(fit_registry))) {
    row <- fit_registry[index, ]
    out <- fit_one(row, smoke = TRUE)
    component <- check_fit(
      out$result, data$markers$marker_ids[seq_len(12L)])
    rows[[index]] <- data.frame(
      fit_id = row$fit_id, route = row$route, schedule = row$schedule,
      seconds = out$seconds, component_iterations = dim(component)[1L],
      component_chains = dim(component)[2L],
      component_markers = dim(component)[3L], passed = TRUE)
  }
  utils::write.csv(do.call(rbind, rows), file.path(output_root, "smoke.csv"),
                   row.names = FALSE)
  message("All ten kernel-composition smoke fits passed.")
  quit(save = "no", status = 0L)
}

summarise_bands <- function(values, fit_id, quantity) {
  labels <- c("0-49", "50-74", "75-99", "100-149", "150-199", "200+")
  rows <- list()
  transitions <- list()
  for (chain in seq_len(ncol(values))) {
    band <- cut(values[, chain], breaks = c(-Inf, 49, 74, 99, 149, 199, Inf),
                labels = labels, right = TRUE)
    tab <- table(factor(band, levels = labels)) / length(band)
    run <- rle(as.character(band))
    matrix_count <- table(
      factor(head(band, -1L), levels = labels),
      factor(tail(band, -1L), levels = labels))
    rows[[chain]] <- data.frame(
      fit_id = fit_id, quantity = quantity, chain = chain, band = labels,
      proportion = as.numeric(tab), band_transitions = sum(diff(as.integer(band)) != 0),
      maximum_run = max(run$lengths),
      low_to_high = sum(head(values[, chain], -1L) <= 74 &
                        tail(values[, chain], -1L) >= 150),
      high_to_low = sum(head(values[, chain], -1L) >= 150 &
                        tail(values[, chain], -1L) <= 74),
      stringsAsFactors = FALSE)
    transitions[[chain]] <- data.frame(
      fit_id = fit_id, quantity = quantity, chain = chain,
      from = rep(labels, each = length(labels)), to = rep(labels, length(labels)),
      count = as.numeric(matrix_count), stringsAsFactors = FALSE)
  }
  list(summary = do.call(rbind, rows), matrix = do.call(rbind, transitions))
}

component_diagnostics <- function(component, marker_truth, A, fit_id) {
  keep <- 3001:9000
  component <- component[keep, , , drop = FALSE]
  nchains <- dim(component)[2L]
  causal <- marker_truth$true_nonnull
  enriched <- A[, 2L] > 0.5
  null_annot <- A[, 4L] > 0.5
  continuous_group <- cut(A[, 3L], stats::quantile(A[, 3L], 0:4 / 4),
                          include.lowest = TRUE, labels = FALSE)
  traces <- vector("list", nchains)
  changes <- vector("list", nchains)
  representative_causal <- which(causal)[1:4]
  representative_noncausal <- which(!causal)[1:4]
  for (chain in seq_len(nchains)) {
    state <- component[, chain, , drop = TRUE]
    counts <- vapply(0:3, function(k) rowSums(state == k), numeric(nrow(state)))
    causal_counts <- vapply(0:3, function(k)
      rowSums(state[, causal, drop = FALSE] == k), numeric(nrow(state)))
    traces[[chain]] <- data.frame(
      fit_id = fit_id, iteration = keep, chain = chain,
      component_0 = counts[, 1L], component_1 = counts[, 2L],
      component_2 = counts[, 3L], component_3 = counts[, 4L],
      active_count = rowSums(counts[, -1L, drop = FALSE]),
      active_enriched = rowSums(state[, enriched, drop = FALSE] > 0),
      active_null_annotation = rowSums(state[, null_annot, drop = FALSE] > 0),
      active_continuous_q1 = rowSums(state[, continuous_group == 1L, drop = FALSE] > 0),
      active_continuous_q4 = rowSums(state[, continuous_group == 4L, drop = FALSE] > 0),
      causal_active = rowSums(causal_counts[, -1L, drop = FALSE]),
      causal_component_1 = causal_counts[, 2L],
      causal_component_2 = causal_counts[, 3L],
      causal_component_3 = causal_counts[, 4L])
    changed <- rowSums(state[-1L, , drop = FALSE] !=
                       state[-nrow(state), , drop = FALSE])
    entering <- rowSums(state[-1L, , drop = FALSE] > 0 &
                        state[-nrow(state), , drop = FALSE] == 0)
    leaving <- rowSums(state[-1L, , drop = FALSE] == 0 &
                       state[-nrow(state), , drop = FALSE] > 0)
    count_change <- abs(counts[-1L, , drop = FALSE] -
                        counts[-nrow(counts), , drop = FALSE])
    changes[[chain]] <- data.frame(
      fit_id = fit_id, iteration = keep[-1L], chain = chain,
      markers_changing = changed, entering_null = leaving, leaving_null = entering,
      largest_coordinated_occupancy_change = apply(count_change, 1L, max))
  }
  list(trace = do.call(rbind, traces), changes = do.call(rbind, changes),
       representative = c(representative_causal, representative_noncausal))
}

stick_prior_draws <- function(alpha, A, fit_id) {
  alpha <- alpha[alpha$iteration > 3000L & alpha$iteration <= 9000L, ]
  alpha_only <- alpha[alpha$parameter == "alpha", ]
  rows <- list()
  groups <- list(enriched = A[, 2L] > 0.5,
                 continuous_high = A[, 3L] >= stats::quantile(A[, 3L], 0.75),
                 continuous_low = A[, 3L] <= stats::quantile(A[, 3L], 0.25),
                 null = A[, 4L] > 0.5)
  for (chain in sort(unique(alpha_only$chain))) {
    z <- alpha_only[alpha_only$chain == chain, ]
    for (iteration in sort(unique(z$iteration))) {
      zi <- z[z$iteration == iteration, ]
      matrix_alpha <- matrix(NA_real_, nrow = ncol(A), ncol = 3L)
      for (stick in seq_len(3L)) {
        zs <- zi[zi$stick == paste0("component_", stick - 1L, "_stick"), ]
        matrix_alpha[, stick] <- zs$value[match(colnames(A), zs$annotation)]
      }
      p <- stats::pnorm(A %*% matrix_alpha)
      rows[[length(rows) + 1L]] <- data.frame(
        fit_id = fit_id, iteration = iteration, chain = chain,
        mean_continuation_stick1 = mean(p[, 1L]),
        mean_continuation_stick2 = mean(p[, 2L]),
        mean_continuation_stick3 = mean(p[, 3L]),
        enriched_contrast_stick1 = mean(p[groups$enriched, 1L]) -
          mean(p[!groups$enriched, 1L]),
        enriched_contrast_stick2 = mean(p[groups$enriched, 2L]) -
          mean(p[!groups$enriched, 2L]),
        enriched_contrast_stick3 = mean(p[groups$enriched, 3L]) -
          mean(p[!groups$enriched, 3L]),
        continuous_contrast_stick1 = mean(p[groups$continuous_high, 1L]) -
          mean(p[groups$continuous_low, 1L]),
        continuous_contrast_stick2 = mean(p[groups$continuous_high, 2L]) -
          mean(p[groups$continuous_low, 2L]),
        continuous_contrast_stick3 = mean(p[groups$continuous_high, 3L]) -
          mean(p[groups$continuous_low, 3L]),
        null_contrast_stick1 = mean(p[groups$null, 1L]) - mean(p[!groups$null, 1L]),
        null_contrast_stick2 = mean(p[groups$null, 2L]) - mean(p[!groups$null, 2L]),
        null_contrast_stick3 = mean(p[groups$null, 3L]) - mean(p[!groups$null, 3L]),
        expected_active_count = sum(p[, 1L]),
        minimum_nonnull_prior = min(p[, 1L]),
        maximum_nonnull_prior = max(p[, 1L]))
    }
  }
  do.call(rbind, rows)
}

if (phase == "fit") {
  if (!condition %in% fit_registry$fit_id)
    stop("Unknown --condition: ", condition)
  row <- fit_registry[fit_registry$fit_id == condition, ]
  prefix <- file.path(output_root, condition)
  result_checkpoint <- paste0(prefix, "_fit_result.rds")
  if (file.exists(result_checkpoint)) {
    out <- readRDS(result_checkpoint)
  } else {
    out <- fit_one(row, smoke = FALSE)
    saveRDS(out, result_checkpoint, compress = FALSE)
  }
  component <- check_fit(out$result, data$markers$marker_ids)
  expected_dimension <- c(9000L, 4L, 1500L)
  if (!identical(dim(component), expected_dimension))
    stop("Full component trace has an invalid dimension.")
  saveRDS(component, paste0(prefix, "_component_trace.rds"), compress = FALSE)

  occupancy <- component_diagnostics(
    component, bundle$marker_truth, annotations, condition)
  alpha <- extract_annotation_coefficient_traces(out$result, ncol(annotations))
  alpha$quantity <- ifelse(alpha$parameter == "alpha",
    paste("alpha", alpha$annotation, alpha$stick, sep = ":"),
    paste("sigmaSqAlpha", alpha$stick, sep = ":"))
  prior <- summarise_drawwise_annotation_prior(
    alpha, annotations, spec$controls$simulation$mixture_var,
    marker_truth = bundle$marker_truth, retain_marker_summary = FALSE)
  if (!identical(prior$status, "ok")) stop(prior$reason, call. = FALSE)
  stick <- stick_prior_draws(alpha, annotations, condition)
  expected <- reshape(
    stick[c("iteration", "chain", "expected_active_count")],
    idvar = "iteration", timevar = "chain", direction = "wide")
  expected_matrix <- as.matrix(expected[, -1L, drop = FALSE])
  active <- reshape(
    occupancy$trace[c("iteration", "chain", "active_count")],
    idvar = "iteration", timevar = "chain", direction = "wide")
  active_matrix <- as.matrix(active[, -1L, drop = FALSE])
  active_bands <- summarise_bands(active_matrix, condition, "realized_active_count")
  expected_bands <- summarise_bands(expected_matrix, condition, "expected_active_count")

  fit_row <- data.frame(fit_id = condition, update_alpha = TRUE)
  basic_occupancy <- study06_occupancy_summary(
    component, bundle$marker_truth, condition)
  draws <- study06_power_required_traces(
    out$result, component, basic_occupancy, annotations,
    bundle$marker_truth, fit_row, spec)
  selected <- study06_selected_diagnostics(draws, spec)

  probabilities <- extract_marker_probabilities(out$result)$posterior_inclusion
  if (is.matrix(probabilities)) probabilities <- probabilities[, 1L]
  pip <- as.numeric(probabilities[match(data$markers$marker_ids,
                                       names(probabilities))])
  effect <- extract_marker_effects(out$result)[, 1L]
  effect <- effect[match(data$markers$marker_ids, names(effect))]
  prediction <- as.numeric(data$scaled$test %*% effect)
  genetic_truth <- as.numeric(simulation$truth$genetic_values[
    data$split$test_ids, 1L])
  phenotype <- as.numeric(simulation$truth$phenotypes[
    data$split$test_ids, 1L])
  metrics <- study06_power_metrics(
    pip, effect, bundle$marker_truth, prediction, genetic_truth, phenotype,
    condition)

  native <- .benchmark_native_fit(out$result)
  trace_quantities <- native$convergence_traces$quantities
  representative <- occupancy$representative
  representative_id <- data$markers$marker_ids[representative]
  representative_index <- which(trace_quantities$marker_id %in% representative_id &
    trace_quantities$group %in% c("selected_b", "selected_d"))
  representative_trace <- list(
    quantities = trace_quantities[representative_index, , drop = FALSE],
    values = native$convergence_traces$values[, , representative_index, drop = FALSE])
  saveRDS(representative_trace, paste0(prefix, "_representative_trace.rds"))

  drift <- native$diagnostics$native$low_rank_residual
  drift_row <- data.frame(
    fit_id = condition, applicable = !is.null(drift),
    maximum_absolute_drift = if (is.null(drift)) NA_real_ else
      max(drift$low_rank_residual_max_abs_drift),
    rebuild_count = if (is.null(drift)) NA_real_ else
      sum(drift$low_rank_residual_rebuild_count))
  runtime <- data.frame(
    fit_id = condition, route = row$route, schedule = row$schedule,
    allocation_updates = row$allocation_updates,
    annotation_updates = row$annotation_updates,
    recorded_iterations = 9000L, retained_iterations = 6000L, chains = 4L,
    kernel_A_applications = 9000L * 4L * row$allocation_updates,
    kernel_H_applications = 9000L * 4L * row$annotation_updates,
    seconds = out$seconds)
  convergence <- selected$rows
  convergence$fit_id <- condition
  convergence$ess_bulk_per_second <- convergence$ess_bulk / out$seconds
  convergence$ess_tail_per_second <- convergence$ess_tail / out$seconds

  utils::write.csv(runtime, paste0(prefix, "_runtime.csv"), row.names = FALSE)
  utils::write.csv(convergence, paste0(prefix, "_convergence.csv"), row.names = FALSE)
  utils::write.csv(occupancy$trace, paste0(prefix, "_occupancy.csv"), row.names = FALSE)
  utils::write.csv(occupancy$changes, paste0(prefix, "_allocation_changes.csv"), row.names = FALSE)
  utils::write.csv(rbind(active_bands$summary, expected_bands$summary),
                   paste0(prefix, "_regime_summary.csv"), row.names = FALSE)
  utils::write.csv(rbind(active_bands$matrix, expected_bands$matrix),
                   paste0(prefix, "_regime_transitions.csv"), row.names = FALSE)
  utils::write.csv(stick, paste0(prefix, "_stick_prior_draws.csv"), row.names = FALSE)
  utils::write.csv(metrics, paste0(prefix, "_power_metrics.csv"), row.names = FALSE)
  utils::write.csv(drift_row, paste0(prefix, "_residual_drift.csv"), row.names = FALSE)
  message("Completed kernel-composition fit: ", condition)
  quit(save = "no", status = 0L)
}

if (phase == "aggregate") {
  suffixes <- c("runtime", "convergence", "regime_summary",
                "regime_transitions", "power_metrics", "residual_drift")
  for (suffix in suffixes) {
    paths <- file.path(output_root, paste0(fit_registry$fit_id, "_", suffix, ".csv"))
    if (any(!file.exists(paths)))
      stop("Missing fit output for aggregate suffix: ", suffix)
    combined <- do.call(rbind, lapply(paths, utils::read.csv,
                                     check.names = FALSE))
    utils::write.csv(combined, file.path(output_root, paste0("all_", suffix, ".csv")),
                     row.names = FALSE)
  }
  variance_rows <- list()
  annotation_range_rows <- list()
  marker_agreement_rows <- list()
  representative_rows <- list()
  eigen_rows <- list()
  representative_index <- c(which(bundle$marker_truth$true_nonnull)[1:4],
                            which(!bundle$marker_truth$true_nonnull)[1:4])
  for (fit_id in fit_registry$fit_id) {
    checkpoint <- readRDS(file.path(output_root, paste0(fit_id, "_fit_result.rds")))
    native <- .benchmark_native_fit(checkpoint$result)
    chains <- native$chains
    for (chain in seq_along(chains)) {
      z <- chains[[chain]]
      keep <- 3001:9000
      for (quantity in c("vbs", "vgs", "ves"))
        variance_rows[[length(variance_rows) + 1L]] <- data.frame(
          fit_id = fit_id, chain = chain, quantity = quantity,
          mean = mean(z[[quantity]][keep]), sd = stats::sd(z[[quantity]][keep]),
          minimum = min(z[[quantity]][keep]), maximum = max(z[[quantity]][keep]))
      h2 <- z$vgs[keep] / (z$vgs[keep] + z$ves[keep])
      variance_rows[[length(variance_rows) + 1L]] <- data.frame(
        fit_id = fit_id, chain = chain, quantity = "heritability",
        mean = mean(h2), sd = stats::sd(h2), minimum = min(h2), maximum = max(h2))
      for (marker in representative_index) {
        representative_rows[[length(representative_rows) + 1L]] <- data.frame(
          fit_id = fit_id, chain = chain,
          marker_id = data$markers$marker_ids[marker],
          causal = bundle$marker_truth$true_nonnull[marker],
          posterior_effect = z$bm[marker], pip = z$dm[marker])
      }
    }
    for (left in seq_len(length(chains) - 1L))
      for (right in (left + 1L):length(chains))
        marker_agreement_rows[[length(marker_agreement_rows) + 1L]] <- data.frame(
          fit_id = fit_id, chain_left = left, chain_right = right,
          pip_correlation = stats::cor(chains[[left]]$dm, chains[[right]]$dm),
          effect_correlation = stats::cor(chains[[left]]$bm, chains[[right]]$bm),
          pip_mean_absolute_difference = mean(abs(
            chains[[left]]$dm - chains[[right]]$dm)))
    alpha <- extract_annotation_coefficient_traces(
      checkpoint$result, ncol(annotations))
    alpha <- alpha[alpha$iteration > 3000L & alpha$iteration <= 9000L, ]
    annotation_range_rows[[length(annotation_range_rows) + 1L]] <- aggregate(
      value ~ parameter + annotation + stick, alpha,
      function(x) c(mean = mean(x), sd = stats::sd(x),
                    minimum = min(x), maximum = max(x)))
    annotation_range_rows[[length(annotation_range_rows)]]$fit_id <- fit_id
    blocks <- native$input$eigen_diagnostics$blocks
    if (!is.null(blocks)) {
      kept <- if ("n_kept" %in% names(blocks)) {
        blocks$n_kept
      } else {
        blocks$retained_rank
      }
      size <- if ("size" %in% names(blocks)) blocks$size else blocks$block_size
      eigen_rows[[length(eigen_rows) + 1L]] <- data.frame(
        fit_id = fit_id, blocks = nrow(blocks), minimum_kept = min(kept),
        maximum_kept = max(kept), all_modes_retained = all(kept == size))
    }
    rm(checkpoint, native, chains, alpha)
    gc(FALSE)
  }
  annotation_ranges <- do.call(rbind, annotation_range_rows)
  annotation_values <- annotation_ranges$value
  annotation_ranges$value <- NULL
  annotation_ranges <- cbind(annotation_ranges,
    mean = annotation_values[, "mean"], sd = annotation_values[, "sd"],
    minimum = annotation_values[, "minimum"],
    maximum = annotation_values[, "maximum"])
  utils::write.csv(do.call(rbind, variance_rows),
                   file.path(output_root, "variance_summary.csv"), row.names = FALSE)
  utils::write.csv(annotation_ranges,
                   file.path(output_root, "annotation_range_summary.csv"), row.names = FALSE)
  utils::write.csv(do.call(rbind, marker_agreement_rows),
                   file.path(output_root, "marker_chain_agreement.csv"), row.names = FALSE)
  utils::write.csv(do.call(rbind, representative_rows),
                   file.path(output_root, "representative_marker_summary.csv"), row.names = FALSE)
  if (length(eigen_rows)) utils::write.csv(do.call(rbind, eigen_rows),
    file.path(output_root, "eigen_retention_summary.csv"), row.names = FALSE)

  convergence <- utils::read.csv(file.path(output_root, "all_convergence.csv"))
  convergence$family <- ifelse(grepl("^alpha:", convergence$quantity), "alpha",
    ifelse(grepl("sigmaSqAlpha", convergence$quantity), "sigmaSqAlpha",
    ifelse(convergence$quantity == "prior_expected_active", "expected_active",
    ifelse(convergence$quantity == "occupancy_traced_active_count", "realized_active",
    ifelse(convergence$quantity == "effect_variance", "effect_variance",
    ifelse(convergence$quantity == "heritability", "heritability", "other"))))))
  efficiency <- do.call(rbind, lapply(
    split(convergence[convergence$family != "other", ],
          list(convergence$fit_id[convergence$family != "other"],
               convergence$family[convergence$family != "other"]), drop = TRUE),
    function(x) data.frame(
      fit_id = x$fit_id[1L], family = x$family[1L],
      minimum_bulk_ess = min(x$ess_bulk), minimum_tail_ess = min(x$ess_tail),
      median_bulk_ess_per_second = median(x$ess_bulk_per_second),
      median_tail_ess_per_second = median(x$ess_tail_per_second),
      maximum_rhat = max(x$rhat), maximum_relative_mcse = max(x$relative_mcse))))
  utils::write.csv(efficiency, file.path(output_root, "efficiency_summary.csv"),
                   row.names = FALSE)
  message("Aggregated all ten kernel-composition fits.")
}
