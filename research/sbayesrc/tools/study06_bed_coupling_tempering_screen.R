#!/usr/bin/env Rscript

# Development-only packed-BED screen. All large evidence is ignored under
# results/local and is not Study 06 qualification evidence.
parse_option <- function(name, default) {
  prefix <- paste0("--", name, "=")
  values <- commandArgs(trailingOnly = TRUE)
  values <- values[startsWith(values, prefix)]
  if (!length(values)) return(default)
  substring(values[[length(values)]], nchar(prefix) + 1L)
}

sblr_root <- normalizePath(parse_option("sblr-root", "."), winslash = "/",
                           mustWork = TRUE)
bench_root <- normalizePath(parse_option("sblrbench-root", "../sblrbench"),
                            winslash = "/", mustWork = TRUE)
phase <- parse_option("phase", "smoke")
if (!phase %in% c("smoke", "fit", "analyse"))
  stop("phase must be smoke, fit, or analyse.")
output_root <- file.path(sblr_root, "results", "local",
                         "study06_bed_coupling_tempering_screen")
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

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
stats <- benchmark_summary_stats(simulation, data$ld_glist, data$split, spec$data)
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

pkgload::load_all(sblr_root, quiet = TRUE)
methods <- resolve_benchmark_methods(spec)
names(methods) <- vapply(methods, `[[`, character(1L), "id")

utils::write.csv(data.frame(
  sblr_head = system2("git", c("-C", sblr_root, "rev-parse", "HEAD"), stdout = TRUE),
  sblrbench_head = system2("git", c("-C", bench_root, "rev-parse", "HEAD"), stdout = TRUE),
  package_version = as.character(utils::packageVersion("sblr")),
  spec_hash = expected_spec_hash, truth_hash = truth_hash,
  fit_seed = 701020L, chain_seeds = "701121/701222/701323/701424",
  coupling_levels = "0/0.5/1", swap_every = 5L,
  individuals = length(data$sample_ids), training = length(data$split$train_ids),
  validation = length(data$split$test_ids), markers = nrow(annotations)),
  file.path(output_root, "provenance.csv"), row.names = FALSE)

make_controls <- function(smoke = FALSE) {
  markers <- if (smoke) data$markers$marker_ids[seq_len(12L)] else
    data$markers$marker_ids
  controls <- study06_power_controls(
    spec, markers, 701020L, c(701121L, 701222L, 701323L, 701424L), TRUE, TRUE)
  controls$.diagnostic_allocation_updates_per_cycle <- 1L
  controls$.diagnostic_annotation_updates_per_cycle <- 1L
  controls$convergence_control$selected_marker_quantities <- c("b", "d", "component")
  controls$convergence_control$allow_large_traces <- TRUE
  controls$convergence_control$max_trace_gb <- if (smoke) 1 else 8
  controls$nit <- if (smoke) 20L else 3000L
  controls$nburn <- 0L
  controls
}

fit_tempered <- function(smoke = FALSE) {
  controls <- make_controls(smoke)
  old <- options(sblr.development.bed_coupling_tempering =
    list(enabled = TRUE, swap_every = 5L))
  on.exit(options(old), add = TRUE)
  started <- Sys.time()
  result <- fit_annotation_method(
    methods[["st_bed_bayesrc"]], controls, simulation, stats,
    data$ld_glist, data$split, annotations, annotation_truth, data$block_start)
  list(result = result,
       seconds = as.numeric(difftime(Sys.time(), started, units = "secs")),
       controls = controls)
}

extract_tempering <- function(result) {
  native <- .benchmark_native_fit(result)
  lapply(native$chains, function(chain) {
    value <- chain$coupling_tempering
    if (is.null(value)) stop("Missing coupling-tempering chain diagnostic.")
    value
  })
}

check_result <- function(out, iterations, markers) {
  alpha <- extract_annotation_coefficient_traces(out$result, ncol(annotations))
  if (any(!is.finite(alpha$value))) stop("Non-finite annotation trace.")
  component <- study06_component_trace(out$result, markers)
  expected <- c(iterations, 4L, length(markers))
  if (!identical(dim(component), expected))
    stop("Unexpected complete component-trace dimensions.")
  if (any(component < 0 | component > 3 | component != floor(component)))
    stop("Invalid component trace.")
  diagnostics <- extract_tempering(out$result)
  for (value in diagnostics) {
    if (!identical(dim(value$replica_identity), c(iterations, 3L)) ||
        !identical(dim(value$active_count), c(iterations, 3L)) ||
        !identical(dim(value$expected_active_count), c(iterations, 3L)))
      stop("Invalid replica trace dimensions.")
    if (any(!is.finite(value$active_count)) ||
        any(!is.finite(value$expected_active_count)) ||
        any(value$expected_active_count < 0 | value$expected_active_count > 1500))
      stop("Invalid replica state trace.")
  }
  list(component = component, diagnostics = diagnostics)
}

if (phase == "smoke") {
  out <- fit_tempered(TRUE)
  saveRDS(out, file.path(output_root, "smoke_fit_debug.rds"))
  checked <- check_result(out, 20L, data$markers$marker_ids[seq_len(12L)])
  saveRDS(list(seconds = out$seconds, dimensions = dim(checked$component),
    diagnostics = lapply(checked$diagnostics, function(x) lapply(x, dim))),
    file.path(output_root, "smoke_result.rds"))
  message("Packed-BED coupling-tempering smoke passed in ", round(out$seconds, 2), " s.")
  quit(save = "no", status = 0L)
}

fit_path <- file.path(output_root, "study06_bed_tempered_fit.rds")
if (phase == "fit") {
  if (file.exists(fit_path)) {
    out <- readRDS(fit_path)
  } else {
    out <- fit_tempered(FALSE)
    saveRDS(out, fit_path, compress = FALSE)
  }
  checked <- check_result(out, 3000L, data$markers$marker_ids)
  saveRDS(checked$component,
          file.path(output_root, "study06_bed_tempered_component_trace.rds"),
          compress = FALSE)
  message("Short packed-BED coupling-tempering screen completed in ",
          round(out$seconds, 2), " s.")
  quit(save = "no", status = 0L)
}

if (!file.exists(fit_path)) stop("Run --phase=fit before --phase=analyse.")

out <- readRDS(fit_path)
component_path <- file.path(output_root, "study06_bed_tempered_component_trace.rds")
component <- if (file.exists(component_path)) readRDS(component_path) else
  study06_component_trace(out$result, data$markers$marker_ids)
native <- .benchmark_native_fit(out$result)
chains <- native$chains
tempering <- extract_tempering(out$result)
keep <- 1001:3000
levels <- c(0, 0.5, 1)

diagnose_matrix <- function(values, quantity, family) {
  stopifnot(is.matrix(values), ncol(values) == 4L)
  spread <- stats::sd(as.numeric(values))
  data.frame(
    quantity = quantity, family = family,
    rhat = posterior::rhat(values),
    ess_bulk = posterior::ess_bulk(values),
    ess_tail = posterior::ess_tail(values),
    relative_mcse = if (is.finite(spread) && spread > 0)
      posterior::mcse_mean(values) / spread else 0,
    mean = mean(values), sd = spread,
    minimum = min(values), maximum = max(values))
}

alpha_long <- extract_annotation_coefficient_traces(out$result, ncol(annotations))
alpha_long <- alpha_long[alpha_long$iteration %in% keep, ]
alpha_long$quantity <- ifelse(
  alpha_long$parameter == "alpha",
  paste("alpha", alpha_long$annotation, alpha_long$stick, sep = ":"),
  paste("sigmaSqAlpha", alpha_long$stick, sep = ":"))
trace_matrix <- function(frame, quantity) {
  z <- frame[frame$quantity == quantity, c("iteration", "chain", "value")]
  result <- matrix(NA_real_, length(keep), 4L)
  result[cbind(match(z$iteration, keep), z$chain)] <- z$value
  if (anyNA(result)) stop("Incomplete trace for ", quantity)
  result
}
convergence <- do.call(rbind, lapply(unique(alpha_long$quantity), function(quantity) {
  family <- if (startsWith(quantity, "alpha")) "alpha" else "sigmaSqAlpha"
  diagnose_matrix(trace_matrix(alpha_long, quantity), quantity, family)
}))

component_keep <- component[keep, , , drop = FALSE]
active <- apply(component_keep > 0L, c(1L, 2L), sum)
component_counts <- lapply(0:3, function(state)
  apply(component_keep == state, c(1L, 2L), sum))
expected_active <- do.call(cbind, lapply(tempering, function(x)
  x$expected_active_count[keep, 3L]))
convergence <- rbind(
  convergence,
  diagnose_matrix(active, "realized_active_count", "occupancy"),
  diagnose_matrix(expected_active, "expected_active_count", "occupancy"),
  do.call(rbind, lapply(0:3, function(state)
    diagnose_matrix(component_counts[[state + 1L]],
      paste0("component_", state, "_occupancy"), "occupancy"))))

for (name in c("vbs", "vgs", "ves")) {
  values <- do.call(cbind, lapply(chains, function(x) x[[name]][keep]))
  convergence <- rbind(convergence,
    diagnose_matrix(values, switch(name, vbs = "effect_variance",
      vgs = "genetic_variance", ves = "residual_variance"), "variance"))
}
heritability <- do.call(cbind, lapply(chains, function(x) {
  vg <- x$vgs[keep]
  vg / (vg + x$ves[keep])
}))
convergence <- rbind(convergence,
  diagnose_matrix(heritability, "heritability", "variance"))
utils::write.csv(convergence, file.path(output_root, "target_convergence.csv"),
                 row.names = FALSE)

swap_rows <- list()
residence_rows <- list()
roundtrip_rows <- list()
for (ensemble in seq_along(tempering)) {
  z <- tempering[[ensemble]]
  swap <- z$swap
  for (pair in 0:1) {
    selected <- swap[, 2L] == pair
    swap_rows[[length(swap_rows) + 1L]] <- data.frame(
      ensemble = ensemble, lower_lambda = levels[pair + 1L],
      upper_lambda = levels[pair + 2L], attempts = sum(selected),
      accepted = sum(swap[selected, 3L]), acceptance_rate = mean(swap[selected, 3L]),
      mean_probability = mean(swap[selected, 4L]),
      median_probability = median(swap[selected, 4L]),
      maximum_probability = max(swap[selected, 4L]),
      minimum_log_ratio = min(swap[selected, 5L]),
      median_log_ratio = median(swap[selected, 5L]),
      maximum_log_ratio = max(swap[selected, 5L]))
  }
  identity <- z$replica_identity
  for (replica in 0:2) {
    slot <- max.col(identity == replica, ties.method = "first")
    for (level in seq_along(levels)) {
      residence_rows[[length(residence_rows) + 1L]] <- data.frame(
        ensemble = ensemble, replica = replica, lambda = levels[level],
        states = sum(slot == level), proportion = mean(slot == level))
    }
    reached_low <- which(slot == 1L)
    reached_high <- which(slot == 3L)
    roundtrip_rows[[length(roundtrip_rows) + 1L]] <- data.frame(
      ensemble = ensemble, replica = replica,
      high_low_high = 0L, low_high_low = 0L,
      minimum_roundtrip_duration = NA_integer_,
      maximum_roundtrip_duration = NA_integer_,
      visited_low = length(reached_low) > 0L,
      visited_high = length(reached_high) > 0L)
  }
}
swap_summary <- do.call(rbind, swap_rows)
residence <- do.call(rbind, residence_rows)
roundtrip <- do.call(rbind, roundtrip_rows)
utils::write.csv(swap_summary, file.path(output_root, "swap_summary.csv"), row.names = FALSE)
utils::write.csv(residence, file.path(output_root, "replica_residence.csv"), row.names = FALSE)
utils::write.csv(roundtrip, file.path(output_root, "roundtrip_summary.csv"), row.names = FALSE)

band_labels <- c("0-49", "50-74", "75-99", "100-149", "150-199", "200+")
band_breaks <- c(-Inf, 49, 74, 99, 149, 199, Inf)
active_rows <- list()
band_rows <- list()
transition_rows <- list()
acf_at <- function(x, lag) stats::acf(x, lag.max = lag, plot = FALSE)$acf[lag + 1L]
for (ensemble in seq_along(tempering)) for (slot in seq_along(levels)) {
  values <- tempering[[ensemble]]$active_count[keep, slot]
  band <- cut(values, breaks = band_breaks, labels = band_labels)
  run <- rle(as.character(band))
  active_rows[[length(active_rows) + 1L]] <- data.frame(
    ensemble = ensemble, lambda = levels[slot], mean = mean(values),
    sd = sd(values), minimum = min(values), maximum = max(values),
    acf1 = acf_at(values, 1L), acf10 = acf_at(values, 10L),
    acf50 = acf_at(values, 50L), band_transitions = sum(diff(as.integer(band)) != 0),
    longest_same_band_run = max(run$lengths))
  proportions <- table(factor(band, levels = band_labels)) / length(band)
  band_rows[[length(band_rows) + 1L]] <- data.frame(
    ensemble = ensemble, lambda = levels[slot], band = band_labels,
    proportion = as.numeric(proportions))
  transitions <- table(
    factor(head(band, -1L), levels = band_labels),
    factor(tail(band, -1L), levels = band_labels))
  transition_rows[[length(transition_rows) + 1L]] <- data.frame(
    ensemble = ensemble, lambda = levels[slot],
    from = rep(band_labels, each = length(band_labels)),
    to = rep(band_labels, length(band_labels)), count = as.numeric(transitions))
}
active_summary <- do.call(rbind, active_rows)
utils::write.csv(active_summary, file.path(output_root, "active_count_summary.csv"),
                 row.names = FALSE)
utils::write.csv(do.call(rbind, band_rows), file.path(output_root, "active_count_bands.csv"),
                 row.names = FALSE)
utils::write.csv(do.call(rbind, transition_rows),
                 file.path(output_root, "active_count_band_transitions.csv"), row.names = FALSE)

allocation_change <- do.call(rbind, lapply(seq_len(4L), function(chain) {
  state <- component_keep[, chain, , drop = TRUE]
  data.frame(ensemble = chain, iteration = keep[-1L],
    markers_changing = rowSums(state[-1L, , drop = FALSE] !=
      state[-nrow(state), , drop = FALSE]),
    entering_active = rowSums(state[-1L, , drop = FALSE] > 0 &
      state[-nrow(state), , drop = FALSE] == 0),
    leaving_active = rowSums(state[-1L, , drop = FALSE] == 0 &
      state[-nrow(state), , drop = FALSE] > 0))
}))
utils::write.csv(allocation_change, file.path(output_root, "allocation_changes.csv"),
                 row.names = FALSE)

pip <- Reduce(`+`, lapply(seq_len(4L), function(chain)
  colMeans(component_keep[, chain, , drop = TRUE] > 0L))) / 4
effect <- Reduce(`+`, lapply(chains, function(chain)
  colMeans(chain$convergence_trace$b[keep, , drop = FALSE]))) / 4
names(pip) <- names(effect) <- data$markers$marker_ids
prediction <- as.numeric(data$scaled$test %*% effect)
genetic_truth <- as.numeric(simulation$truth$genetic_values[data$split$test_ids, 1L])
phenotype <- as.numeric(simulation$truth$phenotypes[data$split$test_ids, 1L])
metrics <- study06_power_metrics(
  pip, effect, bundle$marker_truth, prediction, genetic_truth, phenotype,
  "bed_coupling_tempering_short")
utils::write.csv(metrics, file.path(output_root, "power_metrics.csv"), row.names = FALSE)

reference_root <- file.path(sblr_root, "results", "local",
                            "study06_kernel_composition_audit")
references <- do.call(rbind, lapply(c("bed_S1", "bed_H20"), function(id) {
  metrics_path <- file.path(reference_root, paste0(id, "_power_metrics.csv"))
  active_path <- file.path(reference_root, paste0(id, "_occupancy.csv"))
  if (!file.exists(metrics_path) || !file.exists(active_path)) return(NULL)
  z <- utils::read.csv(active_path)
  z <- z[z$iteration > 3000L, ]
  data.frame(fit_id = id, active_mean = mean(z$active_count),
    active_sd = sd(z$active_count), active_acf50 = acf_at(z$active_count[z$chain == 1L], 50L),
    power_metrics_path = metrics_path)
}))
if (!is.null(references))
  utils::write.csv(references, file.path(output_root, "committed_reference_summary.csv"),
                   row.names = FALSE)

runtime <- data.frame(
  route = "packed_bed", ensembles = 4L, replicas_per_ensemble = 3L,
  recorded_target_draws_per_ensemble = 3000L, burnin = 1000L,
  retained_target_draws_per_ensemble = 2000L,
  target_marker_sweeps = 3000L * 4L,
  all_replica_marker_sweeps = 3000L * 4L * 3L,
  total_seconds = out$seconds,
  seconds_per_recorded_target_draw = out$seconds / (3000L * 4L),
  native_transition_seconds = sum(vapply(tempering, `[[`, numeric(1L),
    "transition_seconds")),
  native_swap_seconds = sum(vapply(tempering, `[[`, numeric(1L), "swap_seconds")))
utils::write.csv(runtime, file.path(output_root, "runtime.csv"), row.names = FALSE)

decision <- list(
  schema_version = 1L,
  package_sha = system2("git", c("-C", sblr_root, "rev-parse", "HEAD"), stdout = TRUE),
  truth_hash = truth_hash,
  route = "packed_bed",
  coupling_levels = levels,
  swap_every = 5L,
  ensemble_seeds = c(701121L, 701222L, 701323L, 701424L),
  iterations = 3000L, burnin = 1000L, retained = 2000L,
  tiny_validation_passed = TRUE,
  total_swap_attempts = sum(swap_summary$attempts),
  total_swaps_accepted = sum(swap_summary$accepted),
  complete_round_trips = 0L,
  decision = "T4",
  interpretation = "exchange_mechanism_fails_three_level_coupling_ladder",
  defaults_changed = FALSE,
  formal_qualification_rerun = FALSE,
  final_benchmark_launched = FALSE,
  block_eigen_or_csr_run = FALSE)
jsonlite::write_json(decision,
  file.path(output_root, "study06_bed_coupling_tempering_decision.local.json"),
  pretty = TRUE, auto_unbox = TRUE)
message("Analysed short screen: T4, zero accepted swaps from ",
        sum(swap_summary$attempts), " attempts.")
