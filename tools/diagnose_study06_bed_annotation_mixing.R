#!/usr/bin/env Rscript

# Read-only analysis of the preserved Study 06 9,000 x 4 BED BayesRC traces.

parse_args <- function(x) {
  out <- list(
    sblr_root = normalizePath(".", winslash = "/", mustWork = TRUE),
    sblrbench_root = normalizePath("../sblrbench", winslash = "/",
      mustWork = TRUE),
    output = file.path(tempdir(), "study06-bed-annotation-mixing"))
  for (arg in x) {
    pair <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1L]]
    if (length(pair) != 2L || !pair[1L] %in% names(out))
      stop("Arguments must use --name=value; unknown argument: ", arg,
        call. = FALSE)
    out[[pair[1L]]] <- pair[2L]
  }
  out <- lapply(out, normalizePath, winslash = "/", mustWork = FALSE)
  out
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
if (startsWith(args$output, paste0(args$sblrbench_root, "/")))
  stop("--output must be outside the read-only sblrbench repository.")
dir.create(args$output, recursive = TRUE, showWarnings = FALSE)

checkpoint_dir <- file.path(args$sblrbench_root,
  "results/local/06_annotation_models/qualification/checkpoints")
checkpoints <- file.path(checkpoint_dir, paste0(c("informative_annotations",
  "uninformative_annotations"), "--r1--st_bed_bayesrc.rds"))
before <- unname(tools::md5sum(checkpoints))

pkgload::load_all(args$sblr_root, quiet = TRUE, recompile = FALSE)
pkgload::load_all(args$sblrbench_root, quiet = TRUE, compile = FALSE)
old_working_directory <- getwd()
on.exit(setwd(old_working_directory), add = TRUE)
setwd(args$sblrbench_root)
study_env <- new.env(parent = globalenv())
sys.source(file.path(args$sblrbench_root,
  "studies/06_annotation_models/spec.R"), envir = study_env)
sys.source(file.path(args$sblrbench_root,
  "studies/06_annotation_models/annotation-design.R"), envir = study_env)
spec <- study_env$spec
data <- prepare_prediction_data(spec, file.path(args$sblrbench_root,
  "results/local/06_annotation_models"))
annotations <- study_env$construct_annotation_design(data$markers$marker_ids,
  spec)
annotation_truth <- study_env$construct_annotation_truth(annotations, spec)
coordinates <- benchmark_annotation_seeds(spec, "benchmark",
  mode = "qualification")
simulation_coordinates <- benchmark_annotation_seeds(spec, "benchmark",
  mode = "final")

all_draws <- list()
separation_rows <- list()
for (scenario_index in seq_along(checkpoints)) {
  scenario <- c("informative_annotations",
    "uninformative_annotations")[[scenario_index]]
  coordinate <- coordinates[coordinates$scenario == scenario &
    coordinates$replicate == 1L &
    coordinates$method == "st_bed_bayesrc", , drop = FALSE]
  simulation_row <- simulation_coordinates[
    simulation_coordinates$scenario == scenario &
      simulation_coordinates$replicate == 1L, , drop = FALSE][1L, ]
  simulation <- study_env$simulate_annotation_architecture(
    as.list(simulation_row), data$scaled$all, data$split$train_rows,
    annotations, annotation_truth, spec)
  bundle <- list(spec = spec, simulation = simulation,
    annotations = annotations, marker_truth = simulation$extras$marker_truth)
  payload <- readRDS(checkpoints[[scenario_index]])
  if (!identical(payload$checkpoint_schema, "sblrbench-semantic-v2") ||
      !inherits(payload$result, "sblrbench_result"))
    stop("Unexpected qualification checkpoint contract: ", checkpoints[[scenario_index]])
  draws <- getFromNamespace(".annotation_required_traces", "sblrbench")(
    payload$result, coordinate, bundle, 4L)
  if (!identical(sort(unique(draws$iteration)), seq_len(9000L)) ||
      !identical(sort(unique(draws$chain)), seq_len(4L)))
    stop("Preserved history is not exactly 9,000 x 4.")
  all_draws[[scenario]] <- draws
  chains <- payload$result$native_fit$chains
  for (chain in seq_along(chains)) {
    component <- as.integer(chains[[chain]]$component)
    for (stick in 0:2) {
      eligible <- stick == 0L | component > (stick - 1L)
      success <- component > stick
      separation_rows[[length(separation_rows) + 1L]] <- data.frame(
        scenario = scenario, chain = chain, stick = stick,
        eligible = sum(eligible), successes = sum(eligible & success),
        failures = sum(eligible & !success),
        complete_separation = sum(eligible & success) == 0L ||
          sum(eligible & !success) == 0L,
        stringsAsFactors = FALSE)
    }
  }
}
draws <- do.call(rbind, all_draws)
separation <- do.call(rbind, separation_rows)

chain_groups <- split(draws, interaction(draws$scenario, draws$quantity,
  draws$chain, drop = TRUE, lex.order = TRUE))
chain_summary <- do.call(rbind, lapply(chain_groups, function(x) data.frame(
  scenario = x$scenario[[1L]], quantity = x$quantity[[1L]],
  chain = x$chain[[1L]], mean = mean(x$value), sd = stats::sd(x$value),
  minimum = min(x$value), maximum = max(x$value),
  acf_1 = stats::acf(x$value, lag.max = 100L, plot = FALSE)$acf[2L],
  acf_10 = stats::acf(x$value, lag.max = 100L, plot = FALSE)$acf[11L],
  acf_50 = stats::acf(x$value, lag.max = 100L, plot = FALSE)$acf[51L],
  acf_100 = stats::acf(x$value, lag.max = 100L, plot = FALSE)$acf[101L],
  stringsAsFactors = FALSE)))

quantity_groups <- split(draws, interaction(draws$scenario, draws$quantity,
  drop = TRUE, lex.order = TRUE))
diagnostics <- do.call(rbind, lapply(quantity_groups, function(x) {
  z <- benchmark_scalar_diagnostics(x, spec$qualification$thresholds)
  cbind(scenario = x$scenario[[1L]], quantity = x$quantity[[1L]], z,
    stringsAsFactors = FALSE)
}))

correlations <- do.call(rbind, lapply(split(draws,
    interaction(draws$scenario, draws$chain, drop = TRUE)), function(x) {
  wide <- reshape(x[c("iteration", "quantity", "value")],
    idvar = "iteration", timevar = "quantity", direction = "wide")
  values <- as.matrix(wide[-1L])
  names <- sub("^value\\.", "", colnames(values))
  correlation <- stats::cor(values)
  index <- which(upper.tri(correlation), arr.ind = TRUE)
  data.frame(scenario = x$scenario[[1L]], chain = x$chain[[1L]],
    quantity_1 = names[index[, 1L]], quantity_2 = names[index[, 2L]],
    correlation = correlation[index], stringsAsFactors = FALSE)
}))

utils::write.csv(chain_summary, file.path(args$output, "chain_summary.csv"),
  row.names = FALSE)
utils::write.csv(diagnostics, file.path(args$output, "diagnostics.csv"),
  row.names = FALSE)
utils::write.csv(correlations,
  file.path(args$output, "cross_parameter_correlations.csv"), row.names = FALSE)
utils::write.csv(separation,
  file.path(args$output, "final_stick_separation.csv"), row.names = FALSE)

grDevices::pdf(file.path(args$output, "registered_quantity_traces.pdf"),
  width = 11, height = 8.5, onefile = TRUE)
for (scenario in names(all_draws)) {
  x <- all_draws[[scenario]]
  for (quantity in unique(x$quantity)) {
    z <- x[x$quantity == quantity, ]
    plot(range(z$iteration), range(z$value), type = "n",
      xlab = "Iteration", ylab = quantity, main = scenario)
    for (chain in seq_len(4L)) {
      value <- z[z$chain == chain, ]
      lines(value$iteration, value$value, col = chain)
    }
  }
}
grDevices::dev.off()

after <- unname(tools::md5sum(checkpoints))
if (!identical(before, after)) stop("A preserved BED checkpoint changed.")
cat("BED mixing outputs:", args$output, "\n")
cat("Worst diagnostics by scenario:\n")
print(do.call(rbind, lapply(split(diagnostics, diagnostics$scenario),
  function(x) x[order(-x$rhat, x$ess_bulk), ][1:8, ])), row.names = FALSE)
cat("Strongest absolute cross-parameter correlations:\n")
print(head(correlations[order(-abs(correlations$correlation)), ], 20L),
  row.names = FALSE)
cat("Final-state stick separation:\n")
print(separation, row.names = FALSE)
cat("Verified: preserved BED checkpoints are byte-identical.\n")
