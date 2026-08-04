#!/usr/bin/env Rscript

# Read-only, non-qualification Study 06 BED pilot for the proper probit-stick
# intercept prior. It writes no benchmark checkpoint or capsule output.

scenario <- Sys.getenv("SBLR_STUDY06_SCENARIO", "informative_annotations")
iterations <- as.integer(Sys.getenv("SBLR_STUDY06_ITERATIONS", "300"))
sblr_root <- normalizePath(Sys.getenv("SBLR_ROOT", "."), winslash = "/")
sblrbench_root <- normalizePath(Sys.getenv("SBLRBENCH_ROOT", "../sblrbench"),
  winslash = "/")
if (!scenario %in% c("informative_annotations", "uninformative_annotations"))
  stop("Invalid SBLR_STUDY06_SCENARIO.")
if (!is.finite(iterations) || iterations < 1L)
  stop("SBLR_STUDY06_ITERATIONS must be positive.")

evidence <- c(
  file.path(sblrbench_root,
    "results/local/06_annotation_models/checkpoints/data/human_glist.rds"),
  file.path(sblrbench_root,
    "results/local/06_annotation_models/checkpoints/ld",
    "training_ld_train1400-test600-seed3101_glist.rds"),
  file.path(sblrbench_root, "studies/06_annotation_models/spec.R"),
  file.path(sblrbench_root,
    "studies/06_annotation_models/annotation-design.R")
)
before <- unname(tools::md5sum(evidence))

pkgload::load_all(sblr_root, quiet = TRUE, recompile = FALSE)
pkgload::load_all(sblrbench_root, quiet = TRUE, compile = FALSE)
old <- getwd()
on.exit(setwd(old), add = TRUE)
setwd(sblrbench_root)
study <- new.env(parent = globalenv())
sys.source(file.path(sblrbench_root,
  "studies/06_annotation_models/spec.R"), study)
sys.source(file.path(sblrbench_root,
  "studies/06_annotation_models/annotation-design.R"), study)
spec <- study$spec
data <- prepare_prediction_data(spec,
  file.path(sblrbench_root, "results/local/06_annotation_models"))
coordinates <- benchmark_annotation_seeds(spec, "benchmark",
  mode = "qualification")
hit <- coordinates$scenario == scenario & coordinates$replicate == 1L &
  coordinates$method == "st_bed_bayesrc"
coordinate <- coordinates[hit, , drop = FALSE]
if (nrow(coordinate) != 1L) stop("Expected one BED qualification coordinate.")
simulation_coordinates <- benchmark_annotation_seeds(spec, "benchmark",
  mode = "final")
simulation_coordinate <- as.list(simulation_coordinates[
  simulation_coordinates$scenario == scenario &
    simulation_coordinates$replicate == 1L, , drop = FALSE][1L, ])
annotations <- study$construct_annotation_design(data$markers$marker_ids, spec)
annotation_truth <- study$construct_annotation_truth(annotations, spec)
simulation <- study$simulate_annotation_architecture(
  simulation_coordinate, data$scaled$all, data$split$train_rows,
  annotations, annotation_truth, spec)
stats <- benchmark_summary_stats(simulation, data$ld_glist, data$split,
  spec$data)
controls <- annotation_method_controls(spec, coordinate, "benchmark",
  mode = "qualification")
controls$nit <- iterations
controls$nburn <- 0L
controls$nthin <- 1L
controls$nchains <- 4L
controls$ncores <- 4L
controls$keep_chains <- TRUE
controls$convergence <- "extended"
controls$convergence_control <- list(warn = FALSE,
  extended_groups = c("annotations", "probability"), keep_traces = TRUE,
  max_trace_gb = 2, allow_large_traces = FALSE)
controls$intercept_flat <- NULL
controls$annotation_intercept_prior <- list(
  distribution = "normal", mean = "initial_mixture", sd = 1)
methods <- resolve_benchmark_methods(spec)
method <- methods[[match("st_bed_bayesrc",
  vapply(methods, `[[`, character(1L), "id"))]]

cat("Study 06 proper-prior BED pilot\n")
cat("scenario:", scenario, "replicate: 1 iterations:", iterations,
  "chains: 4\n")
cat("chain seeds:", paste(controls$chain_seeds, collapse = ","), "\n")
result <- fit_annotation_method(method, controls, simulation, stats,
  data$ld_glist, data$split, annotations, annotation_truth)
draws <- getFromNamespace(".annotation_required_traces", "sblrbench")(
  result, coordinate,
  list(spec = spec, simulation = simulation, stats = stats,
    annotations = annotations, annotation_truth = annotation_truth,
    marker_truth = simulation$extras$marker_truth), 4L)

alpha <- draws[grepl("^(alpha|sigmaSqAlpha)", draws$quantity), ]
ranges <- aggregate(value ~ chain + quantity, alpha,
  function(x) c(min = min(x), max = max(x), mean = mean(x)))
print(ranges, row.names = FALSE)
if (any(!is.finite(alpha$value))) stop("Pilot produced non-finite annotation traces.")
native <- result$native_fit
cat("effect variance range:", range(native$vbs), "\n")
cat("genetic variance range:", range(native$vgs), "\n")
cat("residual variance range:", range(native$ves), "\n")
cat("max absolute alpha:", max(abs(alpha$value[grepl("^alpha", alpha$quantity)])),
  "\n")
cat("sigmaSqAlpha range:",
  range(alpha$value[grepl("^sigmaSqAlpha", alpha$quantity)]), "\n")
after <- unname(tools::md5sum(evidence))
cat("evidence files unchanged:", identical(before, after), "\n")
if (!identical(before, after)) stop("A read-only evidence file changed.")
