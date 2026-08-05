#!/usr/bin/env Rscript

# Read-only Study 06 SBayesRC diagnostic replay. This tool never writes below
# sblrbench; diagnostic state is written to --output.

parse_args <- function(x) {
  out <- list(
    sblr_root = normalizePath(".", winslash = "/", mustWork = TRUE),
    sblrbench_root = normalizePath("../sblrbench", winslash = "/",
      mustWork = TRUE),
    scenario = "informative_annotations",
    chain = 2L,
    route = "csr",
    diagnostic_iterations = 50L,
    output = file.path(tempdir(), "study06-sbayesrc-failure.tsv"),
    replay = TRUE
  )
  for (arg in x) {
    pair <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1L]]
    if (length(pair) != 2L || !pair[1L] %in% names(out))
      stop("Arguments must use --name=value; unknown argument: ", arg,
        call. = FALSE)
    out[[pair[1L]]] <- pair[2L]
  }
  out$chain <- as.integer(out$chain)
  out$diagnostic_iterations <- as.integer(out$diagnostic_iterations)
  out$replay <- !tolower(as.character(out$replay)) %in% c("false", "0", "no")
  out$sblr_root <- normalizePath(out$sblr_root, winslash = "/", mustWork = TRUE)
  out$sblrbench_root <- normalizePath(out$sblrbench_root, winslash = "/",
    mustWork = TRUE)
  out$output <- normalizePath(out$output, winslash = "/", mustWork = FALSE)
  out
}

hash_files <- function(paths) {
  paths <- normalizePath(paths, winslash = "/", mustWork = TRUE)
  stats::setNames(unname(tools::md5sum(paths)), paths)
}

read_failure_state <- function(path) {
  fields <- strsplit(readLines(path, warn = FALSE), "\t", fixed = TRUE)
  width <- max(lengths(fields))
  x <- as.data.frame(do.call(rbind, lapply(fields, function(value)
    c(value, rep("", width - length(value))))), stringsAsFactors = FALSE)
  meta <- x[x[[1L]] == "meta", 2:3, drop = FALSE]
  marker <- x[x[[1L]] == "marker" & x[[2L]] != "index", , drop = FALSE]
  names(marker) <- c("record", "index", "effect", "residual",
    "rebuilt_residual", "score", "diagonal", "component")
  numeric_columns <- setdiff(names(marker), "record")
  marker[numeric_columns] <- lapply(marker[numeric_columns], as.numeric)
  list(meta = stats::setNames(meta[[2L]], meta[[1L]]), marker = marker,
    raw = x)
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
if (!args$scenario %in% c("informative_annotations",
    "uninformative_annotations"))
  stop("--scenario must identify one of the two Study 06 scenarios.",
    call. = FALSE)
if (is.na(args$chain) || args$chain < 1L || args$chain > 4L)
  stop("--chain must be 1, 2, 3, or 4.", call. = FALSE)
if (!args$route %in% c("csr", "retained"))
  stop("--route must be csr or retained.", call. = FALSE)
if (startsWith(args$output, paste0(args$sblrbench_root, "/")))
  stop("--output must be outside the read-only sblrbench repository.",
    call. = FALSE)
dir.create(dirname(args$output), recursive = TRUE, showWarnings = FALSE)

checkpoint_root <- file.path(args$sblrbench_root,
  "results/local/06_annotation_models/checkpoints")
evidence_paths <- c(
  file.path(checkpoint_root, "data/human_glist.rds"),
  file.path(checkpoint_root, "ld",
    "training_ld_train1400-test600-seed3101_glist.rds"),
  file.path(args$sblrbench_root, "studies/06_annotation_models/spec.R"),
  file.path(args$sblrbench_root,
    "studies/06_annotation_models/annotation-design.R"))
before <- hash_files(evidence_paths)

pkgload::load_all(args$sblr_root, quiet = TRUE, recompile = FALSE)
pkgload::load_all(args$sblrbench_root, quiet = TRUE, compile = FALSE)
old_working_directory <- getwd()
on.exit(setwd(old_working_directory), add = TRUE)
# The preserved Glist intentionally retains capsule-relative BED paths.
setwd(args$sblrbench_root)
study_env <- new.env(parent = globalenv())
sys.source(file.path(args$sblrbench_root,
  "studies/06_annotation_models/spec.R"), envir = study_env)
sys.source(file.path(args$sblrbench_root,
  "studies/06_annotation_models/annotation-design.R"), envir = study_env)
spec <- study_env$spec
validate_benchmark_spec(spec)
if (!identical(as.character(utils::packageVersion("sblr")), "0.2.0"))
  stop("Study 06 requires sblr 0.2.0.", call. = FALSE)

data_root <- file.path(args$sblrbench_root,
  "results/local/06_annotation_models")
data <- prepare_prediction_data(spec, data_root)
coordinates <- benchmark_annotation_seeds(spec, "benchmark",
  mode = "qualification")
hit <- coordinates$scenario == args$scenario &
  coordinates$replicate == 1L &
  coordinates$method == "st_csr_sbayesrc"
if (sum(hit) != 1L) stop("Expected one Study 06 CSR qualification coordinate.")
coordinate <- coordinates[hit, , drop = FALSE]
simulation_coordinates <- benchmark_annotation_seeds(spec, "benchmark",
  mode = "final")
simulation_hit <- simulation_coordinates$scenario == args$scenario &
  simulation_coordinates$replicate == 1L
simulation_coordinate <- as.list(simulation_coordinates[
  which(simulation_hit)[1L], , drop = FALSE])
annotations <- study_env$construct_annotation_design(
  data$markers$marker_ids, spec)
annotation_truth <- study_env$construct_annotation_truth(annotations, spec)
simulation <- study_env$simulate_annotation_architecture(
  simulation_coordinate, data$scaled$all, data$split$train_rows,
  annotations, annotation_truth, spec)
summary_stats <- benchmark_summary_stats(simulation, data$ld_glist,
  data$split, spec$data)
bundle <- list(coordinate = simulation_coordinate, simulation = simulation,
  stats = summary_stats, annotations = annotations,
  annotation_truth = annotation_truth)
controls <- annotation_method_controls(spec, coordinate, "benchmark",
  mode = "qualification")
authoritative_seeds <- controls$chain_seeds
controls$nchains <- 1L
controls$ncores <- 1L
controls$chain_seeds <- authoritative_seeds[args$chain]
method <- resolve_benchmark_methods(spec)
method <- method[[match("st_csr_sbayesrc",
  vapply(method, `[[`, character(1L), "id"))]]

stopifnot(
  length(data$sample_ids) == 2000L,
  length(data$split$train_ids) == 1400L,
  length(data$markers$marker_ids) == 37991L,
  identical(rownames(bundle$annotations), data$markers$marker_ids),
  identical(bundle$simulation$data$marker_ids, data$markers$marker_ids),
  identical(as.integer(bundle$stats$n), 1400L),
  identical(as.integer(controls$nit), 9000L),
  identical(as.integer(controls$nburn), 0L),
  identical(as.integer(controls$chain_seeds),
    as.integer(authoritative_seeds[args$chain])))

cat("Study 06 scalar SBayesRC operator diagnostic\n")
cat("sblr HEAD:", system2("git", c("-C", shQuote(args$sblr_root),
  "rev-parse", "HEAD"), stdout = TRUE), "\n")
cat("sblrbench HEAD:", system2("git", c("-C", shQuote(args$sblrbench_root),
  "rev-parse", "HEAD"), stdout = TRUE), "\n")
cat("scenario:", args$scenario, "replicate: 1 logical chain:", args$chain,
  "seed:", controls$chain_seeds, "\n")
cat("samples: 1400 markers: 37991 iterations:",
  if (identical(args$route, "retained")) args$diagnostic_iterations else 9000L,
  "\n")
cat("LD:", normalizePath(evidence_paths[2L], winslash = "/"), "\n")
cat("diagnostic output:", args$output, "\n")

if (identical(args$route, "retained")) {
  # This is a post-identification diagnostic, not an authoritative replay of
  # the historical flat-prior qualification coordinate.
  controls$intercept_flat <- NULL
  controls$annotation_intercept_prior <- list(
    distribution = "normal", mean = "initial_mixture", sd = 1)
  controls$nit <- args$diagnostic_iterations
  controls$nburn <- 0L
  controls$nthin <- 1L
  controls$nchains <- 1L
  controls$ncores <- 1L
  controls$keep_chains <- FALSE
  controls$convergence <- "none"
  controls$convergence_control <- NULL
  controls$gamma <- controls$mixture_var
  controls$mixture_var <- NULL
  controls$alpha_init <- bundle$annotation_truth$uninformative_annotations
  retained <- do.call(stblr_block_eigen, c(list(
    stats = bundle$stats, Glist = data$ld_glist,
    block_start = seq.int(1L, nrow(bundle$annotations), by = 1024L),
    method = "sbayesrc", annotation = bundle$annotations,
    representation = "low_rank", eigen_prop = 0.995,
    low_rank_residual_rebuild_every = 1L), controls))
  diagnostic <- retained$diagnostics$native$low_rank_residual
  cat("retained result: completed\n")
  cat("minimum residual variance:", min(retained$ves), "\n")
  cat("maximum residual variance:", max(retained$ves), "\n")
  cat("minimum genetic variance:", min(retained$vgs), "\n")
  cat("residual rebuild count:",
    diagnostic$low_rank_residual_rebuild_count, "\n")
  cat("maximum reduced-residual drift:",
    diagnostic$low_rank_residual_max_abs_drift, "\n")
  if (!all(is.finite(retained$ves)) || any(retained$ves <= 0))
    stop("Retained diagnostic produced invalid residual variance.")
  after <- hash_files(evidence_paths)
  if (!identical(before, after)) stop("Read-only Study 06 evidence changed.")
  cat("evidence files unchanged: TRUE\n")
  quit(save = "no", status = 0L)
}

if (args$replay) {
  if (file.exists(args$output)) unlink(args$output)
  old <- Sys.getenv("SBLR_SBAYESRC_FAILURE_STATE_PATH", unset = NA_character_)
  on.exit(if (is.na(old)) Sys.unsetenv("SBLR_SBAYESRC_FAILURE_STATE_PATH") else
    Sys.setenv(SBLR_SBAYESRC_FAILURE_STATE_PATH = old), add = TRUE)
  Sys.setenv(SBLR_SBAYESRC_FAILURE_STATE_PATH = args$output)
  error <- tryCatch({
    fit_annotation_method(method, controls, bundle$simulation, bundle$stats,
      data$ld_glist, data$split, bundle$annotations,
      bundle$annotation_truth)
    "NO ERROR"
  }, error = conditionMessage)
  cat("replay result:", error, "\n")
}

if (!file.exists(args$output))
  stop("No failure-state file exists at --output.", call. = FALSE)
state <- read_failure_state(args$output)
print(as.data.frame(as.list(state$meta)), row.names = FALSE)

# Memory-safe exact-reference counterfactual: form Xb, never X'X.
X <- data$scaled$all[data$split$train_rows, , drop = FALSE]
b <- state$marker$effect
wy <- state$marker$score
if (ncol(X) != length(b)) stop("Captured effect vector is misaligned.")
xb <- as.numeric(X %*% b)
q_exact <- sum(xb * xb)
bwy <- sum(b * wy)
yy <- as.numeric(state$meta[["yy"]])
prior <- as.numeric(state$meta[["prior_contribution"]])
q_sparse <- as.numeric(state$meta[["independent_quadratic"]])
scale_exact <- yy - 2 * bwy + q_exact + prior
scale_sparse <- yy - 2 * bwy + q_sparse + prior
cat(sprintf(paste0("exact counterfactual: q_exact=%.17g q_sparse=%.17g ",
  "delta_q=%.17g scale_exact=%.17g scale_sparse=%.17g\n"),
  q_exact, q_sparse, q_exact - q_sparse, scale_exact, scale_sparse))

# Exact corrected scores are computed without materializing genome-wide LD.
corrected_exact <- wy - as.numeric(crossprod(X, xb))
corrected_sparse <- state$marker$rebuilt_residual
score_error <- corrected_sparse - corrected_exact
cat(sprintf(paste0("corrected-score error: max_abs=%.17g rms=%.17g ",
  "relative_l2=%.17g\n"), max(abs(score_error)),
  sqrt(mean(score_error^2)),
  sqrt(sum(score_error^2)) / max(1, sqrt(sum(corrected_exact^2)))))

# Omitted-LD decomposition is restricted to active markers. This state has a
# few thousand active markers; no genome-wide dense matrix is constructed.
active <- which(b != 0)
if (length(active) <= 5000L) {
  X_active <- X[, active, drop = FALSE]
  cross_active <- crossprod(X_active)
  diagonal <- diag(cross_active)
  upper <- which(upper.tri(cross_active), arr.ind = TRUE)
  marker_i <- active[upper[, 1L]]
  marker_j <- active[upper[, 2L]]
  cross_value <- cross_active[upper]
  correlation <- cross_value /
    sqrt(diagonal[upper[, 1L]] * diagonal[upper[, 2L]])
  contribution <- 2 * b[marker_i] * b[marker_j] * cross_value
  distance <- marker_j - marker_i
  within_window <- distance <= 1000L
  above_threshold <- correlation^2 >= 0.001
  category <- ifelse(!within_window, "omitted_outside_1000_markers",
    ifelse(!above_threshold, "omitted_below_r2_0.001", "retained_sparse"))
  bins <- cut(abs(correlation), c(0, 0.01, 0.02, sqrt(0.001), 0.05,
    0.1, 0.2, 0.5, 1, Inf), include.lowest = TRUE, right = FALSE)
  contribution_sign <- ifelse(contribution >= 0, "positive", "negative")
  grouping <- interaction(category, bins, contribution_sign, drop = TRUE,
    lex.order = TRUE)
  groups <- split(seq_along(contribution), grouping)
  decomposition <- do.call(rbind, lapply(groups, function(index) data.frame(
    category = category[index[[1L]]],
    absolute_ld_bin = as.character(bins[index[[1L]]]),
    contribution_sign = contribution_sign[index[[1L]]],
    pair_count = length(index), quadratic = sum(contribution[index]),
    stringsAsFactors = FALSE)))
  decomposition$scenario <- args$scenario
  utils::write.csv(decomposition,
    sub("\\.tsv$", "-omitted-ld-decomposition.csv", args$output),
    row.names = FALSE)
  omitted <- category != "retained_sparse"
  cat(sprintf(paste0("active-pair decomposition: active=%d diagonal=%.17g ",
    "retained_pairs_q=%.17g omitted_pairs_q=%.17g ",
    "omitted_positive=%.17g omitted_negative=%.17g\n"),
    length(active), sum(diagonal * b[active]^2),
    sum(contribution[!omitted]), sum(contribution[omitted]),
    sum(contribution[omitted & contribution >= 0]),
    sum(contribution[omitted & contribution < 0])))
  high_ld_opposite <- abs(correlation) >= 0.8 &
    sign(b[marker_i]) != sign(b[marker_j])
  cat("high-LD opposite-signed active pairs:", sum(high_ld_opposite), "\n")
  rm(X_active, cross_active)
}

after <- hash_files(evidence_paths)
if (!identical(before, after))
  stop("A read-only Study 06 evidence file changed during diagnostics.")
cat("Verified: selected sblrbench evidence files are byte-identical.\n")
