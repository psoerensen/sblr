# Development-only Phase 5D diagnostic for the frozen Study 07 case.
#
# This script reads ../sblrbench but writes only below sblr/results/local.
# It does not alter MCEM scientific settings. Usage examples:
#   Rscript research/sbayesrc/tools/sbayesrc_mcem_study07_diagnosis.R continuous baseline 200
#   Rscript research/sbayesrc/tools/sbayesrc_mcem_study07_diagnosis.R selection baseline 200

args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args)) args[[1L]] else "continuous"
start <- if (length(args) >= 2L) args[[2L]] else "baseline"
max_outer <- if (length(args) >= 3L) as.integer(args[[3L]]) else 200L
replicate_id <- if (length(args) >= 4L) as.integer(args[[4L]]) else 1L

repo_root <- normalizePath(".", winslash = "/", mustWork = TRUE)
bench_root <- normalizePath(file.path(repo_root, "..", "sblrbench"),
  winslash = "/", mustWork = TRUE)
output_dir <- file.path(repo_root, "results", "local", "sbayesrc_mcem",
  "phase5D")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages(devtools::load_all(repo_root, compile = FALSE,
  quiet = TRUE))

load_frozen_study07 <- function() {
  environment <- new.env(parent = globalenv())
  old <- setwd(bench_root)
  on.exit(setwd(old), add = TRUE)
  sys.source(file.path("R", "utils.R"), envir = environment)
  sys.source(file.path("R", "benchmark-data.R"), envir = environment)
  spec <- source(file.path("studies", "07_joint_em_sbayesrc", "spec.R"),
    local = environment)$value
  sys.source(file.path("studies", "07_joint_em_sbayesrc", "data.R"),
    envir = environment)
  # Read-only equivalent of study07_load_inputs(). The benchmark helper can
  # create block-CSR files if absent; Phase 5D must never write sblrbench.
  paths <- unlist(spec$source$paths, use.names = TRUE)
  stopifnot(all(file.exists(paths)))
  observed <- vapply(paths, environment$study07_sha256, character(1))
  stopifnot(identical(unname(observed),
    unname(spec$source$sha256[names(paths)])))
  truth <- readRDS(paths[["truth"]])
  reference <- readRDS(paths[["learned_block"]])
  glist <- readRDS(paths[["human_glist"]])
  for (field in c("bedfiles", "bimfiles", "famfiles")) {
    glist[[field]] <- vapply(glist[[field]], function(path) {
      normalizePath(file.path(bench_root, path), winslash = "/",
        mustWork = TRUE)
    }, character(1))
  }
  ld_glist <- readRDS(paths[["training_ld"]])
  rows <- as.integer(ld_glist$sparseLD$rows)
  marker_ids <- as.character(truth$marker_ids)
  selected_index <- match(marker_ids,
    as.character(ld_glist$sparseLD$marker_names))
  stopifnot(!anyNA(selected_index), !anyDuplicated(selected_index),
    identical(rownames(truth$annotations), marker_ids),
    identical(colnames(truth$annotations), spec$model$annotation_columns))
  y <- matrix(as.numeric(truth$training_y), ncol = 1L,
    dimnames = list(glist$ids[rows], "trait1"))
  working_glist <- environment$benchmark_set_training_af(glist, 1L,
    marker_ids, as.numeric(truth$gwas$freq))
  stats <- sblr::make_summary_stats(Glist = working_glist, y = y, chr = 1L,
    rows = rows, scale = TRUE, nthreads = 1L)
  native_reference <- reference$result$native_fit
  stopifnot(max(abs(as.numeric(stats$wy[[1L]]) -
    as.numeric(native_reference$wy[, 1L]))) < 1e-5)
  study06 <- new.env(parent = environment)
  sys.source(file.path("studies", "06_annotation_models", "spec.R"),
    envir = study06)
  sys.source(file.path("studies", "06_annotation_models",
    "annotation-design.R"), envir = study06)
  alpha_truth <- study06$construct_annotation_truth(truth$annotations,
    study06$spec)$informative_annotations
  prior_truth <- study06$annotation_marker_probabilities(truth$annotations,
    alpha_truth, spec$model$gamma)
  input <- native_reference$input
  intercept_spec <- input$annotation_intercept_prior
  intercept_native <- rbind(
    type = rep(0, ncol(alpha_truth)),
    mean = as.numeric(intercept_spec$mean),
    precision = as.numeric(intercept_spec$precision))
  colnames(intercept_native) <- colnames(alpha_truth)
  raw <- truth$marker_truth$raw_effect
  active <- truth$marker_truth$true_nonnull & is.finite(raw) & raw != 0
  effect_scale <- stats::median(truth$effects[active] / raw[active])
  residual <- truth$training_y - as.numeric(truth$training_x %*% truth$effects)
  truth_variance <- c(B = effect_scale^2, E = stats::var(residual),
    Vg = stats::var(as.numeric(truth$training_x %*% truth$effects)))
  truth_variance <- c(truth_variance,
    h2 = truth_variance[["Vg"]] /
      (truth_variance[["Vg"]] + truth_variance[["E"]]))
  csr_prefix_relative <- file.path(spec$output$local_dir, "inputs",
    "study06_exact_block")
  csr_suffix <- c(".row_ptr.u64.bin", ".col_idx.u32.0based.bin",
    ".values.f32.bin", ".meta.txt")
  stopifnot(all(file.exists(paste0(csr_prefix_relative, csr_suffix))))
  csr_prefix <- file.path(normalizePath(dirname(csr_prefix_relative),
    winslash = "/", mustWork = TRUE), basename(csr_prefix_relative))
  inputs <- list(
    truth = truth, marker_truth = truth$marker_truth,
    alpha_truth = alpha_truth, prior_truth = prior_truth,
    variance_truth = truth_variance, stats = stats, Glist = working_glist,
    training_rows = rows, block_start = which(!duplicated(truth$block)),
    csr_prefix = csr_prefix, alpha_init = input$alpha_init,
    intercept_spec = intercept_spec, intercept_native = intercept_native,
    source_hashes = observed)
  list(spec = spec, inputs = inputs, environment = environment)
}

frozen <- load_frozen_study07()
spec <- frozen$spec
inputs <- frozen$inputs
helper <- frozen$environment
if (!start %in% spec$em$starts) stop("Unknown Study 07 start: ", start)

run_continuous <- function() {
  alpha <- helper$study07_alpha_starts(inputs)[[start]]
  .stblr_mcem_sbayesrc_block_eigen(
    stats = inputs$stats, Glist = inputs$Glist,
    annotation = inputs$truth$annotations,
    block_start = inputs$block_start,
    B = matrix(spec$model$initial_B, 1L, 1L),
    E = matrix(spec$model$initial_E, 1L, 1L),
    ssb_prior = spec$model$ssb_prior,
    sse_prior = spec$model$sse_prior,
    gamma = spec$model$gamma,
    alpha_init = alpha,
    sigmaSqAlpha_init = spec$model$sigmaSqAlpha_fixed_em,
    intercept_prior_resolved = inputs$intercept_native,
    representation = spec$operator$representation,
    eigen_prop = spec$operator$eigen_prop,
    residual_policy = spec$operator$residual_policy,
    block_ve_mode = spec$operator$block_ve_mode,
    updateB = spec$operator$updateB,
    updateE = spec$operator$updateE,
    inner_sweeps = spec$em$inner_sweeps,
    inner_burn = spec$em$inner_burn,
    final_sweeps = spec$em$final_sweeps,
    final_burn = spec$em$final_burn,
    damping = spec$em$damping,
    tol_alpha = spec$em$tol_alpha,
    tol_prior = spec$em$tol_prior,
    min_outer = spec$em$min_outer,
    max_outer = max_outer,
    ncores = 1L,
    seed = unname(spec$em$seeds[[start]]),
    return_responsibilities = TRUE,
    verbose = FALSE,
    .objective_diagnostics = TRUE,
    .diagnostic_responsibility_iterations = c(50L, 100L, 150L, 200L))
}

run_selection <- function() {
  alpha <- helper$study07_alpha_starts(inputs)[[start]]
  delta <- helper$study07_delta_starts()[[start]]
  .stblr_mcem_sbayesrc_s_block_eigen(
    stats = inputs$stats, Glist = inputs$Glist,
    annotation = inputs$truth$annotations,
    block_start = inputs$block_start,
    B = matrix(spec$model$initial_B, 1L, 1L),
    E = matrix(spec$model$initial_E, 1L, 1L),
    ssb_prior = spec$model$ssb_prior,
    sse_prior = spec$model$sse_prior,
    gamma = spec$model$gamma,
    alpha_init = alpha,
    delta_init = delta,
    pi_a = spec$model$selection$pi_A_fixed_em,
    tau2 = spec$model$selection$tau2_fixed_em,
    intercept_prior_resolved = inputs$intercept_native,
    representation = spec$operator$representation,
    eigen_prop = spec$operator$eigen_prop,
    residual_policy = spec$operator$residual_policy,
    block_ve_mode = spec$operator$block_ve_mode,
    updateB = spec$operator$updateB,
    updateE = spec$operator$updateE,
    inner_sweeps = spec$em$inner_sweeps,
    inner_burn = spec$em$inner_burn,
    selection_sweeps = spec$em$selection_sweeps,
    selection_burn = spec$em$selection_burn,
    final_sweeps = spec$em$final_sweeps,
    final_burn = spec$em$final_burn,
    damping = spec$em$damping,
    tol_alpha = spec$em$tol_alpha,
    tol_prior = spec$em$tol_prior,
    min_outer = spec$em$min_outer,
    max_outer = max_outer,
    ncores = 1L,
    seed = unname(spec$em$selection_seeds[[start]]),
    return_responsibilities = TRUE,
    verbose = FALSE,
    .diagnostic_responsibility_iterations = c(50L, 100L, 150L, 200L))
}

run_fixed_estep <- function() {
  effort <- max_outer
  checkpoint_path <- file.path(bench_root, "results", "local",
    "07_joint_em_sbayesrc", "checkpoints",
    paste0("sbayesrc_em--", start, ".rds"))
  checkpoint <- readRDS(checkpoint_path)$result
  alpha <- checkpoint$mcem$alpha_map
  B <- as.matrix(checkpoint$mcem$genomic_hyperparameters$B_final)
  E <- as.matrix(checkpoint$mcem$genomic_hyperparameters$E_final)
  intercept <- .sbayesrc_mcem_intercept_prior(inputs$intercept_native,
    ncol(alpha))
  intercept_spec <- list(distribution = "normal", mean = intercept$mean,
    sd = 1 / sqrt(intercept$precision))
  raw <- .stblr_csr_sbayesrc_block_eigen(
    stats = inputs$stats, Glist = inputs$Glist,
    annotation = inputs$truth$annotations,
    block_start = inputs$block_start,
    representation = spec$operator$representation,
    eigen_prop = spec$operator$eigen_prop,
    residual_policy = spec$operator$residual_policy,
    block_ve_mode = spec$operator$block_ve_mode,
    block_ve_keep_history = FALSE,
    gamma = spec$model$gamma, B = B, E = E,
    ssb_prior = as.numeric(spec$model$ssb_prior),
    sse_prior = as.numeric(spec$model$sse_prior),
    updateAlpha = FALSE, updateB = spec$operator$updateB,
    updateE = spec$operator$updateE, adjE = 0,
    nit = as.integer(spec$em$inner_sweeps * effort),
    nburn = as.integer(spec$em$inner_burn * effort), nthin = 1L,
    ncores = 1L, seed = as.integer(850000L + effort * 100L + replicate_id),
    nchains = 1L, keep_chains = TRUE,
    b_init = list(as.numeric(checkpoint$genomic$marker$b[, 1L])),
    comp_init = list(as.integer(checkpoint$genomic$marker$state[, 1L])),
    use_comp_init = TRUE, use_r_init = FALSE,
    add_intercept = FALSE, standardize_annotations = FALSE,
    center_binary_annotations = FALSE, alpha_init = alpha,
    sigmaSqAlpha_init = spec$model$sigmaSqAlpha_fixed_em,
    annotation_intercept_prior = intercept_spec, pi_floor = 1e-12,
    nub = 4, nue = 4, .diagnostic_updateSigmaSqAlpha = FALSE,
    .information_diagnostics = TRUE, .return_raw = TRUE)
  responsibility <- raw$chains[[1L]][[1L]]$information_flow$rb_comp_prob
  start_alpha <- helper$study07_alpha_starts(inputs)
  mstep <- lapply(start_alpha, function(initial) {
    result <- .sbayesrc_mcem_m_step(inputs$truth$annotations,
      responsibility, initial, inputs$intercept_native,
      spec$model$sigmaSqAlpha_fixed_em)
    prior <- .sbayesrc_mcem_component_prior(inputs$truth$annotations,
      result$alpha)
    list(alpha = result$alpha, objective = sum(result$objective),
      convergence = result$convergence, component_prior = prior,
      expected_active = sum(1 - prior[, 1L]))
  })
  list(responsibility = responsibility, mstep = mstep,
    B = as.numeric(raw$variance$vb[1L, 1L]),
    E = as.numeric(raw$variance$ve[1L, 1L]))
}

started <- Sys.time()
fit <- switch(mode,
  continuous = run_continuous(),
  selection = run_selection(),
  fixed_estep = run_fixed_estep(),
  stop("Mode must be 'continuous', 'selection', or 'fixed_estep'."))
elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))

artifact <- list(
  phase = "5D",
  mode = mode,
  start = start,
  max_outer = max_outer,
  replicate_id = replicate_id,
  elapsed_seconds = elapsed,
  sblr_head = system2("git", c("rev-parse", "HEAD"), stdout = TRUE),
  benchmark_specification_sha256 =
    "241c15afab8fefc571e38e625130de6e4ab58b958c30cff369e26998ee30fa56",
  fit = fit)
path <- if (identical(mode, "fixed_estep")) {
  file.path(output_dir, sprintf("fixed-estep--%s--effort%d--rep%d.rds",
    start, max_outer, replicate_id))
} else {
  file.path(output_dir, sprintf("%s--%s--outer%d.rds", mode, start,
    max_outer))
}
saveRDS(artifact, path, compress = "xz")
message("Wrote ", path, " (", round(elapsed, 1), " seconds).")
