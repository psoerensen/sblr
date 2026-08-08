# Summarize Phase 5D artifacts. Development-only; reads sblrbench frozen
# checkpoints and writes compact CSV evidence below sblr/results/local.

repo_root <- normalizePath(".", winslash = "/", mustWork = TRUE)
bench_root <- normalizePath(file.path(repo_root, "..", "sblrbench"),
  winslash = "/", mustWork = TRUE)
output_dir <- file.path(repo_root, "results", "local", "sbayesrc_mcem",
  "phase5D")
suppressPackageStartupMessages(devtools::load_all(repo_root, compile = FALSE,
  quiet = TRUE))

old <- setwd(bench_root)
spec <- source(file.path("studies", "07_joint_em_sbayesrc", "spec.R"))$value
truth <- readRDS(spec$source$paths$truth)
reference <- readRDS(spec$source$paths$learned_block)
intercept_spec <- reference$result$native_fit$input$annotation_intercept_prior
setwd(old)
A <- as.matrix(truth$annotations)
intercept_native <- rbind(type = rep(0, 3L),
  mean = as.numeric(intercept_spec$mean),
  precision = as.numeric(intercept_spec$precision))
starts <- c("baseline", "mixed", "negative", "positive")
iterations <- c(50L, 100L, 150L, 200L)

pairwise_max <- function(values) {
  pairs <- combn(seq_along(values), 2L)
  max(apply(pairs, 2L, function(index) {
    max(abs(values[[index[[1L]]]] - values[[index[[2L]]]]))
  }))
}

pairwise_min_cor <- function(values) {
  pairs <- combn(seq_along(values), 2L)
  min(apply(pairs, 2L, function(index) {
    stats::cor(as.numeric(values[[index[[1L]]]]),
      as.numeric(values[[index[[2L]]]]))
  }))
}

load_artifact <- function(mode, start) {
  readRDS(file.path(output_dir,
    sprintf("%s--%s--outer200.rds", mode, start)))
}

continuous <- setNames(lapply(starts, function(start) {
  load_artifact("continuous", start)
}), starts)

continuous_state <- function(artifact, iteration) {
  fit <- artifact$fit
  alpha <- fit$mcem$history$alpha[, , iteration + 1L]
  prior <- .sbayesrc_mcem_component_prior(A, alpha)
  list(alpha = alpha, prior = prior, active = 1 - prior[, 1L],
    expected_active = sum(1 - prior[, 1L]))
}

continuous_trajectory <- do.call(rbind, lapply(starts, function(start) {
  artifact <- continuous[[start]]
  summary <- artifact$fit$mcem$history$summary
  rows <- summary[summary$outer %in% iterations, , drop = FALSE]
  rows$start <- start
  rows$converged <- artifact$fit$mcem$converged
  rows$completed_outer <- artifact$fit$mcem$n_outer
  rows$elapsed_seconds <- artifact$elapsed_seconds
  rows
}))
utils::write.csv(continuous_trajectory,
  file.path(output_dir, "continuous_outer_trajectory.csv"), row.names = FALSE)

continuous_comparison <- do.call(rbind, lapply(iterations, function(iteration) {
  state <- lapply(continuous, continuous_state, iteration = iteration)
  data.frame(
    outer = iteration,
    max_alpha_difference = pairwise_max(lapply(state, `[[`, "alpha")),
    max_component_prior_difference = pairwise_max(lapply(state, `[[`, "prior")),
    max_active_prior_difference = pairwise_max(lapply(state, `[[`, "active")),
    min_active_prior_correlation = pairwise_min_cor(lapply(state, `[[`, "active")),
    expected_active_range = diff(range(vapply(state, `[[`, numeric(1L),
      "expected_active"))))
}))

final_genomic <- lapply(continuous, function(artifact) artifact$fit$genomic)
continuous_final <- data.frame(
  metric = c("genomic_pip", "beta"),
  max_difference = c(
    pairwise_max(lapply(final_genomic, function(x) x$marker$dm[, 1L])),
    pairwise_max(lapply(final_genomic, function(x) x$marker$bm[, 1L]))),
  min_correlation = c(
    pairwise_min_cor(lapply(final_genomic, function(x) x$marker$dm[, 1L])),
    pairwise_min_cor(lapply(final_genomic, function(x) x$marker$bm[, 1L])))
)
utils::write.csv(continuous_comparison,
  file.path(output_dir, "continuous_between_start.csv"), row.names = FALSE)
utils::write.csv(continuous_final,
  file.path(output_dir, "continuous_final_genomic.csv"), row.names = FALSE)

reproduction <- do.call(rbind, lapply(starts, function(start) {
  frozen <- readRDS(file.path(bench_root, "results", "local",
    "07_joint_em_sbayesrc", "checkpoints",
    paste0("sbayesrc_em--", start, ".rds")))$result
  long <- continuous[[start]]$fit
  data.frame(
    start = start,
    max_alpha_difference = max(abs(
      frozen$mcem$alpha_map - long$mcem$history$alpha[, , 51L])),
    max_history_difference = max(abs(
      as.matrix(frozen$mcem$history$summary) -
        as.matrix(long$mcem$history$summary[seq_len(50L),
          names(frozen$mcem$history$summary), drop = FALSE]))),
    stringsAsFactors = FALSE)
}))
utils::write.csv(reproduction,
  file.path(output_dir, "continuous_iteration50_reproduction.csv"),
  row.names = FALSE)

fixed_paths <- list.files(output_dir, pattern = "^fixed-estep--.*\\.rds$",
  full.names = TRUE)
if (length(fixed_paths)) {
  fixed <- lapply(fixed_paths, readRDS)
  fixed_summary <- do.call(rbind, lapply(fixed, function(x) {
    baseline <- x$fit$mstep[["baseline"]]
    data.frame(start = x$start, effort = x$max_outer,
      replicate = x$replicate_id, B = x$fit$B, E = x$fit$E,
      expected_active = baseline$expected_active,
      objective = baseline$objective,
      stringsAsFactors = FALSE)
  }))
  utils::write.csv(fixed_summary,
    file.path(output_dir, "fixed_estep_summary.csv"), row.names = FALSE)
}

selection_paths <- file.path(output_dir,
  sprintf("selection--%s--outer200.rds", starts))
if (all(file.exists(selection_paths))) {
  selection <- setNames(lapply(selection_paths, readRDS), starts)
  exact_rows <- list()
  row <- 0L
  models <- as.matrix(expand.grid(rep(list(0:1), ncol(A) - 1L)))
  for (start in starts) {
    fit <- selection[[start]]$fit$mcem
    for (iteration in iterations) {
      responsibility <- fit$history$responsibility_checkpoint[[as.character(iteration)]]
      if (is.null(responsibility)) next
      alpha_start <- fit$history$alpha[[iteration]]
      evaluated <- lapply(seq_len(nrow(models)), function(index) {
        .sbayesrc_s_em_model_laplace(A, responsibility, models[index, ],
          fit$pi_A_fixed, fit$tau2_fixed,
          intercept_native,
          alpha_start)
      })
      log_weight <- vapply(evaluated, `[[`, numeric(1L), "log_weight")
      weight <- exp(log_weight - max(log_weight)); weight <- weight / sum(weight)
      exact_pip <- drop(crossprod(weight, models))
      production_pip <- fit$history$annotation_pip_eb[[iteration]]
      row <- row + 1L
      exact_rows[[row]] <- data.frame(start = start, outer = iteration,
        annotation = colnames(A)[-1L], exact_pip = exact_pip,
        mc3_pip = production_pip,
        absolute_difference = abs(exact_pip - production_pip),
        best_model = paste(models[which.max(weight), ], collapse = ""),
        best_model_probability = max(weight), stringsAsFactors = FALSE)
    }
  }
  utils::write.csv(do.call(rbind, exact_rows),
    file.path(output_dir, "selection_exact_model_audit.csv"), row.names = FALSE)
}

message("Phase 5D summaries written to ", output_dir)
