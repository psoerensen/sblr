#!/usr/bin/env Rscript

# Package-side validation for the official-compatible block residual policy.
# All sibling inputs are read-only; outputs stay under this repository's ignored
# results/local tree. This script is not a Study 06 benchmark entry point.

`%||%` <- function(x, y) if (is.null(x) || !length(x)) y else x

repo_root <- normalizePath(Sys.getenv("SBLR_ROOT", "."), winslash = "/",
                           mustWork = TRUE)
bench_root <- normalizePath(Sys.getenv("SBLRBENCH_ROOT", "../sblrbench"),
                            winslash = "/", mustWork = TRUE)
output_root <- file.path(repo_root, "results", "local",
                         "block_gctb_residual_policy")
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
pkgload::load_all(repo_root, quiet = TRUE)
# The preserved large-study Glist intentionally stores sibling-root-relative
# BED paths. Resolve those read-only paths without rewriting the object.
setwd(bench_root)

write_json <- function(x, name) {
  jsonlite::write_json(x, file.path(output_root, name), auto_unbox = TRUE,
                       pretty = TRUE, digits = 17, null = "null")
}

write_bed <- function(path, dosage) {
  dosage_to_code <- c(`0` = 3L, `1` = 2L, `2` = 0L)
  packed <- unlist(lapply(seq_len(ncol(dosage)), function(marker) {
    codes <- unname(dosage_to_code[as.character(dosage[, marker])])
    codes <- c(codes, rep(0L, (-length(codes)) %% 4L))
    vapply(seq(1L, length(codes), by = 4L), function(i)
      sum(codes[i:(i + 3L)] * c(1L, 4L, 16L, 64L)), integer(1))
  }), use.names = FALSE)
  writeBin(as.raw(c(0x6c, 0x1b, 0x01, packed)), path)
}

summarise_block_fit <- function(fit) {
  list(
    residual_policy = fit$input$residual_policy,
    block_ve_mode = fit$input$block_ve_mode,
    mean_block_ve = mean(fit$ves),
    final_mean_block_ve = mean(fit$block_ve$final_per_chain_block),
    mean_summary_heritability = mean(fit$heritability_summary),
    final_active = sum(fit$component[, 1L] > 0),
    block_resampled = sum(fit$block_ve$resampled_per_chain_block),
    minimum_ratio_resets = sum(fit$block_ve$reset_to_phenotype_per_chain_block),
    finite = all(is.finite(c(fit$b, fit$vbs, fit$vgs, fit$ves,
                             fit$heritability_summary,
                             fit$block_ve$final_per_chain_block))))
}

run_large_b0 <- function() {
  summary_path <- file.path(output_root, "large_b0_trace_validation_summary.json")
  if (file.exists(summary_path))
    return(jsonlite::read_json(summary_path, simplifyVector = TRUE))
  continuation <- file.path(
    bench_root, "results/local/06_annotation_models/large_feasibility/continuation")
  bundle_path <- file.path(
    bench_root, "results/local/06_annotation_models/large_feasibility/prepared_bundle.rds")
  captured_path <- file.path(continuation, "b0_iteration0_no_updateE.rds")
  stopifnot(file.exists(bundle_path), file.exists(captured_path))
  bundle <- readRDS(bundle_path)
  captured <- readRDS(captured_path)
  stopifnot(identical(bundle$identities$specification_hash,
                      "b001bc36a5531e5e6b342286a253fc1fd34dad4265359d89d2feaa026d4533df"),
            length(captured$b[, 1L]) == 37991L)

  st <- sblr:::.mtblr_normalize_stats(bundle$gwas$stats)
  provenance <- st$genotype_provenance[[1L]]
  provenance$cls <- unname(provenance$cls)
  provenance$af <- unname(provenance$af)
  reference <- sblr:::.mtblr_block_eigen_reference(bundle$glist, provenance)
  beta <- as.numeric(captured$b[, 1L])
  contract <- sblr:::stblr_block_low_rank_contract_internal(
    reference$bed_files, reference$n_bed, reference$cls, reference$rows,
    reference$af, as.integer(bundle$blocks$block_start - 1L),
    matrix(st$wy[[1L]], nrow = 1L), beta,
    bundle$spec$block$eigen_prop, st$yy[[1L]], 0)
  vy <- st$yy[[1L]] / (st$n - 1)
  nue <- bundle$spec$prior$nue %||% 4
  prior_scale <- (nue - 2) / nue * vy
  residual_ss <- as.numeric(contract$block_residual_norm_squared)
  draw_scale <- residual_ss + nue * prior_scale
  rank <- vapply(contract$factor, nrow, integer(1))
  vg <- as.numeric(contract$block_genetic_variance)
  effect_ss <- as.numeric(contract$block_effect_ss)
  selected <- vapply(seq_along(vg), function(i)
    sblr:::.stblr_block_ve_decision(
      "allMixVe", 0L, vg[[i]], effect_ss[[i]], 1.1)$selected, logical(1))
  set.seed(bundle$spec$seeds$smoke[[1L]])
  proposal <- rep(vy, length(rank))
  proposal[selected] <- draw_scale[selected] /
    stats::rchisq(sum(selected), df = rank[selected] + nue)
  reset <- selected & !(proposal / vy > 0.7)
  proposal[reset] <- vy

  legacy_json <- jsonlite::read_json(
    file.path(continuation, "b0_residual_scale_audit.json"), simplifyVector = TRUE)
  deterministic <- list(
    specification_hash = bundle$identities$specification_hash,
    truth_hash = bundle$identities$truth_hash,
    n = st$n, markers = length(beta), blocks = length(rank),
    captured_active = sum(captured$component > 0),
    legacy_global_residual_scale = legacy_json$terms$block_residual_scale,
    local_scale = list(minimum = min(draw_scale), median = stats::median(draw_scale),
                       maximum = max(draw_scale), all_finite = all(is.finite(draw_scale)),
                       all_nonnegative = all(draw_scale >= 0)),
    selected_blocks = sum(selected), reset_blocks = sum(reset),
    initial_mean_block_ve = vy, resulting_mean_block_ve = mean(proposal))

  common <- list(
    stats = bundle$gwas$stats, Glist = bundle$glist,
    block_start = bundle$blocks$block_start, representation = "low_rank",
    eigen_policy = bundle$spec$block$eigen_policy,
    eigen_prop = bundle$spec$block$eigen_prop, residual_policy = "gctb_block",
    block_ve_mode = "allMixVe", h2 = bundle$spec$prior$h2,
    nit = 12L, nburn = 4L, nthin = 1L,
    seed = bundle$spec$seeds$fit, nchains = 4L, ncores = 4L,
    chain_seeds = bundle$spec$seeds$smoke, keep_chains = TRUE,
    updateB = TRUE, updateE = TRUE)
  alpha_truth <- bundle$calibrated$alpha
  neutral <- alpha_truth * 0
  target <- bundle$spec$mixture$target_pi
  neutral[1L, ] <- c(stats::qnorm(1 - target[["null"]]),
    stats::qnorm((target[["medium"]] + target[["large"]]) /
      (1 - target[["null"]])),
    stats::qnorm(target[["large"]] /
      (target[["medium"]] + target[["large"]])))
  annotation_common <- list(
    method = "sbayesrc", gamma = bundle$spec$prior$mixture_var,
    annotation = bundle$annotations,
    sigmaSqAlpha_init = bundle$spec$prior$sigmaSqAlpha_init,
    sigmaSqAlpha_a = bundle$spec$prior$sigmaSqAlpha_a,
    sigmaSqAlpha_b = bundle$spec$prior$sigmaSqAlpha_b,
    pi_floor = bundle$spec$prior$pi_floor, alpha_update_every = 1L,
    add_intercept = FALSE, standardize_annotations = FALSE,
    center_binary_annotations = FALSE, convergence = "extended",
    convergence_control = list(extended_groups = "annotations",
                               keep_traces = TRUE))
  smoke_path <- file.path(output_root, "large_b0_trace_smokes.rds")
  if (file.exists(smoke_path)) {
    smoke <- readRDS(smoke_path)
    baseline <- smoke$baseline
    fixed <- smoke$fixed_alpha
    learned <- smoke$learned_alpha
  } else {
    baseline <- do.call(sblr::stblr_block_eigen, c(common, list(
      method = "sbayesr", mixture_var = bundle$spec$prior$mixture_var,
      pi = bundle$spec$prior$pi, convergence = "none")))
    fixed <- do.call(sblr::stblr_block_eigen, c(common, annotation_common,
      list(alpha_init = alpha_truth, updateAlpha = FALSE)))
    learned <- do.call(sblr::stblr_block_eigen, c(common, annotation_common,
      list(alpha_init = neutral, updateAlpha = TRUE)))
    saveRDS(list(baseline = baseline, fixed_alpha = fixed, learned_alpha = learned),
            smoke_path, compress = "xz")
  }
  trace_rows <- function(fit, field) {
    vapply(fit$chains, function(chain) {
      value <- chain$convergence_trace[[field]]
      if (is.null(dim(value))) length(value) else dim(value)[1L]
    }, integer(1L))
  }
  fixed_alpha_rows <- trace_rows(fixed, "alpha")
  learned_alpha_rows <- trace_rows(learned, "alpha")
  fixed_sigma_rows <- trace_rows(fixed, "sigmaSqAlpha")
  learned_sigma_rows <- trace_rows(learned, "sigmaSqAlpha")
  stopifnot(all(fixed_alpha_rows == 12L), all(learned_alpha_rows == 12L),
            all(fixed_sigma_rows == 12L), all(learned_sigma_rows == 12L))
  result <- list(deterministic = deterministic,
       baseline = summarise_block_fit(baseline),
       fixed_alpha = c(summarise_block_fit(fixed), list(
         alpha_trace_rows_per_chain = fixed_alpha_rows,
         sigmaSqAlpha_trace_rows_per_chain = fixed_sigma_rows)),
       learned_alpha = c(summarise_block_fit(learned), list(
         alpha_trace_rows_per_chain = learned_alpha_rows,
         sigmaSqAlpha_trace_rows_per_chain = learned_sigma_rows)))
  write_json(result, "large_b0_trace_validation_summary.json")
  result
}

metric_cor <- function(x, y) unname(stats::cor(as.numeric(x), as.numeric(y)))

run_official_comparison <- function() {
  export <- file.path(
    bench_root, "results/local/06_annotation_models/gctb_parity/export")
  official <- file.path(
    bench_root, "results/local/06_annotation_models/gctb_single_trajectory/runs")
  truth <- readRDS(file.path(export, "truth.rds"))
  stopifnot(nrow(truth$training_x) == 1400L, ncol(truth$training_x) == 1500L)
  p <- truth$gwas$freq
  dosage <- sweep(sweep(truth$training_x, 2L, sqrt(2 * p * (1 - p)), "*"),
                  2L, 2 * p, "+")
  stopifnot(max(abs(dosage - round(dosage))) < 1e-12)
  dosage <- round(dosage)
  bed <- file.path(output_root, "official_export_training.bed")
  write_bed(bed, dosage)
  marker <- truth$marker_ids
  ids <- paste0("train", seq_len(nrow(dosage)))
  glist <- list(n = nrow(dosage), ids = ids, bedfiles = bed,
    rsids = list(marker), rsidsLD = list(marker), chr = list(rep(1L, length(marker))),
    pos = list(seq_along(marker)), af = list(colMeans(dosage) / 2),
    maf = list(pmin(colMeans(dosage) / 2, 1 - colMeans(dosage) / 2)))
  y <- truth$training_y
  stats <- list(
    wy = list(T1 = stats::setNames(drop(crossprod(truth$training_x, y)), marker)),
    ww = list(T1 = stats::setNames(colSums(truth$training_x^2), marker)),
    yy = c(T1 = sum(y^2)), n = nrow(dosage), m = length(marker),
    marker_names = marker, trait_names = "T1", bed_files = bed,
    cls = list(seq_along(marker)), rows = seq_len(nrow(dosage)),
    af = list(colMeans(dosage) / 2))
  common <- list(stats = stats, Glist = glist,
    block_start = seq.int(1L, length(marker), by = 100L),
    representation = "low_rank", eigen_policy = "cumulative_positive_mass",
    eigen_prop = 0.995, residual_policy = "gctb_block",
    block_ve_mode = "allMixVe", h2 = 0.5, nit = 6000L, nburn = 3000L,
    nthin = 1L, nchains = 1L, ncores = 1L, keep_chains = TRUE,
    convergence = "none", updateB = TRUE, updateE = TRUE)
  d0_path <- file.path(output_root, "official_comparison_sblr_D0.rds")
  d1_path <- file.path(output_root, "official_comparison_sblr_D1.rds")
  if (file.exists(d0_path)) d0 <- readRDS(d0_path) else {
    d0 <- do.call(sblr::stblr_block_eigen, c(common,
      list(method = "sbayesr", mixture_var = c(0, 0.01, 0.1, 1),
        pi = c(0.88, 0.06, 0.036, 0.024), seed = 711121L)))
    saveRDS(d0, d0_path, compress = "xz")
  }
  if (file.exists(d1_path)) d1 <- readRDS(d1_path) else {
    d1 <- do.call(sblr::stblr_block_eigen, c(common,
      list(method = "sbayesrc", gamma = c(0, 0.01, 0.1, 1),
        annotation = truth$annotations, pi_init = 0.12,
        active_comp_weights = c(0.5, 0.3, 0.2),
        add_intercept = FALSE, standardize_annotations = FALSE,
        center_binary_annotations = FALSE, alpha_update_every = 1L,
        seed = 721121L)))
    saveRDS(d1, d1_path, compress = "xz")
  }
  saveRDS(list(D0 = d0, D1 = d1),
          file.path(output_root, "official_comparison_sblr_fits.rds"), compress = "xz")
  compare_one <- function(id, fit) {
    prefix <- file.path(official, id, id)
    off_snp <- utils::read.table(paste0(prefix, ".txt"), header = TRUE,
                                 stringsAsFactors = FALSE)
    off <- readRDS(paste0(prefix, ".rds"))
    keep <- 3001:9000
    pip <- as.numeric(fit$dm[, 1L]); effect <- as.numeric(fit$bm[, 1L])
    off_pip <- off_snp$PIP; off_effect <- off_snp$BETA
    pred <- drop(truth$validation_x %*% effect)
    phenotype_variance <- sum(y^2) / (length(y) - 1)
    effect_scale <- sqrt(phenotype_variance)
    off_effect_sblr_units <- off_effect * effect_scale
    list(
      sblr_phenotype_variance = phenotype_variance,
      pip_correlation = metric_cor(pip, off_pip),
      effect_correlation = metric_cor(effect, off_effect),
      effect_rmse_after_scale_conversion = sqrt(mean(
        (effect - off_effect_sblr_units)^2)),
      validation_genetic_agreement = metric_cor(
        pred, drop(truth$validation_x %*% off_effect)),
      posterior_mean_active = sum(fit$ncomp[1L, -1L]),
      official_mean_active = mean(rowSums(off$n_comp_hist[keep, -1L, drop = FALSE])),
      summary_heritability = mean(fit$heritability_summary[keep, 1L]),
      official_heritability = mean(off$hsq_hist[keep]),
      total_block_genetic_variance = mean(fit$vgs[keep, 1L]),
      official_total_block_genetic_variance = mean(rowSums(off$hsq_block_hist)),
      official_total_block_genetic_variance_sblr_units =
        mean(rowSums(off$hsq_block_hist)) * phenotype_variance,
      mean_block_ve = mean(fit$ves[keep, 1L]),
      official_mean_block_ve = mean(off$vare_hist),
      official_mean_block_ve_sblr_units =
        mean(off$vare_hist) * phenotype_variance,
      validation_genetic_correlation = metric_cor(pred, truth$validation_genetic),
      phenotype_prediction_correlation = metric_cor(pred, truth$validation_phenotype))
  }
  list(D0 = compare_one("D0", d0), D1 = compare_one("D1", d1))
}

result <- list(
  schema = "sblr-block-gctb-residual-policy-validation-v1",
  package_sha = system2("git", c("-C", shQuote(repo_root), "rev-parse", "HEAD"),
                        stdout = TRUE),
  sblrbench_sha = system2("git", c("-C", shQuote(bench_root), "rev-parse", "HEAD"),
                          stdout = TRUE),
  created_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
  large_b0 = run_large_b0(),
  official_comparison = run_official_comparison())
write_json(result, "validation_summary.json")
print(result)
