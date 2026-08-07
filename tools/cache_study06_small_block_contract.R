#!/usr/bin/env Rscript

# Reconstruct the immutable 1,500-marker Study 06 block likelihood inputs and
# cache only the compact transformed block factors below ignored local results.

sblr_root <- normalizePath(".", winslash = "/", mustWork = TRUE)
bench_root <- normalizePath("../sblrbench", winslash = "/", mustWork = TRUE)
output_root <- file.path(sblr_root, "results", "local",
                         "sbayesrc_particle_marginal_alpha")
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
cache_path <- file.path(output_root, "study06_small_block_contract.rds")
if (file.exists(cache_path)) quit(save = "no", status = 0L)

isolated_library <- file.path(
  bench_root, "results", "local", "current_benchmark_refresh", "rlib")
stopifnot(dir.exists(isolated_library))
.libPaths(c(isolated_library, .libPaths()))
pkgload::load_all(bench_root, quiet = TRUE)
source(file.path(bench_root, "studies", "06_annotation_models",
                 "power-isolation.R"), local = FALSE)

spec_hash <- "241c15afab8fefc571e38e625130de6e4ab58b958c30cff369e26998ee30fa56"
truth_hash <- "169d52bef390022a9106d7e61b200493869b40cbb76f1c8d911eebbc80fea1eb"
spec <- read_benchmark_spec(file.path(
  bench_root, "studies", "06_annotation_models", "spec.R"))
stopifnot(identical(benchmark_annotation_spec_hash(spec), spec_hash))
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
  seed_grid$scenario == "informative_annotations" & seed_grid$replicate == 1L &
    seed_grid$method == "st_bed_bayesrc", , drop = FALSE]
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
stopifnot(identical(
  benchmark_hash_object(study06_power_truth_identity(data, bundle)), truth_hash))

# The retained operator and summary statistics use the same selected BED rows;
# normalize both provenance paths after hashing the registered truth and before
# the package-side operator identity check.
if (length(stats$bed_files) == length(data$ld_glist$bedfiles))
  stats$bed_files <- data$ld_glist$bedfiles

pkgload::load_all(sblr_root, quiet = TRUE)
st <- sblr:::.mtblr_normalize_stats(stats)
provenance <- st$genotype_provenance[[1L]]
provenance$cls <- unname(provenance$cls)
provenance$af <- unname(provenance$af)
operator_glist <- data$ld_glist
if (!is.null(operator_glist$sparseLD))
  operator_glist$sparseLD$rows <- provenance$rows
if (!is.null(operator_glist$ids) && !is.null(provenance$rows))
  operator_glist$idsLD <- operator_glist$ids[provenance$rows]
old_directory <- setwd(bench_root)
reference <- sblr:::.mtblr_block_eigen_reference(operator_glist, provenance)
contract <- sblr:::stblr_block_low_rank_contract_internal(
  reference$bed_files, reference$n_bed, reference$cls, reference$rows,
  reference$af, as.integer(data$block_start - 1L),
  matrix(st$wy[[1L]], nrow = 1L), rep(0, nrow(annotations)),
  0.995,
  st$yy[[1L]], 0)
setwd(old_directory)

fit <- readRDS(file.path(sblr_root, "results", "local",
  "study06_kernel_composition_audit", "block_eigen_PX_screen_fit_result.rds"))
native <- fit$result$native_fit
cache <- list(
  factor = contract$factor,
  transformed_score = contract$transformed_score,
  annotations = annotations,
  gamma = as.numeric(spec$controls$simulation$mixture_var),
  alpha_truth = annotation_truth$alpha,
  alpha_pooled = native$alpha_final[[1L]],
  sigma_sq_alpha = as.numeric(native$sigmaSqAlpha_final[1L, ]),
  vb = mean(native$vbs),
  block_ve = colMeans(native$block_ve$final_per_chain_block),
  n = st$n, yy = st$yy[[1L]],
  specification_hash = spec_hash, truth_hash = truth_hash)
saveRDS(cache, cache_path, compress = FALSE)
cat("Cached", length(cache$factor), "Study 06 blocks at", cache_path, "\n")
