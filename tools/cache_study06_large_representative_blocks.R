#!/usr/bin/env Rscript

# Cache five evenly spaced retained factors from the frozen 76-block large
# Study 06 design. The sibling repository is read only.

sblr_root <- normalizePath(".", winslash = "/", mustWork = TRUE)
bench_root <- normalizePath("../sblrbench", winslash = "/", mustWork = TRUE)
output_root <- file.path(sblr_root, "results", "local",
                         "sbayesrc_particle_marginal_alpha")
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
cache_path <- file.path(output_root, "study06_large_representative_blocks.rds")
if (file.exists(cache_path)) quit(save = "no", status = 0L)
pkgload::load_all(sblr_root, quiet = TRUE)

bundle <- readRDS(file.path(bench_root, "results", "local",
  "06_annotation_models", "large_feasibility", "prepared_bundle.rds"))
captured <- readRDS(file.path(bench_root, "results", "local",
  "06_annotation_models", "large_feasibility", "continuation",
  "b0_iteration0_no_updateE.rds"))
smokes <- readRDS(file.path(sblr_root, "results", "local",
  "block_gctb_residual_policy", "large_b0_smokes.rds"))
st <- sblr:::.mtblr_normalize_stats(bundle$gwas$stats)
provenance <- st$genotype_provenance[[1L]]
provenance$cls <- unname(provenance$cls)
provenance$af <- unname(provenance$af)
old_directory <- setwd(bench_root)
reference <- sblr:::.mtblr_block_eigen_reference(bundle$glist, provenance)
contract <- sblr:::stblr_block_low_rank_contract_internal(
  reference$bed_files, reference$n_bed, reference$cls, reference$rows,
  reference$af, as.integer(bundle$blocks$block_start - 1L),
  matrix(st$wy[[1L]], nrow = 1L), as.numeric(captured$b[, 1L]),
  bundle$spec$block$eigen_prop, st$yy[[1L]], 0)
setwd(old_directory)

selected <- c(1L, 20L, 39L, 58L, 76L)
starts <- bundle$blocks$block_start
ends <- c(starts[-1L] - 1L, nrow(bundle$annotations))
learned <- smokes$learned_alpha
cache <- list(
  block_id = selected,
  factor = contract$factor[selected],
  transformed_score = contract$transformed_score[selected],
  annotations = lapply(selected, function(block)
    bundle$annotations[starts[[block]]:ends[[block]], , drop = FALSE]),
  gamma = as.numeric(captured$mixture_var),
  alpha_truth = bundle$calibrated$alpha,
  alpha_learned = learned$alpha_final[[1L]],
  sigma_sq_alpha = as.numeric(learned$sigmaSqAlpha_final[1L, ]),
  vb = mean(learned$vbs),
  block_ve = colMeans(learned$block_ve$final_per_chain_block)[selected],
  block_count = length(contract$factor), n = st$n, yy = st$yy[[1L]])
saveRDS(cache, cache_path, compress = FALSE)
cat("Cached representative large blocks:", paste(selected, collapse = ", "), "\n")
