pkgload::load_all(".", compile = FALSE, quiet = TRUE)
source(file.path("tests", "testthat", "helper-blr-unified.R"), local = TRUE)
fixture <- blr_unified_fixture()
on.exit(blr_unified_cleanup(fixture), add = TRUE)
prefix <- blr_unified_exact_ld_prefix(fixture$dosage)
on.exit(blr_unified_cleanup_prefix(prefix), add = TRUE)
metadata <- fixture$stats$marker_metadata
metadata$effect_allele <- "A"; metadata$other_allele <- "C"
fixture$stats$marker_metadata <- metadata
ld_metadata <- list(prefix = prefix, marker_ids = fixture$stats$marker_names,
  marker_metadata = metadata, scale = "standardized_genotype",
  source = "make_summary_stats")
common <- list(mixture_var = c(0, .01, .1, 1),
  models = matrix(c(0L, 1L), 2L, 1L), joint_pi = c(.8, rep(.2 / 3, 3)),
  updateB = FALSE, updateE = FALSE, nit = 20L, nburn = 5L,
  nthin = 1L, convergence = "none", memory_warning_gb = 1e9)
cases <- list(
  csr = list(fun = mtblr_csr, args = list(stats = fixture$stats,
    ld_prefix = prefix, ld_metadata = ld_metadata)),
  block_eigen = list(fun = mtblr_block_eigen, args = list(
    stats = fixture$stats, Glist = fixture$Glist, block_start = 1L,
    eigen_filter = "hard_truncate", eigen_tau = 0)),
  packed_bed = list(fun = mtblr_bed, args = list(y = fixture$y,
    Glist = fixture$Glist, residual_covariance = "diagonal")))
rows <- list(); index <- 1L
for (model in c("bayesr", "sbayesr")) for (operator in names(cases)) {
  entry <- cases[[operator]]
  args <- c(entry$args, common, list(method = model,
    selection_s = if (model == "sbayesr") -1 else NULL))
  elapsed <- system.time(fit <- do.call(entry$fun, args))[["elapsed"]]
  rows[[index]] <- data.frame(operator = operator, model = model,
    traits = ncol(fit$bm), patterns = nrow(fit$model_patterns),
    positive_components = length(common$mixture_var) - 1L,
    chains = 1L, cores = 1L, total_elapsed_seconds = elapsed,
    analytical_total_bytes = fit$memory_estimate$estimated_total_bytes,
    private_worker_bytes = fit$memory_estimate$bayesr_private_worker_state_bytes,
    chain_result_bytes = fit$memory_estimate$bayesr_chain_result_bytes,
    component_output_bytes = fit$memory_estimate$bayesr_component_output_bytes,
    fit_bytes = as.numeric(object.size(fit)))
  index <- index + 1L
}
result <- do.call(rbind, rows)
print(result, row.names = FALSE)
cat("Benchmark values are regression signals, not claims of linear speedup or pure adapter overhead.\n")
