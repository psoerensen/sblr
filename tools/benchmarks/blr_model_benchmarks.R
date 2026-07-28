pkgload::load_all(".", compile = FALSE, quiet = TRUE)
source(file.path("tests", "testthat", "helper-blr-unified.R"), local = TRUE)
fixture <- blr_unified_fixture()
on.exit(blr_unified_cleanup(fixture), add = TRUE)
prefix <- blr_unified_exact_ld_prefix(fixture$dosage)
on.exit(blr_unified_cleanup_prefix(prefix), add = TRUE)
metadata <- fixture$stats$marker_metadata
metadata$effect_allele <- "A"
metadata$other_allele <- "C"
fixture$stats$marker_metadata <- metadata
ld_metadata <- list(prefix = prefix, marker_ids = fixture$stats$marker_names,
  marker_metadata = metadata, scale = "standardized_genotype",
  source = "make_summary_stats")
annotations <- cbind(intercept = 1, binary = c(0, 1, 0),
                     continuous = c(-1, 0, 1))
rownames(annotations) <- fixture$stats$marker_names
common <- list(annotations = annotations, add_intercept = FALSE,
  standardize_annotations = FALSE, mixture_var = c(0, .01, .1, 1),
  models = matrix(c(0L, 1L), 2L, 1L), updateB = FALSE, updateE = FALSE,
  nit = 20L, nburn = 5L, nthin = 1L, convergence = "none",
  memory_warning_gb = 1e9)
cases <- list(
  csr = list(fun = mtblr_csr, method = "sbayesrc", args = list(
    stats = fixture$stats, ld_prefix = prefix, ld_metadata = ld_metadata)),
  block_eigen = list(fun = mtblr_block_eigen, method = "sbayesrc",
    args = list(stats = fixture$stats, Glist = fixture$Glist,
      block_start = 1L, eigen_filter = "hard_truncate", eigen_tau = 0)),
  packed_bed = list(fun = mtblr_bed, method = "bayesrc", args = list(
    y = fixture$y, Glist = fixture$Glist,
    residual_covariance = "diagonal")))
rows <- lapply(names(cases), function(operator) {
  entry <- cases[[operator]]
  annotation_time <- system.time(control <- sblr:::.mtblr_bayesrc_controls(
    "bayesrc", annotations, fixture$stats$marker_names,
    sblr:::.mtblr_models(matrix(c(0L, 1L), 2L, 1L), c(.8, .2), .2, 1L),
    sblr:::.mtblr_bayesr_spec("bayesr",
      sblr:::.mtblr_models(matrix(c(0L, 1L), 2L, 1L), c(.8, .2), .2, 1L),
      NULL, 3L, c(0, .01, .1, 1)), add_intercept = FALSE,
    standardize_annotations = FALSE))[["elapsed"]]
  elapsed <- system.time(fit <- do.call(entry$fun, c(entry$args, common,
    list(method = entry$method))))[["elapsed"]]
  memory <- fit$memory_estimate
  data.frame(operator = operator, model = entry$method, traits = ncol(fit$bm),
    annotation_columns = ncol(control$annotations), positive_components = 3L,
    chains = 1L, cores = 1L,
    annotation_preprocessing_seconds = annotation_time,
    total_elapsed_seconds = elapsed,
    shared_annotation_bytes = memory$bayesrc_shared_annotation_bytes,
    private_annotation_bytes_per_worker =
      memory$bayesrc_private_annotation_bytes_per_worker,
    chain_annotation_result_bytes = memory$bayesrc_chain_annotation_result_bytes,
    prior_probability_output_bytes =
      memory$bayesrc_formatted_annotation_output_bytes,
    fit_bytes = as.numeric(object.size(fit)),
    convergence_bytes = as.numeric(object.size(fit$convergence)))
})
print(do.call(rbind, rows), row.names = FALSE)
cat("model benchmark benchmark values are deterministic regression signals; they are not production-runtime or linear-speedup claims.\n")
