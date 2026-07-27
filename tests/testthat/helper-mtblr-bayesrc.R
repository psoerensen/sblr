.mt_bayesrc_fixture <- function() {
  fixture <- blr_unified_fixture()
  prefix <- blr_unified_exact_ld_prefix(fixture$dosage)
  stats <- fixture$stats
  metadata <- stats$marker_metadata
  metadata$effect_allele <- "A"
  metadata$other_allele <- "C"
  stats$marker_metadata <- metadata
  annotation <- cbind(intercept = 1, coding = c(0, 1, 0))
  rownames(annotation) <- stats$marker_names
  list(
    fixture = fixture,
    prefix = prefix,
    stats = stats,
    ld_metadata = list(
      prefix = prefix,
      marker_ids = stats$marker_names,
      marker_metadata = metadata,
      scale = "standardized_genotype",
      source = "make_summary_stats"
    ),
    annotations = annotation
  )
}

.mt_bayesrc_cleanup <- function(x) {
  blr_unified_cleanup(x$fixture)
  blr_unified_cleanup_prefix(x$prefix)
}

.mt_bayesrc_common <- function(method = "sbayesrc", alpha_init = NULL,
                               selection_s = NULL) {
  list(
    method = method,
    annotations = NULL,
    add_intercept = FALSE,
    standardize_annotations = FALSE,
    mixture_var = c(0, .1, 1),
    models = matrix(c(0L, 1L), 2L, 1L),
    alpha_init = alpha_init,
    selection_s = selection_s,
    vb = matrix(.1),
    ve = matrix(.5),
    updateB = FALSE,
    updateE = FALSE,
    updatePi = FALSE,
    updateAlpha = FALSE,
    nit = 8L,
    nburn = 2L,
    nthin = 1L,
    seed = 42L,
    convergence = "none"
  )
}
