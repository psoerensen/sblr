#' Scalar-trait BLR with sparse-LD CSR operators
#'
#' @param stats Scalar-trait summary statistics.
#' @param Glist Optional genotype/provenance list.
#' @param ld_prefix Optional CSR prefix.
#' @param method `"sbayesc"` or `"sbayesr"`; the `s` prefix denotes
#'   summary-statistics data, not MAF scaling.
#' @param nit,nburn,nthin MCMC controls.
#' @param seed Fit-local base seed.
#' @param nchains Number of logical chains per trait.
#' @param ncores Requested concurrent trait-by-chain workers.
#' @param chain_seeds Optional signed chain seeds.
#' @param keep_chains Retain compact logical-chain records.
#' @param convergence Convergence mode: `"auto"` preserves core-only automatic
#'   behavior, `"none"` disables capture, `"core"` requests the five trait-level
#'   diagnostics, and `"extended"` adds applicable low-dimensional diagnostics.
#' @param convergence_control A uniquely named list of convergence thresholds,
#'   trace retention, extended groups, explicit selected markers and quantities,
#'   and diagnostic trace-memory guards.
#' @param memory_warning_gb Analytical warning threshold.
#' @param verbose Print resolved controls.
#' @param effect_maf Optional MAF aligned to the final summary-marker order
#'   for the independent `maf_effect_s` policy.
#' @param allow_reference_maf_for_maf_effect_s Allow explicit reference-panel
#'   MAF fallback when GWAS-summary or by-construction analysis MAF is absent.
#' @param ... Model-specific controls accepted by the canonical implementation.
#' @return A `stblr_fit` and `blr_fit` object.
#' @export
stblr_csr <- function(
  stats, Glist = NULL, ld_prefix = NULL,
  method = c("sbayesc", "sbayesr"), effect_maf = NULL,
  allow_reference_maf_for_maf_effect_s = FALSE, ...,
  nit = 1000, nburn = 500, nthin = 1, seed = 1,
  nchains = 1L, ncores = 1L, chain_seeds = NULL,
  keep_chains = FALSE, convergence = c("auto", "none", "core", "extended"),
  convergence_control = NULL, memory_warning_gb = 8, verbose = FALSE
) {
  dots <- list(...)
  resolved_model <- .blr_resolve_st_model(
    method, dots, c("sbayesc", "sbayesr"), "csr")
  method <- resolved_model$model
  dots <- resolved_model$dots
  maf_info <- .blr_resolve_st_effect_maf(
    effect_maf, allow_reference_maf_for_maf_effect_s,
    resolved_model$maf_effect_s_active, stats, Glist)
  Glist <- maf_info$Glist
  chain <- .blr_chain_controls(
    nit, nburn, nthin, seed, nchains, ncores, chain_seeds, keep_chains)
  conv <- .blr_convergence_controls(convergence, convergence_control,
                                    chain$nchains)
  marker_ids <- names(stats$ww[[1L]]) %||% stats$marker_id %||%
    paste0("V", seq_along(stats$ww[[1L]]))
  trace_spec <- .blr_st_native_trace_spec(
    conv, marker_ids, resolved_model$prior_kernel,
    component_count = if (resolved_model$prior_kernel == "bayesr")
      length(dots$mixture_var %||% c(0, 0.01, 0.1, 1)) else 0L)
  memory <- .blr_st_preflight_memory(
    stats = stats, operator = "csr", chain = chain, conv = conv,
    memory_warning_gb = memory_warning_gb, trace_spec = trace_spec)
  forbidden <- intersect(names(dots), c(
    "nit", "nburn", "nthin", "seed", "nchains", "ncores",
    "chain_seeds", "keep_chains", "method"))
  if (length(forbidden)) stop("Duplicate canonical controls: ",
                              paste(forbidden, collapse = ", "),
                              call. = FALSE)
  fit <- do.call(.stblr_csr_impl, c(list(
    stats = stats, Glist = Glist, ld_prefix = ld_prefix,
    method = resolved_model$kernel,
    nit = chain$nit, nburn = chain$nburn, nthin = chain$nthin,
    seed = chain$seed, nchains = chain$nchains, ncores = chain$ncores,
    chain_seeds = if (length(chain$chain_seeds_native))
      chain$chain_seeds_native else NULL,
    keep_chains = chain$keep_chains || conv$compute || conv$keep_traces,
    .convergence_spec = trace_spec), dots))
  fit <- .blr_finalize_st_public(
    fit, method, "csr", chain, conv, memory_warning_gb, verbose, memory)
  fit$input$effect_scale <- resolved_model$effect_scale
  fit$input$prior_kernel <- resolved_model$prior_kernel
  fit$input$probability_policy <- resolved_model$probability_policy
  fit$input$effect_maf_source <- maf_info$effect_maf_source
  fit$input$effect_maf_population <- maf_info$effect_maf_population
  fit$input$effect_maf_alignment_status <-
    maf_info$effect_maf_alignment_status
  fit$input$effect_maf_fallback_used <-
    maf_info$effect_maf_fallback_used
  fit$data$effect_scale <- resolved_model$effect_scale
  fit$data$effect_maf_source <- maf_info$effect_maf_source
  fit$data$effect_maf_population <- maf_info$effect_maf_population
  fit$data$effect_maf_alignment_status <-
    maf_info$effect_maf_alignment_status
  fit$data$effect_maf_fallback_used <- maf_info$effect_maf_fallback_used
  fit
}

#' Scalar-trait BLR with packed-BED genotypes
#'
#' @param y Phenotype vector or matrix.
#' @param Glist Genotype/BED provenance.
#' @param method `"bayesc"`, `"bayesr"`, or `"bayesrc"`.
#' @inheritParams stblr_csr
#' @return A `stblr_fit` and `blr_fit` object.
#' @export
stblr_bed <- function(
  y, Glist, method = c("bayesc", "bayesr", "bayesrc"), ...,
  nit = 1000, nburn = 500, nthin = 1, seed = 1,
  nchains = 1L, ncores = 1L, chain_seeds = NULL,
  keep_chains = FALSE, convergence = c("auto", "none", "core", "extended"),
  convergence_control = NULL, memory_warning_gb = 8, verbose = FALSE
) {
  dots <- list(...)
  if (length(method) > 1L) method <- method[[1L]]
  if (!is.character(method) || length(method) != 1L || is.na(method) ||
      !method %in% c("bayesc", "bayesr", "bayesrc")) {
    .blr_public_model_error(
      method, "packed_bed", c("bayesc", "bayesr", "bayesrc"))
  }
  chain <- .blr_chain_controls(
    nit, nburn, nthin, seed, nchains, ncores, chain_seeds, keep_chains)
  conv <- .blr_convergence_controls(convergence, convergence_control,
                                    chain$nchains)
  bed_marker_ids <- unlist(Glist$rsids %||% Glist$rsidsLD, use.names = FALSE)
  if (!length(bed_marker_ids)) bed_marker_ids <- paste0(
    "V", seq_len(sum(lengths(Glist$cls))))
  trace_spec <- .blr_st_native_trace_spec(
    conv, bed_marker_ids, method,
    annotations = identical(method, "bayesrc"),
    component_count = if (method %in% c("bayesr", "bayesrc"))
      length(dots$mixture_var %||% c(0, 0.01, 0.1, 1)) else 0L,
    annotation_quantity_count = if (identical(method, "bayesrc")) {
      annotation_columns <- ncol(dots$A %||% dots$annotations %||% matrix(nrow = 0L, ncol = 0L))
      component_count <- length(dots$mixture_var %||% c(0, 0.01, 0.1, 1))
      annotation_columns * (component_count - 1L) + (component_count - 1L)
    } else 0L)
  memory <- .blr_st_preflight_memory(
    y = y, Glist = Glist, operator = "packed_bed", chain = chain,
    conv = conv, memory_warning_gb = memory_warning_gb,
    trace_spec = trace_spec)
  forbidden <- intersect(names(dots), c(
    "nit", "nburn", "nthin", "seed", "nchains", "ncores",
    "chain_seeds", "keep_chains", "method"))
  if (length(forbidden)) stop("Duplicate canonical controls: ",
                              paste(forbidden, collapse = ", "),
                              call. = FALSE)
  fit <- do.call(.stblr_bed_impl, c(list(
    y = y, Glist = Glist, method = method,
    nit = chain$nit, nburn = chain$nburn, nthin = chain$nthin,
    seed = chain$seed, nchains = chain$nchains, ncores = chain$ncores,
    chain_seeds = if (length(chain$chain_seeds_native))
      chain$chain_seeds_native else NULL,
    keep_chains = chain$keep_chains || conv$compute || conv$keep_traces,
    .convergence_spec = trace_spec), dots))
  fit <- .blr_finalize_st_public(
    fit, method, "packed_bed", chain, conv, memory_warning_gb, verbose,
    memory)
  fit$input$effect_scale <- if (method %in% c("bayesr", "bayesrc"))
    "component" else "unit"
  fit$input$prior_kernel <- method
  fit$input$probability_policy <- if (method == "bayesrc")
    "annotation_probit_stick" else "global"
  fit$data$effect_scale <- fit$input$effect_scale
  fit
}
