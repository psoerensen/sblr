#' Scalar-trait BLR with a block-eigen LD operator
#'
#' Fits one independent scalar model per trait using the canonical retained
#' low-rank block-eigen operator. The reconstructed-dense historical operator is
#' retained as an explicit representation for regression and reproducibility.
#' The retained factor follows the GCTB/SBayesRC eigenspace likelihood strategy
#' in `sblr` cross-product units, but uses `sblr`'s global projected
#' residual-variance contract rather than GCTB's block-specific variance
#' procedure.
#'
#' @param stats Scalar-trait summary statistics.
#' @param Glist Genotype/BED provenance used to construct the operator.
#' @param block_start One-based public block starts.
#' @param method One of `"sbayesc"`, `"sbayesr"`, or `"sbayesrc"`;
#'   the `s` prefix denotes summary-statistics data.
#' @param effect_maf Optional marker-aligned allele frequencies used only
#'   when the independent `maf_effect_s` variance-scaling option is active.
#' @param allow_reference_maf_for_maf_effect_s Whether aligned reference-panel
#'   frequencies may be used explicitly when summary-population frequencies
#'   are unavailable. The default is `FALSE`.
#' @param annotation Annotation matrix required for `"sbayesrc"`.
#' @param representation Operator representation. `"low_rank"` is canonical;
#'   `"dense_reconstructed"` retains the historical packed dense operator.
#' @param eigen_policy Representation-specific eigenvalue policy, or `NULL` for
#'   the representation default.
#' @param eigen_prop Cumulative positive-eigenvalue mass target for low rank.
#' @param eigen_tau,eigen_eta Nonnegative reconstructed-dense filter controls.
#' @param nit,nburn,nthin MCMC iteration controls.
#' @param seed Fit-local base seed.
#' @param nchains Number of logical chains per trait.
#' @param ncores Requested concurrent trait-by-chain workers.
#' @param chain_seeds Optional signed chain seeds.
#' @param keep_chains Retain compact logical-chain records.
#' @param convergence Convergence mode: `"auto"`, `"none"`, `"core"`, or
#'   `"extended"`. Automatic mode remains core-only.
#' @param convergence_control A uniquely named convergence-control list. The
#'   extended controls cover diagnostic groups, explicit marker selection,
#'   retained traces, and hard trace-memory guards.
#' @param memory_warning_gb Analytical warning threshold.
#' @param verbose Print resolved execution information.
#' @param ... Model-specific validated controls.
#' @return A `stblr_fit` and `blr_fit` object.
#' @export
stblr_block_eigen <- function(
  stats, Glist, block_start,
  method = c("sbayesc", "sbayesr", "sbayesrc"),
  effect_maf = NULL, allow_reference_maf_for_maf_effect_s = FALSE,
  annotation = NULL,
  representation = c("low_rank", "dense_reconstructed"),
  eigen_policy = NULL, eigen_prop = 0.995, eigen_tau = 0.01, eigen_eta = 0,
  nit = 1000, nburn = 500, nthin = 1, seed = 1,
  nchains = 1L, ncores = 1L, chain_seeds = NULL,
  keep_chains = FALSE, convergence = c("auto", "none", "core", "extended"),
  convergence_control = NULL, memory_warning_gb = 8,
  verbose = FALSE, ...
) {
  dots <- list(...)
  legacy_filter <- dots$eigen_filter
  dots$eigen_filter <- NULL
  if (!is.null(legacy_filter) && missing(representation)) {
    representation <- "dense_reconstructed"
  }
  representation <- match.arg(representation)
  if (!is.null(legacy_filter)) {
    legacy_filter <- match.arg(
      legacy_filter, c("hard_truncate", "ridge_fixed", "ridge_lw"))
    if (identical(representation, "low_rank")) {
      stop("eigen_filter is only supported by representation = 'dense_reconstructed'.",
           call. = FALSE)
    }
    legacy_policy <- if (identical(legacy_filter, "hard_truncate"))
      "absolute_threshold" else legacy_filter
    if (!is.null(eigen_policy) && !identical(eigen_policy, legacy_policy)) {
      stop("eigen_filter and eigen_policy specify different dense policies.",
           call. = FALSE)
    }
    eigen_policy <- legacy_policy
  }
  if (is.null(eigen_policy)) {
    eigen_policy <- if (identical(representation, "low_rank"))
      "cumulative_positive_mass" else "absolute_threshold"
  }
  if (length(eigen_policy) != 1L || is.na(eigen_policy)) {
    stop("eigen_policy must be NULL or one non-missing string.", call. = FALSE)
  }
  supported <- if (identical(representation, "low_rank")) {
    "cumulative_positive_mass"
  } else {
    c("absolute_threshold", "ridge_fixed", "ridge_lw")
  }
  if (!eigen_policy %in% supported) {
    stop("Unsupported representation/eigen_policy combination.", call. = FALSE)
  }
  if (identical(representation, "low_rank") &&
      (length(eigen_prop) != 1L || !is.finite(eigen_prop) ||
       eigen_prop <= 0 || eigen_prop >= 1)) {
    stop("eigen_prop must be finite and strictly between 0 and 1.", call. = FALSE)
  }
  eigen_filter <- switch(eigen_policy,
    absolute_threshold = "hard_truncate",
    ridge_fixed = "ridge_fixed",
    ridge_lw = "ridge_lw",
    cumulative_positive_mass = "hard_truncate")
  resolved_model <- .blr_resolve_st_model(
    method, dots, c("sbayesc", "sbayesr", "sbayesrc"), "block_eigen")
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
    unlist(Glist$rsids %||% Glist$rsidsLD, use.names = FALSE)
  if (!length(marker_ids)) marker_ids <- paste0("V", seq_along(stats$ww[[1L]]))
  trace_spec <- .blr_st_native_trace_spec(
    conv, marker_ids, resolved_model$prior_kernel,
    annotations = identical(method, "sbayesrc"),
    component_count = if (resolved_model$prior_kernel %in% c("bayesr", "bayesrc"))
      length(dots$mixture_var %||% c(0, 0.01, 0.1, 1)) else 0L,
    annotation_quantity_count = if (identical(method, "sbayesrc")) {
      annotation_columns <- ncol(annotation %||% matrix(nrow = 0L, ncol = 0L))
      component_count <- length(dots$mixture_var %||% c(0, 0.01, 0.1, 1))
      annotation_columns * (component_count - 1L) + (component_count - 1L)
    } else 0L)
  memory <- .blr_st_preflight_memory(
    stats = stats, operator = "block_eigen", chain = chain, conv = conv,
    memory_warning_gb = memory_warning_gb, trace_spec = trace_spec)
  common <- list(
    stats = stats, Glist = Glist, block_start = block_start,
    representation = representation, eigen_policy = eigen_policy,
    eigen_prop = eigen_prop, eigen_filter = eigen_filter,
    eigen_tau = eigen_tau, eigen_eta = eigen_eta,
    nit = chain$nit, nburn = chain$nburn, nthin = chain$nthin,
    seed = chain$seed, nchains = chain$nchains, ncores = chain$ncores,
    chain_seeds = if (length(chain$chain_seeds_native))
      chain$chain_seeds_native else NULL,
    keep_chains = chain$keep_chains || conv$compute || conv$keep_traces,
    .convergence_spec = trace_spec)
  fit <- switch(
    method,
    sbayesc = do.call(.stblr_csr_bayesc_block_eigen, c(common, dots)),
    sbayesr = do.call(.stblr_csr_bayesr_block_eigen, c(common, dots)),
    sbayesrc = {
      if (is.null(annotation)) {
        stop("annotation is required for method = 'sbayesrc'.", call. = FALSE)
      }
      do.call(.stblr_csr_sbayesrc_block_eigen,
              c(common[c("stats", "Glist", "block_start", "representation",
                         "eigen_policy", "eigen_prop", "eigen_filter",
                         "eigen_tau", "eigen_eta")],
                list(annotation = annotation),
                common[c("nit", "nburn", "nthin", "seed", "nchains",
                         "ncores", "chain_seeds", "keep_chains",
                         ".convergence_spec")], dots))
    })
  fit <- .blr_finalize_st_public(
    fit, method, "block_eigen", chain, conv, memory_warning_gb, verbose,
    memory)
  fit$data$operator <- fit$input[c(
    "operator_representation", "operator_contract", "operator_scale_contract",
    "eigen_policy", "eigen_prop", "eigen_tau", "eigen_eta", "eigen_blocks")]
  fit$diagnostics$block_eigen <-
    fit$input$eigen_diagnostics %||% fit$diagnostics$block_eigen %||% NULL
  fit$input$effect_scale <- resolved_model$effect_scale
  fit$input$prior_kernel <- resolved_model$prior_kernel
  fit$input$probability_policy <- resolved_model$probability_policy
  fit$input$effect_maf_source <- maf_info$effect_maf_source
  fit$input$effect_maf_population <- maf_info$effect_maf_population
  fit$input$effect_maf_alignment_status <-
    maf_info$effect_maf_alignment_status
  fit$input$effect_maf_fallback_used <- maf_info$effect_maf_fallback_used
  fit$data$effect_maf_source <- maf_info$effect_maf_source
  fit$data$effect_maf_population <- maf_info$effect_maf_population
  fit$data$effect_maf_alignment_status <-
    maf_info$effect_maf_alignment_status
  fit$data$effect_maf_fallback_used <- maf_info$effect_maf_fallback_used
  fit
}
