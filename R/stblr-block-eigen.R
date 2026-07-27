#' Scalar-trait BLR with a block-eigen LD operator
#'
#' Fits one independent scalar model per trait using the canonical reconstructed
#' block-eigen operator. The operator is built from the supplied BED provenance;
#' filtering changes the represented LD operator and is recorded in the fit.
#'
#' @param stats Scalar-trait summary statistics.
#' @param Glist Genotype/BED provenance used to construct the operator.
#' @param block_start One-based public block starts.
#' @param method One of `"sbayesc"`, `"sbayesr"`, or `"sbayesrc"`;
#'   the `s` prefix denotes summary-statistics data.
#' @param selection_maf Optional marker-aligned allele frequencies used only
#'   when the independent `selection_s` variance-scaling option is active.
#' @param allow_reference_maf_for_selection_s Whether aligned reference-panel
#'   frequencies may be used explicitly when summary-population frequencies
#'   are unavailable. The default is `FALSE`.
#' @param annotation Annotation matrix required for `"sbayesrc"`.
#' @param eigen_filter Block filter: hard truncation, fixed ridge, or LW ridge.
#' @param eigen_tau,eigen_eta Nonnegative filter controls.
#' @param nit,nburn,nthin MCMC iteration controls.
#' @param seed Fit-local base seed.
#' @param nchains Number of logical chains per trait.
#' @param ncores Requested concurrent trait-by-chain workers.
#' @param chain_seeds Optional signed chain seeds.
#' @param keep_chains Retain compact logical-chain records.
#' @param convergence Convergence mode.
#' @param convergence_control Named convergence controls.
#' @param memory_warning_gb Analytical warning threshold.
#' @param verbose Print resolved execution information.
#' @param ... Model-specific validated controls.
#' @return A `stblr_fit` and `blr_fit` object.
#' @export
stblr_block_eigen <- function(
  stats, Glist, block_start,
  method = c("sbayesc", "sbayesr", "sbayesrc"),
  selection_maf = NULL, allow_reference_maf_for_selection_s = FALSE,
  annotation = NULL,
  eigen_filter = c("hard_truncate", "ridge_fixed", "ridge_lw"),
  eigen_tau = 0.01, eigen_eta = 0,
  nit = 1000, nburn = 500, nthin = 1, seed = 1,
  nchains = 1L, ncores = 1L, chain_seeds = NULL,
  keep_chains = FALSE, convergence = c("auto", "none", "core"),
  convergence_control = NULL, memory_warning_gb = 8,
  verbose = FALSE, ...
) {
  dots <- list(...)
  resolved_model <- .blr_resolve_st_model(
    method, dots, c("sbayesc", "sbayesr", "sbayesrc"), "block_eigen")
  method <- resolved_model$model
  dots <- resolved_model$dots
  maf_info <- .blr_resolve_st_selection_maf(
    selection_maf, allow_reference_maf_for_selection_s,
    resolved_model$selection_s_active, stats, Glist)
  Glist <- maf_info$Glist
  eigen_filter <- match.arg(eigen_filter)
  chain <- .blr_chain_controls(
    nit, nburn, nthin, seed, nchains, ncores, chain_seeds, keep_chains)
  conv <- .blr_convergence_controls(convergence, convergence_control,
                                    chain$nchains)
  memory <- .blr_st_preflight_memory(
    stats = stats, operator = "block_eigen", chain = chain, conv = conv,
    memory_warning_gb = memory_warning_gb)
  common <- list(
    stats = stats, Glist = Glist, block_start = block_start,
    eigen_filter = eigen_filter, eigen_tau = eigen_tau, eigen_eta = eigen_eta,
    nit = chain$nit, nburn = chain$nburn, nthin = chain$nthin,
    seed = chain$seed, nchains = chain$nchains, ncores = chain$ncores,
    chain_seeds = if (length(chain$chain_seeds_native))
      chain$chain_seeds_native else NULL,
    keep_chains = chain$keep_chains || conv$compute || conv$keep_traces)
  fit <- switch(
    method,
    sbayesc = do.call(.stblr_csr_bayesc_block_eigen, c(common, dots)),
    sbayesr = do.call(.stblr_csr_bayesr_block_eigen, c(common, dots)),
    sbayesrc = {
      if (is.null(annotation)) {
        stop("annotation is required for method = 'sbayesrc'.", call. = FALSE)
      }
      do.call(.stblr_csr_sbayesrc_block_eigen,
              c(common[c("stats", "Glist", "block_start", "eigen_filter",
                         "eigen_tau", "eigen_eta")],
                list(annotation = annotation),
                common[c("nit", "nburn", "nthin", "seed", "nchains",
                         "ncores", "chain_seeds", "keep_chains")], dots))
    })
  fit <- .blr_finalize_st_public(
    fit, method, "block_eigen", chain, conv, memory_warning_gb, verbose,
    memory)
  fit$data$operator <- fit$input[c(
    "eigen_filter", "eigen_tau", "eigen_eta", "eigen_blocks")]
  fit$diagnostics$block_eigen <-
    fit$input$eigen_diagnostics %||% fit$diagnostics$block_eigen %||% NULL
  fit$input$effect_scale <- resolved_model$effect_scale
  fit$input$prior_kernel <- resolved_model$prior_kernel
  fit$input$probability_policy <- resolved_model$probability_policy
  fit$input$selection_maf_source <- maf_info$selection_maf_source
  fit$input$selection_maf_population <- maf_info$selection_maf_population
  fit$input$selection_maf_alignment_status <-
    maf_info$selection_maf_alignment_status
  fit$input$selection_maf_fallback_used <- maf_info$selection_maf_fallback_used
  fit$data$selection_maf_source <- maf_info$selection_maf_source
  fit$data$selection_maf_population <- maf_info$selection_maf_population
  fit$data$selection_maf_alignment_status <-
    maf_info$selection_maf_alignment_status
  fit$data$selection_maf_fallback_used <- maf_info$selection_maf_fallback_used
  fit
}
