#' Scalar-trait BLR with sparse-LD CSR operators
#'
#' @param stats Scalar-trait summary statistics.
#' @param Glist Optional genotype/provenance list.
#' @param ld_prefix Optional CSR prefix.
#' @param method `"bayesc"`, `"sbayesc"`, `"bayesr"`, or `"sbayesr"`.
#' @param nit,nburn,nthin MCMC controls.
#' @param seed Fit-local base seed.
#' @param nchains Number of logical chains per trait.
#' @param ncores Requested concurrent trait-by-chain workers.
#' @param chain_seeds Optional signed chain seeds.
#' @param keep_chains Retain compact logical-chain records.
#' @param convergence Convergence mode.
#' @param convergence_control Named convergence controls.
#' @param memory_warning_gb Analytical warning threshold.
#' @param verbose Print resolved controls.
#' @param ... Model-specific controls accepted by the canonical implementation.
#' @return A `stblr_fit` and `blr_fit` object.
#' @export
stblr_csr <- function(
  stats, Glist = NULL, ld_prefix = NULL,
  method = c("bayesc", "sbayesc", "bayesr", "sbayesr"), ...,
  nit = 1000, nburn = 500, nthin = 1, seed = 1,
  nchains = 1L, ncores = 1L, chain_seeds = NULL,
  keep_chains = FALSE, convergence = c("auto", "none", "core"),
  convergence_control = NULL, memory_warning_gb = 8, verbose = FALSE
) {
  dots <- list(...)
  resolved_model <- .blr_resolve_st_model(
    method, dots, c("bayesc", "sbayesc", "bayesr", "sbayesr"))
  method <- resolved_model$model
  dots <- resolved_model$dots
  chain <- .blr_chain_controls(
    nit, nburn, nthin, seed, nchains, ncores, chain_seeds, keep_chains)
  conv <- .blr_convergence_controls(convergence, convergence_control,
                                    chain$nchains)
  memory <- .blr_st_preflight_memory(
    stats = stats, operator = "csr", chain = chain, conv = conv,
    memory_warning_gb = memory_warning_gb)
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
    keep_chains = chain$keep_chains || conv$compute || conv$keep_traces), dots))
  fit <- .blr_finalize_st_public(
    fit, method, "csr", chain, conv, memory_warning_gb, verbose, memory)
  fit$input$effect_scale <- resolved_model$effect_scale
  fit$input$probability_policy <- resolved_model$probability_policy
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
  keep_chains = FALSE, convergence = c("auto", "none", "core"),
  convergence_control = NULL, memory_warning_gb = 8, verbose = FALSE
) {
  method <- match.arg(method)
  chain <- .blr_chain_controls(
    nit, nburn, nthin, seed, nchains, ncores, chain_seeds, keep_chains)
  conv <- .blr_convergence_controls(convergence, convergence_control,
                                    chain$nchains)
  memory <- .blr_st_preflight_memory(
    y = y, Glist = Glist, operator = "packed_bed", chain = chain,
    conv = conv, memory_warning_gb = memory_warning_gb)
  dots <- list(...)
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
    keep_chains = chain$keep_chains || conv$compute || conv$keep_traces), dots))
  fit <- .blr_finalize_st_public(
    fit, method, "packed_bed", chain, conv, memory_warning_gb, verbose,
    memory)
  fit$input$effect_scale <- if (method %in% c("bayesr", "bayesrc"))
    "component" else "unit"
  fit$input$probability_policy <- if (method == "bayesrc")
    "annotation_probit_stick" else "global"
  fit$data$effect_scale <- fit$input$effect_scale
  fit
}
