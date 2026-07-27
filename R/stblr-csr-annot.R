#' Fit Annotation-Aware CSR ST-BLR Models
#'
#' `stblr_csr_annot()` is the unified public entry point for the current
#' annotation-aware CSR summary-statistics backends. It dispatches to the
#' internal model adapters and standardizes fit metadata and annotation output.
#'
#' Available annotation models are `annotation_model = "prior"` for fixed
#' annotation-informed prior probabilities and variance multipliers,
#' `annotation_model = "learned"` for learned annotation effects on BayesC-like
#' inclusion and variance priors, `annotation_model = "group"` for grouped
#' annotation architectures, and `annotation_model = "sbayesrc"` for
#' SBayesRC-style annotation-dependent component probabilities. Only
#' `annotation_model = "sbayesrc"` supports fixed global `selection_s` and
#' sampled trait-specific `selection_s`.
#'
#' @param stats Sufficient statistics returned by [bed_xtx_xty()].
#' @param Glist Optional object containing sparse-LD metadata. When `ld_prefix`
#'   is not supplied, `Glist$sparseLD$prefix` is used when available.
#' @param annotations Annotation input. For `annotation_model = "prior"`, use a
#'   marker x annotation matrix or a list with optional `A`, `fixed_pi_marker`
#'   or `pi_marker`, and `fixed_vb_multiplier` or `vb_multiplier` elements. For
#'   `"learned"` and `"sbayesrc"`, use a marker x annotation matrix. For
#'   `"group"`, use a length-`m` group vector, factor, or one-column data frame.
#' @param annotation_model Annotation model. Synonyms are accepted:
#'   `"prior"`, `"fixed_prior"`, `"fixed"`; `"learned"`, `"annot"`,
#'   `"annotation"`; `"group"`, `"groups"`; and `"sbayesrc"`/`"SBayesRC"`.
#' @param method Optional lowercase method override. BayesC-like annotation
#'   models use `"bayesc"`; SBayesRC uses `"sbayesrc"`.
#' @param nit,nburn,nthin MCMC iteration controls.
#' @param ncores Number of OpenMP threads.
#' @param seed Sampler seed.
#' @param nchains Number of independent chains.
#' @param chain_seeds Optional integer seeds, one per chain.
#' @param keep_chains Logical; return compact per-chain summaries when
#'   supported.
#' @param convergence Convergence mode: `"auto"`, `"none"`, or `"core"`.
#' @param convergence_control Optional named convergence-control list.
#' @param memory_warning_gb Analytical memory-warning threshold in GiB.
#' @param verbose Print resolved execution metadata.
#' @param updateLDswap Logical; request optional LD-swap/MH. This is supported
#'   for the current annotation-aware CSR models, including SBayesRC.
#' @param selection_s Optional fixed global BayesS-style MAF-scaling parameter.
#'   Currently supported only for `annotation_model = "sbayesrc"` in this
#'   unified annotation interface. The default `selection_s = NULL` with
#'   `estimate_selection_s = FALSE` fits the ordinary model. A finite numeric
#'   scalar with `estimate_selection_s = FALSE` fits a fixed global-S model.
#'   `selection_s` must remain `NULL` when `estimate_selection_s = TRUE`;
#'   fixed and sampled S cannot both be requested.
#' @param estimate_selection_s Logical; estimate one trait-specific
#'   BayesS-style `selection_s` by Metropolis-Hastings. Currently supported
#'   only for `annotation_model = "sbayesrc"`. Sampled `selection_s` is not
#'   supported for `annotation_model = "prior"`, `"learned"`, or `"group"`.
#' @param selection_s_init Initial value for sampled `selection_s`. Defaults to
#'   0 and is used only when `estimate_selection_s = TRUE`.
#' @param selection_s_prior Numeric length-2 lower and upper bounds for the
#'   uniform sampled-`selection_s` prior. Defaults to `c(-3, 2)` and is used
#'   only when `estimate_selection_s = TRUE`.
#' @param selection_s_proposal_sd Random-walk proposal standard deviation for
#'   sampled `selection_s`. Defaults to 0.35 and is used only when
#'   `estimate_selection_s = TRUE`.
#' @param ld_swap_prob Probability per MCMC iteration of attempting LD-swap
#'   moves when `updateLDswap = TRUE`.
#' @param ld_swap_r2 Minimum LD r-squared for candidate swap partners.
#' @param ld_swap_max_friends Maximum number of high-LD friends stored per
#'   marker for swap proposals.
#' @param ld_swap_moves Number of swap attempts when LD-swap is triggered.
#' @param h2 Initial heritability.
#' @param adjE Residual adjustment factor.
#' @param updateB,updateE,updatePi Logical sampler update controls. `updatePi`
#'   is used by the BayesC-like annotation models and ignored by SBayesRC.
#' @param ld_prefix Optional prefix of the disk-backed CSR LD files.
#' @param ... Additional model-specific arguments passed to the existing
#'   annotation wrapper.
#'
#' @return A formatted ST-BLR fit. All returned fits contain `dm`
#'   (marker-by-trait posterior inclusion or non-null probability), `bm`
#'   (marker-by-trait posterior mean effects), `vbs`, `vgs`, `ves`, `vle`,
#'   `vld`, and `input`, plus standardized metadata fields (`method`, `model`,
#'   `backend`, `data_level`, `annotation_model`, `annotations`, `nchains`,
#'   `keep_chains`). Existing wrapper-specific fields are preserved. The `vle`
#'   and `vld` traces follow the same definitions and formatting conventions as
#'   annotation-unaware CSR fits.
#'
#'   Annotation-aware outputs are model-specific and may include
#'   `annotation_summary`, `annotation_pi`, `annotation_effects`,
#'   `annotation_prior`, `alpha`, and `sigmaSqAlpha`. SBayesRC fits also
#'   include BayesR-style `component_probabilities` marker-by-component posterior
#'   probabilities by trait, `dm_component_mean`, and `ncomp` where available.
#'   The SBayesRC null component column is `gamma_0.00`, so
#'   `dm = 1 - P(gamma_0.00)`.
#'
#'   LD-swap-enabled fits include `fit$diagnostics$ld_swap` and, where chain
#'   summaries are retained, `fit$diagnostics$ld_swap_chains` or chain-level
#'   LD-swap entries. With
#'   `keep_chains = TRUE`, compact per-chain summaries are returned in
#'   `chains` for supported CSR annotation backends.
#'
#'   For sampled SBayesRC `selection_s`, the fit also contains `selection_s`,
#'   `selection_s_chain_mean_sd`, `selection_s_chain_mean_min`,
#'   `selection_s_chain_mean_max`,
#'   `selection_s_trace`, and `selection_s_acceptance`. `selection_s_trace` is
#'   an iteration x trait matrix, `selection_s` is the posterior mean by trait,
#'   and `selection_s_acceptance` is the MH acceptance rate by trait. With
#'   `keep_chains = TRUE`, chain-level sampled-S output is available as
#'   the matching flat trait-by-chain record's `selection_s` and
#'   `selection_s_acceptance` fields.
#'
#'   Fine-mapping diagnostics are available through PIP summaries in `dm` and
#'   LD-swap output when `updateLDswap = TRUE`. Credible-set construction is
#'   performed by helper functions such as [make_credible_sets()] and
#'   [extract_stblr_finemap_loci()] from posterior inclusion probabilities and
#'   LD, rather than being a separate sampler return object.
#'
#' @details
#' This interface is CSR summary-statistics only. The current annotation-aware
#' SBayesRC and the BayesC-like annotation-aware CSR backends support multiple
#' native chains and compact chain summaries. LD-swap/MH is currently available
#' for the fixed-prior, group, and learned annotation-aware CSR BayesC
#' backends and the annotation-aware CSR SBayesRC backend.
#'
#' Annotation-aware BayesC models differ from SBayesRC in what annotations
#' control. The `"prior"`, `"learned"`, and `"group"` models affect BayesC-like
#' inclusion probabilities and/or marker-effect variance priors directly.
#' SBayesRC uses annotations to affect component probabilities; `selection_s`,
#' when requested, affects only marker-specific effect-size prior variance.
#'
#' The model-specific adapters are internal. `stblr_csr_annot()` is the sole
#' public annotation-aware CSR fitting entry.
#'
#' CSR effects are on the standardized-genotype scale. The BayesS-style
#' MAF-dependent prior scale used by fixed and sampled `selection_s` is
#' `q_j(S) = h_j^(S + 1)`, where `h_j = 2 p_j (1 - p_j)` and `p_j` is the
#' minor allele frequency. The `+1` exponent appears because the sampler effects
#' are standardized-genotype-scale effects rather than allele-scale effects.
#'
#' For CSR SBayesRC, annotations affect component probabilities and
#' `selection_s` affects marker-specific effect-size prior variance. Fixed
#' `selection_s` uses
#' `b_j | component_j = m, vb, S ~ N(0, vb * gamma_m * q_j(S))`. The null
#' component column is `gamma_0.00` and `dm = 1 - P(gamma_0.00)`.
#'
#' Sampled CSR SBayesRC uses the trait-specific MH log posterior contribution
#'
#' ```text
#' log p(S_t | b_t, gamma_t, vb_t)
#'   = log p(S_t)
#'     - 0.5 sum_\{j: gamma_jt > 0\} [
#'         log q_j(S_t) + b_jt^2 / (vb_t gamma_jt q_j(S_t))
#'       ]
#' ```
#'
#' `S_t` is sampled separately for each trait and chain. Posterior summaries
#' are averaged or summarized across chains in the returned fit object.
#'
#' @examples
#' \dontrun{
#' fit_sampled_s_sbayesrc <- stblr_csr_annot(
#'   stats = stats,
#'   Glist = Glist,
#'   annotations = annotations,
#'   annotation_model = "sbayesrc",
#'   estimate_selection_s = TRUE,
#'   selection_s_prior = c(-3, 2),
#'   selection_s_proposal_sd = 0.35
#' )
#' }
#'
#' @export
stblr_csr_annot <- function(
  stats,
  Glist = NULL,
  annotations,
  annotation_model = c("prior", "learned", "group", "sbayesrc"),
  method = NULL,
  nit = 1000,
  nburn = 100,
  nthin = 1,
  seed = 1,
  nchains = 1L,
  ncores = 1L,
  chain_seeds = NULL,
  keep_chains = FALSE,
  convergence = c("auto", "none", "core"),
  convergence_control = NULL,
  memory_warning_gb = 8,
  verbose = FALSE,
  h2 = 0.3,
  adjE = 0.9,
  updateB = TRUE,
  updateE = TRUE,
  updatePi = TRUE,
  updateLDswap = FALSE,
  selection_s = NULL,
  estimate_selection_s = FALSE,
  selection_s_init = 0,
  selection_s_prior = c(-3, 2),
  selection_s_proposal_sd = 0.35,
  ld_swap_prob = 0.05,
  ld_swap_r2 = 0.8,
  ld_swap_max_friends = 50L,
  ld_swap_moves = 1L,
  ld_prefix = NULL,
  ...
) {
 if (missing(annotations)) {
  stop("annotations must be supplied.")
 }

 annotation_model <- .stblr_match_annotation_model(annotation_model[[1L]])
 chain <- .blr_chain_controls(
  nit, nburn, nthin, seed, nchains, ncores, chain_seeds, keep_chains)
 conv <- .blr_convergence_controls(
  convergence, convergence_control, chain$nchains)
 memory <- .blr_st_preflight_memory(
  stats = stats, operator = "csr", chain = chain, conv = conv,
  memory_warning_gb = memory_warning_gb)
 .stblr_check_annotation_method(method, annotation_model)
 .stblr_check_annotation_chains(annotation_model, nchains, keep_chains, chain_seeds)
 .validate_ld_swap_args(
  updateLDswap, ld_swap_prob, ld_swap_r2, ld_swap_max_friends, ld_swap_moves
 )
 .stblr_validate_sampled_selection_s(
  estimate_selection_s = estimate_selection_s,
  selection_s = selection_s,
  selection_s_init = selection_s_init,
  selection_s_prior = selection_s_prior,
  selection_s_proposal_sd = selection_s_proposal_sd
 )
 ld_prefix <- .stblr_resolve_csr_annotation_ld_prefix(
  Glist = Glist,
  ld_prefix = ld_prefix
 )

 common <- list(
  stats = stats,
  ld_prefix = ld_prefix,
  h2 = h2,
  adjE = adjE,
  updateB = updateB,
  updateE = updateE,
  nit = chain$nit,
  nburn = chain$nburn,
  nthin = chain$nthin,
  ncores = chain$ncores,
  seed = chain$seed,
  nchains = chain$nchains,
  chain_seeds = if (length(chain$chain_seeds_native))
    chain$chain_seeds_native else NULL,
  keep_chains = chain$keep_chains || conv$compute || conv$keep_traces
 )
 probability_policy <- switch(
  annotation_model, prior = "fixed_marker", learned = "learned_logistic",
  group = "group", sbayesrc = "annotation_probit_stick")
 finish <- function(fit) {
  model <- if (annotation_model == "sbayesrc") "sbayesrc" else "bayesc"
  out <- .blr_finalize_st_public(
   fit, model, "csr", chain, conv, memory_warning_gb, verbose, memory)
  out$input$annotation_policy <- probability_policy
  out$input$probability_policy <- probability_policy
  out$input$effect_scale <- if (model == "sbayesrc") "component_maf_s" else "unit"
  out
 }
 extra <- list(...)

 common$updateLDswap <- updateLDswap
 common$ld_swap_prob <- ld_swap_prob
 common$ld_swap_r2 <- ld_swap_r2
 common$ld_swap_max_friends <- ld_swap_max_friends
 common$ld_swap_moves <- ld_swap_moves

 if (annotation_model %in% c("prior", "learned", "group")) {
  if (isTRUE(estimate_selection_s)) {
   stop(
    "estimate_selection_s is currently supported only for annotation_model = \"sbayesrc\".",
    call. = FALSE
   )
  }
  if (!is.null(selection_s)) {
   stop(
    "selection_s is currently supported only for annotation_model = \"sbayesrc\".",
    call. = FALSE
   )
  }
  common$updatePi <- updatePi
 }

 if (annotation_model == "prior") {
  args <- .stblr_prior_annotation_args(annotations, extra)
  args <- c(common, args)
  return(finish(do.call(stblr_csr_prior_annot, args)))
 }

 if (annotation_model == "learned") {
  if ("A" %in% names(extra)) {
   stop("Supply learned annotations through annotations, not both annotations and A.")
  }
  args <- c(common, list(A = annotations), extra)
  return(finish(do.call(stblr_csr_learn_annot, args)))
 }

 if (annotation_model == "group") {
  if ("group" %in% names(extra)) {
   stop("Supply group annotations through annotations, not both annotations and group.")
  }
  args <- c(common, list(group = annotations), extra)
  return(finish(do.call(stblr_csr_group_annot, args)))
 }

 if ("A" %in% names(extra)) {
  stop("Supply SBayesRC annotations through annotations, not both annotations and A.")
 }
 if (is.null(selection_s) && !isTRUE(estimate_selection_s)) selection_s <- 0
 args <- c(
  common,
  list(
   Glist = Glist,
   A = annotations,
   selection_s = selection_s,
   estimate_selection_s = estimate_selection_s,
   selection_s_init = selection_s_init,
   selection_s_prior = selection_s_prior,
   selection_s_proposal_sd = selection_s_proposal_sd
  ),
  extra
 )
 finish(do.call(stblr_csr_sbayesrc_generic, args))
}

.stblr_check_annotation_chains <- function(annotation_model, nchains, keep_chains,
                                          chain_seeds = NULL) {
 if (!is.numeric(nchains) || length(nchains) != 1L ||
     !is.finite(nchains) || nchains < 1 || nchains != floor(nchains)) {
  stop("nchains must be a positive integer scalar.", call. = FALSE)
 }
 if (!is.logical(keep_chains) || length(keep_chains) != 1L ||
     is.na(keep_chains)) {
  stop("keep_chains must be TRUE or FALSE.", call. = FALSE)
 }
 if (!is.null(chain_seeds) &&
     (!is.numeric(chain_seeds) || length(chain_seeds) != as.integer(nchains) ||
      anyNA(chain_seeds) || any(!is.finite(chain_seeds)) ||
      any(chain_seeds != floor(chain_seeds)))) {
  stop("chain_seeds must be NULL or an integer/numeric vector of length nchains.",
       call. = FALSE)
 }
 invisible(TRUE)
}

.stblr_check_annotation_method <- function(method, annotation_model) {
 if (is.null(method)) return(invisible(TRUE))
 method <- tolower(gsub("-", "", method))

 if (annotation_model %in% c("prior", "learned", "group")) {
  if (!method %in% c("bayesc", "bayescsr")) {
   stop(
    "BayesC-like annotation models require method = NULL or \"bayesc\".",
    call. = FALSE
   )
  }
  return(invisible(TRUE))
 }

 if (!method %in% "sbayesrc") {
  stop(
   "annotation_model = \"sbayesrc\" requires method = NULL or \"sbayesrc\".",
   call. = FALSE
  )
 }
 invisible(TRUE)
}

.stblr_resolve_csr_annotation_ld_prefix <- function(Glist = NULL, ld_prefix = NULL) {
 if (!is.null(ld_prefix)) return(ld_prefix)

 candidates <- list(
  if (!is.null(Glist)) Glist$sparseLD$prefix else NULL,
  if (!is.null(Glist)) Glist$ld_prefix else NULL
 )
 for (candidate in candidates) {
  if (!is.null(candidate) && length(candidate) == 1L && nzchar(candidate)) {
   return(candidate)
  }
 }

 stop(
  "ld_prefix must be supplied, or Glist must contain sparseLD$prefix.",
  call. = FALSE
 )
}

.stblr_prior_annotation_args <- function(annotations, extra) {
 args <- list()

 if (is.list(annotations) && !is.data.frame(annotations)) {
  known <- c(
   "A", "annotations", "fixed_pi_marker", "pi_marker",
   "fixed_vb_multiplier", "vb_multiplier", "beta_pi", "beta_vb",
   "use_pi_marker", "use_vb_multiplier"
  )
  unknown <- setdiff(names(annotations), known)
  if (length(unknown) > 0L) {
   stop(
    "Unknown prior annotation field(s): ",
    paste(unknown, collapse = ", "),
    call. = FALSE
   )
  }

  A <- annotations$A %||% annotations$annotations
  if (!is.null(A) && !"A" %in% names(extra)) args$A <- A

  pi_marker <- annotations$fixed_pi_marker %||% annotations$pi_marker
  if (!is.null(pi_marker) && !"fixed_pi_marker" %in% names(extra)) {
   args$fixed_pi_marker <- pi_marker
   if (!"use_pi_marker" %in% names(extra)) args$use_pi_marker <- TRUE
  }

  vb_multiplier <- annotations$fixed_vb_multiplier %||% annotations$vb_multiplier
  if (!is.null(vb_multiplier) && !"fixed_vb_multiplier" %in% names(extra)) {
   args$fixed_vb_multiplier <- vb_multiplier
   if (!"use_vb_multiplier" %in% names(extra)) args$use_vb_multiplier <- TRUE
  }

  if (!is.null(annotations$beta_pi) && !"beta_pi" %in% names(extra)) {
   args$beta_pi <- annotations$beta_pi
  }
  if (!is.null(annotations$beta_vb) && !"beta_vb" %in% names(extra)) {
   args$beta_vb <- annotations$beta_vb
  }
  if (!is.null(annotations$use_pi_marker) && !"use_pi_marker" %in% names(extra)) {
   args$use_pi_marker <- annotations$use_pi_marker
  }
  if (!is.null(annotations$use_vb_multiplier) &&
      !"use_vb_multiplier" %in% names(extra)) {
   args$use_vb_multiplier <- annotations$use_vb_multiplier
  }

  return(c(args, extra))
 }

 if ("A" %in% names(extra)) {
  stop("Supply prior annotation matrix through annotations, not both annotations and A.")
 }
 c(list(A = annotations), extra)
}

.standardize_stblr_annotation_fit <- function(fit, annotation_model) {
 annotation_model <- .stblr_match_annotation_model(annotation_model)
 if (is.null(fit$input)) fit$input <- list()

 meta <- switch(
  annotation_model,
  prior = list(method = "bayesc", model = "bayesc", backend = "csr_prior_bayesc"),
  learned = list(method = "bayesc", model = "bayesc", backend = "csr_annot_bayesc"),
  group = list(method = "bayesc", model = "bayesc", backend = "csr_group_bayesc"),
  sbayesrc = list(method = "sbayesrc", model = "sbayesrc", backend = "csr_sbayesrc")
 )

 fit$input$method <- meta$method
 fit$input$model <- meta$model
 fit$input$backend <- meta$backend
 fit$input$data_level <- "summary"
 fit$input$annotation_model <- annotation_model
 fit$input$annotations <- TRUE
 if (is.null(fit$input$nchains)) fit$input$nchains <- 1L
 if (is.null(fit$input$keep_chains)) fit$input$keep_chains <- FALSE

 if (annotation_model == "prior") {
  fit <- .stblr_add_prior_annotation_aliases(fit)
 } else if (annotation_model == "learned") {
  fit <- .stblr_add_learned_annotation_aliases(fit)
 } else if (annotation_model == "group") {
  fit <- .stblr_add_group_annotation_aliases(fit)
 } else {
  fit <- .stblr_add_sbayesrc_annotation_aliases(fit)
 }

 fit
}

.stblr_add_prior_annotation_aliases <- function(fit) {
 if (!is.null(fit$input$A)) fit$annotation <- fit$input$A

 fit$annotation_prior <- list(
  pi_marker = fit$input$pi_marker,
  vb_multiplier = fit$input$vb_multiplier,
  use_pi_marker = isTRUE(fit$input$use_pi_marker),
  use_vb_multiplier = isTRUE(fit$input$use_vb_multiplier)
 )

 pieces <- list()
 if (!is.null(fit$input$pi_marker)) {
  pieces$pi_marker <- .stblr_summarize_trait_vectors(
   fit$input$pi_marker,
   value = "pi_marker"
  )
 }
 if (!is.null(fit$input$vb_multiplier)) {
  pieces$vb_multiplier <- .stblr_summarize_trait_vectors(
   fit$input$vb_multiplier,
   value = "vb_multiplier"
  )
 }
 if (length(pieces) > 0L) {
  fit$annotation_summary <- do.call(rbind, pieces)
  rownames(fit$annotation_summary) <- NULL
 }

 fit
}

.stblr_add_learned_annotation_aliases <- function(fit) {
 if (!is.null(fit$input$A)) fit$annotation <- fit$input$A
 fit$annotation_effects <- list(pi = fit$eta_pi, variance = fit$eta_vb)

 if (is.matrix(fit$eta_pi) && is.matrix(fit$eta_vb)) {
  fit$annotation_summary <- .stblr_learned_annotation_summary(
   eta_pi = fit$eta_pi,
   eta_vb = fit$eta_vb
  )
 }

 fit
}

.stblr_add_group_annotation_aliases <- function(fit) {
 fit$annotation <- list(
  group = fit$input$group,
  group_names = fit$input$group_names,
  group_size = fit$input$group_size
 )
 fit$annotation_pi <- fit$group_pi
 fit$annotation_variance <- fit$group_vb_multiplier

 if (is.matrix(fit$group_pi) && is.matrix(fit$group_vb_multiplier)) {
  fit$annotation_summary <- .stblr_group_annotation_summary(fit)
 }

 fit
}

.stblr_add_sbayesrc_annotation_aliases <- function(fit) {
 if (!is.null(fit$input$A)) fit$annotation <- fit$input$A
 fit$annotation_effects <- fit$alpha
 fit$annotation_variance <- fit$sigmaSqAlpha

 if (!is.null(fit$alpha) && !is.null(fit$input$gamma)) {
  fit$annotation_pi <- lapply(
   fit$alpha,
   sbayesrc_annotation_pi,
   gamma = fit$input$gamma
  )
  fit$annotation_summary <- .stblr_sbayesrc_annotation_summary(fit)
 }

 fit
}

.stblr_summarize_trait_vectors <- function(x, value) {
 if (!is.list(x)) x <- list(x)
 trait_names <- names(x)
 if (is.null(trait_names)) trait_names <- paste0("T", seq_along(x))

 out <- do.call(rbind, lapply(seq_along(x), function(i) {
  v <- as.numeric(x[[i]])
  data.frame(
   trait = trait_names[[i]],
   value = value,
   mean = mean(v, na.rm = TRUE),
   min = min(v, na.rm = TRUE),
   max = max(v, na.rm = TRUE),
   stringsAsFactors = FALSE
  )
 }))
 out
}

.stblr_learned_annotation_summary <- function(eta_pi, eta_vb) {
 out <- list()
 k <- 1L
 for (trait in rownames(eta_pi)) {
  for (annotation in colnames(eta_pi)) {
   out[[k]] <- data.frame(
    trait = trait,
    annotation = annotation,
    eta_pi = eta_pi[trait, annotation],
    eta_vb = eta_vb[trait, annotation],
    stringsAsFactors = FALSE
   )
   k <- k + 1L
  }
 }
 do.call(rbind, out)
}

.stblr_group_annotation_summary <- function(fit) {
 out <- list()
 k <- 1L
 for (trait in rownames(fit$group_pi)) {
  for (group in colnames(fit$group_pi)) {
   out[[k]] <- data.frame(
    trait = trait,
    group = group,
    group_pi = fit$group_pi[trait, group],
    group_vb_multiplier = fit$group_vb_multiplier[trait, group],
    group_nincluded = fit$group_nincluded[trait, group],
    group_size = fit$group_size[trait, group],
    stringsAsFactors = FALSE
   )
   k <- k + 1L
  }
 }
 do.call(rbind, out)
}

.stblr_sbayesrc_annotation_summary <- function(fit) {
 out <- list()
 k <- 1L
 for (trait in names(fit$alpha)) {
  gamma_mean <- sbayesrc_annotation_gamma_mean(
   fit$alpha[[trait]],
   gamma = fit$input$gamma
  )
  for (annotation in names(gamma_mean)) {
   out[[k]] <- data.frame(
    trait = trait,
    annotation = annotation,
    expected_gamma = gamma_mean[[annotation]],
    stringsAsFactors = FALSE
   )
   k <- k + 1L
  }
 }
 do.call(rbind, out)
}
