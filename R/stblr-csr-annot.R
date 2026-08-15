.stblr_validate_learned_forwarded_args <- function(extra) {
 if (!length(extra)) return(invisible(TRUE))

 extra_names <- names(extra)
 if (is.null(extra_names)) extra_names <- rep("", length(extra))
 invalid_names <- is.na(extra_names) | !nzchar(extra_names) | extra_names == "NA"
 if (any(invalid_names)) {
  stop(
   paste0(
    "Arguments forwarded to annotation_model = \"learned_logistic\" through ",
    "`...` must have unique, nonempty names; unnamed, empty-name, and ",
    "NA-name arguments are not accepted."
   ),
   call. = FALSE
  )
 }

 duplicate_names <- unique(extra_names[duplicated(extra_names)])
 if (length(duplicate_names)) {
  stop(
   paste0(
    "Arguments forwarded to annotation_model = \"learned_logistic\" through ",
    "`...` must have unique, nonempty names; duplicated argument(s): ",
    paste(sprintf("`%s`", duplicate_names), collapse = ", "), "."
   ),
   call. = FALSE
  )
 }

 learned_probability_controls <- c("pi_min", "pi_max", "pi_marker")
 ambiguous_probability_controls <- vapply(
  extra_names,
  function(arg) {
   if (identical(arg, "pi_marker")) return(FALSE)
   any(startsWith(learned_probability_controls, arg))
  },
  logical(1)
 )
 if (any(ambiguous_probability_controls)) {
  offending <- unique(extra_names[ambiguous_probability_controls])
  stop(
   paste0(
    "Argument(s) ",
    paste(sprintf("`%s`", offending), collapse = ", "),
    " are not accepted for annotation_model = \"learned_logistic\". ",
    "Use the exact name `pi_marker` for the maintained initial global ",
    "probability. Learned annotation probabilities are not controlled by ",
    "`pi_min` or `pi_max`."
   ),
   call. = FALSE
  )
 }

 invisible(TRUE)
}

#' Fit Annotation-Aware CSR ST-BLR Models
#'
#' `stblr_csr_annot()` is the unified public entry point for the current
#' annotation-aware CSR summary-statistics backends. It dispatches to the
#' internal model adapters and standardizes fit metadata and annotation output.
#'
#' Available annotation policies are `annotation_model = "fixed_marker"` for fixed
#' annotation-informed prior probabilities and variance multipliers,
#' `annotation_model = "learned_logistic"` for learned annotation effects on BayesC-like
#' inclusion and variance priors, `annotation_model = "group"` for grouped
#' annotation architectures, `annotation_model = "annotation_probit_stick"` for
#' SBayesRC-style annotation-dependent component probabilities, and
#' `annotation_model = "log_variance"` for BayesC-LV or BayesR-LV. In the LV
#' models annotations rescale non-null prior variances while inclusion or
#' component probabilities remain global. Only
#' `annotation_model = "annotation_probit_stick"` supports fixed global `maf_effect_s` and
#' sampled trait-specific `maf_effect_s`.
#'
#' @param stats Sufficient statistics returned by [bed_xtx_xty()].
#' @param Glist Optional object containing sparse-LD metadata. When `ld_prefix`
#'   is not supplied, `Glist$sparseLD$prefix` is used when available.
#' @param annotations Annotation input. For `annotation_model = "fixed_marker"`, use a
#'   marker x annotation matrix or a list with optional `A`, `fixed_pi_marker`
#'   or `pi_marker`, and `fixed_vb_multiplier` or `vb_multiplier` elements. For
#'   `"learned_logistic"` and `"annotation_probit_stick"`, use a marker x
#'   annotation matrix. For `"log_variance"`, binary columns are centered only
#'   and continuous columns are centered and standardized to SD one; intercept,
#'   constant, non-finite, duplicate, and rank-deficient designs are rejected. For
#'   `"group"`, use a length-`m` group vector, factor, or one-column data frame.
#' @param annotation_model Annotation policy. Must be one of `"fixed_marker"`,
#'   `"group"`, `"learned_logistic"`, `"annotation_probit_stick"`, or
#'   `"log_variance"`.
#' @param method Optional lowercase method override. BayesC-like annotation
#'   models use `"sbayesc"`; SBayesRC uses `"sbayesrc"`. Log variance accepts
#'   `"sbayesc"` or `"sbayesr"`.
#' @param nit,nburn,nthin MCMC iteration controls.
#' @param ncores Number of OpenMP threads.
#' @param seed Sampler seed.
#' @param nchains Number of independent chains.
#' @param chain_seeds Optional integer seeds, one per chain.
#' @param keep_chains Logical; return compact per-chain summaries when
#'   supported.
#' @param convergence Convergence mode: `"auto"`, `"none"`, `"core"`, or
#'   `"extended"`. Automatic mode remains core-only.
#' @param convergence_control Optional uniquely named convergence-control list,
#'   including extended groups and diagnostic trace-memory guards.
#' @param memory_warning_gb Analytical memory-warning threshold in GiB.
#' @param verbose Print resolved execution metadata.
#' @param updateLDswap Logical; request optional LD-swap/MH. This is supported
#'   for the current annotation-aware CSR models, including SBayesRC.
#' @param maf_effect_s Optional fixed global BayesS-style MAF-scaling parameter.
#' @param effect_maf Optional MAF aligned to the final summary-marker order.
#' @param allow_reference_maf_for_maf_effect_s Allow explicit reference-panel
#'   MAF fallback when no GWAS-summary or by-construction analysis MAF exists.
#'   Currently supported only for `annotation_model = "annotation_probit_stick"` in this
#'   unified annotation interface. The default `maf_effect_s = NULL` with
#'   `estimate_maf_effect_s = FALSE` fits the ordinary model. A finite numeric
#'   scalar with `estimate_maf_effect_s = FALSE` fits a fixed global-S model.
#'   `maf_effect_s` must remain `NULL` when `estimate_maf_effect_s = TRUE`;
#'   fixed and sampled S cannot both be requested.
#' @param estimate_maf_effect_s Logical; estimate one trait-specific
#'   BayesS-style `maf_effect_s` by Metropolis-Hastings. Currently supported
#'   only for `annotation_model = "annotation_probit_stick"`. Sampled
#'   `maf_effect_s` is not supported for `annotation_model = "fixed_marker"`,
#'   `"learned_logistic"`, or `"group"`.
#' @param maf_effect_s_init Initial value for sampled `maf_effect_s`. Defaults to
#'   0 and is used only when `estimate_maf_effect_s = TRUE`.
#' @param maf_effect_s_prior Numeric length-2 lower and upper bounds for the
#'   uniform sampled-`maf_effect_s` prior. Defaults to `c(-3, 2)` and is used
#'   only when `estimate_maf_effect_s = TRUE`.
#' @param maf_effect_s_proposal_sd Random-walk proposal standard deviation for
#'   sampled `maf_effect_s`. Defaults to 0.35 and is used only when
#'   `estimate_maf_effect_s = TRUE`.
#' @param ld_swap_prob Probability per MCMC iteration of attempting LD-swap
#'   moves when `updateLDswap = TRUE`.
#' @param ld_swap_r2 Minimum LD r-squared for candidate swap partners.
#' @param ld_swap_max_friends Maximum number of high-LD friends stored per
#'   marker for swap proposals.
#' @param ld_swap_moves Number of swap attempts when LD-swap is triggered.
#' @param theta_prior_sd Positive fixed Gaussian prior SD for log-variance
#'   annotation coefficients. The validated version-1 default is 0.7.
#' @param theta_init Optional annotation-by-trait initial coefficient matrix.
#' @param updateTheta Logical; update log-variance coefficients with elliptical
#'   slice sampling. Set `FALSE` to hold `theta_init` fixed.
#' @param h2 Requested initial expected genetic-variance fraction under the
#'   resolved marker, group, annotation, component, and MAF-S prior weights
#'   used by the selected model.
#' @param adjE Residual adjustment factor.
#' @param updateB,updateE,updatePi Logical sampler update controls. `updatePi`
#'   is used by the BayesC-like annotation models and ignored by SBayesRC.
#' @param ld_prefix Optional prefix of the disk-backed CSR LD files.
#' @param ... Additional model-specific arguments passed to the internal
#'   annotation backend. For `"annotation_probit_stick"`, use `mixture_var`
#'   for the ordered component variance multipliers.
#'
#' @return A formatted ST-BLR fit. All returned fits contain `dm`
#'   (marker-by-trait posterior inclusion or non-null probability), `bm`
#'   (marker-by-trait posterior mean effects), `vbs`, `vgs`, `ves`, `vle`,
#'   `vld`, and `input`, plus standardized metadata fields (`method`, `model`,
#'   `backend`, `data_level`, `annotation_model`, `annotations`, `nchains`,
#'   `keep_chains`). The `vle` and `vld` traces follow the same definitions and
#'   formatting conventions as annotation-unaware CSR fits.
#'
#'   Annotation-aware outputs are model-specific and may include
#'   `annotation_summary`, `annotation_pi`, `annotation_effects`,
#'   `annotation_prior`, `alpha`, and `sigmaSqAlpha`. Log-variance fits include
#'   `theta`, `theta_summary`, `annotation_variance_ratio`,
#'   `annotation_transform`, and `marker_prior_scale`. For a binary annotation,
#'   `exp(theta)` is the conditional annotated/unannotated prior-variance ratio;
#'   for a continuous annotation it is the ratio for a one-SD increase.
#'   SBayesRC fits also
#'   include BayesR-style `component_probabilities`: marker-by-component
#'   posterior probabilities by trait, including component zero. The `dm`
#'   field is the posterior non-null probability.
#'
#'   LD-swap-enabled fits include `fit$diagnostics$ld_swap` and, where chain
#'   summaries are retained, `fit$diagnostics$ld_swap_chains` or chain-level
#'   LD-swap entries. With
#'   `keep_chains = TRUE`, compact per-chain summaries are returned in
#'   `chains` for supported CSR annotation backends.
#'
#'   For sampled SBayesRC `maf_effect_s`, the fit also contains `maf_effect_s`,
#'   `maf_effect_s_chain_mean_sd`, `maf_effect_s_chain_mean_min`,
#'   `maf_effect_s_chain_mean_max`,
#'   `maf_effect_s_trace`, and `maf_effect_s_acceptance`. `maf_effect_s_trace` is
#'   an iteration x trait matrix, `maf_effect_s` is the posterior mean by trait,
#'   and `maf_effect_s_acceptance` is the MH acceptance rate by trait. With
#'   `keep_chains = TRUE`, chain-level sampled-S output is available as
#'   the matching flat trait-by-chain record's `maf_effect_s` and
#'   `maf_effect_s_acceptance` fields.
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
#' control. The `"fixed_marker"`, `"learned_logistic"`, and `"group"` policies
#' affect BayesC-like inclusion probabilities and/or marker-effect variance
#' priors directly.
#' SBayesRC uses annotations to affect component probabilities; `maf_effect_s`,
#' when requested, affects only marker-specific effect-size prior variance.
#'
#' The model-specific adapters are internal. `stblr_csr_annot()` is the sole
#' public annotation-aware CSR fitting entry.
#'
#' CSR effects are on the standardized-genotype scale. The BayesS-style
#' MAF-dependent prior scale used by fixed and sampled `maf_effect_s` is
#' `q_j(S) = h_j^(S + 1)`, where `h_j = 2 p_j (1 - p_j)` and `p_j` is the
#' minor allele frequency. The `+1` exponent appears because the sampler effects
#' are standardized-genotype-scale effects rather than allele-scale effects.
#'
#' For CSR SBayesRC, annotations affect component probabilities and
#' `maf_effect_s` affects marker-specific effect-size prior variance. Fixed
#' `maf_effect_s` uses
#' `b_j | component_j = m, vb, S ~ N(0, vb * mixture_var_m * q_j(S))`. The
#' null component has multiplier zero, and `dm` is the posterior non-null
#' probability.
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
#'   annotation_model = "annotation_probit_stick",
#'   estimate_maf_effect_s = TRUE,
#'   maf_effect_s_prior = c(-3, 2),
#'   maf_effect_s_proposal_sd = 0.35
#' )
#' }
#'
#' @export
stblr_csr_annot <- function(
  stats,
  Glist = NULL,
  annotations,
  annotation_model = c("fixed_marker", "group", "learned_logistic",
                       "annotation_probit_stick", "log_variance"),
  method = NULL,
  nit = 1000,
  nburn = 100,
  nthin = 1,
  seed = 1,
  nchains = 1L,
  ncores = 1L,
  chain_seeds = NULL,
  keep_chains = FALSE,
  convergence = c("auto", "none", "core", "extended"),
  convergence_control = NULL,
  memory_warning_gb = 8,
  verbose = FALSE,
  h2 = 0.3,
  adjE = 0.9,
  updateB = TRUE,
  updateE = TRUE,
  updatePi = TRUE,
  updateLDswap = FALSE,
  maf_effect_s = NULL,
  effect_maf = NULL,
  allow_reference_maf_for_maf_effect_s = FALSE,
  estimate_maf_effect_s = FALSE,
  maf_effect_s_init = 0,
  maf_effect_s_prior = c(-3, 2),
  maf_effect_s_proposal_sd = 0.35,
  ld_swap_prob = 0.05,
  ld_swap_r2 = 0.8,
  ld_swap_max_friends = 50L,
  ld_swap_moves = 1L,
  theta_prior_sd = 0.7,
  theta_init = NULL,
  updateTheta = TRUE,
  ld_prefix = NULL,
  ...
) {
 .blr_validate_exact_public_call(
  sys.call(), sys.function(), "stblr_csr_annot()")
 raw_capture <- .blr_begin_st_raw_capture()
 on.exit(.blr_end_st_raw_capture(raw_capture), add = TRUE)
 if (missing(annotations)) {
  stop("annotations must be supplied.")
 }

 annotation_policy <- .stblr_match_annotation_model(annotation_model[[1L]])
 annotation_model <- .stblr_annotation_backend(annotation_policy)
 chain <- .blr_chain_controls(
  nit, nburn, nthin, seed, nchains, ncores, chain_seeds, keep_chains)
 conv <- .blr_convergence_controls(
  convergence, convergence_control, chain$nchains)
 extra <- list(...)
 if (annotation_model == "learned") {
  .stblr_validate_learned_forwarded_args(extra)
 }
 annotation_targets <- list(
  stblr_csr_prior_annot, stblr_csr_learn_annot,
  stblr_csr_group_annot, stblr_csr_sbayesrc_generic,
  .stblr_csr_logvar_bayesc, .stblr_csr_logvar_bayesr)
 annotation_dot_names <- unique(c(
  "mixture_var",
  unlist(lapply(annotation_targets, function(fun) names(formals(fun))),
         use.names = FALSE)))
 extra <- .blr_capture_forwarded_args(
  extra,
  accepted = setdiff(annotation_dot_names, c(
   "", "...", "stats", "Glist", "ld_prefix", "annotations", "annotation",
   "annotation_info", "A", "group", "method", "nit", "nburn", "nthin",
   "seed", "nchains", "ncores", "chain_seeds", "keep_chains", "h2",
   "adjE", "updateB", "updateE", "updatePi", "updateLDswap",
   "maf_effect_s", "estimate_maf_effect_s", "maf_effect_s_init",
   "maf_effect_s_prior", "maf_effect_s_proposal_sd", "ld_swap_prob",
   "ld_swap_r2", "ld_swap_max_friends", "ld_swap_moves",
   "theta_prior_sd", "theta_init", "updateTheta", ".convergence_spec")),
  what = "stblr_csr_annot(...)"
 )
 marker_ids <- names(stats$ww[[1L]]) %||% stats$marker_id
 if (!length(marker_ids)) marker_ids <- paste0("V", seq_along(stats$ww[[1L]]))
 trait_ids <- stats$trait_names %||% names(stats$wy)
 if (is.null(trait_ids) || length(trait_ids) != length(stats$wy) ||
     anyNA(trait_ids) || any(!nzchar(trait_ids)) || anyDuplicated(trait_ids)) {
  trait_ids <- paste0("T", seq_along(stats$wy))
 }
 if (annotation_model == "logvar") {
  method <- method %||% "sbayesc"
  if (length(method) != 1L || is.na(method) ||
      !method %in% c("sbayesc", "sbayesr")) {
   stop("annotation_model = \"log_variance\" requires method = \"sbayesc\" or \"sbayesr\".",
        call. = FALSE)
  }
  if (!is.numeric(theta_prior_sd) || length(theta_prior_sd) != 1L ||
      !is.finite(theta_prior_sd) || theta_prior_sd <= 0) {
   stop("theta_prior_sd must be a positive finite scalar.", call. = FALSE)
  }
  if (!is.logical(updateTheta) || length(updateTheta) != 1L || is.na(updateTheta)) {
   stop("updateTheta must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.null(maf_effect_s) || isTRUE(estimate_maf_effect_s)) {
   stop("log_variance version 1 does not support maf_effect_s.", call. = FALSE)
  }
  annotation_info <- .stblr_preprocess_logvar_annotations(annotations, marker_ids)
  trace_spec <- .blr_st_native_trace_spec(
   conv, marker_ids, if (method == "sbayesr") "bayesr" else "bayesc",
   annotations = TRUE,
   component_count = if (method == "sbayesr")
    length(extra$mixture_var %||% c(0, 0.01, 0.1, 1)) else 0L,
   annotation_quantity_count = ncol(annotation_info$X))
  memory <- .blr_st_preflight_memory(
   stats = stats, operator = "csr", chain = chain, conv = conv,
   memory_warning_gb = memory_warning_gb, trace_spec = trace_spec)
  resolved_spec <- resolve_blr_spec_from_wrapper(
   method, "csr", trait_ids, marker_ids, chain, sample_sizes = stats$n,
   probability_policy = "global",
   marker_scale_policy = "annotation_log_variance",
   update_flags = list(
    marker_effects = isTRUE(updateB), residual_variance = isTRUE(updateE),
    probability = isTRUE(updatePi), annotation_scale = isTRUE(updateTheta)),
   migration_actions = attr(extra, "migration_actions") %||% character(),
   execution_contract_version = 0L
  )
  .validate_ld_swap_args(
   updateLDswap, ld_swap_prob, ld_swap_r2, ld_swap_max_friends, ld_swap_moves)
  ld_prefix <- .stblr_resolve_csr_annotation_ld_prefix(Glist, ld_prefix)
  args <- c(list(
   stats = stats, ld_prefix = ld_prefix, annotation_info = annotation_info,
   theta_prior_sd = theta_prior_sd, theta_init = theta_init,
   updateTheta = updateTheta, h2 = h2, adjE = adjE,
   nit = chain$nit, nburn = chain$nburn, nthin = chain$nthin,
   ncores = chain$ncores, seed = chain$seed_native, nchains = chain$nchains,
   keep_chains = chain$keep_chains || conv$compute || conv$keep_traces,
   chain_seeds = if (length(chain$chain_seeds_native))
    chain$chain_seeds_native else NULL,
   updateB = updateB, updateE = updateE, updatePi = updatePi,
   updateLDswap = updateLDswap, ld_swap_prob = ld_swap_prob,
   ld_swap_r2 = ld_swap_r2, ld_swap_max_friends = ld_swap_max_friends,
   ld_swap_moves = ld_swap_moves, .convergence_spec = trace_spec), extra)
  fit <- if (method == "sbayesr") {
   do.call(.stblr_csr_logvar_bayesr, args)
  } else {
   do.call(.stblr_csr_logvar_bayesc, args)
  }
  attr(fit, "blr_resolved_spec") <- resolved_spec
  logvar_diagnostics <- fit$diagnostics$logvar
  out <- .blr_finalize_st_public(
   fit, method, "csr", chain, conv, memory_warning_gb, verbose, memory)
  out$diagnostics$logvar <- logvar_diagnostics
  identity <- paste0(method, "_logvar")
  out$model <- identity
  out$annotation_model <- "log_variance"
  out$input$method <- method
  out$input$model <- identity
  out$input$backend <- paste0("csr_logvar_", sub("sbayes", "bayes", method))
  out$input$annotation_model <- "log_variance"
  out$input$annotation_transform <- annotation_info$transform
  out$input$annotation_marker_alignment_status <- annotation_info$marker_alignment
  out <- .stblr_attach_csr_operator_contract(out, Glist, ld_prefix)
  return(out)
 }
 component_count <- if (annotation_model == "sbayesrc")
  length(extra$mixture_var %||% c(0, 0.01, 0.1, 1)) else 0L
 annotation_quantity_count <- switch(
  annotation_model,
  prior = 0L,
  group = 2L * length(unique(as.character(annotations))),
  learned = 2L * (ncol(as.matrix(annotations)) +
    as.integer(isTRUE(extra$add_intercept %||% FALSE))),
  sbayesrc = {
   processed_columns <- ncol(as.matrix(annotations)) +
    as.integer(isTRUE(extra$add_intercept %||% TRUE))
   processed_columns * (component_count - 1L) + (component_count - 1L)
  })
 trace_spec <- .blr_st_native_trace_spec(
  conv, marker_ids,
  if (annotation_model == "sbayesrc") "bayesrc" else "bayesc",
  annotations = annotation_model != "prior",
  component_count = component_count,
  annotation_quantity_count = annotation_quantity_count)
 memory <- .blr_st_preflight_memory(
  stats = stats, operator = "csr", chain = chain, conv = conv,
  memory_warning_gb = memory_warning_gb, trace_spec = trace_spec)
 .stblr_check_annotation_method(method, annotation_model)
 if (annotation_model %in% c("prior", "learned", "group")) {
  if (isTRUE(estimate_maf_effect_s)) {
   stop(
    "estimate_maf_effect_s is currently supported only for annotation_model = \"annotation_probit_stick\".",
    call. = FALSE)
  }
  if (!is.null(maf_effect_s)) {
   stop(
    "maf_effect_s is currently supported only for annotation_model = \"annotation_probit_stick\".",
    call. = FALSE)
  }
 }
 .stblr_validate_sampled_maf_effect_s(
  estimate_maf_effect_s = estimate_maf_effect_s,
  maf_effect_s = maf_effect_s,
  maf_effect_s_init = maf_effect_s_init,
  maf_effect_s_prior = maf_effect_s_prior,
  maf_effect_s_proposal_sd = maf_effect_s_proposal_sd
 )
 maf_info <- .blr_resolve_st_effect_maf(
  effect_maf, allow_reference_maf_for_maf_effect_s,
  !is.null(maf_effect_s) || isTRUE(estimate_maf_effect_s), stats, Glist)
 Glist <- maf_info$Glist
 .stblr_check_annotation_chains(annotation_model, nchains, keep_chains, chain_seeds)
 .validate_ld_swap_args(
  updateLDswap, ld_swap_prob, ld_swap_r2, ld_swap_max_friends, ld_swap_moves
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
  seed = chain$seed_native,
  nchains = chain$nchains,
  chain_seeds = if (length(chain$chain_seeds_native))
    chain$chain_seeds_native else NULL,
  keep_chains = chain$keep_chains || conv$compute || conv$keep_traces,
  .convergence_spec = trace_spec
 )
 probability_policy <- switch(
  annotation_model, prior = "fixed_marker", learned = "learned_logistic",
  group = "group", sbayesrc = "annotation_probit_stick")
 resolved_spec <- resolve_blr_spec_from_wrapper(
  if (annotation_model == "sbayesrc") "sbayesrc" else "sbayesc",
  "csr", trait_ids, marker_ids, chain, sample_sizes = stats$n,
  probability_policy = probability_policy,
  marker_scale_policy = if (annotation_model == "sbayesrc" &&
    (!is.null(maf_effect_s) || isTRUE(estimate_maf_effect_s)))
    "component_maf_s" else if (annotation_model == "sbayesrc")
      "component" else "unit",
  component_multipliers = if (annotation_model == "sbayesrc")
    extra$mixture_var %||% c(0, 0.01, 0.1, 1) else NULL,
  update_flags = list(
   marker_effects = isTRUE(updateB), residual_variance = isTRUE(updateE),
   probability = isTRUE(updatePi)),
  migration_actions = attr(extra, "migration_actions") %||% character(),
  execution_contract_version = if (annotation_model %in% c("prior", "learned"))
   1L else 0L
 )
 if (annotation_model %in% c("prior", "learned")) {
  common$.execution_contract <- .blr_native_execution_contract(resolved_spec)
 }
 finish <- function(fit) {
  attr(fit, "blr_resolved_spec") <- resolved_spec
  model <- if (annotation_model == "sbayesrc") "sbayesrc" else "sbayesc"
  out <- .blr_finalize_st_public(
   fit, model, "csr", chain, conv, memory_warning_gb, verbose, memory)
  out$input$annotation_policy <- probability_policy
  out$input$probability_policy <- probability_policy
  active_s <- !is.null(maf_effect_s) || isTRUE(estimate_maf_effect_s)
  out$input$prior_kernel <- if (model == "sbayesrc") "bayesrc" else "bayesc"
  out$input$effect_scale <- if (model == "sbayesrc") {
    if (active_s) "component_maf_s" else "component"
  } else "unit"
  out$input$effect_scale_policy <- out$input$effect_scale
  out$input$effect_maf_source <- maf_info$effect_maf_source
  out$input$effect_maf_population <- maf_info$effect_maf_population
  out$input$effect_maf_alignment_status <-
   maf_info$effect_maf_alignment_status
  out$input$effect_maf_fallback_used <- maf_info$effect_maf_fallback_used
  out$data$effect_scale <- out$input$effect_scale
  out$data$effect_maf_source <- maf_info$effect_maf_source
  out$data$effect_maf_population <- maf_info$effect_maf_population
  out$data$effect_maf_alignment_status <-
   maf_info$effect_maf_alignment_status
  out$data$effect_maf_fallback_used <- maf_info$effect_maf_fallback_used
  .stblr_attach_csr_operator_contract(out, Glist, ld_prefix)
 }
 common$updateLDswap <- updateLDswap
 common$ld_swap_prob <- ld_swap_prob
 common$ld_swap_r2 <- ld_swap_r2
 common$ld_swap_max_friends <- ld_swap_max_friends
 common$ld_swap_moves <- ld_swap_moves

 if (annotation_model %in% c("prior", "learned", "group")) {
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
 if ("gamma" %in% names(extra)) {
  stop("Use mixture_var, not the retired gamma argument.", call. = FALSE)
 }
 if ("mixture_var" %in% names(extra)) {
  extra$gamma <- extra$mixture_var
  extra$mixture_var <- NULL
 }
 args <- c(
  common,
  list(
   Glist = Glist,
   A = annotations,
   maf_effect_s = maf_effect_s,
   estimate_maf_effect_s = estimate_maf_effect_s,
   maf_effect_s_init = maf_effect_s_init,
   maf_effect_s_prior = maf_effect_s_prior,
   maf_effect_s_proposal_sd = maf_effect_s_proposal_sd
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
 if (length(method) != 1L || is.na(method)) {
  stop("method must be one canonical lowercase model identifier.", call. = FALSE)
 }

 if (annotation_model %in% c("prior", "learned", "group")) {
  if (!method %in% "sbayesc") {
   stop(
    "Summary-statistics BayesC annotation models require method = NULL or \"sbayesc\".",
    call. = FALSE
   )
  }
  return(invisible(TRUE))
 }

 if (!method %in% "sbayesrc") {
  stop(
   paste0(
    "annotation_model = \"annotation_probit_stick\" requires method = ",
    "NULL or \"sbayesrc\"."),
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
 annotation_model <- .stblr_match_annotation_backend(annotation_model)
 if (is.null(fit$input)) fit$input <- list()

 meta <- switch(
  annotation_model,
  prior = list(method = "sbayesc", model = "sbayesc", backend = "csr_prior_bayesc"),
  learned = list(method = "sbayesc", model = "sbayesc", backend = "csr_annot_bayesc"),
  group = list(method = "sbayesc", model = "sbayesc", backend = "csr_group_bayesc"),
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
  fit <- .stblr_add_prior_annotation_metadata(fit)
 } else if (annotation_model == "learned") {
  fit <- .stblr_add_learned_annotation_metadata(fit)
 } else if (annotation_model == "group") {
  fit <- .stblr_add_group_annotation_metadata(fit)
 } else {
  fit <- .stblr_add_sbayesrc_annotation_metadata(fit)
 }

 fit
}

.stblr_add_prior_annotation_metadata <- function(fit) {
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

.stblr_add_learned_annotation_metadata <- function(fit) {
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

.stblr_add_group_annotation_metadata <- function(fit) {
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

.stblr_add_sbayesrc_annotation_metadata <- function(fit) {
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
