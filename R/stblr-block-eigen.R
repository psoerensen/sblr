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
#' @param annotations Annotation matrix for
#'   `annotation_model = "log_variance"`. Binary columns are centered only;
#'   continuous columns are centered and standardized to SD one.
#' @param annotation_model Optional annotation policy. The retained-block LV
#'   models use `"log_variance"`; `NULL` preserves ordinary model dispatch.
#' @param theta_prior_sd Positive finite Gaussian prior SD for log-variance
#'   annotation coefficients. The validated version-1 default is `0.7`.
#' @param theta_init Optional annotation-by-trait initial coefficient matrix.
#' @param updateTheta Whether to update log-variance coefficients by ESS.
#' @param representation Operator representation. `"low_rank"` is canonical;
#'   `"dense_reconstructed"` retains the historical packed dense operator.
#' @param eigen_policy Representation-specific eigenvalue policy, or `NULL` for
#'   the representation default.
#' @param eigen_prop Cumulative positive-eigenvalue mass target for low rank.
#' @param eigen_tau,eigen_eta Nonnegative reconstructed-dense filter controls.
#' @param low_rank_residual_rebuild_every Non-negative integer interval for
#'   rebuilding retained reduced residuals after completed MCMC iterations.
#'   The default is 100; zero disables periodic rebuilding after initialization.
#' @param residual_policy Residual-variance contract. `"gctb_block"` is the
#'   default for retained block-eigen SBayesR/SBayesRC, `"fixed_block"` keeps
#'   every block at the phenotype variance, and
#'   `"global_projected_legacy"` reproduces the historical projected-global
#'   experimental contract.
#' @param block_ve_mode Official block update mode for `"gctb_block"`.
#' @param resam_thresh Official effect-SS to block-Vg eligibility threshold.
#' @param minimum_ve_ratio Official sampled-block-Ve to phenotype-variance
#'   compatibility threshold.
#' @param block_ve_keep_history Retain complete per-chain block-Ve histories.
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
  annotation = NULL, annotations = NULL, annotation_model = NULL,
  theta_prior_sd = 0.7, theta_init = NULL, updateTheta = TRUE,
  representation = c("low_rank", "dense_reconstructed"),
  eigen_policy = NULL, eigen_prop = 0.995, eigen_tau = 0.01, eigen_eta = 0,
  low_rank_residual_rebuild_every = 100L,
  residual_policy = NULL,
  block_ve_mode = c("allMixVe", "mixVe", "samVe", "fixVe"),
  resam_thresh = 1.1, minimum_ve_ratio = 0.7,
  block_ve_keep_history = FALSE,
  nit = 1000, nburn = 500, nthin = 1, seed = 1,
  nchains = 1L, ncores = 1L, chain_seeds = NULL,
  keep_chains = FALSE, convergence = c("auto", "none", "core", "extended"),
  convergence_control = NULL, memory_warning_gb = 8,
  verbose = FALSE, ...
) {
  .blr_validate_exact_public_call(
    sys.call(), sys.function(), "stblr_block_eigen()")
  raw_capture <- .blr_begin_st_raw_capture()
  on.exit(.blr_end_st_raw_capture(raw_capture), add = TRUE)
  rebuild_interval_supplied <- !missing(low_rank_residual_rebuild_every)
  if (!is.null(annotation_model)) {
    if (length(annotation_model) != 1L || is.na(annotation_model) ||
        !identical(annotation_model, "log_variance")) {
      stop("annotation_model must be NULL or 'log_variance'.", call. = FALSE)
    }
    if (!is.null(annotation) && !is.null(annotations)) {
      stop("Supply log-variance annotations through only one of annotations or annotation.",
           call. = FALSE)
    }
    if (!is.numeric(theta_prior_sd) || length(theta_prior_sd) != 1L ||
        !is.finite(theta_prior_sd) || theta_prior_sd <= 0) {
      stop("theta_prior_sd must be a positive finite scalar.", call. = FALSE)
    }
    if (!is.logical(updateTheta) || length(updateTheta) != 1L ||
        is.na(updateTheta)) {
      stop("updateTheta must be TRUE or FALSE.", call. = FALSE)
    }
  }
  dots <- list(...)
  block_targets <- list(
    .stblr_csr_bayesc_block_eigen, .stblr_csr_bayesr_block_eigen,
    .stblr_csr_sbayesrc_block_eigen, .stblr_block_logvar_bayesc,
    .stblr_block_logvar_bayesr)
  block_dot_names <- unique(c(
    "eigen_filter",
    names(formals(.stblr_csr_sbayesrc_generic_impl)),
    unlist(lapply(block_targets, function(fun) names(formals(fun))),
           use.names = FALSE)))
  dots <- .blr_capture_forwarded_args(
    dots,
    accepted = setdiff(block_dot_names, c(
      "", "...", "stats", "Glist", "block_start", "representation",
      "eigen_policy", "eigen_prop", "eigen_tau",
      "eigen_eta", "low_rank_residual_rebuild_every", "residual_policy",
      "block_ve_mode", "resam_thresh", "minimum_ve_ratio",
      "block_ve_keep_history", "nit", "nburn", "nthin", "seed",
      "nchains", "ncores", "chain_seeds", "keep_chains",
      ".convergence_spec", ".native_fun", ".native_args", ".input_extra",
      ".return_raw", ".information_diagnostics",
      ".diagnostic_block_px", ".diagnostic_block_px_log_scale_sd",
      "ld_prefix", "annotation", "annotation_info", "A",
      "theta_prior_sd", "theta_init", "updateTheta")),
    what = "stblr_block_eigen(...)"
  )
  legacy_filter <- dots$eigen_filter
  dots$eigen_filter <- NULL
  if (!is.null(legacy_filter) && missing(representation)) {
    representation <- "dense_reconstructed"
  }
  representation <- match.arg(representation)
  if (!is.numeric(low_rank_residual_rebuild_every) ||
      length(low_rank_residual_rebuild_every) != 1L ||
      !is.finite(low_rank_residual_rebuild_every) ||
      low_rank_residual_rebuild_every < 0 ||
      low_rank_residual_rebuild_every != floor(low_rank_residual_rebuild_every)) {
    stop("low_rank_residual_rebuild_every must be a non-negative integer scalar.",
         call. = FALSE)
  }
  low_rank_residual_rebuild_every <-
    as.integer(low_rank_residual_rebuild_every)
  if (identical(representation, "low_rank")) {
    if (isTRUE(dots$use_r_init)) {
      stop(
        paste0(
          "r_init is not supported for representation = \"retained_low_rank\"; ",
          "a reduced-residual restart contract has not been implemented."
        ),
        call. = FALSE
      )
    }
    if (isTRUE(dots$rebuild_r_before_updateE)) {
      stop(
        paste0(
          "rebuild_r_before_updateE is incompatible with retained low rank; ",
          "use low_rank_residual_rebuild_every instead."
        ),
        call. = FALSE
      )
    }
  } else {
    if (rebuild_interval_supplied && low_rank_residual_rebuild_every != 0L) {
      stop(
        "low_rank_residual_rebuild_every is only supported by representation = 'low_rank'.",
        call. = FALSE
      )
    }
    low_rank_residual_rebuild_every <- 0L
  }
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
  logvar <- identical(annotation_model, "log_variance")
  if (logvar && !method %in% c("sbayesc", "sbayesr")) {
    stop("annotation_model = 'log_variance' requires method = 'sbayesc' or 'sbayesr'.",
         call. = FALSE)
  }
  if (logvar && (!is.null(dots$maf_effect_s) ||
                 isTRUE(dots$estimate_maf_effect_s))) {
    stop("maf_effect_s is not supported jointly with log-variance annotations.",
         call. = FALSE)
  }
  if (is.null(residual_policy)) {
    residual_policy <- if (
      identical(representation, "low_rank") &&
      method %in% c("sbayesr", "sbayesrc")) "gctb_block" else
        "global_projected_legacy"
  }
  residual_policy <- match.arg(
    residual_policy,
    c("gctb_block", "fixed_block", "global_projected_legacy"))
  block_ve_mode <- match.arg(block_ve_mode)
  if (!identical(representation, "low_rank") &&
      !identical(residual_policy, "global_projected_legacy")) {
    stop("Block residual policies require representation = 'low_rank'.",
         call. = FALSE)
  }
  if (identical(method, "sbayesc") &&
      !identical(residual_policy, "global_projected_legacy")) {
    stop("Block residual policies currently apply only to SBayesR/SBayesRC.",
         call. = FALSE)
  }
  if (identical(residual_policy, "fixed_block")) block_ve_mode <- "fixVe"
  if (!is.numeric(resam_thresh) || length(resam_thresh) != 1L ||
      !is.finite(resam_thresh) || resam_thresh <= 0) {
    stop("resam_thresh must be a positive finite scalar.", call. = FALSE)
  }
  if (!is.numeric(minimum_ve_ratio) || length(minimum_ve_ratio) != 1L ||
      !is.finite(minimum_ve_ratio) || minimum_ve_ratio <= 0) {
    stop("minimum_ve_ratio must be a positive finite scalar.", call. = FALSE)
  }
  if (!is.logical(block_ve_keep_history) ||
      length(block_ve_keep_history) != 1L || is.na(block_ve_keep_history)) {
    stop("block_ve_keep_history must be TRUE or FALSE.", call. = FALSE)
  }
  uses_block_ve <- residual_policy %in% c("gctb_block", "fixed_block")
  if (uses_block_ve) {
    if ("adjE" %in% names(dots) && !identical(as.numeric(dots$adjE), 0)) {
      stop("adjE is only available with residual_policy = 'global_projected_legacy'.",
           call. = FALSE)
    }
    if ("updateE_start" %in% names(dots) &&
        !identical(as.integer(dots$updateE_start), 0L)) {
      stop("updateE_start is only configurable for the legacy residual policy.",
           call. = FALSE)
    }
    if ("updateE_every" %in% names(dots) &&
        !identical(as.integer(dots$updateE_every), 1L)) {
      stop("updateE_every is only configurable for the legacy residual policy.",
           call. = FALSE)
    }
    dots$adjE <- 0
    if (identical(dots$updateE, FALSE)) block_ve_mode <- "fixVe"
  }
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
  trait_ids <- stats$trait_names %||% names(stats$wy)
  if (is.null(trait_ids) || length(trait_ids) != length(stats$wy) ||
      anyNA(trait_ids) || any(!nzchar(trait_ids)) || anyDuplicated(trait_ids)) {
    trait_ids <- paste0("T", seq_along(stats$wy))
  }
  annotation_info <- NULL
  if (logvar) {
    logvar_annotations <- annotations %||% annotation
    if (is.null(logvar_annotations)) {
      stop("annotations are required for annotation_model = 'log_variance'.",
           call. = FALSE)
    }
    annotation_info <- .stblr_preprocess_logvar_annotations(
      logvar_annotations, marker_ids)
  }
  trace_spec <- .blr_st_native_trace_spec(
    conv, marker_ids, resolved_model$prior_kernel,
    annotations = identical(method, "sbayesrc") && !logvar,
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
  execution_version <- if (!logvar && method %in% c("sbayesc", "sbayesr")) 1L else 0L
  resolved_spec <- resolve_blr_spec_from_wrapper(
    method, "block_eigen", trait_ids, marker_ids, chain,
    sample_sizes = stats$n,
    probability_policy = if (method == "sbayesrc")
      "annotation_probit_stick" else "global",
    marker_scale_policy = if (logvar) "annotation_log_variance" else
      resolved_model$effect_scale,
    residual_policy = residual_policy,
    component_multipliers = if (
      resolved_model$prior_kernel %in% c("bayesr", "bayesrc"))
      dots$mixture_var %||% c(0, 0.01, 0.1, 1) else NULL,
    update_flags = list(
      marker_effects = isTRUE(dots$updateB %||% TRUE),
      residual_variance = isTRUE(dots$updateE %||% TRUE),
      probability = isTRUE(dots$updatePi %||% TRUE)),
    numerical_controls = list(
      representation = representation, eigen_policy = eigen_policy,
      eigen_prop = eigen_prop, eigen_tau = eigen_tau, eigen_eta = eigen_eta,
      low_rank_residual_rebuild_every = low_rank_residual_rebuild_every),
    migration_actions = attr(dots, "migration_actions") %||% character(),
    execution_contract_version = execution_version
  )
  execution_contract <- .blr_native_execution_contract(resolved_spec)
  common <- list(
    stats = stats, Glist = Glist, block_start = block_start,
    representation = representation, eigen_policy = eigen_policy,
    eigen_prop = eigen_prop, eigen_filter = eigen_filter,
    eigen_tau = eigen_tau, eigen_eta = eigen_eta,
    low_rank_residual_rebuild_every = low_rank_residual_rebuild_every,
    residual_policy = residual_policy, block_ve_mode = block_ve_mode,
    resam_thresh = resam_thresh, minimum_ve_ratio = minimum_ve_ratio,
    block_ve_keep_history = block_ve_keep_history,
    nit = chain$nit, nburn = chain$nburn, nthin = chain$nthin,
    seed = chain$seed_native, nchains = chain$nchains, ncores = chain$ncores,
    chain_seeds = if (length(chain$chain_seeds_native))
      chain$chain_seeds_native else NULL,
    keep_chains = chain$keep_chains || conv$compute || conv$keep_traces,
    .convergence_spec = trace_spec,
    .execution_contract = execution_contract)
  fit <- if (logvar) {
    logvar_common <- c(common[setdiff(names(common), ".execution_contract")], list(
      annotation_info = annotation_info,
      theta_prior_sd = theta_prior_sd,
      theta_init = theta_init,
      updateTheta = updateTheta
    ))
    if (identical(method, "sbayesc")) {
      do.call(.stblr_block_logvar_bayesc, c(
        logvar_common[setdiff(
          names(logvar_common),
          c("residual_policy", "block_ve_mode", "resam_thresh",
            "minimum_ve_ratio", "block_ve_keep_history"))], dots))
    } else {
      do.call(.stblr_block_logvar_bayesr, c(logvar_common, dots))
    }
  } else switch(
    method,
    sbayesc = do.call(
      .stblr_csr_bayesc_block_eigen,
      c(common[setdiff(
        names(common),
        c("residual_policy", "block_ve_mode", "resam_thresh",
          "minimum_ve_ratio", "block_ve_keep_history"))], dots)),
    sbayesr = do.call(.stblr_csr_bayesr_block_eigen, c(common, dots)),
    sbayesrc = {
      if (is.null(annotation)) {
        stop("annotation is required for method = 'sbayesrc'.", call. = FALSE)
      }
      do.call(.stblr_csr_sbayesrc_block_eigen,
              c(common[c("stats", "Glist", "block_start", "representation",
                         "eigen_policy", "eigen_prop", "eigen_filter",
                         "eigen_tau", "eigen_eta",
                         "low_rank_residual_rebuild_every")],
                common[c("residual_policy", "block_ve_mode", "resam_thresh",
                         "minimum_ve_ratio", "block_ve_keep_history")],
                list(annotation = annotation),
                common[c("nit", "nburn", "nthin", "seed", "nchains",
                         "ncores", "chain_seeds", "keep_chains",
                         ".convergence_spec")], dots))
    })
  attr(fit, "blr_resolved_spec") <- resolved_spec
  logvar_diagnostics <- if (logvar) fit$diagnostics$logvar else NULL
  fit <- .blr_finalize_st_public(
    fit, method, "block_eigen", chain, conv, memory_warning_gb, verbose,
    memory)
  fit$data$operator <- fit$input[c(
    "operator_representation", "operator_contract", "operator_scale_contract",
    "eigen_policy", "eigen_prop", "eigen_tau", "eigen_eta", "eigen_blocks")]
  fit$input$operator_role <- if (identical(representation, "low_rank"))
    "canonical_scalable_summary_statistics" else
      "historical_reconstructed_block_operator"
  fit$input$operator_approximate <- TRUE
  fit$input$residual_policy <- residual_policy
  fit$input$block_ve_mode <- block_ve_mode
  fit$input$resam_thresh <- resam_thresh
  fit$input$minimum_ve_ratio <- minimum_ve_ratio
  fit$input$block_ve_definition <- if (uses_block_ve)
    "mean of chain-specific block residual variances" else
      "legacy global projected residual variance"
  fit$input$heritability_definition <- if (uses_block_ve)
    "sum of block genetic variances divided by phenotype variance" else
      "genetic variance divided by genetic plus residual variance"
  fit$data$operator$operator_role <- fit$input$operator_role
  fit$data$operator$operator_approximate <- TRUE
  fit$data$operator$limitation <- if (identical(representation, "low_rank"))
    paste("Projected block likelihood retaining the configured positive",
      "spectral mass; cross-block LD is omitted and residual variance follows",
      if (uses_block_ve) "the selected block-specific policy." else
        "the explicit legacy global-projected policy.")
    else paste("Block-reconstructed approximation retained for regression",
      "and reproducibility, not an exact full-LD reference.")
  fit$diagnostics$block_eigen <-
    fit$input$eigen_diagnostics %||% fit$diagnostics$block_eigen %||% NULL
  if (logvar) fit$diagnostics$logvar <- logvar_diagnostics
  block_diag <- fit$diagnostics$block_eigen$blocks %||% NULL
  if (!is.null(block_diag) && "retained_rank" %in% names(block_diag)) {
    fit$input$block_retained_ranks <- as.integer(block_diag$retained_rank)
    fit$input$block_sample_sizes <- matrix(
      rep(as.integer(stats$n), length.out = length(stats$yy)) %o%
        rep(1L, nrow(block_diag)),
      nrow = length(stats$yy))
  }
  fit$input$effect_scale <- if (logvar) "annotation_log_variance" else
    resolved_model$effect_scale
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
  if (logvar) {
    identity <- paste0(method, "_logvar")
    fit$model <- identity
    fit$annotation_model <- "log_variance"
    fit$input$model <- identity
    fit$input$annotation_model <- "log_variance"
    fit$input$annotation_transform <- annotation_info$transform
    fit$input$annotation_marker_alignment_status <-
      annotation_info$marker_alignment
  }
  fit
}
