#' Fit Annotation-Aware CSR ST-BLR Models
#'
#' `stblr_csr_annot()` is the unified public entry point for the current
#' annotation-aware CSR summary-statistics backends. It dispatches to the
#' existing explicit wrappers and standardizes fit metadata and annotation
#' output aliases while preserving wrapper-specific fields.
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
#' @param method Optional method override. BayesC-like annotation models accept
#'   `"bayesc"`/`"bayesC"`; SBayesRC accepts `"sbayesrc"` or `"bayesr"`.
#' @param nit,nburn,nthin MCMC iteration controls.
#' @param ncores Number of OpenMP threads.
#' @param seed Sampler seed.
#' @param nchains Number of independent chains.
#' @param chain_seeds Optional integer seeds, one per chain.
#' @param keep_chains Logical; return compact per-chain summaries when
#'   supported.
#' @param updateLDswap Logical; request optional LD-swap/MH. This is currently
#'   supported only for `annotation_model = "prior"` among annotation-aware CSR
#'   models.
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
#' @return A formatted ST-BLR fit. All returned fits contain `dm`, `bm`, and
#'   `input`, standardized metadata fields (`method`, `model`, `backend`,
#'   `data_level`, `annotation_model`, `annotations`, `nchains`,
#'   `keep_chains`), and model-specific annotation aliases when they can be
#'   constructed reliably. Existing wrapper-specific fields are preserved.
#'
#' @details
#' This interface is CSR summary-statistics only. The current annotation-aware
#' SBayesRC and the BayesC-like annotation-aware CSR backends support multiple
#' native chains and compact chain summaries. LD-swap/MH is currently available
#' for the fixed-prior annotation-aware CSR BayesC backend only; learned,
#' group, and SBayesRC annotation models error if `updateLDswap = TRUE`.
#'
#' Existing wrappers remain supported:
#' [stblr_csr_prior_annot()], [stblr_csr_learn_annot()],
#' [stblr_csr_group_annot()], and [stblr_csr_sbayesrc_generic()].
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
  ncores = 1,
  seed = 1,
  nchains = 1L,
  chain_seeds = NULL,
  keep_chains = FALSE,
  h2 = 0.3,
  adjE = 0.9,
  updateB = TRUE,
  updateE = TRUE,
  updatePi = TRUE,
  updateLDswap = FALSE,
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
 .stblr_check_annotation_method(method, annotation_model)
 .stblr_check_annotation_chains(annotation_model, nchains, keep_chains, chain_seeds)
 .validate_ld_swap_args(
  updateLDswap, ld_swap_prob, ld_swap_r2, ld_swap_max_friends, ld_swap_moves
 )
 if (annotation_model != "prior" && isTRUE(updateLDswap)) {
  stop(
   "LD-swap/MH is currently implemented only for annotation_model = 'prior' among annotation-aware CSR models.",
   call. = FALSE
  )
 }
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
  nit = nit,
  nburn = nburn,
  nthin = nthin,
  ncores = ncores,
  seed = seed,
  nchains = nchains,
  chain_seeds = chain_seeds,
  keep_chains = keep_chains
 )
 extra <- list(...)

 if (annotation_model %in% c("prior", "learned", "group")) {
 common$updatePi <- updatePi
  common$updateLDswap <- updateLDswap
  common$ld_swap_prob <- ld_swap_prob
  common$ld_swap_r2 <- ld_swap_r2
  common$ld_swap_max_friends <- ld_swap_max_friends
  common$ld_swap_moves <- ld_swap_moves
 }

 if (annotation_model == "prior") {
  args <- .stblr_prior_annotation_args(annotations, extra)
  args <- c(common, args)
  return(do.call(stblr_csr_prior_annot, args))
 }

 if (annotation_model == "learned") {
  if ("A" %in% names(extra)) {
   stop("Supply learned annotations through annotations, not both annotations and A.")
  }
  args <- c(common, list(A = annotations), extra)
  return(do.call(stblr_csr_learn_annot, args))
 }

 if (annotation_model == "group") {
  if ("group" %in% names(extra)) {
   stop("Supply group annotations through annotations, not both annotations and group.")
  }
  args <- c(common, list(group = annotations), extra)
  return(do.call(stblr_csr_group_annot, args))
 }

 if ("A" %in% names(extra)) {
  stop("Supply SBayesRC annotations through annotations, not both annotations and A.")
 }
 args <- c(common, list(A = annotations), extra)
 do.call(stblr_csr_sbayesrc_generic, args)
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

 if (!method %in% c("sbayesrc", "bayesr")) {
  stop(
   "annotation_model = \"sbayesrc\" requires method = NULL, \"sbayesrc\", or \"bayesr\".",
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
