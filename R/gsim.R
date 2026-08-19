#' Simulate genomic data from an in-memory matrix or a qgg Glist
#'
#' `gsim()` is the single public simulation entry point used by sblrbench. In
#' Glist mode, causal variants are selected from marker metadata before any
#' genotypes are read; only those causal columns are requested from
#' `qgg::getG()` to construct genetic values. Genome-wide genotypes are read in
#' chunks only when `compute_sumstats = TRUE`.
#'
#' @param Glist Optional qgg genotype-list object.
#' @param W Optional in-memory genotype matrix. At most one of `Glist` and `W`
#'   may be supplied. If both are NULL, independent binomial genotypes are
#'   simulated.
#' @param ids Optional sample IDs. In Glist mode, `n` randomly selected IDs are
#'   used when `ids` is NULL and `n` is smaller than the available sample size.
#' @param rsids Optional eligible marker IDs.
#' @param A Optional marker-by-annotation matrix. Row names, when present, are
#'   matched to marker IDs. Binary and continuous columns may be mixed.
#' @param architecture One of `"bayesc"`, `"bayesr"`,
#'   `"major_polygenic"`, `"maf_dependent"`, `"clustered"`, or `"fixed"`.
#' @param nt Number of traits.
#' @param h2 Target trait heritabilities.
#' @param vg Target genetic variances. Marker effects are rescaled to these
#'   realized variances when `scale_effects = TRUE`.
#' @param pi Mixture probabilities, with the null component first.
#' @param mixture_variances Relative mixture-component variances, with zero
#'   first.
#' @param marker_multipliers Optional positive marker-specific relative
#'   active-effect variance multipliers, \eqn{q_j}. A supplied vector must have
#'   unique nonempty names exactly matching the simulation markers and is
#'   aligned once to canonical marker order. The multiplier is applied only to
#'   non-null effect draws, separately from component probabilities and the
#'   component multiplier. `NULL` gives unit multipliers. Values are not
#'   clipped or normalized. The aligned vector and compact provenance are
#'   returned. SBayesRV is one possible research motivation for supplying this
#'   generic fixed variance truth; `gsim()` does not implement SBayesRV
#'   inference. Non-unit multipliers are not defined for
#'   `architecture = "fixed"`, where `beta` supplies realized effects directly.
#' @param n_causal Optional exact number of non-null markers. If supplied,
#'   markers are sampled without replacement using their annotation-informed
#'   non-null probabilities.
#' @param beta Fixed marker effects for `architecture = "fixed"`.
#' @param alpha Annotation effects, with dimensions `ncol(A)` by
#'   `length(pi)-1`. Effects act on centered annotations and continuation logits.
#' @param annotation_model `"none"` or `"sbayesrc"`.
#' @param rg,re Genetic-effect and residual correlation matrices (or scalar
#'   off-diagonal correlations for multiple traits).
#' @param maf Optional marker MAFs, named by marker ID or in marker order.
#' @param maf_exponent For `maf_dependent`, effect variance is proportional to
#'   `[2p(1-p)]^maf_exponent` before realized-variance calibration.
#' @param block_id Optional marker block labels for `clustered`.
#' @param n_hot_blocks Number of enriched blocks in the clustered architecture.
#' @param cluster_enrichment Active-odds multiplier in enriched blocks.
#' @param n,m Sample and marker counts when genotypes are simulated. In Glist
#'   mode, `n = NULL` uses all available samples.
#' @param maf_min,maf_max MAF range for simulated genotypes.
#' @param n_annotations Number of annotations to simulate when `A` is NULL.
#' @param annotation_types Simulated annotation types: `"binary"` and/or
#'   `"continuous"`.
#' @param annotation_prob Probability for simulated binary annotations.
#' @param seed Optional seed.
#' @param standardize_W Whether loaded causal genotypes are standardized.
#' @param scale_effects Whether effects are calibrated to `vg`.
#' @param return_genotypes Return causal genotype columns in `W_causal`.
#' @param return_marker_probabilities Return the full marker probability surface.
#' @param compute_sumstats Compute marginal GWAS statistics. Glist data are read
#'   in chunks for this optional step.
#' @param chunk_size Maximum marker columns requested in each Glist read.
#' @param getG_fun Optional replacement for `qgg::getG`, primarily for tests.
#'
#' @return A list containing phenotypes, genetic values, residuals, exact marker
#'   effects and states, annotation truth, probability surfaces, the complete
#'   canonical `marker_multipliers` vector, compact multiplier provenance under
#'   `settings$marker_multipliers`, and optional summary statistics.
#' @export
gsim <- function(
  Glist = NULL,
  W = NULL,
  ids = NULL,
  rsids = NULL,
  A = NULL,
  architecture = c("bayesc", "bayesr", "major_polygenic",
                   "maf_dependent", "clustered", "fixed"),
  nt = 1L,
  h2 = 0.5,
  vg = 1,
  pi = NULL,
  mixture_variances = NULL,
  marker_multipliers = NULL,
  n_causal = NULL,
  beta = NULL,
  alpha = NULL,
  annotation_model = c("none", "sbayesrc"),
  rg = NULL,
  re = 0,
  maf = NULL,
  maf_exponent = -0.5,
  block_id = NULL,
  n_hot_blocks = NULL,
  cluster_enrichment = 20,
  n = NULL,
  m = 1000L,
  maf_min = 0.05,
  maf_max = 0.5,
  n_annotations = 0L,
  annotation_types = NULL,
  annotation_prob = 0.1,
  seed = NULL,
  standardize_W = TRUE,
  scale_effects = TRUE,
  return_genotypes = FALSE,
  return_marker_probabilities = TRUE,
  compute_sumstats = FALSE,
  chunk_size = 5000L,
  getG_fun = NULL
) {
  architecture <- match.arg(architecture)
  annotation_model <- match.arg(annotation_model)
  if (!is.null(seed)) set.seed(seed)

  if (!is.null(Glist) && !is.null(W)) {
    .gsim_stop("Supply at most one of Glist and W.")
  }
  nt <- as.integer(nt)
  .gsim_check_scalar(nt, "nt", 1, Inf, integer = TRUE)
  if (!(length(h2) %in% c(1L, nt))) {
    .gsim_stop("h2 must be scalar or have length nt.")
  }
  if (!(length(vg) %in% c(1L, nt))) {
    .gsim_stop("vg must be scalar or have length nt.")
  }
  h2 <- rep(h2, length.out = nt)
  vg <- rep(vg, length.out = nt)
  if (any(!is.finite(h2) | h2 <= 0 | h2 >= 1)) {
    .gsim_stop("h2 must contain nt values strictly between zero and one.")
  }
  if (length(vg) != nt || any(!is.finite(vg) | vg <= 0)) {
    .gsim_stop("vg must contain nt positive finite values.")
  }
  rg <- .gsim_validate_corr(rg, nt, "rg")
  re <- .gsim_validate_corr(re, nt, "re")
  .gsim_check_scalar(chunk_size, "chunk_size", 1, Inf, integer = TRUE)

  mode <- if (!is.null(Glist)) "Glist" else if (!is.null(W)) "W" else "simulated"
  W_full_raw <- NULL
  chr_map <- NULL
  getG_resolved <- NULL

  if (mode == "Glist") {
    marker_ids <- .gsim_marker_ids_from_glist(Glist, rsids)
    available_ids <- as.character(Glist$ids %||%
      .gsim_stop("Glist must contain sample IDs in Glist$ids."))
    if (is.null(ids)) {
      if (!is.null(n)) {
        .gsim_check_scalar(n, "n", 1, length(available_ids), integer = TRUE)
        sample_ids <- sample(available_ids, n)
      } else {
        sample_ids <- available_ids
      }
    } else {
      sample_ids <- as.character(ids)
      if (anyNA(match(sample_ids, available_ids))) {
        .gsim_stop("Some ids are absent from Glist$ids.")
      }
    }
    chr_map <- .gsim_chr_map_from_glist(Glist, marker_ids)
    getG_resolved <- .gsim_resolve_getG(getG_fun)
  } else {
    if (mode == "simulated") {
      n <- n %||% 1000L
      .gsim_check_scalar(n, "n", 2, Inf, integer = TRUE)
      .gsim_check_scalar(m, "m", 1, Inf, integer = TRUE)
      if (maf_min <= 0 || maf_max > 0.5 || maf_min >= maf_max) {
        .gsim_stop("Require 0 < maf_min < maf_max <= 0.5.")
      }
      p <- stats::runif(m, maf_min, maf_max)
      W <- vapply(p, function(x) stats::rbinom(n, 2L, x), numeric(n))
      W <- matrix(W, nrow = n, ncol = m)
      colnames(W) <- paste0("m", seq_len(m))
      rownames(W) <- paste0("id", seq_len(n))
    }
    W_full_raw <- .gsim_prepare_W(W, ids)
    marker_ids <- colnames(W_full_raw)
    sample_ids <- rownames(W_full_raw)
    if (!is.null(rsids)) {
      pos <- match(as.character(rsids), marker_ids)
      if (anyNA(pos)) .gsim_stop("Some rsids are absent from W.")
      W_full_raw <- W_full_raw[, pos, drop = FALSE]
      marker_ids <- colnames(W_full_raw)
    }
    chr_map <- setNames(rep(NA_character_, length(marker_ids)), marker_ids)
  }

  m_total <- length(marker_ids)
  multiplier <- .gsim_prepare_marker_multipliers(
    marker_multipliers, marker_ids
  )
  marker_multipliers <- multiplier$value
  if (architecture == "fixed" && !multiplier$settings$all_ones) {
    .gsim_stop(
      "non-unit marker_multipliers are not defined for architecture = 'fixed'."
    )
  }
  annotation <- .gsim_prepare_annotations(
    A, marker_ids, n_annotations, annotation_types, annotation_prob
  )
  A <- annotation$A

  maf_all <- if (!is.null(maf)) {
    .gsim_as_named_vector(maf, marker_ids, "maf")
  } else if (!is.null(W_full_raw)) {
    setNames(.gsim_maf_from_W(W_full_raw), marker_ids)
  } else {
    setNames(rep(NA_real_, m_total), marker_ids)
  }
  block_id <- .gsim_as_named_vector(
    block_id, marker_ids, "block_id", numeric_only = FALSE
  )
  if (is.null(block_id) && architecture == "clustered") {
    if (all(is.na(chr_map))) {
      .gsim_stop("clustered architecture requires block_id or chromosome groups in Glist.")
    }
    block_id <- chr_map
  }

  hot_blocks <- NULL
  probability <- continuation <- NULL
  stick_order <- baseline_continuation <- NULL

  if (architecture == "fixed") {
    if (is.null(beta)) .gsim_stop("architecture = 'fixed' requires beta.")
    B <- .gsim_align_fixed_beta(beta, marker_ids, nt)
    component <- ifelse(rowSums(abs(B)) > 0, 2L, 1L)
    pi <- c(mean(component == 1L), mean(component != 1L))
    mixture_variances <- c(0, 1)
    alpha <- matrix(0, if (is.null(A)) 0L else ncol(A), 1L)
    if (!is.null(A)) rownames(alpha) <- colnames(A)
  } else {
    defaults <- .gsim_architecture_defaults(
      architecture, m_total, n_causal, pi, mixture_variances
    )
    pi <- defaults$pi
    mixture_variances <- defaults$mixture_variances
    alpha <- .gsim_prepare_alpha(alpha, A, length(pi) - 1L)

    if (annotation_model == "none" && any(alpha != 0)) {
      .gsim_stop("Non-zero alpha requires annotation_model = 'sbayesrc'.")
    }
    A_for_probability <- if (annotation_model == "sbayesrc") {
      annotation$A_centered
    } else {
      if (is.null(A)) NULL else matrix(0, nrow(A), ncol(A))
    }
    surface <- .gsim_probability_surface(
      pi, mixture_variances, A_for_probability, alpha
    )
    probability <- surface$probability
    continuation <- surface$continuation
    stick_order <- surface$stick_order
    baseline_continuation <- surface$baseline_continuation

    if (nrow(probability) == 1L && m_total > 1L) {
      probability <- probability[rep(1L, m_total), , drop = FALSE]
      continuation <- continuation[rep(1L, m_total), , drop = FALSE]
    }
    rownames(probability) <- marker_ids
    rownames(continuation) <- marker_ids

    if (architecture == "clustered") {
      clustered <- .gsim_apply_cluster_weights(
        probability, block_id, n_hot_blocks, cluster_enrichment
      )
      probability <- clustered$probability
      hot_blocks <- clustered$hot_blocks
      continuation <- .gsim_probability_to_continuation(
        probability, stick_order
      )
    }

    component <- .gsim_draw_components(probability, n_causal)
    causal_tmp <- which(component != 1L)
    if (!length(causal_tmp)) {
      .gsim_stop("The component draw contained no causal marker; increase pi or n_causal.")
    }
    B <- NULL
  }

  causal_idx <- which(component != 1L)
  if (!length(causal_idx)) {
    .gsim_stop("The architecture contains no causal marker.")
  }
  causal_ids <- marker_ids[causal_idx]

  if (mode == "Glist") {
    W_causal_raw <- .gsim_load_glist_markers(
      Glist, causal_ids, sample_ids, chr_map, getG_resolved, chunk_size
    )
  } else {
    W_causal_raw <- W_full_raw[, causal_idx, drop = FALSE]
  }

  maf_causal_observed <- .gsim_maf_from_W(W_causal_raw)
  missing_maf <- !is.finite(maf_all[causal_ids])
  maf_all[causal_ids[missing_maf]] <- maf_causal_observed[missing_maf]

  if (architecture != "fixed") {
    beta_maf <- if (architecture == "maf_dependent") maf_all else NULL
    B <- .gsim_draw_beta(
      component, mixture_variances, nt, rg,
      marker_multipliers = marker_multipliers,
      maf = beta_maf,
      maf_exponent = maf_exponent
    )
    rownames(B) <- marker_ids
  }
  colnames(B) <- paste0("D", seq_len(nt))

  W_causal <- .gsim_impute_and_standardize(
    W_causal_raw, standardize = standardize_W
  )
  B_causal <- B[causal_idx, , drop = FALSE]
  G <- W_causal %*% B_causal
  colnames(G) <- paste0("D", seq_len(nt))
  rownames(G) <- sample_ids

  raw_vg <- apply(G, 2L, stats::var)
  if (any(!is.finite(raw_vg) | raw_vg <= 0)) {
    .gsim_stop("The sampled effects produced zero or non-finite genetic variance.")
  }
  effect_scale <- rep(1, nt)
  if (scale_effects) {
    effect_scale <- sqrt(vg / raw_vg)
    B <- sweep(B, 2L, effect_scale, "*")
    B_causal <- B[causal_idx, , drop = FALSE]
    G <- W_causal %*% B_causal
  }

  realized_vg <- apply(G, 2L, stats::var)
  residual_var <- realized_vg * (1 - h2) / h2
  D_e <- diag(sqrt(residual_var), nrow = nt, ncol = nt)
  Sigma_e <- D_e %*% re %*% D_e
  E <- matrix(stats::rnorm(length(sample_ids) * nt), nrow = length(sample_ids),
              ncol = nt) %*% chol(Sigma_e)
  colnames(E) <- colnames(G)
  rownames(E) <- sample_ids
  Y <- G + E
  colnames(Y) <- colnames(G)
  rownames(Y) <- sample_ids

  sumstats <- NULL
  if (compute_sumstats) {
    if (mode == "Glist") {
      sumstats <- .gsim_stream_sumstats(
        Glist, marker_ids, sample_ids, chr_map, getG_resolved, Y,
        chunk_size, standardize_W
      )
    } else {
      W_stats <- .gsim_impute_and_standardize(
        W_full_raw, standardize = standardize_W
      )
      sumstats <- .gsim_compute_sumstats(W_stats, Y, marker_ids)
    }
  }

  exactness <- .gsim_exactness(G, E, Y, h2, probability)
  causal_table <- data.frame(
    rsid = causal_ids,
    component = component[causal_idx],
    maf = unname(maf_all[causal_ids]),
    stringsAsFactors = FALSE
  )
  for (t in seq_len(nt)) causal_table[[paste0("beta_D", t)]] <- B_causal[, t]

  out <- list(
    y = if (nt == 1L) Y[, 1L] else Y,
    Y = Y,
    G = G,
    E = E,
    B = B,
    B_causal = B_causal,
    causal = causal_table,
    causal_rsids = causal_ids,
    component = setNames(component, marker_ids),
    A = A,
    annotation_types = annotation$types,
    alpha = alpha,
    pi = pi,
    mixture_variances = mixture_variances,
    marker_multipliers = marker_multipliers,
    marker_probabilities = if (return_marker_probabilities) probability else NULL,
    continuation_probabilities = if (return_marker_probabilities) continuation else NULL,
    stick_order = stick_order,
    baseline_continuation = baseline_continuation,
    maf = maf_all,
    block_id = block_id,
    hot_blocks = hot_blocks,
    h2_target = h2,
    h2_observed = exactness$h2_observed,
    vg_target = vg,
    vg_observed = realized_vg,
    rg_target = rg,
    rg_observed = if (nt == 1L) matrix(1, 1L, 1L) else stats::cor(G),
    re_target = re,
    re_observed = if (nt == 1L) matrix(1, 1L, 1L) else stats::cor(E),
    Sigma_e = Sigma_e,
    sumstats = sumstats,
    exactness = exactness,
    settings = list(
      mode = mode,
      architecture = architecture,
      annotation_model = annotation_model,
      seed = seed,
      n = length(sample_ids),
      m = m_total,
      nt = nt,
      n_causal = length(causal_idx),
      standardize_W = standardize_W,
      scale_effects = scale_effects,
      effect_scale = effect_scale,
      marker_multipliers = multiplier$settings,
      chunk_size = chunk_size,
      compute_sumstats = compute_sumstats
    )
  )
  if (return_genotypes) out$W_causal <- W_causal
  class(out) <- c("gsim", "list")
  out
}

#' @export
print.gsim <- function(x, ...) {
  cat("gsim object\n")
  cat("  mode / architecture:", x$settings$mode, "/",
      x$settings$architecture, "\n")
  cat("  samples / markers / causal:", x$settings$n, "/",
      x$settings$m, "/", x$settings$n_causal, "\n")
  cat("  target h2:", paste(format(x$h2_target, digits = 3), collapse = ", "),
      "\n")
  invisible(x)
}
