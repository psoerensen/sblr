#' Summarise Post-Hoc MAF Architecture in an ST-BLR Fit
#'
#' Computes a lightweight post-hoc diagnostic for the relationship between
#' marker heterozygosity and posterior marker signal. This is a descriptive
#' regression on fitted marker summaries, not a sampler-level estimate of the
#' BayesS `S` parameter, not an MCMC estimate of a BayesS selection parameter,
#' and not a sampler-level prior.
#'
#' The diagnostic regresses either `log(PIP * bm^2 + epsilon)` or
#' `log(bm^2 + epsilon)` on `log(h)`, where `h = 2p(1 - p)`.
#'
#' @param fit Fitted ST-BLR object with marker-by-trait `dm` and `bm` fields.
#' @param maf Optional numeric vector of minor allele frequencies. Values in
#'   `(0, 0.5]` are accepted.
#' @param h Optional numeric vector of marker heterozygosities,
#'   `2 * maf * (1 - maf)`. If supplied, `maf` is ignored. If neither `h` nor
#'   `maf` is supplied, the helper looks for recoverable values in `fit`.
#' @param markers Optional marker subset. A character vector is applied to all
#'   traits. A named list can provide trait-specific character vectors.
#' @param min_pip Optional PIP threshold in `[0, 1]`; applied separately by
#'   trait after any explicit `markers` subset.
#' @param top_n Optional positive integer number of highest-PIP markers to keep
#'   per trait after any explicit `markers` and `min_pip` filters.
#' @param response Response used in the post-hoc regression. The default
#'   `"log_pip_weighted_bm2"` regresses `log(PIP * bm^2 + epsilon)`.
#'   `"log_bm2"` regresses `log(bm^2 + epsilon)`.
#' @param epsilon Positive scalar added before taking logs.
#'
#' @return A data frame with one row per trait and columns `trait`,
#'   `maf_effect_s_posthoc`, `se`, `p_value`, `r2`, `intercept`, `n_markers`,
#'   `n_effective_markers`, `method`, `response`, and `marker_filter`.
#'
#' @examples
#' fit <- list(
#'   dm = matrix(seq(0.1, 0.8, length.out = 8), ncol = 1,
#'               dimnames = list(paste0("m", 1:8), "trait1")),
#'   bm = matrix(seq(0.01, 0.08, length.out = 8), ncol = 1,
#'               dimnames = list(paste0("m", 1:8), "trait1"))
#' )
#' maf <- stats::setNames(seq(0.05, 0.4, length.out = 8), paste0("m", 1:8))
#' summarise_architecture(fit, maf = maf)
#' summarise_architecture(fit, maf = maf, min_pip = 0.01)
#' summarise_architecture(fit, maf = maf, top_n = 5)
#' causal_by_trait <- list(trait1 = paste0("m", 1:5))
#' summarise_architecture(fit, maf = maf, markers = causal_by_trait)
#'
#' @export
summarise_architecture <- function(
    fit,
    maf = NULL,
    h = NULL,
    markers = NULL,
    min_pip = NULL,
    top_n = NULL,
    response = c("log_pip_weighted_bm2", "log_bm2"),
    epsilon = 1e-12
) {
  if (!is.list(fit)) {
    stop("fit must be a list-like ST-BLR fit object.")
  }

  dm <- .stblr_maf_effect_s_as_matrix(fit$dm, "fit$dm")
  bm <- .stblr_maf_effect_s_as_matrix(fit$bm, "fit$bm")
  if (!identical(dim(dm), dim(bm))) {
    stop("fit$dm and fit$bm must have matching dimensions.")
  }
  if (!identical(rownames(dm), rownames(bm))) {
    stop("fit$dm and fit$bm must have matching marker names.")
  }
  marker_names <- rownames(dm)
  if (is.null(marker_names) || any(!nzchar(marker_names))) {
    marker_names <- paste0("marker", seq_len(nrow(dm)))
    rownames(dm) <- marker_names
    rownames(bm) <- marker_names
  }

  if (!is.numeric(epsilon) || length(epsilon) != 1L ||
      !is.finite(epsilon) || epsilon <= 0) {
    stop("epsilon must be a positive finite numeric scalar.")
  }
  response <- match.arg(response)
  if (!is.null(min_pip) &&
      (!is.numeric(min_pip) || length(min_pip) != 1L ||
       !is.finite(min_pip) || min_pip < 0 || min_pip > 1)) {
    stop("min_pip must be a finite numeric scalar in [0, 1].")
  }
  if (!is.null(top_n) &&
      (!is.numeric(top_n) || length(top_n) != 1L ||
       !is.finite(top_n) || top_n < 1 || top_n != as.integer(top_n))) {
    stop("top_n must be a positive integer scalar.")
  }
  markers <- .stblr_maf_effect_s_validate_markers(markers, marker_names)

  h <- .stblr_maf_effect_s_prepare_h(
    h = h,
    maf = maf,
    fit = fit,
    markers = marker_names,
    n_markers = nrow(dm)
  )
  log_h <- log(h)
  trait_names <- colnames(dm)
  if (is.null(trait_names) || any(!nzchar(trait_names))) {
    trait_names <- paste0("trait", seq_len(ncol(dm)))
  }

  rows <- lapply(seq_len(ncol(dm)), function(j) {
    trait <- trait_names[j]
    pip <- dm[, j]
    effect <- bm[, j]
    signal <- effect^2
    if (identical(response, "log_pip_weighted_bm2")) {
      signal <- pmax(pip, 0) * signal
    }

    idx <- .stblr_maf_effect_s_filter_indices(
      pip = pip,
      markers = markers,
      trait = trait,
      marker_names = marker_names,
      min_pip = min_pip,
      top_n = top_n
    )
    y <- log(signal + epsilon)
    keep <- is.finite(log_h[idx]) & is.finite(y[idx])
    n_effective <- sum(keep)

    intercept <- NA_real_
    slope <- NA_real_
    se <- NA_real_
    p_value <- NA_real_
    r2 <- NA_real_
    if (length(idx) >= 5L && n_effective >= 5L &&
        length(unique(log_h[idx][keep])) >= 2L) {
      fit_lm <- stats::lm(y[idx][keep] ~ log_h[idx][keep])
      fit_summary <- summary(fit_lm)
      coefs <- stats::coef(fit_summary)
      intercept <- unname(coefs[1L])
      slope <- unname(coefs[2L])
      se <- unname(coefs[2L, "Std. Error"])
      p_value <- unname(coefs[2L, "Pr(>|t|)"])
      r2 <- unname(fit_summary$r.squared)
    }

    data.frame(
      trait = trait,
      maf_effect_s_posthoc = slope,
      se = se,
      p_value = p_value,
      r2 = r2,
      intercept = intercept,
      n_markers = length(idx),
      n_effective_markers = n_effective,
      method = "posthoc_regression",
      response = response,
      marker_filter = .stblr_maf_effect_s_marker_filter_label(
        markers = markers,
        min_pip = min_pip,
        top_n = top_n
      ),
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, rows)
}

.stblr_maf_effect_s_as_matrix <- function(x, name) {
  if (is.null(x)) {
    stop(name, " must be present.")
  }
  x <- as.matrix(x)
  if (!is.numeric(x)) {
    stop(name, " must be numeric.")
  }
  if (nrow(x) == 0L || ncol(x) == 0L) {
    stop(name, " must have at least one marker and one trait.")
  }
  x
}

.stblr_maf_effect_s_prepare_h <- function(h, maf, fit, markers, n_markers) {
  if (is.null(h) && is.null(maf)) {
    recovered <- .stblr_maf_effect_s_recover_h_or_maf(fit)
    h <- recovered$h
    maf <- recovered$maf
  }
  if (is.null(h) && is.null(maf)) {
    stop("Either maf or h must be supplied or recoverable from fit.")
  }

  source <- "h"
  values <- h
  if (is.null(values)) {
    source <- "maf"
    values <- maf
  }

  if (!is.numeric(values)) {
    stop(source, " must be numeric.")
  }
  values <- as.numeric(values)
  names(values) <- names(if (source == "h") h else maf)

  values <- .stblr_maf_effect_s_align_vector(
    values = values,
    markers = markers,
    n_markers = n_markers,
    name = source
  )
  if (any(!is.finite(values))) {
    stop(source, " must contain only finite values.")
  }

  if (identical(source, "maf")) {
    if (any(values <= 0 | values > 0.5)) {
      stop("maf must contain values in (0, 0.5].")
    }
    values <- 2 * values * (1 - values)
  }

  if (any(values <= 0)) {
    stop("h must contain positive finite heterozygosity values.")
  }

  pmax(values, .Machine$double.eps)
}

.stblr_maf_effect_s_align_vector <- function(values, markers, n_markers, name) {
  if (is.null(markers)) {
    if (length(values) == 0L) {
      stop(name, " must not be empty.")
    }
    if (length(values) != n_markers) {
      stop(name, " must have length equal to the number of markers.")
    }
    return(values)
  }

  if (!is.null(names(values)) && any(nzchar(names(values)))) {
    missing_markers <- setdiff(markers, names(values))
    if (length(missing_markers) > 0L) {
      stop(name, " is missing values for markers in fit$dm.")
    }
    return(values[markers])
  }

  if (length(values) != length(markers)) {
    stop(name, " must have length equal to the number of markers.")
  }
  values
}

.stblr_maf_effect_s_recover_h_or_maf <- function(fit) {
  h_candidates <- list(
    fit$h,
    fit$heterozygosity,
    fit$input$h,
    fit$input$heterozygosity
  )
  for (candidate in h_candidates) {
    if (!is.null(candidate)) {
      return(list(h = candidate, maf = NULL))
    }
  }

  maf_candidates <- list(
    fit$maf,
    fit$input$maf
  )
  for (candidate in maf_candidates) {
    if (!is.null(candidate)) {
      return(list(h = NULL, maf = candidate))
    }
  }

  list(h = NULL, maf = NULL)
}

.stblr_maf_effect_s_validate_markers <- function(markers, marker_names) {
  if (is.null(markers)) {
    return(NULL)
  }
  if (is.character(markers)) {
    .stblr_maf_effect_s_check_marker_vector(markers, marker_names)
    return(markers)
  }
  if (is.list(markers) && !is.null(names(markers)) &&
      all(nzchar(names(markers)))) {
    lapply(markers, .stblr_maf_effect_s_check_marker_vector,
           marker_names = marker_names)
    return(markers)
  }
  stop("markers must be NULL, a character vector, or a named list of character vectors.")
}

.stblr_maf_effect_s_check_marker_vector <- function(markers, marker_names) {
  if (!is.character(markers) || anyNA(markers)) {
    stop("markers must contain marker names as character values.")
  }
  missing_markers <- setdiff(markers, marker_names)
  if (length(missing_markers) > 0L) {
    stop("markers contains values not found in fit$dm.")
  }
  invisible(markers)
}

.stblr_maf_effect_s_filter_indices <- function(pip, markers, trait, marker_names,
                                              min_pip, top_n) {
  idx <- seq_along(marker_names)

  if (!is.null(markers)) {
    marker_subset <- markers
    if (is.list(markers)) {
      marker_subset <- markers[[trait]]
      if (is.null(marker_subset)) {
        marker_subset <- character()
      }
    }
    idx <- match(marker_subset, marker_names)
  }

  if (!is.null(min_pip)) {
    idx <- idx[is.finite(pip[idx]) & pip[idx] >= min_pip]
  }

  if (!is.null(top_n) && length(idx) > top_n) {
    ord <- order(pip[idx], decreasing = TRUE, na.last = NA)
    idx <- idx[ord[seq_len(min(top_n, length(ord)))]]
  }

  idx
}

.stblr_maf_effect_s_marker_filter_label <- function(markers, min_pip, top_n) {
  parts <- character()
  if (!is.null(markers)) {
    parts <- c(parts, if (is.list(markers)) "markers=list" else "markers=character")
  }
  if (!is.null(min_pip)) {
    parts <- c(parts, paste0("min_pip=", format(min_pip, scientific = FALSE)))
  }
  if (!is.null(top_n)) {
    parts <- c(parts, paste0("top_n=", top_n))
  }
  if (length(parts) == 0L) {
    return("all")
  }
  paste(parts, collapse = ";")
}
