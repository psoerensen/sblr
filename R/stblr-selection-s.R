#' Summarise Post-Hoc MAF Architecture in an ST-BLR Fit
#'
#' Computes a lightweight post-hoc diagnostic for the relationship between
#' marker heterozygosity and posterior marker signal. This is a descriptive
#' regression on fitted marker summaries, not an MCMC estimate of a BayesS
#' selection parameter and not a sampler-level prior.
#'
#' The diagnostic regresses either `log(PIP * bm^2 + epsilon)` or
#' `log(bm^2 + epsilon)` on `log(h)`, where `h = 2p(1 - p)`.
#'
#' @param fit Fitted ST-BLR object with marker-by-trait `dm` and `bm` fields.
#' @param maf Optional numeric vector of minor allele frequencies. Values in
#'   `[0, 1]` are accepted and converted to minor-allele frequencies with
#'   `pmin(maf, 1 - maf)`.
#' @param h Optional numeric vector of marker heterozygosities,
#'   `2 * maf * (1 - maf)`. If supplied, `maf` is ignored.
#' @param epsilon Positive scalar added before taking logs.
#' @param use_pip_weights Logical. If `TRUE`, regress
#'   `log(PIP * bm^2 + epsilon)`; otherwise regress `log(bm^2 + epsilon)`.
#'
#' @return A data frame with one row per trait and columns `trait`,
#'   `selection_s_posthoc`, `intercept`, `n_markers`, `n_effective_markers`,
#'   `method`, and `response`.
#'
#' @examples
#' fit <- list(
#'   dm = matrix(c(0.2, 0.8, 0.5), ncol = 1,
#'               dimnames = list(paste0("m", 1:3), "trait1")),
#'   bm = matrix(c(0.01, 0.04, -0.02), ncol = 1,
#'               dimnames = list(paste0("m", 1:3), "trait1"))
#' )
#' summarise_stblr_maf_architecture(fit, maf = c(0.05, 0.2, 0.4))
#'
#' @export
summarise_stblr_maf_architecture <- function(
    fit,
    maf = NULL,
    h = NULL,
    epsilon = 1e-12,
    use_pip_weights = TRUE
) {
  if (!is.list(fit)) {
    stop("fit must be a list-like ST-BLR fit object.")
  }

  dm <- .stblr_selection_s_as_matrix(fit$dm, "fit$dm")
  bm <- .stblr_selection_s_as_matrix(fit$bm, "fit$bm")
  if (!identical(dim(dm), dim(bm))) {
    stop("fit$dm and fit$bm must have matching dimensions.")
  }

  if (!is.numeric(epsilon) || length(epsilon) != 1L ||
      !is.finite(epsilon) || epsilon <= 0) {
    stop("epsilon must be a positive finite numeric scalar.")
  }
  if (!is.logical(use_pip_weights) || length(use_pip_weights) != 1L ||
      is.na(use_pip_weights)) {
    stop("use_pip_weights must be TRUE or FALSE.")
  }

  h <- .stblr_selection_s_prepare_h(
    h = h,
    maf = maf,
    markers = rownames(dm),
    n_markers = nrow(dm)
  )
  log_h <- log(h)
  trait_names <- colnames(dm)
  if (is.null(trait_names) || any(!nzchar(trait_names))) {
    trait_names <- paste0("trait", seq_len(ncol(dm)))
  }

  rows <- lapply(seq_len(ncol(dm)), function(j) {
    pip <- dm[, j]
    effect <- bm[, j]
    signal <- effect^2
    response <- "log_bm2"
    if (use_pip_weights) {
      signal <- pmax(pip, 0) * signal
      response <- "log_pip_weighted_bm2"
    }

    y <- log(signal + epsilon)
    keep <- is.finite(log_h) & is.finite(y)
    n_effective <- sum(keep)

    intercept <- NA_real_
    slope <- NA_real_
    if (n_effective >= 2L && length(unique(log_h[keep])) >= 2L) {
      design <- cbind(intercept = 1, log_h = log_h[keep])
      coefs <- stats::lm.fit(design, y[keep])$coefficients
      intercept <- unname(coefs[1L])
      slope <- unname(coefs[2L])
    }

    data.frame(
      trait = trait_names[j],
      selection_s_posthoc = slope,
      intercept = intercept,
      n_markers = length(h),
      n_effective_markers = n_effective,
      method = "posthoc_regression",
      response = response,
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, rows)
}

.stblr_selection_s_as_matrix <- function(x, name) {
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

.stblr_selection_s_prepare_h <- function(h, maf, markers, n_markers) {
  if (is.null(h) && is.null(maf)) {
    stop("Either maf or h must be supplied.")
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

  values <- .stblr_selection_s_align_vector(
    values = values,
    markers = markers,
    n_markers = n_markers,
    name = source
  )
  if (any(!is.finite(values))) {
    stop(source, " must contain only finite values.")
  }

  if (identical(source, "maf")) {
    if (any(values < 0 | values > 1)) {
      stop("maf must contain values in [0, 1].")
    }
    maf_minor <- pmin(values, 1 - values)
    values <- 2 * maf_minor * (1 - maf_minor)
  }

  if (any(values <= 0 | values > 0.5)) {
    stop("h must contain heterozygosity values in (0, 0.5].")
  }

  pmax(values, .Machine$double.eps)
}

.stblr_selection_s_align_vector <- function(values, markers, n_markers, name) {
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
