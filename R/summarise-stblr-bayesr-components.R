#' Summarise BayesR Component Probabilities from an ST-BLR Fit
#'
#' Computes lightweight marker-level diagnostic summaries for formatted ST-BLR
#' BayesR fit objects. In these objects `fit$dm` is expected to be the non-null
#' posterior inclusion probability, `P(component > 0)`, and `fit$component_probabilities` is
#' expected to be a named list of marker-by-component probability matrices, one
#' per trait. When a `component_0` column is present it is treated as the null
#' component, so `fit$dm` should equal `1 - component_0`.
#'
#' This is a diagnostic summary helper, not a convergence diagnostic. The
#' `sum_pip` column is the posterior expected number of active markers for each
#' trait.
#'
#' @param fit Formatted ST-BLR BayesR fit object.
#' @param pip_thresholds Numeric vector of PIP thresholds used to count markers
#'   with `dm` greater than each threshold.
#' @param include_components Logical; include mean component probabilities and
#'   maximum non-null component probabilities from `fit$component_probabilities`.
#' @param include_chain_stability Logical; include chain-stability summaries
#'   from `fit$dm_chain_mean_sd` and `fit$bm_chain_mean_sd` when available.
#' @param top_unstable Integer. If positive and `fit$dm_chain_mean_sd` is available, return
#'   a list with the summary data frame and the top markers by PIP standard
#'   deviation for each trait.
#'
#' @return A data frame with one row per trait, unless `top_unstable > 0` and
#'   `fit$dm_chain_mean_sd` is available, in which case a list with elements `summary` and
#'   `unstable` is returned.
#'
#' @examples
#' comp <- matrix(
#'   c(0.8, 0.1, 0.1,
#'     0.2, 0.3, 0.5,
#'     0.6, 0.3, 0.1),
#'   nrow = 3,
#'   byrow = TRUE,
#'   dimnames = list(paste0("m", 1:3), paste0("component_", 0:2))
#' )
#' fit <- list(
#'   dm = matrix(1 - comp[, "component_0"], ncol = 1,
#'               dimnames = list(rownames(comp), "trait1")),
#'   bm = matrix(c(0.01, -0.03, 0.02), ncol = 1,
#'               dimnames = list(rownames(comp), "trait1")),
#'   component_probabilities = list(trait1 = comp)
#' )
#' summarise_components(fit)
#'
#' @export
summarise_components <- function(
    fit,
    pip_thresholds = c(0.001, 0.01, 0.05, 0.5, 0.95),
    include_components = TRUE,
    include_chain_stability = TRUE,
    top_unstable = 0L
) {
  if (!is.list(fit)) {
    stop("fit must be a list-like ST-BLR BayesR fit object.")
  }

  dm <- .stblr_bayesr_summary_as_matrix(fit$dm, "fit$dm", required = TRUE)
  bm <- .stblr_bayesr_summary_as_matrix(fit$bm, "fit$bm", required = FALSE)
  dm_chain_mean_sd <- .stblr_bayesr_summary_as_matrix(fit$dm_chain_mean_sd, "fit$dm_chain_mean_sd", required = FALSE)
  bm_chain_mean_sd <- .stblr_bayesr_summary_as_matrix(fit$bm_chain_mean_sd, "fit$bm_chain_mean_sd", required = FALSE)
  dm_component_mean <- .stblr_bayesr_summary_as_matrix(
    fit$dm_component_mean,
    "fit$dm_component_mean",
    required = FALSE
  )

  if (is.null(fit$component_probabilities)) {
    stop("fit$component_probabilities must be present.")
  }
  if (!is.list(fit$component_probabilities) || is.null(names(fit$component_probabilities)) ||
      any(!nzchar(names(fit$component_probabilities)))) {
    stop("fit$component_probabilities must be a named list.")
  }

  if (!is.numeric(pip_thresholds) || any(!is.finite(pip_thresholds))) {
    stop("pip_thresholds must be a finite numeric vector.")
  }
  pip_thresholds <- as.numeric(pip_thresholds)

  if (!is.logical(include_components) || length(include_components) != 1L ||
      is.na(include_components)) {
    stop("include_components must be TRUE or FALSE.")
  }
  if (!is.logical(include_chain_stability) ||
      length(include_chain_stability) != 1L ||
      is.na(include_chain_stability)) {
    stop("include_chain_stability must be TRUE or FALSE.")
  }
  if (!is.numeric(top_unstable) || length(top_unstable) != 1L ||
      !is.finite(top_unstable) || top_unstable < 0) {
    stop("top_unstable must be a non-negative integer scalar.")
  }
  top_unstable <- as.integer(top_unstable)

  trait_names <- colnames(dm)
  if (is.null(trait_names) || any(!nzchar(trait_names))) {
    trait_names <- paste0("trait", seq_len(ncol(dm)))
    colnames(dm) <- trait_names
  }

  .stblr_bayesr_summary_check_dims(bm, dm, "fit$bm")
  .stblr_bayesr_summary_check_dims(dm_chain_mean_sd, dm, "fit$dm_chain_mean_sd")
  .stblr_bayesr_summary_check_dims(bm_chain_mean_sd, dm, "fit$bm_chain_mean_sd")
  .stblr_bayesr_summary_check_dims(dm_component_mean, dm, "fit$dm_component_mean")

  rows <- lapply(seq_len(ncol(dm)), function(j) {
    trait <- trait_names[j]
    pip <- dm[, j]
    row <- data.frame(
      trait = trait,
      n_markers = nrow(dm),
      mean_pip = mean(pip, na.rm = TRUE),
      sum_pip = sum(pip, na.rm = TRUE),
      max_pip = max(pip, na.rm = TRUE),
      stringsAsFactors = FALSE
    )

    for (thr in pip_thresholds) {
      row[[paste0("n_pip_gt_", .stblr_bayesr_summary_threshold_name(thr))]] <-
        sum(pip > thr, na.rm = TRUE)
    }

    if (isTRUE(include_chain_stability) && !is.null(dm_chain_mean_sd)) {
      row$mean_pip_sd <- mean(dm_chain_mean_sd[, j], na.rm = TRUE)
      row$max_pip_sd <- max(dm_chain_mean_sd[, j], na.rm = TRUE)
    }
    if (!is.null(bm)) {
      abs_bm <- abs(bm[, j])
      row$mean_abs_effect <- mean(abs_bm, na.rm = TRUE)
      row$max_abs_effect <- max(abs_bm, na.rm = TRUE)
    }
    if (isTRUE(include_chain_stability) && !is.null(bm_chain_mean_sd)) {
      row$mean_effect_sd <- mean(bm_chain_mean_sd[, j], na.rm = TRUE)
      row$max_effect_sd <- max(bm_chain_mean_sd[, j], na.rm = TRUE)
    }
    if (!is.null(dm_component_mean)) {
      row$mean_component_index <- mean(dm_component_mean[, j], na.rm = TRUE)
      row$max_component_index <- max(dm_component_mean[, j], na.rm = TRUE)
    }

    cp <- fit$component_probabilities[[trait]]
    if (is.null(cp)) {
      stop("fit$component_probabilities is missing trait '", trait, "'.")
    }
    cp <- .stblr_bayesr_summary_as_matrix(
      cp,
      paste0("fit$component_probabilities[['", trait, "']]"),
      required = TRUE
    )
    if (nrow(cp) != nrow(dm)) {
      stop("fit$component_probabilities[['", trait, "']] must have nrow(fit$dm) rows.")
    }
    if (is.null(colnames(cp))) {
      colnames(cp) <- paste0("component_", seq.int(0L, ncol(cp) - 1L))
    }

    if ("component_0" %in% colnames(cp)) {
      component0 <- cp[, "component_0"]
      row$max_dm_component0_diff <- max(abs(pip - (1 - component0)), na.rm = TRUE)
    } else {
      row$max_dm_component0_diff <- NA_real_
    }

    if (isTRUE(include_components)) {
      component_means <- as.list(colMeans(cp, na.rm = TRUE))
      names(component_means) <- paste0("mean_", colnames(cp))
      row <- cbind(row, as.data.frame(component_means, check.names = FALSE))

      nonzero <- setdiff(colnames(cp), "component_0")
      if (length(nonzero) > 0L) {
        component_max <- as.list(apply(cp[, nonzero, drop = FALSE], 2L, max, na.rm = TRUE))
        names(component_max) <- paste0("max_", nonzero)
        row <- cbind(row, as.data.frame(component_max, check.names = FALSE))
      }
    }

    row
  })

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  class(out) <- c("stblr_bayesr_component_summary", class(out))

  if (top_unstable > 0L && !is.null(dm_chain_mean_sd)) {
    unstable <- lapply(seq_len(ncol(dm)), function(j) {
      trait <- trait_names[j]
      marker <- rownames(dm)
      if (is.null(marker)) marker <- paste0("marker", seq_len(nrow(dm)))
      ord <- order(dm_chain_mean_sd[, j], decreasing = TRUE, na.last = NA)
      ord <- ord[seq_len(min(top_unstable, length(ord)))]
      data.frame(
        marker = marker[ord],
        trait = trait,
        pip = dm[ord, j],
        pip_sd = dm_chain_mean_sd[ord, j],
        stringsAsFactors = FALSE
      )
    })
    names(unstable) <- trait_names
    return(list(summary = out, unstable = unstable))
  }

  out
}

.stblr_bayesr_summary_as_matrix <- function(x, name, required) {
  if (is.null(x)) {
    if (isTRUE(required)) stop(name, " must be present and matrix-like.")
    return(NULL)
  }
  if (is.matrix(x)) return(x)
  if (is.data.frame(x)) return(as.matrix(x))
  stop(name, " must be matrix-like.")
}

.stblr_bayesr_summary_check_dims <- function(x, dm, name) {
  if (is.null(x)) return(invisible(NULL))
  if (!identical(dim(x), dim(dm))) {
    stop(name, " must have the same dimensions as fit$dm.")
  }
  invisible(NULL)
}

.stblr_bayesr_summary_threshold_name <- function(x) {
  gsub("[^[:alnum:]]+", "_", as.character(x))
}
