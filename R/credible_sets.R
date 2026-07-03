#' Construct Credible Sets from Posterior Inclusion Probabilities and LD
#'
#' Constructs marker credible sets from posterior inclusion probabilities
#' (PIPs) and a dense regional LD/correlation matrix. This is an R-level
#' post-processing utility for fitted ST-BLR models, but it can also be used
#' directly with any named PIP vector and regional LD matrix.
#'
#' Dense LD is intended for regional credible sets. Markers are selected by PIP
#' or by an LD-smoothed PIP score, and each credible set contains remaining
#' markers in LD with the lead marker until the requested cumulative PIP
#' coverage is reached or no further candidate marker is available.
#'
#' @param pip Named numeric vector of posterior inclusion probabilities.
#' @param LD Dense square LD/correlation matrix. If marker names are available
#'   in both `pip` and `LD`, inputs are aligned by common marker names.
#' @param coverage Credible-set target cumulative PIP, in `(0, 1)`.
#' @param min_r2 Minimum squared correlation to the lead marker for set
#'   membership.
#' @param pip_cutoff PIPs below this value are set to zero for set construction.
#' @param method Lead-marker selection method. `"pip"` uses the highest
#'   remaining PIP; `"ld_pip"` uses an LD-smoothed PIP score.
#' @param max_sets Maximum number of credible sets to construct.
#' @param allow_incomplete Logical; if `TRUE`, return sets whose cumulative PIP
#'   does not reach `coverage` when all available LD-neighborhood candidates
#'   have been used. The default `FALSE` reports only credible sets with
#'   `cs_pip >= coverage`.
#' @param remove Marker removal rule after a credible set is constructed.
#'   `"ld_neighborhood"` suppresses all candidate markers in the lead marker's
#'   LD neighborhood from future lead selection. `"credible_set"` suppresses
#'   only markers included in the selected minimal credible set.
#'
#' @return A list with `summary`, `sets`, aligned `pip`, and `parameters`.
#' @export
make_credible_sets_from_ld <- function(
    pip,
    LD,
    coverage = 0.95,
    min_r2 = 0.5,
    pip_cutoff = 0.001,
    method = c("pip", "ld_pip"),
    max_sets = Inf,
    allow_incomplete = FALSE,
    remove = c("ld_neighborhood", "credible_set")
) {
  method <- match.arg(method)
  remove <- match.arg(remove)
  .check_scalar_range(coverage, "coverage", lower = 0, upper = 1,
                      lower_open = TRUE, upper_open = TRUE)
  .check_scalar_range(min_r2, "min_r2", lower = 0, upper = 1)
  .check_scalar_range(pip_cutoff, "pip_cutoff", lower = 0, upper = 1)
  if (!is.logical(allow_incomplete) || length(allow_incomplete) != 1L ||
      is.na(allow_incomplete)) {
    stop("allow_incomplete must be TRUE or FALSE.")
  }
  if (!is.numeric(max_sets) || length(max_sets) != 1L ||
      is.na(max_sets) || max_sets < 0) {
    stop("max_sets must be a non-negative scalar.")
  }

  pip_names <- names(pip)
  pip <- as.numeric(pip)
  names(pip) <- pip_names

  if (!is.matrix(LD)) LD <- as.matrix(LD)
  storage.mode(LD) <- "double"
  if (length(dim(LD)) != 2L || nrow(LD) != ncol(LD)) {
    stop("LD must be a dense square matrix.")
  }

  ld_names <- rownames(LD)
  if (is.null(ld_names)) ld_names <- colnames(LD)

  if (!is.null(pip_names) && !is.null(ld_names)) {
    common <- intersect(pip_names, ld_names)
    if (length(common) < 1L) {
      stop("pip and LD have names but no markers in common.")
    }
    pip <- pip[common]
    LD <- LD[common, common, drop = FALSE]
  } else {
    if (length(pip) != nrow(LD)) {
      stop("LD must have the same number of rows/columns as pip.")
    }
    if (is.null(pip_names) && !is.null(ld_names)) {
      names(pip) <- ld_names
    }
  }

  m <- length(pip)
  if (m < 1L) stop("pip must contain at least one marker.")
  if (nrow(LD) != m) stop("LD must be the same length as pip after matching.")

  marker <- names(pip)
  if (is.null(marker)) marker <- paste0("marker", seq_len(m))
  names(pip) <- marker
  rownames(LD) <- colnames(LD) <- marker

  pip[!is.finite(pip)] <- 0
  pip[pip < pip_cutoff] <- 0

  r2 <- LD^2
  r2[!is.finite(r2)] <- 0
  diag(r2) <- 1

  remaining <- pip > 0
  summaries <- list()
  sets <- list()
  cs_index <- 0L

  while (any(remaining) && cs_index < max_sets) {
    rem_idx <- which(remaining)
    if (method == "pip") {
      lead <- rem_idx[which.max(pip[rem_idx])]
    } else {
      score <- as.numeric(r2 %*% pip)
      score[!remaining] <- -Inf
      lead <- which.max(score)
    }

    candidate_idx <- rem_idx[r2[rem_idx, lead] >= min_r2]
    if (length(candidate_idx) < 1L) candidate_idx <- lead
    candidate_idx <- candidate_idx[order(pip[candidate_idx], decreasing = TRUE)]

    cumulative <- cumsum(pip[candidate_idx])
    take_n <- which(cumulative >= coverage)[1L]
    if (is.na(take_n)) {
      if (!allow_incomplete) {
        remaining[candidate_idx] <- FALSE
        next
      }
      take_n <- length(candidate_idx)
    }
    selected <- candidate_idx[seq_len(take_n)]

    cs_index <- cs_index + 1L
    cs_name <- paste0("CS", cs_index)
    cs_pip <- sum(pip[selected])
    r2_to_lead <- r2[selected, lead]

    sets[[cs_name]] <- data.frame(
      marker = marker[selected],
      pip = pip[selected],
      r2_to_lead = r2_to_lead,
      rank = seq_along(selected),
      stringsAsFactors = FALSE
    )

    summaries[[cs_name]] <- data.frame(
      cs = cs_name,
      lead_marker = marker[lead],
      lead_pip = pip[lead],
      cs_pip = cs_pip,
      n_markers = length(selected),
      min_r2_to_lead = min(r2_to_lead),
      mean_r2_to_lead = mean(r2_to_lead),
      coverage = coverage,
      min_r2 = min_r2,
      pip_cutoff = pip_cutoff,
      method = method,
      stringsAsFactors = FALSE
    )

    if (remove == "ld_neighborhood") {
      remaining[candidate_idx] <- FALSE
    } else {
      remaining[selected] <- FALSE
    }
  }

  summary <- if (length(summaries) > 0L) {
    do.call(rbind, summaries)
  } else {
    data.frame(
      cs = character(), lead_marker = character(), lead_pip = numeric(),
      cs_pip = numeric(), n_markers = integer(), min_r2_to_lead = numeric(),
      mean_r2_to_lead = numeric(), coverage = numeric(), min_r2 = numeric(),
      pip_cutoff = numeric(), method = character(), stringsAsFactors = FALSE
    )
  }
  rownames(summary) <- NULL

  list(
    summary = summary,
    sets = sets,
    pip = pip,
    parameters = list(
      coverage = coverage,
      min_r2 = min_r2,
      pip_cutoff = pip_cutoff,
      method = method,
      max_sets = max_sets,
      allow_incomplete = allow_incomplete,
      remove = remove
    )
  )
}

#' Construct Approximate Multi-Signal Credible Sets from PIPs and LD
#'
#' Constructs multiple credible sets from marginal posterior inclusion
#' probabilities (PIPs) and a dense regional LD/correlation matrix. This is an
#' approximate R-level post-processing tool for local fine-mapping output, most
#' useful when the regional total PIP is much larger than the lead-marker PIP.
#'
#' The procedure repeatedly selects a lead marker from the remaining marginal
#' PIPs, builds a local LD neighborhood around that lead, and accumulates
#' markers by decreasing PIP until the requested `coverage` is reached. These
#' sets are based on marginal PIPs, not posterior inclusion configurations, and
#' are therefore not identical to SuSiE-style per-effect credible sets. Rigorous
#' per-effect credible sets require posterior configurations or a
#' component-specific model.
#'
#' @param pip Named numeric vector of posterior inclusion probabilities.
#' @param LD Dense square LD/correlation matrix. If marker names are available
#'   in both `pip` and `LD`, inputs are aligned by common marker names.
#' @param coverage Credible-set target cumulative PIP, in `(0, 1)`.
#' @param min_r2 Minimum squared correlation to the lead marker for candidate
#'   membership.
#' @param pip_cutoff PIPs below this value are set to zero for set construction.
#' @param signal_coverage Minimum total remaining PIP used to continue searching
#'   for complete additional signals when `allow_incomplete = FALSE`.
#' @param min_signal_pip Minimum remaining lead-marker PIP required to start a
#'   new signal.
#' @param max_signals Maximum number of signals/credible sets to return.
#' @param method Lead-marker selection method. `"pip"` uses the highest
#'   remaining PIP; `"ld_pip"` uses an LD-smoothed PIP score.
#' @param remove Marker removal rule after a signal is constructed.
#'   `"credible_set"` suppresses only selected credible-set markers from future
#'   lead selection. `"ld_neighborhood"` suppresses all candidate markers in the
#'   lead marker's LD neighborhood.
#' @param allow_incomplete Logical; if `TRUE`, return sets whose cumulative PIP
#'   does not reach `coverage` when all available LD-neighborhood candidates
#'   have been used.
#'
#' @return A list with `summary`, `sets`, aligned `pip`, and `parameters`.
#' @export
make_multisignal_credible_sets_from_ld <- function(
    pip,
    LD,
    coverage = 0.95,
    min_r2 = 0.5,
    pip_cutoff = 0.001,
    signal_coverage = 0.95,
    min_signal_pip = 0.05,
    max_signals = Inf,
    method = c("pip", "ld_pip"),
    remove = c("credible_set", "ld_neighborhood"),
    allow_incomplete = FALSE
) {
  method <- match.arg(method)
  remove <- match.arg(remove)
  .check_scalar_range(coverage, "coverage", lower = 0, upper = 1,
                      lower_open = TRUE, upper_open = TRUE)
  .check_scalar_range(min_r2, "min_r2", lower = 0, upper = 1)
  .check_scalar_range(pip_cutoff, "pip_cutoff", lower = 0, upper = 1)
  .check_scalar_range(signal_coverage, "signal_coverage", lower = 0, upper = Inf)
  .check_scalar_range(min_signal_pip, "min_signal_pip", lower = 0, upper = 1)
  if (!is.numeric(max_signals) || length(max_signals) != 1L ||
      is.na(max_signals) || max_signals < 0) {
    stop("max_signals must be a non-negative scalar.")
  }
  if (!is.logical(allow_incomplete) || length(allow_incomplete) != 1L ||
      is.na(allow_incomplete)) {
    stop("allow_incomplete must be TRUE or FALSE.")
  }

  aligned <- .stblr_align_pip_ld(pip, LD)
  pip <- aligned$pip
  LD <- aligned$LD
  marker <- names(pip)

  pip[!is.finite(pip)] <- 0
  pip[pip < pip_cutoff] <- 0

  r2 <- LD^2
  r2[!is.finite(r2)] <- 0
  diag(r2) <- 1

  remaining <- pip > 0
  summaries <- list()
  sets <- list()
  signal_index <- 0L

  while (any(remaining) && signal_index < max_signals) {
    rem_idx <- which(remaining)
    remaining_pip_before <- sum(pip[rem_idx])
    if (!allow_incomplete && remaining_pip_before < signal_coverage) break

    if (method == "pip") {
      lead <- rem_idx[which.max(pip[rem_idx])]
    } else {
      score <- as.numeric(r2 %*% (pip * remaining))
      score[!remaining] <- -Inf
      lead <- which.max(score)
    }
    if (!is.finite(pip[lead]) || pip[lead] < min_signal_pip) break

    candidate_idx <- rem_idx[r2[rem_idx, lead] >= min_r2]
    if (length(candidate_idx) < 1L) candidate_idx <- lead
    candidate_idx <- candidate_idx[order(pip[candidate_idx], decreasing = TRUE)]

    cumulative <- cumsum(pip[candidate_idx])
    take_n <- which(cumulative >= coverage)[1L]
    complete <- !is.na(take_n)
    if (!complete) {
      if (!allow_incomplete) {
        remaining[lead] <- FALSE
        next
      }
      take_n <- length(candidate_idx)
    }
    selected <- candidate_idx[seq_len(take_n)]

    signal_index <- signal_index + 1L
    cs_name <- paste0("CS", signal_index)
    cs_pip <- sum(pip[selected])
    r2_to_lead <- r2[selected, lead]

    sets[[cs_name]] <- data.frame(
      signal = signal_index,
      marker = marker[selected],
      pip = pip[selected],
      r2_to_lead = r2_to_lead,
      rank = seq_along(selected),
      stringsAsFactors = FALSE
    )

    remove_idx <- if (remove == "ld_neighborhood") candidate_idx else selected
    remaining[remove_idx] <- FALSE
    remaining_pip_after <- sum(pip[remaining])

    summaries[[cs_name]] <- data.frame(
      signal = signal_index,
      cs = cs_name,
      lead_marker = marker[lead],
      lead_pip = pip[lead],
      cs_pip = cs_pip,
      remaining_pip_before = remaining_pip_before,
      remaining_pip_after = remaining_pip_after,
      n_markers = length(selected),
      min_r2_to_lead = min(r2_to_lead),
      mean_r2_to_lead = mean(r2_to_lead),
      coverage = coverage,
      min_r2 = min_r2,
      pip_cutoff = pip_cutoff,
      min_signal_pip = min_signal_pip,
      method = method,
      remove = remove,
      complete = complete,
      stringsAsFactors = FALSE
    )
  }

  summary <- if (length(summaries) > 0L) {
    do.call(rbind, summaries)
  } else {
    data.frame(
      signal = integer(), cs = character(), lead_marker = character(),
      lead_pip = numeric(), cs_pip = numeric(),
      remaining_pip_before = numeric(), remaining_pip_after = numeric(),
      n_markers = integer(), min_r2_to_lead = numeric(),
      mean_r2_to_lead = numeric(), coverage = numeric(), min_r2 = numeric(),
      pip_cutoff = numeric(), min_signal_pip = numeric(), method = character(),
      remove = character(), complete = logical(), stringsAsFactors = FALSE
    )
  }
  rownames(summary) <- NULL

  list(
    summary = summary,
    sets = sets,
    pip = pip,
    parameters = list(
      coverage = coverage,
      min_r2 = min_r2,
      pip_cutoff = pip_cutoff,
      signal_coverage = signal_coverage,
      min_signal_pip = min_signal_pip,
      max_signals = max_signals,
      method = method,
      remove = remove,
      allow_incomplete = allow_incomplete
    )
  )
}

#' Construct Credible Sets for Fitted ST-BLR Models
#'
#' Constructs credible sets from marker-level PIPs in `fit$dm` for fitted
#' ST-BLR models. This is an R-level post-processing function for fits from
#' `stblr_csr()` and `stblr_bed_marker()`.
#'
#' Credible sets can be built from predefined marker sets or from automatically
#' defined loci around high-PIP markers. Dense LD can be supplied directly for
#' regional credible sets. Alternatively, for CSR workflows, `Glist$sparseLD`
#' can be used and regional dense LD matrices are extracted internally from the
#' disk-backed sparse CSR.
#'
#' When using sparse LD from `Glist$sparseLD`, the sparse LD should have been
#' computed with a sufficiently permissive `r2_threshold` and local window to
#' support the requested `min_r2`; otherwise, LD pairs needed for set
#' construction may be absent.
#'
#' @param fit Fitted ST-BLR object containing marker PIPs in `fit$dm`.
#' @param Glist Optional genotype list containing marker map and `sparseLD`.
#' @param LD Optional dense LD matrix, or a named list of dense LD matrices for
#'   multiple predefined sets.
#' @param sets Optional named or unnamed list of marker vectors defining
#'   regions to analyze.
#' @param trait Trait column name or index when `fit$dm` is a matrix/data frame.
#' @param coverage,min_r2,pip_cutoff Credible-set construction parameters
#'   passed to [make_credible_sets_from_ld()].
#' @param locus_pip_cutoff Minimum PIP for automatic lead-locus selection.
#' @param max_locus_distance Locus half-window in base pairs for automatic
#'   locus definition.
#' @param max_loci Maximum number of automatic loci.
#' @param method Lead-marker selection method passed to
#'   [make_credible_sets_from_ld()].
#' @param allow_incomplete Logical; if `TRUE`, return credible-set candidates
#'   whose cumulative PIP does not reach `coverage`. The default `FALSE`
#'   reports only valid credible sets with `cs_pip >= coverage`.
#' @param remove Marker removal rule passed to [make_credible_sets_from_ld()].
#'   The default `"ld_neighborhood"` removes the whole lead LD neighborhood from
#'   future credible-set selection; `"credible_set"` removes only selected
#'   credible-set markers.
#'
#' @return A list with combined `summary`, nested credible `sets`, `loci`,
#'   `parameters`, and resolved `trait`.
#' @export
make_credible_sets <- function(
    fit,
    Glist = NULL,
    LD = NULL,
    sets = NULL,
    trait = 1,
    coverage = 0.95,
    min_r2 = 0.5,
    pip_cutoff = 0.001,
    locus_pip_cutoff = 0.01,
    max_locus_distance = 1e6,
    max_loci = Inf,
    method = c("pip", "ld_pip"),
    allow_incomplete = FALSE,
    remove = c("ld_neighborhood", "credible_set")
) {
  method <- match.arg(method)
  remove <- match.arg(remove)
  
  pip_info <- .stblr_extract_pip(fit, trait = trait)
  pip <- pip_info$pip
  
  if (is.null(sets)) {
    map <- .stblr_marker_map_from_Glist(Glist, fit = fit)
    
    loci <- .stblr_define_loci(
      pip = pip,
      map = map,
      locus_pip_cutoff = locus_pip_cutoff,
      max_locus_distance = max_locus_distance,
      max_loci = max_loci
    )
    
    if (nrow(loci) < 1L) {
      stop("No automatic loci found. Lower locus_pip_cutoff or supply sets.")
    }
    
    sets <- stats::setNames(loci$markers, loci$locus)
  } else {
    sets <- .stblr_clean_marker_sets(sets, pip)
    loci <- .stblr_loci_from_sets(sets, pip, Glist = Glist, fit = fit)
  }
  
  locus_sets <- sets
  
  use_sparse <- is.null(LD)
  csr <- NULL
  
  if (use_sparse) {
    if (is.null(Glist) || is.null(Glist$sparseLD$prefix) ||
        !nzchar(Glist$sparseLD$prefix)) {
      stop("Supply LD or Glist$sparseLD$prefix.")
    }
    
    if (!is.null(Glist$sparseLD$r2_threshold) &&
        is.finite(Glist$sparseLD$r2_threshold) &&
        min_r2 < Glist$sparseLD$r2_threshold) {
      warning(
        "Glist$sparseLD was built with r2_threshold = ",
        Glist$sparseLD$r2_threshold,
        ", which is higher than min_r2 = ", min_r2,
        "; some LD pairs needed for credible-set construction may be missing."
      )
    }
    
    csr <- sparseLD_read_CSR(Glist$sparseLD$prefix, one_based = TRUE)
  }
  
  summaries <- list()
  out_sets <- list()
  
  for (locus_name in names(sets)) {
    markers <- sets[[locus_name]]
    regional_pip <- pip[markers]
    
    regional_LD <- if (use_sparse) {
      locus_rows <- loci[loci$locus == locus_name, , drop = FALSE]
      idx <- locus_rows$indices[[1L]]
      
      .extract_sparseLD_region_dense(
        csr,
        idx,
        marker_names = markers
      )
    } else {
      .stblr_resolve_regional_ld(
        LD,
        locus_name,
        markers,
        length(sets)
      )
    }
    
    cs <- make_credible_sets_from_ld(
      pip = regional_pip,
      LD = regional_LD,
      coverage = coverage,
      min_r2 = min_r2,
      pip_cutoff = pip_cutoff,
      method = method,
      allow_incomplete = allow_incomplete,
      remove = remove
    )
    
    out_sets[[locus_name]] <- cs$sets
    
    if (nrow(cs$summary) > 0L) {
      loc <- loci[loci$locus == locus_name, , drop = FALSE]
      
      tab <- cbind(
        data.frame(
          locus = locus_name,
          trait = pip_info$trait,
          chr = loc$chr,
          start = loc$start,
          end = loc$end,
          stringsAsFactors = FALSE
        ),
        cs$summary,
        stringsAsFactors = FALSE
      )
      
      summaries[[locus_name]] <- tab
    }
  }
  
  summary <- if (length(summaries) > 0L) {
    do.call(rbind, summaries)
  } else {
    data.frame()
  }
  
  rownames(summary) <- NULL
  
  list(
    summary = summary,
    
    # Credible-set marker lists.
    sets = out_sets,
    
    # All markers in each detected/supplied locus.
    # Useful as direct input to finemap_stblr_csr(..., sets = ...).
    locus_sets = locus_sets,
    
    loci = loci[, setdiff(names(loci), c("markers", "indices")), drop = FALSE],
    
    parameters = list(
      coverage = coverage,
      min_r2 = min_r2,
      pip_cutoff = pip_cutoff,
      locus_pip_cutoff = locus_pip_cutoff,
      max_locus_distance = max_locus_distance,
      max_loci = max_loci,
      method = method,
      allow_incomplete = allow_incomplete,
      remove = remove
    ),
    
    trait = pip_info$trait
  )
}

.check_scalar_range <- function(x, name, lower, upper,
                                lower_open = FALSE, upper_open = FALSE) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x)) {
    stop(name, " must be a finite numeric scalar.")
  }
  low_ok <- if (lower_open) x > lower else x >= lower
  up_ok <- if (upper_open) x < upper else x <= upper
  if (!low_ok || !up_ok) {
    bracket <- paste0(if (lower_open) "(" else "[", lower, ", ", upper,
                      if (upper_open) ")" else "]")
    stop(name, " must be in ", bracket, ".")
  }
  invisible(TRUE)
}

.stblr_align_pip_ld <- function(pip, LD) {
  pip_names <- names(pip)
  pip <- as.numeric(pip)
  names(pip) <- pip_names

  if (!is.matrix(LD)) LD <- as.matrix(LD)
  storage.mode(LD) <- "double"
  if (length(dim(LD)) != 2L || nrow(LD) != ncol(LD)) {
    stop("LD must be a dense square matrix.")
  }

  ld_names <- rownames(LD)
  if (is.null(ld_names)) ld_names <- colnames(LD)

  if (!is.null(pip_names) && !is.null(ld_names)) {
    common <- intersect(pip_names, ld_names)
    if (length(common) < 1L) {
      stop("pip and LD have names but no markers in common.")
    }
    pip <- pip[common]
    LD <- LD[common, common, drop = FALSE]
  } else {
    if (length(pip) != nrow(LD)) {
      stop("LD must have the same number of rows/columns as pip.")
    }
    if (is.null(pip_names) && !is.null(ld_names)) {
      names(pip) <- ld_names
    }
  }

  m <- length(pip)
  if (m < 1L) stop("pip must contain at least one marker.")
  if (nrow(LD) != m) stop("LD must be the same length as pip after matching.")

  marker <- names(pip)
  if (is.null(marker)) marker <- paste0("marker", seq_len(m))
  names(pip) <- marker
  rownames(LD) <- colnames(LD) <- marker

  list(pip = pip, LD = LD)
}

.stblr_extract_pip <- function(fit, trait = 1) {
  if (is.null(fit$dm)) stop("fit must contain marker PIPs in fit$dm.")
  dm <- fit$dm

  if (is.matrix(dm) || is.data.frame(dm)) {
    dm <- as.matrix(dm)
    if (is.character(trait)) {
      if (is.null(colnames(dm)) || !(trait %in% colnames(dm))) {
        stop("trait was not found in colnames(fit$dm).")
      }
      trait_index <- match(trait, colnames(dm))
      trait_name <- trait
    } else {
      trait_index <- as.integer(trait)
      if (length(trait_index) != 1L || is.na(trait_index) ||
          trait_index < 1L || trait_index > ncol(dm)) {
        stop("trait must be a valid fit$dm column index.")
      }
      trait_name <- colnames(dm)[trait_index]
      if (is.null(trait_name) || is.na(trait_name)) trait_name <- trait_index
    }
    pip <- dm[, trait_index]
    names(pip) <- rownames(dm)
  } else {
    pip <- as.numeric(dm)
    names(pip) <- names(dm)
    trait_name <- if (is.character(trait)) trait else 1L
  }

  if (is.null(names(pip))) {
    names(pip) <- paste0("marker", seq_along(pip))
  }

  list(pip = pip, trait = trait_name)
}

.stblr_clean_marker_sets <- function(sets, pip) {
  if (!is.list(sets)) sets <- list(sets)
  if (is.null(names(sets)) || any(!nzchar(names(sets)))) {
    names(sets) <- paste0("locus", seq_along(sets))
  }

  out <- list()
  dropped <- character()
  for (nm in names(sets)) {
    markers <- unique(as.character(sets[[nm]]))
    markers <- markers[markers %in% names(pip)]
    if (length(markers) < 1L) {
      dropped <- c(dropped, nm)
    } else {
      out[[nm]] <- markers
    }
  }
  if (length(dropped) > 0L) {
    warning("Dropping empty marker sets: ", paste(dropped, collapse = ", "))
  }
  if (length(out) < 1L) stop("No supplied sets contain markers in fit$dm.")
  out
}

.stblr_marker_map_from_Glist <- function(Glist, fit = NULL) {
  if (is.null(Glist)) stop("Glist is required to define marker positions.")
  if (is.null(Glist$rsids)) stop("Glist$rsids is required for marker names.")
  if (is.null(Glist$pos)) stop("Glist$pos is required for marker positions.")

  if (!is.null(Glist$sparseLD$chr) && !is.null(Glist$sparseLD$cls)) {
    component_idx <- Glist$sparseLD$chr
    if (!is.list(component_idx)) component_idx <- as.list(component_idx)
    cls <- Glist$sparseLD$cls
    if (!is.list(cls)) cls <- as.list(cls)
    if (length(component_idx) != length(cls)) {
      stop("Glist$sparseLD$chr and Glist$sparseLD$cls must have the same length.")
    }
  } else {
    component_idx <- as.list(seq_along(Glist$rsids))
    cls <- lapply(Glist$rsids, seq_along)
  }

  rows <- list()
  sparse_index <- 0L
  for (k in seq_along(component_idx)) {
    cc <- as.integer(component_idx[[k]])
    if (length(cc) != 1L || is.na(cc)) {
      stop("Each Glist$sparseLD$chr entry must identify one Glist component.")
    }
    cl <- as.integer(cls[[k]])
    if (length(cl) < 1L) next
    if (is.null(Glist$rsids[[cc]]) || is.null(Glist$pos[[cc]])) {
      stop("Glist$rsids and Glist$pos must contain chromosome/file index ", cc, ".")
    }
    chr <- if (!is.null(Glist$chr) && !is.null(Glist$chr[[cc]])) {
      Glist$chr[[cc]][cl]
    } else {
      rep(cc, length(cl))
    }
    sparse_idx <- seq.int(sparse_index + 1L, sparse_index + length(cl))
    sparse_index <- sparse_index + length(cl)
    rows[[length(rows) + 1L]] <- data.frame(
      marker = as.character(Glist$rsids[[cc]][cl]),
      chr = chr,
      pos = as.numeric(Glist$pos[[cc]][cl]),
      index = sparse_idx,
      stringsAsFactors = FALSE
    )
  }

  map <- do.call(rbind, rows)
  rownames(map) <- NULL

  pip_names <- NULL
  if (!is.null(fit) && !is.null(fit$dm)) {
    if (is.matrix(fit$dm) || is.data.frame(fit$dm)) {
      pip_names <- rownames(fit$dm)
    } else {
      pip_names <- names(fit$dm)
    }
  }

  if (!is.null(pip_names)) {
    keep <- match(pip_names, map$marker)
    keep <- keep[!is.na(keep)]
    map <- map[keep, , drop = FALSE]
    rownames(map) <- NULL
  }

  map
}

.stblr_define_loci <- function(pip, map, locus_pip_cutoff,
                               max_locus_distance, max_loci) {
  .check_scalar_range(locus_pip_cutoff, "locus_pip_cutoff", 0, 1)
  if (!is.numeric(max_locus_distance) || length(max_locus_distance) != 1L ||
      !is.finite(max_locus_distance) || max_locus_distance < 0) {
    stop("max_locus_distance must be a finite non-negative scalar.")
  }
  if (!is.numeric(max_loci) || length(max_loci) != 1L ||
      is.na(max_loci) || max_loci < 0) {
    stop("max_loci must be a non-negative scalar.")
  }

  map <- map[match(names(pip), map$marker), , drop = FALSE]
  map <- map[!is.na(map$marker), , drop = FALSE]
  pip <- pip[map$marker]
  pip[!is.finite(pip)] <- 0

  remaining <- rep(TRUE, length(pip))
  loci <- list()
  locus_index <- 0L

  while (any(remaining) && locus_index < max_loci) {
    rem_idx <- which(remaining)
    lead_rel <- rem_idx[which.max(pip[rem_idx])]
    if (pip[lead_rel] < locus_pip_cutoff) break

    same_chr <- map$chr == map$chr[lead_rel]
    in_window <- abs(map$pos - map$pos[lead_rel]) <= max_locus_distance
    idx <- which(same_chr & in_window)

    locus_index <- locus_index + 1L
    locus_name <- paste0("locus", locus_index)
    loci[[locus_name]] <- data.frame(
      locus = locus_name,
      chr = map$chr[lead_rel],
      start = min(map$pos[idx], na.rm = TRUE),
      end = max(map$pos[idx], na.rm = TRUE),
      lead_marker = map$marker[lead_rel],
      lead_pip = pip[lead_rel],
      n_markers = length(idx),
      stringsAsFactors = FALSE
    )
    loci[[locus_name]]$markers <- list(map$marker[idx])
    loci[[locus_name]]$indices <- list(map$index[idx])

    remaining[idx] <- FALSE
  }

  if (length(loci) < 1L) {
    return(data.frame(
      locus = character(), chr = integer(), start = numeric(), end = numeric(),
      lead_marker = character(), lead_pip = numeric(), n_markers = integer(),
      stringsAsFactors = FALSE
    ))
  }
  out <- do.call(rbind, loci)
  rownames(out) <- NULL
  out
}

.stblr_loci_from_sets <- function(sets, pip, Glist = NULL, fit = NULL) {
  map <- tryCatch(
    .stblr_marker_map_from_Glist(Glist, fit = fit),
    error = function(e) NULL
  )

  rows <- lapply(names(sets), function(nm) {
    markers <- sets[[nm]]
    lead <- markers[which.max(pip[markers])]
    loc_map <- if (!is.null(map)) map[match(markers, map$marker), , drop = FALSE] else NULL
    loc_map <- if (!is.null(loc_map)) loc_map[!is.na(loc_map$marker), , drop = FALSE] else NULL

    chr <- start <- end <- NA
    indices <- match(markers, names(pip))
    if (!is.null(loc_map) && nrow(loc_map) > 0L) {
      chr <- if (length(unique(loc_map$chr)) == 1L) loc_map$chr[1L] else NA
      start <- min(loc_map$pos, na.rm = TRUE)
      end <- max(loc_map$pos, na.rm = TRUE)
      indices <- loc_map$index[match(markers, loc_map$marker)]
    }

    row <- data.frame(
      locus = nm,
      chr = chr,
      start = start,
      end = end,
      lead_marker = lead,
      lead_pip = pip[lead],
      n_markers = length(markers),
      stringsAsFactors = FALSE
    )
    .add_list_columns(row, markers = markers, indices = indices)
  })

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

.add_list_columns <- function(x, markers, indices) {
  x$markers <- list(markers)
  x$indices <- list(indices)
  x
}

.stblr_resolve_regional_ld <- function(LD, locus_name, markers, n_sets) {
  if (is.list(LD) && !is.matrix(LD)) {
    if (!locus_name %in% names(LD)) {
      stop("LD is a list but has no entry named ", locus_name, ".")
    }
    LD <- LD[[locus_name]]
  } else if (n_sets > 1L && is.null(rownames(LD))) {
    stop("A single LD matrix for multiple sets must have marker names.")
  }

  if (!is.matrix(LD)) LD <- as.matrix(LD)
  if (!is.null(rownames(LD))) {
    missing <- setdiff(markers, rownames(LD))
    if (length(missing) > 0L) {
      stop("LD is missing markers for ", locus_name, ": ",
           paste(missing, collapse = ", "))
    }
    LD <- LD[markers, markers, drop = FALSE]
  }
  LD
}

.extract_sparseLD_region_dense <- function(csr, idx, marker_names = NULL) {
  ptr <- .csr_component(csr, c("ptr", "row_ptr", "indptr"))
  col_idx <- .csr_component(csr, c("idx", "col_idx", "indices"))
  values <- .csr_component(csr, c("x", "values", "val"))

  ptr <- as.numeric(ptr)
  col_idx <- as.integer(col_idx)
  values <- as.numeric(values)
  idx <- as.integer(idx)

  if (length(idx) < 1L) stop("idx must contain at least one marker index.")
  if (length(col_idx) != length(values)) {
    stop("CSR column indices and values must have the same length.")
  }
  if (anyNA(idx) || any(idx < 1L)) {
    stop("idx must contain one-based positive marker indices.")
  }

  m <- length(idx)
  LD <- diag(1, m)
  if (!is.null(marker_names)) rownames(LD) <- colnames(LD) <- marker_names

  local <- seq_len(m)
  names(local) <- as.character(idx)

  # Row pointers are CSR offsets into col_idx/values and remain zero-based
  # even when sparseLD_read_CSR(one_based = TRUE) returns one-based columns.
  for (ii in seq_along(idx)) {
    row <- idx[ii]
    if (row + 1L > length(ptr)) next
    start <- ptr[row] + 1
    end <- ptr[row + 1L]
    if (!is.finite(start) || !is.finite(end) || end < start) next
    pos <- seq.int(start, end)
    cols <- col_idx[pos]
    vals <- values[pos]
    keep <- as.character(cols) %in% names(local)
    if (!any(keep)) next
    jj <- unname(local[as.character(cols[keep])])
    LD[ii, jj] <- vals[keep]
    LD[jj, ii] <- vals[keep]
  }

  diag(LD) <- 1
  if (m > 10L && all(LD[row(LD) != col(LD)] == 0)) {
    warning(
      "Regional LD matrix has more than 10 markers but all off-diagonal ",
      "entries are zero; LD extraction or sparseLD indexing may be wrong."
    )
  }
  LD
}

.csr_component <- function(csr, candidates) {
  hit <- intersect(candidates, names(csr))
  if (length(hit) < 1L) {
    stop("CSR object is missing component: ", paste(candidates, collapse = " / "))
  }
  csr[[hit[1L]]]
}
