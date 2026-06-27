#' Local Fine-Mapping for CSR ST-BLR Fits
#'
#' Fine-map predefined local marker regions using the global architecture
#' learned from a previous [stblr_csr()] fit. This is a local fine-mapping
#' step, not a genome-wide refit: every marker supplied in each `sets` element
#' is included in that local model. `coverage`, `min_r2`, and `pip_cutoff` are
#' used only for optional credible-set construction after local fitting.
#'
#' The local CSR sampler is run with fixed global parameters `ve`, `vb`, and
#' `pi` when available. By default these are resolved from the global fit
#' posterior traces after burn-in; user-supplied values override the fit. The
#' implementation calls the lower-level CSR sampler with `updateE = FALSE`,
#' `updateB = FALSE`, and `updatePi = FALSE`, so these parameters are fixed for
#' the local runs. `use_residual = TRUE` residualizes regional summary
#' statistics against posterior mean global effects outside the region.
#'
#' @param fit Global `stblr_csr()` fit.
#' @param Glist Genotype/sparse LD object containing `Glist$sparseLD`.
#' @param stats Summary-statistics object used by `stblr_csr()`.
#' @param sets Named or unnamed list of marker vectors defining local regions.
#' @param trait Trait name or column index. If `NULL`, all traits in `fit$dm`
#'   are fine-mapped.
#' @param nruns Number of independent local runs per region.
#' @param use_residual Logical; if `TRUE`, residualize local `wy` statistics
#'   against global posterior mean effects outside the region.
#' @param nit,nburn Local MCMC iteration controls.
#' @param seeds Optional integer vector of length `nruns`. If `NULL`, seeds are
#'   generated deterministically from `fit$input$seed` when available, otherwise
#'   from `seq_len(nruns)`.
#' @param ve,vb,pi Optional scalar or named numeric global-parameter overrides.
#'   `pi` is the marker inclusion probability, not `c(1 - pi, pi)`.
#' @param credible_sets Logical; construct credible sets from aggregated local
#'   PIPs.
#' @param coverage,min_r2,pip_cutoff Credible-set construction parameters.
#' @param cs_mode Credible-set mode. `"single"` keeps the historical
#'   single-signal behavior. `"multi"` uses approximate marginal-PIP
#'   multi-signal post-processing via
#'   [make_multisignal_credible_sets_from_ld()].
#' @param min_signal_pip Minimum remaining lead-marker PIP required to start a
#'   new signal when `cs_mode = "multi"`.
#' @param max_signals Maximum number of multi-signal credible sets per
#'   locus/trait when `cs_mode = "multi"`.
#'
#' @return A list of class `"stblr_finemap"` with locus summaries,
#'   marker-level aggregated results, optional credible sets, locus metadata,
#'   compact run summaries, and resolved parameters.
#' @export
finemap_stblr_csr <- function(
    fit,
    Glist,
    stats,
    sets,
    trait = NULL,
    nruns = 8,
    use_residual = TRUE,
    nit = 20000,
    nburn = 2000,
    seeds = NULL,
    ve = NULL,
    vb = NULL,
    pi = NULL,
    credible_sets = TRUE,
    coverage = 0.95,
    min_r2 = 0.5,
    pip_cutoff = 0.001,
    cs_mode = c("single", "multi"),
    min_signal_pip = 0.05,
    max_signals = Inf,
    verbose = FALSE
) {
  cs_mode <- match.arg(cs_mode)
  
  if (!is.numeric(nruns) || length(nruns) != 1L || is.na(nruns) || nruns < 1) {
    stop("nruns must be a positive scalar.")
  }
  nruns <- as.integer(nruns)
  
  if (!is.numeric(nit) || length(nit) != 1L || is.na(nit) || nit < 1) {
    stop("nit must be a positive scalar.")
  }
  
  if (!is.numeric(nburn) || length(nburn) != 1L || is.na(nburn) || nburn < 0) {
    stop("nburn must be a non-negative scalar.")
  }
  
  if (!is.logical(use_residual) || length(use_residual) != 1L || is.na(use_residual)) {
    stop("use_residual must be TRUE or FALSE.")
  }
  
  if (!is.logical(credible_sets) || length(credible_sets) != 1L || is.na(credible_sets)) {
    stop("credible_sets must be TRUE or FALSE.")
  }
  
  if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
    stop("verbose must be TRUE or FALSE.")
  }
  
  .check_scalar_range(min_signal_pip, "min_signal_pip", lower = 0, upper = 1)
  
  if (!is.numeric(max_signals) || length(max_signals) != 1L ||
      is.na(max_signals) || max_signals < 0) {
    stop("max_signals must be a non-negative scalar.")
  }
  
  if (is.null(Glist$sparseLD$prefix) || !nzchar(Glist$sparseLD$prefix)) {
    stop("Glist$sparseLD$prefix is required.")
  }
  
  trait_names <- .stblr_resolve_trait_names(fit, stats)
  marker_names <- .stblr_resolve_marker_names(fit, stats, Glist)
  trait_indices <- .stblr_resolve_finemap_traits(fit, stats, trait, trait_names)
  
  sets <- .stblr_clean_finemap_sets(sets, marker_names)
  set_indices <- lapply(sets, match, table = marker_names)
  names(set_indices) <- names(sets)
  
  seeds <- .stblr_finemap_seeds(fit, nruns, seeds)
  
  map <- .stblr_finemap_marker_map(Glist, fit, marker_names)
  loci <- .stblr_finemap_loci_from_sets(sets, set_indices, map)
  
  csr <- sparseLD_read_CSR(Glist$sparseLD$prefix, one_based = TRUE)
  
  temp_dir <- tempfile("stblr_finemap_csr_")
  dir.create(temp_dir, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(temp_dir, recursive = TRUE, force = TRUE), add = TRUE)
  
  marker_rows <- list()
  locus_rows <- list()
  run_summaries <- list()
  cs_summaries <- list()
  cs_sets <- list()
  param_rows <- list()
  
  for (trait_index in trait_indices) {
    trait_name <- trait_names[trait_index]
    
    ve0 <- .stblr_get_global_parameter(fit, "ve", trait_index, override = ve)
    vb0 <- .stblr_get_global_parameter(fit, "vb", trait_index, override = vb)
    pi0 <- .stblr_get_global_parameter(fit, "pi", trait_index, override = pi)
    
    .stblr_validate_finemap_parameters(ve0, vb0, pi0)
    
    param_rows[[trait_name]] <- data.frame(
      trait = trait_name,
      ve = ve0,
      vb = vb0,
      pi = pi0,
      stringsAsFactors = FALSE
    )
    
    for (locus_name in names(sets)) {
      key <- paste(locus_name, trait_name, sep = "::")
      
      idx <- set_indices[[locus_name]]
      markers <- sets[[locus_name]]
      
      local_stats <- .stblr_make_local_stats(
        stats = stats,
        fit = fit,
        Glist = Glist,
        idx = idx,
        trait_index = trait_index,
        trait_name = trait_name,
        marker_names = marker_names,
        use_residual = use_residual,
        csr = csr
      )
      
      prefix <- file.path(
        temp_dir,
        paste0(make.names(locus_name), "_", make.names(trait_name))
      )
      
      .stblr_write_csr_subset(csr, idx, prefix)
      
      fits <- vector("list", nruns)
      
      run_local <- function(seed) {
        .stblr_run_local_csr(
          stats = local_stats,
          ld_prefix = prefix,
          ve = ve0,
          vb = vb0,
          pi = pi0,
          nit = as.integer(nit),
          nburn = as.integer(nburn),
          seed = seed
        )
      }
      
      for (run in seq_len(nruns)) {
        if (isTRUE(verbose)) {
          fits[[run]] <- run_local(seeds[run])
        } else {
          fits[[run]] <- .quiet_value(run_local(seeds[run]))
        }
      }
      
      agg <- .stblr_aggregate_finemap_runs(
        fits,
        locus = locus_name,
        trait = trait_name,
        markers = markers,
        map = map
      )
      
      marker_rows[[key]] <- agg$markers
      
      locus_rows[[key]] <- cbind(
        loci[
          loci$locus == locus_name,
          c("locus", "chr", "start", "end", "n_markers"),
          drop = FALSE
        ],
        data.frame(
          trait = trait_name,
          lead_marker = agg$summary$lead_marker,
          lead_pip = agg$summary$lead_pip,
          lead_pip_sd = agg$summary$lead_pip_sd,
          total_pip = agg$summary$total_pip,
          secondary_pip = agg$summary$total_pip - agg$summary$lead_pip,
          nruns = nruns,
          stringsAsFactors = FALSE
        )
      )[, c(
        "locus", "trait", "chr", "start", "end", "n_markers",
        "lead_marker", "lead_pip", "lead_pip_sd",
        "total_pip", "secondary_pip", "nruns"
      )]
      
      if (is.null(run_summaries[[locus_name]])) {
        run_summaries[[locus_name]] <- list()
      }
      
      run_summaries[[locus_name]][[trait_name]] <- lapply(
        fits,
        .stblr_compact_run_summary
      )
      
      if (credible_sets) {
        regional_LD <- .extract_sparseLD_region_dense(
          csr,
          idx,
          marker_names = markers
        )
        
        cs <- .stblr_finemap_credible_sets_from_ld(
          pip = stats::setNames(agg$markers$pip_mean, agg$markers$marker),
          LD = regional_LD,
          coverage = coverage,
          min_r2 = min_r2,
          pip_cutoff = pip_cutoff,
          cs_mode = cs_mode,
          min_signal_pip = min_signal_pip,
          max_signals = max_signals
        )
        
        if (is.null(cs_sets[[locus_name]])) {
          cs_sets[[locus_name]] <- list()
        }
        
        cs_sets[[locus_name]][[trait_name]] <- cs$sets
        
        if (nrow(cs$summary) > 0L) {
          loc <- loci[loci$locus == locus_name, , drop = FALSE]
          
          cs_summaries[[key]] <- cbind(
            data.frame(
              locus = locus_name,
              trait = trait_name,
              chr = loc$chr,
              start = loc$start,
              end = loc$end,
              stringsAsFactors = FALSE
            ),
            cs$summary,
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }
  
  summary <- if (length(locus_rows)) {
    do.call(rbind, locus_rows)
  } else {
    data.frame()
  }
  rownames(summary) <- NULL
  
  markers <- if (length(marker_rows)) {
    do.call(rbind, marker_rows)
  } else {
    data.frame()
  }
  rownames(markers) <- NULL
  
  cs_out <- if (credible_sets) {
    out_summary <- if (length(cs_summaries)) {
      do.call(rbind, cs_summaries)
    } else {
      data.frame()
    }
    rownames(out_summary) <- NULL
    
    list(
      summary = out_summary,
      sets = cs_sets,
      parameters = list(
        coverage = coverage,
        min_r2 = min_r2,
        pip_cutoff = pip_cutoff,
        cs_mode = cs_mode,
        min_signal_pip = min_signal_pip,
        max_signals = max_signals,
        allow_incomplete = cs_mode == "multi",
        remove = if (cs_mode == "multi") "credible_set" else "ld_neighborhood"
      )
    )
  } else {
    NULL
  }
  
  if (credible_sets && cs_mode == "multi" && nrow(summary) > 0L) {
    signal_counts <- if (!is.null(cs_out) && nrow(cs_out$summary) > 0L) {
      stats::aggregate(
        cs_out$summary$signal,
        by = list(
          locus = cs_out$summary$locus,
          trait = cs_out$summary$trait
        ),
        FUN = function(x) length(unique(x))
      )
    } else {
      data.frame(
        locus = character(),
        trait = character(),
        x = integer()
      )
    }
    
    summary$n_signals <- 0L
    
    for (i in seq_len(nrow(signal_counts))) {
      hit <- summary$locus == signal_counts$locus[i] &
        summary$trait == signal_counts$trait[i]
      
      summary$n_signals[hit] <- signal_counts$x[i]
    }
  }
  
  resolved_parameters <- if (length(param_rows)) {
    do.call(rbind, param_rows)
  } else {
    data.frame()
  }
  rownames(resolved_parameters) <- NULL
  
  out <- list(
    summary = summary,
    markers = markers,
    credible_sets = cs_out,
    loci = loci[
      ,
      setdiff(names(loci), c("markers", "indices")),
      drop = FALSE
    ],
    runs = run_summaries,
    parameters = list(
      nruns = nruns,
      nit = as.integer(nit),
      nburn = as.integer(nburn),
      use_residual = use_residual,
      verbose = verbose,
      resolved = resolved_parameters,
      credible_sets = credible_sets,
      coverage = coverage,
      min_r2 = min_r2,
      pip_cutoff = pip_cutoff,
      cs_mode = cs_mode,
      min_signal_pip = min_signal_pip,
      max_signals = max_signals
    )
  )
  
  class(out) <- "stblr_finemap"
  out
}

.quiet_value <- function(expr) {
  value <- NULL
  
  invisible(capture.output(
    value <- suppressMessages(
      suppressWarnings(
        force(expr)
      )
    ),
    type = "output"
  ))
  
  value
}

.stblr_finemap_credible_sets_from_ld <- function(
    pip,
    LD,
    coverage,
    min_r2,
    pip_cutoff,
    cs_mode = c("single", "multi"),
    min_signal_pip = 0.05,
    max_signals = Inf
) {
  cs_mode <- match.arg(cs_mode)
  if (cs_mode == "single") {
    return(make_credible_sets_from_ld(
      pip = pip,
      LD = LD,
      coverage = coverage,
      min_r2 = min_r2,
      pip_cutoff = pip_cutoff,
      allow_incomplete = FALSE,
      remove = "ld_neighborhood"
    ))
  }

  make_multisignal_credible_sets_from_ld(
    pip = pip,
    LD = LD,
    coverage = coverage,
    min_r2 = min_r2,
    pip_cutoff = pip_cutoff,
    min_signal_pip = min_signal_pip,
    max_signals = max_signals,
    remove = "credible_set",
    allow_incomplete = TRUE
  )
}

.stblr_extract_trait_pip <- function(fit, trait) {
  .stblr_extract_pip(fit, trait = trait)
}

.stblr_resolve_trait_names <- function(fit, stats) {
  if (!is.null(stats$trait_names)) {
    return(as.character(stats$trait_names))
  }
  if (!is.null(names(stats$yy))) {
    return(as.character(names(stats$yy)))
  }
  if (!is.null(fit$dm) && (is.matrix(fit$dm) || is.data.frame(fit$dm)) &&
      !is.null(colnames(fit$dm))) {
    return(as.character(colnames(fit$dm)))
  }
  nt <- if (!is.null(stats$yy)) length(stats$yy) else if (!is.null(fit$dm) &&
      (is.matrix(fit$dm) || is.data.frame(fit$dm))) ncol(fit$dm) else 1L
  paste0("T", seq_len(nt))
}

.stblr_resolve_marker_names <- function(fit, stats, Glist) {
  marker_names <- stats$marker_names
  if (is.null(marker_names) && !is.null(stats$ww)) {
    if (is.list(stats$ww) && length(stats$ww) > 0L) {
      marker_names <- names(stats$ww[[1L]])
    } else if (!is.null(rownames(stats$ww))) {
      marker_names <- rownames(stats$ww)
    } else {
      marker_names <- names(stats$ww)
    }
  }
  if (is.null(marker_names) && !is.null(stats$wy)) {
    if (is.list(stats$wy) && length(stats$wy) > 0L) {
      marker_names <- names(stats$wy[[1L]])
    } else if (!is.null(rownames(stats$wy))) {
      marker_names <- rownames(stats$wy)
    } else {
      marker_names <- names(stats$wy)
    }
  }
  if (is.null(marker_names) && !is.null(fit$dm)) {
    marker_names <- if (is.matrix(fit$dm) || is.data.frame(fit$dm)) {
      rownames(fit$dm)
    } else {
      names(fit$dm)
    }
  }
  if (is.null(marker_names)) {
    map <- tryCatch(.stblr_marker_map_from_Glist(Glist, fit = fit), error = function(e) NULL)
    if (!is.null(map)) marker_names <- map$marker
  }
  if (is.null(marker_names) || anyNA(marker_names) || any(!nzchar(marker_names))) {
    stop("Marker names are required in stats, fit$dm, or Glist.")
  }
  marker_names <- as.character(marker_names)
  if (anyDuplicated(marker_names)) {
    stop("Marker names must be unique.")
  }
  marker_names
}

.stblr_resolve_finemap_traits <- function(fit, stats, trait, trait_names) {
  nt <- length(trait_names)
  if (is.null(trait)) {
    return(seq_len(nt))
  }
  if (is.character(trait)) {
    idx <- match(trait, trait_names)
    if (anyNA(idx)) stop("trait was not found in resolved trait names.")
    return(idx)
  }
  idx <- as.integer(trait)
  if (length(idx) < 1L || anyNA(idx) || any(idx < 1L | idx > nt)) {
    stop("trait must be NULL, valid trait names, or valid trait indices.")
  }
  idx
}

.stblr_clean_finemap_sets <- function(sets, marker_names) {
  if (!is.list(sets)) stop("sets must be a list.")
  if (length(sets) < 1L) stop("sets must contain at least one marker set.")
  nms <- names(sets)
  if (is.null(nms)) nms <- rep("", length(sets))
  missing_names <- !nzchar(nms)
  nms[missing_names] <- paste0("locus", which(missing_names))
  names(sets) <- nms

  marker_rank <- seq_along(marker_names)
  names(marker_rank) <- marker_names
  out <- list()
  dropped <- character()
  for (nm in names(sets)) {
    markers <- unique(as.character(sets[[nm]]))
    markers <- markers[markers %in% marker_names]
    if (length(markers) < 1L) {
      dropped <- c(dropped, nm)
      next
    }
    out[[nm]] <- markers[order(marker_rank[markers])]
  }
  if (length(dropped) > 0L) {
    warning("Dropping empty fine-mapping sets: ", paste(dropped, collapse = ", "))
  }
  if (length(out) < 1L) stop("No non-empty fine-mapping sets remain.")
  out
}

.stblr_get_global_parameter <- function(fit, parameter, trait, override = NULL) {
  parameter <- match.arg(parameter, c("ve", "vb", "pi"))
  if (!is.null(override)) {
    return(.stblr_select_parameter_value(override, trait, fit, parameter))
  }
  trace_name <- switch(parameter, ve = "ves", vb = "vbs", pi = "pis")
  if (!is.null(fit[[trace_name]])) {
    x <- fit[[trace_name]]
    burn <- if (!is.null(fit$input$nburn)) as.integer(fit$input$nburn) else 0L
    val <- .stblr_trace_mean_after_burnin(x, trait, burn)
    if (is.finite(val)) return(val)
  }
  if (!is.null(fit[[parameter]])) {
    val <- .stblr_select_parameter_value(fit[[parameter]], trait, fit, parameter)
    if (is.finite(val)) return(val)
  }
  input_name <- switch(parameter, ve = "E", vb = "B", pi = "pi_prior_mean")
  if (!is.null(fit$input[[input_name]])) {
    val <- .stblr_select_parameter_value(fit$input[[input_name]], trait, fit, parameter)
    if (parameter == "pi" && length(val) == 2L) val <- val[2L]
    if (is.finite(val)) return(val)
  }
  if (parameter == "pi" && !is.null(fit$input$pi) && length(fit$input$pi) >= 2L) {
    return(as.numeric(fit$input$pi[2L]))
  }
  stop("Could not resolve global parameter '", parameter, "' from fit; supply an override.")
}

.stblr_trace_mean_after_burnin <- function(x, trait, burn) {
  if (is.matrix(x) || is.data.frame(x)) {
    x <- as.matrix(x)
    col <- .stblr_parameter_trait_column(x, trait)
    keep <- seq.int(min(nrow(x), burn) + 1L, nrow(x))
    if (length(keep) < 1L) keep <- seq_len(nrow(x))
    return(mean(as.numeric(x[keep, col]), na.rm = TRUE))
  }
  if (is.list(x) && !is.data.frame(x)) {
    xx <- x[[as.integer(trait)]]
    keep <- seq.int(min(length(xx), burn) + 1L, length(xx))
    if (length(keep) < 1L) keep <- seq_along(xx)
    return(mean(as.numeric(xx[keep]), na.rm = TRUE))
  }
  xx <- as.numeric(x)
  keep <- seq.int(min(length(xx), burn) + 1L, length(xx))
  if (length(keep) < 1L) keep <- seq_along(xx)
  mean(xx[keep], na.rm = TRUE)
}

.stblr_select_parameter_value <- function(x, trait, fit, parameter) {
  if (is.matrix(x) || is.data.frame(x)) {
    x <- as.matrix(x)
    if (nrow(x) == ncol(x) && parameter %in% c("ve", "vb")) {
      i <- as.integer(trait)
      return(as.numeric(x[i, i]))
    }
    return(as.numeric(x[.stblr_parameter_trait_column(x, trait)]))
  }
  if (!is.null(names(x)) && is.character(trait)) {
    return(as.numeric(x[[trait]]))
  }
  if (!is.null(names(x)) && !is.null(fit$dm) && is.matrix(fit$dm) &&
      !is.null(colnames(fit$dm))) {
    nm <- colnames(fit$dm)[as.integer(trait)]
    if (!is.na(nm) && nm %in% names(x)) return(as.numeric(x[[nm]]))
  }
  x <- as.numeric(x)
  if (length(x) == 1L) return(x)
  if (parameter == "pi" && length(x) == 2L && all(x >= 0) && abs(sum(x) - 1) < 1e-8) {
    return(x[2L])
  }
  x[as.integer(trait)]
}

.stblr_parameter_trait_column <- function(x, trait) {
  if (is.character(trait)) {
    if (is.null(colnames(x)) || !(trait %in% colnames(x))) {
      stop("Requested trait was not found in parameter trace columns.")
    }
    return(match(trait, colnames(x)))
  }
  col <- as.integer(trait)
  if (length(col) != 1L || is.na(col) || col < 1L || col > ncol(x)) {
    stop("Invalid trait index for parameter trace.")
  }
  col
}

.stblr_make_local_stats <- function(stats, fit, Glist, idx, trait_index,
                                    trait_name, marker_names, use_residual,
                                    csr = NULL) {
  wy <- if (use_residual) {
    .stblr_residualize_wy_csr(stats, fit, Glist, idx, trait_index, csr = csr)
  } else {
    .stblr_extract_stat_vector(stats$wy, trait_index)[idx]
  }
  ww <- .stblr_extract_stat_vector(stats$ww, trait_index)[idx]
  yy <- .stblr_extract_stat_scalar(stats$yy, trait_index)
  markers <- marker_names[idx]
  names(wy) <- markers
  names(ww) <- markers
  yy <- stats::setNames(as.numeric(yy), trait_name)
  list(
    wy = stats::setNames(list(as.numeric(wy)), trait_name),
    ww = stats::setNames(list(as.numeric(ww)), trait_name),
    yy = yy,
    n = stats$n,
    m = length(idx),
    marker_names = markers,
    trait_names = trait_name
  )
}

.stblr_extract_stat_vector <- function(x, trait_index) {
  if (is.list(x) && !is.data.frame(x)) {
    return(as.numeric(x[[trait_index]]))
  }
  if (is.matrix(x) || is.data.frame(x)) {
    x <- as.matrix(x)
    return(as.numeric(x[, trait_index]))
  }
  as.numeric(x)
}

.stblr_extract_stat_scalar <- function(x, trait_index) {
  if (is.matrix(x) || is.data.frame(x)) {
    return(as.numeric(as.matrix(x)[trait_index, trait_index]))
  }
  as.numeric(x)[trait_index]
}

.stblr_residualize_wy_csr <- function(stats, fit, Glist, idx, trait, csr = NULL) {
  wy <- .stblr_extract_stat_vector(stats$wy, trait)
  bm <- .stblr_extract_bm(fit, trait)
  if (length(wy) != length(bm)) {
    stop("stats$wy and fit$bm have different marker counts; cannot residualize.")
  }
  if (is.null(csr)) csr <- sparseLD_read_CSR(Glist$sparseLD$prefix, one_based = TRUE)
  ld_b <- .stblr_csr_matvec(csr, bm)
  LD_LL <- .extract_sparseLD_region_dense(csr, idx)
  as.numeric(wy[idx] - ld_b[idx] + as.numeric(LD_LL %*% bm[idx]))
}

.stblr_extract_bm <- function(fit, trait) {
  if (is.null(fit$bm)) stop("fit$bm is required for residualized fine mapping.")
  bm <- fit$bm
  if (is.matrix(bm) || is.data.frame(bm)) {
    return(as.numeric(as.matrix(bm)[, trait]))
  }
  as.numeric(bm)
}

.stblr_csr_matvec <- function(csr, x) {
  ptr <- as.numeric(.csr_component(csr, c("ptr", "row_ptr", "indptr")))
  col_idx <- as.integer(.csr_component(csr, c("idx", "col_idx", "indices")))
  values <- as.numeric(.csr_component(csr, c("x", "values", "val")))
  m <- length(ptr) - 1L
  if (length(x) != m) stop("Length of x does not match CSR dimension.")
  out <- as.numeric(x)
  for (i in seq_len(m)) {
    start <- ptr[i] + 1
    end <- ptr[i + 1L]
    if (!is.finite(start) || !is.finite(end) || end < start) next
    pos <- seq.int(start, end)
    j <- col_idx[pos]
    v <- values[pos]
    out[i] <- out[i] + sum(v * x[j])
    out[j] <- out[j] + v * x[i]
  }
  out
}

.stblr_write_csr_subset <- function(csr, idx, prefix) {
  idx <- as.integer(idx)
  ptr <- as.numeric(.csr_component(csr, c("ptr", "row_ptr", "indptr")))
  col_idx <- as.integer(.csr_component(csr, c("idx", "col_idx", "indices")))
  values <- as.numeric(.csr_component(csr, c("x", "values", "val")))
  local <- seq_along(idx)
  names(local) <- as.character(idx)
  rows_cols <- vector("list", length(idx))
  rows_vals <- vector("list", length(idx))
  nnz <- 0
  for (a in seq_along(idx)) {
    row <- idx[a]
    start <- ptr[row] + 1
    end <- ptr[row + 1L]
    if (!is.finite(start) || !is.finite(end) || end < start) {
      rows_cols[[a]] <- integer()
      rows_vals[[a]] <- numeric()
      next
    }
    pos <- seq.int(start, end)
    keep <- as.character(col_idx[pos]) %in% names(local)
    jj <- unname(local[as.character(col_idx[pos][keep])])
    keep_upper <- jj > a
    jj <- jj[keep_upper]
    vv <- values[pos][keep][keep_upper]
    rows_cols[[a]] <- as.integer(jj - 1L)
    rows_vals[[a]] <- as.numeric(vv)
    nnz <- nnz + length(jj)
  }
  row_ptr <- c(0, cumsum(vapply(rows_cols, length, integer(1))))
  col0 <- unlist(rows_cols, use.names = FALSE)
  vals <- unlist(rows_vals, use.names = FALSE)
  .stblr_write_uint64_file(paste0(prefix, ".row_ptr.u64.bin"), row_ptr)
  .stblr_write_uint32_file(paste0(prefix, ".col_idx.u32.0based.bin"), col0)
  writeBin(as.numeric(vals), paste0(prefix, ".values.f32.bin"), size = 4, endian = "little")
  meta <- c(
    "format=sparse_ld_csr",
    "storage=streamed_upper_triangle",
    paste0("n_bed=", NA_integer_),
    paste0("n_used=", NA_integer_),
    paste0("n_samples_used=", NA_integer_),
    paste0("n_variants=", length(idx)),
    paste0("nnz=", nnz),
    "triangle=upper",
    "diagonal=implicit_1",
    paste0("row_ptr_file=", prefix, ".row_ptr.u64.bin"),
    paste0("col_idx_file=", prefix, ".col_idx.u32.0based.bin"),
    paste0("values_file=", prefix, ".values.f32.bin"),
    "row_ptr_type=uint64",
    "col_idx_type=uint32",
    "values_type=float32",
    "index_base=0",
    "value=r"
  )
  writeLines(meta, paste0(prefix, ".meta.txt"))
  invisible(prefix)
}

.stblr_write_uint64_file <- function(path, x) {
  con <- file(path, "wb")
  on.exit(close(con), add = TRUE)
  for (value in as.numeric(x)) {
    writeBin(.stblr_uint_le_raw(value, 8L), con)
  }
}

.stblr_write_uint32_file <- function(path, x) {
  con <- file(path, "wb")
  on.exit(close(con), add = TRUE)
  for (value in as.numeric(x)) {
    writeBin(.stblr_uint_le_raw(value, 4L), con)
  }
}

.stblr_uint_le_raw <- function(value, nbytes) {
  if (!is.finite(value) || value < 0 || value != floor(value)) {
    stop("Unsigned integer value must be a finite non-negative integer.")
  }
  out <- raw(nbytes)
  for (i in seq_len(nbytes)) {
    byte <- value %% 256
    out[i] <- as.raw(byte)
    value <- floor(value / 256)
  }
  if (value != 0) stop("Unsigned integer value exceeds target byte width.")
  out
}

.stblr_run_local_csr <- function(stats, ld_prefix, ve, vb, pi, nit, nburn, seed) {
  m <- stats$m
  trait_name <- stats$trait_names[1L]
  B <- matrix(vb, nrow = 1L, ncol = 1L, dimnames = list(trait_name, trait_name))
  E <- matrix(ve, nrow = 1L, ncol = 1L, dimnames = list(trait_name, trait_name))
  ssb_prior <- list(stats::setNames(vb, trait_name))
  sse_prior <- list(stats::setNames(ve, trait_name))
  raw <- stblr_cpg_omp_csr(
    wy = stats$wy,
    ww = stats$ww,
    yy = stats$yy,
    b_init = list(rep(0, m)),
    d_init = list(rep(0, m)),
    use_d_init = FALSE,
    r_init = stats$wy,
    use_r_init = FALSE,
    rebuild_r_before_updateE = FALSE,
    ld_prefix = ld_prefix,
    B = B,
    E = E,
    ssb_prior = ssb_prior,
    sse_prior = sse_prior,
    pi = c(1 - pi, pi),
    nub = 4,
    nue = 4,
    updateB = FALSE,
    updateE = FALSE,
    updatePi = FALSE,
    adjE = 0.9,
    n = as.integer(stats$n),
    nit = as.integer(nit),
    nburn = as.integer(nburn),
    nthin = 1L,
    pi_prior_a = 1,
    pi_prior_b = 1,
    ncores = 1L,
    seed = as.integer(seed)
  )
  fit <- .format_stblr_fit(raw, 1L, m, trait_name, stats$marker_names)
  fit$input <- list(
    n = stats$n, m = m, nt = 1L, B = B, E = E,
    updateB = FALSE, updateE = FALSE, updatePi = FALSE,
    nit = nit, nburn = nburn, seed = seed, ld_prefix = ld_prefix,
    pi = c(1 - pi, pi), pi_prior_mean = pi
  )
  fit
}

.stblr_aggregate_finemap_runs <- function(fits, locus, trait, markers, map = NULL) {
  dm <- do.call(cbind, lapply(fits, function(x) as.numeric(x$dm[, 1L])))
  bm <- do.call(cbind, lapply(fits, function(x) as.numeric(x$bm[, 1L])))
  rownames(dm) <- rownames(bm) <- markers
  marker_map <- if (!is.null(map)) map[match(markers, map$marker), , drop = FALSE] else NULL
  chr <- pos <- rep(NA_real_, length(markers))
  if (!is.null(marker_map) && nrow(marker_map) == length(markers)) {
    chr <- marker_map$chr
    pos <- marker_map$pos
  }
  marker_tab <- data.frame(
    locus = locus,
    trait = trait,
    marker = markers,
    chr = chr,
    pos = pos,
    pip_mean = rowMeans(dm, na.rm = TRUE),
    pip_sd = apply(dm, 1L, stats::sd, na.rm = TRUE),
    pip_min = apply(dm, 1L, min, na.rm = TRUE),
    pip_max = apply(dm, 1L, max, na.rm = TRUE),
    bm_mean = rowMeans(bm, na.rm = TRUE),
    bm_sd = apply(bm, 1L, stats::sd, na.rm = TRUE),
    bm_min = apply(bm, 1L, min, na.rm = TRUE),
    bm_max = apply(bm, 1L, max, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
  lead <- which.max(marker_tab$pip_mean)
  summary <- data.frame(
    locus = locus,
    trait = trait,
    lead_marker = marker_tab$marker[lead],
    lead_pip = marker_tab$pip_mean[lead],
    lead_pip_sd = marker_tab$pip_sd[lead],
    total_pip = sum(marker_tab$pip_mean),
    stringsAsFactors = FALSE
  )
  list(markers = marker_tab, summary = summary)
}

.stblr_finemap_marker_map <- function(Glist, fit, marker_names) {
  map <- tryCatch(.stblr_marker_map_from_Glist(Glist, fit = fit), error = function(e) NULL)
  if (is.null(map)) {
    map <- data.frame(
      marker = marker_names,
      chr = NA_integer_,
      pos = NA_real_,
      index = seq_along(marker_names),
      stringsAsFactors = FALSE
    )
  }
  map <- map[match(marker_names, map$marker), , drop = FALSE]
  missing <- is.na(map$marker)
  if (any(missing)) {
    map[missing, "marker"] <- marker_names[missing]
    map[missing, "index"] <- which(missing)
  }
  rownames(map) <- NULL
  map
}

.stblr_finemap_loci_from_sets <- function(sets, set_indices, map) {
  rows <- lapply(names(sets), function(nm) {
    markers <- sets[[nm]]
    loc_map <- map[set_indices[[nm]], , drop = FALSE]
    chr <- if (length(unique(stats::na.omit(loc_map$chr))) == 1L) loc_map$chr[which(!is.na(loc_map$chr))[1L]] else NA
    start <- if (all(is.na(loc_map$pos))) NA_real_ else min(loc_map$pos, na.rm = TRUE)
    end <- if (all(is.na(loc_map$pos))) NA_real_ else max(loc_map$pos, na.rm = TRUE)
    row <- data.frame(
      locus = nm,
      chr = chr,
      start = start,
      end = end,
      n_markers = length(markers),
      stringsAsFactors = FALSE
    )
    .add_list_columns(row, markers = markers, indices = set_indices[[nm]])
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

.stblr_finemap_seeds <- function(fit, nruns, seeds) {
  if (is.null(seeds)) {
    base <- if (!is.null(fit$input$seed) && is.finite(fit$input$seed)) {
      as.integer(fit$input$seed)
    } else {
      0L
    }
    return(base + seq_len(nruns))
  }
  if (length(seeds) != nruns) stop("seeds must have length nruns.")
  as.integer(seeds)
}

.stblr_compact_run_summary <- function(fit) {
  list(dm = fit$dm, bm = fit$bm, vbs = fit$vbs, ves = fit$ves, pis = fit$pis)
}

.stblr_validate_finemap_parameters <- function(ve, vb, pi) {
  if (!is.numeric(ve) || length(ve) != 1L || !is.finite(ve) || ve <= 0) {
    stop("Resolved ve must be a positive finite scalar.")
  }
  if (!is.numeric(vb) || length(vb) != 1L || !is.finite(vb) || vb <= 0) {
    stop("Resolved vb must be a positive finite scalar.")
  }
  if (!is.numeric(pi) || length(pi) != 1L || !is.finite(pi) || pi <= 0 || pi >= 1) {
    stop("Resolved pi must be a finite scalar in (0, 1).")
  }
  invisible(TRUE)
}
