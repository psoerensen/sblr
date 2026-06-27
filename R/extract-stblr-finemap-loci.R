#' Extract Fine-Mapping Loci from a Genome-Wide CSR ST-BLR Fit
#'
#' Summarizes predefined fine-mapping loci from an existing genome-wide
#' [stblr_csr()] posterior fit. This is a post-processing helper: it does not
#' run local MCMC, does not write temporary CSR subsets, and does not refit the
#' model. Marker-level PIPs and posterior mean effects are read directly from
#' `fit$dm` and `fit$bm`.
#'
#' A typical workflow is to fit the genome-wide CSR model first, optionally with
#' LD-swap proposals, define loci from the fitted PIPs with
#' [make_stblr_credible_sets()], and then extract locus summaries:
#'
#' ```r
#' fitMH <- stblr_csr(
#'   stats = stats,
#'   Glist = Glist,
#'   pi_init = 0.001,
#'   pi_prior_a = 1,
#'   pi_prior_b = 1,
#'   h2 = 0.3,
#'   adjE = 0.9,
#'   nit = 1000,
#'   nburn = 100,
#'   ncores = 3,
#'   seed = 10,
#'   scheduled = FALSE,
#'   updateLDswap = TRUE,
#'   ld_swap_prob = 0.10,
#'   ld_swap_r2 = 0.05,
#'   ld_swap_moves = 5
#' )
#'
#' cs_global <- make_stblr_credible_sets(
#'   fit = fitMH,
#'   Glist = Glist,
#'   trait = "D1",
#'   coverage = 0.95,
#'   min_r2 = 0.5,
#'   pip_cutoff = 0.001,
#'   locus_pip_cutoff = 0.01,
#'   max_locus_distance = 1e6
#' )
#'
#' fm <- extract_stblr_finemap_loci(
#'   fit = fitMH,
#'   Glist = Glist,
#'   locus_sets = cs_global$locus_sets,
#'   trait = "D1",
#'   credible_sets = TRUE,
#'   cs_mode = "multi"
#' )
#' ```
#'
#' @param fit Genome-wide `stblr_csr()` fit containing marker PIPs in `fit$dm`
#'   and posterior mean effects in `fit$bm`.
#' @param Glist Genotype/sparse LD object used to resolve marker chromosome and
#'   position information. If `credible_sets = TRUE`, `Glist$sparseLD$prefix`
#'   must identify a readable sparse LD CSR prefix.
#' @param locus_sets Named or unnamed list of marker vectors defining loci to
#'   summarize, for example `cs_global$locus_sets` from
#'   [make_stblr_credible_sets()].
#' @param trait Trait column name or index in `fit$dm` and `fit$bm`.
#' @param credible_sets Logical; construct credible sets from the genome-wide
#'   PIPs restricted to each locus.
#' @param coverage,min_r2,pip_cutoff Credible-set construction parameters.
#' @param cs_mode Credible-set mode. `"single"` calls
#'   [make_credible_sets_from_ld()]. `"multi"` calls
#'   [make_multisignal_credible_sets_from_ld()].
#' @param min_signal_pip Minimum remaining lead-marker PIP required to start a
#'   new signal when `cs_mode = "multi"`.
#' @param max_signals Maximum number of multi-signal credible sets per locus
#'   when `cs_mode = "multi"`.
#'
#' @return A list of class `"stblr_finemap"` with `summary`, `markers`,
#'   optional `credible_sets`, `loci`, and `parameters`. Unlike
#'   [finemap_stblr_csr()], no `runs` element is returned because no local
#'   chains are run.
#' @export
extract_stblr_finemap_loci <- function(
    fit,
    Glist,
    locus_sets,
    trait = 1,
    credible_sets = TRUE,
    coverage = 0.95,
    min_r2 = 0.5,
    pip_cutoff = 0.001,
    cs_mode = c("single", "multi"),
    min_signal_pip = 0.05,
    max_signals = Inf
) {
  cs_mode <- match.arg(cs_mode)
  if (!is.logical(credible_sets) || length(credible_sets) != 1L ||
      is.na(credible_sets)) {
    stop("credible_sets must be TRUE or FALSE.")
  }
  .check_scalar_range(coverage, "coverage", lower = 0, upper = 1,
                      lower_open = TRUE, upper_open = TRUE)
  .check_scalar_range(min_r2, "min_r2", lower = 0, upper = 1)
  .check_scalar_range(pip_cutoff, "pip_cutoff", lower = 0, upper = 1)
  .check_scalar_range(min_signal_pip, "min_signal_pip", lower = 0, upper = 1)
  if (!is.numeric(max_signals) || length(max_signals) != 1L ||
      is.na(max_signals) || max_signals < 0) {
    stop("max_signals must be a non-negative scalar.")
  }

  pip_info <- .stblr_extract_pip(fit, trait = trait)
  pip <- pip_info$pip
  trait_name <- as.character(pip_info$trait)
  if (is.null(names(pip)) || anyNA(names(pip)) || any(!nzchar(names(pip)))) {
    stop("fit$dm must have marker names for fine-mapping extraction.")
  }

  locus_sets <- .stblr_clean_marker_sets(locus_sets, pip)
  marker_names <- names(pip)

  bm <- .stblr_extract_finemap_fit_stat(
    fit, "bm", trait = trait, marker_names = marker_names, required = TRUE
  )
  pip_sd <- .stblr_extract_finemap_fit_stat(
    fit, "dm_sd", trait = trait, marker_names = marker_names, required = FALSE
  )
  pip_min <- .stblr_extract_finemap_fit_stat(
    fit, "dm_min", trait = trait, marker_names = marker_names, required = FALSE
  )
  pip_max <- .stblr_extract_finemap_fit_stat(
    fit, "dm_max", trait = trait, marker_names = marker_names, required = FALSE
  )
  bm_sd <- .stblr_extract_finemap_fit_stat(
    fit, "bm_sd", trait = trait, marker_names = marker_names, required = FALSE
  )
  bm_min <- .stblr_extract_finemap_fit_stat(
    fit, "bm_min", trait = trait, marker_names = marker_names, required = FALSE
  )
  bm_max <- .stblr_extract_finemap_fit_stat(
    fit, "bm_max", trait = trait, marker_names = marker_names, required = FALSE
  )

  map <- .stblr_marker_map_from_Glist(Glist, fit = fit)
  map <- map[match(marker_names, map$marker), , drop = FALSE]

  csr <- NULL
  if (credible_sets) {
    if (is.null(Glist$sparseLD$prefix) || !nzchar(Glist$sparseLD$prefix)) {
      stop("Glist$sparseLD$prefix is required when credible_sets = TRUE.")
    }
    csr <- sparseLD_read_CSR(Glist$sparseLD$prefix, one_based = TRUE)
  }

  marker_rows <- list()
  locus_rows <- list()
  cs_summaries <- list()
  cs_sets <- list()

  for (locus_name in names(locus_sets)) {
    markers <- locus_sets[[locus_name]]
    idx <- match(markers, marker_names)
    loc_map <- map[idx, , drop = FALSE]

    chr <- if (nrow(loc_map) > 0L &&
        length(unique(stats::na.omit(loc_map$chr))) == 1L) {
      loc_map$chr[which(!is.na(loc_map$chr))[1L]]
    } else {
      NA
    }
    start <- if (nrow(loc_map) > 0L && any(!is.na(loc_map$pos))) {
      min(loc_map$pos, na.rm = TRUE)
    } else {
      NA_real_
    }
    end <- if (nrow(loc_map) > 0L && any(!is.na(loc_map$pos))) {
      max(loc_map$pos, na.rm = TRUE)
    } else {
      NA_real_
    }

    marker_tab <- data.frame(
      locus = locus_name,
      trait = trait_name,
      marker = markers,
      chr = loc_map$chr,
      pos = loc_map$pos,
      pip_mean = as.numeric(pip[markers]),
      pip_sd = as.numeric(pip_sd[markers]),
      pip_min = as.numeric(pip_min[markers]),
      pip_max = as.numeric(pip_max[markers]),
      bm_mean = as.numeric(bm[markers]),
      bm_sd = as.numeric(bm_sd[markers]),
      bm_min = as.numeric(bm_min[markers]),
      bm_max = as.numeric(bm_max[markers]),
      stringsAsFactors = FALSE
    )
    marker_rows[[locus_name]] <- marker_tab

    lead <- which.max(marker_tab$pip_mean)
    lead_pip <- marker_tab$pip_mean[lead]
    lead_pip_sd <- marker_tab$pip_sd[lead]
    total_pip <- sum(marker_tab$pip_mean, na.rm = TRUE)

    locus_rows[[locus_name]] <- data.frame(
      locus = locus_name,
      trait = trait_name,
      chr = chr,
      start = start,
      end = end,
      n_markers = length(markers),
      lead_marker = marker_tab$marker[lead],
      lead_pip = lead_pip,
      lead_pip_sd = lead_pip_sd,
      total_pip = total_pip,
      secondary_pip = total_pip - lead_pip,
      stringsAsFactors = FALSE
    )

    if (credible_sets) {
      sparse_idx <- loc_map$index
      if (anyNA(sparse_idx)) {
        stop("Could not resolve sparse LD indices for locus ", locus_name, ".")
      }
      regional_LD <- .extract_sparseLD_region_dense(
        csr,
        sparse_idx,
        marker_names = markers
      )
      cs <- if (cs_mode == "single") {
        make_credible_sets_from_ld(
          pip = stats::setNames(marker_tab$pip_mean, marker_tab$marker),
          LD = regional_LD,
          coverage = coverage,
          min_r2 = min_r2,
          pip_cutoff = pip_cutoff,
          allow_incomplete = FALSE,
          remove = "ld_neighborhood"
        )
      } else {
        make_multisignal_credible_sets_from_ld(
          pip = stats::setNames(marker_tab$pip_mean, marker_tab$marker),
          LD = regional_LD,
          coverage = coverage,
          min_r2 = min_r2,
          pip_cutoff = pip_cutoff,
          min_signal_pip = min_signal_pip,
          max_signals = max_signals,
          remove = "credible_set",
          allow_incomplete = TRUE
        )
      }

      cs_sets[[locus_name]] <- list()
      cs_sets[[locus_name]][[trait_name]] <- cs$sets

      if (nrow(cs$summary) > 0L) {
        cs_summaries[[locus_name]] <- cbind(
          data.frame(
            locus = locus_name,
            trait = trait_name,
            chr = chr,
            start = start,
            end = end,
            stringsAsFactors = FALSE
          ),
          cs$summary,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  summary <- if (length(locus_rows)) do.call(rbind, locus_rows) else data.frame()
  rownames(summary) <- NULL

  markers <- if (length(marker_rows)) do.call(rbind, marker_rows) else data.frame()
  rownames(markers) <- NULL

  loci <- summary[, c("locus", "trait", "chr", "start", "end", "n_markers"),
                  drop = FALSE]

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

  out <- list(
    summary = summary,
    markers = markers,
    credible_sets = cs_out,
    loci = loci,
    parameters = list(
      trait = trait_name,
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

.stblr_extract_finemap_fit_stat <- function(fit, component, trait,
                                            marker_names, required) {
  x <- fit[[component]]
  if (is.null(x)) {
    if (required) stop("fit must contain marker values in fit$", component, ".")
    return(stats::setNames(rep(NA_real_, length(marker_names)), marker_names))
  }

  if (is.matrix(x) || is.data.frame(x)) {
    x <- as.matrix(x)
    if (is.character(trait)) {
      if (is.null(colnames(x)) || !(trait %in% colnames(x))) {
        stop("trait was not found in colnames(fit$", component, ").")
      }
      values <- x[, match(trait, colnames(x))]
    } else {
      trait_index <- as.integer(trait)
      if (length(trait_index) != 1L || is.na(trait_index) ||
          trait_index < 1L || trait_index > ncol(x)) {
        stop("trait must be a valid fit$", component, " column index.")
      }
      values <- x[, trait_index]
    }
    value_names <- rownames(x)
  } else {
    values <- as.numeric(x)
    value_names <- names(x)
  }

  values <- as.numeric(values)
  if (is.null(value_names)) {
    if (length(values) != length(marker_names)) {
      stop("fit$", component, " must have marker names or match fit$dm length.")
    }
    value_names <- marker_names
  }
  names(values) <- value_names

  missing <- setdiff(marker_names, names(values))
  if (length(missing) > 0L) {
    stop("fit$", component, " is missing markers: ",
         paste(head(missing, 5L), collapse = ", "),
         if (length(missing) > 5L) ", ..." else "")
  }

  values[marker_names]
}
