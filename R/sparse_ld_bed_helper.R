#' Compute BED Marker Sufficient Statistics
#'
#' Computes marker cross-products and marker-trait cross-products directly from
#' a PLINK BED file.
#'
#' @param bed_file Path to a PLINK BED file.
#' @param n Number of individuals in the BED file.
#' @param cls Marker column indices.
#' @param af Allele frequencies corresponding to `cls`.
#' @param y Phenotype vector or matrix. If `y` is a vector, `names(y)` must
#'   contain individual IDs matching `Glist$ids`. If `y` is a matrix,
#'   `rownames(y)` must contain matching IDs unless `nrow(y) == Glist$n`.
#' @param rows Advanced/internal option. Optional 1-based BED/FAM row indices
#'   corresponding to rows of `y`. Normally omitted; rows are inferred from
#'   `names(y)` or `rownames(y)` matched to `Glist$ids`.
#' @param scale Standardize markers using allele frequencies.
#' @param nthreads Number of OpenMP threads.
#' @param MG,JB,TB Internal blocking and tiling controls.
#' @return A list containing sufficient statistics for ST-BLR fitting.
#' @name bed_xtx_xty
#' @export
NULL

#' Stream Sparse LD from PLINK BED Files to CSR Files
#'
#' Computes sparse linkage disequilibrium from PLINK BED files and writes a
#' disk-backed CSR representation.
#'
#' @param bed_files Paths to PLINK BED files.
#' @param n Number of individuals in the BED files.
#' @param cls List of marker column indices, one vector per BED file.
#' @param out_prefix Output prefix for the CSR files.
#' @param rows Optional 1-based individual row indices.
#' @param af Optional list of allele-frequency vectors.
#' @param pos_bp Optional marker base-pair positions.
#' @param max_distance_bp Maximum base-pair distance between retained pairs.
#'   A value of zero disables base-pair distance filtering.
#' @param max_distance_variants Maximum marker-index distance between retained
#'   pairs. A value of zero disables marker-index distance filtering.
#' @param r2_threshold Minimum squared-correlation threshold.
#' @param block_size Marker block size.
#' @param nthreads Number of OpenMP threads.
#' @param allow_full_ld Permit both distance filters to be disabled. When
#'   \code{FALSE}, the default, \code{sparseLD_stream_CSR()} stops if
#'   \code{max_distance_variants <= 0} and no positive base-pair filter can be
#'   applied. Setting both distance limits to zero evaluates all marker pairs
#'   before \code{r2_threshold} filtering and can be very expensive. A local LD
#'   default such as \code{max_distance_variants = 1000} is usually safer.
#' @return The output prefix and sparse-LD writing summary.
#' @name sparseLD_stream_CSR
#' @export
sparseLD_stream_CSR <- function(bed_files, n, cls, out_prefix, rows = NULL,
                                af = NULL, pos_bp = NULL,
                                max_distance_bp = 1000000L,
                                max_distance_variants = 1000L,
                                r2_threshold = 0.01,
                                block_size = 1024L,
                                nthreads = 1L,
                                allow_full_ld = FALSE) {
  pos_bp_empty <- is.null(pos_bp) ||
    length(pos_bp) == 0L ||
    (is.list(pos_bp) && sum(lengths(pos_bp)) == 0L)
  bp_filter_disabled <- pos_bp_empty || max_distance_bp <= 0
  variant_filter_disabled <- max_distance_variants <= 0

  if (!isTRUE(allow_full_ld) &&
      variant_filter_disabled &&
      bp_filter_disabled) {
    stop(
      "Both sparse LD distance filters are disabled. ",
      "This means all marker pairs will be evaluated before r2_threshold ",
      "filtering, which can be very expensive. Use a local LD window such as ",
      "max_distance_variants = 1000, or set allow_full_ld = TRUE to explicitly ",
      "permit full pairwise LD evaluation.",
      call. = FALSE
    )
  }

  .Call(
    `_sblr_sparseLD_stream_CSR`,
    bed_files,
    n,
    cls,
    out_prefix,
    rows,
    af,
    pos_bp,
    max_distance_bp,
    max_distance_variants,
    r2_threshold,
    block_size,
    nthreads
  )
}

#' Read a Disk-Backed Sparse-LD CSR Matrix
#'
#' Reads CSR files previously written by [sparseLD_stream_CSR()].
#'
#' @param prefix Prefix of the disk-backed CSR files.
#' @param one_based Convert returned column indices to one-based indexing.
#' @return A list containing CSR row pointers, column indices, and values.
#' @name sparseLD_read_CSR
#' @export
NULL

.resolve_pi_prior <- function(pi_init = 0.001,
                              pi_vb_init = NULL,
                              pi_prior_mean = NULL,
                              pi_prior_strength = NULL,
                              pi_prior_a = NULL,
                              pi_prior_b = NULL) {
  if (is.null(pi_vb_init)) pi_vb_init <- pi_init
  if (is.null(pi_prior_mean)) pi_prior_mean <- pi_init
  if (is.null(pi_prior_strength)) pi_prior_strength <- 2

  for (nm in c("pi_init", "pi_vb_init", "pi_prior_mean")) {
    val <- get(nm)
    if (!is.numeric(val) || length(val) != 1 || !is.finite(val) ||
        val <= 0 || val >= 1) {
      stop(nm, " must be a finite scalar in (0, 1).")
    }
  }

  if (!is.numeric(pi_prior_strength) || length(pi_prior_strength) != 1 ||
      !is.finite(pi_prior_strength) || pi_prior_strength <= 0) {
    stop("pi_prior_strength must be a finite positive scalar.")
  }

  if (is.null(pi_prior_a)) pi_prior_a <- pi_prior_mean * pi_prior_strength
  if (is.null(pi_prior_b)) pi_prior_b <- (1 - pi_prior_mean) * pi_prior_strength

  if (!is.numeric(pi_prior_a) || length(pi_prior_a) != 1 ||
      !is.finite(pi_prior_a) || pi_prior_a <= 0) {
    stop("pi_prior_a must be a finite positive scalar.")
  }

  if (!is.numeric(pi_prior_b) || length(pi_prior_b) != 1 ||
      !is.finite(pi_prior_b) || pi_prior_b <= 0) {
    stop("pi_prior_b must be a finite positive scalar.")
  }

  list(
    pi_init = pi_init,
    pi_vb_init = pi_vb_init,
    pi_prior_mean = pi_prior_mean,
    pi_prior_strength = pi_prior_strength,
    pi_prior_a = pi_prior_a,
    pi_prior_b = pi_prior_b,
    pi = c(1 - pi_init, pi_init)
  )
}

.stblr_prepare_csr_selection_s <- function(selection_s, Glist, m,
                                           scheduled = FALSE,
                                           backend = c("bayesc", "bayesr", "sbayesrc"),
                                           return_log_h = FALSE) {
  backend <- match.arg(backend)
  if (is.null(selection_s) && !isTRUE(return_log_h)) {
    return(list(
      prior_scale = numeric(),
      log_h = numeric(),
      selection_s = NULL,
      fixed = FALSE,
      exponent = NULL
    ))
  }

  if (isTRUE(scheduled)) {
    stop("selection_s is currently supported only for unscheduled CSR BayesC/BayesR/SBayesRC.")
  }
  if (!is.null(selection_s) &&
      (!is.numeric(selection_s) || length(selection_s) != 1L ||
      !is.finite(selection_s))) {
    stop("selection_s must be NULL or a finite numeric scalar.")
  }
  if (is.null(Glist)) {
    stop("Glist with rsidsLD, rsids, and maf is required for selection_s.")
  }
  if (is.null(Glist$rsidsLD) || is.null(Glist$rsids) || is.null(Glist$maf)) {
    stop("Glist must contain rsidsLD, rsids, and maf for selection_s.")
  }

  maf_ld <- unlist(lapply(seq_along(Glist$rsidsLD), function(chr) {
    if (length(Glist$rsids) < chr || is.null(Glist$rsids[[chr]]) ||
        length(Glist$maf) < chr || is.null(Glist$maf[[chr]])) {
      stop("Glist$rsids and Glist$maf must contain each chromosome in Glist$rsidsLD.")
    }
    idx <- match(Glist$rsidsLD[[chr]], Glist$rsids[[chr]])
    if (anyNA(idx)) {
      stop("Could not align MAF to LD marker order for selection_s.")
    }
    Glist$maf[[chr]][idx]
  }), use.names = FALSE)

  if (length(maf_ld) != m) {
    stop("Aligned MAF length must match the number of CSR markers for selection_s.")
  }
  if (!all(is.finite(maf_ld))) {
    stop("Aligned MAF values for selection_s must be finite.")
  }
  if (any(maf_ld <= 0 | maf_ld > 0.5)) {
    stop("Aligned MAF values for selection_s must be in (0, 0.5].")
  }

  h_ld <- pmax(2 * maf_ld * (1 - maf_ld), 1e-8)
  if (length(h_ld) != m || !all(is.finite(h_ld)) || any(h_ld <= 0)) {
    stop("Could not compute positive finite heterozygosity values for selection_s.")
  }

  log_h <- log(h_ld)
  if (is.null(selection_s)) {
    return(list(
      prior_scale = numeric(),
      log_h = log_h,
      selection_s = NULL,
      fixed = FALSE,
      exponent = NULL
    ))
  } else {
    exponent <- selection_s + 1
    prior_scale <- h_ld ^ exponent
    if (!all(is.finite(prior_scale)) || any(prior_scale <= 0)) {
      stop("selection_s prior scale must contain positive finite values.")
    }
  }

  list(
    prior_scale = prior_scale,
    log_h = log_h,
    selection_s = unname(selection_s),
    fixed = TRUE,
    exponent = unname(exponent)
  )
}

.stblr_prepare_csr_bayesc_selection_s <- function(selection_s, Glist, m,
                                                  scheduled = FALSE,
                                                  return_log_h = FALSE) {
  .stblr_prepare_csr_selection_s(
    selection_s = selection_s,
    Glist = Glist,
    m = m,
    scheduled = scheduled,
    backend = "bayesc",
    return_log_h = return_log_h
  )
}

.stblr_validate_sampled_selection_s <- function(estimate_selection_s,
                                                selection_s,
                                                selection_s_init,
                                                selection_s_prior,
                                                selection_s_proposal_sd) {
  if (!is.logical(estimate_selection_s) || length(estimate_selection_s) != 1L ||
      is.na(estimate_selection_s)) {
    stop("estimate_selection_s must be TRUE or FALSE.", call. = FALSE)
  }
  if (!isTRUE(estimate_selection_s)) {
    return(invisible(NULL))
  }
  if (!is.null(selection_s)) {
    stop("selection_s and estimate_selection_s = TRUE cannot both be requested.", call. = FALSE)
  }
  if (!is.numeric(selection_s_init) || length(selection_s_init) != 1L ||
      !is.finite(selection_s_init)) {
    stop("selection_s_init must be a finite numeric scalar.", call. = FALSE)
  }
  if (!is.numeric(selection_s_prior) || length(selection_s_prior) != 2L ||
      anyNA(selection_s_prior) || any(!is.finite(selection_s_prior))) {
    stop("selection_s_prior must be a finite numeric vector of length 2.", call. = FALSE)
  }
  if (selection_s_prior[1L] >= selection_s_prior[2L]) {
    stop("selection_s_prior lower bound must be less than upper bound.", call. = FALSE)
  }
  if (selection_s_init < selection_s_prior[1L] ||
      selection_s_init > selection_s_prior[2L]) {
    stop("selection_s_init must lie within selection_s_prior.", call. = FALSE)
  }
  if (!is.numeric(selection_s_proposal_sd) || length(selection_s_proposal_sd) != 1L ||
      !is.finite(selection_s_proposal_sd) || selection_s_proposal_sd <= 0) {
    stop("selection_s_proposal_sd must be a finite positive numeric scalar.", call. = FALSE)
  }
  invisible(NULL)
}


.make_stblr_priors <- function(y, m, h2, nub, nue, pi_vb_init,
                               pi_prior_mean, trait_names = NULL) {
 y <- as.matrix(y)
 nt <- ncol(y)
 if (is.null(trait_names)) trait_names <- colnames(y)
 if (is.null(trait_names)) trait_names <- paste0("T", seq_len(nt))

 vy <- colSums(y^2) / (nrow(y) - 1)
 B <- diag((vy * h2) / (m * pi_vb_init), nt, nt)
 E <- diag(vy * (1 - h2), nt, nt)
 ssb_prior <- diag(((nub - 2) / nub) * (vy * h2) / (m * pi_prior_mean), nt, nt)
 sse_prior <- diag(((nue - 2) / nue) * (vy * (1 - h2)), nt, nt)

 for (x in c("B", "E", "ssb_prior", "sse_prior")) {
  obj <- get(x)
  rownames(obj) <- colnames(obj) <- trait_names
  assign(x, obj)
 }

 list(
  vy = vy, B = B, E = E, ssb_prior = ssb_prior, sse_prior = sse_prior,
  ssb_prior_list = split(ssb_prior, rep(seq_len(nt), each = nt)),
  sse_prior_list = split(sse_prior, rep(seq_len(nt), each = nt))
 )
}

.stblr_ensure_ld_swap_fields <- function(fit) {
 if (!("ld_swap" %in% names(fit))) fit["ld_swap"] <- list(NULL)
 if (!("ld_swap_chains" %in% names(fit))) fit["ld_swap_chains"] <- list(NULL)
 fit
}

.stblr_ensure_chain_summary_fields <- function(out, raw, meta, marker_mat) {
 fields <- list(
  bm_sd = list(base = "bm", is_sd = TRUE),
  bm_min = list(base = "bm", is_sd = FALSE),
  bm_max = list(base = "bm", is_sd = FALSE),
  dm_sd = list(base = "dm", is_sd = TRUE),
  dm_min = list(base = "dm", is_sd = FALSE),
  dm_max = list(base = "dm", is_sd = FALSE)
 )
 nchains_meta <- suppressWarnings(as.integer(meta$nchains))
 single_chain <- length(nchains_meta) == 1L && !is.na(nchains_meta) &&
  nchains_meta == 1L

 for (nm in names(fields)) {
  raw_val <- raw$marker[[nm]]
  if (!is.null(raw_val)) {
   out[[nm]] <- marker_mat(raw_val)
   next
  }
  base <- out[[fields[[nm]]$base]]
  if (single_chain && !is.null(base)) {
   out[[nm]] <- if (fields[[nm]]$is_sd) base * 0 else base
  } else {
   out[nm] <- list(NULL)
  }
 }
 out
}

.format_stblr_fit <- function(fit, nt, m, trait_names, variable_names,
                              keep_diagnostics = FALSE) {
 nms <- c(
  "bm", "dm", "wy", "r", "b", "d", "o", "vbs", "vgs", "ves",
  "covb", "covg", "cove", "vb", "vg", "ve", "pi", "pim",
  "pitrait", "pimarker"
 )
 has_vle_vld <- length(fit) >= 22 &&
  all(vapply(fit[21:22], function(x) length(unlist(x, use.names = FALSE)) > 0, logical(1)))
 trace_len <- if (length(fit) >= 8 && length(fit[[8]]) > 0L) {
  length(fit[[8]][[1L]])
 } else {
  NA_integer_
 }
 slot23_len <- if (length(fit) >= 23 && length(fit[[23]]) > 0L) {
  length(fit[[23]][[1L]])
 } else {
  NA_integer_
 }
 slot23_ld_swap <- length(fit) >= 23 && slot23_len == 4L &&
  all(vapply(fit[[23]], function(x) {
   length(x) >= 4L && isTRUE(all.equal(as.numeric(x[[4L]]), 1))
  }, logical(1)))
 has_pis <- length(fit) >= 23 && is.finite(trace_len) &&
  slot23_len == trace_len && !slot23_ld_swap
 has_ld_swap <- slot23_ld_swap
 names(fit)[seq_len(min(length(fit), length(nms)))] <-
  nms[seq_len(min(length(fit), length(nms)))]
 if (has_vle_vld) names(fit)[21:22] <- c("vle", "vld")
 if (has_pis) names(fit)[23] <- "pis"
 if (has_ld_swap) names(fit)[23] <- "ld_swap_raw"
 has_chain_summaries <- length(fit) >= 29 &&
  all(vapply(fit[24:29], function(x) length(x) == nt, logical(1)))
 if (has_chain_summaries) {
  names(fit)[24:29] <- c(
   "bm_sd", "bm_min", "bm_max", "dm_sd", "dm_min", "dm_max"
  )
 }
 has_chain_details <- length(fit) >= 32 &&
  all(vapply(fit[30:32], function(x) length(x) == nt, logical(1)))
 if (has_chain_details) {
  names(fit)[30:32] <- c("chain_dm_raw", "chain_bm_raw", "chain_ld_swap_raw")
 }
 selection_slots <- integer()
 if (length(fit) == 25L) {
  selection_slots <- 24:25
 } else if (length(fit) == 31L) {
  selection_slots <- 30:31
 } else if (length(fit) >= 36L) {
  selection_slots <- 33:36
 }
 has_selection_s_trace <- length(selection_slots) >= 2L &&
  length(fit[[selection_slots[1L]]]) == nt &&
  length(fit[[selection_slots[2L]]]) == nt
 if (has_selection_s_trace) {
  names(fit)[selection_slots[1:2]] <- c(
   "selection_s_trace_raw", "selection_s_acceptance_raw"
  )
  if (length(selection_slots) >= 4L) {
   names(fit)[selection_slots[3:4]] <- c(
    "chain_selection_s_raw", "chain_selection_s_acceptance_raw"
   )
  }
 }

 for (i in 1:7) {
  fit[[i]] <- as.matrix(as.data.frame(fit[[i]]))
  rownames(fit[[i]]) <- variable_names
  colnames(fit[[i]]) <- trait_names
 }
 for (i in 8:10) {
  fit[[i]] <- as.matrix(as.data.frame(fit[[i]]))
  rownames(fit[[i]]) <- paste0("Iter", seq_len(nrow(fit[[i]])))
  colnames(fit[[i]]) <- trait_names
 }
 for (i in 11:16) {
  fit[[i]] <- matrix(unlist(fit[[i]]), ncol = nt, byrow = TRUE)
  rownames(fit[[i]]) <- colnames(fit[[i]]) <- trait_names
 }
 fit[[17]] <- fit[[17]][[1]]
 fit[[18]] <- fit[[18]][[1]]

 out <- if (keep_diagnostics) fit[1:20] else fit[1:18]
 if (keep_diagnostics && !is.null(out$pitrait)) {
  out$pitrait <- matrix(unlist(out$pitrait), ncol = 4, byrow = TRUE)
  rownames(out$pitrait) <- trait_names
  colnames(out$pitrait) <- c("log_cpo", "mean_log_cpo", "seconds_mean", "seconds_max")
  out$diagnostics <- out$pitrait
  out$log_cpo <- out$diagnostics[, "log_cpo"]
  out$mean_log_cpo <- out$diagnostics[, "mean_log_cpo"]
 }
 if (keep_diagnostics && !is.null(out$pimarker)) {
  out$pimarker <- matrix(unlist(out$pimarker), ncol = 2, byrow = TRUE)
  rownames(out$pimarker) <- trait_names
  colnames(out$pimarker) <- c("nsamples", "n_used")
 }
 if (has_vle_vld) {
  out$vle <- as.matrix(as.data.frame(fit[[21]]))
  out$vld <- as.matrix(as.data.frame(fit[[22]]))
  rownames(out$vle) <- paste0("Iter", seq_len(nrow(out$vle)))
  rownames(out$vld) <- paste0("Iter", seq_len(nrow(out$vld)))
  colnames(out$vle) <- colnames(out$vld) <- trait_names
 }
 if (has_pis) {
  out$pis <- as.matrix(as.data.frame(fit[[23]]))
  #out$pis <- lapply(fit[[23]], as.numeric)
  #names(out$pis) <- trait_names
  rownames(out$pis) <- paste0("Iter", seq_len(nrow(out$pis)))
  colnames(out$pis) <- trait_names

 }
 if (has_ld_swap) {
  ld_swap <- matrix(unlist(fit[[23]], use.names = FALSE), ncol = 4, byrow = TRUE)
  ld_swap <- ld_swap[, 1:3, drop = FALSE]
  rownames(ld_swap) <- trait_names
  colnames(ld_swap) <- c("attempted", "accepted", "acceptance_rate")
  out$ld_swap <- as.data.frame(ld_swap)
 }
 if (has_chain_summaries) {
  for (nm in c("bm_sd", "bm_min", "bm_max", "dm_sd", "dm_min", "dm_max")) {
   out[[nm]] <- as.matrix(as.data.frame(fit[[nm]]))
   rownames(out[[nm]]) <- variable_names
   colnames(out[[nm]]) <- trait_names
  }
 }
 if (has_chain_details) {
  chains <- setNames(vector("list", nt), trait_names)
  ld_swap_chains <- list()
  for (tt in seq_len(nt)) {
   dm_raw <- as.numeric(fit$chain_dm_raw[[tt]])
   bm_raw <- as.numeric(fit$chain_bm_raw[[tt]])
   diag_raw <- as.numeric(fit$chain_ld_swap_raw[[tt]])
   nchains <- if (m > 0L) length(dm_raw) %/% m else 0L
   trait_chains <- setNames(vector("list", nchains), paste0("chain", seq_len(nchains)))
   trait_ld <- matrix(NA_real_, nrow = nchains, ncol = 3)
   colnames(trait_ld) <- c("attempted", "accepted", "acceptance_rate")
   rownames(trait_ld) <- names(trait_chains)
   for (cc in seq_len(nchains)) {
    idx <- ((cc - 1L) * m + 1L):(cc * m)
    trait_chains[[cc]] <- list(
     dm = stats::setNames(dm_raw[idx], variable_names),
     bm = stats::setNames(bm_raw[idx], variable_names)
    )
    didx <- ((cc - 1L) * 4L + 1L):((cc - 1L) * 4L + 4L)
    if (length(diag_raw) >= max(didx)) {
     trait_chains[[cc]]$ld_swap <- stats::setNames(
      diag_raw[didx][1:3], c("attempted", "accepted", "acceptance_rate")
     )
     trait_ld[cc, ] <- diag_raw[didx][1:3]
    }
   }
   chains[[tt]] <- trait_chains
   ld_swap_chains[[trait_names[tt]]] <- as.data.frame(trait_ld)
  }
  out$chains <- chains
  out$ld_swap_chains <- ld_swap_chains
 }
 if (has_selection_s_trace) {
  out$selection_s_trace <- as.matrix(as.data.frame(fit$selection_s_trace_raw))
  rownames(out$selection_s_trace) <- paste0("Iter", seq_len(nrow(out$selection_s_trace)))
  colnames(out$selection_s_trace) <- trait_names
  out$selection_s_acceptance <- stats::setNames(
   as.numeric(unlist(fit$selection_s_acceptance_raw, use.names = FALSE)),
   trait_names
  )
  if (!is.null(fit$chain_selection_s_raw) && !is.null(out$chains)) {
   for (tt in seq_len(nt)) {
    raw <- as.numeric(fit$chain_selection_s_raw[[tt]])
    acc_raw <- as.numeric(fit$chain_selection_s_acceptance_raw[[tt]])
    if (!length(raw)) next
    nch <- length(out$chains[[tt]])
    trace_n <- if (nch > 0L) length(raw) %/% nch else 0L
    if (nch > 0L && trace_n > 0L) {
     for (cc in seq_len(nch)) {
      idx <- ((cc - 1L) * trace_n + 1L):(cc * trace_n)
      out$chains[[tt]][[cc]]$selection_s <- raw[idx]
      out$chains[[tt]][[cc]]$selection_s_acceptance <-
       if (length(acc_raw) >= cc) acc_raw[cc] else out$selection_s_acceptance[[tt]]
     }
    }
   }
  }
 }
 if (sum(diag(out$covb)) > 0) out$rb <- cov2cor(out$covb)
 if (sum(diag(out$covg)) > 0) out$rg <- cov2cor(out$covg)
 if (sum(diag(out$cove)) > 0) out$re <- cov2cor(out$cove)
 .stblr_ensure_ld_swap_fields(out)
}

.format_stblr_csr_bayesc_fit <- .format_stblr_fit
.format_stblr_bed_bayesc_fit <- .format_stblr_fit

.is_stblr_raw_v1 <- function(raw) {
 is.list(raw) &&
  is.list(raw$schema) &&
  identical(raw$schema$class, "stblr_raw") &&
  identical(as.integer(raw$schema$version), 1L)
}

.format_stblr_raw_v1 <- function(raw, trait_names = NULL, variable_names = NULL) {
 if (!.is_stblr_raw_v1(raw)) {
  stop(".format_stblr_raw_v1() expects an stblr_raw schema version 1 object.")
 }

 meta <- raw$meta
 nt <- as.integer(meta$nt)
 m <- as.integer(meta$m)
 n_trace <- as.integer(meta$n_trace)
 if (is.null(trait_names)) trait_names <- paste0("T", seq_len(nt))
 if (is.null(variable_names)) variable_names <- paste0("V", seq_len(m))

 marker_mat <- function(x) {
  if (is.null(x)) return(NULL)
  x <- as.matrix(x)
  rownames(x) <- variable_names
  colnames(x) <- trait_names
  x
 }
 trace_mat <- function(x) {
  if (is.null(x)) return(NULL)
  x <- as.matrix(x)
  rownames(x) <- paste0("Iter", seq_len(nrow(x)))
  colnames(x) <- trait_names
  x
 }
 variance_mat <- function(x) {
  if (is.null(x)) return(NULL)
  x <- as.matrix(x)
  rownames(x) <- colnames(x) <- trait_names
  x
 }
 named_trait_vec <- function(x) {
  stats::setNames(as.numeric(x), trait_names)
 }

 out <- list(
  bm = marker_mat(raw$marker$bm),
  dm = marker_mat(raw$marker$dm),
  wy = marker_mat(raw$marker$wy),
  r = marker_mat(raw$marker$r),
  b = marker_mat(raw$marker$b),
  d = marker_mat(raw$marker$state),
  vbs = trace_mat(raw$trace$vbs),
  vgs = trace_mat(raw$trace$vgs),
  ves = trace_mat(raw$trace$ves),
  vle = trace_mat(raw$trace$vle),
  vld = trace_mat(raw$trace$vld),
  pis = trace_mat(raw$trace$pis),
  covb = variance_mat(raw$variance$covb),
  covg = variance_mat(raw$variance$covg),
  cove = variance_mat(raw$variance$cove),
  vb = variance_mat(raw$variance$vb),
  vg = variance_mat(raw$variance$vg),
  ve = variance_mat(raw$variance$ve)
 )
 if (identical(meta$model, "bayesr")) {
  out$component <- out$d
 } else if (identical(meta$model, "sbayesrc")) {
  out$component <- out$d
 }

 pi_final <- as.matrix(raw$pi$final)
 pi_mean <- as.matrix(raw$pi$mean)
 pi_names <- raw$pi$names %||% paste0("component_", seq_len(ncol(pi_final)) - 1L)
 colnames(pi_final) <- colnames(pi_mean) <- pi_names
 rownames(pi_final) <- rownames(pi_mean) <- trait_names
 if (identical(meta$model, "bayesr") || identical(meta$model, "sbayesrc")) {
  out$pi <- pi_final
  out$pim <- pi_mean
  out$final_pi <- pi_final
  out$mean_pi <- pi_mean
 } else {
  out$pi <- if (nt == 1L) as.numeric(pi_final[1L, ]) else pi_final
  out$pim <- if (nt == 1L) as.numeric(pi_mean[1L, ]) else pi_mean
 }

 out <- .stblr_ensure_chain_summary_fields(out, raw, meta, marker_mat)

 if (identical(meta$model, "bayesr") || identical(meta$model, "sbayesrc")) {
  component_names <- raw$component$names %||% pi_names
  comp_prob <- raw$component$prob
  if (is.list(comp_prob) && length(comp_prob) == nt) {
   comp_prob <- lapply(seq_len(nt), function(tt) {
    x <- as.matrix(comp_prob[[tt]])
    rownames(x) <- variable_names
    colnames(x) <- component_names
    x
   })
   names(comp_prob) <- trait_names
   out$comp_prob <- comp_prob
   null_component <- if (identical(meta$model, "sbayesrc")) {
    "gamma_0.00"
   } else {
    "component_0"
   }
   if (!null_component %in% colnames(comp_prob[[1L]])) {
    stop("Raw ", meta$model, " component probabilities are missing null component ", null_component, ".")
   }
   dm_from_components <- vapply(
    comp_prob,
    function(x) 1 - x[, null_component],
    numeric(m)
   )
   if (nt == 1L) dm_from_components <- matrix(dm_from_components, ncol = 1L)
   rownames(dm_from_components) <- variable_names
   colnames(dm_from_components) <- trait_names
   out$dm <- dm_from_components
  }
  if (!is.null(raw$component$dm_component_mean)) {
   out$dm_component_mean <- marker_mat(raw$component$dm_component_mean)
  } else if (!is.null(out$comp_prob)) {
   component_index <- seq_along(component_names) - 1L
   dm_component_mean <- vapply(
    out$comp_prob,
    function(x) as.numeric(x %*% component_index),
    numeric(m)
   )
   if (nt == 1L) dm_component_mean <- matrix(dm_component_mean, ncol = 1L)
   rownames(dm_component_mean) <- variable_names
   colnames(dm_component_mean) <- trait_names
   out$dm_component_mean <- dm_component_mean
  }
  if (!is.null(raw$component$ncomp)) {
   ncomp <- as.matrix(raw$component$ncomp)
   rownames(ncomp) <- trait_names
   colnames(ncomp) <- component_names
   out$ncomp <- ncomp
  }
  if (!is.null(raw$component$mixture_var)) {
   out$mixture_var <- stats::setNames(
    as.numeric(raw$component$mixture_var),
    component_names
   )
  }
  if (identical(meta$model, "sbayesrc")) {
   out$gamma <- out$mixture_var
  }
 }

 if (identical(meta$model, "sbayesrc")) {
  annotation_names <- raw$annotation$annotation_names
  if (is.null(annotation_names)) {
   annotation_names <- paste0("A", seq_len(as.integer(meta$n_annotations)))
  }
  step_names <- paste0("step_", seq_len(max(0L, as.integer(meta$n_components) - 1L)))
  alpha_raw <- raw$annotation$alpha_mean %||% raw$annotation$alpha
  if (is.list(alpha_raw) && length(alpha_raw) == nt) {
   alpha <- lapply(seq_len(nt), function(tt) {
    x <- as.matrix(alpha_raw[[tt]])
    rownames(x) <- annotation_names
    colnames(x) <- step_names
    x
   })
   names(alpha) <- trait_names
   out$alpha <- alpha
   alpha_flat <- do.call(rbind, lapply(alpha, function(x) as.numeric(x)))
   rownames(alpha_flat) <- trait_names
   colnames(alpha_flat) <- as.vector(outer(annotation_names, step_names, paste, sep = ":"))
   out$alpha_flat <- alpha_flat
  }
  sigma_raw <- raw$annotation$sigmaSqAlpha_mean %||% raw$annotation$sigmaSqAlpha
  if (!is.null(sigma_raw)) {
   sigma <- as.matrix(sigma_raw)
   if (identical(dim(sigma), c(length(step_names), nt))) {
    sigma <- t(sigma)
   }
   rownames(sigma) <- trait_names
   colnames(sigma) <- step_names
   out$sigmaSqAlpha <- sigma
  }
 }

 if (identical(meta$backend, "csr_prior_bayesc")) {
  out$prior <- raw$prior
 }

 if (identical(meta$backend, "csr_group_bayesc")) {
  group_names <- raw$group$group_names %||%
   paste0("G", seq_len(as.integer(meta$n_groups)))
  group_mat <- function(x) {
   x <- as.matrix(x)
   if (identical(dim(x), c(length(group_names), nt))) x <- t(x)
   rownames(x) <- trait_names
   colnames(x) <- group_names
   x
  }
  out$group <- raw$group
  if (!is.null(raw$group$pi_mean)) {
   out$group_pi <- group_mat(raw$group$pi_mean)
  }
  if (!is.null(raw$group$vb_multiplier_mean)) {
   out$group_vb_multiplier <- group_mat(raw$group$vb_multiplier_mean)
  }
  if (!is.null(raw$group$n_included_mean)) {
   out$group_nincluded <- group_mat(raw$group$n_included_mean)
  }
  if (!is.null(raw$group$size)) {
   out$group_size <- matrix(
    rep(as.numeric(raw$group$size), each = nt),
    nrow = nt,
    ncol = length(group_names),
    byrow = FALSE
   )
   rownames(out$group_size) <- trait_names
   colnames(out$group_size) <- group_names
  }
 }

 if (identical(meta$backend, "csr_annot_bayesc")) {
  annotation_names <- raw$annotation$annotation_names %||%
   paste0("A", seq_len(as.integer(meta$n_annotations)))
  annotation_mat <- function(x) {
   x <- as.matrix(x)
   if (identical(dim(x), c(length(annotation_names), nt))) x <- t(x)
   rownames(x) <- trait_names
   colnames(x) <- annotation_names
   x
  }
  out$annotation <- raw$annotation
  if (!is.null(raw$annotation$eta_pi_mean)) {
   out$eta_pi <- annotation_mat(raw$annotation$eta_pi_mean)
  }
  if (!is.null(raw$annotation$eta_vb_mean)) {
   out$eta_vb <- annotation_mat(raw$annotation$eta_vb_mean)
  }
 }

 set_updateE_diagnostics <- function(x) {
  x <- as.matrix(x)
  if (ncol(x) == 8L) {
   colnames(x) <- c(
    "trait_index", "chain_index", "n_updateE", "min_sse",
    "min_residual_scale", "max_nonzero_components",
    "max_abs_effect", "max_fitted_quadratic"
   )
  } else if (ncol(x) == 9L) {
   colnames(x) <- c(
    "trait_index", "chain_index", "n_updateE", "min_sse",
    "min_sse_iter", "min_residual_scale", "max_nonzero_components",
    "max_abs_effect", "max_fitted_quadratic"
   )
  }
  x
 }
 if (!is.null(raw$diagnostics$updateE)) {
  out$updateE_diagnostics <- set_updateE_diagnostics(raw$diagnostics$updateE)
 }
 if (!is.null(raw$diagnostics$log_cpo)) {
  out$log_cpo <- named_trait_vec(raw$diagnostics$log_cpo)
 }
 if (!is.null(raw$diagnostics$mean_log_cpo)) {
  out$mean_log_cpo <- named_trait_vec(raw$diagnostics$mean_log_cpo)
 }
 if (!is.null(raw$diagnostics$nsamples)) {
  out$nsamples <- named_trait_vec(raw$diagnostics$nsamples)
 }
 if (!is.null(raw$diagnostics$n_used)) {
  out$n_used <- named_trait_vec(raw$diagnostics$n_used)
 }
 out$diagnostics <- raw$diagnostics

 ld_swap <- raw$diagnostics$ld_swap
 if (!is.null(ld_swap)) {
  ld_swap <- as.matrix(ld_swap)
  rownames(ld_swap) <- trait_names
  colnames(ld_swap) <- c("attempted", "accepted", "acceptance_rate")
  out$ld_swap <- as.data.frame(ld_swap)
 }

 has_chains <- !is.null(raw$chains) &&
  length(raw$chains) > 0L &&
  !isTRUE(identical(raw$chains, list())) &&
  any(!vapply(raw$chains, is.null, logical(1)))

 if (has_chains) {
  chains <- setNames(vector("list", nt), trait_names)
  ld_swap_chains <- list()
  for (tt in seq_len(nt)) {
   if (tt > length(raw$chains) || is.null(raw$chains[[tt]])) {
    next
   }
   trait_raw <- raw$chains[[tt]]
   nch <- length(trait_raw)
   if (nch < 1L) {
    next
   }
   trait_chains <- setNames(vector("list", nch), paste0("chain", seq_len(nch)))
   trait_ld <- matrix(NA_real_, nrow = nch, ncol = 3)
   rownames(trait_ld) <- names(trait_chains)
   colnames(trait_ld) <- c("attempted", "accepted", "acceptance_rate")
   for (cc in seq_len(nch)) {
    ch <- trait_raw[[cc]]
    if (is.null(ch)) {
     next
    }
    trait_chains[[cc]] <- list(
     dm = stats::setNames(as.numeric(ch$marker$dm), variable_names),
     bm = stats::setNames(as.numeric(ch$marker$bm), variable_names)
    )
    if (!is.null(ch$marker$state)) {
     trait_chains[[cc]]$state <- stats::setNames(
      as.numeric(ch$marker$state), variable_names
     )
    }
    if (identical(meta$model, "bayesr") || identical(meta$model, "sbayesrc")) {
     trait_chains[[cc]]$component <- trait_chains[[cc]]$state
     if (!is.null(ch$component$prob)) {
      cp <- as.matrix(ch$component$prob)
      rownames(cp) <- variable_names
      colnames(cp) <- raw$component$names %||% pi_names
      trait_chains[[cc]]$comp_prob <- cp
     }
     if (!is.null(ch$component$dm_component_mean)) {
      trait_chains[[cc]]$dm_component_mean <- stats::setNames(
       as.numeric(ch$component$dm_component_mean), variable_names
      )
     }
     if (!is.null(ch$pi$final)) {
      trait_chains[[cc]]$final_pi <- stats::setNames(
       as.numeric(ch$pi$final), raw$component$names %||% pi_names
      )
     }
     if (!is.null(ch$pi$mean)) {
      trait_chains[[cc]]$mean_pi <- stats::setNames(
       as.numeric(ch$pi$mean), raw$component$names %||% pi_names
      )
     }
     if (!is.null(ch$component$ncomp)) {
      trait_chains[[cc]]$ncomp <- stats::setNames(
       as.numeric(ch$component$ncomp), raw$component$names %||% pi_names
      )
     }
     if (identical(meta$model, "sbayesrc") && !is.null(ch$annotation$alpha)) {
      alpha <- as.matrix(ch$annotation$alpha)
      annotation_names <- raw$annotation$annotation_names %||%
       paste0("A", seq_len(as.integer(meta$n_annotations)))
      step_names <- paste0("step_", seq_len(max(0L, as.integer(meta$n_components) - 1L)))
      rownames(alpha) <- annotation_names
      colnames(alpha) <- step_names
      trait_chains[[cc]]$alpha <- alpha
     }
     if (identical(meta$model, "sbayesrc") && !is.null(ch$annotation$sigmaSqAlpha)) {
      step_names <- paste0("step_", seq_len(max(0L, as.integer(meta$n_components) - 1L)))
      trait_chains[[cc]]$sigmaSqAlpha <- stats::setNames(
       as.numeric(ch$annotation$sigmaSqAlpha), step_names
      )
     }
     for (nm in c("vbs", "vgs", "ves", "vle", "vld", "pis")) {
      if (!is.null(ch$trace[[nm]])) {
       trait_chains[[cc]][[nm]] <- as.numeric(ch$trace[[nm]])
       names(trait_chains[[cc]][[nm]]) <- paste0(
        "Iter", seq_along(trait_chains[[cc]][[nm]])
       )
      }
     }
     if (!is.null(ch$diagnostics$updateE)) {
      trait_chains[[cc]]$updateE_diagnostics <- set_updateE_diagnostics(
       matrix(as.numeric(ch$diagnostics$updateE), nrow = 1L)
      )
     }
    }
    if (identical(meta$backend, "csr_group_bayesc") && !is.null(ch$group)) {
     group_names <- raw$group$group_names %||%
      paste0("G", seq_len(as.integer(meta$n_groups)))
     if (!is.null(ch$group$pi_mean)) {
      trait_chains[[cc]]$group_pi <- stats::setNames(
       as.numeric(ch$group$pi_mean), group_names
      )
     }
     if (!is.null(ch$group$vb_multiplier_mean)) {
      trait_chains[[cc]]$group_vb_multiplier <- stats::setNames(
       as.numeric(ch$group$vb_multiplier_mean), group_names
      )
     }
     if (!is.null(ch$group$n_included_mean)) {
      trait_chains[[cc]]$group_nincluded <- stats::setNames(
       as.numeric(ch$group$n_included_mean), group_names
      )
     }
    }
    if (identical(meta$backend, "csr_annot_bayesc") && !is.null(ch$annotation)) {
     annotation_names <- raw$annotation$annotation_names %||%
      paste0("A", seq_len(as.integer(meta$n_annotations)))
     if (!is.null(ch$annotation$eta_pi_mean)) {
      trait_chains[[cc]]$eta_pi <- stats::setNames(
       as.numeric(ch$annotation$eta_pi_mean), annotation_names
      )
     }
     if (!is.null(ch$annotation$eta_vb_mean)) {
      trait_chains[[cc]]$eta_vb <- stats::setNames(
       as.numeric(ch$annotation$eta_vb_mean), annotation_names
      )
     }
    }
    if (!is.null(ch$diagnostics$ld_swap)) {
     ld <- as.numeric(ch$diagnostics$ld_swap[1L, seq_len(3L)])
     trait_chains[[cc]]$ld_swap <- stats::setNames(
      ld, c("attempted", "accepted", "acceptance_rate")
     )
     trait_ld[cc, ] <- ld
    }
    if (!is.null(ch$selection$trace)) {
     trait_chains[[cc]]$selection_s <- as.numeric(ch$selection$trace)
     trait_chains[[cc]]$selection_s_acceptance <-
      as.numeric(ch$selection$acceptance)
    }
   }
   chains[[tt]] <- trait_chains
   if (any(!is.na(trait_ld))) {
    ld_swap_chains[[trait_names[tt]]] <- as.data.frame(trait_ld)
   }
  }
  out$chains <- chains
  if (length(ld_swap_chains)) out$ld_swap_chains <- ld_swap_chains
 } else {
  out$chains <- NULL
  out$ld_swap_chains <- NULL
 }

 if (!is.null(raw$selection$trace)) {
  out$selection_s_trace <- trace_mat(raw$selection$trace)
  out$selection_s <- named_trait_vec(raw$selection$mean)
  out$selection_s_sd <- named_trait_vec(raw$selection$sd)
  out$selection_s_min <- named_trait_vec(raw$selection$min)
  out$selection_s_max <- named_trait_vec(raw$selection$max)
  out$selection_s_acceptance <- named_trait_vec(raw$selection$acceptance)
 } else if (!is.null(raw$selection$mean)) {
  out$selection_s <- named_trait_vec(raw$selection$mean)
 }

 if (sum(diag(out$covb)) > 0) out$rb <- cov2cor(out$covb)
 if (sum(diag(out$covg)) > 0) out$rg <- cov2cor(out$covg)
 if (sum(diag(out$cove)) > 0) out$re <- cov2cor(out$cove)
 .stblr_ensure_ld_swap_fields(out)
}

.format_stblr_bayesr_fit <- function(fit, nt, m, trait_names, variable_names,
                                     n_components,
                                     keep_diagnostics = FALSE) {
 if (!is.list(fit) || length(fit) < 23L) {
  stop(".format_stblr_bayesr_fit() expects a BayesR raw return list.")
 }
 if (!is.numeric(n_components) || length(n_components) != 1L ||
     !is.finite(n_components) || n_components < 2L ||
     n_components != floor(n_components)) {
  stop("n_components must be a finite integer scalar of at least 2.")
 }
 n_components <- as.integer(n_components)

 comp_prob_raw <- fit[[23L]]
 if (!is.list(comp_prob_raw) || length(comp_prob_raw) != nt ||
     any(lengths(comp_prob_raw) != m * n_components)) {
  stop("BayesR component-probability slot must have length n_components * m per trait.")
 }
 component_mean_raw <- if (length(fit) >= 30L) fit[[30L]] else NULL

 out <- .format_stblr_fit(
  fit, nt, m, trait_names, variable_names,
  keep_diagnostics = keep_diagnostics
 )

 component_names <- paste0("component_", seq.int(0L, n_components - 1L))

 pi_raw <- fit[[17L]]
 pim_raw <- fit[[18L]]
 if (!is.list(pi_raw) || length(pi_raw) != nt ||
     any(lengths(pi_raw) != n_components) ||
     !is.list(pim_raw) || length(pim_raw) != nt ||
     any(lengths(pim_raw) != n_components)) {
  stop("BayesR pi and pim slots must have length n_components per trait.")
 }
 out$pi <- t(vapply(pi_raw, as.numeric, numeric(n_components)))
 out$pim <- t(vapply(pim_raw, as.numeric, numeric(n_components)))
 rownames(out$pi) <- rownames(out$pim) <- trait_names
 colnames(out$pi) <- colnames(out$pim) <- component_names
 out$final_pi <- out$pi
 out$mean_pi <- out$pim

 comp_prob <- lapply(seq_len(nt), function(tt) {
  raw <- as.numeric(comp_prob_raw[[tt]])
  mat <- t(matrix(raw, nrow = n_components, ncol = m, byrow = TRUE))
  rownames(mat) <- variable_names
  colnames(mat) <- component_names
  mat
 })
 names(comp_prob) <- trait_names
 out$comp_prob <- comp_prob

 dm_from_components <- vapply(
  comp_prob,
  function(x) 1 - x[, 1L],
  numeric(m)
 )
 if (nt == 1L) {
  dm_from_components <- matrix(dm_from_components, ncol = 1L)
 }
 rownames(dm_from_components) <- variable_names
 colnames(dm_from_components) <- trait_names
 out$dm <- dm_from_components

 if (!is.null(component_mean_raw)) {
  if (!is.list(component_mean_raw) || length(component_mean_raw) != nt ||
      any(lengths(component_mean_raw) != m)) {
   stop("BayesR component-mean slot must have length m per trait.")
  }
  out$dm_component_mean <- as.matrix(as.data.frame(component_mean_raw))
  rownames(out$dm_component_mean) <- variable_names
  colnames(out$dm_component_mean) <- trait_names
 }

 tol <- 1e-8
 valid_prob <- vapply(comp_prob, function(x) {
  all(is.finite(x)) &&
   all(x >= -tol & x <= 1 + tol) &&
   all(abs(rowSums(x) - 1) <= 1e-6)
 }, logical(1))
 if (!all(valid_prob)) {
  stop("BayesR component probabilities must be finite, within [0, 1], and sum to 1 by marker.")
 }

 out
}

.format_stblr_bed_bayesr_fit <- .format_stblr_bayesr_fit

.fit_stblr_bed_bayesr <- function(
  bed_files,
  n,
  cls,
  y,
  b_init,
  sets,
  rows = NULL,
  af = NULL,
  scale = TRUE,
  B,
  E,
  ssb_prior,
  sse_prior,
  pi,
  mixture_var,
  alpha,
  nub = 4,
  nue = 4,
  updateB = TRUE,
  updateE = TRUE,
  updatePi = TRUE,
  adjE = 0,
  nit = 1000,
  nburn = 100,
  nthin = 1,
  rebuild_every = 25,
  full_sweep_every = 10,
  null_skip_base = 50,
  null_skip_max = 200,
  candidate_threshold = 1e-3,
  candidate_lifetime = 20,
  skip_nulls_burnin_only = FALSE,
  return_wy = FALSE,
  return_r = FALSE,
  read_block_size = 64,
  progress_every = 0,
  nchains = 1,
  ncores = 1,
  seed = 10,
  trait_names = NULL,
  variable_names = NULL,
  keep_diagnostics = TRUE
) {
 if (!is.numeric(nchains) || length(nchains) != 1L ||
     !is.finite(nchains) || nchains < 1 || nchains != floor(nchains)) {
  stop("nchains must be a positive integer scalar.")
 }
 nchains <- as.integer(nchains)
 if (!is.numeric(ncores) || length(ncores) != 1L ||
     !is.finite(ncores) || ncores < 1 || ncores != floor(ncores)) {
  stop("ncores must be a positive integer scalar.")
 }
 ncores <- as.integer(ncores)
 if (!is.numeric(mixture_var) || length(mixture_var) < 2L ||
     any(!is.finite(mixture_var))) {
  stop("mixture_var must be a finite numeric vector with at least two components.")
 }
 if (mixture_var[1L] != 0) stop("mixture_var[1] must be 0 for the null component.")
 if (any(mixture_var[-1L] <= 0)) {
  stop("Non-null mixture_var values must be positive.")
 }
 if (!is.numeric(pi) || length(pi) != length(mixture_var) ||
     any(!is.finite(pi)) || any(pi < 0) || sum(pi) <= 0) {
  stop("pi must be a non-negative finite vector matching mixture_var with positive sum.")
 }
 if (!is.numeric(alpha) || length(alpha) != length(mixture_var) ||
     any(!is.finite(alpha)) || any(alpha <= 0)) {
  stop("alpha must be a positive finite vector matching mixture_var.")
 }

 y_mat <- as.matrix(y)
 nt <- ncol(y_mat)
 m <- length(unlist(cls, use.names = FALSE))
 if (is.null(trait_names)) {
  trait_names <- colnames(y_mat)
  if (is.null(trait_names)) trait_names <- paste0("D", seq_len(nt))
 }
 if (is.null(variable_names)) variable_names <- paste0("m", seq_len(m))
 if (length(trait_names) != nt) stop("trait_names must have length ncol(y).")
 if (length(variable_names) != m) {
  stop("variable_names must match the number of markers implied by cls.")
 }

 raw <- stblr_cpg_omp_bed_marker_scheduled_chains_bayesr(
  bed_files = bed_files,
  n = as.integer(n),
  cls = cls,
  y = y_mat,
  b_init = b_init,
  sets = as.integer(sets),
  rows = rows,
  af = af,
  scale = scale,
  B = B,
  E = E,
  ssb_prior = ssb_prior,
  sse_prior = sse_prior,
  pi = as.numeric(pi),
  c = as.numeric(mixture_var),
  alpha = as.numeric(alpha),
  nub = nub,
  nue = nue,
  updateB = updateB,
  updateE = updateE,
  updatePi = updatePi,
  adjE = adjE,
  nit = as.integer(nit),
  nburn = as.integer(nburn),
  nthin = as.integer(nthin),
  rebuild_every = as.integer(rebuild_every),
  full_sweep_every = as.integer(full_sweep_every),
  null_skip_base = as.integer(null_skip_base),
  null_skip_max = as.integer(null_skip_max),
  candidate_threshold = candidate_threshold,
  candidate_lifetime = as.integer(candidate_lifetime),
  skip_nulls_burnin_only = skip_nulls_burnin_only,
  return_wy = return_wy,
  return_r = return_r,
  read_block_size = as.integer(read_block_size),
  progress_every = as.integer(progress_every),
  nchains = nchains,
  ncores = ncores,
  seed = as.integer(seed)
 )

 if (.is_stblr_raw_v1(raw)) {
  fit <- .format_stblr_raw_v1(raw, trait_names, variable_names)
 } else {
  fit <- .format_stblr_bed_bayesr_fit(
   raw,
   nt = nt,
   m = m,
   trait_names = trait_names,
   variable_names = variable_names,
   n_components = length(mixture_var),
   keep_diagnostics = keep_diagnostics
  )
 }
 fit$input <- list(
  method = "bayesr",
  model = "bayesr",
  backend = "bed_bayesr",
  data_level = "individual",
  annotations = FALSE,
  scheduled = TRUE,
  keep_chains = FALSE,
  nchains = nchains,
  ncores = ncores,
  seed = as.integer(seed),
  n = as.integer(n),
  m = m,
  nt = nt,
  mixture_var = as.numeric(mixture_var),
  pi = as.numeric(pi),
  alpha = as.numeric(alpha),
  nit = as.integer(nit),
  nburn = as.integer(nburn),
  nthin = as.integer(nthin)
 )
 fit
}

# Compatibility alias for scripts/tests that used the old internal helper name.
.stblr_bed_marker_bayesr_experimental <- function(...) {
 .fit_stblr_bed_bayesr(...)
}

.format_stblr_csr_bayesr_fit <- function(fit, nt, m, trait_names,
                                         variable_names, n_components) {
 if (.is_stblr_raw_v1(fit)) {
  return(.format_stblr_raw_v1(fit, trait_names, variable_names))
 }
 if (!is.list(fit) || is.null(names(fit))) {
  stop(".format_stblr_csr_bayesr_fit() expects a named CSR BayesR return list.")
 }
 if (!is.numeric(n_components) || length(n_components) != 1L ||
     !is.finite(n_components) || n_components < 2L ||
     n_components != floor(n_components)) {
  stop("n_components must be a finite integer scalar of at least 2.")
 }
 n_components <- as.integer(n_components)
 component_names <- paste0("component_", seq.int(0L, n_components - 1L))

 set_marker_matrix <- function(x, name) {
  if (is.null(x)) stop("CSR BayesR return is missing field: ", name)
  x <- as.matrix(x)
  if (!identical(dim(x), c(m, nt))) {
   stop("CSR BayesR field ", name, " must have dimension m x nt.")
  }
  rownames(x) <- variable_names
  colnames(x) <- trait_names
  x
 }

 set_trace_matrix <- function(x, name) {
  if (is.null(x)) return(NULL)
  x <- as.matrix(x)
  if (ncol(x) != nt) {
   stop("CSR BayesR field ", name, " must have nt columns.")
  }
  rownames(x) <- paste0("Iter", seq_len(nrow(x)))
  colnames(x) <- trait_names
  x
 }

 set_trait_square <- function(x, name) {
  if (is.null(x)) return(matrix(NA_real_, nt, nt,
                               dimnames = list(trait_names, trait_names)))
  x <- as.matrix(x)
  if (!identical(dim(x), c(nt, nt))) {
   stop("CSR BayesR field ", name, " must have dimension nt x nt.")
  }
  rownames(x) <- colnames(x) <- trait_names
  x
 }

 out <- list(
  bm = set_marker_matrix(fit$bm, "bm"),
  dm = set_marker_matrix(fit$dm, "dm"),
  wy = set_marker_matrix(fit$wy, "wy"),
  r = set_marker_matrix(fit$r, "r"),
  b = set_marker_matrix(fit$b, "b"),
  component = set_marker_matrix(fit$component, "component"),
  vbs = set_trace_matrix(fit$vbs, "vbs"),
  vgs = set_trace_matrix(fit$vgs, "vgs"),
  ves = set_trace_matrix(fit$ves, "ves"),
  covb = set_trait_square(fit$covb, "covb"),
  covg = set_trait_square(fit$covg, "covg"),
  cove = set_trait_square(fit$cove, "cove"),
  vb = set_trait_square(fit$vb, "vb"),
  vg = set_trait_square(fit$vg, "vg"),
  ve = set_trait_square(fit$ve, "ve"),
  pi = as.matrix(fit$pi),
  pim = as.matrix(fit$pim),
  vle = set_trace_matrix(fit$vle, "vle"),
  vld = set_trace_matrix(fit$vld, "vld"),
  bm_sd = set_marker_matrix(fit$bm_sd, "bm_sd"),
  bm_min = set_marker_matrix(fit$bm_min, "bm_min"),
  bm_max = set_marker_matrix(fit$bm_max, "bm_max"),
  dm_sd = set_marker_matrix(fit$dm_sd, "dm_sd"),
  dm_min = set_marker_matrix(fit$dm_min, "dm_min"),
  dm_max = set_marker_matrix(fit$dm_max, "dm_max"),
  dm_component_mean = set_marker_matrix(
   fit$dm_component_mean, "dm_component_mean"
  )
 )

 if (!identical(dim(out$pi), c(nt, n_components))) {
  stop("CSR BayesR field pi must have dimension nt x n_components.")
 }
 if (!identical(dim(out$pim), c(nt, n_components))) {
  stop("CSR BayesR field pim must have dimension nt x n_components.")
 }
 rownames(out$pi) <- rownames(out$pim) <- trait_names
 colnames(out$pi) <- colnames(out$pim) <- component_names

 comp_prob <- fit$comp_prob
 if (!is.list(comp_prob) || length(comp_prob) != nt) {
  stop("CSR BayesR comp_prob must be a list with one matrix per trait.")
 }
 comp_prob <- lapply(seq_len(nt), function(tt) {
  x <- as.matrix(comp_prob[[tt]])
  if (!identical(dim(x), c(m, n_components))) {
   stop("CSR BayesR comp_prob entries must have dimension m x n_components.")
  }
  rownames(x) <- variable_names
  colnames(x) <- component_names
  x
 })
 names(comp_prob) <- trait_names
 out$comp_prob <- comp_prob

 dm_from_components <- vapply(
  comp_prob,
  function(x) 1 - x[, 1L],
  numeric(m)
 )
 if (nt == 1L) dm_from_components <- matrix(dm_from_components, ncol = 1L)
 rownames(dm_from_components) <- variable_names
 colnames(dm_from_components) <- trait_names
 out$dm <- dm_from_components

 if (!is.null(fit$ncomp)) {
  ncomp <- as.matrix(fit$ncomp)
  if (!identical(dim(ncomp), c(nt, n_components))) {
   stop("CSR BayesR ncomp must have dimension nt x n_components.")
  }
  rownames(ncomp) <- trait_names
  colnames(ncomp) <- component_names
  out$ncomp <- ncomp
 }
 if (!is.null(fit$mixture_var)) {
  out$mixture_var <- stats::setNames(as.numeric(fit$mixture_var), component_names)
 }
 set_updateE_diagnostics <- function(x) {
  x <- as.matrix(x)
  if (ncol(x) == 8L) {
   colnames(x) <- c(
    "trait_index", "chain_index", "n_updateE", "min_sse",
    "min_residual_scale", "max_nonzero_components",
    "max_abs_effect", "max_fitted_quadratic"
   )
  } else if (ncol(x) == 9L) {
   colnames(x) <- c(
    "trait_index", "chain_index", "n_updateE", "min_sse",
    "min_sse_iter", "min_residual_scale", "max_nonzero_components",
    "max_abs_effect", "max_fitted_quadratic"
   )
  }
  x
 }

 if (!is.null(fit$updateE_diagnostics)) {
  out$updateE_diagnostics <- set_updateE_diagnostics(fit$updateE_diagnostics)
 }

 if (!is.null(fit$ld_swap)) {
  ld_swap <- as.matrix(fit$ld_swap)
  if (!identical(dim(ld_swap), c(nt, 3L))) {
   stop("CSR BayesR ld_swap must have dimension nt x 3.")
  }
  rownames(ld_swap) <- trait_names
  colnames(ld_swap) <- c("attempted", "accepted", "acceptance_rate")
  out$ld_swap <- as.data.frame(ld_swap)
 }

 if (!is.null(fit$chains)) {
  if (!is.list(fit$chains) || length(fit$chains) != nt) {
   stop("CSR BayesR chains must be a list with one entry per trait.")
  }
  ld_swap_chain_raw <- fit$ld_swap_chains
  if (!is.null(ld_swap_chain_raw)) {
   ld_swap_chain_raw <- as.matrix(ld_swap_chain_raw)
   if (ncol(ld_swap_chain_raw) != 5L) {
    stop("CSR BayesR ld_swap_chains must have five columns.")
   }
  }
  chains <- lapply(seq_len(nt), function(tt) {
   trait_chains <- fit$chains[[tt]]
   if (!is.list(trait_chains)) {
    stop("CSR BayesR chains entries must be lists of chains.")
   }
   trait_chains <- lapply(seq_along(trait_chains), function(cc) {
    ch <- trait_chains[[cc]]
    if (!is.list(ch)) stop("CSR BayesR chain entry must be a list.")
    chain <- list(
     dm = stats::setNames(as.numeric(ch$dm), variable_names),
     bm = stats::setNames(as.numeric(ch$bm), variable_names),
     comp_prob = as.matrix(ch$comp_prob),
     dm_component_mean = stats::setNames(
      as.numeric(ch$dm_component_mean), variable_names
     ),
     final_pi = stats::setNames(as.numeric(ch$final_pi), component_names),
     mean_pi = stats::setNames(as.numeric(ch$mean_pi), component_names)
    )
    if (!identical(dim(chain$comp_prob), c(m, n_components))) {
     stop("CSR BayesR chain comp_prob entries must have dimension m x n_components.")
    }
    rownames(chain$comp_prob) <- variable_names
    colnames(chain$comp_prob) <- component_names
    for (nm in c("vbs", "vgs", "ves", "vle", "vld")) {
     if (!is.null(ch[[nm]])) {
      chain[[nm]] <- as.numeric(ch[[nm]])
      names(chain[[nm]]) <- paste0("Iter", seq_along(chain[[nm]]))
     }
    }
    if (!is.null(ch$updateE_diagnostics)) {
     chain$updateE_diagnostics <- set_updateE_diagnostics(
      matrix(as.numeric(ch$updateE_diagnostics), nrow = 1L)
     )
    }
    if (!is.null(ch$ld_swap)) {
     chain$ld_swap <- stats::setNames(
      as.numeric(ch$ld_swap)[1:3],
      c("attempted", "accepted", "acceptance_rate")
     )
    }
    if (!is.null(ch$selection_s)) {
     chain$selection_s <- as.numeric(ch$selection_s)
     names(chain$selection_s) <- paste0("Iter", seq_along(chain$selection_s))
    }
    if (!is.null(ch$selection_s_acceptance)) {
     chain$selection_s_acceptance <- as.numeric(ch$selection_s_acceptance)[1L]
    }
    chain
   })
   names(trait_chains) <- paste0("chain", seq_along(trait_chains))
   trait_chains
  })
  names(chains) <- trait_names
  out$chains <- chains
  if (!is.null(ld_swap_chain_raw)) {
   ld_swap_chains <- lapply(seq_len(nt), function(tt) {
    rows <- which(ld_swap_chain_raw[, 1L] == (tt - 1L))
    if (!length(rows)) return(NULL)
    trait_ld <- ld_swap_chain_raw[rows, 3:5, drop = FALSE]
    rownames(trait_ld) <- paste0("chain", seq_len(nrow(trait_ld)))
    colnames(trait_ld) <- c("attempted", "accepted", "acceptance_rate")
    as.data.frame(trait_ld)
   })
   names(ld_swap_chains) <- trait_names
   if (any(vapply(ld_swap_chains, Negate(is.null), logical(1)))) {
    out$ld_swap_chains <- ld_swap_chains
   }
  }
 }

 if (!is.null(fit$selection_s_trace)) {
  out$selection_s_trace <- set_trace_matrix(
   fit$selection_s_trace, "selection_s_trace"
  )
 }
 if (!is.null(fit$selection_s_acceptance)) {
  out$selection_s_acceptance <- stats::setNames(
   as.numeric(fit$selection_s_acceptance),
   trait_names
  )
 }

 tol <- 1e-8
 valid_prob <- vapply(comp_prob, function(x) {
  all(is.finite(x)) &&
   all(x >= -tol & x <= 1 + tol) &&
   all(abs(rowSums(x) - 1) <= 1e-6)
 }, logical(1))
 if (!all(valid_prob)) {
  stop("CSR BayesR component probabilities must be finite, within [0, 1], and sum to 1 by marker.")
 }

 .stblr_ensure_ld_swap_fields(out)
}

#' Fit an Exact CSR Summary-Statistics BayesR Model
#'
#' Fits a supported exact CSR BayesR model for GWAS summary statistics and a
#' disk-backed sparse LD matrix. The standard posterior inclusion probability
#' `dm` is `P(component > 0)`, and `comp_prob` stores marker-by-component
#' posterior probabilities with `component_0` as the null component. This
#' explicit convenience wrapper is also available through
#' `stblr_csr(method = "bayesr")`.
#'
#' Exact CSR BayesR supports multiple chains, compact per-chain output with
#' `keep_chains = TRUE`, chain seeds, and `updateE = TRUE` with strict residual
#' SSE validation. Optional BayesR LD-swap/MH moves relocate the full
#' `(component, b)` state from an active marker to a null LD friend.
#' Active/active swaps, annotation-specific swap priors, and scheduled CSR
#' BayesR are not implemented.
#'
#' @param stats Summary-statistics object with `yy`, `ww`, `wy`, and usually
#'   `n`/`m`, as used by [stblr_csr()].
#' @param Glist Optional genotype list containing `Glist$sparseLD$prefix`.
#'   Required when `ld_prefix` is not supplied.
#' @param pi Initial BayesR mixture probabilities. Defaults to total non-null
#'   probability 0.001 distributed evenly across non-null components.
#' @param mixture_var BayesR mixture variance multipliers. The first value must
#'   be zero for the null component and remaining values must be positive.
#' @param alpha Dirichlet prior shapes for mixture probabilities. Defaults to
#'   `pi * 5e5`.
#' @param h2 Initial heritability.
#' @param adjE Residual adjustment factor.
#' @param nit,nburn,nthin MCMC iteration controls.
#' @param ncores Number of OpenMP threads.
#' @param seed Base sampler seed.
#' @param nchains Number of independent chains per trait.
#' @param keep_chains Logical; return compact per-chain BayesR summaries.
#' @param chain_seeds Optional integer seeds, one per chain.
#' @param updateB,updateE,updatePi Logical sampler update controls.
#' @param updateE_start Zero-based first iteration where `updateE` is allowed.
#'   Defaults to 0.
#' @param updateE_every Positive update interval for `updateE`.
#' @param updateLDswap Logical; attempt active/null BayesR LD-swap
#'   Metropolis-Hastings moves. The move relocates `(component, b)` together and
#'   assumes the current plain CSR BayesR model with global `pi` and
#'   `mixture_var`.
#' @param ld_swap_prob Probability per MCMC iteration of attempting LD-swap
#'   moves when `updateLDswap = TRUE`.
#' @param ld_swap_r2 Minimum LD r-squared for candidate swap partners.
#' @param ld_swap_max_friends Maximum number of high-LD friends stored per
#'   marker for swap proposals.
#' @param ld_swap_moves Number of swap attempts when LD-swap is triggered.
#' @param selection_s Optional fixed global BayesS-style MAF-scaling parameter
#'   for ordinary CSR BayesR. For the standardized-genotype CSR effect scale,
#'   non-null component prior variances are scaled by
#'   `h^(selection_s + 1)`, where `h = 2p(1-p)`. `selection_s = NULL`
#'   preserves the ordinary BayesR prior.
#' @param estimate_selection_s Logical; estimate one trait-specific
#'   BayesS-style `selection_s` by Metropolis-Hastings for CSR BayesR. Fixed
#'   `selection_s` and sampled `selection_s` are mutually exclusive.
#' @param selection_s_init Initial value for sampled `selection_s`.
#' @param selection_s_prior Numeric length-2 lower and upper bounds for the
#'   uniform sampled-`selection_s` prior. Only used when
#'   `estimate_selection_s = TRUE`.
#' @param selection_s_proposal_sd Random-walk proposal standard deviation for
#'   sampled `selection_s`. Only used when `estimate_selection_s = TRUE`.
#' @param ld_prefix,n,m Optional low-level CSR LD prefix, sample sizes, and
#'   marker count.
#' @param nub,nue Prior degrees of freedom.
#' @param use_comp_init,comp_init,use_r_init,r_init,rebuild_r_before_updateE
#'   Advanced initialization controls.
#' @param ... Unsupported options. Passing scheduled CSR controls through `...`
#'   errors clearly.
#'
#' @return A formatted ST-BLR BayesR fit.
#' @export
stblr_csr_bayesr <- function(
  stats,
  Glist = NULL,
  ld_prefix = NULL,
  n = NULL,
  m = NULL,
  pi = NULL,
  mixture_var = c(0, 0.01, 0.1, 1),
  alpha = NULL,
  h2 = 0.3,
  adjE = 0.9,
  nit = 1000,
  nburn = 100,
  nthin = 1,
  ncores = 1,
  seed = 1,
  nchains = 1L,
  keep_chains = FALSE,
  chain_seeds = NULL,
  updateB = TRUE,
  updateE = TRUE,
  updatePi = TRUE,
  updateE_start = NULL,
  updateE_every = 1L,
  updateLDswap = FALSE,
  ld_swap_prob = 0.05,
  ld_swap_r2 = 0.8,
  ld_swap_max_friends = 50L,
  ld_swap_moves = 1L,
  selection_s = NULL,
  estimate_selection_s = FALSE,
  selection_s_init = 0,
  selection_s_prior = c(-3, 2),
  selection_s_proposal_sd = 0.35,
  nub = 4,
  nue = 4,
  use_comp_init = FALSE,
  comp_init = NULL,
  use_r_init = FALSE,
  r_init = NULL,
  rebuild_r_before_updateE = FALSE,
  ...
) {
 dots <- list(...)
 if (length(dots)) {
  dot_names <- names(dots)
  if (is.null(dot_names)) dot_names <- rep("", length(dots))
  unsupported <- intersect(
   dot_names,
   "scheduled"
  )
  if ("scheduled" %in% unsupported && isTRUE(dots$scheduled)) {
   stop("scheduled CSR BayesR is not supported by stblr_csr_bayesr().")
  }
  unknown <- setdiff(dot_names[nzchar(dot_names)], unsupported)
  if (length(unknown)) {
   stop("Unsupported argument(s) in ...: ", paste(unknown, collapse = ", "))
  }
 }
 .validate_ld_swap_args(
  updateLDswap, ld_swap_prob, ld_swap_r2, ld_swap_max_friends, ld_swap_moves
 )
 .stblr_validate_sampled_selection_s(
  estimate_selection_s = estimate_selection_s,
  selection_s = selection_s,
  selection_s_init = selection_s_init,
  selection_s_prior = selection_s_prior,
  selection_s_proposal_sd = selection_s_proposal_sd
 )
 if (!is.numeric(nchains) || length(nchains) != 1L ||
     !is.finite(nchains) || nchains < 1 || nchains != floor(nchains)) {
  stop("nchains must be a positive integer scalar.")
 }
 nchains <- as.integer(nchains)
 if (!is.numeric(ncores) || length(ncores) != 1L ||
     !is.finite(ncores) || ncores < 1 || ncores != floor(ncores)) {
  stop("ncores must be a positive integer scalar.")
 }
 ncores <- as.integer(ncores)
 if (!is.logical(keep_chains) || length(keep_chains) != 1L ||
     is.na(keep_chains)) {
  stop("keep_chains must be TRUE or FALSE.")
 }
 if (!is.null(chain_seeds)) {
  if (!is.numeric(chain_seeds) || length(chain_seeds) != nchains ||
      anyNA(chain_seeds) || any(!is.finite(chain_seeds)) ||
      any(chain_seeds != floor(chain_seeds))) {
   stop("chain_seeds must be NULL or an integer/numeric vector of length nchains.")
  }
  chain_seeds <- as.integer(chain_seeds)
 } else {
  chain_seeds <- integer()
 }
 if (is.null(updateE_start)) {
  updateE_start <- 0L
 }
 if (!is.numeric(updateE_start) || length(updateE_start) != 1L ||
     !is.finite(updateE_start) || updateE_start < 0 ||
     updateE_start != floor(updateE_start)) {
  stop("updateE_start must be NULL or a non-negative integer scalar.")
 }
 updateE_start <- as.integer(updateE_start)
 if (!is.numeric(updateE_every) || length(updateE_every) != 1L ||
     !is.finite(updateE_every) || updateE_every < 1 ||
     updateE_every != floor(updateE_every)) {
  stop("updateE_every must be a positive integer scalar.")
 }
 updateE_every <- as.integer(updateE_every)
 if (!is.numeric(mixture_var) || length(mixture_var) < 2L ||
     any(!is.finite(mixture_var)) || mixture_var[1L] != 0 ||
     any(mixture_var[-1L] <= 0)) {
  stop("mixture_var must start with 0 and have positive non-null components.")
 }
 if (is.null(pi)) {
  active <- 0.001
  pi <- c(1 - active, rep(active / (length(mixture_var) - 1L),
                         length(mixture_var) - 1L))
 }
 if (!is.numeric(pi) || length(pi) != length(mixture_var) ||
     any(!is.finite(pi)) || any(pi < 0) || sum(pi) <= 0) {
  stop("pi must be a non-negative finite vector matching mixture_var with positive sum.")
 }
 pi <- pi / sum(pi)
 if (is.null(alpha)) alpha <- pmax(pi * 5e5, .Machine$double.eps)
 if (!is.numeric(alpha) || length(alpha) != length(mixture_var) ||
     any(!is.finite(alpha)) || any(alpha <= 0)) {
  stop("alpha must be a positive finite vector matching mixture_var.")
 }

 nt <- length(stats$yy)
 if (is.null(n)) n <- stats$n
 if (is.null(n)) stop("n must be supplied or available as stats$n.")
 if (length(n) == 1L) n <- rep(n, nt)
 if (is.null(m)) m <- if (!is.null(stats$m)) stats$m else length(stats$ww[[1L]])
 if (is.null(ld_prefix)) {
  if (is.null(Glist) || is.null(Glist$sparseLD$prefix)) {
   stop("Provide ld_prefix or run make_sparse_ld() and supply Glist.")
  }
  ld_prefix <- Glist$sparseLD$prefix
 }

 trait_names <- names(stats$yy)
 if (is.null(trait_names)) trait_names <- paste0("T", seq_len(nt))
 variable_names <- names(stats$ww[[1L]])
 if (is.null(variable_names)) variable_names <- paste0("V", seq_len(m))
 selection_s_info <- .stblr_prepare_csr_selection_s(
  selection_s = selection_s,
  Glist = Glist,
  m = m,
  scheduled = FALSE,
  backend = "bayesr",
  return_log_h = estimate_selection_s
 )

 vy <- as.numeric(stats$yy) / (as.numeric(n) - 1)
 pi_active <- sum(pi[-1L])
 if (!is.finite(pi_active) || pi_active <= 0) {
  stop("pi must allocate positive probability to non-null components.")
 }
 pri <- list(
  vy = vy,
  B = diag((vy * h2) / (m * pi_active), nt, nt),
  E = diag(vy * (1 - h2), nt, nt),
  ssb_prior = diag(((nub - 2) / nub) * (vy * h2) / (m * pi_active), nt, nt),
  sse_prior = diag(((nue - 2) / nue) * (vy * (1 - h2)), nt, nt)
 )
 for (x in c("B", "E", "ssb_prior", "sse_prior")) {
  rownames(pri[[x]]) <- colnames(pri[[x]]) <- trait_names
 }
 pri$ssb_prior_list <- split(pri$ssb_prior, rep(seq_len(nt), each = nt))
 pri$sse_prior_list <- split(pri$sse_prior, rep(seq_len(nt), each = nt))

 if (is.null(comp_init)) comp_init <- lapply(seq_len(nt), function(i) rep(0, m))
 if (is.null(r_init)) r_init <- stats$wy

 raw <- stblr_cpg_omp_csr_bayesr(
  wy = stats$wy,
  ww = stats$ww,
  yy = stats$yy,
  b_init = lapply(seq_len(nt), function(i) rep(0, m)),
  comp_init = comp_init,
  use_comp_init = use_comp_init,
  r_init = r_init,
  use_r_init = use_r_init,
  rebuild_r_before_updateE = rebuild_r_before_updateE,
  ld_prefix = ld_prefix,
  B = pri$B,
  E = pri$E,
  ssb_prior = pri$ssb_prior_list,
  sse_prior = pri$sse_prior_list,
  pi = as.numeric(pi),
  mixture_var = as.numeric(mixture_var),
  alpha = as.numeric(alpha),
  nub = nub,
  nue = nue,
  updateB = updateB,
  updateE = updateE,
  updatePi = updatePi,
  adjE = adjE,
  n = as.integer(n),
  nit = as.integer(nit),
  nburn = as.integer(nburn),
  nthin = as.integer(nthin),
  ncores = ncores,
  seed = as.integer(seed),
  nchains = nchains,
  keep_chains = keep_chains,
  chain_seeds = chain_seeds,
  updateE_start = updateE_start,
  updateE_every = updateE_every,
  updateLDswap = updateLDswap,
  ld_swap_prob = ld_swap_prob,
  ld_swap_r2 = ld_swap_r2,
  ld_swap_max_friends = as.integer(ld_swap_max_friends),
  ld_swap_moves = as.integer(ld_swap_moves),
  selection_s_prior_scale = selection_s_info$prior_scale,
  estimate_selection_s = estimate_selection_s,
  selection_s_init = selection_s_init,
  selection_s_prior = selection_s_prior,
  selection_s_proposal_sd = selection_s_proposal_sd,
  selection_s_log_h = selection_s_info$log_h
 )

 if (.is_stblr_raw_v1(raw)) {
  if (isTRUE(selection_s_info$fixed)) {
   raw$selection$mean <- stats::setNames(rep(selection_s_info$selection_s, nt), trait_names)
  }
  fit <- .format_stblr_raw_v1(raw, trait_names, variable_names)
 } else {
  fit <- .format_stblr_csr_bayesr_fit(
   raw,
   nt = nt,
   m = m,
   trait_names = trait_names,
   variable_names = variable_names,
   n_components = length(mixture_var)
  )
 }
 if (isTRUE(estimate_selection_s)) {
  keep_idx <- seq.int(nburn + 1L, nrow(fit$selection_s_trace))
  s_trace_keep <- fit$selection_s_trace[keep_idx, , drop = FALSE]
  fit$selection_s <- colMeans(s_trace_keep)
  fit$selection_s_sd <- apply(s_trace_keep, 2L, stats::sd)
  fit$selection_s_min <- apply(s_trace_keep, 2L, min)
  fit$selection_s_max <- apply(s_trace_keep, 2L, max)
 }
 fit$input <- list(
  method = "bayesr",
  model = "bayesr",
  backend = "csr_bayesr",
  data_level = "summary",
  annotations = FALSE,
  scheduled = FALSE,
  nchains = nchains,
  keep_chains = keep_chains,
  updateLDswap = updateLDswap,
  ld_swap_prob = ld_swap_prob,
  ld_swap_r2 = ld_swap_r2,
  ld_swap_max_friends = as.integer(ld_swap_max_friends),
  ld_swap_moves = as.integer(ld_swap_moves),
  n = as.integer(n),
  m = m,
  nt = nt,
  h2 = h2,
  nub = nub,
  nue = nue,
  adjE = adjE,
  nit = as.integer(nit),
  nburn = as.integer(nburn),
  nthin = as.integer(nthin),
  mixture_var = as.numeric(mixture_var),
  pi = as.numeric(pi),
  alpha = as.numeric(alpha),
  updateE = updateE,
  updateE_start = updateE_start,
  updateE_every = updateE_every,
  selection_s = selection_s_info$selection_s,
  selection_s_fixed = selection_s_info$fixed,
  selection_s_exponent = selection_s_info$exponent,
  estimate_selection_s = estimate_selection_s,
  selection_s_init = if (isTRUE(estimate_selection_s)) selection_s_init else NULL,
  selection_s_prior = if (isTRUE(estimate_selection_s)) selection_s_prior else NULL,
  selection_s_proposal_sd = if (isTRUE(estimate_selection_s)) selection_s_proposal_sd else NULL,
  selection_s_scale = "standardized_genotype_effect",
  chain_seeds = if (length(chain_seeds)) chain_seeds else NULL
 )
 fit
}

# Internal clearer helper name for the exact CSR BayesR backend.
.fit_stblr_csr_bayesr <- stblr_csr_bayesr

# Compatibility alias for scripts/tests that used the old internal helper name.
.stblr_csr_bayesr_experimental <- function(
  Glist = NULL,
  stats,
  ld_prefix = NULL,
  n = NULL,
  m = NULL,
  h2 = 0.3,
  mixture_var = c(0, 0.01, 0.1, 1),
  pi = NULL,
  alpha = NULL,
  nub = 4,
  nue = 4,
  updateB = TRUE,
  updateE = TRUE,
  updatePi = TRUE,
  adjE = 0.9,
  nit = 1000,
  nburn = 100,
  nthin = 1,
  ncores = 1,
  seed = 10,
  nchains = 1L,
  keep_chains = FALSE,
  chain_seeds = NULL,
  scheduled = FALSE,
  updateLDswap = FALSE,
  ld_swap_prob = 0.05,
  ld_swap_r2 = 0.8,
  ld_swap_max_friends = 50L,
  ld_swap_moves = 1L,
  use_comp_init = FALSE,
  comp_init = NULL,
  use_r_init = FALSE,
  r_init = NULL,
  rebuild_r_before_updateE = FALSE,
  updateE_start = NULL,
  updateE_every = 1L
) {
 .fit_stblr_csr_bayesr(
  stats = stats,
  Glist = Glist,
  ld_prefix = ld_prefix,
  n = n,
  m = m,
  pi = pi,
  mixture_var = mixture_var,
  alpha = alpha,
  h2 = h2,
  adjE = adjE,
  nit = nit,
  nburn = nburn,
  nthin = nthin,
  ncores = ncores,
  seed = seed,
  nchains = nchains,
  keep_chains = keep_chains,
  chain_seeds = chain_seeds,
  updateB = updateB,
  updateE = updateE,
  updatePi = updatePi,
  updateE_start = updateE_start,
  updateE_every = updateE_every,
  nub = nub,
  nue = nue,
  use_comp_init = use_comp_init,
  comp_init = comp_init,
  use_r_init = use_r_init,
  r_init = r_init,
  rebuild_r_before_updateE = rebuild_r_before_updateE,
  scheduled = scheduled,
  updateLDswap = updateLDswap,
  ld_swap_prob = ld_swap_prob,
  ld_swap_r2 = ld_swap_r2,
  ld_swap_max_friends = ld_swap_max_friends,
  ld_swap_moves = ld_swap_moves
 )
}

.resolve_Glist_markers <- function(Glist, chr = NULL, cls = NULL) {
  bedfiles <- as.character(Glist$bedfiles)
  has_bedfile <- !is.na(bedfiles) & nzchar(bedfiles)

  if (is.null(chr)) {
    chr <- which(has_bedfile)
  } else {
    chr <- as.integer(chr)
  }

  if (length(chr) < 1L || anyNA(chr)) {
    stop("chr must contain valid chromosome/file indices.")
  }

  if (any(chr < 1L | chr > length(bedfiles))) {
    stop("chr contains indices outside Glist$bedfiles.")
  }

  missing_bed <- is.na(bedfiles[chr]) | !nzchar(bedfiles[chr])
  if (any(missing_bed)) {
    stop(
      "Glist$bedfiles is missing for chromosome/file index: ",
      paste(chr[missing_bed], collapse = ", ")
    )
  }

  if (is.null(cls)) {
    cls <- lapply(chr, function(cc) {
      match(Glist$rsidsLD[[cc]], Glist$rsids[[cc]])
    })
  } else if (!is.list(cls)) {
    cls <- list(cls)
  }

  if (length(cls) != length(chr)) {
    stop("cls must have one element per chromosome/file in chr.")
  }

  cls <- lapply(cls, as.integer)
  names(cls) <- paste0("chr", chr)

  af <- Map(function(cc, cl) Glist$af[[cc]][cl], chr, cls)

  marker_names <- unlist(
    Map(function(cc, cl) Glist$rsids[[cc]][cl], chr, cls),
    use.names = FALSE
  )

  list(
    chr = chr,
    bed_files = bedfiles[chr],
    cls = cls,
    af = af,
    marker_names = marker_names
  )
}

.validate_ld_swap_args <- function(updateLDswap, ld_swap_prob, ld_swap_r2,
                                   ld_swap_max_friends, ld_swap_moves) {
  if (!is.logical(updateLDswap) || length(updateLDswap) != 1L ||
      is.na(updateLDswap)) {
    stop("updateLDswap must be TRUE or FALSE.")
  }
  if (!is.numeric(ld_swap_prob) || length(ld_swap_prob) != 1L ||
      !is.finite(ld_swap_prob) || ld_swap_prob < 0 || ld_swap_prob > 1) {
    stop("ld_swap_prob must be a finite numeric scalar in [0, 1].")
  }
  if (!is.numeric(ld_swap_r2) || length(ld_swap_r2) != 1L ||
      !is.finite(ld_swap_r2) || ld_swap_r2 < 0 || ld_swap_r2 > 1) {
    stop("ld_swap_r2 must be a finite numeric scalar in [0, 1].")
  }
  if (!is.numeric(ld_swap_max_friends) || length(ld_swap_max_friends) != 1L ||
      !is.finite(ld_swap_max_friends) || ld_swap_max_friends < 1 ||
      ld_swap_max_friends != floor(ld_swap_max_friends)) {
    stop("ld_swap_max_friends must be a positive integer.")
  }
  if (!is.numeric(ld_swap_moves) || length(ld_swap_moves) != 1L ||
      !is.finite(ld_swap_moves) || ld_swap_moves < 0 ||
      ld_swap_moves != floor(ld_swap_moves)) {
    stop("ld_swap_moves must be a non-negative integer.")
  }
  invisible(TRUE)
}

#' Fit ST-BLR from BED Sufficient Statistics and Sparse LD
#'
#' Annotation-unaware summary-statistics CSR interface. Fits ST-BLR samplers
#' using sufficient statistics from [bed_xtx_xty()] and a disk-backed CSR LD
#' prefix from [sparseLD_stream_CSR()].
#'
#' Available methods are `method = "bayesC"` and `method = "bayesR"`. Both
#' support fixed global `selection_s`, sampled trait-specific `selection_s`,
#' LD-swap/fine-mapping diagnostics where `updateLDswap = TRUE`, and native
#' multi-chain posterior summaries.
#'
#'
#' @param Glist Optional qgg genotype list containing `sparseLD`. If
#'   `ld_prefix` is `NULL`, `Glist$sparseLD$prefix` is used.
#' @param stats Sufficient statistics returned by [bed_xtx_xty()].
#' @param ld_prefix Prefix of the disk-backed CSR LD files. If `NULL`, this is
#'   taken from `Glist$sparseLD$prefix`.
#' @param n Sample size. Defaults to `stats$n` when available.
#' @param m Number of markers. Inferred from `stats` when omitted.
#' @param method CSR model to fit. `"bayesc"` keeps the historical BayesC CSR
#'   behavior and is the default. `"bayesr"` dispatches to the exact CSR BayesR
#'   backend exposed by [stblr_csr_bayesr()]. Scheduled CSR BayesR is not
#'   currently implemented.
#' @param pi_init Initial marker inclusion probability. Defaults to 0.001.
#' @param pi_vb_init Inclusion probability used when initializing marker-effect
#'   variance. Defaults to `pi_init`.
#' @param pi_prior_mean Prior mean for the marker inclusion probability.
#'   Defaults to 0.001.
#' @param pi_prior_strength Total Beta prior strength for the marker inclusion
#'   probability. Defaults to 5e5.
#' @param pi_prior_a,pi_prior_b Optional explicit Beta prior shape parameters.
#'   When supplied, these override `pi_prior_mean` and `pi_prior_strength`.
#' @param h2 Initial heritability. Defaults to 0.3.
#' @param selection_s Optional fixed global BayesS-style MAF-scaling parameter
#'   for ordinary annotation-unaware CSR BayesC and BayesR. The default
#'   `selection_s = NULL` with `estimate_selection_s = FALSE` fits the ordinary
#'   model. A finite numeric scalar with `estimate_selection_s = FALSE` fits a
#'   fixed global-S model. `selection_s` must remain `NULL` when
#'   `estimate_selection_s = TRUE`; fixed and sampled S cannot both be
#'   requested.
#' @param estimate_selection_s Logical; for unscheduled annotation-unaware CSR
#'   BayesC and BayesR, estimate one trait-specific BayesS-style `selection_s`
#'   by Metropolis-Hastings. `selection_s = NULL` and
#'   `estimate_selection_s = TRUE` samples S separately for each trait and
#'   chain. Sampled `selection_s` is not supported by scheduled CSR BayesC,
#'   BayesC-like annotation-aware CSR backends, or BED backends.
#' @param selection_s_init Initial value for sampled `selection_s`. Defaults to
#'   0 and is used only when `estimate_selection_s = TRUE`.
#' @param selection_s_prior Numeric length-2 lower and upper bounds for the
#'   uniform sampled-`selection_s` prior. Only used when
#'   `estimate_selection_s = TRUE`; it does not affect ordinary BayesC or fixed
#'   `selection_s`. Defaults to `c(-3, 2)`, i.e. Uniform(-3, 2).
#' @param selection_s_proposal_sd Random-walk proposal standard deviation for
#'   sampled `selection_s`. Only used when `estimate_selection_s = TRUE`; it
#'   does not affect ordinary BayesC or fixed `selection_s`. Defaults to 0.35.
#' @param nub,nue Prior degrees of freedom.
#' @param updateB,updateE,updatePi Logical sampler update controls.
#' @param adjE Residual adjustment factor. For sparse-LD summary-statistic
#'   models, values greater than 0 can stabilize residual-variance updates.
#' @param nit,nburn,nthin MCMC iteration controls.
#' @param ncores Number of OpenMP threads. For multi-chain regular and
#'   scheduled CSR fits, threads are assigned over trait-chain tasks.
#' @param seed Sampler seed.
#' @param nchains Number of independent MCMC chains per trait. Chains are run
#'   as independent genome-wide trait-by-chain tasks in C++ and aggregated
#'   before return for regular and scheduled CSR fits.
#' @param keep_chains Logical; when `TRUE` and `scheduled = FALSE`, return
#'   compact per-chain `dm`, `bm`, and LD-swap diagnostics in `chains`.
#'   Compact per-chain output is not yet supported for scheduled CSR.
#' @param chain_seeds Optional numeric or integer vector of length `nchains`.
#'   When supplied, each value is used as the base seed for that chain with a
#'   deterministic trait offset. Matrix seed inputs are not currently
#'   supported.
#' @param scheduled Use the scheduled sparse-LD sampler.
#' @param full_sweep_every,null_skip_base,null_skip_max,candidate_threshold
#'   Scheduled sampler controls.
#' @param candidate_lifetime Number of scheduled sweeps that candidates remain
#'   active.
#' @param skip_nulls_burnin_only Restrict null-marker skipping to burn-in.
#' @param wakeup_ld_neighbors,wakeup_diff_threshold,wakeup_max_neighbors
#'   Scheduled neighbor wake-up controls.
#' @param use_d_init,use_r_init,rebuild_r_before_updateE Initialization and
#'   residual rebuilding controls.
#' @param updateLDswap Logical; attempt optional LD-swap Metropolis-Hastings
#'   moves among included markers and high-LD excluded neighbors. These moves
#'   can improve mixing in fine-mapping and high-LD regions, are off by
#'   default, and use an exact proposal-ratio correction in the base CSR
#'   sampler.
#' @param ld_swap_prob Probability per MCMC iteration of attempting LD-swap
#'   moves when `updateLDswap = TRUE`.
#' @param ld_swap_r2 Minimum LD r-squared for candidate swap partners.
#' @param ld_swap_max_friends Maximum number of high-LD friends stored per
#'   marker, prioritized by highest r-squared.
#' @param ld_swap_moves Number of swap attempts when LD-swap is triggered.
#' @param mixture_var,pi,alpha BayesR mixture variance multipliers, initial
#'   mixture probabilities, and Dirichlet prior shapes. Used only when
#'   `method = "bayesr"` and passed through to [stblr_csr_bayesr()]. BayesR
#'   defaults are supplied by [stblr_csr_bayesr()] when these are left `NULL`.
#' @param updateE_start,updateE_every BayesR residual-variance update schedule.
#'   Used only when `method = "bayesr"`.
#' @param use_comp_init,comp_init BayesR component-state initialization
#'   controls. Used only when `method = "bayesr"`.
#' @return A formatted ST-BLR fit. Common posterior fields include `dm`
#'   (marker-by-trait posterior inclusion or non-null probability), `bm`
#'   (marker-by-trait posterior mean effects), `vbs`, `vgs`, `ves`, `vle`,
#'   `vld`, and `input` model/backend/settings metadata. Chain-capable CSR fits
#'   provide `dm_sd`, `dm_min`, `dm_max`, `bm_sd`, `bm_min`, and `bm_max`;
#'   standard traces are averaged by iteration across chains. With
#'   `keep_chains = TRUE`, compact per-chain summaries are returned in
#'   `chains` for unscheduled CSR fits.
#'
#'   For scheduled CSR fits, `pis` contains the sampled inclusion-probability
#'   trace for each trait, averaged by iteration across chains when
#'   `nchains > 1`. LD-swap remains available only for `scheduled = FALSE`;
#'   enabled fits include `ld_swap` and, where chain summaries are retained,
#'   `ld_swap_chains` or chain-level LD-swap entries.
#'
#'   For `method = "bayesR"`, `dm` is `P(component > 0)`, `comp_prob` contains
#'   marker-by-component posterior probabilities by trait, and
#'   `dm_component_mean` contains posterior mean component/effect-class
#'   summaries. The BayesR null component column is `component_0`, so
#'   `dm = 1 - P(component_0)`.
#'
#'   When sampled `selection_s` is enabled, the fit also contains
#'   `selection_s`, `selection_s_sd`, `selection_s_min`, `selection_s_max`,
#'   `selection_s_trace`, and `selection_s_acceptance`. `selection_s_trace` is
#'   an iteration x trait matrix, `selection_s` is the posterior mean by trait,
#'   and `selection_s_acceptance` is the MH acceptance rate by trait. With
#'   `keep_chains = TRUE`, chain-level sampled-S output is available as
#'   `fit$chains[[trait]][[chain]]$selection_s` and
#'   `fit$chains[[trait]][[chain]]$selection_s_acceptance`.
#'
#'   Fine-mapping diagnostics are available through PIP summaries in `dm` and
#'   LD-swap output when `updateLDswap = TRUE`. Credible-set construction is
#'   performed by helper functions such as [make_credible_sets()] and
#'   [extract_stblr_finemap_loci()] from posterior inclusion probabilities and
#'   LD, rather than being a separate sampler return object.
#'
#' @details
#' CSR effects are on the standardized-genotype scale. The BayesS-style
#' MAF-dependent prior scale used by fixed and sampled `selection_s` is
#' `q_j(S) = h_j^(S + 1)`, where `h_j = 2 p_j (1 - p_j)` and `p_j` is the
#' minor allele frequency. The `+1` exponent appears because the sampler effects
#' are standardized-genotype-scale effects rather than allele-scale effects.
#'
#' For fixed `selection_s`, CSR BayesC uses
#' `b_j | d_j = 1, vb, S ~ N(0, vb * q_j(S))`. CSR BayesR uses
#' `b_j | component_j = m, vb, S ~ N(0, vb * gamma_m * q_j(S))`, with
#' component 0 as the point-mass null. In BayesR output, the null component
#' column is `component_0` and `dm = 1 - P(component_0)`.
#'
#' For sampled `selection_s`, BayesC uses the trait-specific MH log posterior
#' contribution
#'
#' ```text
#' log p(S_t | b_t, d_t, vb_t)
#'   = log p(S_t)
#'     - 0.5 sum_\{j: d_jt = 1\} [
#'         log q_j(S_t) + b_jt^2 / (vb_t q_j(S_t))
#'       ]
#' ```
#'
#' BayesR uses
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
#' fit_bc <- stblr_csr(stats = stats, Glist = Glist, method = "bayesc")
#' fit_br <- stblr_csr(stats = stats, Glist = Glist, method = "bayesr")
#'
#' fit_fixed_s <- stblr_csr(
#'   stats = stats,
#'   Glist = Glist,
#'   method = "bayesC",
#'   selection_s = -0.5
#' )
#'
#' fit_sampled_s_bc <- stblr_csr(
#'   stats = stats,
#'   Glist = Glist,
#'   method = "bayesC",
#'   estimate_selection_s = TRUE
#' )
#'
#' fit_sampled_s_br <- stblr_csr(
#'   stats = stats,
#'   Glist = Glist,
#'   method = "bayesR",
#'   estimate_selection_s = TRUE,
#'   selection_s_prior = c(-3, 2),
#'   selection_s_proposal_sd = 0.35
#' )
#' }
#' @export
stblr_csr <- function(Glist=NULL, stats, ld_prefix=NULL, n = NULL, m = NULL,
                      method = c("bayesc", "bayesr"),
                      pi_init = 0.001, pi_vb_init = NULL,
                      pi_prior_mean = 0.001, pi_prior_strength = 5e5,
                      pi_prior_a = NULL, pi_prior_b = NULL, h2 = 0.3,
                      selection_s = NULL, estimate_selection_s = FALSE,
                      selection_s_init = 0,
                      selection_s_prior = c(-3, 2),
                      selection_s_proposal_sd = 0.35,
                      nub = 4, nue = 4, updateB = TRUE, updateE = TRUE,
                      updatePi = TRUE, adjE = 0.9, nit = 1000, nburn = 100,
                      nthin = 1, ncores = 1, seed = 10, nchains = 1L,
                      keep_chains = FALSE, chain_seeds = NULL,
                      scheduled = FALSE,
                      full_sweep_every = 1, null_skip_base = 50,
                      null_skip_max = 200, candidate_threshold = 1e-3,
                      candidate_lifetime = 20, skip_nulls_burnin_only = FALSE,
                      wakeup_ld_neighbors = TRUE, wakeup_diff_threshold = 0,
                      wakeup_max_neighbors = 0, use_d_init = FALSE,
                      use_r_init = FALSE, rebuild_r_before_updateE = FALSE,
                      updateLDswap = FALSE, ld_swap_prob = 0.05,
                      ld_swap_r2 = 0.8, ld_swap_max_friends = 50,
                      ld_swap_moves = 1,
                      mixture_var = NULL, pi = NULL, alpha = NULL,
                      updateE_start = NULL, updateE_every = NULL,
                      use_comp_init = FALSE, comp_init = NULL) {
 method <- tolower(method)
 if (length(method) > 1L) method <- method[1L]
 if (length(method) < 1L || is.na(method) ||
     !method %in% c("bayesc", "bayesr")) {
  stop("method must be one of 'bayesc' or 'bayesr'.")
 }
 .stblr_validate_sampled_selection_s(
  estimate_selection_s = estimate_selection_s,
  selection_s = selection_s,
  selection_s_init = selection_s_init,
  selection_s_prior = selection_s_prior,
  selection_s_proposal_sd = selection_s_proposal_sd
 )
 if (method == "bayesr") {
  bayesc_only <- character()
  if (!missing(pi_init) && !identical(pi_init, 0.001)) {
   bayesc_only <- c(bayesc_only, "pi_init")
  }
  if (!missing(pi_vb_init) && !is.null(pi_vb_init)) {
   bayesc_only <- c(bayesc_only, "pi_vb_init")
  }
  if (!missing(pi_prior_mean) && !identical(pi_prior_mean, 0.001)) {
   bayesc_only <- c(bayesc_only, "pi_prior_mean")
  }
  if (!missing(pi_prior_strength) && !identical(pi_prior_strength, 5e5)) {
   bayesc_only <- c(bayesc_only, "pi_prior_strength")
  }
  if (!missing(pi_prior_a) && !is.null(pi_prior_a)) {
   bayesc_only <- c(bayesc_only, "pi_prior_a")
  }
  if (!missing(pi_prior_b) && !is.null(pi_prior_b)) {
   bayesc_only <- c(bayesc_only, "pi_prior_b")
  }
  if (length(bayesc_only)) {
   stop(
    "`", paste(bayesc_only, collapse = "`, `"),
    "` are BayesC-specific. Use `pi` and `alpha` for method = 'bayesr'."
   )
  }
  if (isTRUE(scheduled)) {
   stop(
    "scheduled CSR BayesR is not currently implemented. ",
    "Use method = 'bayesc' for scheduled CSR or scheduled = FALSE for BayesR."
   )
  }
  br_args <- list(
   stats = stats,
   Glist = Glist,
   ld_prefix = ld_prefix,
   n = n,
   m = m,
   h2 = h2,
   adjE = adjE,
   nit = nit,
   nburn = nburn,
   nthin = nthin,
   ncores = ncores,
   seed = seed,
   nchains = nchains,
   keep_chains = keep_chains,
   chain_seeds = chain_seeds,
   updateB = updateB,
   updateE = updateE,
   updatePi = updatePi,
   updateLDswap = updateLDswap,
   ld_swap_prob = ld_swap_prob,
   ld_swap_r2 = ld_swap_r2,
   ld_swap_max_friends = ld_swap_max_friends,
   ld_swap_moves = ld_swap_moves,
   nub = nub,
   nue = nue,
   use_comp_init = use_comp_init,
   comp_init = comp_init,
   use_r_init = use_r_init,
   rebuild_r_before_updateE = rebuild_r_before_updateE,
   selection_s = selection_s,
   estimate_selection_s = estimate_selection_s,
   selection_s_init = selection_s_init,
   selection_s_prior = selection_s_prior,
   selection_s_proposal_sd = selection_s_proposal_sd
  )
  if (!is.null(mixture_var)) br_args$mixture_var <- mixture_var
  if (!is.null(pi)) br_args$pi <- pi
  if (!is.null(alpha)) br_args$alpha <- alpha
  if (!is.null(updateE_start)) br_args$updateE_start <- updateE_start
  if (!is.null(updateE_every)) br_args$updateE_every <- updateE_every
  fit <- do.call(stblr_csr_bayesr, br_args)
  fit$input$method <- "bayesr"
  fit$input$model <- "bayesr"
  fit$input$backend <- "csr_bayesr"
  fit$input$data_level <- "summary"
  return(fit)
 }
 .validate_ld_swap_args(
  updateLDswap, ld_swap_prob, ld_swap_r2, ld_swap_max_friends, ld_swap_moves
 )
 if (!is.numeric(nchains) || length(nchains) != 1L ||
     !is.finite(nchains) || nchains < 1 || nchains != floor(nchains)) {
  stop("nchains must be a positive integer scalar.")
 }
 nchains <- as.integer(nchains)
 if (!is.logical(keep_chains) || length(keep_chains) != 1L ||
     is.na(keep_chains)) {
  stop("keep_chains must be TRUE or FALSE.")
 }
 if (!is.null(chain_seeds)) {
  if (!is.numeric(chain_seeds) || length(chain_seeds) != nchains ||
      anyNA(chain_seeds) || any(!is.finite(chain_seeds)) ||
      any(chain_seeds != floor(chain_seeds))) {
   stop("chain_seeds must be NULL or an integer/numeric vector of length nchains.")
  }
  chain_seeds <- as.integer(chain_seeds)
 } else {
  chain_seeds <- integer()
 }
 if (isTRUE(scheduled) && isTRUE(keep_chains)) {
  stop("scheduled = TRUE, keep_chains = TRUE is not yet supported.")
 }
 if (isTRUE(scheduled) && isTRUE(updateLDswap)) {
  stop("updateLDswap is currently implemented only for scheduled = FALSE.")
 }
 if (isTRUE(scheduled) && isTRUE(estimate_selection_s)) {
  stop("estimate_selection_s is currently supported only for unscheduled CSR BayesC.")
 }
 nt <- length(stats$yy)
 if (is.null(n)) n <- stats$n
 if (is.null(n)) stop("n must be supplied or available as stats$n.")
 if (is.null(m)) m <- if (!is.null(stats$m)) stats$m else length(stats$ww[[1]])
 selection_s_info <- .stblr_prepare_csr_bayesc_selection_s(
  selection_s = selection_s,
  Glist = Glist,
  m = m,
  scheduled = scheduled,
  return_log_h = estimate_selection_s
 )

 if (is.null(ld_prefix)) {
   if (is.null(Glist) || is.null(Glist$sparseLD$prefix)) {
     stop("Provide ld_prefix or run make_sparse_ld() and supply Glist.")
   }
   ld_prefix <- Glist$sparseLD$prefix
 }

 arch <- .resolve_pi_prior(
   pi_init = pi_init,
   pi_vb_init = pi_vb_init,
   pi_prior_mean = pi_prior_mean,
   pi_prior_strength = pi_prior_strength,
   pi_prior_a = pi_prior_a,
   pi_prior_b = pi_prior_b
 )

 trait_names <- names(stats$yy)
 if (is.null(trait_names)) trait_names <- paste0("T", seq_len(nt))
 variable_names <- names(stats$ww[[1]])
 if (is.null(variable_names)) variable_names <- paste0("V", seq_len(m))
 vy <- as.numeric(stats$yy) / (n - 1)
 pri <- list(
  vy = vy,
  B = diag((vy * h2) / (m * arch$pi_vb_init), nt, nt),
  E = diag(vy * (1 - h2), nt, nt),
  ssb_prior = diag(
   ((nub - 2) / nub) * (vy * h2) / (m * arch$pi_prior_mean), nt, nt
  ),
  sse_prior = diag(((nue - 2) / nue) * (vy * (1 - h2)), nt, nt)
 )
 for (x in c("B", "E", "ssb_prior", "sse_prior")) {
  rownames(pri[[x]]) <- colnames(pri[[x]]) <- trait_names
 }
 pri$ssb_prior_list <- split(pri$ssb_prior, rep(seq_len(nt), each = nt))
 pri$sse_prior_list <- split(pri$sse_prior, rep(seq_len(nt), each = nt))
 common <- list(
  wy = stats$wy, ww = stats$ww, yy = stats$yy,
  b_init = lapply(seq_len(nt), function(i) rep(0, m)),
  d_init = lapply(seq_len(nt), function(i) rep(0, m)),
  use_d_init = use_d_init, r_init = stats$wy, use_r_init = use_r_init,
  rebuild_r_before_updateE = rebuild_r_before_updateE, ld_prefix = ld_prefix,
  B = pri$B, E = pri$E, ssb_prior = pri$ssb_prior_list,
  sse_prior = pri$sse_prior_list, pi = arch$pi, nub = nub, nue = nue,
  updateB = updateB, updateE = updateE, updatePi = updatePi, adjE = adjE,
  n = rep(as.integer(n), nt), nit = nit, nburn = nburn, nthin = nthin
 )
 if (scheduled) {
  raw <- do.call(stblr_cpg_omp_csr_scheduled, c(common, list(
   full_sweep_every = full_sweep_every, null_skip_base = null_skip_base,
   null_skip_max = null_skip_max, candidate_threshold = candidate_threshold,
   candidate_lifetime = candidate_lifetime,
   skip_nulls_burnin_only = skip_nulls_burnin_only,
   wakeup_ld_neighbors = wakeup_ld_neighbors,
   wakeup_diff_threshold = wakeup_diff_threshold,
   wakeup_max_neighbors = wakeup_max_neighbors, pi_prior_a = arch$pi_prior_a,
   pi_prior_b = arch$pi_prior_b, ncores = ncores, seed = seed,
   nchains = nchains, keep_chains = keep_chains, chain_seeds = chain_seeds
  )))
 } else {
 raw <- do.call(stblr_cpg_omp_csr, c(common, list(
   pi_prior_a = arch$pi_prior_a, pi_prior_b = arch$pi_prior_b,
   ncores = ncores, seed = seed, nchains = nchains,
   keep_chains = keep_chains, chain_seeds = chain_seeds,
   updateLDswap = updateLDswap,
   ld_swap_prob = ld_swap_prob, ld_swap_r2 = ld_swap_r2,
   ld_swap_max_friends = as.integer(ld_swap_max_friends),
   ld_swap_moves = as.integer(ld_swap_moves),
   selection_s_prior_scale = selection_s_info$prior_scale,
   estimate_selection_s = estimate_selection_s,
   selection_s_init = selection_s_init,
   selection_s_prior = selection_s_prior,
   selection_s_proposal_sd = selection_s_proposal_sd,
   selection_s_log_h = selection_s_info$log_h
  )))
 }
 if (.is_stblr_raw_v1(raw)) {
  if (isTRUE(selection_s_info$fixed)) {
   raw$selection$mean <- stats::setNames(rep(selection_s_info$selection_s, nt), trait_names)
  }
  fit <- .format_stblr_raw_v1(raw, trait_names, variable_names)
 } else {
  fit <- .format_stblr_fit(raw, nt, m, trait_names, variable_names)
 }
 if (isTRUE(estimate_selection_s)) {
  keep_idx <- seq.int(nburn + 1L, nrow(fit$selection_s_trace))
  s_trace_keep <- fit$selection_s_trace[keep_idx, , drop = FALSE]
  fit$selection_s <- colMeans(s_trace_keep)
  fit$selection_s_sd <- apply(s_trace_keep, 2L, stats::sd)
  fit$selection_s_min <- apply(s_trace_keep, 2L, min)
  fit$selection_s_max <- apply(s_trace_keep, 2L, max)
 }
 fit$input <- c(list(
  n = n, m = m, nt = nt, h2 = h2, nub = nub, nue = nue, vy = vy,
  B = pri$B, E = pri$E, ssb_prior = pri$ssb_prior,
  sse_prior = pri$sse_prior, updateB = updateB, updateE = updateE,
  updatePi = updatePi, adjE = adjE, nit = nit, nburn = nburn, nthin = nthin,
  ncores = ncores, seed = seed, nchains = nchains,
  keep_chains = keep_chains,
  chain_seeds = if (length(chain_seeds)) chain_seeds else NULL,
  chain_seed_rule = if (length(chain_seeds)) {
   "chain_seeds[chain] + 1000003 * (trait + 1)"
  } else if (nchains == 1L) {
   "seed + 1000003 * (trait + 1)"
  } else {
   "seed + 1000003 * (trait + 1) + 9176 * (chain + 1)"
  },
  scheduled = scheduled, ld_prefix = ld_prefix,
  method = "bayesc", model = "bayesc",
  backend = if (isTRUE(scheduled)) "csr_scheduled_bayesc" else "csr_bayesc",
  data_level = "summary",
  annotations = FALSE,
  selection_s = selection_s_info$selection_s,
  selection_s_fixed = selection_s_info$fixed,
  selection_s_exponent = selection_s_info$exponent,
  estimate_selection_s = estimate_selection_s,
  selection_s_init = if (isTRUE(estimate_selection_s)) selection_s_init else NULL,
  selection_s_prior = if (isTRUE(estimate_selection_s)) selection_s_prior else NULL,
  selection_s_proposal_sd = if (isTRUE(estimate_selection_s)) selection_s_proposal_sd else NULL,
  selection_s_scale = "standardized_genotype_effect",
  updateLDswap = updateLDswap, ld_swap_prob = ld_swap_prob,
  ld_swap_r2 = ld_swap_r2,
  ld_swap_max_friends = as.integer(ld_swap_max_friends),
  ld_swap_moves = as.integer(ld_swap_moves)
 ), arch)
 fit
}

.fit_stblr_csr_bayesc <- function(...) {
 args <- list(...)
 args$method <- "bayesc"
 do.call(stblr_csr, args)
}

.make_bed_marker_data <- function(Glist, y, chr, cls, block_size,
                                  rows = NULL) {
  chr <- as.integer(chr)
  if (length(chr) < 1 || anyNA(chr)) {
    stop("chr must contain valid chromosome indices.")
  }

  if (is.null(dim(y))) {
    if (is.null(names(y)) && is.null(rows)) {
      stop(
        "When y is a vector, names(y) must contain IDs matching Glist$ids."
      )
    }

    y_names <- names(y)
    y <- matrix(as.numeric(y), ncol = 1)

    if (!is.null(y_names)) {
      rownames(y) <- y_names
    }
    colnames(y) <- "T1"
  } else {
    y <- as.matrix(y)
  }


  if (is.null(Glist$n) || !is.numeric(Glist$n) || length(Glist$n) != 1 ||
      !is.finite(Glist$n) || Glist$n < 1) {
    stop("Glist$n must be a finite positive scalar giving the total BED sample size.")
  }
  Glist$n <- as.integer(Glist$n)

  if (is.null(rows)) {
    if (!is.null(rownames(y))) {
      if (is.null(Glist$ids)) {
        stop("rownames(y) are present but Glist$ids is missing; cannot match y rows to BED individuals.")
      }

      ids_y <- rownames(y)
      ids_y <- as.character(ids_y)

      if (!is.character(Glist$ids)) {
        stop("Glist$ids must be a character vector in PLINK FAM/BED row order.")
      }
      if (length(Glist$ids) != Glist$n) {
        stop("length(Glist$ids) must equal Glist$n.")
      }

      if (anyDuplicated(ids_y)) {
        dup <- ids_y[duplicated(ids_y)]
        stop(
          "Phenotype IDs in rownames(y) must not contain duplicates. First duplicates: ",
          paste(head(unique(dup), 10), collapse = ", ")
        )
      }

      if (anyDuplicated(Glist$ids)) {
        dup <- Glist$ids[duplicated(Glist$ids)]
        stop(
          "Glist$ids must not contain duplicates when matching y to BED individuals. First duplicates: ",
          paste(head(unique(dup), 10), collapse = ", ")
        )
      }

      ids_g <- Glist$ids
      rows <- match(ids_y, ids_g)
      ok <- !is.na(rows)

      if (!all(ok)) {
        warning(
          sum(!ok),
          " phenotype IDs were not found in Glist$ids and will be dropped."
        )

        y <- y[ok, , drop = FALSE]
        ids_y <- ids_y[ok]
        rows <- rows[ok]
        rownames(y) <- ids_y
      }

      if (nrow(y) < 1L) {
        stop("No phenotype IDs matched Glist$ids.")
      }

      if (anyDuplicated(rows)) {
        dup_rows <- rows[duplicated(rows)]
        stop(
          "Duplicated BED rows after matching phenotype IDs. First duplicates: ",
          paste(head(unique(dup_rows), 10), collapse = ", ")
        )
      }
    } else if (nrow(y) == Glist$n) {
      rows <- NULL
    } else {
      stop(
        "nrow(y) != Glist$n, but rownames(y) are missing. ",
        "Set rownames(y) to IDs matching Glist$ids."
      )
    }
  } else {
    rows <- as.integer(rows)

    if (length(rows) != nrow(y)) {
      stop("length(rows) must equal nrow(y).")
    }
    if (anyNA(rows) || any(rows < 1L) || any(rows > Glist$n)) {
      stop("rows must be valid 1-based individual indices into the BED file.")
    }
    if (anyDuplicated(rows)) {
      dup_rows <- rows[duplicated(rows)]
      stop(
        "rows must not contain duplicate individual indices. First duplicates: ",
        paste(head(unique(dup_rows), 10), collapse = ", ")
      )
    }
  }

  if (is.null(Glist$bedfiles)) {
    stop("Glist$bedfiles is missing.")
  }

  bedfiles <- as.character(Glist$bedfiles)

  if (any(chr < 1L | chr > length(bedfiles))) {
    stop("chr contains indices outside Glist$bedfiles.")
  }

  missing_bed <- is.na(bedfiles[chr]) | !nzchar(bedfiles[chr])

  if (any(missing_bed)) {
    stop(
      "Glist$bedfiles is missing for chromosome/file index: ",
      paste(chr[missing_bed], collapse = ", ")
    )
  }

  if (is.null(cls)) {
    if (is.null(Glist$rsidsLD)) {
      stop("cls is NULL, but Glist$rsidsLD is missing.")
    }
    if (is.null(Glist$rsids)) {
      stop("cls is NULL, but Glist$rsids is missing.")
    }

    cls <- lapply(chr, function(cc) {
      if (length(Glist$rsidsLD) < cc || is.null(Glist$rsidsLD[[cc]])) {
        stop("Glist$rsidsLD[[", cc, "]] is missing.")
      }
      if (length(Glist$rsids) < cc || is.null(Glist$rsids[[cc]])) {
        stop("Glist$rsids[[", cc, "]] is missing.")
      }

      out <- match(Glist$rsidsLD[[cc]], Glist$rsids[[cc]])

      if (anyNA(out)) {
        missing <- Glist$rsidsLD[[cc]][is.na(out)]
        stop(
          "Some Glist$rsidsLD markers were not found in Glist$rsids for chromosome ",
          cc, ". First missing markers: ",
          paste(head(missing, 10), collapse = ", ")
        )
      }

      as.integer(out)
    })

    names(cls) <- paste0("chr", chr)
  } else if (!is.list(cls)) {
    cls <- list(cls)
  }

  if (length(cls) != length(chr)) {
    stop("cls must have one element per chromosome in chr.")
  }
  if (any(vapply(cls, anyNA, logical(1)))) {
    stop("cls must contain no missing marker indices.")
  }

  cls <- lapply(cls, as.integer)

  for (k in seq_along(chr)) {
    cc <- chr[k]

    if (length(Glist$rsids) < cc || is.null(Glist$rsids[[cc]])) {
      stop("Glist$rsids[[", cc, "]] is missing.")
    }

    if (any(cls[[k]] < 1L) || any(cls[[k]] > length(Glist$rsids[[cc]]))) {
      stop("cls[[", k, "]] contains marker indices outside Glist$rsids[[", cc, "]].")
    }
  }

  if (is.null(Glist$af)) {
    af <- vector("list", length(chr))
  } else {
    af <- Map(function(cc, cl) {
      if (length(Glist$af) < cc || is.null(Glist$af[[cc]])) {
        stop("Glist$af[[", cc, "]] is missing.")
      }
      if (length(Glist$af[[cc]]) < max(cl)) {
        stop("Glist$af[[", cc, "]] is shorter than the largest marker index in cls.")
      }
      out <- Glist$af[[cc]][cl]
      if (anyNA(out)) {
        stop("Glist$af[[", cc, "]][cls] contains missing values.")
      }
      out
    }, chr, cls)
  }

  m <- sum(lengths(cls))
  if (m < 1L) {
    stop("No markers selected. Check chr, cls, Glist$rsidsLD, and Glist$rsids.")
  }

  trait_names <- colnames(y)
  if (is.null(trait_names)) {
    trait_names <- paste0("T", seq_len(ncol(y)))
  }
  colnames(y) <- trait_names

  variable_names <- unlist(
    Map(function(cc, cl) Glist$rsids[[cc]][cl], chr, cls),
    use.names = FALSE
  )

  sets <- if (length(chr) == 1) {
    rep(seq_len(ceiling(m / block_size)), each = block_size)[seq_len(m)]
  } else {
    rep(seq_along(chr), lengths(cls))
  }

  list(
    y = y,
    chr = chr,
    bed_files = bedfiles[chr],
    cls = cls,
    af = af,
    rows = rows,
    m = m,
    nt = ncol(y),
    n = Glist$n,
    n_total = Glist$n,
    n_used = nrow(y),
    trait_names = trait_names,
    variable_names = variable_names,
    sets = sets,
    b_init = lapply(seq_len(ncol(y)), function(i) rep(0, m))
  )
}


#' Fit ST-BLR Directly from PLINK BED Markers
#'
#' Fits ST-BLR directly from markers stored in PLINK BED files referenced by a
#' qgg genotype list.
#'
#' @param Glist A qgg genotype list.
#' @param y Phenotype vector or matrix. If `y` is a vector, `names(y)` must
#'   contain individual IDs matching `Glist$ids`. If `y` is a matrix,
#'   `rownames(y)` must contain matching IDs unless `nrow(y) == Glist$n`.
#' @param chr Chromosome/file indices to fit. If `NULL`, all available
#'   non-missing, non-empty entries in `Glist$bedfiles` are used.
#' @param cls Optional marker column indices.
#' @param block_size Marker block size.
#' @param pi_init Initial marker inclusion probability. Defaults to 0.001.
#' @param pi_vb_init Inclusion probability used when initializing marker-effect
#'   variance. Defaults to `pi_init`.
#' @param pi_prior_mean Prior mean for the marker inclusion probability.
#'   Defaults to 0.001.
#' @param pi_prior_strength Total Beta prior strength for the marker inclusion
#'   probability. Defaults to 5e5.
#' @param pi_prior_a,pi_prior_b Optional explicit Beta prior shape parameters.
#'   When supplied, these override `pi_prior_mean` and `pi_prior_strength`.
#' @param h2 Initial heritability. Defaults to 0.3.
#' @param nub,nue Prior degrees of freedom.
#' @param scale Standardize BED markers.
#' @param updateB,updateE,updatePi Logical sampler update controls.
#' @param adjE Residual adjustment factor.
#' @param nit,nburn,nthin MCMC iteration controls.
#' @param rebuild_every,full_sweep_every,candidate_threshold,candidate_lifetime
#'   Sampler scheduling and residual-rebuild controls.
#' @param skip_nulls_burnin_only Restrict null-marker skipping to burn-in.
#' @param ncores Number of OpenMP threads.
#' @param seed Sampler seed.
#' @param null_update_prob,null_skip_base,null_skip_max Null-marker scheduling
#'   controls.
#' @param return_wy,return_r Return optional sampler state.
#' @param rows Optional 1-based individual row indices into the BED file. When
#'   omitted, subset phenotypes can instead be matched using `rownames(y)` and
#'   `Glist$ids`.
#' @param backend Sampler backend. `"auto"` selects the scheduled multi-chain
#'   backend for multi-chromosome fits or when `nchains > 1`, and the scheduled
#'   single-chain backend otherwise. `"scheduled"` forces a scheduled backend;
#'   for multi-chromosome fits it uses the scheduled multi-chain backend.
#'   `"sparse"` is only valid for one chromosome.
#' @param nchains Number of MCMC chains. When `nchains > 1`, the scheduled
#'   multi-chain backend is used and returns chain-stability summaries
#'   `dm_sd`, `dm_min`, `dm_max`, `bm_sd`, `bm_min`, and `bm_max`.
#' @param read_block_size,progress_every Multi-chain BED reading and progress
#'   controls.
#' @return A formatted ST-BLR fit. Fits from the scheduled multi-chain backend
#'   include CSR-compatible chain-stability summaries for marker inclusion
#'   probabilities and posterior effects.
#' @export
stblr_bed_marker <- function(
    Glist, y, chr = NULL, cls = NULL, block_size = 1000,
    pi_init = 0.001, pi_vb_init = NULL, pi_prior_mean = 0.001,
    pi_prior_strength = 5e5, pi_prior_a = NULL, pi_prior_b = NULL,
    h2 = 0.3, nub = 4, nue = 4, scale = TRUE, updateB = TRUE,
    updateE = TRUE, updatePi = TRUE, adjE = 0, nit = 1000, nburn = 100,
    nthin = 1, rebuild_every = 25, full_sweep_every = 10,
    candidate_threshold = 1e-3, candidate_lifetime = 20,
    skip_nulls_burnin_only = FALSE, ncores = 1, seed = 10,
    backend = c("auto", "scheduled", "sparse"),
    null_update_prob = 0.02, null_skip_base = 50,
    null_skip_max = 200, return_wy = FALSE, return_r = FALSE, rows = NULL,
    nchains = 1, read_block_size = 64, progress_every = 0
) {

  backend <- match.arg(backend)
  if (!is.numeric(nchains) || length(nchains) != 1 ||
      !is.finite(nchains) || nchains < 1) {
    stop("nchains must be a finite positive scalar.")
  }
  nchains <- as.integer(nchains)

  if (!is.numeric(ncores) || length(ncores) != 1 ||
      !is.finite(ncores) || ncores < 1) {
    stop("ncores must be a finite positive scalar.")
  }
  ncores <- as.integer(ncores)

  if (is.null(Glist$bedfiles)) {
    stop("Glist$bedfiles is missing.")
  }

  bedfiles <- as.character(Glist$bedfiles)
  has_bedfile <- !is.na(bedfiles) & nzchar(bedfiles)

  if (is.null(chr)) {
    chr <- which(has_bedfile)

    if (length(chr) < 1L) {
      stop("No available BED files found in Glist$bedfiles.")
    }
  } else {
    chr <- as.integer(chr)

    if (length(chr) < 1L || anyNA(chr)) {
      stop("chr must contain valid chromosome/file indices.")
    }

    if (any(chr < 1L | chr > length(bedfiles))) {
      stop("chr contains indices outside Glist$bedfiles.")
    }
  }

  arch <- .resolve_pi_prior(
    pi_init = pi_init,
    pi_vb_init = pi_vb_init,
    pi_prior_mean = pi_prior_mean,
    pi_prior_strength = pi_prior_strength,
    pi_prior_a = pi_prior_a,
    pi_prior_b = pi_prior_b
  )

 dat <- .make_bed_marker_data(
   Glist = Glist,
   y = y,
   chr = chr,
   cls = cls,
   block_size = block_size,
   rows = rows
 )

 use_chains_backend <- FALSE

 if (backend == "auto") {
   use_chains_backend <- length(dat$chr) > 1L || nchains > 1L
   use_scheduled_backend <- !use_chains_backend
   use_sparse_backend <- FALSE
 } else if (backend == "scheduled") {
   use_chains_backend <- nchains > 1L || length(dat$chr) > 1L
   use_scheduled_backend <- !use_chains_backend
   use_sparse_backend <- FALSE
 } else {
   use_chains_backend <- FALSE
   use_scheduled_backend <- FALSE
   use_sparse_backend <- TRUE
 }

 if (use_sparse_backend && length(dat$chr) != 1L) {
   stop("backend = 'sparse' only supports one chromosome.")
 }

 pri <- .make_stblr_priors(
  dat$y, dat$m, h2, nub, nue, arch$pi_vb_init, arch$pi_prior_mean,
  dat$trait_names
 )
 # n is the total BED sample size needed for byte layout; rows selects the
 # individuals used, and n_used is the resulting phenotype/sample size.
 common <- list(
   bed_files = dat$bed_files,
   n = dat$n_total,
   cls = dat$cls,
   y = dat$y,
   b_init = dat$b_init,
   sets = dat$sets,
   rows = dat$rows,
   af = dat$af,
   scale = scale,
   B = pri$B,
   E = pri$E,
   ssb_prior = pri$ssb_prior_list,
   sse_prior = pri$sse_prior_list,
   pi = arch$pi,
   nub = nub,
   nue = nue,
   updateB = updateB,
   updateE = updateE,
   updatePi = updatePi,
   adjE = adjE,
   nit = nit,
   nburn = nburn,
   nthin = nthin,
   rebuild_every = rebuild_every,
   full_sweep_every = full_sweep_every,
   candidate_threshold = candidate_threshold,
   candidate_lifetime = candidate_lifetime,
   skip_nulls_burnin_only = skip_nulls_burnin_only
 )
 if (use_chains_backend) {
   raw <- do.call(stblr_cpg_omp_bed_marker_scheduled_chains, c(common, list(
     null_skip_base = null_skip_base,
     null_skip_max = null_skip_max,
     return_wy = return_wy,
     return_r = return_r,
     read_block_size = as.integer(read_block_size),
     progress_every = as.integer(progress_every),
     pi_prior_a = arch$pi_prior_a,
     pi_prior_b = arch$pi_prior_b,
     nchains = as.integer(nchains),
     ncores = ncores,
     seed = seed
   )))
 } else if (use_scheduled_backend) {
   raw <- do.call(stblr_cpg_omp_bed_marker_scheduled, c(common, list(
     null_skip_base = null_skip_base,
     null_skip_max = null_skip_max,
     return_wy = return_wy,
     return_r = return_r,
     pi_prior_a = arch$pi_prior_a,
     pi_prior_b = arch$pi_prior_b,
     ncores = ncores,
     seed = seed
   )))
 } else {
   raw <- do.call(stblr_cpg_omp_bed_marker_sparse, c(common, list(
     null_update_prob = null_update_prob,
     ncores = ncores,
     seed = seed
   )))
 }
 if (.is_stblr_raw_v1(raw)) {
  fit <- .format_stblr_raw_v1(raw, dat$trait_names, dat$variable_names)
 } else {
  fit <- .format_stblr_fit(
   raw, dat$nt, dat$m, dat$trait_names, dat$variable_names, TRUE
  )
 }
 fit$input <- c(list(
  chr = dat$chr, cls = dat$cls, n = dat$n, n_total = dat$n_total,
  n_used = dat$n_used, m = dat$m,
  nt = dat$nt, block_size = block_size, sets = dat$sets, h2 = h2, nub = nub,
  nue = nue, vy = pri$vy, B = pri$B, E = pri$E, ssb_prior = pri$ssb_prior,
  sse_prior = pri$sse_prior, updateB = updateB, updateE = updateE,
  updatePi = updatePi, adjE = adjE, nit = nit, nburn = nburn, nthin = nthin,
  rebuild_every = rebuild_every, full_sweep_every = full_sweep_every,
  candidate_threshold = candidate_threshold, candidate_lifetime = candidate_lifetime,
  skip_nulls_burnin_only = skip_nulls_burnin_only,
  null_update_prob = null_update_prob, null_skip_base = null_skip_base,
  null_skip_max = null_skip_max, ncores = ncores, seed = seed,
  backend = backend,
  use_chains_backend = use_chains_backend,
  use_scheduled_backend = use_scheduled_backend,
  use_sparse_backend = use_sparse_backend,
  nchains = nchains,
  read_block_size = read_block_size, progress_every = progress_every,
  scale = scale, rows = dat$rows
 ), arch)
 fit
}

.normalize_stblr_method <- function(method) {
 method <- tolower(method)
 if (length(method) > 1L) method <- method[1L]
 if (length(method) < 1L || is.na(method) ||
     !method %in% c("bayesc", "bayesr")) {
  stop("method must be one of 'bayesc' or 'bayesr'.")
 }
 method
}

#' Fit individual-level PLINK BED ST-BLR models
#'
#' Fits individual-level ST-BLR models directly from PLINK BED files referenced
#' by a BED-backed `Glist`. This is the individual-level analogue of
#' [stblr_csr()] for the plain BayesC and BayesR model families.
#'
#' `method = "bayesc"` fits the BayesC BED scheduled-chain backend.
#' `method = "bayesr"` fits the BayesR BED scheduled-chain backend. For BayesR,
#' `dm` is `P(component > 0)`, `comp_prob` contains marker-by-component
#' posterior probabilities, and `dm_component_mean` contains posterior mean
#' component indices. BED scheduled-chain fits expose CPO/log-CPO diagnostics.
#'
#' @param y Phenotype vector or matrix. If `y` is a vector, `names(y)` must
#'   contain individual IDs matching `Glist$ids`. If `y` is a matrix,
#'   `rownames(y)` must contain matching IDs unless `nrow(y) == Glist$n`.
#' @param Glist A BED-backed qgg genotype list containing `bedfiles`, `n`,
#'   `ids`, `rsids`, `rsidsLD`, and optionally `af`.
#' @param method Model family. Case-insensitive `"bayesc"` fits BayesC and is
#'   the default; case-insensitive `"bayesr"` fits BayesR.
#' @param covar Optional covariates. Covariates are not currently adjusted by
#'   this BED interface; pass pre-adjusted phenotypes.
#' @param chr Chromosome/file indices to fit. If `NULL`, all available
#'   non-missing, non-empty entries in `Glist$bedfiles` are used.
#' @param cls Optional marker column indices.
#' @param block_size Marker block size used to construct update sets for a
#'   single BED file.
#' @param rows Optional 1-based individual row indices into the BED file. When
#'   omitted, subset phenotypes can instead be matched using `rownames(y)` and
#'   `Glist$ids`.
#' @param scale Standardize BED markers using allele frequencies.
#' @param nit,nburn,nthin MCMC iteration controls.
#' @param ncores Number of OpenMP threads.
#' @param seed Sampler seed.
#' @param nchains Number of scheduled MCMC chains.
#' @param chain_seeds Optional explicit chain seeds. Not currently supported by
#'   the BED scheduled-chain backends.
#' @param pi_init Initial BayesC marker inclusion probability.
#' @param pi_vb_init Inclusion probability used when initializing BayesC
#'   marker-effect variance. Defaults to `pi_init`.
#' @param pi_prior_mean Prior mean for the BayesC marker inclusion probability.
#' @param pi_prior_strength Total Beta prior strength for the BayesC marker
#'   inclusion probability.
#' @param pi_prior_a,pi_prior_b Optional explicit BayesC Beta prior shape
#'   parameters. When supplied, these override `pi_prior_mean` and
#'   `pi_prior_strength`.
#' @param h2 Initial heritability.
#' @param adjE Residual adjustment factor.
#' @param updateB,updateE,updatePi Logical sampler update controls.
#' @param mixture_var BayesR mixture variance multipliers. The first value must
#'   be zero and non-null values must be positive.
#' @param pi Initial BayesR mixture probabilities. Defaults to total non-null
#'   probability 0.001 split equally across non-null components.
#' @param alpha BayesR Dirichlet prior shape vector for mixture probabilities.
#' @param nub,nue Prior degrees of freedom.
#' @param rebuild_every,full_sweep_every,candidate_threshold,candidate_lifetime
#'   Sampler scheduling and residual-rebuild controls.
#' @param skip_nulls_burnin_only Restrict null-marker skipping to burn-in.
#' @param null_skip_base,null_skip_max Null-marker scheduling controls.
#' @param return_wy,return_r Return optional sampler state.
#' @param read_block_size,progress_every Scheduled-chain BED reading and
#'   progress controls.
#' @param ... Unsupported arguments.
#' @return A formatted ST-BLR fit with standard `dm`, `bm`, and `input` fields.
#'   BayesC fits include `log_cpo`, `mean_log_cpo`, `final_pi`, `mean_pi`, and
#'   chain-summary matrices. BayesR fits additionally include `comp_prob` and
#'   `dm_component_mean`.
#'
#' @examples
#' \dontrun{
#' fitC_bed <- stblr_bed(
#'   y = y,
#'   Glist = Glist,
#'   method = "bayesC",
#'   nit = 1000,
#'   nburn = 100
#' )
#'
#' fitR_bed <- stblr_bed(
#'   y = y,
#'   Glist = Glist,
#'   method = "bayesR",
#'   nit = 1000,
#'   nburn = 100
#' )
#' }
#' @export
stblr_bed <- function(
  y,
  Glist,
  method = c("bayesc", "bayesr"),
  covar = NULL,
  chr = NULL,
  cls = NULL,
  block_size = 1000,
  rows = NULL,
  scale = TRUE,
  nit = 1000,
  nburn = 100,
  nthin = 1,
  ncores = 1,
  seed = 1,
  nchains = 1L,
  chain_seeds = NULL,
  pi_init = 0.001,
  pi_vb_init = NULL,
  pi_prior_mean = 0.001,
  pi_prior_strength = 5e5,
  pi_prior_a = NULL,
  pi_prior_b = NULL,
  h2 = 0.3,
  adjE = 0.9,
  updateB = TRUE,
  updateE = TRUE,
  updatePi = TRUE,
  mixture_var = c(0, 0.01, 0.1, 1),
  pi = NULL,
  alpha = NULL,
  nub = 4,
  nue = 4,
  rebuild_every = 25,
  full_sweep_every = 10,
  candidate_threshold = 1e-3,
  candidate_lifetime = 20,
  skip_nulls_burnin_only = FALSE,
  null_skip_base = 50,
  null_skip_max = 200,
  return_wy = FALSE,
  return_r = FALSE,
  read_block_size = 64,
  progress_every = 0,
  ...
) {
 method <- .normalize_stblr_method(method)
 dots <- list(...)
 if (length(dots)) {
  dot_names <- names(dots)
  if (is.null(dot_names)) dot_names <- rep("", length(dots))
  dot_names <- dot_names[nzchar(dot_names)]
  if (length(dot_names)) {
   stop("Unsupported argument(s) in ...: ", paste(dot_names, collapse = ", "))
  }
  stop("Unsupported unnamed argument(s) in ....")
 }
 if (!is.null(covar)) {
  stop("covar is not currently supported by stblr_bed(); pass pre-adjusted phenotypes.")
 }
 if (!is.null(chain_seeds)) {
  stop("chain_seeds is not currently supported by the BED scheduled-chain backends.")
 }
 if (!is.numeric(nchains) || length(nchains) != 1L ||
     !is.finite(nchains) || nchains < 1 || nchains != floor(nchains)) {
  stop("nchains must be a positive integer scalar.")
 }
 nchains <- as.integer(nchains)
 if (!is.numeric(ncores) || length(ncores) != 1L ||
     !is.finite(ncores) || ncores < 1 || ncores != floor(ncores)) {
  stop("ncores must be a positive integer scalar.")
 }
 ncores <- as.integer(ncores)

 if (is.null(Glist$bedfiles)) {
  stop("Glist must contain BED file information in Glist$bedfiles for stblr_bed().")
 }
 bedfiles <- as.character(Glist$bedfiles)
 has_bedfile <- !is.na(bedfiles) & nzchar(bedfiles)
 if (is.null(chr)) {
  chr <- which(has_bedfile)
  if (length(chr) < 1L) {
   stop("No available BED files found in Glist$bedfiles.")
  }
 } else {
  chr <- as.integer(chr)
  if (length(chr) < 1L || anyNA(chr)) {
   stop("chr must contain valid chromosome/file indices.")
  }
  if (any(chr < 1L | chr > length(bedfiles))) {
   stop("chr contains indices outside Glist$bedfiles.")
  }
 }

 bayesc_only <- character()
 if (method == "bayesr") {
  if (!missing(pi_init) && !identical(pi_init, 0.001)) {
   bayesc_only <- c(bayesc_only, "pi_init")
  }
  if (!missing(pi_vb_init) && !is.null(pi_vb_init)) {
   bayesc_only <- c(bayesc_only, "pi_vb_init")
  }
  if (!missing(pi_prior_mean) && !identical(pi_prior_mean, 0.001)) {
   bayesc_only <- c(bayesc_only, "pi_prior_mean")
  }
  if (!missing(pi_prior_strength) && !identical(pi_prior_strength, 5e5)) {
   bayesc_only <- c(bayesc_only, "pi_prior_strength")
  }
  if (!missing(pi_prior_a) && !is.null(pi_prior_a)) {
   bayesc_only <- c(bayesc_only, "pi_prior_a")
  }
  if (!missing(pi_prior_b) && !is.null(pi_prior_b)) {
   bayesc_only <- c(bayesc_only, "pi_prior_b")
  }
  if (length(bayesc_only)) {
   stop(
    "`", paste(bayesc_only, collapse = "`, `"),
    "` are BayesC-specific. Use `pi` and `alpha` for method = 'bayesr'."
   )
  }
 } else {
  bayesr_only <- character()
  if (!missing(mixture_var)) bayesr_only <- c(bayesr_only, "mixture_var")
  if (!missing(pi) && !is.null(pi)) bayesr_only <- c(bayesr_only, "pi")
  if (!missing(alpha) && !is.null(alpha)) bayesr_only <- c(bayesr_only, "alpha")
  if (length(bayesr_only)) {
   stop(
    "`", paste(bayesr_only, collapse = "`, `"),
    "` are BayesR-specific and are only supported with method = 'bayesr'."
   )
  }
 }

 dat <- .make_bed_marker_data(
  Glist = Glist,
  y = y,
  chr = chr,
  cls = cls,
  block_size = block_size,
  rows = rows
 )

 if (method == "bayesr") {
  if (!is.numeric(mixture_var) || length(mixture_var) < 2L ||
      any(!is.finite(mixture_var)) || mixture_var[1L] != 0 ||
      any(mixture_var[-1L] <= 0)) {
   stop("mixture_var must start with 0 and have positive non-null components.")
  }
  if (is.null(pi)) {
   active <- 0.001
   pi <- c(1 - active, rep(active / (length(mixture_var) - 1L),
                          length(mixture_var) - 1L))
  }
  if (!is.numeric(pi) || length(pi) != length(mixture_var) ||
      any(!is.finite(pi)) || any(pi < 0) || sum(pi) <= 0) {
   stop("pi must be a non-negative finite vector matching mixture_var with positive sum.")
  }
  pi <- pi / sum(pi)
  if (is.null(alpha)) alpha <- pmax(pi * 5e5, .Machine$double.eps)
  if (!is.numeric(alpha) || length(alpha) != length(mixture_var) ||
      any(!is.finite(alpha)) || any(alpha <= 0)) {
   stop("alpha must be a positive finite vector matching mixture_var.")
  }
  pi_active <- sum(pi[-1L])
  if (!is.finite(pi_active) || pi_active <= 0) {
   stop("pi must allocate positive probability to non-null components.")
  }
  pri <- .make_stblr_priors(
   dat$y, dat$m, h2, nub, nue, pi_active, pi_active, dat$trait_names
  )
  fit <- .fit_stblr_bed_bayesr(
   bed_files = dat$bed_files,
   n = dat$n_total,
   cls = dat$cls,
   y = dat$y,
   b_init = dat$b_init,
   sets = dat$sets,
   rows = dat$rows,
   af = dat$af,
   scale = scale,
   B = pri$B,
   E = pri$E,
   ssb_prior = pri$ssb_prior_list,
   sse_prior = pri$sse_prior_list,
   pi = pi,
   mixture_var = mixture_var,
   alpha = alpha,
   nub = nub,
   nue = nue,
   updateB = updateB,
   updateE = updateE,
   updatePi = updatePi,
   adjE = adjE,
   nit = nit,
   nburn = nburn,
   nthin = nthin,
   rebuild_every = rebuild_every,
   full_sweep_every = full_sweep_every,
   null_skip_base = null_skip_base,
   null_skip_max = null_skip_max,
   candidate_threshold = candidate_threshold,
   candidate_lifetime = candidate_lifetime,
   skip_nulls_burnin_only = skip_nulls_burnin_only,
   return_wy = return_wy,
   return_r = return_r,
   read_block_size = read_block_size,
   progress_every = progress_every,
   nchains = nchains,
   ncores = ncores,
   seed = seed,
   trait_names = dat$trait_names,
   variable_names = dat$variable_names,
   keep_diagnostics = TRUE
  )
  fit$input <- c(list(
   method = "bayesr",
   model = "bayesr",
   backend = "bed_bayesr",
   data_level = "individual",
   annotations = FALSE,
   chr = dat$chr,
   cls = dat$cls,
   n = dat$n,
   n_total = dat$n_total,
   n_used = dat$n_used,
   m = dat$m,
   nt = dat$nt,
   block_size = block_size,
   sets = dat$sets,
   h2 = h2,
   nub = nub,
   nue = nue,
   vy = pri$vy,
   B = pri$B,
   E = pri$E,
   ssb_prior = pri$ssb_prior,
   sse_prior = pri$sse_prior,
   updateB = updateB,
   updateE = updateE,
   updatePi = updatePi,
   adjE = adjE,
   nit = as.integer(nit),
   nburn = as.integer(nburn),
   nthin = as.integer(nthin),
   rebuild_every = rebuild_every,
   full_sweep_every = full_sweep_every,
   candidate_threshold = candidate_threshold,
   candidate_lifetime = candidate_lifetime,
   skip_nulls_burnin_only = skip_nulls_burnin_only,
   null_skip_base = null_skip_base,
   null_skip_max = null_skip_max,
   ncores = ncores,
   seed = as.integer(seed),
   nchains = nchains,
   read_block_size = read_block_size,
   progress_every = progress_every,
   scale = scale,
   rows = dat$rows
  ), fit$input)
  fit$input$method <- "bayesr"
  fit$input$model <- "bayesr"
  fit$input$backend <- "bed_bayesr"
  fit$input$data_level <- "individual"
  fit$input$nchains <- nchains
  fit$input <- fit$input[!duplicated(names(fit$input))]
  return(fit)
 }

 arch <- .resolve_pi_prior(
  pi_init = pi_init,
  pi_vb_init = pi_vb_init,
  pi_prior_mean = pi_prior_mean,
  pi_prior_strength = pi_prior_strength,
  pi_prior_a = pi_prior_a,
  pi_prior_b = pi_prior_b
 )
 pri <- .make_stblr_priors(
  dat$y, dat$m, h2, nub, nue, arch$pi_vb_init, arch$pi_prior_mean,
  dat$trait_names
 )
 common <- list(
  bed_files = dat$bed_files,
  n = dat$n_total,
  cls = dat$cls,
  y = dat$y,
  b_init = dat$b_init,
  sets = dat$sets,
  rows = dat$rows,
  af = dat$af,
  scale = scale,
  B = pri$B,
  E = pri$E,
  ssb_prior = pri$ssb_prior_list,
  sse_prior = pri$sse_prior_list,
  pi = arch$pi,
  nub = nub,
  nue = nue,
  updateB = updateB,
  updateE = updateE,
  updatePi = updatePi,
  adjE = adjE,
  nit = nit,
  nburn = nburn,
  nthin = nthin,
  rebuild_every = rebuild_every,
  full_sweep_every = full_sweep_every,
  null_skip_base = null_skip_base,
  null_skip_max = null_skip_max,
  candidate_threshold = candidate_threshold,
  candidate_lifetime = candidate_lifetime,
  skip_nulls_burnin_only = skip_nulls_burnin_only,
  return_wy = return_wy,
  return_r = return_r,
  read_block_size = as.integer(read_block_size),
  progress_every = as.integer(progress_every),
  pi_prior_a = arch$pi_prior_a,
  pi_prior_b = arch$pi_prior_b,
  nchains = nchains,
  ncores = ncores,
  seed = seed
 )
 raw <- do.call(stblr_cpg_omp_bed_marker_scheduled_chains, common)
 if (.is_stblr_raw_v1(raw)) {
  fit <- .format_stblr_raw_v1(raw, dat$trait_names, dat$variable_names)
 } else {
  fit <- .format_stblr_fit(
   raw, dat$nt, dat$m, dat$trait_names, dat$variable_names, TRUE
  )
 }
 if (is.null(fit$final_pi)) fit$final_pi <- fit$pi
 if (is.null(fit$mean_pi)) fit$mean_pi <- fit$pim
 fit$input <- c(list(
  method = "bayesc",
  model = "bayesc",
  backend = "bed_bayesc",
  data_level = "individual",
  annotations = FALSE,
  scheduled = TRUE,
  keep_chains = FALSE,
  chr = dat$chr,
  cls = dat$cls,
  n = dat$n,
  n_total = dat$n_total,
  n_used = dat$n_used,
  m = dat$m,
  nt = dat$nt,
  block_size = block_size,
  sets = dat$sets,
  h2 = h2,
  nub = nub,
  nue = nue,
  vy = pri$vy,
  B = pri$B,
  E = pri$E,
  ssb_prior = pri$ssb_prior,
  sse_prior = pri$sse_prior,
  updateB = updateB,
  updateE = updateE,
  updatePi = updatePi,
  adjE = adjE,
  nit = as.integer(nit),
  nburn = as.integer(nburn),
  nthin = as.integer(nthin),
  rebuild_every = rebuild_every,
  full_sweep_every = full_sweep_every,
  candidate_threshold = candidate_threshold,
  candidate_lifetime = candidate_lifetime,
  skip_nulls_burnin_only = skip_nulls_burnin_only,
  null_skip_base = null_skip_base,
  null_skip_max = null_skip_max,
  ncores = ncores,
  seed = as.integer(seed),
  nchains = nchains,
  read_block_size = read_block_size,
  progress_every = progress_every,
  scale = scale,
  rows = dat$rows
 ), arch)
 fit
}

.fit_stblr_bed_bayesc <- function(...) {
 args <- list(...)
 args$method <- "bayesc"
 do.call(stblr_bed, args)
}



#' Make Sufficient Statistics for ST-BLR
#'
#' Computes genome-wide marker sufficient statistics from PLINK BED files
#' referenced by a `Glist` object. The returned marker order is consistent with
#' the sparse-LD CSR matrix produced by [sparseLD_stream_CSR()] when using the
#' same `chr` and `cls`.
#'
#' @param Glist A qgg genotype list containing `bedfiles`, `n`, `ids`,
#'   `rsids`, `rsidsLD`, and optionally `af`.
#' @param y Phenotype vector or matrix. If `y` is a vector, `names(y)` must
#'   contain individual IDs matching `Glist$ids`. If `y` is a matrix,
#'   `rownames(y)` must contain individual IDs matching `Glist$ids`, unless
#'   `rows` is supplied.
#' @param chr Chromosome/file indices to use. If `NULL`, all available
#'   non-missing, non-empty entries in `Glist$bedfiles` are used.
#' @param cls Optional marker column indices. If `NULL`, these are inferred as
#'   `match(Glist$rsidsLD[[cc]], Glist$rsids[[cc]])` for each chromosome/file.
#'   If supplied as a vector, it is treated as the marker index vector for a
#'   single chromosome/file. If supplied as a list, it must have one element per
#'   chromosome/file in `chr`.
#' @param rows Optional 1-based BED/FAM row indices corresponding to rows of
#'   `y`. If `NULL`, rows are inferred by matching `rownames(y)` or `names(y)`
#'   to `Glist$ids`.
#' @param scale Logical. Standardize BED markers using allele frequencies.
#' @param nthreads Number of OpenMP threads used by [bed_xtx_xty()].
#'
#' @return A list with genome-wide sufficient statistics:
#' \describe{
#'   \item{wy}{List of marker-trait cross-products, one vector per trait.}
#'   \item{ww}{List of marker cross-products, one vector per trait.}
#'   \item{yy}{Trait sums of squares.}
#'   \item{n}{Analysis sample size after phenotype-to-BED matching.}
#'   \item{m}{Total number of markers.}
#'   \item{chr}{Chromosome/file indices used.}
#'   \item{bed_files}{BED files used, in the same order as `chr`.}
#'   \item{cls}{Marker indices used for each chromosome/file.}
#'   \item{af}{Allele frequencies corresponding to `cls`.}
#'   \item{rows}{BED/FAM row indices used.}
#'   \item{marker_names}{Genome-wide marker names in the same order as `wy`
#'     and `ww`.}
#'   \item{trait_names}{Trait names.}
#'   \item{stats_by_chr}{Per-chromosome sufficient statistics returned by
#'     [bed_xtx_xty()].}
#' }
#'
#' @export
make_summary_stats <- function(Glist, y, chr = NULL, cls = NULL, rows = NULL,
                       scale = TRUE, nthreads = 1) {
  bedfiles <- as.character(Glist$bedfiles)
  has_bedfile <- !is.na(bedfiles) & nzchar(bedfiles)

  if (is.null(chr)) {
    chr <- which(has_bedfile)
  } else {
    chr <- as.integer(chr)
  }

  if (length(chr) < 1L || anyNA(chr)) {
    stop("chr must contain valid chromosome/file indices.")
  }

  if (any(chr < 1L | chr > length(bedfiles))) {
    stop("chr contains indices outside Glist$bedfiles.")
  }

  missing_bed <- is.na(bedfiles[chr]) | !nzchar(bedfiles[chr])

  if (any(missing_bed)) {
    stop(
      "Glist$bedfiles is missing for chromosome/file index: ",
      paste(chr[missing_bed], collapse = ", ")
    )
  }

  if (is.null(dim(y))) {
    ids_y <- names(y)

    y <- matrix(as.numeric(y), ncol = 1)

    if (!is.null(ids_y)) {
      rownames(y) <- ids_y
    }

    colnames(y) <- "T1"
  } else {
    y <- as.matrix(y)
  }

  if (is.null(rows) && is.null(rownames(y))) {
    stop(
      "When rows is NULL, y must have names(y) or rownames(y) ",
      "matching Glist$ids."
    )
  }

  if (is.null(rows)) {
    rows <- match(rownames(y), Glist$ids)
    ok <- !is.na(rows)

    if (!all(ok)) {
      warning(
        sum(!ok),
        " phenotype IDs were not found in Glist$ids and will be dropped."
      )

      y <- y[ok, , drop = FALSE]
      rows <- rows[ok]
    }

    if (nrow(y) < 1L) {
      stop("No phenotype IDs matched Glist$ids.")
    }
  } else {
    rows <- as.integer(rows)

    if (length(rows) != nrow(y)) {
      stop("length(rows) must equal nrow(y).")
    }
  }

  rows <- as.integer(rows)

  ord <- order(rows)
  rows <- rows[ord]
  y <- y[ord, , drop = FALSE]

  trait_names <- colnames(y)
  if (is.null(trait_names)) {
    trait_names <- paste0("T", seq_len(ncol(y)))
  }
  colnames(y) <- trait_names

  if (is.null(cls)) {
    cls <- lapply(chr, function(cc) {
      match(Glist$rsidsLD[[cc]], Glist$rsids[[cc]])
    })
  } else if (!is.list(cls)) {
    cls <- list(cls)
  }

  if (length(cls) != length(chr)) {
    stop("cls must have one element per chromosome/file in chr.")
  }

  cls <- lapply(cls, as.integer)
  names(cls) <- paste0("chr", chr)

  af <- Map(function(cc, cl) {
    Glist$af[[cc]][cl]
  }, chr, cls)

  stats_by_chr <- vector("list", length(chr))
  names(stats_by_chr) <- paste0("chr", chr)

  for (k in seq_along(chr)) {
    cc <- chr[k]

    message("Computing sufficient statistics for chromosome ", cc)

    stats_by_chr[[k]] <- bed_xtx_xty(
      bed_file = bedfiles[cc],
      n = Glist$n,
      cls = cls[[k]],
      af = af[[k]],
      y = y,
      rows = rows,
      scale = scale,
      nthreads = nthreads
    )
  }

  nt <- ncol(y)

  marker_names <- unlist(
    Map(function(cc, cl) Glist$rsids[[cc]][cl], chr, cls),
    use.names = FALSE
  )

  wy <- lapply(seq_len(nt), function(t) {
    out <- unlist(
      lapply(stats_by_chr, function(s) s$wy[[t]]),
      use.names = FALSE
    )
    names(out) <- marker_names
    out
  })

  ww <- lapply(seq_len(nt), function(t) {
    out <- unlist(
      lapply(stats_by_chr, function(s) s$ww[[t]]),
      use.names = FALSE
    )
    names(out) <- marker_names
    out
  })

  names(wy) <- trait_names
  names(ww) <- trait_names

  yy <- stats_by_chr[[1]]$yy
  names(yy) <- trait_names

  list(
    wy = wy,
    ww = ww,
    yy = yy,
    n = nrow(y),
    m = length(marker_names),
    chr = chr,
    bed_files = bedfiles[chr],
    cls = cls,
    af = af,
    rows = rows,
    marker_names = marker_names,
    trait_names = trait_names,
    stats_by_chr = stats_by_chr
  )
}


#' Make Sparse LD for ST-BLR
#'
#' Computes a disk-backed sparse-LD CSR matrix from PLINK BED files referenced
#' by a `Glist` object and stores the sparse-LD prefix and resolved marker
#' structure in `Glist$sparseLD`. The returned `Glist` can be used directly
#' with [make_summary_stats()] and [stblr_csr()].
#'
#' @param Glist A qgg genotype list containing `bedfiles`, `n`, `ids`,
#'   `rsids`, `rsidsLD`, `af`, and optionally `idsLD`.
#' @param rows Optional 1-based BED/FAM row indices used as the LD reference
#'   sample. If `NULL` and `Glist$idsLD` is available, rows are inferred as
#'   `match(Glist$idsLD, Glist$ids)`. If still `NULL`, all BED individuals are
#'   used.
#' @param out_prefix Output prefix for the disk-backed sparse-LD CSR files. If
#'   `NULL`, defaults to `"sparseLD"` in the directory of the first selected
#'   BED file.
#' @param chr Chromosome/file indices to use. If `NULL`, all available
#'   non-missing, non-empty entries in `Glist$bedfiles` are used.
#' @param cls Optional marker column indices. If `NULL`, these are inferred as
#'   `match(Glist$rsidsLD[[cc]], Glist$rsids[[cc]])` for each chromosome/file.
#'   If supplied as a vector, it is treated as the marker index vector for a
#'   single chromosome/file. If supplied as a list, it must have one element per
#'   chromosome/file in `chr`.
#' @param pos_bp Optional marker base-pair positions passed to
#'   [sparseLD_stream_CSR()].
#' @param max_distance_bp Maximum base-pair distance between retained marker
#'   pairs. Use 0 to disable base-pair distance filtering.
#' @param max_distance_variants Maximum marker-index distance between retained
#'   marker pairs. Use 0 to disable marker-index distance filtering. At least
#'   one of `max_distance_bp` or `max_distance_variants` should usually be
#'   positive.
#' @param r2_threshold Minimum squared-correlation threshold for retaining LD
#'   entries.
#' @param block_size Marker block size used by [sparseLD_stream_CSR()].
#' @param nthreads Number of OpenMP threads used by [sparseLD_stream_CSR()].
#' @param allow_full_ld Logical. If `FALSE`, the function stops when both
#'   distance filters are disabled, because this evaluates all marker pairs and
#'   may be extremely slow. Set to `TRUE` to explicitly allow full pairwise LD.
#'
#' @return The input `Glist` with an added `sparseLD` element containing:
#' \describe{
#'   \item{prefix}{Output prefix for the sparse-LD CSR files.}
#'   \item{chr}{Chromosome/file indices used.}
#'   \item{bed_files}{BED files used, in the same order as `chr`.}
#'   \item{cls}{Marker indices used for each chromosome/file.}
#'   \item{af}{Allele frequencies corresponding to `cls`.}
#'   \item{rows}{BED/FAM row indices used as the LD reference sample, or
#'     `NULL` if all BED individuals were used.}
#' }
#'
#' @export
make_sparse_ld <- function(Glist,
                          rows = NULL,
                          out_prefix = NULL,
                          chr = NULL,
                          cls = NULL,
                          pos_bp = NULL,
                          max_distance_bp = 0,
                          max_distance_variants = 1000,
                          r2_threshold = 0.001,
                          block_size = 1024,
                          nthreads = 1,
                          allow_full_ld = FALSE) {
  bedfiles <- as.character(Glist$bedfiles)
  chr <- if (is.null(chr)) which(!is.na(bedfiles) & nzchar(bedfiles)) else as.integer(chr)

  if (is.null(cls)) {
    cls <- lapply(chr, function(cc) {
      match(Glist$rsidsLD[[cc]], Glist$rsids[[cc]])
    })
  } else if (!is.list(cls)) {
    cls <- list(cls)
  }

  if (length(cls) != length(chr)) {
    stop("cls must have one element per chromosome/file in chr.")
  }

  cls <- lapply(cls, as.integer)
  names(cls) <- paste0("chr", chr)

  af <- Map(function(cc, cl) Glist$af[[cc]][cl], chr, cls)

  if (is.null(rows) && !is.null(Glist$idsLD)) {
    rows <- match(Glist$idsLD, Glist$ids)
  }
  if (!is.null(rows)) {
    rows <- as.integer(rows)
  }

  if (is.null(out_prefix)) {
    out_prefix <- file.path(dirname(bedfiles[chr[1]]), "sparseLD")
  }

  sparseLD_stream_CSR(
    bed_files = bedfiles[chr],
    n = Glist$n,
    cls = cls,
    out_prefix = out_prefix,
    rows = rows,
    af = af,
    pos_bp = pos_bp,
    max_distance_bp = max_distance_bp,
    max_distance_variants = max_distance_variants,
    r2_threshold = r2_threshold,
    block_size = block_size,
    nthreads = nthreads,
    allow_full_ld = allow_full_ld
  )

  Glist$sparseLD <- list(
    prefix = out_prefix,
    chr = chr,
    bed_files = bedfiles[chr],
    cls = cls,
    af = af,
    rows = rows,
    max_distance_bp = max_distance_bp,
    max_distance_variants = max_distance_variants,
    r2_threshold = r2_threshold,
    block_size = block_size
  )

  Glist
}
