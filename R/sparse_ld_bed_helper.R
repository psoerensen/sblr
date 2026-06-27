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
 if (sum(diag(out$covb)) > 0) out$rb <- cov2cor(out$covb)
 if (sum(diag(out$covg)) > 0) out$rg <- cov2cor(out$covg)
 if (sum(diag(out$cove)) > 0) out$re <- cov2cor(out$cove)
 out
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
#' Fits the single-trait ST-BLR sampler using sufficient statistics from
#' [bed_xtx_xty()] and a disk-backed CSR LD prefix from
#' [sparseLD_stream_CSR()].
#'
#'
#' @param Glist Optional qgg genotype list containing `sparseLD`. If
#'   `ld_prefix` is `NULL`, `Glist$sparseLD$prefix` is used.
#' @param stats Sufficient statistics returned by [bed_xtx_xty()].
#' @param ld_prefix Prefix of the disk-backed CSR LD files. If `NULL`, this is
#'   taken from `Glist$sparseLD$prefix`.
#' @param n Sample size. Defaults to `stats$n` when available.
#' @param m Number of markers. Inferred from `stats` when omitted.
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
#' @param updateB,updateE,updatePi Logical sampler update controls.
#' @param adjE Residual adjustment factor. For sparse-LD summary-statistic
#'   models, values greater than 0 can stabilize residual-variance updates.
#' @param nit,nburn,nthin MCMC iteration controls.
#' @param ncores Number of OpenMP threads.
#' @param seed Sampler seed.
#' @param nchains Number of independent MCMC chains per trait. For
#'   `scheduled = FALSE`, chains are run as independent genome-wide
#'   trait-by-chain tasks in C++ and aggregated before return. `scheduled =
#'   TRUE` currently supports only `nchains = 1`.
#' @param keep_chains Logical; when `TRUE` and `scheduled = FALSE`, return
#'   compact per-chain `dm`, `bm`, and LD-swap diagnostics in `chains`.
#' @param chain_seeds Optional numeric or integer vector of length `nchains`.
#'   When supplied for `scheduled = FALSE`, each value is used as the base seed
#'   for that chain with a deterministic trait offset. Matrix seed inputs are
#'   not currently supported.
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
#' @return A formatted ST-BLR fit. For scheduled CSR fits, `pis` contains the
#'   full sampled inclusion-probability trace for each trait. Multi-chain
#'   regular CSR fits additionally provide `dm_sd`, `dm_min`, `dm_max`,
#'   `bm_sd`, `bm_min`, and `bm_max`; standard traces are averaged by iteration
#'   across chains.
#' @export
stblr_csr <- function(Glist=NULL, stats, ld_prefix=NULL, n = NULL, m = NULL,
                      pi_init = 0.001, pi_vb_init = NULL,
                      pi_prior_mean = 0.001, pi_prior_strength = 5e5,
                      pi_prior_a = NULL, pi_prior_b = NULL, h2 = 0.3,
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
                      ld_swap_moves = 1) {
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
 if (isTRUE(scheduled) && nchains != 1L) {
  stop("nchains > 1 is currently supported only for scheduled = FALSE.")
 }
 if (isTRUE(scheduled) && isTRUE(updateLDswap)) {
  stop("updateLDswap is currently implemented only for scheduled = FALSE.")
 }
 nt <- length(stats$yy)
 if (is.null(n)) n <- stats$n
 if (is.null(n)) stop("n must be supplied or available as stats$n.")
 if (is.null(m)) m <- if (!is.null(stats$m)) stats$m else length(stats$ww[[1]])

 if (is.null(ld_prefix)) {
   if (is.null(Glist) || is.null(Glist$sparseLD$prefix)) {
     stop("Provide ld_prefix or run make_sparseLD() and supply Glist.")
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
   pi_prior_b = arch$pi_prior_b, ncores = ncores, seed = seed
  )))
 } else {
 raw <- do.call(stblr_cpg_omp_csr, c(common, list(
   pi_prior_a = arch$pi_prior_a, pi_prior_b = arch$pi_prior_b,
   ncores = ncores, seed = seed, nchains = nchains,
   keep_chains = keep_chains, chain_seeds = chain_seeds,
   updateLDswap = updateLDswap,
   ld_swap_prob = ld_swap_prob, ld_swap_r2 = ld_swap_r2,
   ld_swap_max_friends = as.integer(ld_swap_max_friends),
   ld_swap_moves = as.integer(ld_swap_moves)
  )))
 }
 fit <- .format_stblr_fit(raw, nt, m, trait_names, variable_names)
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
  updateLDswap = updateLDswap, ld_swap_prob = ld_swap_prob,
  ld_swap_r2 = ld_swap_r2,
  ld_swap_max_friends = as.integer(ld_swap_max_friends),
  ld_swap_moves = as.integer(ld_swap_moves)
 ), arch)
 fit
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
 fit <- .format_stblr_fit(
  raw, dat$nt, dat$m, dat$trait_names, dat$variable_names, TRUE
 )
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
make_stats <- function(Glist, y, chr = NULL, cls = NULL, rows = NULL,
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
#' with [make_stats()] and [stblr_csr()].
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
make_sparseLD <- function(Glist,
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
