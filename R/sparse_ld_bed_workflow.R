#' Compute BED Marker Sufficient Statistics
#'
#' Computes marker cross-products and marker-trait cross-products directly from
#' a PLINK BED file.
#'
#' @param bed_file Path to a PLINK BED file.
#' @param n Number of individuals in the BED file.
#' @param cls Marker column indices.
#' @param af Allele frequencies corresponding to `cls`.
#' @param y Phenotype matrix.
#' @param rows Optional 1-based individual row indices.
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
#' @param max_distance_variants Maximum marker-index distance between retained
#'   pairs.
#' @param r2_threshold Minimum squared-correlation threshold.
#' @param block_size Marker block size.
#' @param nthreads Number of OpenMP threads.
#' @return The output prefix and sparse-LD writing summary.
#' @name sparseLD_stream_CSR
#' @export
NULL

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

.resolve_pi_prior <- function(pi_marker = 0.001, pi_init = NULL,
                              pi_vb_init = NULL, pi_prior_mean = NULL,
                              pi_prior_strength = NULL, pi_prior_a = NULL,
                              pi_prior_b = NULL) {
 if (is.null(pi_init)) pi_init <- pi_marker
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
  pi_marker = pi_marker, pi_init = pi_init, pi_vb_init = pi_vb_init,
  pi_prior_mean = pi_prior_mean, pi_prior_strength = pi_prior_strength,
  pi_prior_a = pi_prior_a, pi_prior_b = pi_prior_b,
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
 has_vle_vld <- length(fit) >= 22
 names(fit)[seq_len(min(length(fit), length(nms)))] <-
  nms[seq_len(min(length(fit), length(nms)))]
 if (has_vle_vld) names(fit)[21:22] <- c("vle", "vld")

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
 if (sum(diag(out$covb)) > 0) out$rb <- cov2cor(out$covb)
 if (sum(diag(out$covg)) > 0) out$rg <- cov2cor(out$covg)
 if (sum(diag(out$cove)) > 0) out$re <- cov2cor(out$cove)
 out
}

#' Fit ST-BLR from BED Sufficient Statistics and Sparse LD
#'
#' Fits the single-trait ST-BLR sampler using sufficient statistics from
#' [bed_xtx_xty()] and a disk-backed CSR LD prefix from
#' [sparseLD_stream_CSR()].
#'
#' @param stats Sufficient statistics returned by [bed_xtx_xty()].
#' @param ld_prefix Prefix of the disk-backed CSR LD files.
#' @param n Sample size. Defaults to `stats$n` when available.
#' @param m Number of markers. Inferred from `stats` when omitted.
#' @param pi_marker Backward-compatible initial inclusion probability.
#' @param pi_init,pi_vb_init,pi_prior_mean,pi_prior_strength,pi_prior_a,pi_prior_b
#'   Inclusion-probability and marker-variance prior controls.
#' @param h2 Initial heritability.
#' @param nub,nue Prior degrees of freedom.
#' @param updateB,updateE,updatePi Logical sampler update controls.
#' @param adjE Residual adjustment factor.
#' @param nit,nburn,nthin MCMC iteration controls.
#' @param ncores Number of OpenMP threads.
#' @param seed Sampler seed.
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
#' @return A formatted ST-BLR fit.
#' @export
stblr_csr <- function(stats, ld_prefix, n = NULL, m = NULL,
                      pi_marker = 0.001, pi_init = NULL, pi_vb_init = NULL,
                      pi_prior_mean = NULL, pi_prior_strength = NULL,
                      pi_prior_a = NULL, pi_prior_b = NULL, h2 = 0.5,
                      nub = 4, nue = 4, updateB = TRUE, updateE = TRUE,
                      updatePi = TRUE, adjE = 0.9, nit = 1000, nburn = 100,
                      nthin = 1, ncores = 3, seed = 10, scheduled = FALSE,
                      full_sweep_every = 10, null_skip_base = 50,
                      null_skip_max = 200, candidate_threshold = 1e-3,
                      candidate_lifetime = 20, skip_nulls_burnin_only = FALSE,
                      wakeup_ld_neighbors = TRUE, wakeup_diff_threshold = 0,
                      wakeup_max_neighbors = 0, use_d_init = FALSE,
                      use_r_init = FALSE, rebuild_r_before_updateE = FALSE) {
 nt <- length(stats$yy)
 if (is.null(n)) n <- stats$n
 if (is.null(n)) stop("n must be supplied or available as stats$n.")
 if (is.null(m)) m <- if (!is.null(stats$m)) stats$m else length(stats$ww[[1]])

 arch <- .resolve_pi_prior(
  pi_marker, pi_init, pi_vb_init, pi_prior_mean, pi_prior_strength,
  pi_prior_a, pi_prior_b
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
   ncores = ncores, seed = seed
  )))
 }
 fit <- .format_stblr_fit(raw, nt, m, trait_names, variable_names)
 fit$input <- c(list(
  n = n, m = m, nt = nt, h2 = h2, nub = nub, nue = nue, vy = vy,
  B = pri$B, E = pri$E, ssb_prior = pri$ssb_prior,
  sse_prior = pri$sse_prior, updateB = updateB, updateE = updateE,
  updatePi = updatePi, adjE = adjE, nit = nit, nburn = nburn, nthin = nthin,
  ncores = ncores, seed = seed, scheduled = scheduled, ld_prefix = ld_prefix
 ), arch)
 fit
}

.make_bed_marker_data <- function(Glist, y, chr, cls, block_size, chains,
                                  rows = NULL) {
 chr <- as.integer(chr)
 if (length(chr) < 1 || anyNA(chr)) stop("chr must contain valid chromosome indices.")
 y <- as.matrix(y)
 if (!chains && length(chr) != 1) stop("Use chains = TRUE for multi-chromosome BED sampling.")

 if (is.null(rows)) {
  if (!is.null(rownames(y))) {
   if (is.null(Glist$ids)) {
    stop("rownames(y) are present but Glist$ids is missing; cannot match y rows to BED individuals.")
   }
   if (anyDuplicated(rownames(y))) stop("rownames(y) must not contain duplicates.")
   if (anyDuplicated(Glist$ids)) {
    stop("Glist$ids must not contain duplicates when matching rownames(y).")
   }
   rows <- match(rownames(y), Glist$ids)
   if (anyNA(rows)) {
    missing_ids <- rownames(y)[is.na(rows)]
    stop(
     "Some rownames(y) were not found in Glist$ids. First missing IDs: ",
     paste(head(missing_ids, 10), collapse = ", ")
    )
   }
  } else if (nrow(y) == Glist$n) {
   rows <- NULL
  } else {
   stop(
    "nrow(y) != Glist$n, but neither rows nor rownames(y) were supplied. ",
    "Provide rows as 1-based BED individual indices, or set rownames(y) ",
    "to IDs matching Glist$ids."
   )
  }
 } else {
  rows <- as.integer(rows)
  if (length(rows) != nrow(y)) stop("length(rows) must equal nrow(y).")
  if (anyNA(rows) || any(rows < 1L) || any(rows > Glist$n)) {
   stop("rows must be valid 1-based individual indices into the BED file.")
  }
  if (anyDuplicated(rows)) stop("rows must not contain duplicate individual indices.")
 }

 if (is.null(cls)) {
  cls <- lapply(chr, function(cc) match(Glist$rsidsLD[[cc]], Glist$rsids[[cc]]))
 } else if (!is.list(cls)) {
  cls <- list(cls)
 }
 if (length(cls) != length(chr) || any(vapply(cls, anyNA, logical(1)))) {
  stop("cls must match chr and contain no missing marker indices.")
 }
 cls <- lapply(cls, as.integer)
 af <- Map(function(cc, cl) Glist$af[[cc]][cl], chr, cls)
 m <- sum(lengths(cls))
 trait_names <- colnames(y)
 if (is.null(trait_names)) trait_names <- paste0("T", seq_len(ncol(y)))
 colnames(y) <- trait_names
 variable_names <- unlist(Map(function(cc, cl) Glist$rsids[[cc]][cl], chr, cls))
 sets <- if (length(chr) == 1) {
  rep(seq_len(ceiling(m / block_size)), each = block_size)[seq_len(m)]
 } else {
  rep(seq_along(chr), lengths(cls))
 }
 list(
  y = y, chr = chr, bed_files = Glist$bedfiles[chr], cls = cls, af = af,
  rows = rows, m = m, nt = ncol(y), n = Glist$n, n_total = Glist$n,
  n_used = nrow(y), trait_names = trait_names, variable_names = variable_names,
  sets = sets, b_init = lapply(seq_len(ncol(y)), function(i) rep(0, m))
 )
}
#' Fit ST-BLR Directly from PLINK BED Markers
#'
#' Fits ST-BLR directly from markers stored in PLINK BED files referenced by a
#' qgg genotype list.
#'
#' @param Glist A qgg genotype list.
#' @param y Phenotype matrix.
#' @param chr Chromosome indices to fit.
#' @param cls Optional marker column indices.
#' @param block_size Marker block size.
#' @param pi_marker Backward-compatible initial inclusion probability.
#' @param pi_init,pi_vb_init,pi_prior_mean,pi_prior_strength,pi_prior_a,pi_prior_b
#'   Inclusion-probability and marker-variance prior controls.
#' @param h2 Initial heritability.
#' @param nub,nue Prior degrees of freedom.
#' @param scale Standardize BED markers.
#' @param updateB,updateE,updatePi Logical sampler update controls.
#' @param adjE Residual adjustment factor.
#' @param nit,nburn,nthin MCMC iteration controls.
#' @param rebuild_every,full_sweep_every,candidate_threshold,candidate_lifetime
#'   Sparse sampler controls.
#' @param skip_nulls_burnin_only Restrict null-marker skipping to burn-in.
#' @param ncores Number of OpenMP threads.
#' @param seed Sampler seed.
#' @param scheduled Use the scheduled BED sampler.
#' @param null_update_prob,null_skip_base,null_skip_max Null-marker scheduling
#'   controls.
#' @param return_wy,return_r Return optional sampler state.
#' @param rows Optional 1-based individual row indices into the BED file. When
#'   omitted, subset phenotypes can instead be matched using `rownames(y)` and
#'   `Glist$ids`.
#' @param chains,nchains Multi-chain controls.
#' @param read_block_size,progress_every Multi-chain BED reading and progress
#'   controls.
#' @return A formatted ST-BLR fit.
#' @export
stblr_bed_marker <- function(
  Glist, y, chr = 1, cls = NULL, block_size = 1000, pi_marker = 0.001,
  pi_init = NULL, pi_vb_init = NULL, pi_prior_mean = NULL,
  pi_prior_strength = NULL, pi_prior_a = NULL, pi_prior_b = NULL,
  h2 = 0.5, nub = 4, nue = 4, scale = TRUE, updateB = TRUE,
  updateE = TRUE, updatePi = TRUE, adjE = 0, nit = 1000, nburn = 100,
  nthin = 1, rebuild_every = 25, full_sweep_every = 10,
  candidate_threshold = 1e-3, candidate_lifetime = 20,
  skip_nulls_burnin_only = FALSE, ncores = 3, seed = 10,
  scheduled = TRUE, null_update_prob = 0.02, null_skip_base = 50,
  null_skip_max = 200, return_wy = FALSE, return_r = FALSE, rows = NULL,
  chains = FALSE, nchains = 1, read_block_size = 64, progress_every = 0
) {
 arch <- .resolve_pi_prior(
  pi_marker, pi_init, pi_vb_init, pi_prior_mean, pi_prior_strength,
  pi_prior_a, pi_prior_b
 )
 dat <- .make_bed_marker_data(
   Glist = Glist,
   y = y,
   chr = chr,
   cls = cls,
   block_size = block_size,
   chains = chains,
   rows = rows
 )
 if (chains && !scheduled) {
  warning("chains=TRUE uses the scheduled multi-chain sampler; setting scheduled=TRUE.")
  scheduled <- TRUE
 }
 if (!chains && !scheduled && length(dat$chr) != 1) {
  stop("The sparse non-scheduled sampler path only supports one chromosome.")
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
 if (chains) {
  raw <- do.call(stblr_cpg_omp_bed_marker_scheduled_chains, c(common, list(
   null_skip_base = null_skip_base, null_skip_max = null_skip_max,
   return_wy = return_wy, return_r = return_r,
   read_block_size = as.integer(read_block_size),
   progress_every = as.integer(progress_every), pi_prior_a = arch$pi_prior_a,
   pi_prior_b = arch$pi_prior_b, nchains = as.integer(nchains),
   ncores = ncores, seed = seed
  )))
 } else if (scheduled) {
  raw <- do.call(stblr_cpg_omp_bed_marker_scheduled, c(common, list(
   null_skip_base = null_skip_base, null_skip_max = null_skip_max,
   return_wy = return_wy, return_r = return_r, pi_prior_a = arch$pi_prior_a,
   pi_prior_b = arch$pi_prior_b, ncores = ncores, seed = seed
  )))
 } else {
  raw <- do.call(stblr_cpg_omp_bed_marker_sparse, c(common, list(
   null_update_prob = null_update_prob, ncores = ncores, seed = seed
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
  scheduled = scheduled, chains = chains, nchains = nchains,
  read_block_size = read_block_size, progress_every = progress_every,
  scale = scale, rows = dat$rows
 ), arch)
 fit
}
