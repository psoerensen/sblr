# Development notebook containing candidate R wrappers, simulations, and
# benchmarks. It is not intended to run during package checks.
#
# Simplified R interface for individual-level BED STBLR marker sampler.
#
# This version reduces duplication by separating the wrapper into small helpers:
#   - resolve_pi_prior()
#   - make_stblr_priors()
#   - make_bed_marker_data()
#   - format_bed_marker_stblr_fit2()
#   - stblr_bed_marker()
#
# It supports:
#   - single-chain scheduled sampler with pi_prior_a/pi_prior_b
#   - multi-chain scheduled sampler
#   - legacy non-scheduled sparse sampler for one chromosome
#   - VLE/VLD outputs from updated C++ samplers


resolve_pi_prior <- function(
  pi_marker = 0.001,
  pi_init = NULL,
  pi_vb_init = NULL,
  pi_prior_mean = NULL,
  pi_prior_strength = NULL,
  pi_prior_a = NULL,
  pi_prior_b = NULL
) {
 if (is.null(pi_init)) pi_init <- pi_marker
 if (is.null(pi_vb_init)) pi_vb_init <- pi_init
 if (is.null(pi_prior_mean)) pi_prior_mean <- pi_init
 if (is.null(pi_prior_strength)) pi_prior_strength <- 2

 for (nm in c("pi_init", "pi_vb_init", "pi_prior_mean")) {
  val <- get(nm)
  if (!is.numeric(val) || length(val) != 1 || !is.finite(val) || val <= 0 || val >= 1) {
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
  pi_marker = pi_marker,
  pi_init = pi_init,
  pi_vb_init = pi_vb_init,
  pi_prior_mean = pi_prior_mean,
  pi_prior_strength = pi_prior_strength,
  pi_prior_a = pi_prior_a,
  pi_prior_b = pi_prior_b,
  pi = c(1 - pi_init, pi_init)
 )
}


make_stblr_priors <- function(
  y,
  m,
  h2 = 0.5,
  nub = 4,
  nue = 4,
  pi_vb_init = 0.001,
  pi_prior_mean = 0.001,
  trait_names = NULL
) {
 y <- as.matrix(y)
 nt <- ncol(y)
 n_used <- nrow(y)

 if (is.null(trait_names)) {
  trait_names <- colnames(y)
 }
 if (is.null(trait_names)) {
  trait_names <- paste0("T", seq_len(nt))
 }

 yy <- colSums(y^2)
 vy <- yy / (n_used - 1)

 B <- diag((vy * h2) / (m * pi_vb_init), nrow = nt, ncol = nt)
 E <- diag(vy * (1 - h2), nrow = nt, ncol = nt)

 ssb_prior <- diag(
  ((nub - 2) / nub) * (vy * h2) / (m * pi_prior_mean),
  nrow = nt,
  ncol = nt
 )

 sse_prior <- diag(
  ((nue - 2) / nue) * (vy * (1 - h2)),
  nrow = nt,
  ncol = nt
 )

 rownames(B) <- colnames(B) <- trait_names
 rownames(E) <- colnames(E) <- trait_names
 rownames(ssb_prior) <- colnames(ssb_prior) <- trait_names
 rownames(sse_prior) <- colnames(sse_prior) <- trait_names

 list(
  vy = vy,
  B = B,
  E = E,
  ssb_prior = ssb_prior,
  sse_prior = sse_prior,
  ssb_prior_list = split(ssb_prior, rep(seq_len(ncol(ssb_prior)), each = nrow(ssb_prior))),
  sse_prior_list = split(sse_prior, rep(seq_len(ncol(sse_prior)), each = nrow(sse_prior)))
 )
}


make_bed_marker_data <- function(
  Glist,
  y,
  chr = 1,
  cls = NULL,
  block_size = 1000,
  chains = FALSE
) {
 chr <- as.integer(chr)
 if (length(chr) < 1 || anyNA(chr)) {
  stop("chr must contain at least one valid chromosome index.")
 }

 y <- as.matrix(y)
 if (nrow(y) < 1 || ncol(y) < 1) {
  stop("y must be a matrix-like object with at least one row and one column.")
 }

 if (!chains && length(chr) != 1) {
  stop("Use chains = TRUE for multi-chromosome BED sampling.")
 }

 if (is.null(cls)) {
  cls_list <- lapply(chr, function(cc) {
   out <- match(Glist$rsidsLD[[cc]], Glist$rsids[[cc]])
   if (anyNA(out)) {
    stop("Missing rsids for chromosome ", cc, ": ", sum(is.na(out)))
   }
   out
  })
 } else if (is.list(cls)) {
  if (length(cls) != length(chr)) {
   stop("When cls is a list, length(cls) must equal length(chr).")
  }
  cls_list <- cls
 } else {
  if (length(chr) != 1) {
   stop("When chr has length > 1, cls must be NULL or a list with one element per chromosome.")
  }
  cls_list <- list(cls)
 }

 cls_list <- lapply(seq_along(cls_list), function(k) {
  x <- as.integer(cls_list[[k]])
  if (anyNA(x)) stop("cls contains NA for chromosome ", chr[k], ".")
  x
 })

 af_list <- Map(function(cc, cl) {
  out <- Glist$af[[cc]][cl]
  if (anyNA(out)) stop("af contains NA after subsetting for chromosome ", cc, ".")
  out
 }, chr, cls_list)

 bed_files <- Glist$bedfiles[chr]
 m <- sum(lengths(cls_list))
 nt <- ncol(y)

 trait_names <- colnames(y)
 if (is.null(trait_names)) {
  trait_names <- paste0("T", seq_len(nt))
  colnames(y) <- trait_names
 }

 variable_names <- unlist(
  Map(function(cc, cl) Glist$rsids[[cc]][cl], chr, cls_list),
  use.names = FALSE
 )
 if (is.null(variable_names) || anyNA(variable_names)) {
  variable_names <- paste0("V", seq_len(m))
 }

 if (length(chr) == 1) {
  sets <- rep(seq_len(ceiling(m / block_size)), each = block_size)[seq_len(m)]
 } else {
  sets <- rep(seq_along(chr), lengths(cls_list))
 }

 list(
  y = y,
  chr = chr,
  bed_files = bed_files,
  cls_list = cls_list,
  af_list = af_list,
  m = m,
  nt = nt,
  n = Glist$n,
  n_used = nrow(y),
  trait_names = trait_names,
  variable_names = variable_names,
  sets = sets,
  b_init = lapply(seq_len(nt), function(t) rep(0, m))
 )
}


format_bed_marker_stblr_fit2 <- function(
  fit,
  nt,
  m,
  trait_names = NULL,
  variable_names = NULL
) {
 nms20 <- c(
  "bm", "dm", "wy", "r", "b", "d", "o",
  "vbs", "vgs", "ves",
  "covb", "covg", "cove",
  "vb", "vg", "ve",
  "pi", "pim", "pitrait", "pimarker"
 )

 has_vle_vld <- length(fit) >= 22

 if (has_vle_vld) {
  names(fit)[seq_along(nms20)] <- nms20
  names(fit)[21:22] <- c("vle", "vld")
 } else {
  names(fit) <- nms20
 }

 if (is.null(trait_names)) trait_names <- paste0("T", seq_len(nt))
 if (is.null(variable_names)) variable_names <- paste0("V", seq_len(m))

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
  rownames(fit[[i]]) <- trait_names
  colnames(fit[[i]]) <- trait_names
 }

 fit[[17]] <- fit[[17]][[1]]
 fit[[18]] <- fit[[18]][[1]]

 ## Keep pitrait and pimarker.
 out <- fit[1:20]

 ## Format diagnostics slots.
 if (!is.null(out$pitrait)) {
  out$pitrait <- matrix(unlist(out$pitrait), ncol = 4, byrow = TRUE)
  rownames(out$pitrait) <- trait_names
  colnames(out$pitrait) <- c(
   "log_cpo",
   "mean_log_cpo",
   "seconds_mean",
   "seconds_max"
  )

  out$diagnostics <- out$pitrait
  out$log_cpo <- out$diagnostics[, "log_cpo"]
  out$mean_log_cpo <- out$diagnostics[, "mean_log_cpo"]
 }

 if (!is.null(out$pimarker)) {
  out$pimarker <- matrix(unlist(out$pimarker), ncol = 2, byrow = TRUE)
  rownames(out$pimarker) <- trait_names
  colnames(out$pimarker) <- c("nsamples", "n_used")
 }

 if (has_vle_vld) {
  out$vle <- as.matrix(as.data.frame(fit[[21]]))
  out$vld <- as.matrix(as.data.frame(fit[[22]]))

  rownames(out$vle) <- paste0("Iter", seq_len(nrow(out$vle)))
  rownames(out$vld) <- paste0("Iter", seq_len(nrow(out$vld)))

  colnames(out$vle) <- trait_names
  colnames(out$vld) <- trait_names
 }

 if (sum(diag(out$covb)) > 0) out$rb <- cov2cor(out$covb)
 if (sum(diag(out$covg)) > 0) out$rg <- cov2cor(out$covg)
 if (sum(diag(out$cove)) > 0) out$re <- cov2cor(out$cove)

 out
}

stblr_bed_marker <- function(
  Glist,
  y,
  chr = 1,
  cls = NULL,
  block_size = 1000,

  pi_marker = 0.001,
  pi_init = NULL,
  pi_vb_init = NULL,
  pi_prior_mean = NULL,
  pi_prior_strength = NULL,
  pi_prior_a = NULL,
  pi_prior_b = NULL,

  h2 = 0.5,
  nub = 4,
  nue = 4,
  scale = TRUE,
  updateB = TRUE,
  updateE = TRUE,
  updatePi = TRUE,
  adjE = 0,
  nit = 1000,
  nburn = 100,
  nthin = 1,
  rebuild_every = 25,
  full_sweep_every = 10,
  candidate_threshold = 1e-3,
  candidate_lifetime = 20,
  skip_nulls_burnin_only = FALSE,
  ncores = 3,
  seed = 10,
  scheduled = TRUE,
  null_update_prob = 0.02,
  null_skip_base = 50,
  null_skip_max = 200,
  return_wy = FALSE,
  return_r = FALSE,
  rows = NULL,

  chains = FALSE,
  nchains = 1,
  read_block_size = 64,
  progress_every = 0
) {
 arch <- resolve_pi_prior(
  pi_marker = pi_marker,
  pi_init = pi_init,
  pi_vb_init = pi_vb_init,
  pi_prior_mean = pi_prior_mean,
  pi_prior_strength = pi_prior_strength,
  pi_prior_a = pi_prior_a,
  pi_prior_b = pi_prior_b
 )

 dat <- make_bed_marker_data(
  Glist = Glist,
  y = y,
  chr = chr,
  cls = cls,
  block_size = block_size,
  chains = chains
 )

 pri <- make_stblr_priors(
  y = dat$y,
  m = dat$m,
  h2 = h2,
  nub = nub,
  nue = nue,
  pi_vb_init = arch$pi_vb_init,
  pi_prior_mean = arch$pi_prior_mean,
  trait_names = dat$trait_names
 )

 base_args <- list(
  bed_files = dat$bed_files,
  n = dat$n,
  cls = dat$cls_list,
  y = dat$y,
  b_init = dat$b_init,
  sets = dat$sets,
  rows = rows,
  af = dat$af_list,
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
  if (!scheduled) {
   warning("chains=TRUE uses the scheduled multi-chain sampler; setting scheduled=TRUE.")
   scheduled <- TRUE
  }

  raw_fit <- do.call(
   stblr_cpg_omp_bed_marker_scheduled_chains,
   c(
    base_args,
    list(
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
    )
   )
  )
 } else if (scheduled) {
  raw_fit <- do.call(
   stblr_cpg_omp_bed_marker_scheduled,
   c(
    base_args,
    list(
     null_skip_base = null_skip_base,
     null_skip_max = null_skip_max,
     return_wy = return_wy,
     return_r = return_r,
     pi_prior_a = arch$pi_prior_a,
     pi_prior_b = arch$pi_prior_b,
     ncores = ncores,
     seed = seed
    )
   )
  )
 } else {
  if (length(dat$chr) != 1) {
   stop("The sparse non-scheduled sampler path only supports one chromosome in this wrapper.")
  }

  raw_fit <- do.call(
   stblr_cpg_omp_bed_marker_sparse,
   c(
    base_args,
    list(
     null_update_prob = null_update_prob,
     ncores = ncores,
     seed = seed
    )
   )
  )
 }

 fit <- format_bed_marker_stblr_fit2(
  fit = raw_fit,
  nt = dat$nt,
  m = dat$m,
  trait_names = dat$trait_names,
  variable_names = dat$variable_names
 )

 fit$input <- c(
  list(
   chr = dat$chr,
   cls = dat$cls_list,
   n = dat$n,
   n_used = dat$n_used,
   m = dat$m,
   nt = dat$nt,
   block_size = block_size,
   n_blocks = length(unique(dat$sets)),
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
   nit = nit,
   nburn = nburn,
   nthin = nthin,
   rebuild_every = rebuild_every,
   full_sweep_every = full_sweep_every,
   candidate_threshold = candidate_threshold,
   candidate_lifetime = candidate_lifetime,
   skip_nulls_burnin_only = skip_nulls_burnin_only,
   null_update_prob = null_update_prob,
   null_skip_base = null_skip_base,
   null_skip_max = null_skip_max,
   ncores = ncores,
   seed = seed,
   scheduled = scheduled,
   chains = chains,
   nchains = nchains,
   read_block_size = read_block_size,
   progress_every = progress_every,
   scale = scale,
   rows = rows
  ),
  arch
 )

 fit
}


# format_stblr_fit <- function(
#   fit,
#   nt,
#   m,
#   trait_names = NULL,
#   variable_names = NULL
# ) {
#  names(fit) <- c(
#   "bm", "dm", "wy", "r", "b", "d", "o",
#   "vbs", "vgs", "ves",
#   "covb", "covg", "cove",
#   "vb", "vg", "ve",
#   "pi", "pim", "pitrait", "pimarker"
#  )
#
#  if (is.null(trait_names)) {
#   trait_names <- paste0("T", seq_len(nt))
#  }
#
#  if (is.null(variable_names)) {
#   variable_names <- paste0("V", seq_len(m))
#  }
#
#  for (i in 1:7) {
#   fit[[i]] <- as.matrix(as.data.frame(fit[[i]]))
#   rownames(fit[[i]]) <- variable_names
#   colnames(fit[[i]]) <- trait_names
#  }
#
#  for (i in 8:10) {
#   fit[[i]] <- as.matrix(as.data.frame(fit[[i]]))
#   rownames(fit[[i]]) <- paste0("Iter", seq_len(nrow(fit[[i]])))
#   colnames(fit[[i]]) <- trait_names
#  }
#
#  for (i in 11:16) {
#   fit[[i]] <- matrix(unlist(fit[[i]]), ncol = nt, byrow = TRUE)
#   rownames(fit[[i]]) <- trait_names
#   colnames(fit[[i]]) <- trait_names
#  }
#
#  fit[[17]] <- fit[[17]][[1]]
#  fit[[18]] <- fit[[18]][[1]]
#
#  fit <- fit[1:18]
#
#  if (sum(diag(fit$covb)) > 0) {
#   fit$rb <- cov2cor(fit$covb)
#  }
#
#  if (sum(diag(fit$covg)) > 0) {
#   fit$rg <- cov2cor(fit$covg)
#  }
#
#  if (sum(diag(fit$cove)) > 0) {
#   fit$re <- cov2cor(fit$cove)
#  }
#
#  fit
# }

format_stblr_fit <- function(fit, nt, m, trait_names = NULL, variable_names = NULL) {
 has_vle_vld <- length(fit) >= 22

 if (has_vle_vld) {
  names(fit)[1:22] <- c(
   "bm", "dm", "wy", "r", "b", "d", "o",
   "vbs", "vgs", "ves",
   "covb", "covg", "cove",
   "vb", "vg", "ve",
   "pi", "pim", "pitrait", "pimarker",
   "vle", "vld"
  )
 } else {
  names(fit) <- c(
   "bm", "dm", "wy", "r", "b", "d", "o",
   "vbs", "vgs", "ves",
   "covb", "covg", "cove",
   "vb", "vg", "ve",
   "pi", "pim", "pitrait", "pimarker"
  )
 }

 if (is.null(trait_names)) {
  trait_names <- paste0("T", seq_len(nt))
 }

 if (is.null(variable_names)) {
  variable_names <- paste0("V", seq_len(m))
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
  colnames(fit[[i]]) <- rownames(fit[[i]]) <- trait_names
 }

 fit[[17]] <- fit[[17]][[1]]
 fit[[18]] <- fit[[18]][[1]]

 out <- fit[1:18]

 if (has_vle_vld) {
  out$vle <- as.matrix(as.data.frame(fit[[21]]))
  out$vld <- as.matrix(as.data.frame(fit[[22]]))

  rownames(out$vle) <- paste0("Iter", seq_len(nrow(out$vle)))
  rownames(out$vld) <- paste0("Iter", seq_len(nrow(out$vld)))

  colnames(out$vle) <- trait_names
  colnames(out$vld) <- trait_names
 }

 if (sum(diag(out$covb)) > 0) out$rb <- cov2cor(out$covb)
 if (sum(diag(out$covg)) > 0) out$rg <- cov2cor(out$covg)
 if (sum(diag(out$cove)) > 0) out$re <- cov2cor(out$cove)

 out
}
stblr_csr <- function(
  stats,
  ld_prefix,
  n = NULL,
  m = NULL,

  # Backward-compatible architecture argument.
  # If the explicit arguments below are NULL, they inherit from pi_marker.
  pi_marker = 0.001,

  # Explicit architecture controls.
  # pi_init: initial inclusion probability used by the chain.
  # pi_vb_init: controls the initial marker-effect variance B.
  # pi_prior_mean: controls ssb_prior and, by default, the Beta prior mean for pi.
  # pi_prior_strength: a + b for the Beta prior on pi.
  pi_init = NULL,
  pi_vb_init = NULL,
  pi_prior_mean = NULL,
  pi_prior_strength = NULL,
  pi_prior_a = NULL,
  pi_prior_b = NULL,

  h2 = 0.5,
  nub = 4,
  nue = 4,
  updateB = TRUE,
  updateE = TRUE,
  updatePi = TRUE,
  adjE = 0.9,
  nit = 1000,
  nburn = 100,
  nthin = 1,
  ncores = 3,
  seed = 10,
  scheduled = FALSE,
  full_sweep_every = 10,
  null_skip_base = 50,
  null_skip_max = 200,
  candidate_threshold = 1e-3,
  candidate_lifetime = 20,
  skip_nulls_burnin_only = FALSE,
  wakeup_ld_neighbors = TRUE,
  wakeup_diff_threshold = 0.0,
  wakeup_max_neighbors = 0,
  use_d_init = FALSE,
  use_r_init = FALSE,
  rebuild_r_before_updateE = FALSE
) {
 nt <- length(stats$yy)

 if (is.null(n)) {
  if (!is.null(stats$n)) {
   n <- stats$n
  } else {
   stop("n must be supplied or available as stats$n.")
  }
 }

 if (is.null(m)) {
  if (!is.null(stats$m)) {
   m <- stats$m
  } else {
   m <- length(stats$ww[[1]])
  }
 }

 ## -------------------------------------------------------------------------
 ## Resolve architecture controls
 ## -------------------------------------------------------------------------
 if (is.null(pi_init)) {
  pi_init <- pi_marker
 }
 if (is.null(pi_vb_init)) {
  pi_vb_init <- pi_init
 }
 if (is.null(pi_prior_mean)) {
  pi_prior_mean <- pi_init
 }
 if (is.null(pi_prior_strength)) {
  # Weak default equivalent to Beta(1,1) if pi_prior_mean = 0.5,
  # but for sparse pi_prior_mean this becomes a very weak sparse prior.
  pi_prior_strength <- 2
 }

 for (nm in c("pi_init", "pi_vb_init", "pi_prior_mean")) {
  val <- get(nm)
  if (!is.numeric(val) || length(val) != 1 || !is.finite(val) || val <= 0 || val >= 1) {
   stop(nm, " must be a finite scalar in (0, 1).")
  }
 }

 if (!is.numeric(pi_prior_strength) || length(pi_prior_strength) != 1 ||
     !is.finite(pi_prior_strength) || pi_prior_strength <= 0) {
  stop("pi_prior_strength must be a finite positive scalar.")
 }

 if (is.null(pi_prior_a)) {
  pi_prior_a <- pi_prior_mean * pi_prior_strength
 }
 if (is.null(pi_prior_b)) {
  pi_prior_b <- (1 - pi_prior_mean) * pi_prior_strength
 }

 if (!is.numeric(pi_prior_a) || length(pi_prior_a) != 1 ||
     !is.finite(pi_prior_a) || pi_prior_a <= 0) {
  stop("pi_prior_a must be a finite positive scalar.")
 }
 if (!is.numeric(pi_prior_b) || length(pi_prior_b) != 1 ||
     !is.finite(pi_prior_b) || pi_prior_b <= 0) {
  stop("pi_prior_b must be a finite positive scalar.")
 }

 trait_names <- names(stats$yy)
 if (is.null(trait_names)) {
  trait_names <- paste0("T", seq_len(nt))
 }

 variable_names <- names(stats$ww[[1]])
 if (is.null(variable_names)) {
  variable_names <- paste0("V", seq_len(m))
 }

 b <- lapply(seq_len(nt), function(i) rep(0, m))
 d <- lapply(seq_len(nt), function(i) rep(0, m))

 vy <- as.numeric(stats$yy) / (n - 1)

 ## Initial marker-effect variance scale.
 B <- diag(
  (vy * h2) / (m * pi_vb_init),
  nrow = nt,
  ncol = nt
 )

 E <- diag(
  vy * (1 - h2),
  nrow = nt,
  ncol = nt
 )

 ## Long-run prior scale for marker-effect variance.
 ssb_prior <- diag(
  ((nub - 2) / nub) * (vy * h2) / (m * pi_prior_mean),
  nrow = nt,
  ncol = nt
 )

 sse_prior <- diag(
  ((nue - 2) / nue) * (vy * (1 - h2)),
  nrow = nt,
  ncol = nt
 )

 rownames(B) <- colnames(B) <- trait_names
 rownames(E) <- colnames(E) <- trait_names
 rownames(ssb_prior) <- colnames(ssb_prior) <- trait_names
 rownames(sse_prior) <- colnames(sse_prior) <- trait_names

 ssb_prior_list <- split(
  ssb_prior,
  rep(seq_len(ncol(ssb_prior)), each = nrow(ssb_prior))
 )

 sse_prior_list <- split(
  sse_prior,
  rep(seq_len(ncol(sse_prior)), each = nrow(sse_prior))
 )

 common_args <- list(
  wy = stats$wy,
  ww = stats$ww,
  yy = stats$yy,
  b_init = b,
  d_init = d,
  use_d_init = use_d_init,
  r_init = stats$wy,
  use_r_init = use_r_init,
  rebuild_r_before_updateE = rebuild_r_before_updateE,
  ld_prefix = ld_prefix,
  B = B,
  E = E,
  ssb_prior = ssb_prior_list,
  sse_prior = sse_prior_list,
  pi = c(1 - pi_init, pi_init),
  nub = nub,
  nue = nue,
  updateB = updateB,
  updateE = updateE,
  updatePi = updatePi,
  adjE = adjE,
  n = rep(as.integer(n), nt),
  nit = nit,
  nburn = nburn,
  nthin = nthin
 )

 if (scheduled) {
  ## Note: this assumes stblr_cpg_omp_csr_scheduled has also been updated to
  ## accept pi_prior_a and pi_prior_b. If not, remove those two arguments below
  ## or update the scheduled C++ sampler similarly.
  raw_fit <- do.call(
   stblr_cpg_omp_csr_scheduled,
   c(
    common_args,
    list(
     full_sweep_every = full_sweep_every,
     null_skip_base = null_skip_base,
     null_skip_max = null_skip_max,
     candidate_threshold = candidate_threshold,
     candidate_lifetime = candidate_lifetime,
     skip_nulls_burnin_only = skip_nulls_burnin_only,
     wakeup_ld_neighbors = wakeup_ld_neighbors,
     wakeup_diff_threshold = wakeup_diff_threshold,
     wakeup_max_neighbors = wakeup_max_neighbors,
     pi_prior_a = pi_prior_a,
     pi_prior_b = pi_prior_b,
     ncores = ncores,
     seed = seed
    )
   )
  )
 } else {
  raw_fit <- do.call(
   stblr_cpg_omp_csr,
   c(
    common_args,
    list(
     pi_prior_a = pi_prior_a,
     pi_prior_b = pi_prior_b,
     ncores = ncores,
     seed = seed
    )
   )
  )
 }

 fit <- format_stblr_fit(
  fit = raw_fit,
  nt = nt,
  m = m,
  trait_names = trait_names,
  variable_names = variable_names
 )

 fit$input <- list(
  n = n,
  m = m,
  nt = nt,

  # Backward-compatible and explicit architecture settings
  pi_marker = pi_marker,
  pi_init = pi_init,
  pi_vb_init = pi_vb_init,
  pi_prior_mean = pi_prior_mean,
  pi_prior_strength = pi_prior_strength,
  pi_prior_a = pi_prior_a,
  pi_prior_b = pi_prior_b,

  h2 = h2,
  nub = nub,
  nue = nue,
  vy = vy,
  B = B,
  E = E,
  ssb_prior = ssb_prior,
  sse_prior = sse_prior,
  updateB = updateB,
  updateE = updateE,
  updatePi = updatePi,
  adjE = adjE,
  nit = nit,
  nburn = nburn,
  nthin = nthin,
  ncores = ncores,
  seed = seed,
  scheduled = scheduled,
  full_sweep_every = full_sweep_every,
  null_skip_base = null_skip_base,
  null_skip_max = null_skip_max,
  candidate_threshold = candidate_threshold,
  candidate_lifetime = candidate_lifetime,
  skip_nulls_burnin_only = skip_nulls_burnin_only,
  wakeup_ld_neighbors = wakeup_ld_neighbors,
  wakeup_diff_threshold = wakeup_diff_threshold,
  wakeup_max_neighbors = wakeup_max_neighbors,
  use_d_init = use_d_init,
  use_r_init = use_r_init,
  rebuild_r_before_updateE = rebuild_r_before_updateE,
  ld_prefix = ld_prefix
 )

 fit
}

#
# prepare_bed_marker_inputs <- function(
#   Glist,
#   y,
#   chr = 1,
#   cls = NULL,
#   block_size = 1000,
#   pi_marker = 0.001,
#   h2 = 0.5,
#   nub = 4,
#   nue = 4,
#   scale = TRUE,
#   nthreads_stats = 4,
#   compute_stats = TRUE
# ) {
#  if (is.null(cls)) {
#   cls <- match(Glist$rsidsLD[[chr]], Glist$rsids[[chr]])
#
#   if (anyNA(cls)) {
#    stop("Some LD rsids were not found in Glist$rsids for this chromosome.")
#   }
#  }
#
#  m <- length(cls)
#  n <- Glist$n
#  nt <- ncol(y)
#
#  trait_names <- colnames(y)
#  if (is.null(trait_names)) {
#   trait_names <- paste0("T", seq_len(nt))
#  }
#
#  variable_names <- Glist$rsids[[chr]][cls]
#  if (is.null(variable_names)) {
#   variable_names <- paste0("V", seq_len(m))
#  }
#
#  sets <- rep(
#   seq_len(ceiling(m / block_size)),
#   each = block_size
#  )[seq_len(m)]
#
#  stats <- NULL
#
#  if (compute_stats) {
#   stats <- bed_xtx_xty(
#    bed_file = Glist$bedfiles[chr],
#    n = n,
#    cls = cls,
#    af = Glist$af[[chr]][cls],
#    y = y,
#    scale = scale,
#    nthreads = nthreads_stats
#   )
#
#   yy_diag <- if (is.matrix(stats$yy)) diag(stats$yy) else stats$yy
#  } else {
#   yy_diag <- colSums(scale(y, center = TRUE, scale = FALSE)^2)
#  }
#
#  vy <- yy_diag / (n - 1)
#
#  B <- diag((vy * h2) / (m * pi_marker), nrow = nt, ncol = nt)
#  E <- diag(vy * (1 - h2), nrow = nt, ncol = nt)
#
#  ssb_prior <- diag(
#   ((nub - 2) / nub) * (vy * h2) / (m * pi_marker),
#   nrow = nt,
#   ncol = nt
#  )
#
#  sse_prior <- diag(
#   ((nue - 2) / nue) * (vy * (1 - h2)),
#   nrow = nt,
#   ncol = nt
#  )
#
#  rownames(B) <- colnames(B) <- trait_names
#  rownames(E) <- colnames(E) <- trait_names
#  rownames(ssb_prior) <- colnames(ssb_prior) <- trait_names
#  rownames(sse_prior) <- colnames(sse_prior) <- trait_names
#
#  ssb_prior_list <- split(
#   ssb_prior,
#   rep(seq_len(ncol(ssb_prior)), each = nrow(ssb_prior))
#  )
#
#  sse_prior_list <- split(
#   sse_prior,
#   rep(seq_len(ncol(sse_prior)), each = nrow(sse_prior))
#  )
#
#  b_init <- lapply(seq_len(nt), function(t) rep(0, m))
#
#  list(
#   chr = chr,
#   cls = cls,
#   m = m,
#   n = n,
#   nt = nt,
#   y = y,
#   scale = scale,
#   block_size = block_size,
#   sets = sets,
#   stats = stats,
#   trait_names = trait_names,
#   variable_names = variable_names,
#   b_init = b_init,
#   pi = c(1 - pi_marker, pi_marker),
#   pi_marker = pi_marker,
#   h2 = h2,
#   nub = nub,
#   nue = nue,
#   vy = vy,
#   B = B,
#   E = E,
#   ssb_prior = ssb_prior,
#   sse_prior = sse_prior,
#   ssb_prior_list = ssb_prior_list,
#   sse_prior_list = sse_prior_list
#  )
# }
#
# format_bed_marker_stblr_fit <- function(
#   fit,
#   nt,
#   m,
#   trait_names = NULL,
#   variable_names = NULL
# ) {
#  names(fit) <- c(
#   "bm", "dm", "wy", "r", "b", "d", "o",
#   "vbs", "vgs", "ves",
#   "covb", "covg", "cove",
#   "vb", "vg", "ve",
#   "pi", "pim", "pitrait", "pimarker"
#  )
#
#  if (is.null(trait_names)) {
#   trait_names <- paste0("T", seq_len(nt))
#  }
#
#  if (is.null(variable_names)) {
#   variable_names <- paste0("V", seq_len(m))
#  }
#
#  for (i in 1:7) {
#   fit[[i]] <- as.matrix(as.data.frame(fit[[i]]))
#   rownames(fit[[i]]) <- variable_names
#   colnames(fit[[i]]) <- trait_names
#  }
#
#  for (i in 8:10) {
#   fit[[i]] <- as.matrix(as.data.frame(fit[[i]]))
#   rownames(fit[[i]]) <- paste0("Iter", seq_len(nrow(fit[[i]])))
#   colnames(fit[[i]]) <- trait_names
#  }
#
#  for (i in 11:16) {
#   fit[[i]] <- matrix(unlist(fit[[i]]), ncol = nt, byrow = TRUE)
#   rownames(fit[[i]]) <- trait_names
#   colnames(fit[[i]]) <- trait_names
#  }
#
#  fit[[17]] <- fit[[17]][[1]]
#  fit[[18]] <- fit[[18]][[1]]
#
#  fit <- fit[1:18]
#
#  if (sum(diag(fit$covb)) > 0) {
#   fit$rb <- cov2cor(fit$covb)
#  }
#
#  if (sum(diag(fit$covg)) > 0) {
#   fit$rg <- cov2cor(fit$covg)
#  }
#
#  if (sum(diag(fit$cove)) > 0) {
#   fit$re <- cov2cor(fit$cove)
#  }
#
#  fit
# }
#
# stblr_bed_marker <- function(
#   Glist,
#   y,
#   chr = 1,
#   cls = NULL,
#   block_size = 1000,
#
#   # Backward-compatible architecture argument.
#   # If the more explicit arguments below are NULL, they inherit from pi_marker.
#   pi_marker = 0.001,
#
#   # Explicit architecture controls.
#   # pi_init: initial inclusion probability used by the chain.
#   # pi_vb_init: controls the initial marker-effect variance B.
#   # pi_prior_mean: controls ssb_prior and, by default, the Beta prior mean for pi.
#   # pi_prior_strength: a + b for the Beta prior on pi.
#   pi_init = NULL,
#   pi_vb_init = NULL,
#   pi_prior_mean = NULL,
#   pi_prior_strength = NULL,
#   pi_prior_a = NULL,
#   pi_prior_b = NULL,
#
#   h2 = 0.5,
#   nub = 4,
#   nue = 4,
#   scale = TRUE,
#   compute_stats = TRUE,
#   nthreads_stats = 4,
#   updateB = TRUE,
#   updateE = TRUE,
#   updatePi = TRUE,
#   adjE = 0,
#   nit = 1000,
#   nburn = 100,
#   nthin = 1,
#   rebuild_every = 25,
#   full_sweep_every = 10,
#   candidate_threshold = 1e-3,
#   candidate_lifetime = 20,
#   skip_nulls_burnin_only = FALSE,
#   ncores = 3,
#   seed = 10,
#   scheduled = FALSE,
#   null_update_prob = 0.02,
#   null_skip_base = 50,
#   null_skip_max = 200,
#   return_wy = FALSE,
#   return_r = FALSE,
#   rows = NULL,
#
#   # Multi-chain packed BED scheduled sampler arguments.
#   chains = FALSE,
#   nchains = 1,
#   read_block_size = 64,
#   progress_every = 0
# ) {
#  ## -------------------------------------------------------------------------
#  ## Basic checks
#  ## -------------------------------------------------------------------------
#  chr <- as.integer(chr)
#  if (length(chr) < 1 || anyNA(chr)) {
#   stop("chr must contain at least one valid chromosome index.")
#  }
#
#  y <- as.matrix(y)
#  if (nrow(y) < 1 || ncol(y) < 1) {
#   stop("y must be a matrix-like object with at least one row and one column.")
#  }
#
#  if (!chains && length(chr) != 1) {
#   stop(
#    "The non-chain wrapper path currently expects length(chr) == 1. ",
#    "Use chains = TRUE for multi-chromosome BED sampling."
#   )
#  }
#
#  ## -------------------------------------------------------------------------
#  ## Validate and resolve architecture controls
#  ## -------------------------------------------------------------------------
#  if (is.null(pi_init)) {
#   pi_init <- pi_marker
#  }
#  if (is.null(pi_vb_init)) {
#   pi_vb_init <- pi_init
#  }
#  if (is.null(pi_prior_mean)) {
#   pi_prior_mean <- pi_init
#  }
#  if (is.null(pi_prior_strength)) {
#   # Beta(1,1)-like default when pi_prior_a/pi_prior_b are not supplied.
#   # This preserves old weak-prior behavior as closely as possible.
#   pi_prior_strength <- 2
#  }
#
#  for (nm in c("pi_init", "pi_vb_init", "pi_prior_mean")) {
#   val <- get(nm)
#   if (!is.numeric(val) || length(val) != 1 || !is.finite(val) || val <= 0 || val >= 1) {
#    stop(nm, " must be a finite scalar in (0, 1).")
#   }
#  }
#
#  if (!is.numeric(pi_prior_strength) || length(pi_prior_strength) != 1 ||
#      !is.finite(pi_prior_strength) || pi_prior_strength <= 0) {
#   stop("pi_prior_strength must be a finite positive scalar.")
#  }
#
#  if (is.null(pi_prior_a)) {
#   pi_prior_a <- pi_prior_mean * pi_prior_strength
#  }
#  if (is.null(pi_prior_b)) {
#   pi_prior_b <- (1 - pi_prior_mean) * pi_prior_strength
#  }
#
#  if (!is.numeric(pi_prior_a) || length(pi_prior_a) != 1 ||
#      !is.finite(pi_prior_a) || pi_prior_a <= 0) {
#   stop("pi_prior_a must be a finite positive scalar.")
#  }
#  if (!is.numeric(pi_prior_b) || length(pi_prior_b) != 1 ||
#      !is.finite(pi_prior_b) || pi_prior_b <= 0) {
#   stop("pi_prior_b must be a finite positive scalar.")
#  }
#
#  if (!is.numeric(nchains) || length(nchains) != 1 || !is.finite(nchains) || nchains < 1) {
#   stop("nchains must be a positive integer.")
#  }
#  nchains <- as.integer(nchains)
#
#  if (!is.numeric(read_block_size) || length(read_block_size) != 1 ||
#      !is.finite(read_block_size) || read_block_size < 1) {
#   stop("read_block_size must be a positive integer.")
#  }
#  read_block_size <- as.integer(read_block_size)
#
#  if (!is.numeric(progress_every) || length(progress_every) != 1 ||
#      !is.finite(progress_every) || progress_every < 0) {
#   stop("progress_every must be a non-negative integer.")
#  }
#  progress_every <- as.integer(progress_every)
#
#  ## -------------------------------------------------------------------------
#  ## Build chromosome/file-specific cls and af inputs.
#  ## For chains=TRUE, this supports chr = 1:22 directly.
#  ## For chains=FALSE, this is restricted to one chromosome above.
#  ## -------------------------------------------------------------------------
#  if (is.null(cls)) {
#   cls_list <- lapply(chr, function(cc) {
#    out <- match(Glist$rsidsLD[[cc]], Glist$rsids[[cc]])
#    if (anyNA(out)) {
#     stop("Missing rsids for chromosome ", cc, ": ", sum(is.na(out)))
#    }
#    out
#   })
#  } else if (is.list(cls)) {
#   if (length(cls) != length(chr)) {
#    stop("When cls is a list, length(cls) must equal length(chr).")
#   }
#   cls_list <- cls
#  } else {
#   if (length(chr) != 1) {
#    stop("When chr has length > 1, cls must be NULL or a list with one element per chromosome.")
#   }
#   cls_list <- list(cls)
#  }
#
#  cls_list <- lapply(seq_along(cls_list), function(k) {
#   x <- as.integer(cls_list[[k]])
#   if (anyNA(x)) stop("cls contains NA for chromosome index ", chr[k], ".")
#   x
#  })
#
#  af_list <- Map(function(cc, cl) {
#   out <- Glist$af[[cc]][cl]
#   if (anyNA(out)) {
#    stop("af contains NA after subsetting for chromosome ", cc, ".")
#   }
#   out
#  }, chr, cls_list)
#
#  bed_files <- Glist$bedfiles[chr]
#  m_fit <- sum(lengths(cls_list))
#  nt <- ncol(y)
#  n_used <- nrow(y)
#
#  trait_names <- colnames(y)
#  if (is.null(trait_names)) {
#   trait_names <- paste0("T", seq_len(nt))
#   colnames(y) <- trait_names
#  }
#
#  variable_names <- unlist(
#   Map(function(cc, cl) Glist$rsids[[cc]][cl], chr, cls_list),
#   use.names = FALSE
#  )
#  if (is.null(variable_names) || anyNA(variable_names)) {
#   variable_names <- paste0("V", seq_len(m_fit))
#  }
#
#  ## -------------------------------------------------------------------------
#  ## Prepare single-chromosome helper inputs when possible.
#  ## This preserves compatibility with the existing helper and formatter.
#  ## For multi-chromosome chain runs, we build the needed objects directly.
#  ## -------------------------------------------------------------------------
#  inp <- NULL
#  if (length(chr) == 1) {
#   inp <- prepare_bed_marker_inputs(
#    Glist = Glist,
#    y = y,
#    chr = chr,
#    cls = cls_list[[1]],
#    block_size = block_size,
#    pi_marker = pi_marker,
#    h2 = h2,
#    nub = nub,
#    nue = nue,
#    scale = scale,
#    nthreads_stats = nthreads_stats,
#    compute_stats = compute_stats
#   )
#
#   # Prefer helper-derived values for backward compatibility.
#   nt <- inp$nt
#   m_fit <- inp$m
#   trait_names <- inp$trait_names
#   variable_names <- inp$variable_names
#  }
#
#  ## -------------------------------------------------------------------------
#  ## Variance and prior setup.
#  ## This decouples chain initialization from the long-run architecture prior.
#  ## -------------------------------------------------------------------------
#  yy <- colSums(y^2)
#  vy <- yy / (n_used - 1)
#
#  B <- diag(
#   (vy * h2) / (m_fit * pi_vb_init),
#   nrow = nt,
#   ncol = nt
#  )
#
#  E <- diag(
#   vy * (1 - h2),
#   nrow = nt,
#   ncol = nt
#  )
#
#  ssb_prior <- diag(
#   ((nub - 2) / nub) * (vy * h2) / (m_fit * pi_prior_mean),
#   nrow = nt,
#   ncol = nt
#  )
#
#  sse_prior <- diag(
#   ((nue - 2) / nue) * (vy * (1 - h2)),
#   nrow = nt,
#   ncol = nt
#  )
#
#  rownames(B) <- colnames(B) <- trait_names
#  rownames(E) <- colnames(E) <- trait_names
#  rownames(ssb_prior) <- colnames(ssb_prior) <- trait_names
#  rownames(sse_prior) <- colnames(sse_prior) <- trait_names
#
#  ssb_prior_list <- split(
#   ssb_prior,
#   rep(seq_len(ncol(ssb_prior)), each = nrow(ssb_prior))
#  )
#
#  sse_prior_list <- split(
#   sse_prior,
#   rep(seq_len(ncol(sse_prior)), each = nrow(sse_prior))
#  )
#
#  pi <- c(1 - pi_init, pi_init)
#  b_init <- lapply(seq_len(nt), function(t) rep(0, m_fit))
#
#  ## -------------------------------------------------------------------------
#  ## Marker update sets
#  ## -------------------------------------------------------------------------
#  if (length(chr) == 1 && !is.null(inp)) {
#   sets <- inp$sets
#  } else {
#   # One update block per chromosome by default for multi-chromosome runs.
#   sets <- rep(seq_along(chr), lengths(cls_list))
#  }
#
#  ## -------------------------------------------------------------------------
#  ## Run sampler
#  ## -------------------------------------------------------------------------
#  if (chains) {
#   if (!scheduled) {
#    warning("chains=TRUE uses the scheduled multi-chain sampler; setting scheduled=TRUE.")
#    scheduled <- TRUE
#   }
#
#   raw_fit <- do.call(
#    stblr_cpg_omp_bed_marker_scheduled_chains,
#    list(
#     bed_files = bed_files,
#     n = Glist$n,
#     cls = cls_list,
#     y = y,
#     b_init = b_init,
#     sets = sets,
#     rows = rows,
#     af = af_list,
#     scale = scale,
#     B = B,
#     E = E,
#     ssb_prior = ssb_prior_list,
#     sse_prior = sse_prior_list,
#     pi = pi,
#     nub = nub,
#     nue = nue,
#     updateB = updateB,
#     updateE = updateE,
#     updatePi = updatePi,
#     adjE = adjE,
#     nit = nit,
#     nburn = nburn,
#     nthin = nthin,
#     rebuild_every = rebuild_every,
#     full_sweep_every = full_sweep_every,
#     null_skip_base = null_skip_base,
#     null_skip_max = null_skip_max,
#     candidate_threshold = candidate_threshold,
#     candidate_lifetime = candidate_lifetime,
#     skip_nulls_burnin_only = skip_nulls_burnin_only,
#     return_wy = return_wy,
#     return_r = return_r,
#     read_block_size = read_block_size,
#     progress_every = progress_every,
#     pi_prior_a = pi_prior_a,
#     pi_prior_b = pi_prior_b,
#     nchains = nchains,
#     ncores = ncores,
#     seed = seed
#    )
#   )
#  } else if (scheduled) {
#   # Updated single-chain scheduled C++ signature includes pi_prior_a/pi_prior_b.
#   raw_fit <- do.call(
#    stblr_cpg_omp_bed_marker_scheduled,
#    list(
#     bed_files = bed_files,
#     n = Glist$n,
#     cls = cls_list,
#     y = y,
#     b_init = b_init,
#     sets = sets,
#     rows = rows,
#     af = af_list,
#     scale = scale,
#     B = B,
#     E = E,
#     ssb_prior = ssb_prior_list,
#     sse_prior = sse_prior_list,
#     pi = pi,
#     nub = nub,
#     nue = nue,
#     updateB = updateB,
#     updateE = updateE,
#     updatePi = updatePi,
#     adjE = adjE,
#     nit = nit,
#     nburn = nburn,
#     nthin = nthin,
#     rebuild_every = rebuild_every,
#     full_sweep_every = full_sweep_every,
#     null_skip_base = null_skip_base,
#     null_skip_max = null_skip_max,
#     candidate_threshold = candidate_threshold,
#     candidate_lifetime = candidate_lifetime,
#     skip_nulls_burnin_only = skip_nulls_burnin_only,
#     return_wy = return_wy,
#     return_r = return_r,
#     pi_prior_a = pi_prior_a,
#     pi_prior_b = pi_prior_b,
#     ncores = ncores,
#     seed = seed
#    )
#   )
#  } else {
#   if (length(chr) != 1) {
#    stop("The sparse non-scheduled sampler path only supports one chromosome in this wrapper.")
#   }
#
#   # Keep the original sparse sampler call unchanged because it does not yet use
#   # pi_prior_a/pi_prior_b.
#   raw_fit <- do.call(
#    stblr_cpg_omp_bed_marker_sparse,
#    list(
#     bed_files = bed_files,
#     n = Glist$n,
#     cls = cls_list,
#     y = y,
#     b_init = b_init,
#     sets = sets,
#     rows = rows,
#     af = af_list,
#     scale = scale,
#     B = B,
#     E = E,
#     ssb_prior = ssb_prior_list,
#     sse_prior = sse_prior_list,
#     pi = pi,
#     nub = nub,
#     nue = nue,
#     updateB = updateB,
#     updateE = updateE,
#     updatePi = updatePi,
#     adjE = adjE,
#     nit = nit,
#     nburn = nburn,
#     nthin = nthin,
#     rebuild_every = rebuild_every,
#     full_sweep_every = full_sweep_every,
#     candidate_threshold = candidate_threshold,
#     candidate_lifetime = candidate_lifetime,
#     skip_nulls_burnin_only = skip_nulls_burnin_only,
#     ncores = ncores,
#     seed = seed,
#     null_update_prob = null_update_prob
#    )
#   )
#  }
#
#  ## -------------------------------------------------------------------------
#  ## Format output
#  ## -------------------------------------------------------------------------
#  fit <- format_bed_marker_stblr_fit(
#   fit = raw_fit,
#   nt = nt,
#   m = m_fit,
#   trait_names = trait_names,
#   variable_names = variable_names
#  )
#
#  # Updated scheduled samplers return:
#  # raw_fit[[21]] = VLE trace
#  # raw_fit[[22]] = VLD trace = VG - VLE
#  if (length(raw_fit) >= 22) {
#   fit$vle <- raw_fit[[21]]
#   fit$vld <- raw_fit[[22]]
#  }
#
#  fit$input <- list(
#   chr = chr,
#   cls = cls_list,
#   n = Glist$n,
#   n_used = n_used,
#   m = m_fit,
#   nt = nt,
#   block_size = block_size,
#   n_blocks = length(unique(sets)),
#   sets = sets,
#
#   # Backward-compatible and explicit architecture settings
#   pi_marker = pi_marker,
#   pi_init = pi_init,
#   pi_vb_init = pi_vb_init,
#   pi_prior_mean = pi_prior_mean,
#   pi_prior_strength = pi_prior_strength,
#   pi_prior_a = pi_prior_a,
#   pi_prior_b = pi_prior_b,
#
#   h2 = h2,
#   nub = nub,
#   nue = nue,
#   vy = vy,
#   B = B,
#   E = E,
#   ssb_prior = ssb_prior,
#   sse_prior = sse_prior,
#   updateB = updateB,
#   updateE = updateE,
#   updatePi = updatePi,
#   adjE = adjE,
#   nit = nit,
#   nburn = nburn,
#   nthin = nthin,
#   rebuild_every = rebuild_every,
#   full_sweep_every = full_sweep_every,
#   candidate_threshold = candidate_threshold,
#   candidate_lifetime = candidate_lifetime,
#   skip_nulls_burnin_only = skip_nulls_burnin_only,
#   null_update_prob = null_update_prob,
#   null_skip_base = null_skip_base,
#   null_skip_max = null_skip_max,
#   ncores = ncores,
#   seed = seed,
#   scheduled = scheduled,
#   chains = chains,
#   nchains = nchains,
#   read_block_size = read_block_size,
#   progress_every = progress_every,
#   scale = scale,
#   rows = rows
#  )
#
#  if (compute_stats && !is.null(inp) && !is.null(inp$stats)) {
#   fit$stats <- inp$stats
#  }
#
#  fit
# }



library(qgg)
library(sblr)

# data_dir <- "C:/Users/au223366/Documents/GitHub/examples/human"
# dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
#
# files <- c("bed", "bim", "fam", "pheno", "covar")
#
# for (f in files) {
#  url <- sprintf("https://github.com/psoerensen/qgdata/raw/main/simulated_human_data/human.%s",f)
#  download.file(url, destfile = file.path(data_dir, paste0("human.", f)), mode = "wb")
# }
#
# Glist <- gprep(
#  study = "Example",
#  bedfiles = "C:/Users/au223366/Documents/GitHub/examples/human/human.bed",
#  bimfiles = "C:/Users/au223366/Documents/GitHub/examples/human/human.bim",
#  famfiles = "C:/Users/au223366/Documents/GitHub/examples/human/human.fam"
# )
#
# rsids <- gfilter(Glist = Glist, excludeMAF = 0.05, excludeMISS = 0.05,
#                  excludeCGAT = TRUE, excludeINDEL = TRUE, excludeDUPS = TRUE, excludeHWE = 1e-12,
#                  excludeMHC = FALSE)
#
# ldfiles <- "C:/Users/au223366/Documents/GitHub/examples/human/human.ld"
# Glist <- gprep(Glist, task = "sparseld", msize = 1000, rsids = rsids, ldfiles = ldfiles,
#                overwrite = TRUE)
# saveRDS(Glist, file = "C:/Users/au223366/Documents/GitHub/examples/human/Glist_sparseLD_1k.RDS",
#         compress = FALSE)

Glist <- readRDS(file = "C:/Users/au223366/Documents/GitHub/examples/human/Glist_sparseLD_1k.RDS")

chr <- 1
rsids <- Glist$rsidsLD[[chr]]
h2 <- c(0.4, 0.5, 0.3)
rg <- matrix(
 c(
  1.0, 0.7, 0.3,
  0.7, 1.0, 0.5,
  0.3, 0.5, 1.0
 ),
 nrow = 3,
 byrow = TRUE
)

sim <- mtsim(
 Glist = Glist,
 chr=chr,
 rsids = rsids,
 nt = 3,
 n_shared = 30,
 n_specific = 10,
 h2 = h2,
 rg = rg,
 re = 0,
 seed = 1
)

stat <- glma(y = scale(sim$y[,1]), rsids=Glist$rsidsLD[[1]], Glist = Glist)
system.time(fitC <- gbayes(stat = stat, Glist = Glist, method = "bayesC", nit = 1000))


# Compute sumstats
y <- scale(sim$y)
cls <- match(Glist$rsidsLD[[chr]], Glist$rsids[[chr]])
chr <- 1
system.time(stats <- bed_xtx_xty(
 bed_file = Glist$bedfiles[chr],
 n = Glist$n,
 cls = cls,
 af = Glist$af[[chr]][cls],
 y = y,
 scale = TRUE,
 nthreads = 4
))

system.time(stats <- bed_xtx_xty(
 bed_file = Glist$bedfiles[chr],
 n = Glist$n,              # total number of individuals in BED
 cls = cls,
 af = NULL,                # computed over selected rows
 y = y,                    # must have nrow = length(rows)
 rows = NULL,        # 1-based individual rows in BED
 scale = TRUE,
 nthreads = 4
))


library(RhpcBLASctl)

blas_set_num_threads(1)
omp_set_num_threads(4)

system.time(out <- sparseLD_stream_CSR(
 bed_files = Glist$bedfiles[1],
 n = Glist$n,
 cls = list(cls),
 out_prefix = file.path(
  "C:/Users/au223366/Documents/GitHub/examples/human",
  "ld_test"
 ),
 rows = NULL,
 af = list(Glist$af[[1]][cls]),
 pos_bp = NULL,                 # optional; not needed if no bp window
 max_distance_bp = 0,            # no bp restriction
 max_distance_variants = 0,      # no variant-count restriction
 r2_threshold = 0.001,
 block_size = 1024,
 nthreads = 4
))

# ld <- sparseLD_read_CSR(
#  file.path("C:/Users/au223366/Documents/GitHub/examples/human", "ld_test"),
#  one_based = FALSE
# )
# summary(diff(ld$row_ptr))
# mean(diff(ld$row_ptr))
# max(diff(ld$row_ptr))
# out$nnz
# sprsrisk_cohl <- sparse_ld_risk(ld, transform = "coherence_log_abs")
# sprsrisk_coh <- sparse_ld_risk(ld, transform = "coherence")
# sprsrisk <- sparse_ld_risk(ld, transform = "signed_sum_abs")



ld_prefix <- file.path(
 "C:/Users/au223366/Documents/GitHub/examples/human",
 "ld_test"
)

system.time(fit_st <- stblr_csr(
 stats = stats,
 ld_prefix = ld_prefix,
 n = Glist$n,
 pi_marker = 0.001,

 pi_init = 0.001,
 pi_vb_init = 0.001,
 pi_prior_mean = 0.001,
 pi_prior_a = 1,
 pi_prior_b = 1,

 h2 = 0.5,
 adjE = 0.9,
 nit = 1000,
 nburn = 100,
 nthin = 1,
 ncores = 3,
 seed = 10,
 rebuild_r_before_updateE = FALSE,
 scheduled = FALSE
))


plot(sim$B[,1],fit_st$dm[,1])
plot(sim$B[,1],fitC[[1]]$dm)

fit_sched <- stblr_csr(
 stats = stats,
 ld_prefix = ld_prefix,
 n = Glist$n,
 pi_marker = 0.001,
 h2 = 0.5,
 adjE = 0.9,
 nit = 1000,
 nburn = 100,
 nthin = 1,
 ncores = 3,
 seed = 10,
 rebuild_r_before_updateE = TRUE,
 scheduled = TRUE,
 full_sweep_every = 10,
 null_skip_base = 50,
 null_skip_max = 200,
 candidate_threshold = 1e-3,
 candidate_lifetime = 20,
 skip_nulls_burnin_only = FALSE,
 wakeup_ld_neighbors = TRUE,
 wakeup_diff_threshold = 0.0,
 wakeup_max_neighbors = 0
)


fit_bed_sparse <- stblr_bed_marker(
 Glist = Glist,
 y = y,
 chr = 1,
 block_size = 1000,
 pi_marker = 0.001,
 h2 = 0.5,
 nit = 1000,
 nburn = 100,
 nthin = 1,
 rebuild_every = 25,
 full_sweep_every = 10,
 null_update_prob = 0.02,
 candidate_threshold = 1e-3,
 candidate_lifetime = 20,
 skip_nulls_burnin_only = FALSE,
 ncores = 3,
 seed = 10,
 scheduled = FALSE
)

fit_bed_sched <- stblr_bed_marker(
 Glist = Glist,
 y = y,
 chr = 1,
 block_size = 1000,

 # Old-style sparse architecture
 pi_marker = 0.001,

 # New explicit architecture controls
 pi_init = 0.001,
 pi_vb_init = 0.001,
 pi_prior_mean = 0.001,
 pi_prior_strength = 1e4,
 pi_prior_a = 1,
 pi_prior_b = 1,

 h2 = 0.5,
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
 ncores = 3,
 seed = 10,
 scheduled = TRUE,
 chains = FALSE
)

fit_bed_sched <- stblr_bed_marker(
 Glist = Glist,
 y = y,
 chr = 1,
 block_size = 1000,

 # Old-style sparse architecture
 pi_marker = 0.001,

 # New explicit architecture controls
 pi_init = 0.001,
 pi_vb_init = 0.001,
 pi_prior_mean = 0.001,
 pi_prior_strength = 1e4,
 pi_prior_a = 1,
 pi_prior_b = 1,

 h2 = 0.5,
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

 ncores = 9,
 seed = 10,
 scheduled = TRUE,
 chains = TRUE,
 nchains = 3,

 read_block_size = 64,
 progress_every = 200
)

# fit_bed_sched <- stblr_bed_marker(
#  Glist = Glist,
#  y = y,
#  chr = 1,
#  block_size = 1000,
#  pi_marker = 0.001,
#  h2 = 0.5,
#  nit = 1000,
#  nburn = 100,
#  nthin = 1,
#  rebuild_every = 25,
#  full_sweep_every = 10,
#  null_skip_base = 50,
#  null_skip_max = 200,
#  candidate_threshold = 1e-3,
#  candidate_lifetime = 20,
#  skip_nulls_burnin_only = FALSE,
#  return_wy = FALSE,
#  return_r = FALSE,
#  ncores = 3,
#  seed = 10,
#  scheduled = TRUE
# )

cor(fit_st$dm[, 1], fit_bed_sched$dm[, 1])
tail(sort(fit_st$dm[, 1] - fit_bed_sched$dm[, 1]))

rec_fit_st <- recovery_summary(fit_st, sim)
rec_fit_sched <- recovery_summary(fit_sched, sim)
rec_fit_bed_sparse <- recovery_summary(fit_bed_sparse, sim)
rec_fit_bed_sched <- recovery_summary(fit_bed_sched, sim)

library(dplyr)
library(ggplot2)
library(tidyr)

rec_list <- list(
 ST = rec_fit_st,
 ST_scheduled = rec_fit_sched,
 BED_sparse = rec_fit_bed_sparse,
 BED_scheduled = rec_fit_bed_sched
)

topk_all <- bind_rows(lapply(names(rec_list), function(x) {
 rec_list[[x]]$topk |>
  mutate(method = x)
}))

global_all <- bind_rows(lapply(names(rec_list), function(x) {
 rec_list[[x]]$global |>
  mutate(method = x)
}))

threshold_all <- bind_rows(lapply(names(rec_list), function(x) {
 rec_list[[x]]$thresholds |>
  mutate(method = x)
}))

ggplot(topk_all, aes(x = K, y = recall_dm, color = method, group = method)) +
 geom_line() +
 geom_point() +
 facet_wrap(~ trait_name) +
 labs(
  x = "Top K markers",
  y = "Recall among true causal variants",
  title = "Top-K causal variant recovery"
 ) +
 theme_bw()

ggplot(topk_all, aes(x = K, y = precision_dm, color = method, group = method)) +
 geom_line() +
 geom_point() +
 facet_wrap(~ trait_name) +
 labs(
  x = "Top K markers",
  y = "Precision",
  title = "Top-K precision"
 ) +
 theme_bw()

global_long <- global_all |>
 select(method, trait_name, auc_dm, ap_dm) |>
 pivot_longer(
  cols = c(auc_dm, ap_dm),
  names_to = "metric",
  values_to = "value"
 )

ggplot(global_long, aes(x = method, y = value, fill = method)) +
 geom_col() +
 facet_grid(metric ~ trait_name) +
 labs(
  x = NULL,
  y = "Value",
  title = "Global ranking performance"
 ) +
 theme_bw() +
 theme(axis.text.x = element_text(angle = 45, hjust = 1),
       legend.position = "none")


ggplot(threshold_all, aes(x = recall, y = precision, color = method)) +
 geom_path(aes(group = method)) +
 geom_point(aes(shape = factor(threshold)), size = 2.5) +
 facet_wrap(~ trait_name) +
 labs(
  x = "Recall",
  y = "Precision",
  shape = "PIP threshold",
  title = "PIP threshold precision-recall tradeoff"
 ) +
 theme_bw()

ggplot(threshold_all, aes(x = factor(threshold), y = n_selected, fill = method)) +
 geom_col(position = position_dodge()) +
 facet_wrap(~ trait_name) +
 labs(
  x = "PIP threshold",
  y = "Number of selected markers",
  title = "Number of selected markers by PIP threshold"
 ) +
 theme_bw()

ggplot(global_all, aes(x = method, y = pip_sum, fill = method)) +
 geom_col() +
 geom_hline(
  aes(yintercept = n_causal),
  linetype = 2
 ) +
 facet_wrap(~ trait_name) +
 labs(
  x = NULL,
  y = "Sum of PIPs",
  title = "Posterior expected number of active markers"
 ) +
 theme_bw() +
 theme(axis.text.x = element_text(angle = 45, hjust = 1),
       legend.position = "none")



rank_long <- global_all |>
 select(method, trait_name, median_rank_dm, worst_rank_dm) |>
 pivot_longer(
  cols = c(median_rank_dm, worst_rank_dm),
  names_to = "rank_type",
  values_to = "rank"
 )

ggplot(rank_long, aes(x = method, y = rank, fill = method)) +
 geom_col() +
 facet_grid(rank_type ~ trait_name, scales = "free_y") +
 labs(
  x = NULL,
  y = "Rank of causal variants",
  title = "Median and worst causal ranks"
 ) +
 theme_bw() +
 theme(axis.text.x = element_text(angle = 45, hjust = 1),
       legend.position = "none")


false_discovery_summary <- function(rec, threshold = 0.01) {
 x <- rec$thresholds[rec$thresholds$threshold == threshold, ]
 data.frame(
  trait = x$trait_name,
  threshold = threshold,
  selected = x$n_selected,
  recovered = x$recovered,
  false_selected = x$n_selected - x$recovered,
  false_per_true = (x$n_selected - x$recovered) / pmax(x$recovered, 1),
  precision = x$precision,
  recall = x$recall
 )
}

lapply(
 list(
  ST = rec_fit_st,
  ST_scheduled = rec_fit_sched,
  BED_sparse = rec_fit_bed_sparse,
  BED_scheduled = rec_fit_bed_sched
 ),
 false_discovery_summary,
 threshold = 0.01
)


fd_all <- do.call(rbind, lapply(names(rec_list), function(method) {
 x <- false_discovery_summary(rec_list[[method]], threshold = 0.01)
 cbind(method = method, x)
}))

ggplot(fd_all, aes(x = method, y = false_per_true, fill = method)) +
 geom_col() +
 facet_wrap(~ trait) +
 labs(
  x = NULL,
  y = "False selected per recovered causal",
  title = "Low-PIP false discovery burden at PIP >= 0.01"
 ) +
 theme_bw() +
 theme(
  axis.text.x = element_text(angle = 45, hjust = 1),
  legend.position = "none"
 )


threshold_all <- bind_rows(lapply(names(rec_list), function(method) {
 rec_list[[method]]$thresholds |>
  mutate(method = method)
}))

ggplot(threshold_all, aes(x = threshold, y = n_selected, color = method)) +
 geom_line() +
 geom_point() +
 scale_x_log10() +
 facet_wrap(~ trait_name) +
 labs(
  x = "PIP threshold",
  y = "Number selected",
  title = "Selected markers across PIP thresholds"
 ) +
 theme_bw()


# high-confidence calls
rec_fit_st$thresholds[rec_fit_st$thresholds$threshold == 0.5, ]

# top-K recovery around the expected causal count
rec_fit_st$topk[rec_fit_st$topk$K == 40, ]


credible_set_summary <- function(fit, sim, coverage = 0.95) {
 dm <- as.matrix(fit$dm)
 Btrue <- as.matrix(sim$B)

 nt <- ncol(Btrue)
 trait_names <- colnames(Btrue)
 if (is.null(trait_names)) trait_names <- paste0("T", seq_len(nt))

 out <- lapply(seq_len(nt), function(t) {
  pip <- dm[, t]
  causal <- Btrue[, t] != 0

  ord <- order(pip, decreasing = TRUE)
  cum_pip <- cumsum(pip[ord])

  k_cs <- which(cum_pip >= coverage)[1]
  cs <- ord[seq_len(k_cs)]

  data.frame(
   trait = t,
   trait_name = trait_names[t],
   coverage = coverage,
   cs_size = length(cs),
   recovered = sum(causal[cs]),
   n_causal = sum(causal),
   precision = sum(causal[cs]) / length(cs),
   recall = sum(causal[cs]) / sum(causal),
   pip_sum = sum(pip)
  )
 })

 do.call(rbind, out)
}

credible_set_summary(fit_st, sim, coverage = 0.95)
credible_set_summary(fit_bed_sparse, sim, coverage = 0.95)


credible_mass_summary <- function(fit, Btrue, coverage = 0.95) {
 if (is.list(Btrue) && !is.null(Btrue$B)) {
  Btrue <- Btrue$B
 }

 dm <- as.matrix(fit$dm)
 Btrue <- as.matrix(Btrue)

 nt <- ncol(Btrue)
 trait_names <- colnames(Btrue)
 if (is.null(trait_names)) trait_names <- paste0("T", seq_len(nt))

 out <- lapply(seq_len(nt), function(t) {
  pip <- dm[, t]
  causal <- Btrue[, t] != 0

  ord <- order(pip, decreasing = TRUE)
  target <- coverage * sum(pip)
  cum_pip <- cumsum(pip[ord])

  k_cs <- which(cum_pip >= target)[1]
  cs <- ord[seq_len(k_cs)]

  data.frame(
   trait = t,
   trait_name = trait_names[t],
   coverage = coverage,
   target_pip_mass = target,
   pip_sum = sum(pip),
   cs_size = length(cs),
   recovered = sum(causal[cs]),
   n_causal = sum(causal),
   precision = sum(causal[cs]) / length(cs),
   recall = sum(causal[cs]) / sum(causal),
   mean_pip_in_set = mean(pip[cs]),
   min_pip_in_set = min(pip[cs])
  )
 })

 do.call(rbind, out)
}


credible_mass_summary(fit_st, sim, coverage = 0.95)
credible_mass_summary(fit_sched, sim, coverage = 0.95)
credible_mass_summary(fit_bed_sparse, sim, coverage = 0.95)
credible_mass_summary(fit_bed_sched, sim, coverage = 0.95)

cs_compare <- rbind(
 cbind(method = "ST", credible_mass_summary(fit_st, sim, 0.95)),
 cbind(method = "ST_scheduled", credible_mass_summary(fit_sched, sim, 0.95)),
 cbind(method = "BED_sparse", credible_mass_summary(fit_bed_sparse, sim, 0.95)),
 cbind(method = "BED_scheduled", credible_mass_summary(fit_bed_sched, sim, 0.95))
)

cs_compare[, c(
 "method", "trait_name", "pip_sum", "cs_size",
 "recovered", "precision", "recall", "min_pip_in_set"
)]

library(ggplot2)

ggplot(cs_compare, aes(x = method, y = cs_size, fill = method)) +
 geom_col() +
 facet_wrap(~ trait_name, scales = "free_y") +
 labs(
  x = NULL,
  y = "Markers needed for 95% total PIP mass",
  title = "Posterior mass concentration differs strongly between CSR and BED"
 ) +
 theme_bw() +
 theme(
  axis.text.x = element_text(angle = 45, hjust = 1),
  legend.position = "none"
 )


pip_floor_summary <- function(fit, thresholds = c(0.001, 0.002, 0.005, 0.01, 0.05, 0.1, 0.5)) {
 dm <- as.matrix(fit$dm)

 do.call(rbind, lapply(seq_len(ncol(dm)), function(t) {
  pip <- dm[, t]
  data.frame(
   trait = colnames(dm)[t],
   threshold = thresholds,
   n_above = sapply(thresholds, function(th) sum(pip >= th)),
   mass_above = sapply(thresholds, function(th) sum(pip[pip >= th])),
   total_mass = sum(pip),
   frac_mass_above = sapply(thresholds, function(th) sum(pip[pip >= th]) / sum(pip))
  )
 }))
}

pip_floor_summary(fit_st)
pip_floor_summary(fit_bed_sparse)

credible_mass_summary_trimmed <- function(fit, sim, coverage = 0.95, min_pip = 0.005) {
 Btrue <- if (is.list(sim) && !is.null(sim$B)) sim$B else sim
 dm <- as.matrix(fit$dm)
 Btrue <- as.matrix(Btrue)

 nt <- ncol(Btrue)
 trait_names <- colnames(Btrue)
 if (is.null(trait_names)) trait_names <- paste0("T", seq_len(nt))

 out <- lapply(seq_len(nt), function(t) {
  pip <- dm[, t]
  causal <- Btrue[, t] != 0

  keep <- pip >= min_pip
  pip2 <- pip
  pip2[!keep] <- 0

  ord <- order(pip2, decreasing = TRUE)
  target <- coverage * sum(pip2)
  cum_pip <- cumsum(pip2[ord])

  k_cs <- if (target > 0) which(cum_pip >= target)[1] else 0
  cs <- if (k_cs > 0) ord[seq_len(k_cs)] else integer(0)

  data.frame(
   trait = t,
   trait_name = trait_names[t],
   coverage = coverage,
   min_pip = min_pip,
   pip_sum_raw = sum(pip),
   pip_sum_trimmed = sum(pip2),
   cs_size = length(cs),
   recovered = sum(causal[cs]),
   n_causal = sum(causal),
   precision = if (length(cs) > 0) sum(causal[cs]) / length(cs) else NA_real_,
   recall = sum(causal[cs]) / sum(causal),
   min_pip_in_set = if (length(cs) > 0) min(pip[cs]) else NA_real_
  )
 })

 do.call(rbind, out)
}
credible_mass_summary_trimmed(fit_st, sim, 0.95, min_pip = 0.005)
credible_mass_summary_trimmed(fit_st, sim, 0.95, min_pip = 0.01)
credible_mass_summary_trimmed(fit_bed_sparse, sim, 0.95, min_pip = 0.005)

background_pip_summary <- function(fit, thresholds = c(0.002, 0.005, 0.01, 0.05)) {
 dm <- as.matrix(fit$dm)

 do.call(rbind, lapply(seq_len(ncol(dm)), function(t) {
  pip <- dm[, t]
  trait <- colnames(dm)[t]
  if (is.null(trait)) trait <- paste0("T", t)

  do.call(rbind, lapply(thresholds, function(th) {
   data.frame(
    trait = trait,
    threshold = th,
    total_mass = sum(pip),
    n_below = sum(pip < th),
    mass_below = sum(pip[pip < th]),
    frac_mass_below = sum(pip[pip < th]) / sum(pip),
    n_above = sum(pip >= th),
    mass_above = sum(pip[pip >= th]),
    frac_mass_above = sum(pip[pip >= th]) / sum(pip)
   )
  }))
 }))
}
background_pip_summary(fit_st, 0.01)
background_pip_summary(fit_bed_sparse, 0.01)

credible_mass_summary_trimmed(fit_st, sim, coverage = 0.95, min_pip = 0.01)
credible_mass_summary_trimmed(fit_bed_sparse, sim, coverage = 0.95, min_pip = 0.01)

beta_estimation_summary <- function(fit, Btrue) {
 if (is.list(Btrue) && !is.null(Btrue$B)) {
  Btrue <- Btrue$B
 }

 bm <- as.matrix(fit$bm)
 Btrue <- as.matrix(Btrue)

 if (!all(dim(bm) == dim(Btrue))) {
  stop("dim(fit$bm) must match dim(Btrue).")
 }

 nt <- ncol(Btrue)
 trait_names <- colnames(Btrue)
 if (is.null(trait_names)) trait_names <- paste0("T", seq_len(nt))

 out <- lapply(seq_len(nt), function(t) {
  causal <- Btrue[, t] != 0

  err_all <- bm[, t] - Btrue[, t]
  err_causal <- bm[causal, t] - Btrue[causal, t]

  data.frame(
   trait = t,
   trait_name = trait_names[t],
   cor_all = cor(bm[, t], Btrue[, t]),
   cor_causal = cor(bm[causal, t], Btrue[causal, t]),
   rmse_all = sqrt(mean(err_all^2)),
   rmse_causal = sqrt(mean(err_causal^2)),
   mae_all = mean(abs(err_all)),
   mae_causal = mean(abs(err_causal)),
   bias_all = mean(err_all),
   bias_causal = mean(err_causal),
   mean_abs_beta_hat_causal = mean(abs(bm[causal, t])),
   mean_abs_beta_true_causal = mean(abs(Btrue[causal, t])),
   shrinkage_ratio = mean(abs(bm[causal, t])) / mean(abs(Btrue[causal, t]))
  )
 })

 do.call(rbind, out)
}

fit_list <- list(
 ST = fit_st,
 ST_scheduled = fit_sched,
 BED_sparse = fit_bed_sparse,
 BED_scheduled = fit_bed_sched
)

beta_cmp <- do.call(rbind, lapply(names(fit_list), function(method) {
 x <- beta_estimation_summary(fit_list[[method]], sim)
 x$method <- method
 x
}))

beta_cmp[, c(
 "method", "trait_name",
 "cor_all", "cor_causal",
 "rmse_all", "rmse_causal",
 "shrinkage_ratio"
)]

plot_beta_true_vs_est <- function(fit_list, Btrue) {
 if (is.list(Btrue) && !is.null(Btrue$B)) {
  Btrue <- Btrue$B
 }

 out <- do.call(rbind, lapply(names(fit_list), function(method) {
  fit <- fit_list[[method]]
  bm <- as.matrix(fit$bm)

  do.call(rbind, lapply(seq_len(ncol(Btrue)), function(t) {
   causal <- Btrue[, t] != 0

   data.frame(
    method = method,
    trait = colnames(Btrue)[t],
    beta_true = Btrue[causal, t],
    beta_hat = bm[causal, t]
   )
  }))
 }))

 ggplot(out, aes(x = beta_true, y = beta_hat, color = method)) +
  geom_point(alpha = 0.7) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  facet_wrap(~ trait, scales = "free") +
  labs(
   x = "True causal effect",
   y = "Posterior mean effect",
   title = "Estimated vs true effects at causal markers"
  ) +
  theme_bw()
}

plot_beta_true_vs_est(fit_list, sim)


dm_estimation_summary <- function(
  fit,
  Btrue,
  pip_thresholds = c(0.01, 0.05, 0.1, 0.5)
) {
 if (is.list(Btrue) && !is.null(Btrue$B)) {
  Btrue <- Btrue$B
 }

 dm <- as.matrix(fit$dm)
 Btrue <- as.matrix(Btrue)

 if (!all(dim(dm) == dim(Btrue))) {
  stop("dim(fit$dm) must match dim(Btrue).")
 }

 auc_binary <- function(score, label) {
  label <- as.logical(label)
  n_pos <- sum(label)
  n_neg <- sum(!label)

  if (n_pos == 0 || n_neg == 0) return(NA_real_)

  r <- rank(score, ties.method = "average")
  (sum(r[label]) - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
 }

 average_precision <- function(score, label) {
  label <- as.logical(label)
  n_pos <- sum(label)

  if (n_pos == 0) return(NA_real_)

  ord <- order(score, decreasing = TRUE)
  y <- label[ord]

  mean(cumsum(y)[y] / which(y))
 }

 safe_log_loss <- function(p, y, eps = 1e-12) {
  p <- pmin(pmax(p, eps), 1 - eps)
  -mean(y * log(p) + (1 - y) * log(1 - p))
 }

 nt <- ncol(Btrue)
 trait_names <- colnames(Btrue)
 if (is.null(trait_names)) trait_names <- paste0("T", seq_len(nt))

 global <- lapply(seq_len(nt), function(t) {
  pip <- dm[, t]
  causal <- Btrue[, t] != 0
  y <- as.numeric(causal)

  mean_pip_causal <- mean(pip[causal], na.rm = TRUE)
  mean_pip_noncausal <- mean(pip[!causal], na.rm = TRUE)

  data.frame(
   trait = t,
   trait_name = trait_names[t],
   n_markers = length(pip),
   n_causal = sum(causal),
   pip_sum = sum(pip, na.rm = TRUE),
   pip_sum_over_n_causal = sum(pip, na.rm = TRUE) / sum(causal),
   mean_pip_causal = mean_pip_causal,
   mean_pip_noncausal = mean_pip_noncausal,
   median_pip_causal = median(pip[causal], na.rm = TRUE),
   median_pip_noncausal = median(pip[!causal], na.rm = TRUE),
   pip_enrichment = mean_pip_causal / mean_pip_noncausal,
   auc = auc_binary(pip, causal),
   average_precision = average_precision(pip, causal),
   brier_score = mean((pip - y)^2, na.rm = TRUE),
   log_loss = safe_log_loss(pip, y),
   max_pip_noncausal = max(pip[!causal], na.rm = TRUE),
   min_pip_causal = min(pip[causal], na.rm = TRUE)
  )
 })

 thresholds <- lapply(seq_len(nt), function(t) {
  pip <- dm[, t]
  causal <- Btrue[, t] != 0
  n_causal <- sum(causal)

  do.call(rbind, lapply(pip_thresholds, function(thr) {
   selected <- pip >= thr
   n_selected <- sum(selected, na.rm = TRUE)
   recovered <- sum(causal[selected], na.rm = TRUE)
   false_selected <- n_selected - recovered

   precision <- if (n_selected > 0) recovered / n_selected else NA_real_
   recall <- recovered / n_causal

   data.frame(
    trait = t,
    trait_name = trait_names[t],
    threshold = thr,
    n_selected = n_selected,
    recovered = recovered,
    false_selected = false_selected,
    false_per_true = false_selected / pmax(recovered, 1),
    precision = precision,
    recall = recall,
    f1 = if (!is.na(precision) && (precision + recall) > 0) {
     2 * precision * recall / (precision + recall)
    } else {
     NA_real_
    }
   )
  }))
 })

 out <- list(
  global = do.call(rbind, global),
  thresholds = do.call(rbind, thresholds)
 )

 rownames(out$global) <- NULL
 rownames(out$thresholds) <- NULL

 out
}

fit_list <- list(
 ST = fit_st,
 ST_scheduled = fit_sched,
 BED_sparse = fit_bed_sparse,
 BED_scheduled = fit_bed_sched
)

dm_cmp <- do.call(rbind, lapply(names(fit_list), function(method) {
 x <- dm_estimation_summary(fit_list[[method]], sim)$global
 x$method <- method
 x
}))

dm_cmp[, c(
 "method", "trait_name",
 "pip_sum", "pip_sum_over_n_causal",
 "mean_pip_causal", "mean_pip_noncausal",
 "auc", "average_precision",
 "brier_score", "log_loss",
 "max_pip_noncausal", "min_pip_causal"
)]


dm_thr_cmp <- do.call(rbind, lapply(names(fit_list), function(method) {
 x <- dm_estimation_summary(fit_list[[method]], sim)$thresholds
 x$method <- method
 x
}))

dm_thr_cmp[dm_thr_cmp$threshold %in% c(0.01, 0.05, 0.1, 0.5), ]


library(ggplot2)

ggplot(dm_cmp, aes(x = method, y = pip_sum, fill = method)) +
 geom_col() +
 geom_hline(aes(yintercept = n_causal), linetype = 2) +
 facet_wrap(~ trait_name) +
 labs(
  x = NULL,
  y = "Sum of PIPs",
  title = "Posterior expected number of included markers"
 ) +
 theme_bw() +
 theme(
  axis.text.x = element_text(angle = 45, hjust = 1),
  legend.position = "none"
 )


ggplot(
 dm_thr_cmp[dm_thr_cmp$threshold == 0.01, ],
 aes(x = method, y = false_per_true, fill = method)
) +
 geom_col() +
 facet_wrap(~ trait_name) +
 labs(
  x = NULL,
  y = "False selected per recovered causal",
  title = "Low-threshold PIP false discovery burden"
 ) +
 theme_bw() +
 theme(
  axis.text.x = element_text(angle = 45, hjust = 1),
  legend.position = "none"
 )



pip_df <- bind_rows(lapply(names(fit_list), function(method) {
 fit <- fit_list[[method]]
 data.frame(
  method = method,
  marker = seq_len(nrow(fit$dm)),
  fit$dm
 ) |>
  tidyr::pivot_longer(
   cols = -c(method, marker),
   names_to = "trait_name",
   values_to = "pip"
  )
}))

ggplot(pip_df, aes(x = pip + 1e-6, color = method)) +
 stat_ecdf() +
 scale_x_log10() +
 facet_wrap(~ trait_name) +
 labs(
  x = "PIP + 1e-6",
  y = "ECDF",
  title = "PIP distribution: CSR has heavier low-PIP tail"
 ) +
 theme_bw()



ggplot(pip_df, aes(x = method, y = pip, fill = method)) +
 geom_violin(scale = "width") +
 scale_y_log10() +
 facet_wrap(~ trait_name) +
 labs(
  x = NULL,
  y = "PIP",
  title = "Low-PIP tail comparison"
 ) +
 theme_bw() +
 theme(
  axis.text.x = element_text(angle = 45, hjust = 1),
  legend.position = "none"
 )






m <- length(cls)
n <- Glist$n
nt <- length(stats$yy)

b <- lapply(seq_len(nt), function(x) rep(0, m))

models <- rep(list(0:1), nt)
models <- t(do.call(expand.grid, models))
models <- split(models, rep(seq_len(ncol(models)), each = nrow(models)))

pi_marker <- 0.001
pimodels <- c(
 1 - pi_marker,
 rep(pi_marker / (length(models) - 1), length(models) - 1)
)


nub <- nue <- 4

vy <- diag(stats$yy) / (n - 1)

h2 <- 0.5
vg <- diag(diag(vy) * h2)
ve <- diag(diag(vy) * (1 - h2))

expected_active <- m * pi_marker

vb <- diag(diag(vg) / expected_active)

ssb_prior <- diag(((nub - 2.0) / nub) * diag(vg) / expected_active)
sse_prior <- diag(((nue - 2.0) / nue) * diag(ve))
trait_names <- names(stats$yy)

if(is.null(trait_names)) trait_names <- paste0("T", seq_len(length(stats$yy)))

rownames(vy) <- colnames(vy) <- trait_names
rownames(vg) <- colnames(vg) <- trait_names
rownames(ve) <- colnames(ve) <- trait_names
rownames(vb) <- colnames(vb) <- trait_names
rownames(ssb_prior) <- colnames(ssb_prior) <- trait_names
rownames(sse_prior) <- colnames(sse_prior) <- trait_names



system.time(fit <- mtblr_cpg_omp_csr(
 wy = stats$wy,
 ww = stats$ww,
 yy = split(diag(stats$yy), rep(1:length(stats$yy), each = length(stats$yy))),
 b_init = b,
 ld_prefix = file.path(
  "C:/Users/au223366/Documents/GitHub/examples/human",
  "ld_test"
 ),
 B = vb,
 E = ve,
 ssb_prior=split(ssb_prior, rep(1:ncol(ssb_prior), each = nrow(ssb_prior))),
 sse_prior=split(sse_prior, rep(1:ncol(sse_prior), each = nrow(sse_prior))),
 models = models,
 pi = pimodels,
 nub = 4,
 nue = 4,
 updateB = TRUE,
 updateE = TRUE,
 updatePi = TRUE,
 adjE = 0.9,
 n = rep(as.integer(n), nt),
 nit = 500,
 nburn = 100,
 nthin = 1,
 seed = 10,
 method = 4
))

names(fit) <- c("bm","dm","wy","r","b","d","o",
                "vbs","vgs","ves",
                "covb","covg","cove",
                "vb","vg","ve",
                "pi","pim","pitrait","pimarker")

trait_names <- names(stats$yy)
if(is.null(trait_names)) trait_names <- paste0("T",1:nt)
variable_names <- names(stats$ww[[1]])
if(is.null(variable_names)) variable_names <- paste0("V",1:length(stats$ww[[1]]))

for(i in 1:7){
 fit[[i]] <- as.matrix(as.data.frame(fit[[i]]))
 rownames(fit[[i]]) <- variable_names
 colnames(fit[[i]]) <- trait_names
}
for(i in 8:10){
 fit[[i]] <- as.matrix(as.data.frame(fit[[i]]))
 rownames(fit[[i]]) <- paste0("Iter",1:nrow(fit[[i]]))
 colnames(fit[[i]]) <- trait_names
}

for(i in 11:16){
 fit[[i]] <- matrix(unlist(fit[[i]]), ncol = nt, byrow = TRUE)
 colnames(fit[[i]]) <- rownames(fit[[i]]) <- trait_names
}
fit[[17]] <- fit[[17]][[1]]
fit[[18]] <- fit[[18]][[1]]
names(fit[[17]]) <- sapply(models,paste,collapse="_")
names(fit[[18]]) <- sapply(models,paste,collapse="_")
fit <- fit[1:18]
if(sum(diag(fit$covb))>0) fit$rb <- cov2cor(fit$covb)
if(sum(diag(fit$covg))>0) fit$rg <- cov2cor(fit$covg)
if(sum(diag(fit$cove))>0) fit$re <- cov2cor(fit$cove)


colSums(fit$dm)

sapply(c(0.01, 0.05, 0.1, 0.2, 0.5),
       function(th) sum(fit$dm[,1] > th))

cor(fit$bm[,1], sim$B[,1])
cor(fit$dm[,1], fitC$stat$dm)
cor(fitC$stat$bm, sim$B[,1])
plot(fit$dm[,1], fitC$stat$dm)
plot(fit$bm[,1], sim$B[,1])
plot(fitC$stat$bm, sim$B[,1])
plot(abs(fit$bm[,1]), fit$dm[,1])

colMeans(fit$vbs)
colMeans(fit$vgs)
colMeans(fit$ves)
fit$rg
fit$rb
fit$re

for (t in 1:3) {
 causal_t <- sim$B[,t] != 0
 cat("\nTrait", t, "\n")

 for (K in c(40, 60, 100, 150, 250)) {
  top_dm <- order(fit$dm[,t], decreasing = TRUE)[1:K]
  top_bm <- order(abs(fit$bm[,t]), decreasing = TRUE)[1:K]

  cat("K =", K,
      " dm:", sum(causal_t[top_dm]), "/", K,
      " |bm|:", sum(causal_t[top_bm]), "/", K, "\n")
 }
}


colMeans(fit$vgs) + colMeans(fit$ves)

causal_any <- rowSums(sim$B != 0) > 0
cor(sim$B[causal_any, ])


m <- length(cls)
n <- Glist$n
nt <- length(stats$yy)

b <- lapply(seq_len(nt), function(x) rep(0, m))
d <- lapply(seq_len(nt), function(x) rep(0, m))

nub <- nue <- 4

vy <- stats$yy / (n - 1)
h2 <- 0.5
pi <- 0.001

vb <- diag((vy * h2) / (m * pi))
ve <- diag(vy * (1 - h2))

ssb_prior <- diag(((nub - 2) / nub) * (vy * h2) / (m * pi))
sse_prior <- diag(((nue - 2) / nue) * (vy * (1 - h2)))

Sys.setenv(
 MKL_NUM_THREADS = "1",
 MKL_DYNAMIC = "FALSE",
 OMP_NUM_THREADS = "3",
 OMP_DYNAMIC = "FALSE"
)

system.time(fit <- stblr_cpg_omp_csr(
 wy = stats$wy,
 ww = stats$ww,
 yy = stats$yy,
 b_init = b,
 d_init = d,
 use_d_init = FALSE,
 r_init = stats$wy,
 use_r_init = FALSE,
 rebuild_r_before_updateE = FALSE,
 #ld_row_ptr = ld$row_ptr,
 #ld_col_idx = ld$col_idx,
 #ld_values = ld$values,
 #ld_col_idx_one_based = TRUE,
 ld_prefix = file.path(
  "C:/Users/au223366/Documents/GitHub/examples/human",
  "ld_test"
 ),
 B = vb,
 E = ve,
 ssb_prior = split(ssb_prior, rep(1:ncol(ssb_prior), each = nrow(ssb_prior))),
 sse_prior = split(sse_prior, rep(1:ncol(sse_prior), each = nrow(sse_prior))),
 pi = c(1 - pi, pi),
 nub = 4,
 nue = 4,
 updateB = TRUE,
 updateE = TRUE,
 updatePi = TRUE,
 adjE = 0.9,
 n = rep(as.integer(n), nt),
 nit = 1000,
 nburn = 100,
 nthin = 1,
 ncores=3,
 seed = 10
))

system.time(fit <- stblr_cpg_omp_csr_scheduled(
 wy = stats$wy,
 ww = stats$ww,
 yy = stats$yy,
 b_init = b,
 d_init = d,
 use_d_init = FALSE,
 r_init = stats$wy,
 use_r_init = FALSE,
 rebuild_r_before_updateE = FALSE,
 ld_prefix = file.path(
  "C:/Users/au223366/Documents/GitHub/examples/human",
  "ld_test"
 ),
 B = vb,
 E = ve,
 ssb_prior = split(ssb_prior, rep(1:ncol(ssb_prior), each = nrow(ssb_prior))),
 sse_prior = split(sse_prior, rep(1:ncol(sse_prior), each = nrow(sse_prior))),
 pi = c(1 - pi_marker, pi_marker),
 nub = 4,
 nue = 4,
 updateB = TRUE,
 updateE = TRUE,
 updatePi = TRUE,
 adjE = 0.9,
 n = rep(as.integer(n), nt),
 nit = 1000,
 nburn = 100,
 nthin = 1,
 full_sweep_every = 10,
 null_skip_base = 50,
 null_skip_max = 200,
 candidate_threshold = 1e-3,
 candidate_lifetime = 20,
 skip_nulls_burnin_only = FALSE,
 wakeup_ld_neighbors = TRUE,
 wakeup_diff_threshold = 0.0,
 wakeup_max_neighbors = 0,
 ncores = 3,
 seed = 10
))

names(fit) <- c("bm","dm","wy","r","b","d","o",
                "vbs","vgs","ves",
                "covb","covg","cove",
                "vb","vg","ve",
                "pi","pim","pitrait","pimarker")

trait_names <- names(stats$yy)
if(is.null(trait_names)) trait_names <- paste0("T",1:nt)
variable_names <- names(stats$ww[[1]])
if(is.null(variable_names)) variable_names <- paste0("V",1:length(stats$ww[[1]]))

for(i in 1:7){
 fit[[i]] <- as.matrix(as.data.frame(fit[[i]]))
 rownames(fit[[i]]) <- variable_names
 colnames(fit[[i]]) <- trait_names
}
for(i in 8:10){
 fit[[i]] <- as.matrix(as.data.frame(fit[[i]]))
 rownames(fit[[i]]) <- paste0("Iter",1:nrow(fit[[i]]))
 colnames(fit[[i]]) <- trait_names
}

for(i in 11:16){
 fit[[i]] <- matrix(unlist(fit[[i]]), ncol = nt, byrow = TRUE)
 colnames(fit[[i]]) <- rownames(fit[[i]]) <- trait_names
}
fit[[17]] <- fit[[17]][[1]]
fit[[18]] <- fit[[18]][[1]]
fit <- fit[1:18]
if(sum(diag(fit$covb))>0) fit$rb <- cov2cor(fit$covb)
if(sum(diag(fit$covg))>0) fit$rg <- cov2cor(fit$covg)
if(sum(diag(fit$cove))>0) fit$re <- cov2cor(fit$cove)

colSums(fit$dm)

sapply(c(0.01, 0.05, 0.1, 0.2, 0.5),
       function(th) sum(fit$dm[,1] > th))

cor(fit$bm[,1], sim$B[,1])
cor(fit$dm[,1], fitC$stat$dm)
cor(fitC$stat$bm, sim$B[,1])
plot(fit$dm[,1], fitC$stat$dm)
plot(fit$bm[,1], sim$B[,1])
plot(fitC$stat$bm, sim$B[,1])

colMeans(fit$vbs)
colMeans(fit$vgs)
colMeans(fit$ves)
fit$rg
fit$rb
fit$re

for (t in 1:3) {
 causal_t <- sim$B[,t] != 0
 cat("\nTrait", t, "\n")

 for (K in c(40, 60, 100, 150, 250)) {
  top_dm <- order(fit$dm[,t], decreasing = TRUE)[1:K]
  top_bm <- order(abs(fit$bm[,t]), decreasing = TRUE)[1:K]

  cat("K =", K,
      " dm:", sum(causal_t[top_dm]), "/", K,
      " |bm|:", sum(causal_t[top_bm]), "/", K, "\n")
 }
}


causal_any <- rowSums(sim$B != 0) > 0
cor(sim$B[causal_any, ])

colMeans(fit$vgs) + colMeans(fit$ves)
colSums(fit$dm)
cor(fit$bm[,1], sim$B[,1])
cor(fit$dm[,1], fitC$stat$dm)


vg_set <- estimate_vg_by_set(fitST, sets, n)
prop_set <- sweep(vg_set, 2, colSums(vg_set), "/")
prop_of_vg <- sweep(vg_set, 2, colMeans(fitST$vgs), "/")
prop_of_vp <- vg_set  # if phenotypes are standardized to variance ~1

vg_set_bm <- estimate_vg_by_set(fitST, sets, n)
vg_set_dm_weighted <- estimate_weighted_vg_by_set(fitST, sets, n)
pip_enrichment <- estimate_pip_enrichment(fitST$dm, sets)

cov2cor_safe <- function(G) {
 d <- diag(G)
 if (any(d <= 0)) return(matrix(NA_real_, nrow(G), ncol(G)))
 cov2cor(G)
}

G_sets <- estimate_G_by_set(fit, sets, n)
R_sets <- lapply(G_sets, cov2cor_safe)


shared_dm <- apply(fit$dm, 1, min)       # evidence shared across all traits
any_dm    <- apply(fit$dm, 1, max)       # evidence in any trait
sum_dm    <- rowSums(fit$dm)             # total inclusion burden

# Marker-level ST summaries
b_st  <- fitST$bm
dm_st <- fitST$dm
q_st  <- fitST$wy - fitST$r

# Set-level variance
vg_sets <- estimate_vg_by_set(fitST, sets, n)

# Set-level PIP enrichment
pip_enrich <- estimate_pip_enrichment(dm_st, sets)

# Optional: MT covariance decomposition
G_sets <- estimate_G_by_set(fitMT, sets, n)
R_sets <- lapply(G_sets, cov2cor_safe)


##################################################################


# format_stblr_fit <- function(fit, nt, m, trait_names = NULL, variable_names = NULL) {
#  names(fit) <- c(
#   "bm", "dm", "wy", "r", "b", "d", "o",
#   "vbs", "vgs", "ves",
#   "covb", "covg", "cove",
#   "vb", "vg", "ve",
#   "pi", "pim", "pitrait", "pimarker"
#  )
#
#  if (is.null(trait_names)) {
#   trait_names <- paste0("T", seq_len(nt))
#  }
#
#  if (is.null(variable_names)) {
#   variable_names <- paste0("V", seq_len(m))
#  }
#
#  for (i in 1:7) {
#   fit[[i]] <- as.matrix(as.data.frame(fit[[i]]))
#   rownames(fit[[i]]) <- variable_names
#   colnames(fit[[i]]) <- trait_names
#  }
#
#  for (i in 8:10) {
#   fit[[i]] <- as.matrix(as.data.frame(fit[[i]]))
#   rownames(fit[[i]]) <- paste0("Iter", seq_len(nrow(fit[[i]])))
#   colnames(fit[[i]]) <- trait_names
#  }
#
#  for (i in 11:16) {
#   fit[[i]] <- matrix(unlist(fit[[i]]), ncol = nt, byrow = TRUE)
#   colnames(fit[[i]]) <- rownames(fit[[i]]) <- trait_names
#  }
#
#  fit[[17]] <- fit[[17]][[1]]
#  fit[[18]] <- fit[[18]][[1]]
#
#  fit <- fit[1:18]
#
#  if (sum(diag(fit$covb)) > 0) fit$rb <- cov2cor(fit$covb)
#  if (sum(diag(fit$covg)) > 0) fit$rg <- cov2cor(fit$covg)
#  if (sum(diag(fit$cove)) > 0) fit$re <- cov2cor(fit$cove)
#
#  fit
# }
format_stblr_fit <- function(fit, nt, m, trait_names = NULL, variable_names = NULL) {
 has_vle_vld <- length(fit) >= 22

 if (has_vle_vld) {
  names(fit)[1:22] <- c(
   "bm", "dm", "wy", "r", "b", "d", "o",
   "vbs", "vgs", "ves",
   "covb", "covg", "cove",
   "vb", "vg", "ve",
   "pi", "pim", "pitrait", "pimarker",
   "vle", "vld"
  )
 } else {
  names(fit) <- c(
   "bm", "dm", "wy", "r", "b", "d", "o",
   "vbs", "vgs", "ves",
   "covb", "covg", "cove",
   "vb", "vg", "ve",
   "pi", "pim", "pitrait", "pimarker"
  )
 }

 if (is.null(trait_names)) {
  trait_names <- paste0("T", seq_len(nt))
 }

 if (is.null(variable_names)) {
  variable_names <- paste0("V", seq_len(m))
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
  colnames(fit[[i]]) <- rownames(fit[[i]]) <- trait_names
 }

 fit[[17]] <- fit[[17]][[1]]
 fit[[18]] <- fit[[18]][[1]]

 out <- fit[1:18]

 if (has_vle_vld) {
  out$vle <- as.matrix(as.data.frame(fit[[21]]))
  out$vld <- as.matrix(as.data.frame(fit[[22]]))

  rownames(out$vle) <- paste0("Iter", seq_len(nrow(out$vle)))
  rownames(out$vld) <- paste0("Iter", seq_len(nrow(out$vld)))

  colnames(out$vle) <- trait_names
  colnames(out$vld) <- trait_names
 }

 if (!is.null(out$pitrait) && is.matrix(out$pitrait) && ncol(out$pitrait) >= 4) {
  out$diagnostics <- out$pitrait

  colnames(out$diagnostics) <- c(
   "log_cpo",
   "mean_log_cpo",
   "seconds_mean",
   "seconds_max"
  )

  rownames(out$diagnostics) <- trait_names

  out$log_cpo <- out$diagnostics[, "log_cpo"]
  out$mean_log_cpo <- out$diagnostics[, "mean_log_cpo"]
 }

 if (sum(diag(out$covb)) > 0) out$rb <- cov2cor(out$covb)
 if (sum(diag(out$covg)) > 0) out$rg <- cov2cor(out$covg)
 if (sum(diag(out$cove)) > 0) out$re <- cov2cor(out$cove)

 out
}

library(qgg)
library(sblr)

data_dir <- "C:/Users/au223366/Documents/GitHub/examples/human"

Glist <- readRDS(file.path(data_dir, "Glist_sparseLD_1k.RDS"))

# -----------------------------------------------------------------------------
# Simulate traits
# -----------------------------------------------------------------------------

chr <- 1

rg <- matrix(
 c(
  1.0, 0.7, 0.3,
  0.7, 1.0, 0.5,
  0.3, 0.5, 1.0
 ),
 nrow = 3,
 byrow = TRUE
)

sim <- mtsim(
 Glist = Glist,
 chr = chr,
 rsids = Glist$rsidsLD[[chr]],
 nt = 3,
 n_shared = 30,
 n_specific = 10,
 h2 = c(0.4, 0.5, 0.3),
 rg = rg,
 re = 0,
 seed = 1
)

y <- as.matrix(scale(sim$y))

stat <- glma(y = scale(sim$y[,1]), rsids=Glist$rsidsLD[[1]], Glist = Glist)
system.time(fitC <- gbayes(stat = stat, Glist = Glist, method = "bayesC", nit = 1000))


# -----------------------------------------------------------------------------
# Marker subset and sufficient statistics for diagnostics only
# -----------------------------------------------------------------------------

cls <- match(Glist$rsidsLD[[chr]], Glist$rsids[[chr]])


cls <- match(Glist$rsidsLD[[chr]], Glist$rsids[[chr]])
chr <- 1
system.time(stats <- bed_xtx_xty(
 bed_file = Glist$bedfiles[chr],
 n = Glist$n,
 cls = cls,
 af = Glist$af[[chr]][cls],
 y = y,
 scale = TRUE,
 nthreads = 4
))


m <- length(cls)
n <- Glist$n
nt <- ncol(y)

# -----------------------------------------------------------------------------
# Block definition for new BED-backed sampler
# -----------------------------------------------------------------------------

block_size <- 1000

sets <- rep(
 seq_len(ceiling(m / block_size)),
 each = block_size
)[seq_len(m)]

table(sets)

# -----------------------------------------------------------------------------
# Prior setup for ST-BLR
# -----------------------------------------------------------------------------

b <- lapply(seq_len(nt), function(t) rep(0, m))

nub <- 4
nue <- 4

vy <- stats$yy / (n - 1)

h2 <- 0.5
pi_marker <- 0.001

vb <- diag((vy * h2) / (m * pi_marker), nt, nt)
ve <- diag(vy * (1 - h2), nt, nt)

ssb_prior <- diag(((nub - 2) / nub) * (vy * h2) / (m * pi_marker), nt, nt)
sse_prior <- diag(((nue - 2) / nue) * (vy * (1 - h2)), nt, nt)

trait_names <- colnames(y)
if (is.null(trait_names)) trait_names <- paste0("T", seq_len(nt))

rownames(vb) <- colnames(vb) <- trait_names
rownames(ve) <- colnames(ve) <- trait_names
rownames(ssb_prior) <- colnames(ssb_prior) <- trait_names
rownames(sse_prior) <- colnames(sse_prior) <- trait_names

ssb_prior_list <- split(
 ssb_prior,
 rep(seq_len(ncol(ssb_prior)), each = nrow(ssb_prior))
)

sse_prior_list <- split(
 sse_prior,
 rep(seq_len(ncol(sse_prior)), each = nrow(sse_prior))
)

# -----------------------------------------------------------------------------
# Run new BED-backed block ST-BLR
# -----------------------------------------------------------------------------

Sys.setenv(OMP_NUM_THREADS = "3")
Sys.setenv(OMP_DYNAMIC = "FALSE")

fit <- stblr_cpg_omp_bed_marker_sparse(
 bed_files = Glist$bedfiles[chr],
 n = Glist$n,
 cls = list(cls),
 y = y,
 b_init = b,
 sets = sets,
 rows = NULL,
 af = list(Glist$af[[chr]][cls]),
 scale = TRUE,
 B = vb,
 E = ve,
 ssb_prior = split(ssb_prior, rep(1:ncol(ssb_prior), each = nrow(ssb_prior))),
 sse_prior = split(sse_prior, rep(1:ncol(sse_prior), each = nrow(sse_prior))),
 pi = c(1 - pi, pi),
 nub = 4,
 nue = 4,
 updateB = TRUE,
 updateE = TRUE,
 updatePi = TRUE,
 adjE = 0.9,
 nit = 1000,
 nburn = 100,
 nthin = 1,
 rebuild_every = 25,
 full_sweep_every = 10,
 null_update_prob = 0.02,
 candidate_threshold = 1e-3,
 candidate_lifetime = 20,
 skip_nulls_burnin_only = TRUE,
 ncores = 3,
 seed = 10
)

system.time(
 fitBED_raw <- stblr_cpg_omp_bed_blocks(
  bed_files = Glist$bedfiles[chr],
  n = Glist$n,
  cls = list(cls),
  y = y,
  b_init = b,
  sets = sets,
  rows = NULL,
  af = list(Glist$af[[chr]][cls]),
  scale = TRUE,
  B = vb,
  E = ve,
  ssb_prior = ssb_prior_list,
  sse_prior = sse_prior_list,
  pi = c(1 - pi_marker, pi_marker),
  nub = nub,
  nue = nue,
  updateB = TRUE,
  updateE = TRUE,
  updatePi = TRUE,
  adjE = 0,
  nit = 100,
  nburn = 50,
  nthin = 10,
  block_nit = 1,
  rebuild_every = 0,
  ncores = 3,
  seed = 10
 )
)

system.time(
 fitBED_raw_resid <- stblr_cpg_omp_bed_residual_blocks(
 bed_files = Glist$bedfiles[chr],
 n = Glist$n,
 cls = list(cls),
 y = y,
 b_init = b,
 sets = sets,
 rows = NULL,
 af = list(Glist$af[[chr]][cls]),
 scale = TRUE,
 B = vb,
 E = ve,
 ssb_prior = ssb_prior_list,
 sse_prior = sse_prior_list,
 pi = c(1 - pi_marker, pi_marker),
 nub = nub,
 nue = nue,
 updateB = TRUE,
 updateE = TRUE,
 updatePi = TRUE,
 adjE = 0,
 nit = 100,
 nburn = 50,
 nthin = 10,
 block_nit = 1,
 rebuild_every = 0,
 ncores = 3,
 seed = 10
))

# -----------------------------------------------------------------------------
# Basic diagnostics
# -----------------------------------------------------------------------------

fitBED <- format_stblr_fit(fitBED_raw, nt=3, m=m)

colSums(fitBED$dm)

sapply(
 c(0.01, 0.05, 0.1, 0.2, 0.5),
 function(th) sum(fitBED$dm[, 1] > th)
)

cor(fitBED$bm[, 1], sim$B[, 1])
cor(fitBED$dm[, 1], fitC$stat$dm)
cor(fitC$stat$bm, sim$B[, 1])

plot(fitBED$dm[, 1], fitC$stat$dm,
     xlab = "BED-block STBLR dm",
     ylab = "gbayes dm")

plot(fitBED$bm[, 1], sim$B[, 1],
     xlab = "BED-block STBLR bm",
     ylab = "True simulated B")

plot(fitC$stat$bm, sim$B[, 1],
     xlab = "gbayes bm",
     ylab = "True simulated B")

colMeans(fitBED$vbs)
colMeans(fitBED$vgs)
colMeans(fitBED$ves)

fitBED$rg
fitBED$rb
fitBED$re

colMeans(fitBED$vgs) + colMeans(fitBED$ves)

for (t in 1:3) {
 causal_t <- sim$B[, t] != 0
 cat("\nTrait", t, "\n")

 for (K in c(40, 60, 100, 150, 250)) {
  top_dm <- order(fitBED$dm[, t], decreasing = TRUE)[1:K]
  top_bm <- order(abs(fitBED$bm[, t]), decreasing = TRUE)[1:K]

  cat(
   "K =", K,
   " dm:", sum(causal_t[top_dm]), "/", K,
   " |bm|:", sum(causal_t[top_bm]), "/", K,
   "\n"
  )
 }
}

recovery_summary <- function(fit, Btrue, Ks = c(40, 60, 100, 150, 250)) {
 nt <- ncol(Btrue)

 out <- do.call(rbind, lapply(seq_len(nt), function(t) {
  causal <- Btrue[, t] != 0
  n_causal <- sum(causal)

  do.call(rbind, lapply(Ks, function(K) {
   top_dm <- order(fit$dm[, t], decreasing = TRUE)[seq_len(K)]
   top_bm <- order(abs(fit$bm[, t]), decreasing = TRUE)[seq_len(K)]

   data.frame(
    trait = t,
    K = K,
    n_causal = n_causal,

    recovered_dm = sum(causal[top_dm]),
    precision_dm = sum(causal[top_dm]) / K,
    recall_dm = sum(causal[top_dm]) / n_causal,

    recovered_bm = sum(causal[top_bm]),
    precision_bm = sum(causal[top_bm]) / K,
    recall_bm = sum(causal[top_bm]) / n_causal
   )
  }))
 }))

 out
}

rec <- recovery_summary(fitBED, sim$B)
rec

proxy_recovery <- function(fit, Btrue, pos, window = 50000, pip_threshold = 0.1) {
 nt <- ncol(Btrue)

 out <- lapply(seq_len(nt), function(t) {
  causal <- which(Btrue[, t] != 0)
  selected <- which(fit$dm[, t] > pip_threshold)

  recovered_exact <- causal %in% selected

  recovered_proxy <- vapply(causal, function(i) {
   any(abs(pos[selected] - pos[i]) <= window)
  }, logical(1))

  data.frame(
   trait = t,
   n_causal = length(causal),
   exact = sum(recovered_exact),
   proxy = sum(recovered_proxy),
   exact_recall = mean(recovered_exact),
   proxy_recall = mean(recovered_proxy)
  )
 })

 do.call(rbind, out)
}

proxy_recovery(
 fit = fitBED,
 Btrue = sim$B,
 pos = Glist$pos[[chr]][cls],
 window = 50000,
 pip_threshold = 0.1
)



# -----------------------------------------------------------------------------
# Marker subset and sufficient statistics for diagnostics only
# -----------------------------------------------------------------------------

chr <- 1

cls <- match(Glist$rsidsLD[[chr]], Glist$rsids[[chr]])

if (anyNA(cls)) {
 stop("Some LD rsids were not found in Glist$rsids for this chromosome.")
}

system.time(stats <- bed_xtx_xty(
 bed_file = Glist$bedfiles[chr],
 n = Glist$n,
 cls = cls,
 af = Glist$af[[chr]][cls],
 y = y,
 scale = TRUE,
 nthreads = 4
))

m <- length(cls)
n <- Glist$n
nt <- ncol(y)

# -----------------------------------------------------------------------------
# Block/order definition for new BED marker-wise sampler
# -----------------------------------------------------------------------------

block_size <- 1000

sets <- rep(
 seq_len(ceiling(m / block_size)),
 each = block_size
)[seq_len(m)]

table(sets)

# -----------------------------------------------------------------------------
# Prior setup for ST-BLR
# -----------------------------------------------------------------------------

b <- lapply(seq_len(nt), function(t) rep(0, m))

nub <- 4
nue <- 4

h2 <- 0.5
pi_marker <- 0.001

yy_diag <- if (is.matrix(stats$yy)) diag(stats$yy) else stats$yy

vy <- yy_diag / (n - 1)

vb <- diag((vy * h2) / (m * pi_marker), nt, nt)
ve <- diag(vy * (1 - h2), nt, nt)

ssb_prior <- diag(
 ((nub - 2) / nub) * (vy * h2) / (m * pi_marker),
 nt,
 nt
)

sse_prior <- diag(
 ((nue - 2) / nue) * (vy * (1 - h2)),
 nt,
 nt
)

trait_names <- colnames(y)
if (is.null(trait_names)) trait_names <- paste0("T", seq_len(nt))

rownames(vb) <- colnames(vb) <- trait_names
rownames(ve) <- colnames(ve) <- trait_names
rownames(ssb_prior) <- colnames(ssb_prior) <- trait_names
rownames(sse_prior) <- colnames(sse_prior) <- trait_names

ssb_prior_list <- split(
 ssb_prior,
 rep(seq_len(ncol(ssb_prior)), each = nrow(ssb_prior))
)

sse_prior_list <- split(
 sse_prior,
 rep(seq_len(ncol(sse_prior)), each = nrow(sse_prior))
)

# -----------------------------------------------------------------------------
# Run new packed-BED marker-wise ST-BLR
# -----------------------------------------------------------------------------

Sys.setenv(
 OMP_NUM_THREADS = "3",
 OMP_DYNAMIC = "FALSE",
 MKL_NUM_THREADS = "1",
 MKL_DYNAMIC = "FALSE",
 OPENBLAS_NUM_THREADS = "1"
)

system.time(fit <- stblr_cpg_omp_bed_marker_sparse(
 bed_files = Glist$bedfiles[chr],
 n = Glist$n,
 cls = list(cls),
 y = y,
 b_init = b,
 sets = sets,
 rows = NULL,
 af = list(Glist$af[[chr]][cls]),
 scale = TRUE,
 B = vb,
 E = ve,
 ssb_prior = ssb_prior_list,
 sse_prior = sse_prior_list,
 pi = c(1 - pi_marker, pi_marker),
 nub = nub,
 nue = nue,
 updateB = TRUE,
 updateE = TRUE,
 updatePi = TRUE,
 adjE = 0.9,
 nit = 1000,
 nburn = 100,
 nthin = 1,
 rebuild_every = 25,
 full_sweep_every = 10,
 null_update_prob = 0.02,
 candidate_threshold = 1e-3,
 candidate_lifetime = 20,
 skip_nulls_burnin_only = FALSE,
 ncores = 3,
 seed = 10
))

system.time(fit <- stblr_cpg_omp_bed_marker_scheduled(
 bed_files = Glist$bedfiles[chr],
 n = Glist$n,
 cls = list(cls),
 y = y,
 b_init = b,
 sets = sets,
 rows = NULL,
 af = list(Glist$af[[chr]][cls]),
 scale = TRUE,
 B = vb,
 E = ve,
 ssb_prior = ssb_prior_list,
 sse_prior = sse_prior_list,
 pi = c(1 - pi_marker, pi_marker),
 nub = nub,
 nue = nue,
 updateB = TRUE,
 updateE = TRUE,
 updatePi = TRUE,
 adjE = 0.9,
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
 ncores = 3,
 seed = 10
))


# -----------------------------------------------------------------------------
# Basic diagnostics
# -----------------------------------------------------------------------------

names(fit) <- c("bm","dm","wy","r","b","d","o",
                "vbs","vgs","ves",
                "covb","covg","cove",
                "vb","vg","ve",
                "pi","pim","pitrait","pimarker")

trait_names <- colnames(stats$yy)
if(is.null(trait_names)) trait_names <- paste0("T",1:nt)
variable_names <- names(stats$ww[[1]])
if(is.null(variable_names)) variable_names <- paste0("V",1:length(stats$ww[[1]]))

for(i in 1:7){
 fit[[i]] <- as.matrix(as.data.frame(fit[[i]]))
 rownames(fit[[i]]) <- variable_names
 colnames(fit[[i]]) <- trait_names
}
for(i in 8:10){
 fit[[i]] <- as.matrix(as.data.frame(fit[[i]]))
 rownames(fit[[i]]) <- paste0("Iter",1:nrow(fit[[i]]))
 colnames(fit[[i]]) <- trait_names
}

for(i in 11:16){
 fit[[i]] <- matrix(unlist(fit[[i]]), ncol = nt, byrow = TRUE)
 colnames(fit[[i]]) <- rownames(fit[[i]]) <- trait_names
}
fit[[17]] <- fit[[17]][[1]]
fit[[18]] <- fit[[18]][[1]]
fit <- fit[1:18]
if(sum(diag(fit$covb))>0) fit$rb <- cov2cor(fit$covb)
if(sum(diag(fit$covg))>0) fit$rg <- cov2cor(fit$covg)
if(sum(diag(fit$cove))>0) fit$re <- cov2cor(fit$cove)

colSums(fit$dm)

sapply(c(0.01, 0.05, 0.1, 0.2, 0.5),
       function(th) sum(fit$dm[,1] > th))

cor(fit$bm[,1], sim$B[,1])
cor(fit$dm[,1], fitC$stat$dm)
cor(fitC$stat$bm, sim$B[,1])
plot(fit$dm[,1], fitC$stat$dm)
plot(fit$bm[,1], sim$B[,1])
plot(fitC$stat$bm, sim$B[,1])

colMeans(fit$vbs)
colMeans(fit$vgs)
colMeans(fit$ves)
fit$rg
fit$rb
fit$re

for (t in 1:3) {
 causal_t <- sim$B[,t] != 0
 cat("\nTrait", t, "\n")

 for (K in c(40, 60, 100, 150, 250)) {
  top_dm <- order(fit$dm[,t], decreasing = TRUE)[1:K]
  top_bm <- order(abs(fit$bm[,t]), decreasing = TRUE)[1:K]

  cat("K =", K,
      " dm:", sum(causal_t[top_dm]), "/", K,
      " |bm|:", sum(causal_t[top_bm]), "/", K, "\n")
 }
}


causal_any <- rowSums(sim$B != 0) > 0
cor(sim$B[causal_any, ])

colMeans(fit$vgs) + colMeans(fit$ves)
colSums(fit$dm)
cor(fit$bm[,1], sim$B[,1])
cor(fit$dm[,1], fitC$stat$dm)


compare_stblr_priors <- function(fit_csr, fit_bed, tolerance = 1e-10) {
 x <- fit_csr$input
 y <- fit_bed$input

 checks <- list(
  n_equal = identical(x$n, y$n),
  m_equal = identical(x$m, y$m),
  nt_equal = identical(x$nt, y$nt),
  pi_marker_equal = isTRUE(all.equal(x$pi_marker, y$pi_marker, tolerance = tolerance)),
  h2_equal = isTRUE(all.equal(x$h2, y$h2, tolerance = tolerance)),
  nub_equal = isTRUE(all.equal(x$nub, y$nub, tolerance = tolerance)),
  nue_equal = isTRUE(all.equal(x$nue, y$nue, tolerance = tolerance)),
  B_equal = isTRUE(all.equal(x$B, y$B, tolerance = tolerance)),
  E_equal = isTRUE(all.equal(x$E, y$E, tolerance = tolerance)),
  ssb_prior_equal = isTRUE(all.equal(x$ssb_prior, y$ssb_prior, tolerance = tolerance)),
  sse_prior_equal = isTRUE(all.equal(x$sse_prior, y$sse_prior, tolerance = tolerance)),
  updateB_equal = identical(x$updateB, y$updateB),
  updateE_equal = identical(x$updateE, y$updateE),
  updatePi_equal = identical(x$updatePi, y$updatePi),
  adjE_equal = isTRUE(all.equal(x$adjE, y$adjE, tolerance = tolerance)),
  nit_equal = identical(x$nit, y$nit),
  nburn_equal = identical(x$nburn, y$nburn),
  nthin_equal = identical(x$nthin, y$nthin),
  seed_equal = identical(x$seed, y$seed)
 )

 data.frame(
  field = names(checks),
  equal = unlist(checks),
  row.names = NULL
 )
}

recovery_summary <- function(
  fit,
  sim,
  Ks = c(40, 60, 100, 150, 250),
  trait_names = NULL,
  pip_thresholds = c(0.01, 0.05, 0.1, 0.5, 0.7, 0.9)
) {
 if (is.null(fit$dm)) stop("fit$dm is missing.")
 if (is.null(fit$bm)) stop("fit$bm is missing.")

 dm <- as.matrix(fit$dm)
 bm <- as.matrix(fit$bm)
 Btrue <- as.matrix(sim$B)

 if (!all(dim(dm) == dim(Btrue))) {
  stop("dim(fit$dm) must match dim(Btrue).")
 }

 if (!all(dim(bm) == dim(Btrue))) {
  stop("dim(fit$bm) must match dim(Btrue).")
 }

 m <- nrow(Btrue)
 nt <- ncol(Btrue)

 Ks <- sort(unique(Ks))
 Ks <- Ks[Ks <= m]

 if (length(Ks) == 0) {
  stop("No valid K values. All K values exceed the number of markers.")
 }

 if (is.null(trait_names)) {
  trait_names <- colnames(Btrue)
  if (is.null(trait_names)) {
   trait_names <- paste0("T", seq_len(nt))
  }
 }

 auc_binary <- function(score, label) {
  label <- as.logical(label)
  n_pos <- sum(label)
  n_neg <- sum(!label)

  if (n_pos == 0 || n_neg == 0) {
   return(NA_real_)
  }

  r <- rank(score, ties.method = "average")
  (sum(r[label]) - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
 }

 average_precision <- function(score, label) {
  label <- as.logical(label)
  n_pos <- sum(label)

  if (n_pos == 0) {
   return(NA_real_)
  }

  ord <- order(score, decreasing = TRUE)
  y <- label[ord]

  mean(cumsum(y)[y] / which(y))
 }

 topk <- list()
 global <- list()
 thresholds <- list()

 for (t in seq_len(nt)) {
  causal <- Btrue[, t] != 0
  n_causal <- sum(causal)

  if (n_causal == 0) {
   warning("Trait ", t, " has no causal markers in Btrue.")
   next
  }

  score_dm <- dm[, t]
  score_bm <- abs(bm[, t])

  rank_dm <- rank(-score_dm, ties.method = "average")
  rank_bm <- rank(-score_bm, ties.method = "average")

  causal_rank_dm <- rank_dm[causal]
  causal_rank_bm <- rank_bm[causal]

  mean_pip_causal <- mean(score_dm[causal], na.rm = TRUE)
  mean_pip_noncausal <- mean(score_dm[!causal], na.rm = TRUE)

  global[[t]] <- data.frame(
   trait = t,
   trait_name = trait_names[t],
   n_markers = m,
   n_causal = n_causal,
   pip_sum = sum(score_dm, na.rm = TRUE),
   mean_pip_causal = mean_pip_causal,
   mean_pip_noncausal = mean_pip_noncausal,
   pip_enrichment = mean_pip_causal / mean_pip_noncausal,
   auc_dm = auc_binary(score_dm, causal),
   auc_bm = auc_binary(score_bm, causal),
   ap_dm = average_precision(score_dm, causal),
   ap_bm = average_precision(score_bm, causal),
   mean_rank_dm = mean(causal_rank_dm),
   median_rank_dm = median(causal_rank_dm),
   best_rank_dm = min(causal_rank_dm),
   worst_rank_dm = max(causal_rank_dm),
   mean_rank_bm = mean(causal_rank_bm),
   median_rank_bm = median(causal_rank_bm),
   best_rank_bm = min(causal_rank_bm),
   worst_rank_bm = max(causal_rank_bm)
  )

  topk[[t]] <- do.call(rbind, lapply(Ks, function(K) {
   top_dm <- order(score_dm, decreasing = TRUE)[seq_len(K)]
   top_bm <- order(score_bm, decreasing = TRUE)[seq_len(K)]

   recovered_dm <- sum(causal[top_dm])
   recovered_bm <- sum(causal[top_bm])

   precision_dm <- recovered_dm / K
   recall_dm <- recovered_dm / n_causal

   precision_bm <- recovered_bm / K
   recall_bm <- recovered_bm / n_causal

   f1_dm <- if ((precision_dm + recall_dm) > 0) {
    2 * precision_dm * recall_dm / (precision_dm + recall_dm)
   } else {
    0
   }

   f1_bm <- if ((precision_bm + recall_bm) > 0) {
    2 * precision_bm * recall_bm / (precision_bm + recall_bm)
   } else {
    0
   }

   data.frame(
    trait = t,
    trait_name = trait_names[t],
    K = K,
    K_eff = min(K, m),
    n_causal = n_causal,
    recovered_dm = recovered_dm,
    precision_dm = precision_dm,
    recall_dm = recall_dm,
    f1_dm = f1_dm,
    recovered_bm = recovered_bm,
    precision_bm = precision_bm,
    recall_bm = recall_bm,
    f1_bm = f1_bm
   )
  }))

  thresholds[[t]] <- do.call(rbind, lapply(pip_thresholds, function(thr) {
   selected <- score_dm >= thr
   n_selected <- sum(selected, na.rm = TRUE)
   recovered <- sum(causal[selected], na.rm = TRUE)

   precision <- if (n_selected > 0) recovered / n_selected else NA_real_
   recall <- recovered / n_causal

   data.frame(
    trait = t,
    trait_name = trait_names[t],
    threshold = thr,
    n_selected = n_selected,
    recovered = recovered,
    precision = precision,
    recall = recall,
    f1 = if (!is.na(precision) && (precision + recall) > 0) {
     2 * precision * recall / (precision + recall)
    } else {
     NA_real_
    }
   )
  }))
 }

 out <- list(
  topk = do.call(rbind, topk),
  global = do.call(rbind, global),
  thresholds = do.call(rbind, thresholds)
 )

 rownames(out$topk) <- NULL
 rownames(out$global) <- NULL
 rownames(out$thresholds) <- NULL

 out
}

sparse_ld_risk <- function(ld, transform = c("coherence_log_abs", "coherence", "signed_sum")) {
 transform <- match.arg(transform)

 m <- ld$nrow
 row_ptr <- ld$row_ptr
 col_idx <- ld$col_idx
 values <- ld$values

 col <- if (!is.null(ld$index_base) && ld$index_base == 1) {
  col_idx
 } else {
  col_idx + 1
 }

 signed_sum <- numeric(m)
 abs_sum <- numeric(m)
 ldscore <- numeric(m)
 degree <- numeric(m)

 for (i in seq_len(m)) {
  start <- row_ptr[i] + 1
  end <- row_ptr[i + 1]

  if (end >= start) {
   idx <- start:end
   j <- col[idx]
   r <- values[idx]

   signed_sum[i] <- signed_sum[i] + sum(r)
   abs_sum[i] <- abs_sum[i] + sum(abs(r))
   ldscore[i] <- ldscore[i] + sum(r^2)
   degree[i] <- degree[i] + length(r)

   signed_sum[j] <- signed_sum[j] + r
   abs_sum[j] <- abs_sum[j] + abs(r)
   ldscore[j] <- ldscore[j] + r^2
   degree[j] <- degree[j] + 1
  }
 }

 coherence <- abs(signed_sum) / pmax(abs_sum, 1e-300)

 risk <- switch(
  transform,
  coherence_log_abs = coherence * log1p(abs_sum),
  coherence = coherence,
  signed_sum = abs(signed_sum)
 )

 # Scale to a stable, interpretable range.
 risk_scaled <- risk / median(risk[risk > 0], na.rm = TRUE)

 data.frame(
  ld_risk = risk_scaled,
  ld_risk_raw = risk,
  signed_sum = signed_sum,
  abs_sum = abs_sum,
  coherence = coherence,
  ldscore = ldscore,
  degree = degree
 )
}


