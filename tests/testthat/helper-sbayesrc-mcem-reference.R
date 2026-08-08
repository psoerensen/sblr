.mcem_ref_component_prior <- function(A, alpha, floor = 1e-12) {
  q <- matrix(pmin(pmax(pnorm(drop(A %*% alpha)), floor), 1 - floor),
              nrow(A), ncol(alpha))
  out <- matrix(0, nrow(A), ncol(alpha) + 1L)
  remaining <- rep(1, nrow(A))
  for (stick in seq_len(ncol(alpha))) {
    out[, stick] <- remaining * (1 - q[, stick])
    remaining <- remaining * q[, stick]
  }
  out[, ncol(out)] <- remaining
  out <- pmax(out, floor)
  out / rowSums(out)
}

.mcem_ref_soft_stick <- function(responsibility) {
  component_count <- ncol(responsibility)
  eligible <- success <- matrix(0, nrow(responsibility), component_count - 1L)
  for (stick in seq_len(component_count - 1L)) {
    eligible[, stick] <- rowSums(
      responsibility[, stick:component_count, drop = FALSE]
    )
    success[, stick] <- rowSums(
      responsibility[, (stick + 1L):component_count, drop = FALSE]
    )
  }
  list(eligible = eligible, success = success)
}

.mcem_ref_m_step <- function(A, responsibility, alpha_start,
                             intercept_mean, intercept_sd, slope_sd) {
  soft <- .mcem_ref_soft_stick(responsibility)
  alpha <- matrix(0, ncol(A), ncol(responsibility) - 1L)
  for (stick in seq_len(ncol(alpha))) {
    e <- soft$eligible[, stick]
    y <- soft$success[, stick]
    mean <- c(intercept_mean[stick], rep(0, ncol(A) - 1L))
    sd <- c(intercept_sd[stick], rep(slope_sd[stick], ncol(A) - 1L))
    fn <- function(parameter) {
      probability <- pmin(pmax(pnorm(drop(A %*% parameter)), 1e-12), 1 - 1e-12)
      -(
        sum(y * log(probability) + (e - y) * log1p(-probability)) -
          0.5 * sum(((parameter - mean) / sd)^2)
      )
    }
    gr <- function(parameter) {
      eta <- drop(A %*% parameter)
      probability <- pmin(pmax(pnorm(eta), 1e-12), 1 - 1e-12)
      d_eta <- dnorm(eta) * (y / probability - (e - y) / (1 - probability))
      -(drop(crossprod(A, d_eta)) - (parameter - mean) / sd^2)
    }
    alpha[, stick] <- optim(
      alpha_start[, stick], fn, gr, method = "BFGS",
      control = list(maxit = 500L, reltol = 1e-12)
    )$par
  }
  alpha
}

.mcem_write_csr_correlation <- function(correlation) {
  prefix <- tempfile("mcem_sbayesrc_csr_")
  marker_count <- nrow(correlation)
  edge <- which(upper.tri(correlation) & correlation != 0, arr.ind = TRUE)
  edge <- edge[order(edge[, 1L], edge[, 2L]), , drop = FALSE]
  row_ptr <- c(0, cumsum(tabulate(edge[, 1L], nbins = marker_count)))
  sblr:::.stblr_write_uint64_file(
    paste0(prefix, ".row_ptr.u64.bin"), row_ptr
  )
  sblr:::.stblr_write_uint32_file(
    paste0(prefix, ".col_idx.u32.0based.bin"), edge[, 2L] - 1L
  )
  writeBin(
    as.numeric(correlation[edge]), paste0(prefix, ".values.f32.bin"),
    size = 4L, endian = .Platform$endian
  )
  writeLines(c(
    "format=sparse_ld_csr", "storage=streamed_upper_triangle",
    "n_bed=NA", "n_used=NA", "n_samples_used=NA",
    paste0("n_variants=", marker_count), paste0("nnz=", nrow(edge)),
    "triangle=upper", "diagonal=implicit_1",
    paste0("row_ptr_file=", prefix, ".row_ptr.u64.bin"),
    paste0("col_idx_file=", prefix, ".col_idx.u32.0based.bin"),
    paste0("values_file=", prefix, ".values.f32.bin"),
    "row_ptr_type=uint64", "col_idx_type=uint32", "values_type=float32",
    "index_base=0", "value=r"
  ), paste0(prefix, ".meta.txt"))
  prefix
}

.mcem_write_bed <- function(path, dosage) {
 dosage_to_code <- c(`0` = 3L, `1` = 2L, `2` = 0L)
 packed <- unlist(lapply(seq_len(nrow(dosage)), function(marker) {
  codes <- unname(dosage_to_code[as.character(dosage[marker, ])])
  codes <- c(codes, rep(0L, (-length(codes)) %% 4L))
  vapply(seq(1L, length(codes), by = 4L), function(index) {
   sum(codes[index:(index + 3L)] * c(1L, 4L, 16L, 64L))
  }, integer(1))
 }))
 writeBin(as.raw(c(0x6c, 0x1b, 0x01, packed)), path)
}

.mcem_block_fixture <- function(seed = 9182L, marker_count = 12L,
                                sample_count = 80L) {
 set.seed(seed)
 allele_frequency <- seq(0.15, 0.42, length.out = marker_count)
 dosage <- vapply(allele_frequency, function(p) {
  stats::rbinom(sample_count, 2L, p)
 }, numeric(sample_count))
 dosage <- t(dosage)
 observed_frequency <- rowMeans(dosage) / 2
 scale <- sqrt(2 * observed_frequency * (1 - observed_frequency))
 X <- t((dosage - 2 * observed_frequency) / scale)
 annotation_signal <- as.numeric(scale(seq_len(marker_count)))
 A <- cbind(intercept = 1, annotation_signal = annotation_signal)
 gamma <- c(0, 0.05, 0.2)
 alpha_truth <- matrix(c(-0.7, -0.55, 0.35, -0.2), 2L, 2L)
 component_prior <- .mcem_ref_component_prior(A, alpha_truth)
 component <- apply(component_prior, 1L, function(probability) {
  sample.int(length(probability), 1L, prob = probability) - 1L
 })
 beta <- numeric(marker_count)
 active <- component > 0L
 beta[active] <- stats::rnorm(
  sum(active), sd = sqrt(gamma[component[active] + 1L] * 0.08)
 )
 y <- drop(X %*% beta + stats::rnorm(sample_count, sd = sqrt(0.7)))
 diagonal <- colSums(X^2)
 XtX <- crossprod(X)
 correlation <- XtX / sqrt(base::outer(diagonal, diagonal))
 diag(correlation) <- 1
 prefix <- .mcem_write_csr_correlation(correlation)
 bed_file <- tempfile("sbayesrc_em_", fileext = ".bed")
 .mcem_write_bed(bed_file, dosage)
 marker_names <- paste0("rs", seq_len(marker_count))
 sample_names <- paste0("id", seq_len(sample_count))
 stats <- list(
  wy = list(trait1 = stats::setNames(drop(crossprod(X, y)), marker_names)),
  ww = list(trait1 = stats::setNames(diagonal, marker_names)),
  yy = stats::setNames(sum(y^2), "trait1"),
  n = sample_count, m = marker_count,
  bed_files = bed_file, cls = list(seq_len(marker_count)),
  rows = seq_len(sample_count), af = list(observed_frequency),
  marker_names = marker_names, trait_names = "trait1"
 )
 Glist <- list(
  n = sample_count, ids = sample_names, bedfiles = bed_file,
  rsids = list(marker_names), rsidsLD = list(marker_names),
  chr = list(rep(1L, marker_count)), pos = list(seq_len(marker_count) * 100),
  af = list(observed_frequency), maf = list(pmin(observed_frequency,
                                                1 - observed_frequency))
 )
 phenotype_variance <- sum(y^2) / (sample_count - 1)
 list(
  stats = stats, Glist = Glist, A = A, gamma = gamma,
  alpha_truth = alpha_truth, prefix = prefix, bed_file = bed_file,
  B = matrix(0.08, 1L, 1L), E = matrix(phenotype_variance, 1L, 1L),
  ssb_prior = list(0.08), sse_prior = list(phenotype_variance),
  intercept_prior = rbind(
   type = c(0, 0), mean = c(-0.6, -0.5), precision = c(1, 1)
  )
 )
}

.mcem_cleanup_block_fixture <- function(fixture) {
 .mcem_cleanup_csr(fixture$prefix)
 unlink(fixture$bed_file)
 invisible(NULL)
}

.mcem_run_block <- function(fixture, alpha_start, seed,
                            updateB = FALSE, updateE = FALSE,
                            inner_sweeps = 600L, inner_burn = 220L,
                            max_outer = 25L) {
 sblr:::.stblr_mcem_sbayesrc_block_eigen(
  stats = fixture$stats, Glist = fixture$Glist,
  annotation = fixture$A, block_start = 1L,
  B = fixture$B, E = fixture$E,
  ssb_prior = fixture$ssb_prior, sse_prior = fixture$sse_prior,
  gamma = fixture$gamma, alpha_init = alpha_start,
  sigmaSqAlpha_init = c(1, 1),
  intercept_prior_resolved = fixture$intercept_prior,
  representation = "low_rank", eigen_prop = 0.999999,
  residual_policy = "gctb_block",
  block_ve_mode = if (isTRUE(updateE)) "allMixVe" else "fixVe",
  updateB = updateB, updateE = updateE,
  inner_sweeps = inner_sweeps, inner_burn = inner_burn,
  final_sweeps = inner_sweeps + 100L, final_burn = inner_burn,
  max_outer = max_outer, seed = seed, ncores = 1L
 )
}

.mcem_cleanup_csr <- function(prefix) {
  unlink(paste0(prefix, c(
    ".row_ptr.u64.bin", ".col_idx.u32.0based.bin", ".values.f32.bin",
    ".meta.txt"
  )))
  invisible(NULL)
}

.mcem_exact_fixture <- function(correlated = FALSE, marker_count = 6L,
                                seed = 20271280L) {
  set.seed(seed)
  sample_size <- 100L
  gamma <- c(0, 0.1, 1)
  correlation <- if (correlated) {
    toeplitz(0.45^(0:(marker_count - 1L)))
  } else diag(marker_count)
  annotation <- as.numeric(scale(rnorm(marker_count)))
  A <- cbind(Intercept = 1, annotation = annotation)
  baseline <- c(0.72, 0.2, 0.08)
  intercept_prior <- sblr:::.sbayesrc_resolve_intercept_prior(
    baseline, list(mean = "initial_mixture", sd = 1)
  )
  alpha_start <- matrix(0, 2L, 2L,
                        dimnames = list(colnames(A), c("stick_1", "stick_2")))
  alpha_start[1L, ] <- intercept_prior$mean
  score <- c(-1.4, 1.8, 0.7, -0.5, 2.1, -1.0)[seq_len(marker_count)] *
    sqrt(sample_size)
  prefix <- .mcem_write_csr_correlation(correlation)
  list(
    A = A, gamma = gamma, alpha_start = alpha_start,
    sigmaSqAlpha = rep(1, 2L), intercept_prior = intercept_prior,
    correlation = correlation, XtX = sample_size * correlation,
    score = score, sample_size = sample_size, sigma2_beta = 0.08,
    sigma2_e = 1, prefix = prefix,
    wy = list(score), ww = list(rep(sample_size, marker_count)),
    yy = sum(score^2) / sample_size + sample_size,
    b = list(rep(0, marker_count)), component = list(rep(0, marker_count)),
    r = list(score), B = matrix(0.08, 1L, 1L), E = matrix(1, 1L, 1L),
    ssb = list(0.08), sse = list(1)
  )
}

.mcem_log_sum_exp <- function(x) {
  maximum <- max(x)
  maximum + log(sum(exp(x - maximum)))
}

.mcem_allocation_space <- function(marker_count, component_count) {
  as.matrix(expand.grid(
    rep(list(0:(component_count - 1L)), marker_count),
    KEEP.OUT.ATTRS = FALSE
  ))
}

.mcem_configuration_log_bf <- function(component, fixture) {
  active <- which(component > 0L)
  if (!length(active)) return(0)
  variance <- fixture$sigma2_beta * fixture$gamma[component[active] + 1L]
  precision <- fixture$XtX[active, active, drop = FALSE] / fixture$sigma2_e +
    diag(1 / variance, length(active))
  h <- fixture$score[active] / fixture$sigma2_e
  R <- chol(precision)
  solved <- backsolve(R, forwardsolve(t(R), h))
  0.5 * (-2 * sum(log(diag(R))) - sum(log(variance)) + sum(h * solved))
}

.mcem_exact_target <- function(fixture) {
  configuration <- .mcem_allocation_space(
    nrow(fixture$A), length(fixture$gamma)
  )
  log_bf <- apply(configuration, 1L, .mcem_configuration_log_bf,
                  fixture = fixture)
  objective <- function(parameter) {
    alpha <- matrix(parameter, ncol(fixture$A), ncol(fixture$alpha_start))
    prior <- .mcem_ref_component_prior(fixture$A, alpha)
    log_weight <- log_bf
    for (marker in seq_len(nrow(fixture$A))) {
      log_weight <- log_weight + log(prior[cbind(
        rep.int(marker, nrow(configuration)), configuration[, marker] + 1L
      )])
    }
    intercept_mean <- fixture$intercept_prior$mean
    log_alpha_prior <- 0
    for (stick in seq_len(ncol(alpha))) {
      mean <- c(intercept_mean[stick], 0)
      log_alpha_prior <- log_alpha_prior - 0.5 * sum((alpha[, stick] - mean)^2)
    }
    .mcem_log_sum_exp(log_weight) + log_alpha_prior
  }
  fit <- optim(
    as.numeric(fixture$alpha_start), function(x) -objective(x),
    method = "BFGS", control = list(maxit = 1000L, reltol = 1e-12)
  )
  list(
    alpha = matrix(fit$par, ncol(fixture$A), ncol(fixture$alpha_start),
                   dimnames = dimnames(fixture$alpha_start)),
    objective = -fit$value, convergence = fit$convergence
  )
}

.mcem_run_fixture <- function(fixture, seed = 20271290L,
                              sweeps = 500L, burn = 150L,
                              outer = 12L) {
  sblr:::.stblr_mcem_sbayesrc_csr(
    fixture$wy, fixture$ww, fixture$yy, fixture$b, fixture$component,
    fixture$r, fixture$prefix, fixture$B, fixture$E, fixture$ssb,
    fixture$sse, fixture$A, fixture$gamma, fixture$alpha_start,
    fixture$sigmaSqAlpha, fixture$intercept_prior$native,
    fixture$sample_size, inner_sweeps = sweeps, inner_burn = burn,
    final_sweeps = sweeps, final_burn = burn, damping = 0.5,
    min_outer = 3L, max_outer = outer, ncores = 1L, seed = seed
  )
}
