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
