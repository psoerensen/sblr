# Isolated standard SBayesRC alpha-recovery reference
#
# This standalone base-R experiment deliberately adapts the validated
# mathematics in ignored local input
# `local_reference/test_sbayesrc_alpha_recovery_v2.R` (SHA-256
# 9f1edd471d019aaf1e65a63be289abe47f5e162460c58ef5d4639730bcf082e3).
# The compact fixed-z reference was
# `local_reference/test_sbayesrc_continuous_reference.R` (SHA-256
# bba54addbec93a2219f01bd0269714fa2a0032a1327092e5ee58732637b41c16).
# Jian Zeng's 2024 R file is provenance only (SHA-256
# 50c393dd579459ce611fbc58f884d56a032531437b1d190411a5d7db842ab683).
#
# This is an isolated probit hierarchy. It does not read Study 06 fits, invoke
# an sblr production sampler, regenerate frozen truth, or implement SBayesRC-S.

set.seed(20260807)

N <- 1500L
K <- 3L
annotation_names <- c("enriched_binary", "continuous_signal", "null_annotation")
alpha_nonint_true <- rbind(
  enriched_binary = c(1.60, 0.30, 0.20),
  continuous_signal = c(0.30, 0.15, 0.10),
  null_annotation = c(0, 0, 0)
)
colnames(alpha_nonint_true) <- paste0("stick", seq_len(K))
target_continue <- c(171 / 1500, 87 / 171, 37 / 87)
tau2 <- rep(1, K)
n_chains <- 4L
n_iter <- 8000L
burn <- 2000L
n_rep <- 20L

rtruncnorm_probit <- function(mu, d, eps = 1e-12) {
  p0 <- pmin(pmax(pnorm(-mu), eps), 1 - eps)
  u <- numeric(length(mu))
  is1 <- d == 1L
  if (any(is1)) u[is1] <- runif(sum(is1), p0[is1], 1 - eps)
  if (any(!is1)) u[!is1] <- runif(sum(!is1), eps, p0[!is1])
  mu + qnorm(u)
}

rhat_basic <- function(x) {
  x <- as.matrix(x)
  n <- nrow(x)
  m <- ncol(x)
  if (m < 2L || n < 2L) return(NA_real_)
  W <- mean(apply(x, 2, var))
  if (!is.finite(W) || W <= 0) return(NA_real_)
  B <- n * var(colMeans(x))
  sqrt((((n - 1) / n) * W + B / n) / W)
}

calibrate_intercept <- function(A, slopes, target) {
  f <- function(a0) mean(pnorm(a0 + as.numeric(A %*% slopes))) - target
  uniroot(f, c(-12, 12), tol = 1e-12)$root
}

exact_alpha_posterior <- function(z, A, tau2) {
  X <- cbind(Intercept = 1, A)
  B <- diag(c(0, rep(1 / tau2, ncol(X) - 1L)))
  precision <- crossprod(X) + B
  covariance <- solve(precision)
  mean <- as.numeric(covariance %*% crossprod(X, z))
  names(mean) <- colnames(X)
  dimnames(covariance) <- list(colnames(X), colnames(X))
  list(mean = mean, cov = covariance, precision = precision)
}

blocked_alpha_draw <- function(z, A, tau2) {
  post <- exact_alpha_posterior(z, A, tau2)
  as.numeric(post$mean + backsolve(chol(post$precision), rnorm(length(post$mean))))
}

scalar_alpha_sweep <- function(z, A, alpha, tau2, intercept_n = nrow(A)) {
  X <- cbind(Intercept = 1, A)
  residual <- as.numeric(z - X %*% alpha)

  old <- alpha[1L]
  v <- 1 / intercept_n
  alpha[1L] <- rnorm(1L, v * (sum(residual) + intercept_n * old), sqrt(v))
  residual <- residual + old - alpha[1L]

  for (j in 2:ncol(X)) {
    x <- X[, j]
    old <- alpha[j]
    x2 <- sum(x * x)
    v <- 1 / (x2 + 1 / tau2)
    alpha[j] <- rnorm(1L, v * (sum(x * residual) + x2 * old), sqrt(v))
    residual <- residual + x * (old - alpha[j])
  }
  alpha
}

simulate_fixture <- function() {
  A <- cbind(
    enriched_binary = rbinom(N, 1L, 0.10),
    continuous_signal = as.numeric(scale(rnorm(N))),
    null_annotation = as.numeric(scale(rnorm(N)))
  )
  eligible <- d_list <- z_true_list <- vector("list", K)
  alpha_true <- matrix(NA_real_, 4L, K,
    dimnames = list(c("Intercept", annotation_names), colnames(alpha_nonint_true)))
  alpha_true[-1L, ] <- alpha_nonint_true
  idx <- seq_len(N)

  for (k in seq_len(K)) {
    Ak <- A[idx, , drop = FALSE]
    a0 <- calibrate_intercept(Ak, alpha_nonint_true[, k], target_continue[k])
    alpha_true[1L, k] <- a0
    z <- rnorm(length(idx), a0 + as.numeric(Ak %*% alpha_nonint_true[, k]), 1)
    d <- as.integer(z > 0)
    eligible[[k]] <- idx
    z_true_list[[k]] <- z
    d_list[[k]] <- d
    if (k < K) idx <- idx[d == 1L]
  }
  list(A = A, alpha_true = alpha_true, eligible = eligible,
    d_list = d_list, z_true_list = z_true_list)
}

fixed_z_chain <- function(z, A, tau2, iterations = 30000L, discard = 5000L) {
  alpha <- rep(0, ncol(A) + 1L)
  draws <- matrix(NA_real_, iterations - discard, length(alpha),
    dimnames = list(NULL, c("Intercept", colnames(A))))
  pos <- 0L
  for (iter in seq_len(iterations)) {
    alpha <- scalar_alpha_sweep(z, A, alpha, tau2)
    if (iter > discard) {
      pos <- pos + 1L
      draws[pos, ] <- alpha
    }
  }
  draws
}

run_fixed_z <- function(fixture) {
  rows <- lapply(seq_len(K), function(k) {
    idx <- fixture$eligible[[k]]
    A <- fixture$A[idx, , drop = FALSE]
    exact <- exact_alpha_posterior(fixture$z_true_list[[k]], A, tau2[k])
    draws <- fixed_z_chain(fixture$z_true_list[[k]], A, tau2[k])
    data.frame(
      stick = k,
      eligible = length(idx),
      enriched_eligible = sum(A[, "enriched_binary"] == 1),
      max_mean_error = max(abs(colMeans(draws) - exact$mean)),
      max_covariance_error = max(abs(cov(draws) - exact$cov))
    )
  })
  do.call(rbind, rows)
}

run_probit_chain <- function(A, d, tau2, update, intercept_contract,
                             global_m = nrow(A), init = NULL) {
  alpha <- if (is.null(init)) rnorm(ncol(A) + 1L, 0, 0.5) else as.numeric(init)
  kept <- matrix(NA_real_, n_iter - burn, length(alpha),
    dimnames = list(NULL, c("Intercept", colnames(A))))
  X <- cbind(Intercept = 1, A)
  pos <- 0L
  for (iter in seq_len(n_iter)) {
    z <- rtruncnorm_probit(as.numeric(X %*% alpha), d)
    if (update == "blocked") {
      alpha <- blocked_alpha_draw(z, A, tau2)
    } else {
      ni <- if (intercept_contract == "eligible") nrow(A) else global_m
      alpha <- scalar_alpha_sweep(z, A, alpha, tau2, ni)
    }
    if (iter > burn) {
      pos <- pos + 1L
      kept[pos, ] <- alpha
    }
  }
  kept
}

fit_multichain <- function(A, d, tau2, update,
                           intercept_contract = "eligible", global_m = nrow(A)) {
  chains <- lapply(seq_len(n_chains), function(chain) {
    run_probit_chain(A, d, tau2, update, intercept_contract, global_m,
      init = rnorm(ncol(A) + 1L, 0, 1))
  })
  pooled <- do.call(rbind, chains)
  out <- data.frame(
    parameter = colnames(pooled),
    mean = colMeans(pooled),
    sd = apply(pooled, 2, sd),
    lower = apply(pooled, 2, quantile, 0.025),
    upper = apply(pooled, 2, quantile, 0.975),
    rhat = NA_real_
  )
  for (j in seq_len(ncol(pooled))) {
    out$rhat[j] <- rhat_basic(do.call(cbind, lapply(chains, function(x) x[, j])))
  }
  out
}

run_single_fixture <- function(fixture, update) {
  do.call(rbind, lapply(seq_len(K), function(k) {
    idx <- fixture$eligible[[k]]
    A <- fixture$A[idx, , drop = FALSE]
    out <- fit_multichain(A, fixture$d_list[[k]], tau2[k], update)
    out$truth <- fixture$alpha_true[out$parameter, k]
    out$bias <- out$mean - out$truth
    out$covered <- out$lower <= out$truth & out$upper >= out$truth
    out$stick <- k
    out$update <- update
    out
  }))
}

run_intercept_contract_diagnostic <- function(fixture) {
  do.call(rbind, lapply(2:K, function(k) {
    idx <- fixture$eligible[[k]]
    A <- fixture$A[idx, , drop = FALSE]
    eligible <- fit_multichain(A, fixture$d_list[[k]], tau2[k], "scalar",
      "eligible", N)
    total_m <- fit_multichain(A, fixture$d_list[[k]], tau2[k], "scalar",
      "global_m", N)
    data.frame(
      stick = k,
      parameter = rownames(fixture$alpha_true),
      truth = fixture$alpha_true[, k],
      eligible_mean = eligible$mean,
      total_m_mean = total_m$mean,
      eligible_rhat = eligible$rhat,
      total_m_rhat = total_m$rhat
    )
  }))
}

run_repeated_recovery <- function() {
  records <- list()
  pos <- 0L
  for (replicate in seq_len(n_rep)) {
    fixture <- simulate_fixture()
    for (k in seq_len(K)) {
      idx <- fixture$eligible[[k]]
      A <- fixture$A[idx, , drop = FALSE]
      fit <- fit_multichain(A, fixture$d_list[[k]], tau2[k], "blocked")
      truth <- fixture$alpha_true[fit$parameter, k]
      for (j in seq_len(nrow(fit))) {
        pos <- pos + 1L
        records[[pos]] <- data.frame(
          replicate = replicate,
          stick = k,
          parameter = fit$parameter[j],
          eligible = length(idx),
          enriched_eligible = sum(A[, "enriched_binary"] == 1),
          truth = truth[j],
          post_mean = fit$mean[j],
          post_sd = fit$sd[j],
          lower = fit$lower[j],
          upper = fit$upper[j],
          rhat = fit$rhat[j],
          covered = fit$lower[j] <= truth[j] && fit$upper[j] >= truth[j]
        )
      }
    }
  }
  records <- do.call(rbind, records)
  records$bias <- records$post_mean - records$truth
  groups <- split(records, list(records$stick, records$parameter), drop = TRUE)
  summary <- do.call(rbind, lapply(groups, function(x) data.frame(
    stick = x$stick[1L],
    parameter = x$parameter[1L],
    mean_eligible = mean(x$eligible),
    mean_enriched_eligible = mean(x$enriched_eligible),
    truth_mean = mean(x$truth),
    posterior_mean = mean(x$post_mean),
    bias = mean(x$bias),
    rmse = sqrt(mean(x$bias^2)),
    mean_post_sd = mean(x$post_sd),
    coverage_95 = mean(x$covered),
    mean_rhat = mean(x$rhat),
    max_rhat = max(x$rhat)
  )))
  rownames(summary) <- NULL
  summary[order(summary$stick,
    match(summary$parameter, c("Intercept", annotation_names))), ]
}

output_dir <- "results/local/06_annotation_models/alpha_recovery_addendum"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

fixture <- simulate_fixture()
fixed_z <- run_fixed_z(fixture)
single_blocked <- run_single_fixture(fixture, "blocked")
single_scalar <- run_single_fixture(fixture, "scalar")
intercept_contract <- run_intercept_contract_diagnostic(fixture)
repeated <- run_repeated_recovery()

write.csv(fixed_z, file.path(output_dir, "fixed_z_summary.csv"), row.names = FALSE)
write.csv(rbind(single_blocked, single_scalar),
  file.path(output_dir, "single_fixture_summary.csv"), row.names = FALSE)
write.csv(intercept_contract,
  file.path(output_dir, "intercept_contract_summary.csv"), row.names = FALSE)
write.csv(repeated, file.path(output_dir, "alpha_recovery_summary.csv"), row.names = FALSE)

stopifnot(
  max(fixed_z$max_mean_error) <= 0.04,
  max(fixed_z$max_covariance_error) <= 0.05,
  max(repeated$max_rhat) <= 1.05
)

cat("ALPHA-R1 PASS — isolated standard SBayesRC alpha hierarchy validated\n")
cat("Compact output:", normalizePath(output_dir, winslash = "/"), "\n")
