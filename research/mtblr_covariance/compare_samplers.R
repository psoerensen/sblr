script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
research_dir <- if (length(script_arg)) {
  dirname(normalizePath(sub("^--file=", "", script_arg[1L]), mustWork = TRUE))
} else {
  normalizePath("research/mtblr_covariance", mustWork = TRUE)
}

for (file in c("mtblr_reference_model.R", "mtblr_exact_reference.R",
               "mtblr_full_latent.R", "mtblr_active_only.R",
               "mtblr_current_hybrid.R")) {
  source(file.path(research_dir, file), local = FALSE)
}

fixture <- function(M = 2L) {
  X <- cbind(
    c(-1.3, -0.8, -0.2, 0.4, 0.9, 1.0),
    c(0.7, -1.1, 0.5, -0.6, 1.2, -0.7)
  )[, seq_len(M), drop = FALSE]
  alpha <- matrix(c(0.45, -0.25, -0.18, 0.35), 2L, 2L, byrow = TRUE)[seq_len(M), , drop = FALSE]
  Y <- X %*% alpha + matrix(c(
    0.10, -0.05, -0.08, 0.03, 0.04, 0.08,
    -0.06, 0.09, 0.03, -0.07, 0.05, -0.04
  ), 6L, 2L)
  list(X = X, Y = Y, Ve = matrix(c(0.35, 0.07, 0.07, 0.28), 2L),
       Vb = matrix(c(0.32, 0.10, 0.10, 0.26), 2L), alpha = alpha)
}

summarize_sampler <- function(fit) {
  vb11 <- fit$Vb[1L, 1L, ]
  vb12 <- fit$Vb[1L, 2L, ]
  c(pip1 = mean(fit$delta[, 1L]),
    pip2 = if (ncol(fit$delta) > 1L) mean(fit$delta[, 2L]) else NA_real_,
    alpha11 = mean(fit$alpha[, 1L, 1L]),
    alpha12 = mean(fit$alpha[, 1L, 2L]),
    prediction11 = mean(fit$fitted[, 1L, 1L]),
    Vb11 = mean(vb11), Vb12 = mean(vb12),
    Vb11_lag1 = mt_chain_diagnostics(vb11)$lag1,
    Vb11_ess = mt_chain_diagnostics(vb11)$ess)
}

run_sparse_inactive_comparison <- function(marker_counts = c(2L, 12L),
                                           pi = 0.03,
                                           n_iter = 5000L,
                                           burn = 1000L,
                                           seed = 931L) {
  Ve <- matrix(c(0.4, 0.05, 0.05, 0.35), 2L)
  nu0 <- 7
  prior_mean <- matrix(c(0.30, 0.06, 0.06, 0.24), 2L)
  Psi0 <- mt_iw_scale_from_mean(prior_mean, nu0)
  rows <- list()
  q <- 0L
  for (M in marker_counts) {
    X <- matrix(0, 4L, M)
    Y <- matrix(c(0.1, -0.1, 0.05, -0.05,
                  -0.05, 0.08, -0.04, 0.02), 4L, 2L)
    for (kind in c("full", "active")) {
      fit <- if (kind == "full") {
        mtblr_full_latent(X, Y, Ve, pi, nu0, Psi0,
                          n_iter, burn, seed + M)
      } else {
        mtblr_active_only(X, Y, Ve, pi, nu0, Psi0,
                          n_iter, burn, seed + M + 100L)
      }
      diagnostic <- mt_chain_diagnostics(fit$Vb[1L, 1L, ])
      q <- q + 1L
      rows[[q]] <- data.frame(
        sampler = kind, markers = M,
        mean_active = mean(rowSums(fit$delta)),
        Vb11 = mean(fit$Vb[1L, 1L, ]),
        Vb11_lag1 = diagnostic$lag1,
        Vb11_ess = diagnostic$ess
      )
    }
  }
  do.call(rbind, rows)
}

run_comparison <- function(n_iter = 12000L, burn = 2000L,
                           seeds = c(2901L, 2902L, 2903L),
                           grid_n = 15L) {
  f <- fixture(2L)
  pi <- 0.18
  nu0 <- 7
  prior_mean <- matrix(c(0.30, 0.06, 0.06, 0.24), 2L)
  Psi0 <- mt_iw_scale_from_mean(prior_mean, nu0)
  exact_known <- mt_exact_bayesc_known_vb(f$X, f$Y, f$Ve, f$Vb, pi)
  grid_fixture <- fixture(1L)
  grid <- mt_vb_grid_reference(grid_fixture$X, grid_fixture$Y,
                               grid_fixture$Ve, pi, nu0, Psi0,
                               grid_n = grid_n)

  rows <- list()
  fits <- list()
  q <- 0L
  for (seed in seeds) {
    for (kind in c("full", "active", "hybrid")) {
      q <- q + 1L
      elapsed <- system.time({
        fit <- switch(
          kind,
          full = mtblr_full_latent(f$X, f$Y, f$Ve, pi, nu0, Psi0,
                                   n_iter, burn, seed),
          active = mtblr_active_only(f$X, f$Y, f$Ve, pi, nu0, Psi0,
                                     n_iter, burn, seed),
          hybrid = mtblr_current_hybrid(f$X, f$Y, f$Ve, pi, nu0, Psi0,
                                        n_iter, burn, seed)
        )
      })[["elapsed"]]
      rows[[q]] <- c(sampler = kind, seed = seed,
                     summarize_sampler(fit), seconds = elapsed)
      fits[[paste(kind, seed, sep = "_")]] <- fit
    }
  }
  rows <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE)
  numeric_columns <- setdiff(names(rows), "sampler")
  rows[numeric_columns] <- lapply(rows[numeric_columns], as.numeric)

  known_full <- mtblr_full_latent(f$X, f$Y, f$Ve, pi, nu0, Psi0,
                                  16000L, 3000L, 813L,
                                  Vb_init = f$Vb, update_vb = FALSE)
  known_active <- mtblr_active_only(f$X, f$Y, f$Ve, pi, nu0, Psi0,
                                    16000L, 3000L, 814L,
                                    Vb_init = f$Vb, update_vb = FALSE)

  list(
    configuration = list(pi = pi, nu0 = nu0, Psi0 = Psi0,
                         n_iter = n_iter, burn = burn, seeds = seeds),
    exact_known_vb = exact_known,
    known_vb_monte_carlo = rbind(
      exact = c(exact_known$pip, exact_known$alpha_mean[1L, ]),
      full_latent = c(colMeans(known_full$delta),
                      apply(known_full$alpha[, 1L, , drop = FALSE], 3L, mean)),
      active_only = c(colMeans(known_active$delta),
                      apply(known_active$alpha[, 1L, , drop = FALSE], 3L, mean))
    ),
    Vb_grid = grid,
    chain_summary = rows,
    sampler_mean = aggregate(. ~ sampler, rows[, setdiff(names(rows), "seed")], mean),
    sparse_inactive = run_sparse_inactive_comparison(),
    fits = fits
  )
}

if (sys.nframe() == 0L) {
  result <- run_comparison()
  cat("Known-Vb exact and Monte Carlo comparison\n")
  print(signif(result$known_vb_monte_carlo, 5))
  cat("\nNumerical grid posterior mean of Vb\n")
  print(signif(result$Vb_grid$mean, 5))
  cat("Grid boundary-mass upper bound:",
      signif(result$Vb_grid$boundary_mass_upper_bound, 4), "\n")
  cat("\nAcross-chain sampler summaries\n")
  print(result$sampler_mean, digits = 5, row.names = FALSE)
  cat("\nSparse-inclusion/inactive-marker sensitivity\n")
  print(result$sparse_inactive, digits = 5, row.names = FALSE)
}
