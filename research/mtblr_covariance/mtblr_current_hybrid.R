mt_current_covariance_heuristic <- function(alpha, delta, nu0, ssb_prior,
                                            correlation_shrink = 20,
                                            rho_limit = 0.95,
                                            eigen_floor = 1e-10) {
  T <- ncol(alpha)
  out <- matrix(0, T, T)
  for (t in seq_len(T)) {
    active <- delta != 0L & alpha[, t] != 0
    shape <- nu0 + sum(active)
    scale <- nu0 * ssb_prior[t, t] + sum(alpha[active, t]^2)
    out[t, t] <- scale / rchisq(1L, shape)
    out[t, t] <- max(out[t, t], eigen_floor)
  }
  if (T > 1L) {
    for (t in seq_len(T - 1L)) for (s in (t + 1L):T) {
      shared <- delta != 0L & alpha[, t] != 0 & alpha[, s] != 0
      n_shared <- sum(shared)
      empirical <- if (n_shared) sum(alpha[shared, t] * alpha[shared, s]) / n_shared else 0
      prior_weight <- n_shared / (n_shared + nu0)
      shrink_weight <- n_shared / (n_shared + correlation_shrink)
      cov_ts <- prior_weight * shrink_weight * empirical
      limit <- rho_limit * sqrt(out[t, t] * out[s, s])
      out[t, s] <- out[s, t] <- max(-limit, min(limit, cov_ts))
    }
  }
  eig <- eigen(mt_symmetrize(out), symmetric = TRUE)
  eig$values <- pmax(eig$values, eigen_floor)
  mt_symmetrize(eig$vectors %*% (eig$values * t(eig$vectors)))
}

mtblr_current_hybrid <- function(X, Y, Ve, pi, nu0, ssb_prior,
                                 n_iter = 10000L, burn = 2000L,
                                 seed = 1L) {
  dat <- mt_validate_data(X, Y, Ve)
  mt_validate_iw(nu0, ssb_prior, require_mean = TRUE)
  iw <- list(nu = nu0, Psi = ssb_prior)
  n_iter <- as.integer(n_iter)
  burn <- as.integer(burn)
  if (n_iter <= burn || burn < 0L) {
    stop("n_iter must be greater than burn >= 0.", call. = FALSE)
  }
  set.seed(seed)
  M <- dat$M
  T <- dat$T
  beta <- matrix(0, M, T)
  delta <- integer(M)
  alpha <- matrix(0, M, T)
  residual <- dat$Y
  Vb <- mt_iw_mean(iw$nu, iw$Psi)
  keep <- n_iter - burn
  delta_draws <- matrix(0L, keep, M)
  alpha_draws <- array(0, c(keep, M, T))
  fitted_draws <- array(0, c(keep, dat$N, T))
  Vb_latent_draws <- Vb_output_draws <- array(0, c(T, T, keep))
  out_i <- 0L

  for (iter in seq_len(n_iter)) {
    Vb <- mt_current_covariance_heuristic(alpha, delta, nu0, ssb_prior)
    Vb_latent <- mt_rinvwishart(nu0 + M, ssb_prior + crossprod(beta))

    for (j in seq_len(M)) {
      partial <- residual + tcrossprod(dat$X[, j], alpha[j, ])
      cond <- mt_marker_conditional(
        partial, dat$X[, j], Vb_latent, dat$Ve, pi
      )
      delta[j] <- rbinom(1L, 1L, cond$probability)
      if (delta[j] == 1L) {
        beta[j, ] <- mt_rmvnorm(cond$mean, cond$covariance)
        alpha[j, ] <- beta[j, ]
      } else {
        beta[j, ] <- mt_rmvnorm(rep(0, T), Vb_latent)
        alpha[j, ] <- 0
      }
      residual <- partial - tcrossprod(dat$X[, j], alpha[j, ])
    }

    Vb <- mt_current_covariance_heuristic(alpha, delta, nu0, ssb_prior)
    if (iter > burn) {
      out_i <- out_i + 1L
      delta_draws[out_i, ] <- delta
      alpha_draws[out_i, , ] <- alpha
      fitted_draws[out_i, , ] <- dat$X %*% alpha
      Vb_latent_draws[, , out_i] <- Vb_latent
      Vb_output_draws[, , out_i] <- Vb
    }
  }
  ans <- list(sampler = "current_hybrid", delta = delta_draws,
              alpha = alpha_draws, fitted = fitted_draws,
              Vb_latent = Vb_latent_draws,
              Vb = Vb_output_draws, seed = seed,
              note = "Vb drives the returned trace; Vb_latent drives marker updates.")
  mt_validate_draws(ans$Vb)
  mt_validate_draws(ans$Vb_latent)
  ans
}
