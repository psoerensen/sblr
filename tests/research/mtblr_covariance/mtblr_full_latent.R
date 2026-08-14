mtblr_full_latent <- function(X, Y, Ve, pi, nu0, Psi0,
                              n_iter = 10000L, burn = 2000L,
                              seed = 1L, Vb_init = NULL,
                              update_vb = TRUE) {
  dat <- mt_validate_data(X, Y, Ve)
  mt_validate_iw(nu0, Psi0, require_mean = TRUE)
  iw <- list(nu = nu0, Psi = Psi0)
  if (length(pi) != 1L || !is.finite(pi) || pi <= 0 || pi >= 1) {
    stop("pi must be one finite number in (0, 1).", call. = FALSE)
  }
  n_iter <- as.integer(n_iter)
  burn <- as.integer(burn)
  if (n_iter <= burn || burn < 0L) {
    stop("n_iter must be greater than burn >= 0.", call. = FALSE)
  }
  if (is.null(Vb_init)) Vb_init <- mt_iw_mean(iw$nu, iw$Psi)
  mt_assert_spd(Vb_init, "Vb_init")

  set.seed(seed)
  M <- dat$M
  T <- dat$T
  beta <- matrix(0, M, T)
  for (j in seq_len(M)) beta[j, ] <- mt_rmvnorm(rep(0, T), Vb_init)
  delta <- integer(M)
  alpha <- beta * delta
  residual <- dat$Y - dat$X %*% alpha
  Vb <- Vb_init

  keep <- n_iter - burn
  delta_draws <- matrix(0L, keep, M)
  beta_draws <- array(0, c(keep, M, T))
  alpha_draws <- array(0, c(keep, M, T))
  Vb_draws <- array(0, c(T, T, keep))
  fitted_draws <- array(0, c(keep, dat$N, T))

  out_i <- 0L
  for (iter in seq_len(n_iter)) {
    for (j in seq_len(M)) {
      partial <- residual + tcrossprod(dat$X[, j], alpha[j, ])
      cond <- mt_marker_conditional(
        partial, dat$X[, j], Vb, dat$Ve, pi
      )
      delta[j] <- rbinom(1L, 1L, cond$probability)
      if (delta[j] == 1L) {
        beta[j, ] <- mt_rmvnorm(cond$mean, cond$covariance)
        alpha[j, ] <- beta[j, ]
      } else {
        beta[j, ] <- mt_rmvnorm(rep(0, T), Vb)
        alpha[j, ] <- 0
      }
      residual <- partial - tcrossprod(dat$X[, j], alpha[j, ])
    }

    if (isTRUE(update_vb)) {
      Vb <- mt_rinvwishart(iw$nu + M, iw$Psi + crossprod(beta))
    }

    if (iter > burn) {
      out_i <- out_i + 1L
      delta_draws[out_i, ] <- delta
      beta_draws[out_i, , ] <- beta
      alpha_draws[out_i, , ] <- alpha
      Vb_draws[, , out_i] <- Vb
      fitted_draws[out_i, , ] <- dat$X %*% alpha
    }
  }

  ans <- list(
    sampler = "full_latent",
    delta = delta_draws,
    beta = beta_draws,
    alpha = alpha_draws,
    Vb = Vb_draws,
    fitted = fitted_draws,
    pi = pi,
    prior = iw,
    update_vb = isTRUE(update_vb),
    seed = seed
  )
  mt_validate_draws(ans$Vb)
  ans
}
