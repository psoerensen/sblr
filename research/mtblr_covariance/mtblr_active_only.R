mtblr_active_only <- function(X, Y, Ve, pi, nu0, Psi0,
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
  delta <- integer(M)
  alpha <- matrix(0, M, T)
  residual <- dat$Y
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
        beta[j, ] <- 0
        alpha[j, ] <- 0
      }
      residual <- partial - tcrossprod(dat$X[, j], alpha[j, ])
    }

    if (isTRUE(update_vb)) {
      active <- which(delta == 1L)
      S <- if (length(active)) crossprod(beta[active, , drop = FALSE]) else matrix(0, T, T)
      Vb <- mt_rinvwishart(iw$nu + length(active), iw$Psi + S)
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
    sampler = "active_only",
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

mtblr_bayesr_active_only <- function(X, Y, Ve, component_prob, gamma,
                                     nu0, Psi0, n_iter = 10000L,
                                     burn = 2000L, seed = 1L,
                                     Vb_init = NULL, update_vb = TRUE,
                                     q = NULL) {
  dat <- mt_validate_data(X, Y, Ve)
  mt_validate_iw(nu0, Psi0, require_mean = TRUE)
  iw <- list(nu = nu0, Psi = Psi0)
  if (length(component_prob) != length(gamma) || gamma[1L] != 0 ||
      any(gamma[-1L] <= 0) || any(!is.finite(gamma)) ||
      any(component_prob <= 0) || abs(sum(component_prob) - 1) > 1e-12) {
    stop("BayesR requires probabilities summing to one and gamma = c(0, positive values).", call. = FALSE)
  }
  n_iter <- as.integer(n_iter)
  burn <- as.integer(burn)
  if (n_iter <= burn || burn < 0L) {
    stop("n_iter must be greater than burn >= 0.", call. = FALSE)
  }
  set.seed(seed)
  M <- dat$M
  T <- dat$T
  q <- mt_validate_marker_scale(q, M)
  state <- integer(M)
  theta <- alpha <- matrix(0, M, T)
  residual <- dat$Y
  if (is.null(Vb_init)) Vb_init <- mt_iw_mean(iw$nu, iw$Psi)
  mt_assert_spd(Vb_init, "Vb_init")
  Vb <- Vb_init
  keep <- n_iter - burn
  state_draws <- matrix(0L, keep, M)
  Vb_draws <- array(0, c(T, T, keep))
  alpha_draws <- array(0, c(keep, M, T))
  out_i <- 0L
  for (iter in seq_len(n_iter)) {
    for (j in seq_len(M)) {
      partial <- residual + tcrossprod(dat$X[, j], alpha[j, ])
      cond <- mt_component_conditional(
        partial, dat$X[, j], Vb, dat$Ve, component_prob, gamma, q[j]
      )
      state[j] <- sample.int(length(gamma), 1L, prob = cond$probability) - 1L
      if (state[j] > 0L) {
        k <- state[j] + 1L
        theta[j, ] <- mt_rmvnorm(cond$candidates[[k]]$mean,
                                 cond$candidates[[k]]$covariance)
        alpha[j, ] <- theta[j, ]
      } else {
        theta[j, ] <- alpha[j, ] <- 0
      }
      residual <- partial - tcrossprod(dat$X[, j], alpha[j, ])
    }
    if (isTRUE(update_vb)) {
      sufficient <- mt_covariance_sufficient(theta, state, gamma, q,
                                              storage = "scaled")
      Vb <- mt_rinvwishart(
        iw$nu + sum(state > 0L),
        iw$Psi + sufficient$statistic
      )
    }
    if (iter > burn) {
      out_i <- out_i + 1L
      state_draws[out_i, ] <- state
      Vb_draws[, , out_i] <- Vb
      alpha_draws[out_i, , ] <- alpha
    }
  }
  ans <- list(sampler = "bayesr_active_only", state = state_draws,
              theta = alpha_draws, alpha = alpha_draws, Vb = Vb_draws,
              gamma = gamma, q = q, effect_storage = "scaled_theta",
              component_prob = component_prob, seed = seed,
              update_vb = isTRUE(update_vb))
  mt_validate_draws(ans$Vb)
  ans
}
