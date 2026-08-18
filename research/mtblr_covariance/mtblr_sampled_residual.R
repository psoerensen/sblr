# Independent Phase 4b residual-covariance qualification reference.

mt_residual_covariance_posterior <- function(residual, df, scale) {
  residual <- as.matrix(residual)
  if (!is.numeric(residual) || nrow(residual) < 1L || ncol(residual) != 2L ||
      any(!is.finite(residual))) {
    stop("residual must be a finite N x 2 matrix.", call. = FALSE)
  }
  mt_validate_iw(df, scale)
  statistic <- crossprod(residual)
  list(
    degrees_of_freedom = df + nrow(residual),
    scale = mt_symmetrize(scale + statistic),
    statistic = mt_symmetrize(statistic))
}

mt_rsample_residual_covariance <- function(residual, df, scale, draws) {
  posterior <- mt_residual_covariance_posterior(residual, df, scale)
  draws <- as.integer(draws)
  if (length(draws) != 1L || is.na(draws) || draws < 1L) {
    stop("draws must be one positive integer.", call. = FALSE)
  }
  value <- array(0, c(2L, 2L, draws))
  for (draw in seq_len(draws)) {
    value[, , draw] <- mt_rinvwishart(
      posterior$degrees_of_freedom, posterior$scale)
  }
  c(posterior, list(draws = value))
}

mtblr_pattern_sampler_sampled_residual <- function(
    X, Y, patterns = mt_pattern_space(),
    dirichlet_prior = rep(1, nrow(patterns)),
    marker_df, marker_scale, residual_df, residual_scale,
    Vb_init, Ve_init, Pi_init = dirichlet_prior / sum(dirichlet_prior),
    n_iter = 4000L, burn = 1000L, seed = 1L,
    update_vb = TRUE, update_pi = TRUE, update_ve = TRUE) {
  dat <- mt_validate_data(X, Y, Ve_init)
  patterns <- mt_validate_patterns(patterns, Pi_init)
  if (ncol(patterns) != 2L || length(dirichlet_prior) != nrow(patterns) ||
      any(!is.finite(dirichlet_prior)) || any(dirichlet_prior <= 0)) {
    stop("The sampled-residual reference requires valid two-trait patterns.",
         call. = FALSE)
  }
  mt_validate_iw(marker_df, marker_scale)
  mt_validate_iw(residual_df, residual_scale)
  mt_assert_spd(Vb_init, "Vb_init")
  n_iter <- as.integer(n_iter); burn <- as.integer(burn)
  if (n_iter <= burn || burn < 0L) {
    stop("n_iter must exceed burn >= 0.", call. = FALSE)
  }

  set.seed(seed)
  null_index <- which(rowSums(patterns) == 0L)
  if (length(null_index) != 1L) stop("Exactly one null pattern is required.")
  state <- rep(null_index, dat$M)
  latent <- matrix(NA_real_, dat$M, 2L)
  realised <- matrix(0, dat$M, 2L)
  residual <- dat$Y
  Vb <- Vb_init; Ve <- Ve_init; Pi <- Pi_init
  keep <- n_iter - burn
  state_draws <- matrix(0L, keep, dat$M)
  realised_draws <- array(0, c(keep, dat$M, 2L))
  Vb_draws <- Ve_draws <- array(0, c(2L, 2L, keep))
  Pi_draws <- matrix(0, keep, nrow(patterns))
  fitted_draws <- array(0, c(keep, dat$N, 2L))
  out_index <- 0L

  for (iteration in seq_len(n_iter)) {
    for (marker in seq_len(dat$M)) {
      partial <- residual + tcrossprod(dat$X[, marker], realised[marker, ])
      conditional <- mt_pattern_conditional(
        partial, dat$X[, marker], Vb, Ve, patterns, Pi)
      state[marker] <- mt_draw_pattern_joint(conditional)
      active <- which(patterns[state[marker], ] == 1L)
      if (!length(active)) {
        latent[marker, ] <- NA_real_
        realised[marker, ] <- 0
      } else {
        active_value <- mt_rmvnorm(
          conditional$active_mean[[state[marker]]],
          conditional$active_covariance[[state[marker]]])
        latent[marker, ] <- mt_complete_latent(
          patterns[state[marker], ], active_value, Vb)
        realised[marker, ] <- patterns[state[marker], ] * latent[marker, ]
      }
      residual <- partial - tcrossprod(dat$X[, marker], realised[marker, ])
    }

    counts <- tabulate(state, nbins = nrow(patterns))
    if (isTRUE(update_pi)) Pi <- mt_rdirichlet(dirichlet_prior + counts)
    included <- which(state != null_index)
    marker_statistic <- if (length(included)) {
      crossprod(latent[included, , drop = FALSE])
    } else matrix(0, 2L, 2L)
    if (isTRUE(update_vb)) {
      Vb <- mt_rinvwishart(
        marker_df + length(included), marker_scale + marker_statistic)
    }
    if (isTRUE(update_ve)) {
      residual_posterior <- mt_residual_covariance_posterior(
        residual, residual_df, residual_scale)
      Ve <- mt_rinvwishart(
        residual_posterior$degrees_of_freedom,
        residual_posterior$scale)
    }

    if (iteration > burn) {
      out_index <- out_index + 1L
      state_draws[out_index, ] <- state
      realised_draws[out_index, , ] <- realised
      Vb_draws[, , out_index] <- Vb
      Ve_draws[, , out_index] <- Ve
      Pi_draws[out_index, ] <- Pi
      fitted_draws[out_index, , ] <- dat$Y - residual
    }
  }
  list(
    patterns = patterns, state = state_draws, alpha = realised_draws,
    Vb = Vb_draws, Ve = Ve_draws, Pi = Pi_draws, fitted = fitted_draws,
    final = list(state = state, latent = latent, realised = realised,
                 residual = residual, Vb = Vb, Ve = Ve, Pi = Pi),
    update_vb = update_vb, update_pi = update_pi, update_ve = update_ve)
}
