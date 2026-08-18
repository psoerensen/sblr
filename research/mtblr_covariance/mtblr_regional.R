mt_current_pattern_heuristic <- function(alpha, mask, nu0, ssb_prior,
                                         subset = seq_len(nrow(alpha))) {
  T <- ncol(alpha)
  out <- matrix(0, T, T)
  for (t in seq_len(T)) {
    active <- subset[mask[subset, t] == 1L]
    scale <- nu0 * ssb_prior[t, t] + sum(alpha[active, t]^2)
    out[t, t] <- max(scale / rchisq(1L, nu0 + length(active)), 1e-8)
  }
  if (T > 1L) for (t in seq_len(T - 1L)) for (s in (t + 1L):T) {
    shared <- subset[mask[subset, t] == 1L & mask[subset, s] == 1L]
    n_shared <- length(shared)
    empirical <- if (n_shared) sum(alpha[shared, t] * alpha[shared, s]) /
      n_shared else 0
    w <- n_shared / (n_shared + nu0)
    covariance <- w * empirical + (1 - w) * ssb_prior[t, s]
    rho <- (n_shared / (n_shared + 20)) * covariance /
      sqrt(out[t, t] * out[s, s])
    rho <- max(-0.95, min(0.95, rho))
    out[t, s] <- out[s, t] <- rho * sqrt(out[t, t] * out[s, s])
  }
  eig <- eigen(mt_symmetrize(out), symmetric = TRUE)
  eig$values <- pmax(eig$values, 1e-8)
  mt_symmetrize(eig$vectors %*% (eig$values * t(eig$vectors)))
}

mt_current_pattern_conditional <- function(partial, x, Vb, Ve, patterns, Pi) {
  patterns <- mt_validate_patterns(patterns, Pi)
  Bi <- solve(Vb)
  Ei <- solve(Ve)
  score_trait <- diag(Ei) * as.numeric(crossprod(x, partial))
  log_weight <- log(Pi)
  mean <- covariance <- vector("list", nrow(patterns))
  for (s in seq_len(nrow(patterns))) {
    rhs <- patterns[s, ] * score_trait
    precision <- Bi
    diag(precision) <- diag(precision) + patterns[s, ] * sum(x^2) * diag(Ei)
    covariance[[s]] <- mt_symmetrize(solve(precision))
    mean[[s]] <- as.numeric(covariance[[s]] %*% rhs)
    log_weight[s] <- log_weight[s] - 0.5 * mt_logdet_spd(precision) +
      0.5 * sum(rhs * mean[[s]])
  }
  probability <- exp(log_weight - max(log_weight))
  list(probability = probability / sum(probability), mean = mean,
       covariance = covariance)
}

mtblr_pattern_current_hybrid <- function(X, Y, Ve,
                                         patterns = mt_pattern_space(),
                                         Pi_init = rep(1 / nrow(patterns), nrow(patterns)),
                                         dirichlet_prior = rep(1, nrow(patterns)),
                                         sets = list(seq_len(ncol(X))),
                                         nu0, ssb_prior,
                                         n_iter = 5000L, burn = 1000L,
                                         seed = 1L) {
  dat <- mt_validate_data(X, Y, Ve)
  patterns <- mt_validate_patterns(patterns, Pi_init)
  flat <- unlist(sets, use.names = FALSE)
  if (!is.list(sets) || any(lengths(sets) == 0L) || anyDuplicated(flat) ||
      !identical(sort(as.integer(flat)), seq_len(dat$M))) {
    stop("sets must be a disjoint, complete, nonempty marker partition.", call. = FALSE)
  }
  mt_validate_iw(nu0, ssb_prior, require_mean = TRUE)
  set.seed(seed)
  T <- dat$T
  M <- dat$M
  state <- rep(which(rowSums(patterns) == 0L), M)
  beta <- alpha <- matrix(0, M, T)
  mask <- patterns[state, , drop = FALSE]
  residual <- dat$Y
  Pi <- Pi_init
  Vb <- mt_iw_mean(nu0, ssb_prior)
  keep <- n_iter - burn
  state_draws <- matrix(0L, keep, M)
  Vb_output <- Vb_latent <- array(0, c(T, T, keep))
  Pi_draws <- matrix(0, keep, nrow(patterns))
  out_i <- 0L
  for (iter in seq_len(n_iter)) {
    counts <- dirichlet_prior
    for (set in sets) {
      invisible(mt_current_pattern_heuristic(alpha, mask, nu0, ssb_prior, set))
      # This all-marker draw immediately overwrites the set heuristic in source.
      Vb <- mt_rinvwishart(nu0 + M, ssb_prior + crossprod(beta))
      for (j in set) {
        partial <- residual + tcrossprod(dat$X[, j], alpha[j, ])
        conditional <- mt_current_pattern_conditional(partial, dat$X[, j], Vb,
                                                      Ve, patterns, Pi)
        state[j] <- sample.int(nrow(patterns), 1L, prob = conditional$probability)
        counts[state[j]] <- counts[state[j]] + 1
        beta[j, ] <- mt_rmvnorm(conditional$mean[[state[j]]],
                                conditional$covariance[[state[j]]])
        alpha[j, ] <- patterns[state[j], ] * beta[j, ]
        mask[j, ] <- patterns[state[j], ]
        residual <- partial - tcrossprod(dat$X[, j], alpha[j, ])
      }
    }
    Pi <- mt_rdirichlet(counts)
    Vb_reported <- mt_current_pattern_heuristic(alpha, mask, nu0, ssb_prior)
    if (iter > burn) {
      out_i <- out_i + 1L
      state_draws[out_i, ] <- state
      Vb_latent[, , out_i] <- Vb
      Vb_output[, , out_i] <- Vb_reported
      Pi_draws[out_i, ] <- Pi
    }
    Vb <- Vb_reported
  }
  list(sampler = "current_pattern_hybrid", patterns = patterns,
       state = state_draws, Vb_latent = Vb_latent, Vb = Vb_output,
       Pi = Pi_draws, sets = sets, seed = seed)
}

mtblr_regional_sampler <- function(X, Y, Ve, region,
                                   region_levels = unique(region),
                                   patterns = mt_pattern_space(),
                                   dirichlet_prior = rep(1, nrow(patterns)),
                                   nu0, Psi0,
                                   n_iter = 8000L, burn = 1500L, seed = 1L,
                                   augmentation = c("completed_active", "full"),
                                   covariance_mode = c("regional", "global", "fixed"),
                                   Vb_fixed = NULL,
                                   regional_pi = TRUE,
                                   update_pi = TRUE,
                                   Pi_init = NULL,
                                   global_dirichlet_prior = NULL) {
  augmentation <- match.arg(augmentation)
  covariance_mode <- match.arg(covariance_mode)
  controls <- mt_validate_mcmc_controls(n_iter, burn, seed)
  n_iter <- controls$n_iter
  burn <- controls$burn
  seed <- controls$seed
  for (flag in c("regional_pi", "update_pi")) {
    value <- get(flag)
    if (length(value) != 1L || !is.logical(value) || is.na(value)) {
      stop(flag, " must be one non-missing logical value.", call. = FALSE)
    }
  }
  dat <- mt_validate_data(X, Y, Ve)
  patterns <- mt_validate_patterns(patterns)
  if (ncol(patterns) != dat$T) stop("patterns must have T columns.", call. = FALSE)
  if (!is.atomic(region_levels) || !length(region_levels) || anyNA(region_levels) ||
      anyDuplicated(region_levels) || length(region) != dat$M || anyNA(region) ||
      any(!region %in% region_levels)) {
    stop("region_levels must be unique, non-missing declared levels covering every marker.",
         call. = FALSE)
  }
  mt_validate_iw(nu0, Psi0, require_mean = TRUE)
  R <- length(region_levels)
  K <- nrow(patterns)
  T <- dat$T
  regional_dirichlet_prior <- NULL
  if (regional_pi) {
    if (is.numeric(dirichlet_prior) && is.null(dim(dirichlet_prior)) &&
        length(dirichlet_prior) == K) {
      dirichlet_prior <- matrix(dirichlet_prior, R, K, byrow = TRUE)
    }
    if (!is.numeric(dirichlet_prior) ||
        !identical(dim(dirichlet_prior), c(R, K)) ||
        any(!is.finite(dirichlet_prior)) || any(dirichlet_prior <= 0)) {
      stop("regional Dirichlet prior must be a finite, strictly positive R by K matrix.",
           call. = FALSE)
    }
    regional_dirichlet_prior <- dirichlet_prior
  } else {
    if (!is.numeric(global_dirichlet_prior) ||
        !is.null(dim(global_dirichlet_prior)) ||
        length(global_dirichlet_prior) != K ||
        any(!is.finite(global_dirichlet_prior)) ||
        any(global_dirichlet_prior <= 0)) {
      stop("shared-global Pi requires one finite, strictly positive global_dirichlet_prior vector of length K.",
           call. = FALSE)
    }
    global_dirichlet_prior <- as.numeric(global_dirichlet_prior)
  }
  if (covariance_mode == "fixed") {
    if (is.null(Vb_fixed)) stop("fixed mode requires Vb_fixed.", call. = FALSE)
    if (is.matrix(Vb_fixed)) Vb_fixed <- replicate(R, Vb_fixed, simplify = FALSE)
    if (!is.list(Vb_fixed) || length(Vb_fixed) != R) {
      stop("Vb_fixed must provide one matrix per region.", call. = FALSE)
    }
    for (r in seq_len(R)) {
      if (!identical(dim(as.matrix(Vb_fixed[[r]])), c(T, T))) {
        stop("each Vb_fixed matrix must be T by T.", call. = FALSE)
      }
      mt_assert_spd(Vb_fixed[[r]], paste0("Vb_fixed[[", r, "]]"))
    }
  }
  set.seed(seed)
  region_index <- match(region, region_levels)
  V0 <- mt_iw_mean(nu0, Psi0)
  Vb <- if (covariance_mode == "fixed") Vb_fixed else replicate(R, V0, simplify = FALSE)
  if (is.null(Pi_init)) {
    Pi <- if (regional_pi) {
      regional_dirichlet_prior / rowSums(regional_dirichlet_prior)
    } else {
      matrix(global_dirichlet_prior / sum(global_dirichlet_prior), R, K,
             byrow = TRUE)
    }
  } else {
    if (!regional_pi && is.numeric(Pi_init) && length(Pi_init) == K) {
      Pi_init <- matrix(Pi_init, R, K, byrow = TRUE)
    }
    Pi <- as.matrix(Pi_init)
  }
  if (!identical(dim(Pi), c(R, K)) || any(!is.finite(Pi)) || any(Pi <= 0) ||
      any(abs(rowSums(Pi) - 1) > 1e-10)) stop("Pi_init must be positive R by K probabilities.")
  if (!regional_pi &&
      any(abs(Pi - matrix(Pi[1L, ], R, K, byrow = TRUE)) > 1e-12)) {
    stop("shared-global Pi requires identical Pi_init rows.", call. = FALSE)
  }
  null <- which(rowSums(patterns) == 0L)
  state <- rep(null, dat$M)
  beta <- alpha <- matrix(0, dat$M, T)
  if (augmentation == "full") for (j in seq_len(dat$M)) {
    beta[j, ] <- mt_rmvnorm(rep(0, T), Vb[[region_index[j]]])
  }
  residual <- dat$Y
  keep <- n_iter - burn
  state_draws <- matrix(0L, keep, dat$M)
  Vb_draws <- array(0, c(T, T, R, keep))
  Pi_draws <- array(0, c(keep, R, K))
  out_i <- 0L
  for (iter in seq_len(n_iter)) {
    for (j in seq_len(dat$M)) {
      r <- region_index[j]
      partial <- residual + tcrossprod(dat$X[, j], alpha[j, ])
      conditional <- mt_pattern_conditional(partial, dat$X[, j], Vb[[r]], Ve,
                                            patterns, Pi[r, ])
      state[j] <- mt_draw_pattern_joint(conditional)
      active <- which(patterns[state[j], ] == 1L)
      if (length(active)) {
        value <- mt_rmvnorm(conditional$active_mean[[state[j]]],
                            conditional$active_covariance[[state[j]]])
        beta[j, ] <- mt_complete_latent(patterns[state[j], ], value, Vb[[r]])
        alpha[j, ] <- patterns[state[j], ] * beta[j, ]
      } else {
        alpha[j, ] <- 0
        beta[j, ] <- if (augmentation == "full")
          mt_rmvnorm(rep(0, T), Vb[[r]]) else numeric(T)
      }
      residual <- partial - tcrossprod(dat$X[, j], alpha[j, ])
    }
    if (isTRUE(update_pi) && regional_pi) {
      for (r in seq_len(R)) {
        idx <- which(region_index == r)
        Pi[r, ] <- mt_rdirichlet(regional_dirichlet_prior[r, ] +
                                   tabulate(state[idx], nbins = K))
      }
    } else if (isTRUE(update_pi)) {
      common <- mt_rdirichlet(global_dirichlet_prior + tabulate(state, nbins = K))
      Pi <- matrix(common, R, K, byrow = TRUE)
    }
    if (covariance_mode == "regional") for (r in seq_len(R)) {
      idx <- which(region_index == r)
      included <- if (augmentation == "full") idx else idx[state[idx] != null]
      S <- if (length(included)) crossprod(beta[included, , drop = FALSE]) else
        matrix(0, T, T)
      Vb[[r]] <- mt_rinvwishart(nu0 + length(included), Psi0 + S)
    } else if (covariance_mode == "global") {
      included <- if (augmentation == "full") seq_len(dat$M) else which(state != null)
      S <- if (length(included)) crossprod(beta[included, , drop = FALSE]) else
        matrix(0, T, T)
      common <- mt_rinvwishart(nu0 + length(included), Psi0 + S)
      Vb <- replicate(R, common, simplify = FALSE)
    }
    if (iter > burn) {
      out_i <- out_i + 1L
      state_draws[out_i, ] <- state
      for (r in seq_len(R)) Vb_draws[, , r, out_i] <- Vb[[r]]
      Pi_draws[out_i, , ] <- Pi
    }
  }
  for (r in seq_len(R)) {
    regional_draws <- array(Vb_draws[, , r, , drop = FALSE], c(T, T, keep))
    mt_validate_draws(regional_draws)
  }
  list(sampler = paste("regional", augmentation, covariance_mode, sep = "_"),
       patterns = patterns, region = region, region_levels = region_levels,
       state = state_draws, Vb = Vb_draws, Pi = Pi_draws,
       alpha = alpha, beta = beta, seed = seed,
       regional_pi = regional_pi,
       regional_dirichlet_prior = regional_dirichlet_prior,
       global_dirichlet_prior = global_dirichlet_prior)
}
