mt_pattern_index <- function(pattern, patterns) {
  hit <- which(rowSums(abs(patterns - matrix(pattern, nrow(patterns),
                                             ncol(patterns), byrow = TRUE))) == 0L)
  if (length(hit) != 1L) stop("pattern is absent or duplicated.", call. = FALSE)
  hit
}

mt_draw_pattern_joint <- function(conditional) {
  sample.int(length(conditional$probability), 1L,
             prob = conditional$probability)
}

mt_pattern_graph_connected <- function(patterns) {
  patterns <- mt_validate_patterns(patterns)
  if (nrow(patterns) == 1L) return(TRUE)
  adjacency <- as.matrix(stats::dist(patterns, method = "manhattan")) == 1
  visited <- rep(FALSE, nrow(patterns))
  frontier <- 1L
  while (length(frontier)) {
    node <- frontier[1L]
    frontier <- frontier[-1L]
    if (visited[node]) next
    visited[node] <- TRUE
    frontier <- c(frontier, which(adjacency[node, ] & !visited))
  }
  all(visited)
}

mt_draw_pattern_coordinate <- function(current_index, conditional, patterns) {
  current <- patterns[current_index, ]
  for (trait in seq_len(ncol(patterns))) {
    candidates <- which(vapply(seq_len(nrow(patterns)), function(s) {
      all(patterns[s, -trait, drop = FALSE] == current[-trait])
    }, logical(1L)))
    probability <- conditional$probability[candidates]
    probability <- probability / sum(probability)
    current_index <- sample(candidates, 1L, prob = probability)
    current <- patterns[current_index, ]
  }
  current_index
}

mtblr_pattern_sampler <- function(X, Y, Ve, patterns = mt_pattern_space(),
                                  dirichlet_prior = rep(1, nrow(patterns)),
                                  nu0, Psi0, n_iter = 10000L, burn = 2000L,
                                  seed = 1L,
                                  augmentation = c("full", "completed_active"),
                                  state_transition = c("joint", "coordinate"),
                                  Vb_init = NULL, Pi_init = NULL,
                                  update_vb = TRUE, update_pi = TRUE) {
  augmentation <- match.arg(augmentation)
  state_transition <- match.arg(state_transition)
  dat <- mt_validate_data(X, Y, Ve)
  patterns <- mt_validate_patterns(patterns)
  if (ncol(patterns) != dat$T) stop("patterns must have T columns.", call. = FALSE)
  if (state_transition == "coordinate" && !mt_pattern_graph_connected(patterns)) {
    stop("coordinate updating requires a pattern graph connected by one-coordinate moves; use joint updating for this restricted pattern set.",
         call. = FALSE)
  }
  if (length(dirichlet_prior) != nrow(patterns) ||
      any(!is.finite(dirichlet_prior)) || any(dirichlet_prior <= 0)) {
    stop("Dirichlet prior must be positive and match patterns.", call. = FALSE)
  }
  mt_validate_iw(nu0, Psi0, require_mean = TRUE)
  n_iter <- as.integer(n_iter)
  burn <- as.integer(burn)
  if (n_iter <= burn || burn < 0L) stop("n_iter must exceed burn >= 0.", call. = FALSE)
  if (is.null(Vb_init)) Vb_init <- mt_iw_mean(nu0, Psi0)
  mt_assert_spd(Vb_init, "Vb_init")
  if (is.null(Pi_init)) Pi_init <- dirichlet_prior / sum(dirichlet_prior)
  mt_validate_patterns(patterns, Pi_init)

  set.seed(seed)
  M <- dat$M
  T <- dat$T
  null_index <- which(rowSums(patterns) == 0L)
  state <- rep(null_index, M)
  beta <- alpha <- matrix(0, M, T)
  if (augmentation == "full") {
    for (j in seq_len(M)) beta[j, ] <- mt_rmvnorm(rep(0, T), Vb_init)
  }
  residual <- dat$Y
  Vb <- Vb_init
  Pi <- Pi_init
  keep <- n_iter - burn
  state_draws <- matrix(0L, keep, M)
  beta_draws <- alpha_draws <- array(0, c(keep, M, T))
  Vb_draws <- array(0, c(T, T, keep))
  Pi_draws <- matrix(0, keep, nrow(patterns))
  fitted_draws <- array(0, c(keep, dat$N, T))
  transitions <- matrix(0L, nrow(patterns), nrow(patterns),
                        dimnames = list(rownames(patterns), rownames(patterns)))
  out_i <- 0L

  for (iter in seq_len(n_iter)) {
    for (j in seq_len(M)) {
      previous <- state[j]
      partial <- residual + tcrossprod(dat$X[, j], alpha[j, ])
      conditional <- mt_pattern_conditional(partial, dat$X[, j], Vb, Ve,
                                            patterns, Pi)
      state[j] <- if (state_transition == "joint") {
        mt_draw_pattern_joint(conditional)
      } else {
        mt_draw_pattern_coordinate(state[j], conditional, patterns)
      }
      transitions[previous, state[j]] <- transitions[previous, state[j]] + 1L
      active <- which(patterns[state[j], ] == 1L)
      if (!length(active)) {
        alpha[j, ] <- 0
        beta[j, ] <- if (augmentation == "full") {
          mt_rmvnorm(rep(0, T), Vb)
        } else numeric(T)
      } else {
        active_value <- mt_rmvnorm(conditional$active_mean[[state[j]]],
                                   conditional$active_covariance[[state[j]]])
        beta[j, ] <- mt_complete_latent(patterns[state[j], ], active_value, Vb)
        alpha[j, ] <- patterns[state[j], ] * beta[j, ]
      }
      residual <- partial - tcrossprod(dat$X[, j], alpha[j, ])
    }

    counts <- tabulate(state, nbins = nrow(patterns))
    if (isTRUE(update_pi)) Pi <- mt_rdirichlet(dirichlet_prior + counts)
    if (isTRUE(update_vb)) {
      included <- if (augmentation == "full") seq_len(M) else
        which(state != null_index)
      statistic <- if (length(included)) {
        crossprod(beta[included, , drop = FALSE])
      } else matrix(0, T, T)
      Vb <- mt_rinvwishart(nu0 + length(included), Psi0 + statistic)
    }

    if (iter > burn) {
      out_i <- out_i + 1L
      state_draws[out_i, ] <- state
      beta_draws[out_i, , ] <- beta
      alpha_draws[out_i, , ] <- alpha
      Vb_draws[, , out_i] <- Vb
      Pi_draws[out_i, ] <- Pi
      fitted_draws[out_i, , ] <- dat$X %*% alpha
    }
  }
  ans <- list(sampler = paste(augmentation, state_transition, sep = "_"),
              patterns = patterns, state = state_draws, beta = beta_draws,
              alpha = alpha_draws, Vb = Vb_draws, Pi = Pi_draws,
              fitted = fitted_draws, transitions = transitions,
              seed = seed, update_vb = update_vb, update_pi = update_pi)
  mt_validate_draws(ans$Vb)
  if (any(abs(rowSums(ans$Pi) - 1) > 1e-12)) stop("stored Pi is not normalized.")
  ans
}

mt_pattern_summary <- function(fit) {
  patterns <- fit$patterns
  K <- nrow(patterns)
  M <- ncol(fit$state)
  pattern_probability <- matrix(0, M, K)
  for (j in seq_len(M)) for (s in seq_len(K)) {
    pattern_probability[j, s] <- mean(fit$state[, j] == s)
  }
  list(pattern_probability = pattern_probability,
       trait_pip = pattern_probability %*% patterns,
       pleiotropic_probability = pattern_probability[, rowSums(patterns) == ncol(patterns)],
       alpha_mean = apply(fit$alpha, c(2L, 3L), mean),
       Vb_mean = apply(fit$Vb, c(1L, 2L), mean),
       Pi_mean = colMeans(fit$Pi),
       fitted_mean = apply(fit$fitted, c(2L, 3L), mean))
}
