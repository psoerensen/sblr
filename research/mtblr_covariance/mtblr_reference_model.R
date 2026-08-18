# Transparent utilities for the standalone MTBLR covariance prototype.

mt_symmetrize <- function(x) (x + t(x)) / 2

mt_assert_spd <- function(x, name = deparse(substitute(x))) {
  x <- as.matrix(x)
  if (!is.numeric(x) || nrow(x) < 1L || nrow(x) != ncol(x) ||
      any(!is.finite(x))) {
    stop(name, " must be a finite nonempty square matrix.", call. = FALSE)
  }
  if (max(abs(x - t(x))) > 1e-10) {
    stop(name, " must be symmetric.", call. = FALSE)
  }
  ch <- tryCatch(chol(mt_symmetrize(x)), error = function(e) NULL)
  if (is.null(ch)) stop(name, " must be positive definite.", call. = FALSE)
  invisible(x)
}

mt_logdet_spd <- function(x) {
  mt_assert_spd(x)
  2 * sum(log(diag(chol(x))))
}

mt_validate_iw <- function(df, scale, require_mean = FALSE) {
  mt_assert_spd(scale, "inverse-Wishart scale")
  p <- nrow(scale)
  if (length(df) != 1L || !is.finite(df) || df <= p - 1) {
    stop("inverse-Wishart df must be finite and greater than T - 1.",
         call. = FALSE)
  }
  if (require_mean && df <= p + 1) {
    stop("inverse-Wishart mean requires df greater than T + 1.",
         call. = FALSE)
  }
  invisible(TRUE)
}

mt_iw_scale_from_mean <- function(mean, df) {
  mt_assert_spd(mean, "inverse-Wishart prior mean")
  mt_validate_iw(df, mean, require_mean = TRUE)
  (df - nrow(mean) - 1) * mean
}

mt_iw_mean <- function(df, scale) {
  mt_validate_iw(df, scale, require_mean = TRUE)
  scale / (df - nrow(scale) - 1)
}

mt_rinvwishart <- function(df, scale) {
  mt_validate_iw(df, scale)
  out <- solve(stats::rWishart(1L, df, solve(scale))[, , 1L])
  out <- mt_symmetrize(out)
  mt_assert_spd(out, "inverse-Wishart draw")
  out
}

mt_multigamma_log <- function(a, p) {
  p * (p - 1) * log(pi) / 4 + sum(lgamma(a + (1 - seq_len(p)) / 2))
}

mt_log_iwishart <- function(V, df, scale) {
  mt_validate_iw(df, scale)
  mt_assert_spd(V, "V")
  p <- nrow(V)
  0.5 * df * mt_logdet_spd(scale) -
    0.5 * df * p * log(2) - mt_multigamma_log(df / 2, p) -
    0.5 * (df + p + 1) * mt_logdet_spd(V) -
    0.5 * sum(diag(scale %*% solve(V)))
}

mt_rmvnorm <- function(mean, covariance) {
  mean <- as.numeric(mean)
  mt_assert_spd(covariance, "normal covariance")
  if (length(mean) != nrow(covariance) || any(!is.finite(mean))) {
    stop("normal mean and covariance dimensions differ.", call. = FALSE)
  }
  as.numeric(mean + t(chol(covariance)) %*% stats::rnorm(length(mean)))
}

mt_log_mvn_zero <- function(x, covariance) {
  x <- as.numeric(x)
  mt_assert_spd(covariance, "marginal covariance")
  if (length(x) != nrow(covariance) || any(!is.finite(x))) {
    stop("normal observation and covariance dimensions differ.", call. = FALSE)
  }
  -0.5 * (length(x) * log(2 * pi) + mt_logdet_spd(covariance) +
            drop(crossprod(x, solve(covariance, x))))
}

mt_validate_data <- function(X, Y, Ve) {
  X <- as.matrix(X); Y <- as.matrix(Y); Ve <- as.matrix(Ve)
  if (!is.numeric(X) || !is.numeric(Y) || nrow(X) != nrow(Y) ||
      nrow(X) < 1L || ncol(X) < 1L || ncol(Y) != nrow(Ve) ||
      any(!is.finite(X)) || any(!is.finite(Y))) {
    stop("X, Y, and Ve have incompatible or nonfinite dimensions.",
         call. = FALSE)
  }
  mt_assert_spd(Ve, "Ve")
  list(X = X, Y = Y, Ve = Ve, N = nrow(X), M = ncol(X), T = ncol(Y))
}

mt_validate_mcmc_controls <- function(n_iter, burn, seed) {
  controls <- list(n_iter = n_iter, burn = burn, seed = seed)
  for (name in names(controls)) {
    value <- controls[[name]]
    if (length(value) != 1L || !is.numeric(value) || !is.finite(value) ||
        value != floor(value)) {
      stop(name, " must be one finite integer-valued number.", call. = FALSE)
    }
  }
  if (burn < 0 || n_iter <= burn) {
    stop("n_iter must exceed burn >= 0.", call. = FALSE)
  }
  if (seed < 0 || seed > .Machine$integer.max) {
    stop("seed must be between 0 and .Machine$integer.max.", call. = FALSE)
  }
  list(n_iter = as.integer(n_iter), burn = as.integer(burn), seed = as.integer(seed))
}

mt_validate_marker_scale <- function(q, M) {
  if (is.null(q)) q <- rep(1, M)
  if (!is.numeric(q) || !is.null(dim(q)) || length(q) != M ||
      any(!is.finite(q)) || any(q <= 0)) {
    stop("q must be a finite, strictly positive numeric vector of length M.",
         call. = FALSE)
  }
  as.numeric(q)
}

mt_marker_conditional <- function(partial, x, Vb, Ve, prior_probability) {
  if (length(prior_probability) != 1L || !is.finite(prior_probability) ||
      prior_probability <= 0 || prior_probability >= 1) {
    stop("prior_probability must lie strictly inside (0, 1).", call. = FALSE)
  }
  mt_assert_spd(Vb, "Vb"); mt_assert_spd(Ve, "Ve")
  Ei <- solve(Ve); Bi <- solve(Vb)
  precision <- Bi + sum(x^2) * Ei
  score <- as.numeric(Ei %*% as.numeric(crossprod(x, partial)))
  covariance <- solve(precision)
  mean <- as.numeric(covariance %*% score)
  log_bf <- 0.5 * (mt_logdet_spd(Bi) - mt_logdet_spd(precision) +
                     sum(score * mean))
  probability <- stats::plogis(stats::qlogis(prior_probability) + log_bf)
  list(probability = probability, mean = mean,
       covariance = mt_symmetrize(covariance), log_bf = log_bf)
}

mt_component_conditional <- function(partial, x, Vb, Ve,
                                     component_probability, gamma, q = 1) {
  if (length(component_probability) != length(gamma) || gamma[1L] != 0 ||
      any(!is.finite(component_probability)) ||
      any(component_probability < 0) || sum(component_probability) <= 0 ||
      any(!is.finite(gamma)) || any(gamma[-1L] <= 0)) {
    stop("invalid BayesR component probabilities or multipliers.",
         call. = FALSE)
  }
  if (length(q) != 1L || !is.numeric(q) || !is.finite(q) || q <= 0) {
    stop("q must be one finite, strictly positive marker multiplier.",
         call. = FALSE)
  }
  component_probability <- component_probability / sum(component_probability)
  Ei <- solve(Ve); Bi <- solve(Vb)
  score <- as.numeric(Ei %*% as.numeric(crossprod(x, partial)))
  log_weight <- rep(-Inf, length(gamma)); candidates <- vector("list", length(gamma))
  if (component_probability[1L] > 0) log_weight[1L] <- log(component_probability[1L])
  for (k in 2:length(gamma)) {
    if (component_probability[k] == 0) next
    prior_precision <- Bi / (gamma[k] * q)
    precision <- prior_precision + sum(x^2) * Ei
    covariance <- solve(precision)
    mean <- as.numeric(covariance %*% score)
    log_weight[k] <- log(component_probability[k]) + 0.5 * (
      mt_logdet_spd(prior_precision) - mt_logdet_spd(precision) +
        sum(score * mean))
    candidates[[k]] <- list(mean = mean, covariance = mt_symmetrize(covariance))
  }
  weight <- exp(log_weight - max(log_weight)); weight <- weight / sum(weight)
  list(probability = weight, candidates = candidates)
}

mt_covariance_sufficient <- function(effect, state, gamma = NULL, q = NULL,
                                     storage = c("base", "scaled")) {
  effect <- as.matrix(effect); state <- as.integer(state)
  storage <- match.arg(storage)
  if (nrow(effect) != length(state) || any(!is.finite(effect))) {
    stop("effect rows must align with finite marker states.", call. = FALSE)
  }
  active <- state > 0L
  if (storage == "base") {
    weights <- rep(1, sum(active))
  } else {
    if (is.null(gamma)) stop("scaled storage requires gamma.", call. = FALSE)
    if (any(state < 0L) || any(state + 1L > length(gamma)) || gamma[1L] != 0 ||
        any(!is.finite(gamma)) || any(gamma[-1L] <= 0)) {
      stop("invalid BayesR state or gamma.", call. = FALSE)
    }
    q <- mt_validate_marker_scale(q, nrow(effect))
    weights <- gamma[state[active] + 1L] * q[active]
  }
  statistic <- matrix(0, ncol(effect), ncol(effect))
  if (any(active)) for (i in seq_along(which(active))) {
    value <- effect[which(active)[i], ]
    statistic <- statistic + tcrossprod(value) / weights[i]
  }
  list(statistic = mt_symmetrize(statistic), active = sum(active),
       storage = storage)
}

mt_chain_diagnostics <- function(x) {
  x <- as.numeric(x)
  if (length(x) < 3L || any(!is.finite(x))) stop("diagnostic input is invalid.", call. = FALSE)
  ac <- stats::acf(x, plot = FALSE, lag.max = min(1000L, length(x) - 1L))$acf[-1L]
  first_nonpositive <- which(ac <= 0)[1L]
  positive <- if (is.na(first_nonpositive)) ac else ac[seq_len(first_nonpositive - 1L)]
  tau <- 1 + 2 * sum(positive)
  list(lag1 = unname(stats::acf(x, plot = FALSE, lag.max = 1L)$acf[2L]),
       ess = length(x) / max(tau, 1), mean = mean(x), variance = stats::var(x))
}

mt_validate_draws <- function(draws) {
  if (length(dim(draws)) != 3L || dim(draws)[1L] != dim(draws)[2L]) {
    stop("covariance draws must be T by T by draws.", call. = FALSE)
  }
  for (i in seq_len(dim(draws)[3L])) mt_assert_spd(draws[, , i], "Vb draw")
  invisible(TRUE)
}
