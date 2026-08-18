source(file.path("tests", "research", "mtblr_covariance",
                 "mtblr_pattern_scale_reference.R"))

patterns <- rbind(`0_0` = c(0L, 0L), `1_0` = c(1L, 0L),
                  `0_1` = c(0L, 1L), `1_1` = c(1L, 1L))
Vb <- matrix(c(0.7, 0.18, 0.18, 0.45), 2L)
h <- c(0.8, -0.35)
L <- matrix(c(1.4, 0.22, 0.22, 1.1), 2L)
pi <- c(0.55, 0.15, 0.12, 0.18)
gamma <- c(0.1, 1)
omega <- c(0.7, 0.3)
q <- 1.25

fit <- mtblr_pattern_scale_conditional(h, L, Vb, patterns, pi,
                                        gamma, omega, q)
stopifnot(nrow(fit$states) == 7L, sum(is.na(fit$states$scale)) == 1L,
          abs(sum(fit$probability) - 1) < 1e-14,
          all(is.finite(fit$probability)), all(fit$probability > 0))

# Independently check one partial-pattern state by direct Gaussian algebra.
state <- which(fit$states$pattern == 2L & fit$states$scale == 2L)
prior <- gamma[2L] * q * Vb[1L, 1L, drop = FALSE]
Q <- solve(prior) + L[1L, 1L, drop = FALSE]
stopifnot(max(abs(fit$active_covariance[[state]] - solve(Q))) < 1e-14,
          max(abs(fit$active_mean[[state]] - solve(Q, h[1L]))) < 1e-14)

# One positive scale at one reduces exactly to the Cheng active-pattern target.
cheng <- mtblr_pattern_scale_conditional(h, L, Vb, patterns, pi, 1, 1, 1)
manual <- numeric(nrow(patterns))
manual[1L] <- log(pi[1L])
for (pattern in 2:nrow(patterns)) {
  A <- which(patterns[pattern, ] == 1L)
  prior <- Vb[A, A, drop = FALSE]
  Q <- solve(prior) + L[A, A, drop = FALSE]
  m <- solve(Q, h[A])
  manual[pattern] <- log(pi[pattern]) - determinant(prior, TRUE)$modulus / 2 -
    determinant(Q, TRUE)$modulus / 2 + drop(crossprod(h[A], m)) / 2
}
manual <- exp(manual - max(manual)); manual <- manual / sum(manual)
stopifnot(max(abs(cheng$probability - manual)) < 1e-14)

# T=1 is the usual null plus positive normal-scale mixture.
univariate <- mtblr_pattern_scale_conditional(
  h = 0.6, likelihood_precision = matrix(1.7), Vb = matrix(0.4),
  patterns = rbind(0L, 1L), pattern_probability = c(0.8, 0.2),
  scales = c(0.01, 0.1, 1), scale_probability = c(0.6, 0.3, 0.1), q = 1.3)
stopifnot(length(univariate$probability) == 4L,
          abs(sum(univariate$probability) - 1) < 1e-14)

# Null plus all-active patterns is the multivariate BayesR scale mixture.
all_active <- mtblr_pattern_scale_conditional(
  h, L, Vb, patterns = rbind(c(0L, 0L), c(1L, 1L)),
  pattern_probability = c(0.55, 0.45), scales = gamma,
  scale_probability = omega, q = q)
manual_all <- c(log(0.55), rep(NA_real_, length(gamma)))
for (scale_index in seq_along(gamma)) {
  prior <- gamma[scale_index] * q * Vb
  Q <- solve(prior) + L
  mean <- solve(Q, h)
  manual_all[scale_index + 1L] <- log(0.45) + log(omega[scale_index]) -
    determinant(prior, TRUE)$modulus / 2 -
    determinant(Q, TRUE)$modulus / 2 + drop(crossprod(h, mean)) / 2
}
manual_all <- exp(manual_all - max(manual_all))
manual_all <- manual_all / sum(manual_all)
stopifnot(max(abs(all_active$probability - manual_all)) < 1e-14)

# Scaled theta outer products remove gamma*q exactly once.
base <- rbind(c(0.4, -0.2), c(-0.1, 0.3))
scale <- c(0.1, 1.0); marker_q <- c(1.25, 0.8)
theta <- base * sqrt(scale * marker_q)
statistic <- mtblr_pattern_scale_vb_statistic(
  theta, pattern = c(2L, 4L), scale = scale, q = marker_q)
stopifnot(max(abs(statistic - crossprod(base))) < 1e-14)

cat("Phase 7A pattern-by-scale reference: PASS\n")
