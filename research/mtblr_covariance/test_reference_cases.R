script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
research_dir <- if (length(script_arg)) {
  dirname(normalizePath(sub("^--file=", "", script_arg[1L]), mustWork = TRUE))
} else {
  normalizePath("research/mtblr_covariance", mustWork = TRUE)
}
for (file in c("mtblr_reference_model.R", "mtblr_exact_reference.R",
               "mtblr_full_latent.R", "mtblr_active_only.R",
               "mtblr_current_hybrid.R", "compare_samplers.R")) {
  source(file.path(research_dir, file), local = FALSE)
}

expect_close <- function(actual, expected, tolerance, label) {
  error <- max(abs(actual - expected))
  if (!is.finite(error) || error > tolerance) {
    stop(sprintf("%s: max error %.6g exceeds %.6g", label, error, tolerance),
         call. = FALSE)
  }
  invisible(error)
}

f <- fixture(2L)
nu0 <- 7
prior_mean <- matrix(c(0.30, 0.06, 0.06, 0.24), 2L)
Psi0 <- mt_iw_scale_from_mean(prior_mean, nu0)
stopifnot(max(abs(mt_iw_mean(nu0, Psi0) - prior_mean)) < 1e-14)
stopifnot(inherits(try(mt_rinvwishart(1, diag(2)), silent = TRUE), "try-error"))

# Known-Vb BayesC: exact enumeration is independent of either transition.
pi <- 0.18
exact <- mt_exact_bayesc_known_vb(f$X, f$Y, f$Ve, f$Vb, pi)
full_known <- mtblr_full_latent(f$X, f$Y, f$Ve, pi, nu0, Psi0,
                                26000L, 4000L, 2021L, f$Vb, FALSE)
active_known <- mtblr_active_only(f$X, f$Y, f$Ve, pi, nu0, Psi0,
                                  26000L, 4000L, 2022L, f$Vb, FALSE)
expect_close(colMeans(full_known$delta), exact$pip, 0.025,
             "full-latent known-Vb PIP")
expect_close(colMeans(active_known$delta), exact$pip, 0.025,
             "active-only known-Vb PIP")
expect_close(apply(full_known$alpha, c(2L, 3L), mean), exact$alpha_mean, 0.035,
             "full-latent known-Vb alpha mean")
expect_close(apply(active_known$alpha, c(2L, 3L), mean), exact$alpha_mean, 0.035,
             "active-only known-Vb alpha mean")

# Reproducibility and positive-definite covariance draws.
full_a <- mtblr_full_latent(f$X, f$Y, f$Ve, pi, nu0, Psi0, 2000L, 500L, 71L)
full_b <- mtblr_full_latent(f$X, f$Y, f$Ve, pi, nu0, Psi0, 2000L, 500L, 71L)
stopifnot(identical(full_a$delta, full_b$delta), identical(full_a$Vb, full_b$Vb))
mt_validate_draws(full_a$Vb)
active_a <- mtblr_active_only(f$X, f$Y, f$Ve, pi, nu0, Psi0,
                              1200L, 200L, 72L)
active_b <- mtblr_active_only(f$X, f$Y, f$Ve, pi, nu0, Psi0,
                              1200L, 200L, 72L)
stopifnot(identical(active_a$delta, active_b$delta),
          identical(active_a$Vb, active_b$Vb))
hybrid_a <- mtblr_current_hybrid(f$X, f$Y, f$Ve, pi, nu0, Psi0,
                                 600L, 100L, 73L)
hybrid_b <- mtblr_current_hybrid(f$X, f$Y, f$Ve, pi, nu0, Psi0,
                                 600L, 100L, 73L)
stopifnot(identical(hybrid_a$delta, hybrid_b$delta),
          identical(hybrid_a$Vb, hybrid_b$Vb),
          identical(hybrid_a$Vb_latent, hybrid_b$Vb_latent))

# BayesR component scaling and exact state probabilities at known Vb.  The
# prototype stores scaled component effects theta_j with covariance
# gamma[k_j] * q[j] * Vb; alpha_j equals theta_j for an active shared component.
gamma <- c(0, 0.1, 1)
component_probability <- c(0.55, 0.25, 0.20)
exact_r <- mt_exact_bayesr_known_vb(f$X[, 1L, drop = FALSE], f$Y,
                                    f$Ve, f$Vb,
                                    component_probability, gamma)
exact_r_q1 <- mt_exact_bayesr_known_vb(f$X[, 1L, drop = FALSE], f$Y,
                                       f$Ve, f$Vb,
                                       component_probability, gamma, q = 1)
stopifnot(identical(exact_r$state_probability, exact_r_q1$state_probability),
          identical(exact_r$alpha_mean, exact_r_q1$alpha_mean))
bayesr_default_q <- mtblr_bayesr_active_only(
  f$X[, 1L, drop = FALSE], f$Y, f$Ve,
  component_probability, gamma, nu0, Psi0,
  300L, 50L, 3030L, f$Vb, FALSE
)
bayesr_explicit_q1 <- mtblr_bayesr_active_only(
  f$X[, 1L, drop = FALSE], f$Y, f$Ve,
  component_probability, gamma, nu0, Psi0,
  300L, 50L, 3030L, f$Vb, FALSE, q = 1
)
stopifnot(identical(bayesr_default_q$state, bayesr_explicit_q1$state),
          identical(bayesr_default_q$theta, bayesr_explicit_q1$theta),
          identical(bayesr_default_q$Vb, bayesr_explicit_q1$Vb))
stopifnot(abs(sum(exact_r$state_probability) - 1) < 1e-14)
bayesr_known <- mtblr_bayesr_active_only(
  f$X[, 1L, drop = FALSE], f$Y, f$Ve,
  component_probability, gamma, nu0, Psi0,
  26000L, 4000L, 3031L, f$Vb, FALSE
)
state_frequency <- tabulate(bayesr_known$state[, 1L] + 1L,
                            nbins = length(gamma)) / nrow(bayesr_known$state)
expect_close(state_frequency, exact_r$state_probability, 0.025,
             "BayesR known-Vb state probabilities")

# Unequal marker scales enter both the exact likelihood and sampler once.
q_unequal <- c(0.35, 1.8)
exact_r_q <- mt_exact_bayesr_known_vb(
  f$X, f$Y, f$Ve, f$Vb, component_probability, gamma, q_unequal
)
bayesr_q <- mtblr_bayesr_active_only(
  f$X, f$Y, f$Ve, component_probability, gamma, nu0, Psi0,
  10000L, 2000L, 3032L, f$Vb, FALSE, q_unequal
)
component_frequency <- matrix(0, ncol(f$X), length(gamma))
for (j in seq_len(ncol(f$X))) for (k in 0:(length(gamma) - 1L)) {
  component_frequency[j, k + 1L] <- mean(bayesr_q$state[, j] == k)
}
expect_close(component_frequency, exact_r_q$component_probability, 0.05,
             "BayesR unequal-q component probabilities")
expect_close(apply(bayesr_q$alpha, c(2L, 3L), mean), exact_r_q$alpha_mean,
             0.07, "BayesR unequal-q realised-effect means")
stopifnot(identical(bayesr_q$q, q_unequal),
          identical(bayesr_q$effect_storage, "scaled_theta"))

# The covariance statistic removes gamma[k_j] * q_j exactly once from theta.
theta_demo <- rbind(c(0, 0), c(0.2, -0.1), c(0.3, 0.4))
q_demo <- c(0.7, 0.4, 2.5)
stat <- mt_covariance_sufficient(theta_demo, c(0L, 1L, 2L), gamma,
                                 q_demo, storage = "scaled")$statistic
manual <- tcrossprod(theta_demo[2L, ]) / (gamma[2L] * q_demo[2L]) +
  tcrossprod(theta_demo[3L, ]) / (gamma[3L] * q_demo[3L])
gamma_only <- tcrossprod(theta_demo[2L, ]) / gamma[2L] +
  tcrossprod(theta_demo[3L, ]) / gamma[3L]
double_removed <- tcrossprod(theta_demo[2L, ]) / (gamma[2L] * q_demo[2L])^2 +
  tcrossprod(theta_demo[3L, ]) / (gamma[3L] * q_demo[3L])^2
stopifnot(max(abs(stat - manual)) < 1e-14,
          max(abs(stat - gamma_only)) > 1e-3,
          max(abs(stat - double_removed)) > 1e-3)
stopifnot(inherits(try(mt_validate_marker_scale(c(1, 0), 2L), silent = TRUE),
                   "try-error"),
          inherits(try(mt_validate_marker_scale(1, 2L), silent = TRUE),
                   "try-error"),
          inherits(try(mt_validate_marker_scale(c("1", "2"), 2L), silent = TRUE),
                   "try-error"),
          inherits(try(mt_validate_marker_scale(c(1, Inf), 2L), silent = TRUE),
                   "try-error"))

# A low-dimensional numerical reference, evaluated at two grid resolutions.
g <- fixture(1L)
grid11 <- mt_vb_grid_reference(g$X, g$Y, g$Ve, pi, nu0, Psi0, grid_n = 11L)
grid15 <- mt_vb_grid_reference(g$X, g$Y, g$Ve, pi, nu0, Psi0, grid_n = 15L)
expect_close(grid11$mean, grid15$mean, 0.05, "Vb grid refinement")
if (grid15$boundary_mass_upper_bound > 0.12) {
  stop("Vb grid has excessive boundary mass; widen the numerical domain.", call. = FALSE)
}

# Unknown-Vb samplers target the same marginal model; use independent chains.
full_unknown <- mtblr_full_latent(g$X, g$Y, g$Ve, pi, nu0, Psi0,
                                  16000L, 3000L, 4101L)
active_unknown <- mtblr_active_only(g$X, g$Y, g$Ve, pi, nu0, Psi0,
                                    16000L, 3000L, 4102L)
full_vb_mean <- apply(full_unknown$Vb, c(1L, 2L), mean)
active_vb_mean <- apply(active_unknown$Vb, c(1L, 2L), mean)
expect_close(full_vb_mean, grid15$mean, 0.09, "full-latent unknown-Vb mean")
expect_close(active_vb_mean, grid15$mean, 0.09, "active-only unknown-Vb mean")
expect_close(full_vb_mean, active_vb_mean, 0.08, "sampler Vb agreement")
expect_close(colMeans(full_unknown$delta), colMeans(active_unknown$delta), 0.04,
             "sampler PIP agreement")

# Inactive markers affect augmentation mixing, not the active-only IW statistic.
inactive_stat <- mt_covariance_sufficient(matrix(0, 8L, 2L), integer(8L))
stopifnot(inactive_stat$active == 0L, all(inactive_stat$statistic == 0))

cat("All MTBLR covariance reference cases passed.\n")
