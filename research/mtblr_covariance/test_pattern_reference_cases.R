script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
research_dir <- if (length(script_arg)) {
  dirname(normalizePath(sub("^--file=", "", script_arg[1L]), mustWork = TRUE))
} else normalizePath("research/mtblr_covariance", mustWork = TRUE)
for (file in c("mtblr_reference_model.R", "mtblr_exact_reference.R",
               "mtblr_pattern_reference.R", "mtblr_pattern_samplers.R",
               "mtblr_regional.R", "compare_samplers.R")) {
  source(file.path(research_dir, file), local = FALSE)
}

expect_close <- function(actual, expected, tolerance, label) {
  error <- max(abs(actual - expected))
  if (!is.finite(error) || error > tolerance) {
    stop(sprintf("%s: max error %.6g exceeds %.6g", label, error, tolerance),
         call. = FALSE)
  }
}

expect_error_pattern <- function(expression, pattern, label) {
  message <- tryCatch({
    force(expression)
    NA_character_
  }, error = conditionMessage)
  if (is.na(message) || !grepl(pattern, message, perl = TRUE)) {
    stop(label, ": expected error matching `", pattern, "`, got `",
         message, "`.", call. = FALSE)
  }
  invisible(message)
}

patterns <- mt_pattern_space()
stopifnot(identical(unname(patterns), matrix(c(0L, 1L, 0L, 1L,
                                                  0L, 0L, 1L, 1L), 4L, 2L)))
Vb <- matrix(c(.40, .16, .16, .30), 2L)

# Conditional completion is checked against its independent Gaussian formula.
set.seed(1401)
active <- .55
draw <- replicate(12000L, mt_complete_latent(c(1L, 0L), active, Vb)[2L])
expected_mean <- Vb[2L, 1L] / Vb[1L, 1L] * active
expected_variance <- Vb[2L, 2L] - Vb[2L, 1L]^2 / Vb[1L, 1L]
expect_close(mean(draw), expected_mean, .012, "conditional completion mean")
expect_close(var(draw), expected_variance, .015, "conditional completion variance")

f <- fixture(2L)
nu0 <- 7
Psi0 <- mt_iw_scale_from_mean(matrix(c(.30, .08, .08, .25), 2L), nu0)
prior <- c(6, 1, 1, 3)
exact <- mt_pattern_configuration_reference(f$X, f$Y, f$Ve, f$Vb,
                                            patterns, dirichlet_prior = prior)
stopifnot(abs(sum(exact$configuration_probability) - 1) < 1e-13,
          all(abs(rowSums(exact$pattern_probability) - 1) < 1e-13))

# Known-Vb joint and coordinate kernels both agree with exact enumeration.
fits <- list(
  joint = mtblr_pattern_sampler(f$X, f$Y, f$Ve, patterns, prior, nu0, Psi0,
                                8000L, 1200L, 1501L, "completed_active", "joint",
                                f$Vb, prior / sum(prior), FALSE, TRUE),
  coordinate = mtblr_pattern_sampler(f$X, f$Y, f$Ve, patterns, prior, nu0, Psi0,
                                     10000L, 1500L, 1502L, "completed_active", "coordinate",
                                     f$Vb, prior / sum(prior), FALSE, TRUE)
)
for (name in names(fits)) {
  summary <- mt_pattern_summary(fits[[name]])
  expect_close(summary$pattern_probability, exact$pattern_probability, .055,
               paste(name, "pattern probability"))
  expect_close(summary$trait_pip, exact$trait_pip, .055,
               paste(name, "trait PIP"))
  expect_close(summary$Pi_mean, exact$pi_mean, .04,
               paste(name, "learned Pi"))
  expect_close(summary$alpha_mean, exact$alpha_mean, .07,
               paste(name, "realised effect"))
}

# Full latent and completed-active samplers agree with unknown Vb.
full <- mtblr_pattern_sampler(f$X, f$Y, f$Ve, patterns, prior, nu0, Psi0,
                              7000L, 1000L, 1601L, "full", "joint")
completed <- mtblr_pattern_sampler(f$X, f$Y, f$Ve, patterns, prior, nu0, Psi0,
                                   7000L, 1000L, 1602L, "completed_active", "joint")
sf <- mt_pattern_summary(full)
sc <- mt_pattern_summary(completed)
expect_close(sf$pattern_probability, sc$pattern_probability, .06,
             "full versus completed pattern probabilities")
expect_close(sf$Vb_mean, sc$Vb_mean, .11, "full versus completed Vb")
expect_close(sf$fitted_mean, sc$fitted_mean, .07, "full versus completed predictions")
mt_validate_draws(full$Vb)
mt_validate_draws(completed$Vb)

# Same-seed reproducibility.
a <- mtblr_pattern_sampler(f$X, f$Y, f$Ve, patterns, prior, nu0, Psi0,
                           900L, 200L, 1701L, "completed_active", "coordinate")
b <- mtblr_pattern_sampler(f$X, f$Y, f$Ve, patterns, prior, nu0, Psi0,
                           900L, 200L, 1701L, "completed_active", "coordinate")
stopifnot(identical(a$state, b$state), identical(a$Vb, b$Vb), identical(a$Pi, b$Pi))

# Coordinate updating is valid only when one-coordinate moves connect the
# declared pattern graph. Joint updating remains valid for restricted spaces.
stopifnot(mt_pattern_graph_connected(patterns))
restricted_patterns <- rbind(null = c(0L, 0L), shared = c(1L, 1L))
stopifnot(!mt_pattern_graph_connected(restricted_patterns))
expect_error_pattern(
  mtblr_pattern_sampler(f$X, f$Y, f$Ve, restricted_patterns, c(2, 1),
                        nu0, Psi0, 20L, 5L, 1702L,
                        "completed_active", "coordinate"),
  "connected by one-coordinate moves", "disconnected coordinate graph"
)
restricted_joint <- mtblr_pattern_sampler(
  f$X, f$Y, f$Ve, restricted_patterns, c(2, 1), nu0, Psi0,
  120L, 20L, 1703L, "completed_active", "joint"
)
stopifnot(all(restricted_joint$state %in% 1:2))

# Two fixed regional covariances agree with exact enumeration.
region <- c("A", "B")
region_levels <- c("A", "B")
regional_vb <- list(matrix(c(.42, .12, .12, .26), 2L),
                    matrix(c(.20, -.04, -.04, .36), 2L))
regional_pi <- rbind(c(.55, .20, .10, .15), c(.35, .10, .20, .35))
regional_exact <- mt_regional_fixed_reference(
  f$X, f$Y, f$Ve, patterns, region, region_levels, regional_vb, regional_pi
)
regional_fit <- mtblr_regional_sampler(
  f$X, f$Y, f$Ve, region, region_levels, patterns, rep(1, 4), nu0, Psi0,
  8000L, 1200L, 1801L, "completed_active", "fixed", regional_vb,
  TRUE, FALSE, regional_pi
)
regional_marginal <- matrix(0, 2L, 4L)
for (j in 1:2) for (s in 1:4) regional_marginal[j, s] <-
  mean(regional_fit$state[, j] == s)
expect_close(regional_marginal, regional_exact$pattern_probability, .05,
             "fixed regional pattern probabilities")

# Regional full and completed-active augmentations agree marginally.
regional_full <- mtblr_regional_sampler(
  f$X, f$Y, f$Ve, region, region_levels, patterns, rep(1, 4), nu0, Psi0,
  4000L, 700L, 1811L, "full", "regional"
)
regional_completed <- mtblr_regional_sampler(
  f$X, f$Y, f$Ve, region, region_levels, patterns, rep(1, 4), nu0, Psi0,
  4000L, 700L, 1812L, "completed_active", "regional"
)
regional_probability <- function(fit) {
  out <- matrix(0, 2L, 4L)
  for (j in 1:2) for (s in 1:4) out[j, s] <- mean(fit$state[, j] == s)
  out
}
expect_close(regional_probability(regional_full),
             regional_probability(regional_completed), .08,
             "regional augmentation agreement")
expect_close(apply(regional_full$Vb, c(1L, 2L, 3L), mean),
             apply(regional_completed$Vb, c(1L, 2L, 3L), mean), .14,
             "regional covariance agreement")

# Global mode has one persistent covariance shared by all declared regions.
global_fit <- mtblr_regional_sampler(
  f$X, f$Y, f$Ve, region, region_levels, patterns, rep(1, 4), nu0, Psi0,
  800L, 150L, 1821L, "completed_active", "global"
)
stopifnot(identical(global_fit$Vb[, , 1L, ], global_fit$Vb[, , 2L, ]))

# Shared-global Pi uses one declared prior, not the sum of replicated regional
# priors. With a global covariance, identical schedule, and identical seed, the
# regional wrapper reduces exactly to the ordinary Cheng sampler.
global_prior <- c(4, 1, 2, 3)
ordinary_global <- mtblr_pattern_sampler(
  f$X, f$Y, f$Ve, patterns, global_prior, nu0, Psi0,
  900L, 150L, 1831L, "completed_active", "joint"
)
regional_global <- mtblr_regional_sampler(
  f$X, f$Y, f$Ve, region, region_levels, patterns,
  dirichlet_prior = matrix(rep(1, 8), 2, 4),
  nu0 = nu0, Psi0 = Psi0, n_iter = 900L, burn = 150L, seed = 1831L,
  augmentation = "completed_active", covariance_mode = "global",
  regional_pi = FALSE, global_dirichlet_prior = global_prior
)
stopifnot(identical(ordinary_global$state, regional_global$state),
          identical(ordinary_global$Pi, regional_global$Pi[, 1L, ]),
          identical(ordinary_global$Vb, regional_global$Vb[, , 1L, ]),
          identical(regional_global$Pi[, 1L, ], regional_global$Pi[, 2L, ]),
          identical(regional_global$Vb[, , 1L, ], regional_global$Vb[, , 2L, ]),
          identical(regional_global$global_dirichlet_prior, global_prior))
ordinary_summary <- mt_pattern_summary(ordinary_global)
regional_pattern_probability <- matrix(0, ncol(f$X), nrow(patterns))
for (j in seq_len(ncol(f$X))) for (s in seq_len(nrow(patterns))) {
  regional_pattern_probability[j, s] <- mean(regional_global$state[, j] == s)
}
stopifnot(identical(ordinary_summary$pattern_probability,
                    regional_pattern_probability),
          identical(ordinary_summary$trait_pip,
                    regional_pattern_probability %*% patterns))

# Regional controls and dimensions fail before allocation or sampling. A
# single retained covariance draw preserves its T x T x 1 array shape.
expect_error_pattern(
  mtblr_regional_sampler(f$X, f$Y, f$Ve, region, region_levels, patterns,
                         rep(1, 4), nu0, Psi0, 10.5, 2L, 1841L),
  "n_iter must be one finite integer", "fractional n_iter"
)
expect_error_pattern(
  mtblr_regional_sampler(f$X, f$Y, f$Ve, region, region_levels, patterns,
                         rep(1, 4), nu0, Psi0, 10L, 10L, 1841L),
  "n_iter must exceed burn", "invalid burn"
)
expect_error_pattern(
  mtblr_regional_sampler(f$X, f$Y, f$Ve, region, region_levels, patterns,
                         rep(1, 4), nu0, Psi0, 10L, 2L, NA_real_),
  "seed must be one finite integer", "nonfinite seed"
)
expect_error_pattern(
  mtblr_regional_sampler(f$X, f$Y, f$Ve, region, c("A", "A"), patterns,
                         rep(1, 4), nu0, Psi0, 10L, 2L, 1841L),
  "region_levels must be unique", "duplicated region levels"
)
expect_error_pattern(
  mtblr_regional_sampler(f$X, f$Y, f$Ve, region, region_levels, patterns,
                         rep(1, 4), nu0, Psi0, 10L, 2L, 1841L,
                         covariance_mode = "fixed", Vb_fixed = diag(3)),
  "T by T", "fixed covariance dimensions"
)
expect_error_pattern(
  mtblr_regional_sampler(f$X, f$Y, f$Ve, region, region_levels, patterns,
                         rep(1, 4), nu0, Psi0, 10L, 2L, 1841L,
                         covariance_mode = "global", regional_pi = FALSE),
  "global_dirichlet_prior", "missing global Pi prior"
)
expect_error_pattern(
  mtblr_regional_sampler(f$X, f$Y, f$Ve, region, region_levels, patterns,
                         rep(1, 4), nu0, Psi0, 10L, 2L, 1841L,
                         covariance_mode = "global", regional_pi = FALSE,
                         global_dirichlet_prior = c(1, 1)),
  "global_dirichlet_prior", "wrong global Pi prior length"
)
expect_error_pattern(
  mtblr_regional_sampler(f$X, f$Y, f$Ve, region, region_levels, patterns,
                         matrix(1, 1, 4), nu0, Psi0, 10L, 2L, 1841L),
  "regional Dirichlet prior", "wrong regional Pi prior dimensions"
)
single_draw <- mtblr_regional_sampler(
  f$X, f$Y, f$Ve, region, region_levels, patterns, rep(1, 4), nu0, Psi0,
  1L, 0L, 1842L, "completed_active", "fixed", regional_vb,
  TRUE, FALSE, regional_pi
)
stopifnot(identical(dim(single_draw$Vb), c(2L, 2L, 2L, 1L)))

# Empty and low-information regions retain proper IW/Dirichlet updates.
empty_fit <- mtblr_regional_sampler(
  f$X[, 1L, drop = FALSE], f$Y, f$Ve, "A", c("A", "empty"),
  patterns, rep(1, 4), nu0, Psi0, 1600L, 300L, 1901L,
  "completed_active", "regional"
)
empty_mean <- apply(empty_fit$Vb[, , 2L, ], c(1L, 2L), mean)
expect_close(empty_mean, mt_iw_mean(nu0, Psi0), .06, "empty-region prior mean")
stopifnot(all(is.finite(empty_fit$Pi)),
          all(abs(apply(empty_fit$Pi, c(1L, 2L), sum) - 1) < 1e-12))

# Source reconstruction: set order affects the hybrid transition.
hybrid_ab <- mtblr_pattern_current_hybrid(
  f$X, f$Y, f$Ve, patterns, rep(.25, 4), rep(1, 4),
  list(1L, 2L), nu0, Psi0, 240L, 40L, 2001L
)
hybrid_ba <- mtblr_pattern_current_hybrid(
  f$X, f$Y, f$Ve, patterns, rep(.25, 4), rep(1, 4),
  list(2L, 1L), nu0, Psi0, 240L, 40L, 2001L
)
stopifnot(!identical(hybrid_ab$state, hybrid_ba$state),
          all(is.finite(hybrid_ab$Vb)), all(is.finite(hybrid_ab$Vb_latent)))

cat("All four-pattern and regional MTBLR reference cases passed.\n")
