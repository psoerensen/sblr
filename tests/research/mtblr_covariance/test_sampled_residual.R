script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
research_dir <- if (length(script_arg)) {
  dirname(normalizePath(sub("^--file=", "", script_arg[1L]), mustWork = TRUE))
} else normalizePath("tests/research/mtblr_covariance", mustWork = TRUE)
for (file in c("mtblr_reference_model.R", "mtblr_pattern_reference.R",
               "mtblr_pattern_samplers.R", "mtblr_sampled_residual.R")) {
  source(file.path(research_dir, file), local = FALSE)
}

residual <- cbind(c(-.8, .2, .5, -.1, .4, -.3),
                  c(.3, -.4, .6, .2, -.2, .5))
prior_df <- 1.5
prior_scale <- matrix(c(.7, .14, .14, .6), 2L)
posterior <- mt_residual_covariance_posterior(
  residual, prior_df, prior_scale)
stopifnot(
  identical(posterior$degrees_of_freedom, prior_df + nrow(residual)),
  isTRUE(all.equal(posterior$statistic, crossprod(residual), tolerance = 0)),
  isTRUE(all.equal(posterior$scale,
                   prior_scale + crossprod(residual), tolerance = 0)))

# The prior is proper but has no finite mean; six observations make the
# posterior mean finite under the same degrees-of-freedom/scale convention.
set.seed(741L)
conditional <- mt_rsample_residual_covariance(
  residual, prior_df, prior_scale, 12000L)
expected_mean <- posterior$scale / (posterior$degrees_of_freedom - 3)
observed_mean <- apply(conditional$draws, c(1L, 2L), mean)
if (max(abs(observed_mean - expected_mean)) > .035) {
  stop("Independent residual inverse-Wishart mean exceeded tolerance.",
       call. = FALSE)
}

patterns <- mt_pattern_space()
X <- cbind(
  c(-1.2, -.4, .5, 1.1, .2, -.7),
  c(.8, -.5, -.2, .4, 1.0, -1.1))
Y <- cbind(
  .6 * X[, 1L] + c(.2, -.1, .3, -.2, .1, -.3),
  -.35 * X[, 2L] + c(.1, .2, -.15, .25, -.2, .05))
Vb0 <- matrix(c(.35, .08, .08, .30), 2L)
Ve0 <- matrix(c(.80, .18, .18, .70), 2L)
fit <- mtblr_pattern_sampler_sampled_residual(
  X, Y, patterns, rep(1, 4L), 4.5, diag(c(.5, .45)),
  1.5, matrix(c(.7, .1, .1, .65), 2L), Vb0, Ve0,
  n_iter = 1200L, burn = 200L, seed = 742L)
stopifnot(
  all(is.finite(fit$Ve)),
  all(is.finite(fit$Vb)),
  max(abs(rowSums(fit$Pi) - 1)) < 1e-12,
  isTRUE(all.equal(crossprod(fit$final$residual),
                   crossprod(Y - X %*% fit$final$realised),
                   tolerance = 1e-12)))
for (draw in seq_len(dim(fit$Ve)[3L])) {
  mt_assert_spd(fit$Ve[, , draw], "sampled residual covariance")
}

cat("Sampled residual-covariance research checks passed.\n")
