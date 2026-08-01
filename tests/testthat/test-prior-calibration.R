test_that("global BayesC calibration preserves the legacy formula", {
 m <- 200L; pi0 <- 0.013; vy <- c(T1 = 2.5, T2 = 4); h2 <- c(.2, .6)
 got <- sblr:::.stblr_scalar_prior_calibration(
  vy, h2, 4, 4, rep(pi0, m), rep(pi0, m))
 expect_equal(diag(got$B), vy * h2 / (m * pi0))
 expect_equal(diag(got$B) * got$calibration$prior_weight_initial, vy * h2)
})

test_that("BayesR separates initial pi from the Dirichlet prior mean", {
 grids <- list(c(0, .01, .1, 1), c(0, .2, 2, 5))
 for (gamma in grids) {
  pi <- c(.7, .1, .1, .1)
  alpha <- c(20, 5, 3, 2)
  got <- sblr:::.make_stblr_bayesr_priors(
   vy = 3, m = 50, h2 = .4, nub = 4, nue = 4, pi = pi,
   mixture_var = gamma, alpha = alpha)
  w0 <- 50 * sum((pi / sum(pi)) * gamma)
  wp <- 50 * sum((alpha / sum(alpha)) * gamma)
  expect_equal(unname(got$calibration$prior_weight_initial), w0)
  expect_equal(unname(got$calibration$prior_weight_prior_mean), wp)
  expect_equal(unname(diag(got$B) * w0), 3 * .4)
  expect_equal(unname(diag(got$ssb_prior) * wp), .5 * 3 * .4)
 }
})

test_that("resolved marker probabilities, variance multipliers, and MAF scales combine", {
 p <- c(.01, .2, .4, .8)
 v <- c(.5, 2, 1, 3)
 q <- c(.2, .4, .8, 1.6)
 got <- sblr:::.stblr_scalar_prior_calibration(
  5, .3, 4, 4, p, p, marker_scale = q, variance_multiplier = v)
 direct <- sum(p * v * q)
 expect_equal(unname(got$calibration$prior_weight_initial), direct)
 expect_equal(unname(diag(got$B) * direct), 5 * .3)

 ordinary <- sblr:::.stblr_scalar_prior_calibration(5, .3, 4, 4, p, p)
 ones <- sblr:::.stblr_scalar_prior_calibration(
  5, .3, 4, 4, p, p, marker_scale = rep(1, 4))
 expect_equal(ordinary$B, ones$B)
})

test_that("sampled MAF-S initialization uses the requested initial S once", {
 h <- c(.1, .2, .4)
 info <- list(fixed = FALSE, log_h = log(h), prior_scale = numeric())
 expect_equal(
  sblr:::.stblr_calibration_maf_scale(info, TRUE, -.5, 3), h^.5)
 expect_equal(
  sblr:::.stblr_calibration_maf_scale(info, FALSE, 0, 3), rep(1, 3))
})

test_that("BayesRC resolved component probabilities include gamma", {
 A <- cbind(intercept = 1, x = c(-1, 0, 1))
 gamma <- c(0, .1, 1)
 alpha <- matrix(c(0, 0, -.5, .7), 2, 2)
 probability <- sbayesrc_marker_pi(A, alpha, gamma)
 expected_gamma <- as.numeric(probability %*% gamma)
 got <- sblr:::.stblr_scalar_prior_calibration(
  2, .25, 4, 4, expected_gamma, expected_gamma)
 expect_equal(unname(got$calibration$prior_weight_initial), sum(expected_gamma))
 expect_equal(unname(diag(got$B) * sum(expected_gamma)), .5)
})

test_that("group and learned annotation weights use resolved marker values", {
 group <- c(1L, 1L, 2L, 2L)
 group_pi <- c(.1, .6); group_v <- c(2, .5)
 expected <- group_pi[group]
 multiplier <- group_v[group]
 group_fit <- sblr:::.stblr_scalar_prior_calibration(
  4, .5, 4, 4, expected, expected, variance_multiplier = multiplier)
 expect_equal(unname(group_fit$calibration$prior_weight_initial),
              sum(expected * multiplier))

 A <- cbind(1, c(-1, 0, 1, 2))
 resolved <- sblr:::.stblr_make_prior_from_annotations(
  A, 1, .2, beta_pi = c(0, .5), beta_vb = c(0, -.3))
 learned <- sblr:::.stblr_scalar_prior_calibration(
  4, .5, 4, 4, resolved$pi_marker, resolved$pi_marker,
  variance_multiplier = resolved$vb_multiplier)
 expect_equal(unname(learned$calibration$prior_weight_initial),
  sum(resolved$pi_marker[[1]] * resolved$vb_multiplier[[1]]))
})

test_that("MT BayesC uses trait-pattern-specific expected weights", {
 patterns <- rbind(null = c(0, 0), t1 = c(1, 0), both = c(1, 1))
 probability <- c(.5, .3, .2)
 got <- sblr:::.mtblr_prior_calibration(
  c(2, 5), c(.25, .4), 4, 4, patterns, probability,
  marker_scale = rep(1, 10))
 direct <- matrix(0, 2, 2)
 for (j in 1:10) for (s in seq_along(probability))
  direct <- direct + probability[s] * tcrossprod(patterns[s, ])
 expect_equal(got$calibration$prior_weight_initial, direct)
 expect_equal(diag(got$vb) * diag(direct), c(.5, 2))
})

test_that("MT BayesR weights joint states, components, and marker scales", {
 patterns <- rbind(c(0, 0), c(1, 0), c(1, 1), c(0, 1))
 probability <- c(.6, .1, .2, .1)
 gamma <- c(0, .1, 1, 2)
 q <- c(.2, .5, 1)
 got <- sblr:::.mtblr_prior_calibration(
  c(3, 4), c(.2, .5), 4, 4, patterns, probability,
  gamma = gamma, marker_scale = q)
 direct <- matrix(0, 2, 2)
 for (j in seq_along(q)) for (s in seq_along(probability))
  direct <- direct + probability[s] * gamma[s] * q[j] *
   tcrossprod(patterns[s, ])
 expect_equal(got$calibration$prior_weight_initial, direct)
 expect_equal(diag(got$vb) * diag(direct), c(.6, 2))
})

test_that("fixed BayesRC probabilities reduce algebraically to BayesR", {
 patterns <- rbind(c(0, 0), c(1, 0), c(1, 1), c(0, 1))
 p <- c(.7, .1, .1, .1); gamma <- c(0, .1, 1, .5)
 bayesr <- sblr:::.mtblr_prior_calibration(
  c(2, 2), c(.3, .4), 4, 4, patterns, p, gamma = gamma,
  marker_scale = rep(1, 5))
 bayesrc <- sblr:::.mtblr_prior_calibration(
  c(2, 2), c(.3, .4), 4, 4, patterns,
  matrix(rep(p, each = 5), 5, 4), gamma = gamma,
  marker_scale = rep(1, 5))
 expect_equal(bayesrc$vb, bayesr$vb)
 expect_equal(bayesrc$ssb_prior, bayesr$ssb_prior)
})

test_that("explicit variance overrides remain authoritative", {
 B <- matrix(2); E <- matrix(3); ssb <- matrix(4); sse <- matrix(5)
 got <- sblr:::.stblr_scalar_prior_calibration(
  2, .3, 4, 4, rep(.1, 4), B = B, E = E,
  ssb_prior = ssb, sse_prior = sse)
 expect_equal(unname(got$B), B); expect_equal(unname(got$E), E)
 expect_equal(unname(got$ssb_prior), ssb); expect_equal(unname(got$sse_prior), sse)

 vb <- diag(c(.2, .3)); prior <- diag(c(.1, .2))
 mt <- sblr:::.mtblr_prior_calibration(
  c(2, 3), c(.2, .3), 4, 4,
  rbind(c(0, 0), c(1, 0), c(0, 1)), c(.8, .1, .1),
  marker_scale = rep(1, 5), vb = vb, ssb_prior = prior)
 expect_equal(mt$vb, vb); expect_equal(mt$ssb_prior, prior)
})

test_that("unsafe automatic full MT covariance requires explicit vb", {
 patterns <- rbind(c(0, 0), c(1, 0), c(0, 1), c(1, 1))
 expect_error(sblr:::.mtblr_prior_calibration(
  c(2, 3), c(.2, .3), 4, 4, patterns, rep(.25, 4),
  marker_scale = rep(1, 3), vg = matrix(c(1, .2, .2, 1), 2)),
  "supply an explicit")
})
