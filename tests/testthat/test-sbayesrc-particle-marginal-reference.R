test_that("particle likelihood estimator is unbiased on a tiny exact block", {
  Q <- matrix(c(1, .2, -.1, .7, .8, .3), nrow = 2)
  w <- c(.5, -.2)
  probability <- rbind(c(.75, .25), c(.65, .35), c(.8, .2))
  exact <- sblr:::.sbayesrc_exact_block_allocation(
    Q, w, probability, c(0, 1), vb = .6, ve = 1.1)
  set.seed(8172)
  estimate <- replicate(1500, exp(sblr:::.sbayesrc_sis_block_likelihood(
    Q, w, probability, c(0, 1), .6, 1.1, particles = 8L,
    retain_paths = FALSE)$log_likelihood - exact$log_normalizer))
  standard_error <- stats::sd(estimate) / sqrt(length(estimate))
  expect_lt(abs(mean(estimate) - 1), 4 * standard_error + 0.005)
})

test_that("selected paths recover exact marginals under the extended target", {
  Q <- matrix(c(1, .2, -.1, .7, .8, .3), nrow = 2)
  w <- c(.5, -.2)
  probability <- rbind(c(.75, .25), c(.65, .35), c(.8, .2))
  exact <- sblr:::.sbayesrc_exact_block_allocation(
    Q, w, probability, c(0, 1), vb = .6, ve = 1.1)
  set.seed(9821)
  out <- replicate(2000, {
    value <- sblr:::.sbayesrc_sis_block_likelihood(
      Q, w, probability, c(0, 1), .6, 1.1, particles = 8L)
    c(weight = exp(value$log_likelihood - exact$log_normalizer),
      value$component > 0L)
  })
  empirical <- rowSums(out[-1L, , drop = FALSE] * out[1L, ]) / sum(out[1L, ])
  expect_lt(max(abs(as.numeric(empirical) - unname(exact$pip))), 0.025)
})

test_that("correlated auxiliary proposal preserves standard normal marginals", {
  set.seed(771)
  auxiliary <- sblr:::.sbayesrc_particle_auxiliary(32L, 20L)
  proposed <- sblr:::.sbayesrc_correlate_auxiliary(auxiliary, .9)
  expect_equal(dim(proposed$component), c(32L, 20L))
  expect_true(all(is.finite(unlist(proposed))))
  expect_error(sblr:::.sbayesrc_correlate_auxiliary(auxiliary, 1),
               "strictly between")
})

test_that("particle marginal reference is reproducible and diagnostics-neutral", {
  Q <- diag(4)
  w <- c(.1, -.2, .3, .05)
  probability <- matrix(c(.8, .2), 4, 2, byrow = TRUE)
  set.seed(54)
  auxiliary <- sblr:::.sbayesrc_particle_auxiliary(8L, 4L)
  first <- sblr:::.sbayesrc_sis_block_likelihood(
    Q, w, probability, c(0, 1), 1, 1, 8L, auxiliary)
  second <- sblr:::.sbayesrc_sis_block_likelihood(
    Q, w, probability, c(0, 1), 1, 1, 8L, auxiliary)
  expect_identical(first, second)
})

test_that("particle marginal alpha ratio includes the complete proper prior", {
  alpha <- matrix(c(-1, .2, -.1, -.5, .3, .1), 3, 2)
  proposed <- alpha + .05
  value <- sblr:::.sbayesrc_particle_marginal_log_ratio(
    alpha, proposed, -10, -9.5, c(-1, -.4), c(1, 1), c(.8, 1.2))
  reference <- sblr:::.sbayesrc_alpha_log_prior(
    proposed, c(-1, -.4), c(1, 1), c(.8, 1.2)) -
    sblr:::.sbayesrc_alpha_log_prior(
      alpha, c(-1, -.4), c(1, 1), c(.8, 1.2)) + .5
  expect_equal(value, reference, tolerance = 1e-14)
})

test_that("independent pseudo-marginal MH has the exact tiny alpha marginal", {
  Q <- matrix(c(1, .2, -.1, .7, .8, .3), nrow = 2)
  w <- c(.5, -.2)
  probabilities <- list(
    rbind(c(.85, .15), c(.8, .2), c(.9, .1)),
    rbind(c(.65, .35), c(.55, .45), c(.75, .25)))
  exact_log <- vapply(probabilities, function(probability)
    sblr:::.sbayesrc_exact_block_allocation(
      Q, w, probability, c(0, 1), .6, 1.1)$log_normalizer, numeric(1L))
  exact_second <- stats::plogis(exact_log[[2L]] - exact_log[[1L]])
  set.seed(1771)
  state <- 1L
  current <- sblr:::.sbayesrc_sis_block_likelihood(
    Q, w, probabilities[[state]], c(0, 1), .6, 1.1, 8L)
  trace <- integer(9000L)
  for (iteration in seq_along(trace)) {
    proposed_state <- 3L - state
    proposed <- sblr:::.sbayesrc_sis_block_likelihood(
      Q, w, probabilities[[proposed_state]], c(0, 1), .6, 1.1, 8L)
    if (log(stats::runif(1L)) < min(0, proposed$log_likelihood -
                                    current$log_likelihood)) {
      state <- proposed_state
      current <- proposed
    }
    trace[[iteration]] <- state
  }
  expect_lt(abs(mean(trace[-seq_len(1000L)] == 2L) - exact_second), 0.04)
})
