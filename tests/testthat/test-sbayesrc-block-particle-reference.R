test_that("particle marker proposal is the exact incremental conditional", {
  set.seed(41)
  score <- 0.7
  diagonal <- 2.3
  ve <- 1.2
  vb <- 0.8
  gamma <- c(0, 0.2, 1)
  probability <- c(0.75, 0.2, 0.05)
  value <- sblr:::.sbayesrc_particle_marker_parameters(
    score, diagonal, ve, vb, gamma, probability)
  expect_equal(sum(exp(value$log_weight - value$log_normalizer)), 1,
               tolerance = 1e-14)
  expect_equal(value$variance[2], 1 / (diagonal / ve + 1 / (vb * gamma[2])))
  expect_equal(value$mean[3], value$variance[3] * score / ve)
})

test_that("conditional particle Gibbs preserves a tiny exact block target", {
  skip_on_cran()
  Q <- matrix(c(1, 0.2, -0.1, 0.7, 0.8, 0.3), nrow = 2)
  w <- c(0.5, -0.2)
  probability <- rbind(c(.75, .25), c(.65, .35), c(.8, .2))
  gamma <- c(0, 1)
  exact <- sblr:::.sbayesrc_exact_block_allocation(
    Q, w, probability, gamma, vb = 0.6, ve = 1.1)
  expect_equal(sum(exact$probability), 1, tolerance = 1e-12)

  set.seed(20260807)
  component <- integer(3)
  beta <- numeric(3)
  draws <- matrix(0L, 12000, 3)
  for (iteration in seq_len(nrow(draws))) {
    update <- sblr:::.sbayesrc_particle_block_step(
      Q, w, component, beta, probability, gamma, vb = 0.6, ve = 1.1,
      particles = 16L, resampling_threshold = 0.5)
    component <- update$component
    beta <- update$beta
    draws[iteration, ] <- component
  }
  empirical_pip <- colMeans(draws[-seq_len(2000), , drop = FALSE] > 0L)
  exact_pip <- vapply(seq_len(3), function(marker)
    sum(exact$probability[exact$states[, marker] > 0L]), numeric(1L))
  expect_equal(empirical_pip, exact_pip, tolerance = 0.025)
})

test_that("particle diagnostics do not consume random numbers", {
  Q <- diag(3)
  w <- c(.2, -.1, .3)
  probability <- matrix(c(.8, .2), 3, 2, byrow = TRUE)
  run <- function(diagnostics) {
    set.seed(91)
    sblr:::.sbayesrc_particle_block_step(
      Q, w, integer(3), numeric(3), probability, c(0, 1), 1, 1,
      particles = 8L, retain_diagnostics = diagnostics)
  }
  expect_identical(run(FALSE), run(TRUE)[c("component", "beta")])
})
