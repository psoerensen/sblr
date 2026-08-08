test_that("tiny global selection-model block preserves the exact fixed-z target", {
  fixture <- .sbs_fixture()
  exact <- .sbs_exact_posterior(
    fixture$z, fixture$annotation, fixture$eligible,
    fixture$pi_a, fixture$tau2, fixture$intercept_prior
  )
  transition <- .sbs4c_global_model_transition(exact$model_probability)
  expect_equal(rowSums(transition), rep(1, nrow(transition)), tolerance = 1e-14)
  expect_equal(
    drop(exact$model_probability %*% transition),
    as.numeric(exact$model_probability), tolerance = 1e-14
  )
  expect_lte(
    .sbs4c_detailed_balance_error(exact$model_probability, transition),
    1e-14
  )

  set.seed(20271301L)
  draws <- sample.int(
    length(exact$model_probability), 25000L, replace = TRUE,
    prob = exact$model_probability
  )
  frequency <- tabulate(draws, nbins = length(exact$model_probability)) / 25000
  expect_lte(max(abs(frequency - exact$model_probability)), 0.01)
})
