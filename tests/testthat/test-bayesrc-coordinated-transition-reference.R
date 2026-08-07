test_that("coordinated reference reproduces ordered probit sticks", {
  annotation <- cbind(1, c(-1, 0, 1))
  alpha <- cbind(c(-0.4, 0.7), c(0.2, -0.3), c(0.8, 0.1))
  probability <- sblr:::.bayesrc_coordinated_component_probabilities(
    annotation, alpha)
  continuation <- pnorm(annotation %*% alpha)
  expected <- cbind(
    1 - continuation[, 1L],
    continuation[, 1L] * (1 - continuation[, 2L]),
    continuation[, 1L] * continuation[, 2L] *
      (1 - continuation[, 3L]),
    apply(continuation, 1L, prod))
  expect_equal(probability, expected, tolerance = 1e-14)
  expect_equal(rowSums(probability), rep(1, nrow(annotation)))

  component <- 0:3
  indicators <- outer(component, 0:2, `>`)
  eligible <- sapply(0:2, function(stick) stick == 0L | component >= stick)
  expect_equal(indicators[, 1L], component > 0L)
  expect_equal(eligible[, 2L], component >= 1L)
  expect_equal(eligible[, 3L], component >= 2L)
})

test_that("beta-fixed allocation refresh cannot cross spike support", {
  support <- sblr:::.bayesrc_coordinated_beta_fixed_support
  expect_true(support(c(0L, 1L, 2L), c(0, 0.2, -0.1)))
  expect_false(support(c(0L, 1L), c(0.2, 0.1)))
  expect_false(support(c(1L, 0L), c(0, 0)))
  expect_true(support(c(1L, 2L), c(0.1, -0.3)))
})

test_that("collapsed subset weights match independent one-dimensional integration", {
  score <- 1.3
  diagonal <- 2.4
  ve <- 0.8
  vb <- 0.6
  gamma <- c(0, 0.2, 1.5)
  marker_probability <- matrix(c(0.7, 0.2, 0.1), nrow = 1L)
  state <- sblr:::.bayesrc_coordinated_subset_states(
    score, matrix(diagonal, 1L, 1L), ve, vb, gamma, marker_probability)

  numerical <- numeric(length(gamma))
  numerical[[1L]] <- marker_probability[[1L]]
  for (component in 2:length(gamma)) {
    prior_sd <- sqrt(vb * gamma[[component]])
    integrand <- function(beta) {
      exp(-(diagonal * beta^2 - 2 * score * beta) / (2 * ve)) *
        dnorm(beta, 0, prior_sd)
    }
    numerical[[component]] <- marker_probability[[component]] *
      integrate(integrand, -Inf, Inf, rel.tol = 1e-12)$value
  }
  numerical <- numerical / sum(numerical)
  expect_equal(state$probability, numerical, tolerance = 1e-10)
  expect_equal(state$component[, 1L], 0:2)
})

test_that("coordinated collapsed MH ratio satisfies detailed balance", {
  annotation <- cbind(1, c(-1, 1, -0.5, 0.7))
  alpha_old <- cbind(c(-0.5, 0.4), c(0.1, -0.2))
  alpha_new <- cbind(c(-0.2, 0.8), c(0.1, -0.2))
  probability_old <- sblr:::.bayesrc_coordinated_component_probabilities(
    annotation, alpha_old)
  probability_new <- sblr:::.bayesrc_coordinated_component_probabilities(
    annotation, alpha_new)
  score <- c(0.9, -0.6)
  operator <- matrix(c(2.0, 0.35, 0.35, 1.7), 2L)
  gamma <- c(0, 0.3, 1.2)
  ve <- 0.9
  vb <- 0.7
  old_subset <- sblr:::.bayesrc_coordinated_subset_states(
    score, operator, ve, vb, gamma, probability_old[1:2, ])
  new_subset <- sblr:::.bayesrc_coordinated_subset_states(
    score, operator, ve, vb, gamma, probability_new[1:2, ])
  outside_component <- c(0L, 2L)
  log_ratio <- sblr:::.bayesrc_coordinated_log_mh(
    alpha_old[, 1L], alpha_new[, 1L], outside_component,
    annotation[3:4, , drop = FALSE], alpha_old, alpha_new,
    old_subset, new_subset, intercept_mean = -0.4, intercept_sd = 1,
    sigma_sq_alpha = 1.3)

  old_component <- c(0L, 1L, outside_component)
  new_component <- c(2L, 0L, outside_component)
  old_beta <- c(0, 0.25, 0, -0.15)
  new_beta <- c(-0.2, 0, 0, -0.15)
  full_operator <- matrix(0, 4L, 4L)
  full_operator[1:2, 1:2] <- operator
  diag(full_operator)[3:4] <- c(1.5, 1.8)
  full_score <- c(score, 0.2, -0.4)

  log_target <- function(alpha, component, beta, probability) {
    active <- which(component > 0L)
    if (any(component == 0L & beta != 0)) return(-Inf)
    likelihood <- -(drop(crossprod(beta, full_operator %*% beta)) -
      2 * drop(crossprod(full_score, beta))) / (2 * ve)
    effect_prior <- sum(vapply(active, function(i) {
      dnorm(beta[[i]], 0, sqrt(vb * gamma[component[[i]] + 1L]), log = TRUE)
    }, numeric(1L)))
    allocation <- sum(log(probability[
      cbind(seq_along(component), component + 1L)]))
    alpha_prior <- dnorm(alpha[[1L]], -0.4, 1, log = TRUE) +
      dnorm(alpha[[2L]], 0, sqrt(1.3), log = TRUE)
    likelihood + effect_prior + allocation + alpha_prior
  }
  log_mvn <- function(value, mean, covariance, active) {
    if (!length(active)) return(0)
    difference <- value[active] - mean[active]
    chol_covariance <- chol(covariance[active, active, drop = FALSE])
    -0.5 * (length(active) * log(2 * pi) +
      2 * sum(log(diag(chol_covariance))) +
      sum(backsolve(chol_covariance, difference, transpose = TRUE)^2))
  }
  proposal_log_density <- function(subset, component, beta) {
    row <- which(rowSums(abs(sweep(subset$component, 2L,
                                  as.integer(component), `-`))) == 0)
    expect_length(row, 1L)
    active <- which(component > 0L)
    log(subset$probability[[row]]) +
      log_mvn(beta, subset$mean[[row]], subset$covariance[[row]], active)
  }
  log_q_forward <- proposal_log_density(
    new_subset, new_component[1:2], new_beta[1:2])
  log_q_reverse <- proposal_log_density(
    old_subset, old_component[1:2], old_beta[1:2])
  direct_ratio <- log_target(alpha_new[, 1L], new_component, new_beta,
                             probability_new) -
    log_target(alpha_old[, 1L], old_component, old_beta, probability_old) +
    log_q_reverse - log_q_forward
  expect_equal(log_ratio, direct_ratio, tolerance = 1e-10)

  old_log_target <- log_target(alpha_old[, 1L], old_component, old_beta,
                               probability_old)
  new_log_target <- log_target(alpha_new[, 1L], new_component, new_beta,
                               probability_new)
  log_flow_forward <- old_log_target + log_q_forward + min(0, log_ratio)
  log_flow_reverse <- new_log_target + log_q_reverse + min(0, -log_ratio)
  expect_equal(log_flow_forward, log_flow_reverse, tolerance = 1e-10)
})

test_that("coordinated reference covers boundary states and is RNG neutral", {
  probability <- matrix(rep(c(0.8, 0.15, 0.05), 3L), 3L, byrow = TRUE)
  seed_before <- { set.seed(9182); .Random.seed }
  state <- sblr:::.bayesrc_coordinated_subset_states(
    c(0, 0, 0), diag(3), 1, 0.5, c(0, 0.1, 1), probability)
  expect_identical(.Random.seed, seed_before)
  expect_true(any(rowSums(state$component > 0L) == 0L))
  expect_true(any(rowSums(state$component > 0L) == 3L))
  expect_true(all(is.finite(state$log_weight)))
  expect_equal(sum(state$probability), 1, tolerance = 1e-14)
})
