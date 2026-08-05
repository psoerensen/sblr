pairwise_reference <- function(prior_i, prior_j, gamma, vb, scale_i, scale_j,
                               ve, G, score) {
  states <- expand.grid(component_j = seq_along(gamma) - 1L,
                        component_i = seq_along(gamma) - 1L)
  states <- states[c("component_i", "component_j")]
  states$log_weight <- mapply(function(ci, cj) {
    active <- which(c(ci, cj) > 0L)
    out <- log(prior_i[ci + 1L]) + log(prior_j[cj + 1L])
    if (!length(active)) return(out)
    component <- c(ci, cj)[active]
    variance <- vb * c(scale_i, scale_j)[active] * gamma[component + 1L]
    precision <- G[active, active, drop = FALSE] / ve +
      diag(as.numeric(1 / variance), nrow = length(active))
    covariance <- solve(precision)
    h <- score[active] / ve
    out - 0.5 * (sum(log(variance)) + determinant(precision, logarithm = TRUE)$modulus) +
      0.5 * drop(crossprod(h, covariance %*% h))
  }, states$component_i, states$component_j)
  states$probability <- exp(states$log_weight - max(states$log_weight))
  states$probability <- states$probability / sum(states$probability)
  states
}

pairwise_native <- function(prior_i, prior_j, gamma, vb = 0.4,
                            scale_i = 1, scale_j = 1, ve = 0.7,
                            G = matrix(c(3, 0.8, 0.8, 2), 2),
                            score = c(1.2, -0.4)) {
  getFromNamespace(".st_bayesrc_pairwise_conditional", "sblr")(
    prior_i, prior_j, gamma, vb, scale_i, scale_j, ve,
    G[1, 1], G[2, 2], G[1, 2], score[1], score[2])
}

test_that("pairwise collapsed weights match dense Gaussian reference", {
  cases <- list(
    list(G = matrix(c(3, 0, 0, 2), 2), score = c(1.2, -0.4)),
    list(G = matrix(c(3, 0.8, 0.8, 2), 2), score = c(1.2, -0.4)),
    list(G = matrix(c(3, 2.2, 2.2, 2), 2), score = c(2, 1.5)),
    list(G = matrix(c(3, -1.1, -1.1, 2), 2), score = c(1.2, -0.4))
  )
  gamma <- c(0, 0.05, 0.5)
  prior_i <- c(0.8, 0.15, 0.05)
  prior_j <- c(0.6, 0.1, 0.3)
  for (case in cases) {
    native <- pairwise_native(prior_i, prior_j, gamma, vb = 0.4,
      scale_i = 0.7, scale_j = 1.8, G = case$G, score = case$score)
    reference <- pairwise_reference(prior_i, prior_j, gamma, 0.4, 0.7, 1.8,
      0.7, case$G, case$score)
    expect_equal(unname(native$component),
      unname(as.matrix(reference[c("component_i", "component_j")])))
    expect_equal(as.numeric(native$probability), reference$probability,
                 tolerance = 1e-12)
  }
})

test_that("zero cross-product reduces to independent marker conditionals", {
  gamma <- c(0, 0.1, 1)
  prior_i <- c(0.7, 0.2, 0.1)
  prior_j <- c(0.5, 0.3, 0.2)
  joint <- pairwise_native(prior_i, prior_j, gamma,
    G = diag(c(2.5, 4)), score = c(0.7, -1.1))
  matrix_probability <- matrix(joint$probability, nrow = length(gamma),
                               ncol = length(gamma))
  expect_equal(matrix_probability,
    tcrossprod(rowSums(matrix_probability), colSums(matrix_probability)),
    tolerance = 1e-12)
})

test_that("pairwise conditional has marker-exchange symmetry", {
  gamma <- c(0, 0.1, 1)
  prior <- c(0.8, 0.15, 0.05)
  out <- pairwise_native(prior, prior, gamma, scale_i = 1.3, scale_j = 1.3,
    G = matrix(c(3, 1.2, 1.2, 3), 2), score = c(0.8, 0.8))
  probability <- matrix(out$probability, nrow = 3, ncol = 3)
  expect_equal(probability, t(probability), tolerance = 1e-14)
})

test_that("pairwise conditional rejects indefinite operator pairs", {
  expect_error(pairwise_native(c(0.8, 0.2), c(0.8, 0.2), c(0, 1),
    G = matrix(c(1, 2, 2, 1), 2)), "positive semidefinite")
})
