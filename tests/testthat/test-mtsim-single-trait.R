make_mtsim_test_matrix <- function() {
 set.seed(1)
 W <- matrix(
  rbinom(80L * 40L, size = 2L, prob = 0.3),
  nrow = 80L,
  ncol = 40L
 )
 rownames(W) <- paste0("id", seq_len(nrow(W)))
 colnames(W) <- paste0("m", seq_len(ncol(W)))
 scale(W)
}

test_that("mtsim preserves single-trait covariance dimensions", {
 W <- make_mtsim_test_matrix()

 sim <- mtsim(
  W = W,
  standardize_W = FALSE,
  nt = 1L,
  n_shared = 10L,
  n_specific = 0L,
  h2 = 0.2,
  re = 0,
  seed = 1001L
 )

 expect_equal(dim(sim$B), c(ncol(W), 1L))
 expect_equal(dim(sim$G), c(nrow(W), 1L))
 expect_equal(dim(sim$E), c(nrow(W), 1L))
 expect_equal(dim(sim$Sigma_e), c(1L, 1L))
 expect_true(all(is.finite(sim$B)))
 expect_true(all(is.finite(sim$G)))
 expect_true(all(is.finite(sim$E)))
 expect_true(all(is.finite(sim$Sigma_e)))
 expect_equal(unname(sim$G), unname(W %*% sim$B), tolerance = 1e-12)
 expect_identical(rownames(sim$B), colnames(W))
 expect_identical(rownames(sim$G), rownames(W))
 expect_identical(rownames(sim$E), rownames(W))
 expect_identical(sim$rsids, colnames(W))
 expect_identical(sim$ids, rownames(W))
})

test_that("mtsim supports multiple single-trait heritabilities", {
 W <- make_mtsim_test_matrix()
 h2_values <- c(0.2, 0.3, 0.5, 0.8)

 for (i in seq_along(h2_values)) {
  h2 <- h2_values[i]
  sim <- mtsim(
   W = W,
   standardize_W = FALSE,
   nt = 1L,
   n_shared = 10L,
   n_specific = 0L,
   h2 = h2,
   re = 0,
   seed = 2000L + i
  )

  expect_equal(dim(sim$Sigma_e), c(1L, 1L), info = paste("h2 =", h2))
  expect_true(all(is.finite(c(sim$B, sim$G, sim$E, sim$Sigma_e))))
  expect_equal(unname(sim$G), unname(W %*% sim$B), tolerance = 1e-12)
  expect_equal(sim$Sigma_e[1, 1], (1 - h2) / h2, tolerance = 1e-12)
 }
})

test_that("mtsim preserves multi-trait residual covariance behaviour", {
 W <- make_mtsim_test_matrix()

 sim <- mtsim(
  W = W,
  standardize_W = FALSE,
  nt = 2L,
  n_shared = 10L,
  n_specific = 2L,
  h2 = c(0.3, 0.6),
  rg = 0.2,
  re = 0.25,
  seed = 3001L
 )

 expect_equal(dim(sim$B), c(ncol(W), 2L))
 expect_equal(dim(sim$G), c(nrow(W), 2L))
 expect_equal(dim(sim$E), c(nrow(W), 2L))
 expect_equal(dim(sim$Sigma_e), c(2L, 2L))
 expect_true(all(is.finite(c(sim$B, sim$G, sim$E, sim$Sigma_e))))
 expect_equal(sim$Sigma_e, t(sim$Sigma_e), tolerance = 1e-12)
 expect_true(all(eigen(sim$Sigma_e, symmetric = TRUE, only.values = TRUE)$values > 0))
 expect_equal(unname(sim$G), unname(W %*% sim$B), tolerance = 1e-12)
})

test_that("mtsim retains scalar re semantics for one trait", {
 W <- make_mtsim_test_matrix()

 sim <- mtsim(
  W = W,
  standardize_W = FALSE,
  nt = 1L,
  n_shared = 5L,
  n_specific = 0L,
  h2 = 0.5,
  re = 0,
  seed = 4001L
 )

 expect_equal(dim(sim$re_target), c(1L, 1L))
 expect_equal(sim$re_target[1, 1], 1)
 expect_equal(dim(sim$Sigma_e), c(1L, 1L))
})
