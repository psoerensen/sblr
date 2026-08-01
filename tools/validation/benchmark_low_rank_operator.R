#!/usr/bin/env Rscript

sizes <- c(250L, 500L, 1000L, 2000L)
set.seed(5119)

elapsed <- function(fun, sweeps = 20L) {
  system.time(for (iteration in seq_len(sweeps)) fun())[["elapsed"]] / sweeps
}

rows <- lapply(sizes, function(m) {
  k <- ceiling(0.25 * m)
  delta <- rnorm(m, sd = 0.002)
  dense <- diag(m)
  dense[row(dense) == col(dense) + 1L] <- 0.15
  dense[col(dense) == row(dense) + 1L] <- 0.15
  Q_full <- matrix(rnorm(m * m, sd = 1 / sqrt(m)), m, m)
  Q_low <- Q_full[seq_len(k), , drop = FALSE]

  dense_sweep <- function() {
    residual <- numeric(m)
    for (marker in seq_len(m)) residual <- residual - dense[, marker] * delta[marker]
    residual
  }
  low_sweep <- function(Q) {
    residual <- numeric(nrow(Q))
    for (marker in seq_len(m)) residual <- residual - Q[, marker] * delta[marker]
    residual
  }
  csr_sweep <- function() {
    residual <- numeric(m)
    for (marker in seq_len(m)) {
      index <- max(1L, marker - 1L):min(m, marker + 1L)
      residual[index] <- residual[index] - dense[index, marker] * delta[marker]
    }
    residual
  }
  eigen_time <- system.time(eigen(dense, symmetric = TRUE, only.values = TRUE))[["elapsed"]]
  data.frame(
    block_size = m, retained_rank = k,
    csr_seconds_per_sweep = elapsed(csr_sweep),
    dense_seconds_per_sweep = elapsed(dense_sweep),
    low_rank_full_seconds_per_sweep = elapsed(function() low_sweep(Q_full)),
    low_rank_25pct_seconds_per_sweep = elapsed(function() low_sweep(Q_low)),
    eigendecomposition_seconds = eigen_time,
    csr_storage_bytes = 8 * (m + 1) + 8 * m + 16 * (m - 1),
    dense_packed_storage_bytes = 4 * m * (m + 1) / 2,
    low_rank_full_storage_bytes = 4 * m * m + 8 * m + 8 * m,
    low_rank_25pct_storage_bytes = 4 * m * k + 8 * k + 8 * m,
    dense_chain_state_bytes = 8 * m,
    low_rank_full_chain_state_bytes = 8 * m,
    low_rank_25pct_chain_state_bytes = 8 * k
  )
})
result <- do.call(rbind, rows)
print(result, row.names = FALSE)
cat("Performance microbenchmark complete; interpret full-rank and truncated ranks separately.\n")
