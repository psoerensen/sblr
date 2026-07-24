source("tests/testthat/helper-mtblr-bed-contract.R")

bench <- function(n, m, nt, repetitions = 20L) {
  set.seed(17014L + n + m + nt)
  dosage <- matrix(sample(c(0, 1, 2, NA), n * m, replace = TRUE,
                          prob = c(.28, .4, .28, .04)), n, m)
  af <- runif(m, .1, .9)
  X <- phase17n_transform_genotypes(dosage, af, TRUE)
  R <- matrix(rnorm(n * nt), n, nt)
  j <- max(1L, m %/% 2L)
  packed_bytes <- m * ceiling(n / 4)
  dense_bytes <- as.numeric(object.size(X))
  map_bytes <- 40 * m

  decode <- system.time(for (i in seq_len(repetitions))
    phase17n_transform_genotypes(dosage[, j, drop = FALSE], af[j], TRUE))[["elapsed"]]
  score_workspace <- system.time(for (i in seq_len(repetitions))
    crossprod(X[, j], R))[["elapsed"]]
  update_workspace <- system.time(for (i in seq_len(repetitions))
    phase17n_update_residual(R, X[, j], rep(.01, nt)))[["elapsed"]]
  direct_lookup_estimate <- n * nt * 2

  data.frame(
    n = n, m = m, nt = nt, repetitions = repetitions,
    packed_bytes = packed_bytes, dense_X_bytes = dense_bytes,
    marker_map_bytes = map_bytes, decoded_workspace_bytes = 8 * n,
    marker_decode_seconds = decode,
    marker_score_seconds = score_workspace,
    marker_update_seconds = update_workspace,
    direct_packed_lookup_count_per_score_update = direct_lookup_estimate
  )
}

results <- rbind(
  small_nt1 = bench(250L, 500L, 1L),
  small_nt4 = bench(250L, 500L, 4L),
  moderate_nt2 = bench(2000L, 2000L, 2L),
  moderate_nt6 = bench(2000L, 2000L, 6L)
)
print(results)
cat("BENCHMARK_SCOPE=audit_oracle_only\n")
cat("TIMING_CLAIM=none_tiny_synthetic_signal_only\n")
