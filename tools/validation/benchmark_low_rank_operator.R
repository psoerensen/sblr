#!/usr/bin/env Rscript

benchmark_library <- Sys.getenv("SBLR_BENCHMARK_LIBRARY")
if (nzchar(benchmark_library)) {
  .libPaths(c(benchmark_library, .libPaths()))
  suppressPackageStartupMessages(library(sblr))
} else {
  pkgload::load_all(".", compile = FALSE, quiet = TRUE)
}

markers <- 1000L
ranks <- c(1000L, 750L, 500L, 250L)
repetitions <- 7L
warmup <- 2L
seed <- 7301L

flatten_timings <- function(model, rank, benchmark) {
  do.call(rbind, lapply(names(benchmark$timings), function(operation) {
    timing <- benchmark$timings[[operation]]
    data.frame(
      model = model,
      markers = markers,
      retained_rank = rank,
      operation = operation,
      median_seconds = unname(timing[["median_seconds"]]),
      minimum_seconds = unname(timing[["minimum_seconds"]]),
      repetitions = repetitions,
      warmup = warmup,
      stringsAsFactors = FALSE
    )
  }))
}

native_rows <- lapply(ranks, function(rank) {
  bayesc <- sblr:::stblr_low_rank_bayesc_hot_path_benchmark_internal(
    markers, rank, repetitions, warmup, seed
  )
  bayesr <- sblr:::stblr_low_rank_bayesr_hot_path_benchmark_internal(
    markers, rank, repetitions, warmup, seed
  )
  rbind(
    flatten_timings("BayesC", rank, bayesc),
    flatten_timings("BayesR", rank, bayesr)
  )
})
native <- do.call(rbind, native_rows)

reference_speedup <- function(model, optimized, reference) {
  out <- merge(
    native[native$model == model & native$operation == optimized, ],
    native[native$model == model & native$operation == reference,
           c("retained_rank", "median_seconds")],
    by = "retained_rank", suffixes = c("", "_reference")
  )
  out$relative_speedup <- out$median_seconds_reference / out$median_seconds
  out
}

cat("Native retained-low-rank hot-path benchmark\n")
cat("R:", R.version.string, "\n")
cat("OS:", paste(Sys.info()[c("sysname", "release", "machine")],
                 collapse = " "), "\n")
cat("CPU logical cores:", parallel::detectCores(logical = TRUE), "\n")
cat("Compiler configuration:\n")
print(utils::capture.output(system("R CMD config CXX17FLAGS", intern = TRUE)))
cat("OpenMP/thread information:\n")
print(sblr:::sparseLD_thread_info(1L))
cat("Updates enabled in representative Gibbs iterations:\n")
cat("  BayesC: marker effects/states, effect variance, residual variance, pi, genetic variance\n")
cat("  BayesR: marker effects/components, effect variance, residual variance, pi, genetic variance\n")
print(native, row.names = FALSE)
cat("\nBayesC Gibbs speedup versus direct-quadratic reference:\n")
print(reference_speedup(
  "BayesC", "bayesc_gibbs_iteration",
  "bayesc_gibbs_iteration_direct_reference"
), row.names = FALSE)
cat("\nBayesR Gibbs speedup versus pre-patch arithmetic reference:\n")
print(reference_speedup(
  "BayesR", "bayesr_gibbs_iteration",
  "bayesr_gibbs_iteration_direct_reference"
), row.names = FALSE)

if (identical(tolower(Sys.getenv("SBLR_RUN_ILLUSTRATIVE_R_BENCHMARK")),
              "true")) {
  cat("\nIllustrative R matrix-loop benchmark (not native throughput)\n")
  set.seed(5119)
  elapsed <- function(fun, sweeps = 20L) {
    system.time(for (iteration in seq_len(sweeps)) fun())[["elapsed"]] / sweeps
  }
  illustrative <- lapply(c(250L, 500L, 1000L, 2000L), function(m) {
    k <- ceiling(0.25 * m)
    delta <- rnorm(m, sd = 0.002)
    Q <- matrix(rnorm(k * m, sd = 1 / sqrt(m)), k, m)
    data.frame(
      block_size = m,
      retained_rank = k,
      seconds_per_sweep = elapsed(function() {
        residual <- numeric(k)
        for (marker in seq_len(m)) {
          residual <- residual - Q[, marker] * delta[marker]
        }
        residual
      })
    )
  })
  print(do.call(rbind, illustrative), row.names = FALSE)
}
