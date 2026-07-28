root <- normalizePath(if (file.exists("DESCRIPTION")) "." else "../..",
                      winslash = "/", mustWork = TRUE)
`%||%` <- function(x, y) if (is.null(x)) y else x
source(file.path(root, "R", "blr-extended-convergence.R"), local = TRUE)
source(file.path(root, "R", "mtblr-convergence.R"), local = TRUE)

grid <- expand.grid(
  family = c("stblr", "mtblr"),
  operator = c("csr", "block_eigen", "packed_bed"),
  traits = c(1L, 2L, 3L, 5L, 10L), chains = c(1L, 2L, 4L),
  draws = c(100L, 1000L), selected_markers = c(0L, 1L, 10L, 100L),
  stringsAsFactors = FALSE)
grid$covariance_quantities <- with(grid,
  ifelse(family == "mtblr", 3L * traits * (traits - 1L) / 2L, 0L))
grid$probability_quantities <- 8L
grid$annotation_quantities <- 20L
grid$selected_b <- with(grid, selected_markers * traits)
grid$selected_d <- grid$selected_b
grid$selected_component <- grid$selected_markers

started <- proc.time()[["elapsed"]]
memory <- lapply(seq_len(nrow(grid)), function(i) {
  x <- grid[i, ]
  .blr_extended_trace_memory(
    x$chains, x$draws,
    x$covariance_quantities + x$probability_quantities +
      x$annotation_quantities + x$selected_b,
    x$selected_d + x$selected_component, keep_traces = TRUE)
})
resolution_elapsed <- proc.time()[["elapsed"]] - started
grid$capture_bytes <- vapply(memory, `[[`, numeric(1), "trace_capture_bytes")
grid$workspace_bytes <- vapply(memory, `[[`, numeric(1), "workspace_bytes")
grid$retained_bytes <- vapply(memory, `[[`, numeric(1), "retained_trace_bytes")
grid$summary_bytes <- vapply(memory, `[[`, numeric(1), "summary_output_bytes")

set.seed(21)
draws <- matrix(stats::rnorm(1000L * 4L), 1000L, 4L)
diagnostic_started <- proc.time()[["elapsed"]]
invisible(.blr_convergence_scalar(draws))
diagnostic_elapsed <- proc.time()[["elapsed"]] - diagnostic_started

cat("PHASE21_BENCHMARK_ROWS=", nrow(grid), "\n", sep = "")
cat(sprintf("TRACE_PLAN_RESOLUTION_SECONDS=%.6f\n", resolution_elapsed))
cat(sprintf("SCALAR_DIAGNOSTIC_SECONDS=%.6f\n", diagnostic_elapsed))
cat("MAX_CAPTURE_BYTES=", max(grid$capture_bytes), "\n", sep = "")
cat("MAX_WORKSPACE_BYTES=", max(grid$workspace_bytes), "\n", sep = "")
cat("MAX_RETAINED_BYTES=", max(grid$retained_bytes), "\n", sep = "")
cat("MAX_SUMMARY_BYTES=", max(grid$summary_bytes), "\n", sep = "")
cat("LINEAR_SPEEDUP_CLAIM=FALSE\n")
