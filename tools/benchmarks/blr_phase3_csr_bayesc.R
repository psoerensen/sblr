#!/usr/bin/env Rscript

# Phase 3 deliberately reuses the Phase 2 fixture, warm-up, timing, and sampled
# whole-process RSS implementation so the measurements remain comparable. The
# Phase 2 driver is read only through its function definitions; its main block
# is not executed here.
phase2_script <- file.path("tools", "benchmarks", "blr_phase2_csr_bayesc.R")
phase2_lines <- readLines(phase2_script, warn = FALSE)
phase2_main <- grep("^options <- phase2_benchmark_options\\(\\)$", phase2_lines)
if (length(phase2_main) != 1L) {
  stop("Could not locate the Phase 2 benchmark entry point.")
}
eval(parse(text = phase2_lines[seq_len(phase2_main - 1L)]), envir = environment())

options <- phase2_benchmark_options(commandArgs(trailingOnly = TRUE))
if (!"sblr" %in% loadedNamespaces()) pkgload::load_all(".", compile = FALSE)
results <- phase2_run_benchmarks(options)

if (isTRUE(options$child)) {
  if (is.null(options$result)) stop("--result is required for child mode")
  saveRDS(results, options$result)
} else {
  cat("Phase 3 canonical CSR BayesC benchmark\n")
  print(results, row.names = FALSE)
  cat("\nWarm-up excluded; timing summaries:\n")
  print(aggregate(
    elapsed_seconds ~ fixture + markers + traits + nit + nburn + chains +
      cores + output,
    results,
    function(x) c(mean = mean(x), median = stats::median(x), sd = stats::sd(x),
                  min = min(x), max = max(x))
  ), row.names = FALSE)
  cat("\nEnvironment:\n")
  print(sessionInfo())
  cat("\nMeasurement note: elapsed time covers the complete public call; ",
      "sampled peak RSS uses the unchanged Phase 2 child-process method. ",
      "Initialization and result conversion are not separately instrumented.\n",
      sep = "")
}
