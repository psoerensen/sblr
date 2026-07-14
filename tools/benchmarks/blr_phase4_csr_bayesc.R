#!/usr/bin/env Rscript

# Phase 4 reuses the Phase 3 fixture and measurement implementation verbatim
# so task/seed infrastructure extraction is compared on the same workloads.
phase3_script <- file.path("tools", "benchmarks", "blr_phase3_csr_bayesc.R")
phase3_lines <- readLines(phase3_script, warn = FALSE)
phase3_main <- grep("^options <- phase2_benchmark_options", phase3_lines)
if (length(phase3_main) != 1L) {
  stop("Could not locate the Phase 3 benchmark entry point.")
}
eval(parse(text = phase3_lines[seq_len(phase3_main - 1L)]),
     envir = environment())

options <- phase2_benchmark_options(commandArgs(trailingOnly = TRUE))
if (!"sblr" %in% loadedNamespaces()) pkgload::load_all(".", compile = FALSE)
results <- phase2_run_benchmarks(options)

if (isTRUE(options$child)) {
  if (is.null(options$result)) stop("--result is required for child mode")
  saveRDS(results, options$result)
} else {
  cat("Phase 4 shared-infrastructure CSR BayesC benchmark\n")
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
