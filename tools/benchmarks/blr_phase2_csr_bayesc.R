# Reproducible post-migration benchmark for unscheduled CSR BayesC.
# Run from the package root, for example:
# Rscript tools/benchmarks/blr_phase2_csr_bayesc.R --fixture=moderate --nrep=5

phase1_script <- file.path("tools", "benchmarks", "blr_phase1_csr_bayesc.R")
phase1_lines <- readLines(phase1_script, warn = FALSE)
phase1_entry <- grep("^options <- phase1_benchmark_options", phase1_lines)[1L]
eval(parse(text = phase1_lines[seq_len(phase1_entry - 1L)]), envir = environment())

phase2_benchmark_options <- function(args = commandArgs(trailingOnly = TRUE)) {
  values <- list(
    fixture = "all", nrep = 3L, chains = "all", cores = "all",
    output = "all", peak = TRUE, child = FALSE, result = NULL
  )
  for (arg in args) {
    pieces <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1L]]
    key <- pieces[[1L]]
    value <- if (length(pieces) == 1L) "TRUE" else pieces[[2L]]
    if (!key %in% names(values)) stop("Unknown option: --", key)
    values[[key]] <- value
  }
  values$nrep <- as.integer(values$nrep)
  values$peak <- tolower(as.character(values$peak)) %in% c("true", "1", "yes")
  values$child <- tolower(as.character(values$child)) %in% c("true", "1", "yes")
  values
}

phase2_peak_rss <- function(script, fixture, chains, cores, output) {
  if (!requireNamespace("processx", quietly = TRUE) ||
      !requireNamespace("ps", quietly = TRUE)) {
    return(list(bytes = NA_real_, method =
      "unavailable: optional processx and ps packages are required"))
  }
  result_file <- tempfile("blr_phase2_peak_", fileext = ".rds")
  args <- c(
    script, "--child=TRUE", paste0("--fixture=", fixture), "--nrep=1",
    paste0("--chains=", chains), paste0("--cores=", cores),
    paste0("--output=", output), "--peak=FALSE",
    paste0("--result=", result_file)
  )
  process <- processx::process$new(
    file.path(R.home("bin"), "Rscript.exe"), args,
    stdout = "|", stderr = "|", cleanup_tree = TRUE
  )
  peak <- 0
  child_output <- character()
  deadline <- Sys.time() + 120
  while (process$is_alive()) {
    info <- tryCatch(process$get_memory_info(), error = function(e) NULL)
    if (!is.null(info) && is.finite(info[["rss"]])) peak <- max(peak, info[["rss"]])
    poll <- process$poll_io(10)
    if (identical(unname(poll[["output"]]), "ready")) {
      child_output <- c(child_output, process$read_output_lines())
    }
    if (identical(unname(poll[["error"]]), "ready")) {
      child_output <- c(child_output, process$read_error_lines())
    }
    if (Sys.time() > deadline) {
      process$kill_tree()
      stop("Peak-memory child exceeded 120 seconds. Output:\n",
           paste(child_output, collapse = "\n"))
    }
  }
  process$wait()
  if (!identical(process$get_exit_status(), 0L)) {
    stop("Peak-memory child failed:\n", paste(c(
      child_output, process$read_all_output_lines(),
      process$read_all_error_lines()
    ), collapse = "\n"))
  }
  list(
    bytes = as.numeric(peak),
    method = "maximum child RSS sampled every 10 ms via processx/ps"
  )
}

phase2_run_benchmarks <- function(options) {
  fixtures <- if (identical(options$fixture, "all")) c("tiny", "moderate") else {
    match.arg(options$fixture, c("tiny", "moderate"))
  }
  chains <- if (identical(options$chains, "all")) c(1L, 2L) else as.integer(options$chains)
  cores <- if (identical(options$cores, "all")) c(1L, 2L) else as.integer(options$cores)
  outputs <- if (identical(options$output, "all")) c("minimal", "ordinary") else {
    match.arg(options$output, c("minimal", "ordinary"))
  }
  configurations <- expand.grid(
    fixture = fixtures, chains = chains, cores = cores, output = outputs,
    stringsAsFactors = FALSE
  )
  configurations <- configurations[
    configurations$cores <= configurations$chains, , drop = FALSE
  ]
  script <- normalizePath(
    file.path("tools", "benchmarks", "blr_phase2_csr_bayesc.R"),
    winslash = "/", mustWork = TRUE
  )
  rows <- vector("list", nrow(configurations) * options$nrep)
  row_index <- 0L
  for (index in seq_len(nrow(configurations))) {
    config <- configurations[index, ]
    fixture <- phase1_benchmark_fixture(config$fixture)
    # One unrecorded run warms package dispatch, CSR loading, and native code.
    invisible(phase1_benchmark_once(
      fixture, as.integer(config$chains), as.integer(config$cores), config$output
    ))
    peak <- if (isTRUE(options$peak)) {
      phase2_peak_rss(script, config$fixture, config$chains,
                      config$cores, config$output)
    } else list(bytes = NA_real_, method = "not requested")
    for (repetition in seq_len(options$nrep)) {
      gc()
      run <- phase1_benchmark_once(
        fixture, as.integer(config$chains), as.integer(config$cores), config$output
      )
      row_index <- row_index + 1L
      rows[[row_index]] <- data.frame(
        fixture = config$fixture, markers = fixture$marker_count,
        traits = fixture$trait_count, nit = fixture$nit, nburn = fixture$nburn,
        nthin = fixture$nthin, retained_samples = run$retained_samples,
        trace_rows = run$trace_rows, chains = config$chains, cores = config$cores,
        output = config$output, repetition = repetition,
        elapsed_seconds = run$elapsed, output_bytes = run$output_bytes,
        peak_rss_bytes = peak$bytes, peak_memory_method = peak$method,
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}

options <- phase2_benchmark_options()
if (!"sblr" %in% loadedNamespaces()) pkgload::load_all(".", compile = FALSE)
results <- phase2_run_benchmarks(options)
if (isTRUE(options$child)) {
  if (is.null(options$result)) stop("--result is required for child mode")
  saveRDS(results, options$result)
} else {
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
  cat("\nConversion note: elapsed time covers the complete public call, including ",
      "validation, CSR loading, typed conversion, execution, and result conversion; ",
      "the conversion boundary is not separately instrumented.\n", sep = "")
}
