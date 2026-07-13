# Reproducible Phase 1 baseline for the existing unscheduled CSR BayesC route.
# Run from the package root, for example:
# Rscript tools/benchmarks/blr_phase1_csr_bayesc.R --fixture=moderate --nrep=5

phase1_benchmark_options <- function(args = commandArgs(trailingOnly = TRUE)) {
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

phase1_write_csr <- function(marker_count) {
  prefix <- tempfile("blr_phase1_benchmark_csr_")
  if (marker_count == 1L) {
    row_ptr <- c(0, 0)
    col_idx <- integer()
    values <- numeric()
  } else {
    row_ptr <- c(0:(marker_count - 1L), marker_count - 1L)
    col_idx <- seq_len(marker_count - 1L)
    values <- rep(0.1, marker_count - 1L)
  }
  sblr:::.stblr_write_uint64_file(
    paste0(prefix, ".row_ptr.u64.bin"), row_ptr
  )
  sblr:::.stblr_write_uint32_file(
    paste0(prefix, ".col_idx.u32.0based.bin"), col_idx
  )
  writeBin(values, paste0(prefix, ".values.f32.bin"),
           size = 4, endian = "little")
  writeLines(c(
    "format=sparse_ld_csr", "storage=streamed_upper_triangle",
    "n_bed=NA", "n_used=NA", "n_samples_used=NA",
    paste0("n_variants=", marker_count),
    paste0("nnz=", length(values)),
    "triangle=upper", "diagonal=implicit_1",
    paste0("row_ptr_file=", prefix, ".row_ptr.u64.bin"),
    paste0("col_idx_file=", prefix, ".col_idx.u32.0based.bin"),
    paste0("values_file=", prefix, ".values.f32.bin"),
    "row_ptr_type=uint64", "col_idx_type=uint32", "values_type=float32",
    "index_base=0", "value=r"
  ), paste0(prefix, ".meta.txt"))
  prefix
}

phase1_benchmark_fixture <- function(kind = c("tiny", "moderate")) {
  kind <- match.arg(kind)
  if (kind == "tiny") {
    marker_count <- 3L
    trait_count <- 1L
    nit <- 8L
    nburn <- 2L
  } else {
    marker_count <- 2000L
    trait_count <- 2L
    nit <- 80L
    nburn <- 20L
  }
  marker_ids <- paste0("m", seq_len(marker_count))
  trait_ids <- paste0("trait", seq_len(trait_count))
  sample_size <- rep(10000L, trait_count)
  wy <- lapply(seq_len(trait_count), function(trait) {
    stats::setNames(
      4 * sin(seq_len(marker_count) / (13 + trait)) +
        2 * cos(seq_len(marker_count) / (7 + trait)),
      marker_ids
    )
  })
  ww <- lapply(seq_len(trait_count), function(trait) {
    stats::setNames(rep(sample_size[[trait]], marker_count), marker_ids)
  })
  names(wy) <- names(ww) <- trait_ids
  list(
    kind = kind,
    stats = list(
      wy = wy, ww = ww,
      yy = stats::setNames(as.numeric(sample_size), trait_ids),
      n = sample_size[[1L]], m = marker_count,
      marker_names = marker_ids, trait_names = trait_ids
    ),
    prefix = phase1_write_csr(marker_count),
    marker_count = marker_count,
    trait_count = trait_count,
    nit = nit,
    nburn = nburn,
    nthin = 1L
  )
}

phase1_benchmark_once <- function(fixture, chains, cores, output) {
  keep_chains <- identical(output, "ordinary") && chains > 1L
  started <- proc.time()[["elapsed"]]
  fit <- stblr_csr(
    stats = fixture$stats,
    ld_prefix = fixture$prefix,
    pi_init = 0.01,
    pi_prior_mean = 0.01,
    pi_prior_strength = 100,
    updateB = TRUE,
    updateE = TRUE,
    updatePi = TRUE,
    nit = fixture$nit,
    nburn = fixture$nburn,
    nthin = fixture$nthin,
    seed = 1701L,
    nchains = chains,
    keep_chains = keep_chains,
    ncores = cores,
    updateLDswap = FALSE,
    scheduled = FALSE
  )
  elapsed <- proc.time()[["elapsed"]] - started
  list(
    elapsed = unname(elapsed),
    trace_rows = nrow(fit$vbs),
    retained_samples = length(seq.int(
      fixture$nburn + 1L,
      fixture$nburn + fixture$nit,
      by = fixture$nthin
    )),
    output_bytes = as.numeric(object.size(fit)),
    fit = fit
  )
}

phase1_peak_rss <- function(script, fixture, chains, cores, output) {
  if (!requireNamespace("processx", quietly = TRUE) ||
      !requireNamespace("ps", quietly = TRUE)) {
    return(list(bytes = NA_real_, method =
      "unavailable: optional processx and ps packages are required"))
  }
  result_file <- tempfile("blr_phase1_peak_", fileext = ".rds")
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
  while (process$is_alive()) {
    info <- tryCatch(process$get_memory_info(), error = function(e) NULL)
    if (!is.null(info) && is.finite(info[["rss"]])) {
      peak <- max(peak, info[["rss"]])
    }
    Sys.sleep(0.01)
  }
  process$wait()
  if (!identical(process$get_exit_status(), 0L)) {
    stop(
      "Peak-memory child failed:\n",
      paste(c(process$read_all_output_lines(),
              process$read_all_error_lines()), collapse = "\n")
    )
  }
  list(
    bytes = as.numeric(peak),
    method = "maximum child RSS sampled every 10 ms via processx/ps"
  )
}

phase1_run_benchmarks <- function(options) {
  fixtures <- if (identical(options$fixture, "all")) {
    c("tiny", "moderate")
  } else {
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
  configurations <- configurations[configurations$cores <= configurations$chains, , drop = FALSE]
  rows <- vector("list", nrow(configurations) * options$nrep)
  row_index <- 0L
  script <- normalizePath(
    file.path("tools", "benchmarks", "blr_phase1_csr_bayesc.R"),
    winslash = "/", mustWork = TRUE
  )
  for (index in seq_len(nrow(configurations))) {
    config <- configurations[index, ]
    fixture <- phase1_benchmark_fixture(config$fixture)
    peak <- if (isTRUE(options$peak)) {
      phase1_peak_rss(script, config$fixture, config$chains,
                      config$cores, config$output)
    } else {
      list(bytes = NA_real_, method = "not requested")
    }
    for (repetition in seq_len(options$nrep)) {
      gc()
      run <- phase1_benchmark_once(
        fixture, as.integer(config$chains), as.integer(config$cores),
        config$output
      )
      row_index <- row_index + 1L
      rows[[row_index]] <- data.frame(
        fixture = config$fixture,
        markers = fixture$marker_count,
        traits = fixture$trait_count,
        nit = fixture$nit,
        nburn = fixture$nburn,
        nthin = fixture$nthin,
        retained_samples = run$retained_samples,
        trace_rows = run$trace_rows,
        chains = config$chains,
        cores = config$cores,
        output = config$output,
        repetition = repetition,
        elapsed_seconds = run$elapsed,
        output_bytes = run$output_bytes,
        peak_rss_bytes = peak$bytes,
        peak_memory_method = peak$method,
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}

options <- phase1_benchmark_options()
if (!"sblr" %in% loadedNamespaces()) {
  pkgload::load_all(".", compile = FALSE)
}
results <- phase1_run_benchmarks(options)
if (isTRUE(options$child)) {
  if (is.null(options$result)) stop("--result is required for child mode")
  saveRDS(results, options$result)
} else {
  print(results, row.names = FALSE)
  cat("\nTiming summaries:\n")
  print(aggregate(
    elapsed_seconds ~ fixture + markers + traits + nit + nburn + chains +
      cores + output,
    results,
    function(x) c(mean = mean(x), sd = stats::sd(x), min = min(x), max = max(x))
  ), row.names = FALSE)
  cat("\nEnvironment:\n")
  print(sessionInfo())
  cat("\nOutput note: 'minimal' is the least-retained current public mode ",
      "(keep_chains = FALSE); 'ordinary' retains compact chains when ",
      "chains > 1. The production backend has no selective marker/trace ",
      "output switch.\n", sep = "")
}
