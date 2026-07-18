measure_peak_rss <- function(command, args = character(), interval_seconds = 0.05,
                             env = character()) {
  stopifnot(interval_seconds > 0)
  if (!requireNamespace("processx", quietly = TRUE) ||
      !requireNamespace("ps", quietly = TRUE)) {
    stop("measure_peak_rss requires the development packages processx and ps")
  }
  started <- Sys.time()
  process <- processx::process$new(command, args, env = env,
                                   stdout = "|", stderr = "|")
  peak <- 0
  final <- NA_real_
  samples <- 0L
  sample_tree <- function(pid) {
    root_rss <- tryCatch({
      info <- process$get_memory_info()
      unname(info[[if ("rss" %in% names(info)) "rss" else "wset"]])
    }, error = function(e) 0)
    root <- tryCatch(ps::ps_handle(pid), error = function(e) NULL)
    if (is.null(root)) return(root_rss)
    handles <- tryCatch(ps::ps_children(root, recursive = TRUE),
                        error = function(e) list())
    rss <- vapply(handles, function(handle) {
      tryCatch({
        info <- ps::ps_memory_info(handle)
        field <- if ("rss" %in% names(info)) "rss" else "wset"
        unname(info[[field]])
      }, error = function(e) 0)
    }, numeric(1))
    root_rss + sum(rss)
  }
  repeat {
    rss <- sample_tree(process$get_pid())
    if (is.finite(rss)) {
      final <- rss
      peak <- max(peak, rss)
      samples <- samples + 1L
    }
    if (!process$is_alive()) break
    Sys.sleep(interval_seconds)
  }
  process$wait()
  list(
    exit_status = process$get_exit_status(),
    peak_rss_bytes = peak,
    final_sampled_rss_bytes = final,
    sampling_interval_seconds = interval_seconds,
    sample_count = samples,
    process_tree_method = "ps root plus recursive descendants",
    platform = paste(Sys.info()[c("sysname", "release")], collapse = " "),
    elapsed_seconds = as.numeric(difftime(Sys.time(), started, units = "secs")),
    stdout = paste(process$read_all_output_lines(), collapse = "\n"),
    stderr = paste(process$read_all_error_lines(), collapse = "\n")
  )
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (!identical(args, "--smoke")) {
    stop("usage: Rscript tools/benchmarks/measure_peak_rss.R --smoke")
  }
  child <- file.path(R.home("bin"), "Rscript")
  expression <- paste0(
    "x <- raw(16 * 1024^2); ",
    "cat(length(x)); flush.console(); Sys.sleep(1.5)"
  )
  result <- measure_peak_rss(child, c("-e", expression), 0.02)
  dput(result)
  if (result$exit_status != 0L || result$peak_rss_bytes <= 0) quit(status = 1L)
}
