pkgload::load_all(".", compile = FALSE, quiet = TRUE)
source("tests/testthat/helper-mtblr-bed-contract.R", local = TRUE)
source("tests/testthat/helper-mtblr-bed-public.R", local = TRUE)

benchmark_case <- function(size) {
  if (size == "small") return(phase17p_case(nt = 2L))
  n <- 28L
  m <- 24L
  dosage <- outer(seq_len(n), seq_len(m), function(i, j) (i + 2 * j) %% 3)
  dosage[(row(dosage) + 3 * col(dosage)) %% 17 == 0] <- NA_real_
  path <- tempfile("phase17s-moderate-", fileext = ".bed")
  phase17n_write_bed(path, dosage)
  ids <- paste0("id", seq_len(n))
  Y <- vapply(seq_len(3L), function(t) {
    y <- sin(seq_len(n) * (.11 + t / 19)) +
      cos(seq_len(n) * (.07 + t / 23))
    y - mean(y)
  }, numeric(n))
  colnames(Y) <- paste0("T", 1:3)
  rownames(Y) <- ids
  Glist <- list(
    n = n, ids = ids, bedfiles = path,
    rsids = list(paste0("rs", seq_len(m))),
    rsidsLD = list(paste0("rs", seq_len(m))),
    af = list(seq(.18, .42, length.out = m)))
  list(fixture = list(paths = path), Glist = Glist, Y = Y, rows = NULL)
}

rows <- list()
index <- 0L
for (size in c("small", "moderate"))
  for (mode in c("diagonal", "full"))
    for (nchains in c(1L, 2L, 4L))
      for (ncores in c(1L, 2L, 4L))
        for (keep in c(FALSE, TRUE)) {
          case <- benchmark_case(size)
          args <- phase17p_public_args(
            case, mode, updates = TRUE, center = FALSE,
            nchains = nchains, ncores = ncores, keep_chains = keep,
            nit = 8L, nburn = 2L, memory_warning_gb = Inf)
          elapsed <- system.time(
            fit <- suppressWarnings(do.call(mtblr_bed, args)))["elapsed"]
          d <- fit$chain_diagnostics
          index <- index + 1L
          rows[[index]] <- data.frame(
            size = size, n = fit$input$n, m = fit$input$m, nt = fit$input$nt,
            mode = mode, nchains = nchains, requested_cores = ncores,
            used_workers = d$used_workers, keep_chains = keep,
            public_preparation_plus_format_seconds = max(
              0, unname(elapsed) - d$dispatch_seconds),
            native_dispatch_seconds = d$dispatch_seconds,
            total_public_seconds = unname(elapsed),
            chain_seconds_mean = d$seconds_mean,
            chain_seconds_max = d$seconds_max,
            requested_memory_gib = fit$memory_estimate$estimated_total_gib,
            execution_memory_gib =
              fit$memory_estimate$execution_estimated_total_gib,
            fit_bytes = as.numeric(object.size(fit)),
            retained_chain_bytes = as.numeric(object.size(fit$chains)))
          unlink(case$fixture$paths)
        }
print(do.call(rbind, rows), row.names = FALSE)
cat("Phase 17S regression signals only; no pure-adapter or linear-speedup claim.\n")
