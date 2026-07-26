source(file.path("tests", "testthat",
                 "helper-mtblr-bed-extended-convergence-contract.R"))
pkgload::load_all(".", compile = FALSE, quiet = TRUE)

grid <- expand.grid(
  nchains = c(2L, 4L, 8L), nit = c(100L, 1000L, 5000L),
  nt = c(2L, 5L, 10L, 20L, 50L),
  nmodels = c(2L, 16L, 256L, 4096L, 16384L),
  selected_patterns = c(0L, 10L, 100L, 1000L),
  selected_markers = c(0L, 10L, 100L, 1000L, 10000L),
  marker_quantities = c("b", "d", "both"),
  keep_traces = c(FALSE, TRUE), KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE)

bytes <- lapply(seq_len(nrow(grid)), function(i) {
  x <- grid[i, ]
  C <- x$nchains; N <- x$nit; T <- x$nt
  qoff <- T * (T - 1) / 2
  psel <- min(x$selected_patterns, x$nmodels)
  b <- 8 * C * N * qoff
  g <- 8 * C * N * qoff
  e <- 8 * C * N * qoff
  mass <- 8 * C * N
  selected_pi <- 8 * C * N * psel
  all_pi <- 8 * C * N * x$nmodels
  marker_b <- if (x$marker_quantities %in% c("b", "both"))
    8 * C * N * x$selected_markers * T else 0
  marker_d <- if (x$marker_quantities %in% c("d", "both"))
    4 * C * N * x$selected_markers * T else 0
  captured <- b + g + e + mass + selected_pi + marker_b + marker_d
  data.frame(
    tier1_bytes = 8 * C * N * 3 * T,
    B_cov_bytes = b, G_cov_bytes = g, E_cov_bytes = e,
    probability_mass_bytes = mass, selected_pattern_bytes = selected_pi,
    all_pattern_bytes = all_pi, selected_b_bytes = marker_b,
    selected_d_bytes = marker_d, captured_total_bytes = captured,
    retained_total_bytes = if (x$keep_traces) captured else 0,
    workspace_bytes = 8 * C * N * 8,
    summary_rows = 3 * T + 3 * qoff + 1 + psel +
      as.integer(x$marker_quantities %in% c("b", "both")) *
        x$selected_markers * T +
      as.integer(x$marker_quantities %in% c("d", "both")) *
        x$selected_markers * T)
})
grid <- cbind(grid, do.call(rbind, bytes))

timings <- do.call(rbind, lapply(c(2L, 4L, 8L), function(C) {
  do.call(rbind, lapply(c(100L, 1000L, 5000L), function(N) {
    base <- outer(seq_len(N), seq_len(C), function(i, j)
      sin(i / 11 + j / 7) + cos(i / 29 - j / 5))
    traces <- list(
      covariance = base,
      probability = plogis(base),
      zero_inflated_b = base * (base > .3),
      binary_d = (base > .3) * 1)
    do.call(rbind, lapply(names(traces), function(kind) {
      elapsed <- system.time(
        sblr:::.mtblr_convergence_scalar(traces[[kind]]))["elapsed"]
      data.frame(nchains = C, nit = N, kind = kind,
                 scalar_seconds = unname(elapsed))
    }))
  }))
}))

cat("PHASE17W_SYNTHETIC_FEASIBILITY_ONLY=TRUE\n")
cat("GRID_ROWS=", nrow(grid), "\n", sep = "")
cat("DEFAULT_CAPTURE_LIMIT_GIB=2\n")
cat("MAX_CAPTURE_GIB=",
    format(max(grid$captured_total_bytes) / 1024^3, digits = 6), "\n", sep = "")
cat("MAX_ALL_PATTERN_GIB=",
    format(max(grid$all_pattern_bytes) / 1024^3, digits = 6), "\n", sep = "")
cat("MAX_RETAINED_GIB=",
    format(max(grid$retained_total_bytes) / 1024^3, digits = 6), "\n", sep = "")
cat("CAPTURE_CASES_ABOVE_LIMIT=",
    sum(grid$captured_total_bytes > 2 * 1024^3), "\n", sep = "")
closest <- order(abs(grid$captured_total_bytes - 2 * 1024^3))[seq_len(6L)]
print(grid[closest, c("nchains", "nit", "nt", "nmodels",
  "selected_patterns", "selected_markers", "marker_quantities",
  "keep_traces", "captured_total_bytes", "retained_total_bytes",
  "workspace_bytes", "summary_rows")], row.names = FALSE)
print(aggregate(scalar_seconds ~ kind + nchains + nit, timings, max),
      row.names = FALSE)
cat("No production-runtime or linear-scaling claim is made.\n")
