root <- normalizePath(".", mustWork = TRUE)
suppressPackageStartupMessages(pkgload::load_all(root, compile = FALSE,
                                                 quiet = TRUE))
environment <- new.env(parent = globalenv())
sys.source(file.path(root, "tests/testthat/helper-mtblr-bed-contract.R"),
           environment)
sys.source(file.path(root, "tests/testthat/helper-mtblr-bed-internal.R"),
           environment)

run <- function(label, nt, mode, iterations) {
  case <- environment$phase17o_case(
    nt = nt, residual_covariance = mode, updates = TRUE,
    multiple_sets = nt > 2L)
  on.exit(environment$phase17o_cleanup(case), add = TRUE)
  case$nit <- as.integer(iterations)
  case$nburn <- 2L
  elapsed <- system.time(raw <- environment$phase17o_call(case))[["elapsed"]]
  n <- nrow(case$fixture$X)
  m <- ncol(case$fixture$X)
  stride <- 64 * ceiling(ceiling(n / 4) / 64)
  data.frame(
    label = label, n = n, m = m, nt = nt,
    models = length(case$models), sets = length(case$sets),
    residual_mode = mode, iterations = case$nit, burnin = case$nburn,
    thinning = case$nthin,
    bed_read_seconds = NA_real_,
    map_wy_order_preparation_seconds = NA_real_,
    mcmc_seconds = NA_real_,
    final_marker_r_seconds = NA_real_,
    total_seconds = elapsed,
    timing_separation = "not exposed by internal binding",
    packed_bytes = m * stride,
    sample_residual_bytes = 8 * n * nt,
    workspace_bytes = 8 * n,
    fit_object_bytes = as.numeric(object.size(raw)),
    completed_fit_rss = NA_real_
  )
}

result <- rbind(
  run("small_diagonal", 2L, "diagonal", 5L),
  run("small_full", 2L, "full", 5L),
  run("moderate_diagonal", 3L, "diagonal", 20L),
  run("moderate_full", 3L, "full", 20L)
)
print(result, row.names = FALSE)
cat("Phase 17O timings are regression signals, not performance claims.\n")
