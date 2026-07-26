pkgload::load_all(".", compile = FALSE, quiet = TRUE)
source("tests/testthat/helper-mtblr-bed-contract.R", local = TRUE)
source("tests/testthat/helper-mtblr-bed-internal.R", local = TRUE)
source("tests/testthat/helper-mtblr-bed-multichain-contract.R", local = TRUE)
source("tests/testthat/helper-mtblr-bed-chains-internal.R", local = TRUE)
source("tests/testthat/helper-mtblr-bed-convergence-contract.R", local = TRUE)
source("tests/testthat/helper-mtblr-bed-convergence-engine.R", local = TRUE)

elapsed <- function(expr) unname(system.time(force(expr))[["elapsed"]])
grid <- expand.grid(
  nchains = c(2L, 4L, 8L), nit = c(100L, 1000L, 5000L),
  nt = c(1L, 5L, 20L))
synthetic <- do.call(rbind, lapply(seq_len(nrow(grid)), function(i) {
  nchains <- grid$nchains[i]; nit <- grid$nit[i]; nt <- grid$nt[i]
  draws <- vapply(seq_len(nchains), function(chain) {
    sin(seq_len(nit) * (.013 + chain / 1000)) +
      cos(seq_len(nit) * (.007 + chain / 1200))
  }, numeric(nit))
  scalar_seconds <- elapsed(
    sblr:::.mtblr_convergence_scalar(draws))
  values <- array(rep(draws, 3L * nt), c(nit, nchains, 3L * nt))
  quantities <- data.frame(
    quantity_index = seq_len(3L * nt),
    group = rep(c("B_diag", "G_diag", "E_diag"), each = nt),
    trait_index = rep(seq_len(nt), 3L), updated = TRUE)
  bundle <- list(
    schema = list(class = "mtblr_convergence_trace_bundle", version = 1L),
    scope = "core", nchains = nchains,
    postburn_draws_per_chain = nit,
    quantities = quantities, values = values)
  tier1_seconds <- elapsed(sblr:::.mtblr_convergence_tier1(
    bundle, paste0("T", seq_len(nt))))
  memory <- sblr:::.mtblr_convergence_memory_estimate(
    nchains, nit, nt, keep_traces = TRUE)
  data.frame(
    nchains = nchains, nit = nit, nt = nt,
    scalar_diagnostic_seconds = scalar_seconds,
    tier1_diagnostic_seconds = tier1_seconds,
    trace_capture_bytes = memory$trace_capture_bytes,
    workspace_bytes = memory$maximum_workspace_bytes,
    retained_trace_bytes = memory$retained_trace_bytes,
    estimated_total_gib = memory$estimated_total_gib)
}))
print(synthetic, row.names = FALSE)

case <- phase17o_case(nt = 2L, updates = TRUE)
case$nit <- 20L
case$nburn <- 5L
on.exit(phase17o_cleanup(case), add = TRUE)
route_seconds <- elapsed(native <- phase17u_native_call(
  case, nchains = 2L, ncores = 1L))
diagnostic_seconds <- elapsed(diagnosed <- phase17u_diagnose(
  native, colnames(case$Y), TRUE, TRUE, keep_traces = TRUE))
cat("ACTUAL_TRACE_EXTRACTION_PLUS_ROUTE_SECONDS=", route_seconds, "\n", sep = "")
cat("ACTUAL_TIER1_ENGINE_SECONDS=", diagnostic_seconds, "\n", sep = "")
cat("ACTUAL_TRACE_BYTES=", as.numeric(object.size(native$trace_bundle)), "\n",
    sep = "")
cat("ACTUAL_RETAINED_TRACE_BYTES=",
    as.numeric(object.size(diagnosed$convergence_traces)), "\n", sep = "")
cat("Phase 17U internal regression signals only; no public-runtime claim.\n")
