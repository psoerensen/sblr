source("tests/testthat/helper-mtblr-bed-convergence-contract.R", local = TRUE)
grid <- expand.grid(nchains = c(2L, 4L, 8L), nit = c(100L, 1000L, 5000L),
                    nt = c(1L, 5L, 20L), nmodels = c(2L, 16L, 256L, 4096L),
                    selected_markers = c(0L, 10L, 100L, 1000L))
memory <- do.call(rbind, lapply(seq_len(nrow(grid)), function(i) {
  cbind(grid[i, ], as.data.frame(as.list(phase17t_memory(
    grid$nchains[i], grid$nit[i], grid$nt[i], grid$nmodels[i],
    grid$selected_markers[i]))))
}))
fixture <- phase17t_fixtures()$well_mixed
elapsed <- function(expr) unname(system.time(force(expr))[["elapsed"]])
timing <- data.frame(
  synthetic_rhat_seconds = elapsed(replicate(100, phase17t_rhat(fixture))),
  synthetic_ess_seconds = elapsed(replicate(100, phase17t_ess_bulk(fixture))),
  synthetic_mcse_seconds = elapsed(replicate(100, phase17t_mcse_mean(fixture))))
print(memory, row.names = FALSE)
print(timing, row.names = FALSE)
cat("Phase 17T contract-oracle signals only; no production-runtime claim.\n")
