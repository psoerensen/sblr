source("tests/testthat/helper-mtblr-bed-multichain-contract.R", local = TRUE)
cases <- expand.grid(size = c("small", "moderate"), nchains = c(1L, 2L, 4L, 8L),
                     ncores = c(1L, 2L, 4L), keep_chains = c(FALSE, TRUE),
                     stringsAsFactors = FALSE)
rows <- lapply(seq_len(nrow(cases)), function(i) {
  z <- cases[i, ]
  dims <- if (z$size == "small") c(n = 100, m = 500, nt = 2, models = 4) else
    c(n = 1000, m = 5000, nt = 4, models = 16)
  mem <- phase17q_memory(dims["n"], dims["m"], dims["nt"], dims["models"],
                         40, z$nchains, z$ncores, z$keep_chains)
  chains <- lapply(seq_len(z$nchains), function(ch) phase17q_chain(
    ch, ch + 2L, matrix(ch, 20, dims["nt"]), matrix(ch / 10, 20, dims["nt"]),
    matrix(as.integer(ch %% 2), 20, dims["nt"])))
  elapsed <- system.time(phase17q_aggregate(chains, z$keep_chains))["elapsed"]
  data.frame(z, n = dims["n"], m = dims["m"], nt = dims["nt"],
             model_count = dims["models"], worker_count = mem$worker_count,
             shared_bytes = mem$shared_bytes,
             private_bytes_per_live_worker = mem$private_state_bytes_per_chain,
             result_bytes_per_chain = mem$result_bytes_per_chain,
             retained_chain_bytes = mem$estimated_retained_output_bytes,
             estimated_total_gib = mem$estimated_total_gib,
             synthetic_aggregation_seconds = unname(elapsed))
})
out <- do.call(rbind, rows)
print(out, row.names = FALSE)
cat("Audit-only synthetic calculations; no OpenMP speedup claim.\n")
