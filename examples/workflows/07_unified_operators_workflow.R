# Phase 18 unified operator workflow
#
# This script assumes `stats`, `Glist`, `y`, `ld_prefix`, and suitable
# `block_start` objects have already been prepared. Small iteration counts are
# illustrative only and are not convergence recommendations.

common <- list(
  nit = 1000, nburn = 500, nthin = 1, seed = 1,
  nchains = 4L, ncores = 2L, keep_chains = FALSE,
  convergence = "auto",
  convergence_control = list(warn = TRUE, keep_traces = FALSE)
)

# Scalar-trait family: traits are independent logical tasks.
fit_st_csr <- do.call(stblr_csr, c(
  list(stats = stats, ld_prefix = ld_prefix, method = "sbayesc"), common))
fit_st_block <- do.call(stblr_block_eigen, c(
  list(stats = stats, Glist = Glist, block_start = block_start,
       method = "sbayesc"), common))
fit_st_bed <- do.call(stblr_bed, c(
  list(y = y, Glist = Glist, method = "bayesc"), common))

# Joint small-T family: one logical task is one complete joint chain.
fit_mt_csr <- do.call(mtblr_csr, c(
  list(stats = stats, ld_prefix = ld_prefix, method = "sbayesc"), common))
fit_mt_block <- do.call(mtblr_block_eigen, c(
  list(stats = stats, Glist = Glist, block_start = block_start,
       method = "sbayesc"), common))
fit_mt_bed <- do.call(mtblr_bed, c(
  list(y = y, Glist = Glist, method = "bayesc"), common))

# Every canonical fit shares these top-level contracts.
lapply(list(fit_st_csr, fit_st_block, fit_st_bed,
            fit_mt_csr, fit_mt_block, fit_mt_bed), function(fit) {
  fit[c("family", "model", "operator", "input", "data", "diagnostics",
        "convergence", "convergence_traces", "chains", "memory_estimate")]
})
