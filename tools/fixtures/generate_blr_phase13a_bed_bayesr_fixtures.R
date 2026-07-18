# Manual maintenance tool; Phase 13A tests never regenerate references.
pkgload::load_all(".", compile = FALSE, quiet = TRUE)
source(file.path("tests", "testthat", "fixtures", "blr-phase13a-bed-bayesr-reference.R"))
out <- file.path("tests", "testthat", "fixtures", "blr_phase13a_bed_bayesr")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
configs <- list(one_chain_one_core = c(1L, 1L, 71L),
  two_chains_one_core = c(1L, 2L, 73L),
  two_chains_two_cores = c(2L, 2L, 73L))
for (name in names(configs)) {
  z <- configs[[name]]
  value <- phase13a_capture(z[1], z[2], z[3])
  value$metadata <- list(source_commit = "8c41d18", R = R.version.string,
    route = "stblr_bed(method = bayesr)", reference_mode = "fresh deterministic",
    samples = 6L, markers = 2L, components = 4L,
    scales = c(0, .01, .1, 1), initial_pi = c(.95, .03, .015, .005),
    seed = z[3], cores = z[1], chains = z[2], nit = 6L, nburn = 2L,
    nthin = 1L, full_sweep_every = 10L, null_skip_base = 50L,
    null_skip_max = 200L, schema = "stblr_raw_v1",
    rng_ownership = "logical-chain-owned")
  saveRDS(value, file.path(out, paste0(name, ".rds")), version = 3)
}
