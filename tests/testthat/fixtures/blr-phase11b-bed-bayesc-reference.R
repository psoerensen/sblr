source(testthat::test_path("fixtures", "blr-phase11a-bed-reference.R"))

phase11b_capture <- function(route = c("multichain", "single"), ncores = 1L,
    nchains = 2L, seed = 71L) {
  route <- match.arg(route)
  if (route == "multichain") {
    return(phase11a_capture("bayesc", ncores, nchains, seed))
  }
  x <- phase11a_fixture()
  ns <- asNamespace("sblr")
  assign(".phase11b_raw", NULL, envir = .GlobalEnv)
  suppressMessages(trace(".as_stblr_fit",
    tracer = quote(assign(".phase11b_raw", raw, envir = .GlobalEnv)),
    where = ns, print = FALSE))
  on.exit(suppressMessages(try(untrace(".as_stblr_fit", where = ns), silent = TRUE)))
  fit <- sblr::stblr_bed_marker(x$Glist, x$y, pi_init = .5,
    pi_prior_mean = .5, pi_prior_strength = 4, nit = 6L, nburn = 2L,
    nthin = 1L, seed = seed, ncores = ncores, nchains = 1L,
    backend = "scheduled", updateB = FALSE, updateE = FALSE,
    rebuild_every = 2L, read_block_size = 2L)
  list(raw = get(".phase11b_raw", envir = .GlobalEnv), fit = fit)
}

phase11b_metadata <- function(route, nchains, ncores) list(
  starting_commit = "e2035ab",
  rng_ownership = "bed_scheduled_bayesc_chain_rng_v1",
  reference_mode = "post-correction deterministic",
  route = route, samples = 6L, markers = 2L, seed = 71L,
  nchains = nchains, ncores = ncores, nit = 6L, nburn = 2L, nthin = 1L,
  scheduler = list(full_sweep_every = 10L, null_skip_base = 50L,
    null_skip_max = 200L, candidate_threshold = .001,
    candidate_lifetime = 20L, skip_nulls_burnin_only = FALSE),
  schema = "stblr_raw_v1", R = R.version.string)
