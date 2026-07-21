source(testthat::test_path("fixtures", "blr-phase11a-bed-reference.R"))

phase13a_capture <- function(ncores = 1L, nchains = 1L, seed = 71L,
                             full_sweep_every = 10L,
                             null_skip_base = 50L,
                             null_skip_max = 200L,
                             updatePi = FALSE,
                             progress_every = 0) {
  x <- phase11a_fixture()
  ns <- asNamespace("sblr")
  assign(".phase13a_raw", NULL, envir = .GlobalEnv)
  suppressMessages(trace(".as_stblr_fit",
    tracer = quote(assign(".phase13a_raw", raw, envir = .GlobalEnv)),
    where = ns, print = FALSE))
  on.exit(suppressMessages(try(untrace(".as_stblr_fit", where = ns), silent = TRUE)))
  fit <- sblr::stblr_bed(y = x$y, Glist = x$Glist, method = "bayesr",
    nit = 6L, nburn = 2L, nthin = 1L, seed = seed, ncores = ncores,
    nchains = nchains, updateB = FALSE, updateE = FALSE,
    updatePi = updatePi, mixture_var = c(0, .01, .1, 1),
    pi = c(.95, .03, .015, .005), rebuild_every = 2L,
    read_block_size = 2L, full_sweep_every = full_sweep_every,
    null_skip_base = null_skip_base, null_skip_max = null_skip_max,
    progress_every = progress_every)
  list(raw = get(".phase13a_raw", envir = .GlobalEnv), fit = fit)
}

phase13a_normalize <- phase11a_normalize
