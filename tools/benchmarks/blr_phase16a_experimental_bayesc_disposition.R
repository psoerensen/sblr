source("tests/testthat/fixtures/blr-phase11a-bed-reference.R")

cat("Phase 16A experimental packed-BED BayesC disposition evidence\n")
cat("These workloads are experimental capability checks, not canonical benchmarks.\n")
cat("completed-fit RSS is not peak RSS; page-cache effects apply.\n")

x <- phase11a_fixture()
run_route <- function(backend) {
  gc()
  before <- sum(gc()[, 2])
  elapsed <- system.time({
    fit <- sblr::stblr_bed_marker(x$Glist, x$y, backend = backend,
      pi_init = .5, pi_prior_mean = .5, pi_prior_strength = 4,
      nit = 20L, nburn = 5L, nthin = 1L, seed = 71L, ncores = 1L,
      updateB = FALSE, updateE = FALSE, rebuild_every = 2L,
      full_sweep_every = 2L, null_update_prob = .5)
  })[["elapsed"]]
  list(seconds = unname(elapsed), completed_fit_rss_mib = sum(gc()[, 2]),
       rss_delta_mib = sum(gc()[, 2]) - before, backend = fit$input$backend,
       samples = x$Glist$n, markers = sum(lengths(x$Glist$rsids)),
       bed_bytes = unname(file.info(x$Glist$bedfiles)$size))
}

for (backend in c("scheduled", "sparse")) {
  invisible(run_route(backend))
  values <- replicate(5L, run_route(backend), simplify = FALSE)
  times <- vapply(values, `[[`, numeric(1), "seconds")
  cat("\n", backend, " route\n", sep = "")
  print(list(disposition = "retain as explicitly experimental",
    reason = if (backend == "scheduled") "deterministic Phase 11B oracle"
      else "distinct null_update_prob scheduler research policy",
    repetitions = times, mean = mean(times), median = median(times),
    minimum = min(times), maximum = max(times), range = diff(range(times)),
    completed_fit_rss_mib = values[[length(values)]]$completed_fit_rss_mib,
    samples = values[[1]]$samples, markers = values[[1]]$markers,
    bed_bytes = values[[1]]$bed_bytes,
    R = R.version.string, package = as.character(utils::packageVersion("sblr"))))
}
