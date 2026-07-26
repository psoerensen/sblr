root <- normalizePath(".", mustWork = TRUE)
suppressPackageStartupMessages(pkgload::load_all(root, compile = FALSE,
                                                 quiet = TRUE))
environment <- new.env(parent = globalenv())
for (helper in c("helper-mtblr-bed-contract.R",
                 "helper-mtblr-bed-internal.R",
                 "helper-mtblr-bed-public.R")) {
  sys.source(file.path(root, "tests/testthat", helper), environment)
}

run <- function(label, mode, nt, iterations) {
  case <- environment$phase17p_case(nt = nt)
  on.exit(environment$phase17p_cleanup(case), add = TRUE)
  args <- environment$phase17p_public_args(case, mode, TRUE)
  args$nit <- as.integer(iterations)
  args$nburn <- 2L
  public_seconds <- system.time(
    fit <- suppressWarnings(do.call(mtblr_bed, args)))[["elapsed"]]
  internal_args <- environment$phase17p_native_args(args)
  internal_seconds <- system.time(
    raw <- do.call(getFromNamespace("mtblr_bed_internal", "sblr"),
                   internal_args))[["elapsed"]]
  data.frame(
    label = label, n = fit$input$n, m = fit$input$m, nt = nt,
    models = nrow(fit$input$models), sets = length(fit$input$sets),
    residual_mode = mode, iterations = iterations,
    public_preparation_seconds = NA_real_,
    native_call_seconds_separable = NA_real_,
    total_public_seconds = public_seconds,
    total_internal_seconds = internal_seconds,
    analytical_memory_gib = fit$memory_estimate$estimated_total_gib,
    fit_object_bytes = as.numeric(object.size(fit)),
    completed_fit_rss = NA_real_, raw_object_bytes = as.numeric(object.size(raw))
  )
}

print(rbind(
  run("small_full", "full", 2L, 5L),
  run("small_diagonal", "diagonal", 2L, 5L),
  run("moderate_full", "full", 3L, 20L),
  run("moderate_diagonal", "diagonal", 3L, 20L)
), row.names = FALSE)
cat("Phase 17P timings are regression signals; the public/internal difference is not asserted to be pure adapter overhead.\n")
