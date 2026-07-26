phase17s_public_args <- function(case, nchains = 1L, ncores = 1L,
                                 chain_seeds = NULL, keep_chains = FALSE,
                                 residual_covariance = "full",
                                 updates = FALSE,
                                 convergence = "none", ...) {
  args <- phase17p_public_args(
    case, residual_covariance, updates, center = FALSE,
    nchains = nchains, ncores = ncores, chain_seeds = chain_seeds,
    keep_chains = keep_chains, convergence = convergence, ...)
  args
}

phase17s_internal_args <- function(public_args) {
  args <- phase17p_native_args(public_args)
  args$nchains <- as.integer(phase17p_or(public_args$nchains, 1L))
  args$ncores <- as.integer(phase17p_or(public_args$ncores, 1L))
  seeds <- phase17p_or(public_args$chain_seeds, NULL)
  args$chain_seeds <- if (is.null(seeds)) integer() else seeds
  args$keep_chains <- phase17p_or(public_args$keep_chains, FALSE)
  args
}

phase17s_internal <- function(public_args) {
  do.call(sblr:::mtblr_bed_chains_internal,
          phase17s_internal_args(public_args))
}

phase17s_fit_numerics <- function(fit) {
  fields <- c("bm", "dm", "wy", "r", "b", "d", "marker_order",
              "vbs", "vgs", "ves", "covb", "covg", "cove",
              "vb", "vg", "ve", "pi", "pim",
              "bm_sd", "bm_min", "bm_max", "dm_sd", "dm_min", "dm_max")
  setNames(lapply(fields, function(field) unname(fit[[field]])), fields)
}

phase17s_raw_numerics <- function(raw) {
  values <- list(
    bm = raw$marker$bm, dm = raw$marker$dm, wy = raw$marker$wy,
    r = raw$marker$r, b = raw$marker$b, d = raw$marker$state,
    marker_order = raw$marker$order,
    vbs = raw$trace$vbs, vgs = raw$trace$vgs, ves = raw$trace$ves,
    covb = raw$variance$covb, covg = raw$variance$covg,
    cove = raw$variance$cove, vb = raw$variance$vb,
    vg = raw$variance$vg, ve = raw$variance$ve,
    pi = raw$pi$final, pim = raw$pi$mean,
    bm_sd = raw$marker$bm_sd, bm_min = raw$marker$bm_min,
    bm_max = raw$marker$bm_max, dm_sd = raw$marker$dm_sd,
    dm_min = raw$marker$dm_min, dm_max = raw$marker$dm_max)
  lapply(values, unname)
}

phase17s_compare_public_internal <- function(args, tolerance = 1e-12) {
  fit <- do.call(mtblr_bed, args)
  raw <- phase17s_internal(args)
  testthat::expect_equal(phase17s_fit_numerics(fit),
                         phase17s_raw_numerics(raw), tolerance = tolerance)
  testthat::expect_identical(fit$nchains, raw$meta$nchains)
  testthat::expect_identical(fit$chain_seeds,
                             raw$diagnostics$mt_bed$chain_seeds)
  testthat::expect_identical(fit$chain_diagnostics$used_workers,
                             raw$diagnostics$mt_bed$used_workers)
  testthat::expect_identical(is.null(fit$chains), is.null(raw$chains))
  if (!is.null(fit$chains)) {
    testthat::expect_identical(names(fit$chains), names(raw$chains))
  }
  invisible(list(fit = fit, raw = raw))
}

phase17s_without_timing <- function(fit) {
  fit$bed_diagnostics[c("chain_seconds", "seconds_mean", "seconds_max",
                        "dispatch_seconds", "requested_cores",
                        "used_workers")] <- NULL
  fit$chain_diagnostics[c("chain_seconds", "seconds_mean", "seconds_max",
                          "dispatch_seconds", "requested_cores",
                          "used_workers")] <- NULL
  if (!is.null(fit$chains)) {
    fit$chains <- lapply(fit$chains, function(chain) {
      chain$diagnostics$seconds <- NULL
      chain
    })
  }
  fit$input$used_workers <- NULL
  fit$input$ncores <- NULL
  fit$input$ncores_requested <- NULL
  fit$input$memory_estimate <- NULL
  fit$input$memory_warning_gb <- NULL
  fit$memory_estimate <- NULL
  fit
}

phase17s_preexisting_fields <- function(fit) {
  fields <- c("bm", "dm", "wy", "r", "b", "d", "marker_order",
              "vbs", "vgs", "ves", "covb", "covg", "cove",
              "vb", "vg", "ve", "pi", "pim", "rb", "rg", "re",
              "bed_diagnostics", "phenotype_preprocessing")
  fit[fields]
}
