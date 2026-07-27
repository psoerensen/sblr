phase17p_case <- function(nt = 2L, use_all_rows = FALSE,
                          matched_ids = FALSE, uncentered = FALSE) {
  fixture <- phase17n_fixture()
  full_ids <- c(paste0("a", 1:4), paste0("b", 1:3))
  Glist <- list(
    n = fixture$n_bed, ids = paste0("id", seq_len(fixture$n_bed)),
    bedfiles = fixture$paths,
    rsids = list(full_ids[1:4], full_ids[5:7]),
    rsidsLD = list(c("a4", "a2", "a1"), c("b3", "b1")),
    af = list(c(.27, .43, .20, .31), c(.46, .20, .38))
  )
  rows <- if (use_all_rows) NULL else fixture$rows
  n <- if (use_all_rows) fixture$n_bed else length(fixture$rows)
  Y <- vapply(seq_len(nt), function(trait) {
    x <- sin(seq_len(n) * (.31 + trait / 11)) +
      cos(seq_len(n) * (.19 + trait / 17))
    if (!uncentered) x <- x - mean(x)
    x
  }, numeric(n))
  colnames(Y) <- paste0("T", seq_len(nt))
  if (matched_ids) {
    stopifnot(!use_all_rows)
    rownames(Y) <- Glist$ids[fixture$rows]
    rows <- NULL
  }
  list(fixture = fixture, Glist = Glist, Y = Y, rows = rows)
}

phase17p_cleanup <- function(case) {
  unlink(case$fixture$paths)
  invisible()
}

phase17p_or <- function(x, y) if (is.null(x)) y else x

phase17p_public_args <- function(case, residual_covariance = "full",
                                  updates = FALSE, center = FALSE, ...) {
  args <- list(
    y = case$Y, Glist = case$Glist, rows = case$rows,
    residual_covariance = residual_covariance, center = center,
    updateB = updates, updateE = updates, updatePi = updates,
    nit = 5L, nburn = 2L, nthin = 1L, seed = 17016L,
    memory_warning_gb = Inf)
  extra <- list(...)
  args[names(extra)] <- extra
  args
}

phase17p_native_args <- function(public_args) {
  y <- public_args$y
  dat <- sblr:::.make_bed_marker_data(
    public_args$Glist, y, phase17p_or(public_args$chr, NULL),
    phase17p_or(public_args$cls, NULL),
    phase17p_or(public_args$block_size, 1000L),
    phase17p_or(public_args$rows, NULL))
  Y <- as.matrix(dat$y)
  if (isTRUE(phase17p_or(public_args$center, TRUE))) {
    Y <- sweep(Y, 2L, colMeans(Y), "-")
  }
  nt <- ncol(Y); m <- dat$m
  mod <- sblr:::.mtblr_models(
    phase17p_or(public_args$models, NULL),
    phase17p_or(public_args$pimodels, NULL),
    phase17p_or(public_args$pi, .001), nt)
  null <- which(rowSums(mod$matrix) == 0L)
  p_active <- 1 - sum(mod$probabilities[null])
  labels <- unique(dat$sets)
  defaults <- lapply(labels, function(label) which(dat$sets == label))
  sets <- sblr:::.mtblr_sets(phase17p_or(public_args$sets, defaults), m)
  h2 <- phase17p_or(public_args$h2, .5)
  if (length(h2) == 1L) h2 <- rep(h2, nt)
  vy <- colSums(Y^2) / (nrow(Y) - 1)
  vg <- phase17p_or(public_args$vg, diag(vy * h2, nt))
  ve <- phase17p_or(public_args$ve, diag(vy * (1 - h2), nt))
  vb <- phase17p_or(public_args$vb, vg / (m * p_active))
  nub <- phase17p_or(public_args$nub, 4)
  nue <- phase17p_or(public_args$nue, 4)
  ssb <- phase17p_or(public_args$ssb_prior,
                     ((nub - 2) / nub) * vg / (m * p_active))
  sse <- phase17p_or(public_args$sse_prior,
                     ((nue - 2) / nue) * ve)
  init <- sblr:::.mtblr_bed_initialization(
    phase17p_or(public_args$beta, NULL),
    phase17p_or(public_args$b, NULL),
    phase17p_or(public_args$state, NULL), mod$matrix, m, nt)
  list(
    bed_files = dat$bed_files, n_bed = dat$n_total, cls = dat$cls,
    rows = dat$rows, af = unlist(dat$af, use.names = FALSE), Y = Y,
    beta_init = lapply(seq_len(nt), function(t) init$beta[, t]),
    b_init = lapply(seq_len(nt), function(t) init$b[, t]),
    state_init = lapply(seq_len(nt), function(t) init$state[, t]),
    sets = sets$native, B = vb, E = ve,
    ssb_prior = lapply(seq_len(nt), function(t) ssb[t, ]),
    sse_prior = lapply(seq_len(nt), function(t) sse[t, ]),
    models = mod$native, pi = mod$probabilities, nub = nub, nue = nue,
    updateB = phase17p_or(public_args$updateB, TRUE),
    updateE = phase17p_or(public_args$updateE, TRUE),
    updatePi = phase17p_or(public_args$updatePi, TRUE),
    residual_covariance = phase17p_or(public_args$residual_covariance, "full"),
    nit = phase17p_or(public_args$nit, 1000L),
    nburn = phase17p_or(public_args$nburn, 500L),
    nthin = phase17p_or(public_args$nthin, 1L),
    seed = phase17p_or(public_args$seed, 1L),
    method = 4L)
}

phase17p_compare_public_internal <- function(args, tolerance = 1e-12) {
  fit <- do.call(mtblr_bed, args)
  raw <- do.call(sblr:::mtblr_bed_internal, phase17p_native_args(args))
  marker_fields <- c("bm", "dm", "wy", "r", "b")
  trace_fields <- c("vbs", "vgs", "ves")
  variance_fields <- c(
    covb = "cov_b_mean", covg = "cov_g_mean", cove = "cov_e_mean",
    vb = "cov_b_final", vg = "cov_g_final", ve = "cov_e_final")
  actual <- c(
    lapply(marker_fields, function(field) unname(fit[[field]])),
    lapply(trace_fields, function(field) unname(fit[[field]])),
    lapply(unname(variance_fields), function(field) unname(fit[[field]])),
    list(pi = unname(fit$pi_final), pim = unname(fit$pi_mean)))
  names(actual) <- c(marker_fields, trace_fields, names(variance_fields),
                     "pi", "pim")
  expected <- c(
    raw$marker[marker_fields], raw$trace[trace_fields],
    raw$variance[names(variance_fields)],
    list(pi = raw$pi$final, pim = raw$pi$mean))
  testthat::expect_equal(actual, expected, tolerance = tolerance)
  testthat::expect_identical(unname(fit$d), raw$marker$state)
  testthat::expect_identical(fit$marker_order, raw$marker$order)
  testthat::expect_equal(
    list(fit$raw_schema_version, fit$input$backend, fit$input$data_level),
    list(1L, "mt_bed_bayesc", "individual"))
  serial_diagnostics <- names(raw$diagnostics$mt_bed)
  testthat::expect_equal(fit$diagnostics$mt_bed[serial_diagnostics],
                         raw$diagnostics$mt_bed)
  invisible(fit)
}
