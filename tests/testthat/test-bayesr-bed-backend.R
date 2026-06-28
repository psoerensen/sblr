source_sblr_test_file <- function(path) {
  candidates <- c(path, file.path("..", "..", path))
  path <- candidates[file.exists(candidates)][1L]
  if (is.na(path)) stop("Could not find ", path, call. = FALSE)
  source(path)
}

if (!exists(".format_stblr_bayesr_fit", mode = "function")) {
  source_sblr_test_file("R/sparse_ld_bed_helper.R")
}
if (!exists("check_stblr_backend_consistency", mode = "function")) {
  source_sblr_test_file("R/check-stblr-backend-consistency.R")
}
if (!exists("extract_stblr_finemap_loci", mode = "function")) {
  source_sblr_test_file("R/credible_sets.R")
  source_sblr_test_file("R/extract-stblr-finemap-loci.R")
}

make_bayesr_bed_raw <- function(bm, comp_prob, component_mean,
                                bm_chain = NULL, dm_chain = NULL) {
  m <- nrow(comp_prob)
  k <- ncol(comp_prob)
  nt <- 1L
  trace_len <- 5L
  dm <- 1 - comp_prob[, 1L]

  if (is.null(bm_chain)) bm_chain <- list(as.numeric(bm))
  if (is.null(dm_chain)) dm_chain <- list(as.numeric(dm))
  nchains <- length(bm_chain)

  bm_sd <- if (nchains > 1L) {
    apply(do.call(cbind, bm_chain), 1L, stats::sd)
  } else {
    rep(0, m)
  }
  dm_sd <- if (nchains > 1L) {
    apply(do.call(cbind, dm_chain), 1L, stats::sd)
  } else {
    rep(0, m)
  }
  bm_min <- apply(do.call(cbind, bm_chain), 1L, min)
  bm_max <- apply(do.call(cbind, bm_chain), 1L, max)
  dm_min <- apply(do.call(cbind, dm_chain), 1L, min)
  dm_max <- apply(do.call(cbind, dm_chain), 1L, max)

  raw <- vector("list", 30L)
  raw[[1L]] <- list(as.numeric(bm))
  raw[[2L]] <- list(as.numeric(dm))
  raw[[3L]] <- list(rep(0, m))
  raw[[4L]] <- list(rep(0, m))
  raw[[5L]] <- list(as.numeric(bm))
  raw[[6L]] <- list(rep(0, m))
  raw[[7L]] <- list(seq.int(0L, m - 1L))
  raw[[8L]] <- list(seq_len(trace_len) / 10)
  raw[[9L]] <- list(seq_len(trace_len) / 9)
  raw[[10L]] <- list(seq_len(trace_len) / 8)
  for (i in 11:16) raw[[i]] <- list(1)
  raw[[17L]] <- list(rep(1 / k, k))
  raw[[18L]] <- list(rep(1 / k, k))
  raw[[19L]] <- list(c(-10, -10, 0, 0))
  raw[[20L]] <- list(c(2, 6))
  raw[[21L]] <- list(seq_len(trace_len) / 7)
  raw[[22L]] <- list(seq_len(trace_len) / 6)
  raw[[23L]] <- list(as.vector(comp_prob))
  raw[[24L]] <- list(bm_sd)
  raw[[25L]] <- list(bm_min)
  raw[[26L]] <- list(bm_max)
  raw[[27L]] <- list(dm_sd)
  raw[[28L]] <- list(dm_min)
  raw[[29L]] <- list(dm_max)
  raw[[30L]] <- list(as.numeric(component_mean))
  raw
}

format_bayesr_bed_test_fit <- function(raw, nchains) {
  format_fun <- if (exists(".format_stblr_bayesr_fit", mode = "function")) {
    .format_stblr_bayesr_fit
  } else {
    getFromNamespace(".format_stblr_bayesr_fit", "sblr")
  }
  fit <- format_fun(
    raw,
    nt = 1L,
    m = 3L,
    trait_names = "D1",
    variable_names = paste0("m", 1:3),
    n_components = 3L,
    keep_diagnostics = TRUE
  )
  fit$input <- list(
    model = "bayesr",
    backend = "bed_scheduled_chains_bayesr",
    scheduled = TRUE,
    keep_chains = FALSE,
    nchains = nchains
  )
  fit
}

test_that("BED BayesR experimental helper remains internal", {
  expect_true(exists(".stblr_bed_marker_bayesr_experimental", mode = "function"))
  helper_args <- names(formals(.stblr_bed_marker_bayesr_experimental))
  expect_true(all(c("mixture_var", "pi", "alpha", "nchains") %in% helper_args))
})

test_that("BED BayesR formatter exposes non-null PIP as standard dm", {
  comp_prob <- matrix(
    c(
      0.70, 0.20, 0.10,
      0.10, 0.30, 0.60,
      0.40, 0.40, 0.20
    ),
    nrow = 3,
    byrow = TRUE
  )
  raw <- make_bayesr_bed_raw(
    bm = c(0.01, -0.02, 0.03),
    comp_prob = comp_prob,
    component_mean = c(0.4, 1.5, 0.8)
  )

  fit <- format_bayesr_bed_test_fit(raw, nchains = 1L)

  expect_equal(as.numeric(fit$dm[, "D1"]), 1 - comp_prob[, 1], tolerance = 1e-12)
  expect_true(all(fit$dm >= -1e-12 & fit$dm <= 1 + 1e-12))
  expect_named(fit$comp_prob, "D1")
  # comp_prob is marker x component; component_0 is the null component.
  expect_identical(rownames(fit$comp_prob$D1), paste0("m", 1:3))
  expect_identical(colnames(fit$comp_prob$D1), paste0("component_", 0:2))
  expect_equal(unname(fit$comp_prob$D1), comp_prob, tolerance = 1e-12)
  expect_equal(unname(rowSums(fit$comp_prob$D1)), rep(1, 3), tolerance = 1e-12)
  expect_true(all(fit$comp_prob$D1 >= -1e-12 & fit$comp_prob$D1 <= 1 + 1e-12))
  expect_equal(unname(fit$dm_component_mean[, "D1"]), c(0.4, 1.5, 0.8))
})

test_that("BED BayesR formatter exposes single-chain summary convention", {
  comp_prob <- matrix(
    c(0.80, 0.10, 0.10, 0.25, 0.50, 0.25, 0.60, 0.30, 0.10),
    nrow = 3,
    byrow = TRUE
  )
  bm <- c(0.02, 0.10, -0.05)
  raw <- make_bayesr_bed_raw(bm, comp_prob, component_mean = c(0.2, 1, 0.5))
  fit <- format_bayesr_bed_test_fit(raw, nchains = 1L)

  for (nm in c("bm_sd", "bm_min", "bm_max", "dm_sd", "dm_min", "dm_max")) {
    expect_true(nm %in% names(fit))
    expect_equal(dim(fit[[nm]]), dim(fit$dm))
    expect_identical(rownames(fit[[nm]]), rownames(fit$dm))
    expect_identical(colnames(fit[[nm]]), colnames(fit$dm))
  }
  expect_equal(fit$bm_sd, fit$bm * 0, tolerance = 1e-12)
  expect_equal(fit$dm_sd, fit$dm * 0, tolerance = 1e-12)
  expect_equal(fit$bm_min, fit$bm, tolerance = 1e-12)
  expect_equal(fit$bm_max, fit$bm, tolerance = 1e-12)
  expect_equal(fit$dm_min, fit$dm, tolerance = 1e-12)
  expect_equal(fit$dm_max, fit$dm, tolerance = 1e-12)
})

test_that("BED BayesR formatter exposes finite multi-chain summaries", {
  comp_prob <- matrix(
    c(0.65, 0.25, 0.10, 0.20, 0.30, 0.50, 0.50, 0.25, 0.25),
    nrow = 3,
    byrow = TRUE
  )
  bm_chain <- list(c(0.00, 0.10, -0.10), c(0.04, 0.20, -0.04))
  dm_chain <- list(c(0.30, 0.70, 0.40), c(0.40, 0.90, 0.60))
  raw <- make_bayesr_bed_raw(
    bm = rowMeans(do.call(cbind, bm_chain)),
    comp_prob = comp_prob,
    component_mean = c(0.45, 1.4, 0.75),
    bm_chain = bm_chain,
    dm_chain = dm_chain
  )
  fit <- format_bayesr_bed_test_fit(raw, nchains = 2L)

  expect_true(all(is.finite(fit$bm_sd)))
  expect_true(all(is.finite(fit$dm_sd)))
  expect_true(all(fit$bm_sd >= -1e-12))
  expect_true(all(fit$dm_sd >= -1e-12))
  expect_true(all(fit$bm_min <= fit$bm + 1e-12))
  expect_true(all(fit$bm <= fit$bm_max + 1e-12))
  expect_true(all(fit$dm_min <= fit$dm + 1e-12))
  expect_true(all(fit$dm <= fit$dm_max + 1e-12))

  chk <- check_stblr_backend_consistency(
    fit,
    require_chain_summaries = TRUE,
    verbose = FALSE
  )
  expect_true(chk$ok)
})

test_that("BED BayesR formatted fit is compatible with fine-mapping extractor", {
  comp_prob <- matrix(
    c(0.80, 0.10, 0.10, 0.25, 0.50, 0.25, 0.60, 0.30, 0.10),
    nrow = 3,
    byrow = TRUE
  )
  bm <- c(0.02, 0.10, -0.05)
  raw <- make_bayesr_bed_raw(bm, comp_prob, component_mean = c(0.2, 1, 0.5))
  fit <- format_bayesr_bed_test_fit(raw, nchains = 1L)
  glist <- list(
    rsids = list(paste0("m", 1:3)),
    chr = list(rep(1L, 3)),
    pos = list(c(100, 200, 300)),
    sparseLD = list(chr = 1L, cls = list(1:3), prefix = NULL)
  )

  fm <- extract_stblr_finemap_loci(
    fit = fit,
    Glist = glist,
    locus_sets = list(regionA = paste0("m", 1:3)),
    trait = "D1",
    credible_sets = FALSE
  )

  expect_equal(fm$markers$pip_sd, as.numeric(fit$dm_sd[, "D1"]))
  expect_equal(fm$markers$pip_min, as.numeric(fit$dm_min[, "D1"]))
  expect_equal(fm$markers$pip_max, as.numeric(fit$dm_max[, "D1"]))
  expect_equal(fm$markers$bm_sd, as.numeric(fit$bm_sd[, "D1"]))
  expect_equal(fm$markers$bm_min, as.numeric(fit$bm_min[, "D1"]))
  expect_equal(fm$markers$bm_max, as.numeric(fit$bm_max[, "D1"]))
})
