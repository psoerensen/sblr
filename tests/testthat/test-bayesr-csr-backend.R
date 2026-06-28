source_sblr_test_file <- function(path) {
  candidates <- c(path, file.path("..", "..", path))
  path <- candidates[file.exists(candidates)][1L]
  if (is.na(path)) stop("Could not find ", path, call. = FALSE)
  source(path)
}

if (!exists(".format_stblr_csr_bayesr_fit", mode = "function")) {
  source_sblr_test_file("R/sparse_ld_bed_helper.R")
}
if (!exists("check_stblr_backend_consistency", mode = "function")) {
  source_sblr_test_file("R/check-stblr-backend-consistency.R")
}
if (!exists("extract_stblr_finemap_loci", mode = "function")) {
  source_sblr_test_file("R/credible_sets.R")
  source_sblr_test_file("R/extract-stblr-finemap-loci.R")
}

make_bayesr_csr_raw <- function(nchains = 1L) {
  comp_prob <- matrix(
    c(
      0.70, 0.20, 0.10,
      0.10, 0.30, 0.60,
      0.40, 0.40, 0.20
    ),
    nrow = 3,
    byrow = TRUE
  )
  dm <- 1 - comp_prob[, 1L]
  bm <- c(0.01, -0.02, 0.03)
  bm_chain <- if (nchains == 1L) {
    list(bm)
  } else {
    list(c(0.00, -0.03, 0.02), c(0.02, -0.01, 0.04))
  }
  dm_chain <- if (nchains == 1L) {
    list(dm)
  } else {
    list(c(0.20, 0.80, 0.50), c(0.40, 1.00, 0.70))
  }
  bm_mat <- do.call(cbind, bm_chain)
  dm_mat <- do.call(cbind, dm_chain)
  trace_len <- 4L

  list(
    bm = matrix(rowMeans(bm_mat), ncol = 1),
    dm = matrix(rowMeans(dm_mat), ncol = 1),
    wy = matrix(rep(0, 3), ncol = 1),
    r = matrix(rep(0, 3), ncol = 1),
    b = matrix(bm, ncol = 1),
    component = matrix(c(0, 2, 1), ncol = 1),
    vbs = matrix(seq_len(trace_len) / 10, ncol = 1),
    vgs = matrix(seq_len(trace_len) / 9, ncol = 1),
    ves = matrix(seq_len(trace_len) / 8, ncol = 1),
    covb = matrix(1),
    covg = matrix(1),
    cove = matrix(1),
    vb = matrix(1),
    vg = matrix(1),
    ve = matrix(1),
    pi = matrix(c(0.40, 0.30, 0.30), nrow = 1),
    pim = matrix(c(0.45, 0.25, 0.30), nrow = 1),
    vle = matrix(seq_len(trace_len) / 7, ncol = 1),
    vld = matrix(seq_len(trace_len) / 6, ncol = 1),
    bm_sd = matrix(if (nchains == 1L) rep(0, 3) else apply(bm_mat, 1, stats::sd), ncol = 1),
    bm_min = matrix(apply(bm_mat, 1, min), ncol = 1),
    bm_max = matrix(apply(bm_mat, 1, max), ncol = 1),
    dm_sd = matrix(if (nchains == 1L) rep(0, 3) else apply(dm_mat, 1, stats::sd), ncol = 1),
    dm_min = matrix(apply(dm_mat, 1, min), ncol = 1),
    dm_max = matrix(apply(dm_mat, 1, max), ncol = 1),
    comp_prob = list(comp_prob),
    dm_component_mean = matrix(c(0.4, 1.5, 0.8), ncol = 1),
    ncomp = matrix(colSums(comp_prob), nrow = 1),
    mixture_var = c(0, 0.01, 0.1)
  )
}

format_bayesr_csr_test_fit <- function(raw, nchains) {
  fit <- .format_stblr_csr_bayesr_fit(
    raw,
    nt = 1L,
    m = 3L,
    trait_names = "D1",
    variable_names = paste0("m", 1:3),
    n_components = 3L
  )
  fit$input <- list(
    model = "bayesr",
    backend = "csr_bayesr",
    scheduled = FALSE,
    keep_chains = FALSE,
    updateLDswap = FALSE,
    nchains = nchains
  )
  fit
}

test_that("CSR BayesR formatter exposes standard and component fields", {
  raw <- make_bayesr_csr_raw(nchains = 1L)
  fit <- format_bayesr_csr_test_fit(raw, nchains = 1L)

  for (nm in c(
    "dm", "bm", "dm_sd", "dm_min", "dm_max",
    "bm_sd", "bm_min", "bm_max", "comp_prob", "dm_component_mean"
  )) {
    expect_true(nm %in% names(fit))
  }

  expect_equal(as.numeric(fit$dm[, "D1"]), 1 - raw$comp_prob[[1]][, 1], tolerance = 1e-12)
  expect_true(all(fit$dm >= -1e-12 & fit$dm <= 1 + 1e-12))
  expect_identical(rownames(fit$comp_prob$D1), paste0("m", 1:3))
  expect_identical(colnames(fit$comp_prob$D1), paste0("component_", 0:2))
  expect_equal(unname(rowSums(fit$comp_prob$D1)), rep(1, 3), tolerance = 1e-12)
  expect_true(all(fit$comp_prob$D1 >= -1e-12 & fit$comp_prob$D1 <= 1 + 1e-12))
  expect_equal(unname(fit$dm_component_mean[, "D1"]), c(0.4, 1.5, 0.8))
})

test_that("CSR BayesR formatter exposes single-chain summary convention", {
  fit <- format_bayesr_csr_test_fit(make_bayesr_csr_raw(nchains = 1L), nchains = 1L)

  for (nm in c("bm_sd", "bm_min", "bm_max", "dm_sd", "dm_min", "dm_max")) {
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

test_that("CSR BayesR formatter exposes finite multi-chain summaries", {
  fit <- format_bayesr_csr_test_fit(make_bayesr_csr_raw(nchains = 2L), nchains = 2L)

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

test_that("CSR BayesR formatted fit is compatible with fine-mapping extractor", {
  fit <- format_bayesr_csr_test_fit(make_bayesr_csr_raw(nchains = 1L), nchains = 1L)
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

test_that("CSR BayesR experimental helper rejects unsupported modes early", {
  expect_error(
    .stblr_csr_bayesr_experimental(stats = list(), scheduled = TRUE),
    "scheduled CSR BayesR"
  )
  expect_error(
    .stblr_csr_bayesr_experimental(stats = list(), updateLDswap = TRUE),
    "LD-swap"
  )
  expect_error(
    .stblr_csr_bayesr_experimental(stats = list(), keep_chains = TRUE),
    "keep_chains"
  )
})
