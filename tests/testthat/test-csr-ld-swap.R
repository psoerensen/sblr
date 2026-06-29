make_tiny_csr_prefix <- function() {
  prefix <- tempfile("tiny_ld_swap_csr_")
  row_ptr <- c(0, 1, 1, 1)
  col_idx <- 1L
  values <- 0.95

  sblr:::.stblr_write_uint64_file(paste0(prefix, ".row_ptr.u64.bin"), row_ptr)
  sblr:::.stblr_write_uint32_file(paste0(prefix, ".col_idx.u32.0based.bin"), col_idx)
  writeBin(as.numeric(values), paste0(prefix, ".values.f32.bin"),
           size = 4, endian = "little")

  writeLines(
    c(
      "format=sparse_ld_csr",
      "storage=streamed_upper_triangle",
      "n_bed=NA",
      "n_used=NA",
      "n_samples_used=NA",
      "n_variants=3",
      "nnz=1",
      "triangle=upper",
      "diagonal=implicit_1",
      paste0("row_ptr_file=", prefix, ".row_ptr.u64.bin"),
      paste0("col_idx_file=", prefix, ".col_idx.u32.0based.bin"),
      paste0("values_file=", prefix, ".values.f32.bin"),
      "row_ptr_type=uint64",
      "col_idx_type=uint32",
      "values_type=float32",
      "index_base=0",
      "value=r"
    ),
    paste0(prefix, ".meta.txt")
  )

  prefix
}

tiny_csr_stats <- function() {
  markers <- paste0("m", 1:3)
  list(
    wy = list(trait1 = stats::setNames(c(8, 7.5, 0.2), markers)),
    ww = list(trait1 = stats::setNames(rep(100, 3), markers)),
    yy = stats::setNames(100, "trait1"),
    n = 100L,
    m = 3L,
    marker_names = markers,
    trait_names = "trait1"
  )
}

test_that("CSR LD-swap arguments are validated", {
  stats <- tiny_csr_stats()
  prefix <- make_tiny_csr_prefix()

  expect_error(
    stblr_csr(stats = stats, ld_prefix = prefix, updateLDswap = NA),
    "updateLDswap"
  )
  expect_error(
    stblr_csr(stats = stats, ld_prefix = prefix, ld_swap_prob = 2),
    "ld_swap_prob"
  )
  expect_error(
    stblr_csr(stats = stats, ld_prefix = prefix, ld_swap_r2 = -0.1),
    "ld_swap_r2"
  )
  expect_error(
    stblr_csr(stats = stats, ld_prefix = prefix, ld_swap_max_friends = 0),
    "ld_swap_max_friends"
  )
  expect_error(
    stblr_csr(stats = stats, ld_prefix = prefix, ld_swap_moves = -1),
    "ld_swap_moves"
  )
  expect_error(
    stblr_csr(stats = stats, ld_prefix = prefix, nchains = 0),
    "nchains"
  )
  expect_error(
    stblr_csr(stats = stats, ld_prefix = prefix, nchains = -1),
    "nchains"
  )
  expect_error(
    stblr_csr(stats = stats, ld_prefix = prefix, keep_chains = NA),
    "keep_chains"
  )
  expect_error(
    stblr_csr(stats = stats, ld_prefix = prefix, nchains = 2, chain_seeds = 1),
    "chain_seeds"
  )
  expect_error(
    stblr_csr(stats = stats, ld_prefix = prefix, scheduled = TRUE, keep_chains = TRUE),
    "keep_chains"
  )
  expect_error(
    stblr_csr(stats = stats, ld_prefix = prefix, scheduled = TRUE, updateLDswap = TRUE),
    "updateLDswap"
  )
  expect_error(
    stblr_csr(stats = stats, ld_prefix = prefix, method = "invalid"),
    "method"
  )
})

expect_csr_chain_summary_shape <- function(fit) {
  for (nm in c("dm_sd", "dm_min", "dm_max", "bm_sd", "bm_min", "bm_max")) {
    expect_true(nm %in% names(fit))
    expect_equal(dim(fit[[nm]]), dim(fit$dm))
    expect_identical(rownames(fit[[nm]]), rownames(fit$dm))
    expect_identical(colnames(fit[[nm]]), colnames(fit$dm))
    expect_false(anyNA(fit[[nm]]))
    expect_true(all(is.finite(fit[[nm]])))
  }
}

tiny_csr_glist <- function() {
  list(
    rsids = list(paste0("m", 1:3)),
    chr = list(rep(1L, 3)),
    pos = list(c(100, 200, 300)),
    sparseLD = list(chr = 1L, cls = list(1:3), prefix = make_tiny_csr_prefix())
  )
}

test_that("CSR sampler returns zero LD-swap diagnostics by default", {
  fit <- stblr_csr(
    stats = tiny_csr_stats(),
    ld_prefix = make_tiny_csr_prefix(),
    pi_init = 0.5,
    pi_prior_mean = 0.5,
    pi_prior_strength = 2,
    updateB = FALSE,
    updateE = FALSE,
    updatePi = FALSE,
    nit = 2,
    nburn = 0,
    seed = 11
  )

  expect_true("ld_swap" %in% names(fit))
  expect_equal(fit$ld_swap$attempted, 0)
  expect_equal(fit$ld_swap$accepted, 0)
  expect_equal(fit$ld_swap$acceptance_rate, 0)
  expect_equal(fit$input$method, "bayesc")
  expect_equal(fit$input$model, "bayesc")
  expect_equal(fit$input$backend, "csr_bayesc")
  expect_equal(fit$input$data_level, "summary")
  expect_equal(fit$input$scheduled, FALSE)
  expect_equal(fit$input$nchains, 1L)
})

test_that("CSR method defaults to BayesC and explicit BayesC preserves metadata", {
  args <- list(
    stats = tiny_csr_stats(),
    ld_prefix = make_tiny_csr_prefix(),
    pi_init = 0.5,
    pi_prior_mean = 0.5,
    pi_prior_strength = 2,
    updateB = FALSE,
    updateE = FALSE,
    updatePi = FALSE,
    nit = 2,
    nburn = 0,
    seed = 11
  )

  fit_default <- do.call(stblr_csr, args)
  fit_method <- do.call(stblr_csr, c(args, list(method = "BayesC")))

  expect_equal(class(fit_method), class(fit_default))
  expect_equal(names(fit_method), names(fit_default))
  expect_equal(dim(fit_method$dm), dim(fit_default$dm))
  expect_equal(fit_default$input$method, "bayesc")
  expect_equal(fit_method$input$method, "bayesc")
  expect_equal(fit_method$input$model, "bayesc")
  expect_equal(fit_method$input$backend, "csr_bayesc")
  expect_equal(fit_method$input$data_level, "summary")
})

test_that("CSR sampler runs with LD-swap enabled on tiny CSR LD", {
  fit <- stblr_csr(
    stats = tiny_csr_stats(),
    ld_prefix = make_tiny_csr_prefix(),
    pi_init = 0.5,
    pi_prior_mean = 0.5,
    pi_prior_strength = 2,
    updateB = FALSE,
    updateE = FALSE,
    updatePi = FALSE,
    nit = 3,
    nburn = 0,
    seed = 12,
    updateLDswap = TRUE,
    ld_swap_prob = 1,
    ld_swap_r2 = 0.8,
    ld_swap_max_friends = 2,
    ld_swap_moves = 2
  )

  expect_true("ld_swap" %in% names(fit))
  expect_equal(fit$ld_swap$attempted, 6)
  expect_true(fit$ld_swap$accepted >= 0)
  expect_true(fit$ld_swap$acceptance_rate >= 0)
  expect_true(fit$ld_swap$acceptance_rate <= 1)
})

test_that("CSR sampler returns multi-chain posterior summaries", {
  fit <- stblr_csr(
    stats = tiny_csr_stats(),
    ld_prefix = make_tiny_csr_prefix(),
    pi_init = 0.5,
    pi_prior_mean = 0.5,
    pi_prior_strength = 2,
    updateB = FALSE,
    updateE = FALSE,
    updatePi = FALSE,
    nit = 3,
    nburn = 0,
    seed = 14,
    nchains = 2,
    ncores = 2
  )

  expect_equal(fit$input$nchains, 2L)
  expect_csr_chain_summary_shape(fit)
})

test_that("scheduled CSR returns one-chain posterior summaries", {
  fit <- stblr_csr(
    stats = tiny_csr_stats(),
    ld_prefix = make_tiny_csr_prefix(),
    scheduled = TRUE,
    pi_init = 0.5,
    pi_prior_mean = 0.5,
    pi_prior_strength = 2,
    updateB = FALSE,
    updateE = FALSE,
    updatePi = FALSE,
    full_sweep_every = 1,
    null_skip_base = 1,
    nit = 2,
    nburn = 0,
    seed = 17,
    nchains = 1,
    ncores = 1
  )

  expect_equal(fit$input$nchains, 1L)
  expect_equal(fit$input$backend, "csr_scheduled_bayesc")
  expect_equal(fit$input$data_level, "summary")
  expect_equal(fit$input$scheduled, TRUE)
  expect_true("pis" %in% names(fit))
  expect_csr_chain_summary_shape(fit)
  expect_equal(fit$dm_sd, fit$dm * 0, tolerance = 1e-12)
  expect_equal(fit$bm_sd, fit$bm * 0, tolerance = 1e-12)
  expect_equal(fit$dm_min, fit$dm, tolerance = 1e-12)
  expect_equal(fit$dm_max, fit$dm, tolerance = 1e-12)
  expect_equal(fit$bm_min, fit$bm, tolerance = 1e-12)
  expect_equal(fit$bm_max, fit$bm, tolerance = 1e-12)
})

test_that("scheduled CSR returns multi-chain posterior summaries", {
  fit <- stblr_csr(
    stats = tiny_csr_stats(),
    ld_prefix = make_tiny_csr_prefix(),
    scheduled = TRUE,
    pi_init = 0.5,
    pi_prior_mean = 0.5,
    pi_prior_strength = 2,
    updateB = FALSE,
    updateE = FALSE,
    updatePi = FALSE,
    full_sweep_every = 1,
    null_skip_base = 1,
    nit = 3,
    nburn = 0,
    seed = 18,
    nchains = 2,
    ncores = 2
  )

  expect_equal(fit$input$nchains, 2L)
  expect_equal(fit$input$backend, "csr_scheduled_bayesc")
  expect_equal(fit$input$data_level, "summary")
  expect_equal(fit$input$scheduled, TRUE)
  expect_true("pis" %in% names(fit))
  expect_csr_chain_summary_shape(fit)
  expect_true(all(fit$dm_sd >= -1e-12))
  expect_true(all(fit$bm_sd >= -1e-12))
  expect_true(all(fit$dm_min <= fit$dm + 1e-12))
  expect_true(all(fit$dm <= fit$dm_max + 1e-12))
  expect_true(all(fit$bm_min <= fit$bm + 1e-12))
  expect_true(all(fit$bm <= fit$bm_max + 1e-12))

  fm <- extract_stblr_finemap_loci(
    fit = fit,
    Glist = tiny_csr_glist(),
    locus_sets = list(regionA = rownames(fit$dm)),
    trait = "trait1",
    credible_sets = FALSE
  )
  expect_equal(fm$markers$pip_sd, as.numeric(fit$dm_sd[, "trait1"]))
  expect_equal(fm$markers$pip_min, as.numeric(fit$dm_min[, "trait1"]))
  expect_equal(fm$markers$pip_max, as.numeric(fit$dm_max[, "trait1"]))
  expect_equal(fm$markers$bm_sd, as.numeric(fit$bm_sd[, "trait1"]))
  expect_equal(fm$markers$bm_min, as.numeric(fit$bm_min[, "trait1"]))
  expect_equal(fm$markers$bm_max, as.numeric(fit$bm_max[, "trait1"]))
})

test_that("CSR sampler can keep compact per-chain summaries", {
  fit <- stblr_csr(
    stats = tiny_csr_stats(),
    ld_prefix = make_tiny_csr_prefix(),
    pi_init = 0.5,
    pi_prior_mean = 0.5,
    pi_prior_strength = 2,
    updateB = FALSE,
    updateE = FALSE,
    updatePi = FALSE,
    nit = 3,
    nburn = 0,
    seed = 15,
    nchains = 2,
    keep_chains = TRUE,
    ncores = 2
  )

  expect_true("chains" %in% names(fit))
  expect_length(fit$chains$trait1, 2)
  dm_chain <- vapply(fit$chains$trait1, function(x) x$dm, numeric(3))
  bm_chain <- vapply(fit$chains$trait1, function(x) x$bm, numeric(3))
  expect_equal(as.numeric(fit$dm[, "trait1"]), unname(rowMeans(dm_chain)), tolerance = 1e-12)
  expect_equal(as.numeric(fit$bm[, "trait1"]), unname(rowMeans(bm_chain)), tolerance = 1e-12)
})

test_that("CSR LD-swap diagnostics aggregate across chains", {
  fit <- stblr_csr(
    stats = tiny_csr_stats(),
    ld_prefix = make_tiny_csr_prefix(),
    pi_init = 0.5,
    pi_prior_mean = 0.5,
    pi_prior_strength = 2,
    updateB = FALSE,
    updateE = FALSE,
    updatePi = FALSE,
    nit = 3,
    nburn = 0,
    seed = 16,
    nchains = 2,
    keep_chains = TRUE,
    ncores = 2,
    updateLDswap = TRUE,
    ld_swap_prob = 1,
    ld_swap_r2 = 0.8,
    ld_swap_max_friends = 2,
    ld_swap_moves = 2
  )

  expect_true("ld_swap" %in% names(fit))
  expect_equal(fit$ld_swap$attempted, 12)
  expect_true(fit$ld_swap$accepted >= 0)
  expect_true(fit$ld_swap$acceptance_rate >= 0)
  expect_true(fit$ld_swap$acceptance_rate <= 1)
  expect_true("ld_swap_chains" %in% names(fit))
  expect_equal(sum(fit$ld_swap_chains$trait1$attempted), fit$ld_swap$attempted)
})

test_that("local CSR runner accepts LD-swap pass-through arguments", {
  fit <- sblr:::.stblr_run_local_csr(
    stats = tiny_csr_stats(),
    ld_prefix = make_tiny_csr_prefix(),
    ve = 1,
    vb = 0.1,
    pi = 0.5,
    nit = 2,
    nburn = 0,
    seed = 13,
    updateLDswap = TRUE,
    ld_swap_prob = 1,
    ld_swap_r2 = 0.8,
    ld_swap_max_friends = 2,
    ld_swap_moves = 1
  )

  expect_true("ld_swap" %in% names(fit))
  expect_true(isTRUE(fit$input$updateLDswap))
  expect_equal(fit$input$ld_swap_moves, 1L)
})
