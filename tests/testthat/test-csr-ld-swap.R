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
})

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
