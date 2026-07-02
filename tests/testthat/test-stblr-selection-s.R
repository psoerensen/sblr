make_selection_s_csr_prefix <- function() {
  prefix <- tempfile("tiny_selection_s_csr_")
  row_ptr <- c(0, 1, 1, 1)
  col_idx <- 1L
  values <- 0.4

  sblr:::.stblr_write_uint64_file(paste0(prefix, ".row_ptr.u64.bin"), row_ptr)
  sblr:::.stblr_write_uint32_file(
    paste0(prefix, ".col_idx.u32.0based.bin"),
    col_idx
  )
  writeBin(
    as.numeric(values),
    paste0(prefix, ".values.f32.bin"),
    size = 4,
    endian = "little"
  )

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

selection_s_csr_stats <- function() {
  markers <- paste0("m", 1:3)
  list(
    wy = list(trait1 = stats::setNames(c(8, 3, 0.4), markers)),
    ww = list(trait1 = stats::setNames(rep(100, 3), markers)),
    yy = stats::setNames(100, "trait1"),
    n = 100L,
    m = 3L,
    marker_names = markers,
    trait_names = "trait1"
  )
}

make_selection_s_csr_fit <- function() {
  stblr_csr(
    stats = selection_s_csr_stats(),
    ld_prefix = make_selection_s_csr_prefix(),
    nit = 8,
    nburn = 2,
    nthin = 1,
    ncores = 1,
    seed = 11,
    pi_init = 0.5,
    pi_prior_mean = 0.5,
    pi_prior_strength = 2,
    updateB = FALSE,
    updateE = FALSE,
    updatePi = FALSE
  )
}

test_that("summarise_stblr_maf_architecture works for a small CSR fit", {
  fit <- make_selection_s_csr_fit()

  out <- summarise_stblr_maf_architecture(
    fit,
    maf = c(m1 = 0.05, m2 = 0.2, m3 = 0.4)
  )

  expect_s3_class(out, "data.frame")
  expect_equal(out$trait, "trait1")
  expect_equal(out$method, "posthoc_regression")
  expect_equal(out$response, "log_pip_weighted_bm2")
  expect_equal(out$n_markers, 3)
  expect_equal(out$n_effective_markers, 3)
  expect_true(all(is.finite(out$selection_s_posthoc)))
  expect_true(all(is.finite(out$intercept)))
})

test_that("summarise_stblr_maf_architecture handles h and unweighted response", {
  fit <- list(
    dm = matrix(c(0.1, 0.5, 0.9), ncol = 1,
                dimnames = list(paste0("m", 1:3), "trait1")),
    bm = matrix(c(0.01, -0.02, 0.03), ncol = 1,
                dimnames = list(paste0("m", 1:3), "trait1"))
  )
  fit_before <- fit

  out <- summarise_stblr_maf_architecture(
    fit,
    h = c(m1 = 0.095, m2 = 0.32, m3 = 0.48),
    use_pip_weights = FALSE
  )

  expect_equal(out$method, "posthoc_regression")
  expect_equal(out$response, "log_bm2")
  expect_true(all(is.finite(out$selection_s_posthoc)))
  expect_identical(fit, fit_before)
})

test_that("summarise_stblr_maf_architecture handles missing MAF or h clearly", {
  fit <- list(
    dm = matrix(c(0.1, 0.5, 0.9), ncol = 1),
    bm = matrix(c(0.01, -0.02, 0.03), ncol = 1)
  )

  expect_error(
    summarise_stblr_maf_architecture(fit),
    "Either maf or h must be supplied"
  )
  expect_error(
    summarise_stblr_maf_architecture(fit, maf = c(0.1, 0.2)),
    "length equal to the number of markers"
  )
  expect_error(
    summarise_stblr_maf_architecture(fit, h = c(0.1, 0.2, 0.7)),
    "h must contain heterozygosity values"
  )
})
