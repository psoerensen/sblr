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

make_selection_s_subset_fit <- function(n_markers = 20) {
  markers <- paste0("m", seq_len(n_markers))
  trait_names <- c("D1", "D2")
  pip_d1 <- seq(0.001, 0.2, length.out = n_markers)
  pip_d2 <- rev(pip_d1)
  effect_d1 <- seq(0.01, 0.08, length.out = n_markers)
  effect_d2 <- rev(effect_d1)
  dm <- cbind(D1 = pip_d1, D2 = pip_d2)
  bm <- cbind(D1 = effect_d1, D2 = effect_d2)
  rownames(dm) <- markers
  rownames(bm) <- markers
  list(
    dm = dm,
    bm = bm,
    maf = stats::setNames(seq(0.05, 0.45, length.out = n_markers), markers),
    markers = markers,
    trait_names = trait_names
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
  expect_true(is.na(out$selection_s_posthoc))
  expect_true(is.na(out$intercept))
})

test_that("summarise_stblr_maf_architecture handles h and unweighted response", {
  fit <- list(
    dm = matrix(seq(0.1, 0.9, length.out = 6), ncol = 1,
                dimnames = list(paste0("m", 1:6), "trait1")),
    bm = matrix(seq(0.01, 0.06, length.out = 6), ncol = 1,
                dimnames = list(paste0("m", 1:6), "trait1"))
  )
  fit_before <- fit

  out <- summarise_stblr_maf_architecture(
    fit,
    h = stats::setNames(seq(0.095, 0.48, length.out = 6), paste0("m", 1:6)),
    use_pip_weights = FALSE
  )

  expect_equal(out$method, "posthoc_regression")
  expect_equal(out$response, "log_bm2")
  expect_true(all(is.finite(out$selection_s_posthoc)))
  expect_identical(fit, fit_before)
})

test_that("summarise_stblr_maf_architecture supports all-marker and filtered subsets", {
  fit <- make_selection_s_subset_fit()

  out_all <- summarise_stblr_maf_architecture(fit, maf = fit$maf)
  expect_equal(out_all$n_markers, c(20L, 20L))
  expect_equal(out_all$marker_filter, c("all", "all"))
  expect_true(all(c("se", "p_value", "r2", "marker_filter") %in% names(out_all)))
  expect_true(all(is.finite(out_all$selection_s_posthoc)))

  out_min <- summarise_stblr_maf_architecture(fit, maf = fit$maf, min_pip = 0.01)
  expect_equal(out_min$marker_filter, rep("min_pip=0.01", 2))
  expect_true(all(out_min$n_markers < out_all$n_markers))
  expect_true(all(out_min$n_markers >= 5))

  out_top <- summarise_stblr_maf_architecture(fit, maf = fit$maf, top_n = 10)
  expect_equal(out_top$n_markers, c(10L, 10L))
  expect_equal(out_top$marker_filter, rep("top_n=10", 2))
})

test_that("summarise_stblr_maf_architecture supports character marker subsets", {
  fit <- make_selection_s_subset_fit()
  marker_subset <- paste0("m", 1:8)

  out <- summarise_stblr_maf_architecture(
    fit,
    maf = unname(fit$maf),
    markers = marker_subset
  )

  expect_equal(out$n_markers, c(8L, 8L))
  expect_equal(out$n_effective_markers, c(8L, 8L))
  expect_equal(out$marker_filter, rep("markers=character", 2))
  expect_true(all(is.finite(out$selection_s_posthoc)))
})

test_that("summarise_stblr_maf_architecture supports trait-specific marker subsets", {
  fit <- make_selection_s_subset_fit()
  causal_by_trait <- list(
    D1 = paste0("m", 1:6),
    D2 = paste0("m", 7:14)
  )

  out <- summarise_stblr_maf_architecture(
    fit,
    maf = fit$maf,
    markers = causal_by_trait
  )

  expect_equal(out$n_markers, c(6L, 8L))
  expect_equal(out$marker_filter, rep("markers=list", 2))
  expect_true(all(is.finite(out$selection_s_posthoc)))
})

test_that("summarise_stblr_maf_architecture combines marker and min_pip filters", {
  fit <- make_selection_s_subset_fit()
  marker_subset <- paste0("m", 1:12)

  out <- summarise_stblr_maf_architecture(
    fit,
    maf = fit$maf,
    markers = marker_subset,
    min_pip = 0.05
  )

  expected_d1 <- sum(fit$dm[marker_subset, "D1"] >= 0.05)
  expected_d2 <- sum(fit$dm[marker_subset, "D2"] >= 0.05)
  expect_equal(out$n_markers, c(expected_d1, expected_d2))
  expect_equal(out$marker_filter, rep("markers=character;min_pip=0.05", 2))
})

test_that("summarise_stblr_maf_architecture returns NA statistics with fewer than five markers", {
  fit <- make_selection_s_subset_fit()

  out <- summarise_stblr_maf_architecture(
    fit,
    maf = fit$maf,
    markers = paste0("m", 1:4)
  )

  expect_equal(out$n_markers, c(4L, 4L))
  expect_true(all(is.na(out$selection_s_posthoc)))
  expect_true(all(is.na(out$se)))
  expect_true(all(is.na(out$p_value)))
  expect_true(all(is.na(out$r2)))
})

test_that("summarise_stblr_maf_architecture validates marker filter inputs", {
  fit <- make_selection_s_subset_fit()

  expect_error(
    summarise_stblr_maf_architecture(fit, maf = fit$maf, min_pip = -0.1),
    "min_pip"
  )
  expect_error(
    summarise_stblr_maf_architecture(fit, maf = fit$maf, min_pip = 1.1),
    "min_pip"
  )
  expect_error(
    summarise_stblr_maf_architecture(fit, maf = fit$maf, top_n = 0),
    "top_n"
  )
  expect_error(
    summarise_stblr_maf_architecture(fit, maf = fit$maf, top_n = 1.5),
    "top_n"
  )
})

test_that("summarise_stblr_maf_architecture handles missing MAF or h clearly", {
  fit <- list(
    dm = matrix(c(0.1, 0.5, 0.9), ncol = 1),
    bm = matrix(c(0.01, -0.02, 0.03), ncol = 1)
  )

  expect_error(
    summarise_stblr_maf_architecture(fit),
    "Either maf or h must be supplied or recoverable from fit"
  )
  expect_error(
    summarise_stblr_maf_architecture(fit, maf = c(0.1, 0.2)),
    "length equal to the number of markers"
  )
  expect_error(
    summarise_stblr_maf_architecture(fit, maf = c(0.1, 0.2, 0.7)),
    "maf must contain values"
  )
})
