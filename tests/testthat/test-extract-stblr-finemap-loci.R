make_extract_finemap_fit <- function() {
  markers <- paste0("m", 1:4)
  list(
    dm = matrix(
      c(0.10, 0.60, 0.35, 0.02),
      ncol = 1,
      dimnames = list(markers, "D1")
    ),
    bm = matrix(
      c(0.01, 0.20, -0.10, 0.00),
      ncol = 1,
      dimnames = list(markers, "D1")
    ),
    dm_sd = matrix(
      c(0.01, 0.06, 0.03, 0.002),
      ncol = 1,
      dimnames = list(markers, "D1")
    ),
    dm_min = matrix(
      c(0.08, 0.50, 0.30, 0.01),
      ncol = 1,
      dimnames = list(markers, "D1")
    ),
    dm_max = matrix(
      c(0.12, 0.70, 0.40, 0.03),
      ncol = 1,
      dimnames = list(markers, "D1")
    ),
    bm_sd = matrix(
      c(0.001, 0.02, 0.01, 0.000),
      ncol = 1,
      dimnames = list(markers, "D1")
    ),
    bm_min = matrix(
      c(0.00, 0.15, -0.12, 0.00),
      ncol = 1,
      dimnames = list(markers, "D1")
    ),
    bm_max = matrix(
      c(0.02, 0.25, -0.08, 0.00),
      ncol = 1,
      dimnames = list(markers, "D1")
    )
  )
}

make_extract_finemap_glist <- function(prefix = NULL) {
  list(
    rsids = list(paste0("m", 1:4)),
    chr = list(rep(1L, 4)),
    pos = list(c(100, 200, 300, 400)),
    sparseLD = list(
      chr = 1L,
      cls = list(1:4),
      prefix = prefix
    )
  )
}

make_extract_finemap_csr_prefix <- function() {
  prefix <- tempfile("extract_finemap_csr_")
  row_ptr <- c(0, 2, 3, 3, 3)
  col_idx <- c(1L, 2L, 2L)
  values <- c(sqrt(0.8), sqrt(0.6), sqrt(0.7))

  sblr:::.stblr_write_uint64_file(paste0(prefix, ".row_ptr.u64.bin"), row_ptr)
  sblr:::.stblr_write_uint32_file(paste0(prefix, ".col_idx.u32.0based.bin"), col_idx)
  writeBin(as.numeric(values), paste0(prefix, ".values.f32.bin"),
           size = 4, endian = "little")

  writeLines(
    c(
      "format=sparse_ld_csr",
      "storage=streamed_upper_triangle",
      "n_bed=4",
      "n_used=4",
      "n_samples_used=100",
      "n_variants=4",
      "nnz=3",
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

test_that("extract_stblr_finemap_loci summarizes fitted posterior without sparse LD", {
  fit <- make_extract_finemap_fit()
  Glist <- make_extract_finemap_glist()
  locus_sets <- list(regionA = c("m1", "m2", "m3"))

  out <- extract_stblr_finemap_loci(
    fit = fit,
    Glist = Glist,
    locus_sets = locus_sets,
    trait = "D1",
    credible_sets = FALSE
  )

  expect_s3_class(out, "stblr_finemap")
  expect_null(out$credible_sets)
  expect_false("runs" %in% names(out))

  expect_named(
    out$markers,
    c(
      "locus", "trait", "marker", "chr", "pos", "pip_mean", "pip_sd",
      "pip_min", "pip_max", "bm_mean", "bm_sd", "bm_min", "bm_max"
    )
  )
  expect_named(
    out$summary,
    c(
      "locus", "trait", "chr", "start", "end", "n_markers",
      "lead_marker", "lead_pip", "lead_pip_sd", "total_pip",
      "secondary_pip"
    )
  )

  expect_equal(out$summary$lead_marker, "m2")
  expect_equal(out$summary$lead_pip, 0.60)
  expect_equal(out$summary$lead_pip_sd, 0.06)
  expect_equal(out$summary$total_pip, 1.05)
  expect_equal(out$summary$secondary_pip, 0.45)
  expect_equal(out$markers$pip_min, c(0.08, 0.50, 0.30))
  expect_equal(out$markers$pip_max, c(0.12, 0.70, 0.40))
  expect_equal(out$markers$bm_sd, c(0.001, 0.02, 0.01))
  expect_equal(out$markers$bm_min, c(0.00, 0.15, -0.12))
  expect_equal(out$markers$bm_max, c(0.02, 0.25, -0.08))
})

test_that("extract_stblr_finemap_loci builds credible sets from sparse LD", {
  fit <- make_extract_finemap_fit()
  prefix <- make_extract_finemap_csr_prefix()
  Glist <- make_extract_finemap_glist(prefix)
  locus_sets <- list(regionA = c("m1", "m2", "m3"))

  out <- extract_stblr_finemap_loci(
    fit = fit,
    Glist = Glist,
    locus_sets = locus_sets,
    trait = "D1",
    credible_sets = TRUE,
    coverage = 0.70,
    min_r2 = 0.5,
    cs_mode = "single"
  )

  expect_s3_class(out, "stblr_finemap")
  expect_s3_class(out$credible_sets$summary, "data.frame")
  expect_gt(nrow(out$credible_sets$summary), 0)
  expect_named(out$credible_sets$sets, "regionA")
  expect_named(out$credible_sets$sets$regionA, "D1")
})
