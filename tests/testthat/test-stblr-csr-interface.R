make_stblr_csr_interface_prefix <- function(m = 3L) {
  prefix <- tempfile("stblr_csr_interface_ld_")
  row_ptr <- rep(0, m + 1L)

  sblr:::.stblr_write_uint64_file(paste0(prefix, ".row_ptr.u64.bin"), row_ptr)
  file.create(paste0(prefix, ".col_idx.u32.0based.bin"))
  file.create(paste0(prefix, ".values.f32.bin"))

  writeLines(
    c(
      "format=sparse_ld_csr",
      "storage=streamed_upper_triangle",
      "n_bed=NA",
      "n_used=NA",
      "n_samples_used=NA",
      paste0("n_variants=", m),
      "nnz=0",
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

make_stblr_csr_interface_stats <- function() {
  markers <- paste0("m", 1:3)
  list(
    wy = list(trait1 = stats::setNames(c(2, -1, 0.5), markers)),
    ww = list(trait1 = stats::setNames(rep(50, 3), markers)),
    yy = stats::setNames(50, "trait1"),
    n = 50L,
    m = 3L,
    marker_names = markers,
    trait_names = "trait1"
  )
}

test_that("stblr_csr fits BayesC through public method interface", {
  skip_if_not(
    exists("stblr_cpg_omp_csr", mode = "function"),
    "native BayesC CSR symbol is not loaded"
  )

  fit <- stblr_csr(
    stats = make_stblr_csr_interface_stats(),
    ld_prefix = make_stblr_csr_interface_prefix(),
    method = "bayesC",
    pi_init = 0.5,
    pi_prior_mean = 0.5,
    pi_prior_strength = 2,
    updateB = FALSE,
    updateE = FALSE,
    updatePi = FALSE,
    nit = 2,
    nburn = 0,
    nchains = 1L,
    ncores = 1L,
    seed = 101L
  )

  expect_true(all(c("dm", "bm", "input") %in% names(fit)))
  expect_equal(fit$input$method, "bayesc")
  expect_equal(fit$input$model, "bayesc")
  expect_equal(fit$input$backend, "csr_bayesc")
  expect_equal(fit$input$data_level, "summary")
  expect_equal(fit$input$scheduled, FALSE)
  expect_equal(fit$input$nchains, 1L)
})

test_that("stblr_csr fits BayesR through public method interface", {
  skip_if_not(
    exists("stblr_cpg_omp_csr_bayesr", mode = "function"),
    "native BayesR CSR symbol is not loaded"
  )

  fit <- stblr_csr(
    stats = make_stblr_csr_interface_stats(),
    ld_prefix = make_stblr_csr_interface_prefix(),
    method = "bayesR",
    updateB = FALSE,
    updateE = FALSE,
    updatePi = FALSE,
    nit = 2,
    nburn = 0,
    nchains = 1L,
    ncores = 1L,
    seed = 102L
  )

  expect_true(all(c("dm", "bm", "comp_prob", "dm_component_mean", "input") %in% names(fit)))
  expect_equal(fit$input$method, "bayesr")
  expect_equal(fit$input$model, "bayesr")
  expect_equal(fit$input$backend, "csr_bayesr")
  expect_equal(fit$input$data_level, "summary")
  expect_equal(fit$input$scheduled, FALSE)
  expect_equal(fit$input$nchains, 1L)
  expect_equal(
    unname(as.numeric(fit$dm[, "trait1"])),
    unname(1 - fit$comp_prob$trait1[, "component_0"]),
    tolerance = 1e-8
  )
})
