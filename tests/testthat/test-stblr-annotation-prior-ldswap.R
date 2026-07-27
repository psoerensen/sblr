make_tiny_annotation_prior_ldswap_csr_prefix <- function() {
  prefix <- tempfile("tiny_annotation_prior_ldswap_csr_")
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

tiny_annotation_prior_ldswap_stats <- function() {
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

tiny_annotation_prior_ldswap_A <- function() {
  matrix(
    c(1, 0, 1, 1, 0, 1),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(paste0("m", 1:3), c("coding", "qtl"))
  )
}

tiny_annotation_prior_ldswap_annotations <- function() {
  list(
    A = tiny_annotation_prior_ldswap_A(),
    fixed_pi_marker = list(c(0.55, 0.25, 0.15)),
    fixed_vb_multiplier = list(c(1.0, 1.4, 0.8)),
    use_pi_marker = TRUE,
    use_vb_multiplier = TRUE
  )
}

fit_tiny_annotation_prior_ldswap <- function(updateLDswap = TRUE,
                                             nchains = 1L,
                                             keep_chains = FALSE,
                                             chain_seeds = NULL,
                                             seed = 31L) {
  stblr_csr_annot(
    stats = tiny_annotation_prior_ldswap_stats(),
    ld_prefix = make_tiny_annotation_prior_ldswap_csr_prefix(),
    annotations = tiny_annotation_prior_ldswap_annotations(),
    annotation_model = "prior",
    pi_init = 0.35,
    pi_prior_mean = 0.35,
    pi_prior_strength = 2,
    updateB = FALSE,
    updateE = FALSE,
    updatePi = FALSE,
    nit = 3,
    nburn = 0,
    ncores = 1L,
    seed = seed,
    nchains = nchains,
    keep_chains = keep_chains,
    chain_seeds = chain_seeds,
    updateLDswap = updateLDswap,
    ld_swap_prob = 1,
    ld_swap_r2 = 0.8,
    ld_swap_max_friends = 2L,
    ld_swap_moves = 2L
  )
}

expect_annotation_prior_ldswap_diagnostics <- function(fit, expected_attempted) {
  ld_swap <- fit$diagnostics$ld_swap
  if (is.null(expected_attempted)) {
    expect_null(ld_swap)
    return(invisible())
  }
  expect_identical(colnames(ld_swap),
                   c("attempted", "accepted", "acceptance_rate"))
  expect_equal(ld_swap$attempted, expected_attempted)
  expect_true(ld_swap$accepted >= 0)
  expect_true(ld_swap$acceptance_rate >= 0)
  expect_true(ld_swap$acceptance_rate <= 1)
}

test_that("fixed-prior annotation CSR BayesC is backward compatible without LD-swap", {
  fit <- fit_tiny_annotation_prior_ldswap(updateLDswap = FALSE)

  expect_equal(fit$input$backend, "csr_prior_bayesc")
  expect_false(isTRUE(fit$input$updateLDswap))
  expect_equal(dim(fit$dm), c(3L, 1L))
  expect_equal(dim(fit$bm), c(3L, 1L))
  expect_annotation_prior_ldswap_diagnostics(fit, expected_attempted = NULL)
})

test_that("fixed-prior annotation CSR BayesC supports LD-swap diagnostics", {
  fit <- fit_tiny_annotation_prior_ldswap(updateLDswap = TRUE)

  expect_true(isTRUE(fit$input$updateLDswap))
  expect_equal(fit$input$backend, "csr_prior_bayesc")
  expect_equal(dim(fit$dm), c(3L, 1L))
  expect_equal(dim(fit$bm), c(3L, 1L))
  expect_annotation_prior_ldswap_diagnostics(fit, expected_attempted = 6)
})

test_that("fixed-prior annotation CSR BayesC supports LD-swap with kept chains", {
  fit <- fit_tiny_annotation_prior_ldswap(
    updateLDswap = TRUE,
    nchains = 2L,
    keep_chains = TRUE,
    chain_seeds = c(41L, 42L)
  )

  expect_equal(fit$input$nchains, 2L)
  expect_true(isTRUE(fit$input$keep_chains))
  expect_annotation_prior_ldswap_diagnostics(fit, expected_attempted = 12)
  expect_equal(sum(fit$diagnostics$ld_swap_chains$trait1$attempted),
               fit$diagnostics$ld_swap$attempted)
  expect_true("chains" %in% names(fit))
  expect_length(fit$chains, 2L)
  expect_true(all(vapply(fit$chains, function(ch) {
    all(c("dm", "bm", "ld_swap") %in% names(ch)) &&
      length(ch$dm) == 3L &&
      length(ch$bm) == 3L &&
      identical(names(ch$dm), tiny_annotation_prior_ldswap_stats()$marker_names) &&
      identical(names(ch$bm), tiny_annotation_prior_ldswap_stats()$marker_names)
  }, logical(1))))
})

test_that("fixed-prior annotation CSR BayesC LD-swap is reproducible with chain seeds", {
  args <- list(
    updateLDswap = TRUE,
    nchains = 2L,
    keep_chains = TRUE,
    chain_seeds = c(51L, 52L)
  )
  fit_a <- do.call(fit_tiny_annotation_prior_ldswap, args)
  fit_b <- do.call(fit_tiny_annotation_prior_ldswap, args)

  expect_equal(fit_a$dm, fit_b$dm, tolerance = 1e-12)
  expect_equal(fit_a$bm, fit_b$bm, tolerance = 1e-12)
  expect_equal(fit_a$diagnostics$ld_swap, fit_b$diagnostics$ld_swap)
})

test_that("old fixed-prior annotation wrapper is not exported", {
  expect_false("stblr_csr_prior_annot" %in% getNamespaceExports("sblr"))
})
