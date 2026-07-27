stblr_csr_annot <- getFromNamespace("stblr_csr_annot", "sblr")

make_tiny_annotation_sbayesrc_ldswap_csr_prefix <- function() {
  prefix <- tempfile("tiny_annotation_sbayesrc_ldswap_csr_")
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

tiny_annotation_sbayesrc_ldswap_stats <- function() {
  markers <- paste0("m", 1:3)
  list(
    wy = list(trait1 = stats::setNames(c(8, 0.1, 0.2), markers)),
    ww = list(trait1 = stats::setNames(rep(100, 3), markers)),
    yy = stats::setNames(100, "trait1"),
    n = 100L,
    m = 3L,
    marker_names = markers,
    trait_names = "trait1"
  )
}

tiny_annotation_sbayesrc_ldswap_A <- function() {
  matrix(
    c(1, 0, 1, 1, 0, 1),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(paste0("m", 1:3), c("coding", "qtl"))
  )
}

fit_tiny_annotation_sbayesrc_ldswap <- function(updateLDswap = TRUE,
                                                nchains = 1L,
                                                keep_chains = FALSE,
                                                chain_seeds = NULL,
                                                seed = 71L) {
  stats <- tiny_annotation_sbayesrc_ldswap_stats()
  stblr_csr_annot(
    stats = stats,
    Glist = list(
      rsidsLD = list(stats$marker_names),
      rsids = list(stats$marker_names),
      maf = list(rep(0.2, stats$m))
    ),
    ld_prefix = make_tiny_annotation_sbayesrc_ldswap_csr_prefix(),
    annotations = tiny_annotation_sbayesrc_ldswap_A(),
    annotation_model = "sbayesrc",
    gamma = c(0, 0.1, 1),
    pi_init = 0.45,
    pi_prior_mean = 0.45,
    pi_prior_strength = 2,
    updateAlpha = FALSE,
    updateB = FALSE,
    updateE = FALSE,
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

expect_sbayesrc_ldswap_diagnostics <- function(fit, available = TRUE) {
  if (!available) {
    expect_null(fit$diagnostics$ld_swap)
    return(invisible())
  }
  expect_identical(colnames(fit$diagnostics$ld_swap),
                   c("attempted", "accepted", "acceptance_rate"))
  expect_true(all(fit$diagnostics$ld_swap$attempted >=
                  fit$diagnostics$ld_swap$accepted))
  expect_true(all(fit$diagnostics$ld_swap$accepted >= 0))
  expect_true(all(fit$diagnostics$ld_swap$acceptance_rate >= 0))
  expect_true(all(fit$diagnostics$ld_swap$acceptance_rate <= 1))
}

test_that("SBayesRC annotation fit is backward compatible without LD-swap", {
  skip_if_not(
    exists("stblr_cpg_omp_csr_sbayesrc", mode = "function"),
    "native SBayesRC CSR symbol is not loaded"
  )

  fit <- fit_tiny_annotation_sbayesrc_ldswap(updateLDswap = FALSE)

  expect_equal(fit$input$backend, "csr_sbayesrc")
  expect_false(isTRUE(fit$input$updateLDswap))
  expect_equal(dim(fit$dm), c(3L, 1L))
  expect_equal(dim(fit$bm), c(3L, 1L))
  expect_sbayesrc_ldswap_diagnostics(fit, available = FALSE)
  expect_true("component_probabilities" %in% names(fit))
})

test_that("SBayesRC annotation fit supports LD-swap diagnostics", {
  skip_if_not(
    exists("stblr_cpg_omp_csr_sbayesrc", mode = "function"),
    "native SBayesRC CSR symbol is not loaded"
  )

  fit <- fit_tiny_annotation_sbayesrc_ldswap(updateLDswap = TRUE)

  expect_true(isTRUE(fit$input$updateLDswap))
  expect_equal(fit$input$backend, "csr_sbayesrc")
  expect_equal(fit$input$ld_swap_moves, 2L)
  expect_equal(dim(fit$dm), c(3L, 1L))
  expect_equal(dim(fit$bm), c(3L, 1L))
  expect_sbayesrc_ldswap_diagnostics(fit)
  expect_true("alpha" %in% names(fit))
  expect_true("sigmaSqAlpha" %in% names(fit))
  expect_true("ncomp" %in% names(fit))
  expect_equal(dim(fit$component_probabilities$trait1), c(3L, 3L))
  expect_equal(unname(rowSums(fit$component_probabilities$trait1)), rep(1, 3), tolerance = 1e-8)
  expect_true(all(is.finite(fit$component_probabilities$trait1)))
  expect_true(all(fit$component_probabilities$trait1 >= 0 & fit$component_probabilities$trait1 <= 1))
})

test_that("SBayesRC annotation LD-swap works with kept chains", {
  skip_if_not(
    exists("stblr_cpg_omp_csr_sbayesrc", mode = "function"),
    "native SBayesRC CSR symbol is not loaded"
  )

  fit <- fit_tiny_annotation_sbayesrc_ldswap(
    updateLDswap = TRUE,
    nchains = 2L,
    keep_chains = TRUE,
    chain_seeds = c(81L, 82L)
  )

  expect_equal(fit$input$nchains, 2L)
  expect_true(isTRUE(fit$input$keep_chains))
  expect_sbayesrc_ldswap_diagnostics(fit)
  expect_true("ld_swap_chains" %in% names(fit$diagnostics))
  expect_equal(sum(fit$diagnostics$ld_swap_chains$trait1$attempted),
               fit$diagnostics$ld_swap$attempted)
  expect_length(fit$chains, 2L)
  expect_true(all(vapply(fit$chains, function(ch) {
    all(c("dm", "bm", "ld_swap", "component_probabilities",
          "alpha", "sigmaSqAlpha") %in% names(ch)) &&
      length(ch$dm) == 3L &&
      length(ch$bm) == 3L &&
      identical(names(ch$dm), tiny_annotation_sbayesrc_ldswap_stats()$marker_names) &&
      identical(names(ch$bm), tiny_annotation_sbayesrc_ldswap_stats()$marker_names) &&
      identical(dim(ch$component_probabilities), c(3L, 3L))
  }, logical(1))))
})

test_that("SBayesRC annotation LD-swap is reproducible with chain seeds", {
  skip_if_not(
    exists("stblr_cpg_omp_csr_sbayesrc", mode = "function"),
    "native SBayesRC CSR symbol is not loaded"
  )

  args <- list(
    updateLDswap = TRUE,
    nchains = 2L,
    keep_chains = TRUE,
    chain_seeds = c(91L, 92L)
  )
  fit_a <- do.call(fit_tiny_annotation_sbayesrc_ldswap, args)
  fit_b <- do.call(fit_tiny_annotation_sbayesrc_ldswap, args)

  expect_equal(fit_a$dm, fit_b$dm, tolerance = 1e-12)
  expect_equal(fit_a$bm, fit_b$bm, tolerance = 1e-12)
  expect_equal(fit_a$diagnostics$ld_swap, fit_b$diagnostics$ld_swap)
})

test_that("old SBayesRC annotation wrapper is not exported", {
  expect_false("stblr_csr_sbayesrc_generic" %in% getNamespaceExports("sblr"))
})
