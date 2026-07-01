source_sblr_test_file <- function(path) {
  candidates <- c(path, file.path("..", "..", path))
  path <- candidates[file.exists(candidates)][1L]
  if (is.na(path)) stop("Could not find ", path, call. = FALSE)
  source(path)
}

if (!exists("stblr_csr_annot", mode = "function") ||
    !"nchains" %in% names(formals(stblr_csr_annot))) {
  source_sblr_test_file("R/annotation-helpers.R")
  source_sblr_test_file("R/stblr-csr-annot.R")
}

make_tiny_sbayesrc_chain_csr_prefix <- function(m = 4L) {
  prefix <- tempfile("tiny_sbayesrc_chain_csr_")
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

tiny_sbayesrc_chain_stats <- function(m = 4L) {
  markers <- paste0("m", seq_len(m))
  list(
    wy = list(trait1 = stats::setNames(c(2, -1, 0.5, 1)[seq_len(m)], markers)),
    ww = list(trait1 = stats::setNames(rep(50, m), markers)),
    yy = stats::setNames(50, "trait1"),
    n = 50L,
    m = m,
    marker_names = markers,
    trait_names = "trait1"
  )
}

tiny_sbayesrc_chain_matrix <- function(m = 4L) {
  markers <- paste0("m", seq_len(m))
  matrix(
    c(
      1, 0,
      0, 1,
      1, 1,
      0, 0
    ),
    nrow = m,
    byrow = TRUE,
    dimnames = list(markers, c("coding", "qtl"))
  )
}

fit_tiny_sbayesrc_chains <- function(
    nchains = 1L, chain_seeds = NULL, keep_chains = FALSE) {
  stats <- tiny_sbayesrc_chain_stats()
  sblr::stblr_csr_annot(
    stats = stats,
    ld_prefix = make_tiny_sbayesrc_chain_csr_prefix(stats$m),
    annotations = tiny_sbayesrc_chain_matrix(),
    annotation_model = "sbayesrc",
    gamma = c(0, 0.1, 1),
    pi_init = 0.35,
    pi_prior_mean = 0.35,
    pi_prior_strength = 2,
    updateAlpha = FALSE,
    updateB = FALSE,
    updateE = FALSE,
    nit = 3,
    nburn = 0,
    ncores = 1L,
    seed = 101L,
    nchains = nchains,
    chain_seeds = chain_seeds,
    keep_chains = keep_chains
  )
}

expect_marker_summary_dims <- function(fit, stats) {
  expect_equal(dim(fit$dm), c(stats$m, length(stats$yy)))
  expect_equal(dim(fit$bm), c(stats$m, length(stats$yy)))
  for (nm in c("dm_sd", "dm_min", "dm_max", "bm_sd", "bm_min", "bm_max")) {
    expect_true(nm %in% names(fit), info = nm)
    expect_equal(dim(fit[[nm]]), dim(fit$dm), info = nm)
  }
}

test_that("single-chain SBayesRC annotation fit remains backward compatible", {
  skip_if_not(
    exists("stblr_cpg_omp_csr_sbayesrc", mode = "function"),
    "native SBayesRC CSR symbol is not loaded"
  )

  stats <- tiny_sbayesrc_chain_stats()
  fit <- fit_tiny_sbayesrc_chains(nchains = 1L)

  expect_equal(dim(fit$dm), c(stats$m, length(stats$yy)))
  expect_equal(dim(fit$bm), c(stats$m, length(stats$yy)))
  expect_equal(fit$input$nchains, 1L)
  expect_false(fit$input$keep_chains)
  expect_false("chains" %in% names(fit))
  expect_true("comp_prob" %in% names(fit))
  expect_equal(dim(fit$comp_prob$trait1), c(stats$m, 3L))
})

test_that("native SBayesRC annotation fit supports multiple chains", {
  skip_if_not(
    exists("stblr_cpg_omp_csr_sbayesrc", mode = "function"),
    "native SBayesRC CSR symbol is not loaded"
  )

  stats <- tiny_sbayesrc_chain_stats()
  fit <- fit_tiny_sbayesrc_chains(nchains = 2L, chain_seeds = c(11L, 12L))

  expect_equal(fit$input$nchains, 2L)
  expect_false(fit$input$keep_chains)
  expect_marker_summary_dims(fit, stats)
  expect_true(is.list(fit$comp_prob))
  expect_equal(dim(fit$comp_prob$trait1), c(stats$m, 3L))
  expect_equal(unname(rowSums(fit$comp_prob$trait1)), rep(1, stats$m), tolerance = 1e-8)
  expect_true(all(fit$comp_prob$trait1 >= 0 & fit$comp_prob$trait1 <= 1))
})

test_that("native SBayesRC annotation fit can keep compact chains", {
  skip_if_not(
    exists("stblr_cpg_omp_csr_sbayesrc", mode = "function"),
    "native SBayesRC CSR symbol is not loaded"
  )

  stats <- tiny_sbayesrc_chain_stats()
  fit <- fit_tiny_sbayesrc_chains(
    nchains = 2L,
    chain_seeds = c(21L, 22L),
    keep_chains = TRUE
  )

  expect_true(isTRUE(fit$input$keep_chains))
  expect_true(is.list(fit$chains))
  expect_length(fit$chains, length(stats$yy))
  expect_length(fit$chains$trait1, 2L)
  for (chain in fit$chains$trait1) {
    expect_length(chain$dm, stats$m)
    expect_length(chain$bm, stats$m)
    expect_equal(dim(chain$comp_prob), c(stats$m, 3L))
    expect_equal(dim(chain$alpha), dim(fit$alpha$trait1))
    expect_length(chain$sigmaSqAlpha, ncol(fit$alpha$trait1))
  }
})

test_that("SBayesRC annotation chain seeds are reproducible", {
  skip_if_not(
    exists("stblr_cpg_omp_csr_sbayesrc", mode = "function"),
    "native SBayesRC CSR symbol is not loaded"
  )

  fit1 <- fit_tiny_sbayesrc_chains(nchains = 2L, chain_seeds = c(31L, 32L))
  fit2 <- fit_tiny_sbayesrc_chains(nchains = 2L, chain_seeds = c(31L, 32L))

  expect_equal(fit1$dm, fit2$dm, tolerance = 1e-12)
  expect_equal(fit1$bm, fit2$bm, tolerance = 1e-12)
})

test_that("annotation chain argument validation is shared", {
  stats <- tiny_sbayesrc_chain_stats()
  A <- tiny_sbayesrc_chain_matrix()
  prefix <- make_tiny_sbayesrc_chain_csr_prefix(stats$m)

  expect_error(
    stblr_csr_annot(
      stats = stats, ld_prefix = prefix, annotations = A,
      annotation_model = "prior", nchains = 0L
    ),
    "nchains"
  )
  expect_error(
    stblr_csr_annot(
      stats = stats, ld_prefix = prefix, annotations = A,
      annotation_model = "learned", nchains = 2L, chain_seeds = 1L
    ),
    "chain_seeds"
  )
  expect_error(
    stblr_csr_annot(
      stats = stats, ld_prefix = prefix,
      annotations = c("a", "b", "a", "b"),
      annotation_model = "group", keep_chains = NA
    ),
    "keep_chains"
  )
})
