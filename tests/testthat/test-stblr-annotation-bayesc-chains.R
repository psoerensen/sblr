make_tiny_annotation_bayesc_chain_csr_prefix <- function(m = 4L) {
  prefix <- tempfile("tiny_annotation_bayesc_chain_csr_")
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

tiny_annotation_bayesc_chain_stats <- function(m = 4L) {
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

tiny_annotation_bayesc_chain_matrix <- function(m = 4L) {
  markers <- paste0("m", seq_len(m))
  matrix(
    c(1, 0, 0, 1, 1, 1, 0, 0),
    nrow = m,
    byrow = TRUE,
    dimnames = list(markers, c("coding", "qtl"))
  )
}

expect_bayesc_chain_fit <- function(fit, stats, nchains = 2L, keep_chains = TRUE) {
  expect_equal(fit$input$nchains, nchains)
  expect_identical(fit$input$keep_chains, keep_chains)
  expect_equal(dim(fit$dm), c(stats$m, length(stats$yy)))
  expect_equal(dim(fit$bm), c(stats$m, length(stats$yy)))
  for (nm in c("dm_sd", "dm_min", "dm_max", "bm_sd", "bm_min", "bm_max")) {
    expect_true(nm %in% names(fit), info = nm)
    expect_equal(dim(fit[[nm]]), dim(fit$dm), info = nm)
  }
  if (keep_chains) {
    expect_true(is.list(fit$chains))
    expect_identical(names(fit$chains), names(stats$yy))
    expect_length(fit$chains[[1]], nchains)
    for (chain in fit$chains[[1]]) {
      expect_length(chain$dm, stats$m)
      expect_length(chain$bm, stats$m)
      expect_identical(names(chain$dm), stats$marker_names)
      expect_identical(names(chain$bm), stats$marker_names)
    }
  }
}

fit_prior_bayesc_chains <- function(nchains = 2L, keep_chains = TRUE,
                                    chain_seeds = c(10L, 20L)) {
  stats <- tiny_annotation_bayesc_chain_stats()
  A <- tiny_annotation_bayesc_chain_matrix()
  sblr::stblr_csr_annot(
    stats = stats,
    ld_prefix = make_tiny_annotation_bayesc_chain_csr_prefix(stats$m),
    annotations = list(
      A = A,
      fixed_pi_marker = list(rep(0.35, stats$m)),
      fixed_vb_multiplier = list(c(1, 1.2, 0.8, 1))
    ),
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
    seed = 100L,
    nchains = nchains,
    keep_chains = keep_chains,
    chain_seeds = chain_seeds
  )
}

test_that("fixed-prior BayesC annotations support native chains", {
  skip_if_not(exists("stblr_cpg_omp_csr_prior", mode = "function"))
  stats <- tiny_annotation_bayesc_chain_stats()

  fit1 <- fit_prior_bayesc_chains(nchains = 1L, keep_chains = FALSE, chain_seeds = NULL)
  expect_equal(fit1$input$nchains, 1L)
  expect_false("chains" %in% names(fit1))

  fit <- fit_prior_bayesc_chains()
  expect_bayesc_chain_fit(fit, stats)

  fit_b <- fit_prior_bayesc_chains()
  expect_equal(fit$dm, fit_b$dm, tolerance = 1e-12)
  expect_equal(fit$bm, fit_b$bm, tolerance = 1e-12)
})

test_that("learned BayesC annotations support native chains", {
  skip_if_not(exists("stblr_cpg_omp_csr_annot", mode = "function"))
  stats <- tiny_annotation_bayesc_chain_stats()
  A <- tiny_annotation_bayesc_chain_matrix()
  fit <- sblr::stblr_csr_learn_annot(
    stats = stats,
    ld_prefix = make_tiny_annotation_bayesc_chain_csr_prefix(stats$m),
    A = A,
    pi_init = 0.35,
    pi_prior_mean = 0.35,
    pi_prior_strength = 2,
    learn_pi_annot = TRUE,
    learn_vb_annot = TRUE,
    updateB = FALSE,
    updateE = FALSE,
    updatePi = FALSE,
    annot_update_every = 1L,
    rw_sd_eta_pi = 0.01,
    rw_sd_eta_vb = 0.01,
    nit = 3,
    nburn = 0,
    ncores = 1L,
    seed = 101L,
    nchains = 2L,
    keep_chains = TRUE,
    chain_seeds = c(11L, 12L)
  )

  expect_bayesc_chain_fit(fit, stats)
  expect_true(all(c("eta_pi", "eta_vb") %in% names(fit$chains[[1]][[1]])))
})

test_that("group BayesC annotations support native chains", {
  skip_if_not(exists("stblr_cpg_omp_csr_group_annot", mode = "function"))
  stats <- tiny_annotation_bayesc_chain_stats()
  group <- stats::setNames(c("coding", "background", "coding", "background"),
                           stats$marker_names)
  fit <- sblr::stblr_csr_group_annot(
    stats = stats,
    ld_prefix = make_tiny_annotation_bayesc_chain_csr_prefix(stats$m),
    group = group,
    group_names = c("coding", "background"),
    group_pi_init = c(0.35, 0.25),
    group_vb_multiplier_init = c(1.2, 0.8),
    pi_init = 0.3,
    pi_prior_mean = 0.3,
    pi_prior_strength = 2,
    updateB = FALSE,
    updateE = FALSE,
    updatePi = FALSE,
    updateGroupVb = FALSE,
    nit = 3,
    nburn = 0,
    ncores = 1L,
    seed = 102L,
    nchains = 2L,
    keep_chains = TRUE,
    chain_seeds = c(21L, 22L)
  )

  expect_bayesc_chain_fit(fit, stats)
  expect_true(all(c(
    "group_pi", "group_vb_multiplier", "group_nincluded"
  ) %in% names(fit$chains[[1]][[1]])))
})
