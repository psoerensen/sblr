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

tiny_annotation_bayesc_parallel_stats <- function(m = 4L) {
  stats <- tiny_annotation_bayesc_chain_stats(m)
  markers <- stats$marker_names
  stats$wy <- list(
    trait1 = stats::setNames(c(20, -10, 5, 10)[seq_len(m)], markers),
    trait2 = stats::setNames(c(-15, 7.5, 12.5, -5)[seq_len(m)], markers)
  )
  stats$ww <- list(
    trait1 = stats::setNames(rep(50, m), markers),
    trait2 = stats::setNames(rep(50, m), markers)
  )
  stats$yy <- stats::setNames(c(50, 50), c("trait1", "trait2"))
  stats$trait_names <- names(stats$yy)
  stats
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
  for (nm in c("dm_chain_mean_sd", "dm_chain_mean_min", "dm_chain_mean_max",
               "bm_chain_mean_sd", "bm_chain_mean_min", "bm_chain_mean_max")) {
    expect_true(nm %in% names(fit), info = nm)
    expect_equal(dim(fit[[nm]]), dim(fit$dm), info = nm)
  }
  if (keep_chains) {
    expect_true(is.list(fit$chains))
    expect_identical(names(fit$chains), paste0("task", seq_len(nchains)))
    expect_length(fit$chains, nchains)
    for (chain in fit$chains) {
      expect_identical(chain$trait_index, 1L)
      expect_identical(chain$trait_name, names(stats$yy))
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
    annotation_model = "fixed_marker",
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
    convergence = "none",
    chain_seeds = chain_seeds
  )
}

test_that("fixed-prior BayesC annotations support native chains", {
  skip_if_not(exists("stblr_cpg_omp_csr_prior", mode = "function"))
  stats <- tiny_annotation_bayesc_chain_stats()

  fit1 <- fit_prior_bayesc_chains(nchains = 1L, keep_chains = FALSE, chain_seeds = NULL)
  expect_equal(fit1$input$nchains, 1L)
  expect_true("chains" %in% names(fit1))
  expect_null(fit1$chains)

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
  fit <- sblr::stblr_csr_annot(
    stats = stats,
    ld_prefix = make_tiny_annotation_bayesc_chain_csr_prefix(stats$m),
    annotations = A,
    annotation_model = "learned_logistic",
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
    convergence = "none",
    chain_seeds = c(11L, 12L)
  )

  expect_bayesc_chain_fit(fit, stats)
  expect_true(all(c("eta_pi", "eta_vb") %in% names(fit$chains[[1]])))
})

test_that("two-trait learned-logistic fits are reproducible across one and two workers", {
  thread_info <- sblr:::sparseLD_thread_info(2L)
  skip_if_not(isTRUE(thread_info$openmp), "OpenMP is unavailable")
  skip_if_not(
    as.integer(thread_info$actual_threads_requested_region) >= 2L,
    "The OpenMP runtime cannot provide two workers"
  )
  stats <- tiny_annotation_bayesc_parallel_stats()
  A <- tiny_annotation_bayesc_chain_matrix()
  prefix <- make_tiny_annotation_bayesc_chain_csr_prefix(stats$m)
  on.exit(unlink(paste0(prefix, c(
    ".row_ptr.u64.bin", ".col_idx.u32.0based.bin",
    ".values.f32.bin", ".meta.txt"
  ))), add = TRUE)
  common <- list(
    stats = stats,
    ld_prefix = prefix,
    annotations = A,
    annotation_model = "learned_logistic",
    pi_init = 0.35,
    pi_prior_a = 0.7,
    pi_prior_b = 1.3,
    learn_pi_annot = TRUE,
    learn_vb_annot = FALSE,
    eta_pi_init = matrix(c(1, -0.8, 0.6, -0.4),
                         nrow = ncol(A), ncol = length(stats$yy)),
    rw_sd_eta_pi = 0,
    annot_update_every = 1L,
    d_init = list(c(1, 0, 0, 1), c(0, 1, 1, 0)),
    use_d_init = TRUE,
    updateB = FALSE,
    updateE = FALSE,
    updatePi = TRUE,
    nit = 12L,
    nburn = 3L,
    nchains = 2L,
    keep_chains = TRUE,
    convergence = "none",
    chain_seeds = c(7311L, 7312L),
    seed = 7310L
  )
  run <- function(ncores) {
    do.call(sblr::stblr_csr_annot, c(common, list(ncores = ncores)))
  }
  centered_offsets <- sweep(A %*% common$eta_pi_init, 2L,
                            colMeans(A %*% common$eta_pi_init), "-")
  expect_true(all(apply(centered_offsets, 2L, function(x) any(x != 0))))
  serial <- expect_no_warning(run(1L))
  serial_repeat <- expect_no_warning(run(1L))
  parallel <- expect_no_warning(run(2L))
  serial_workers <- serial$diagnostics$native$parallel
  parallel_workers <- parallel$diagnostics$native$parallel
  expect_identical(
    parallel_workers$scope,
    "learned_sampler_trait_parallel_region"
  )
  expect_true(isTRUE(parallel_workers$openmp))
  expect_equal(parallel_workers$requested_thread_count, rep(2L, 2L))
  expect_equal(parallel_workers$configured_thread_count, rep(2L, 2L))
  expect_true(all(parallel_workers$actual_team_size >= 2L))
  expect_true(
    length(parallel_workers$runtime_max_threads_before_request) == 1L &&
      is.finite(parallel_workers$runtime_max_threads_before_request) &&
      parallel_workers$runtime_max_threads_before_request >= 1L
  )
  expect_equal(dim(parallel_workers$trait_worker_id), c(2L, 2L))
  expect_true(all(apply(
    parallel_workers$trait_worker_id,
    2L,
    function(worker) length(unique(worker)) >= 2L
  )))
  expect_equal(serial_workers$requested_thread_count, rep(1L, 2L))
  expect_equal(serial_workers$configured_thread_count, rep(1L, 2L))
  expect_equal(serial_workers$actual_team_size, rep(1L, 2L))
  expect_equal(serial_workers$trait_worker_id, matrix(0L, 2L, 2L))
  expect_true(all(diag(parallel$cov_g_mean) > 0))
  for (field in c(
    "bm", "dm", "wy", "r", "b", "d", "b_final", "d_final",
    "pi_trace", "pi_final",
    "pi_mean", "eta_pi", "eta_vb", "vbs", "vgs", "ves", "vle", "vld",
    "chains", "convergence_traces"
  )) {
    expect_identical(serial[[field]], serial_repeat[[field]], info = field)
    expect_identical(serial[[field]], parallel[[field]], info = field)
  }
})

test_that("learned-logistic updatePi FALSE preserves the global probability", {
  stats <- tiny_annotation_bayesc_chain_stats()
  A <- tiny_annotation_bayesc_chain_matrix()
  prefix <- make_tiny_annotation_bayesc_chain_csr_prefix(stats$m)
  on.exit(unlink(paste0(prefix, c(
    ".row_ptr.u64.bin", ".col_idx.u32.0based.bin",
    ".values.f32.bin", ".meta.txt"
  ))), add = TRUE)
  fit <- sblr::stblr_csr_annot(
    stats = stats,
    ld_prefix = prefix,
    annotations = A,
    annotation_model = "learned_logistic",
    pi_init = 0.35,
    pi_prior_mean = 0.35,
    pi_prior_strength = 2,
    learn_pi_annot = TRUE,
    eta_pi_init = matrix(c(1, -0.8), nrow = ncol(A), ncol = 1L),
    rw_sd_eta_pi = 0,
    annot_update_every = 1L,
    updateB = FALSE,
    updateE = FALSE,
    updatePi = FALSE,
    nit = 5L,
    nburn = 0L,
    ncores = 1L,
    seed = 7313L,
    convergence = "none"
  )
  expect_equal(as.numeric(fit$pi_trace), rep(0.35, 5L), tolerance = 0)
  expect_equal(unname(fit$pi_final), c(0.65, 0.35), tolerance = 0)
})

test_that("group BayesC annotations support native chains", {
  skip_if_not(exists("stblr_cpg_omp_csr_group_annot", mode = "function"))
  stats <- tiny_annotation_bayesc_chain_stats()
  group <- stats::setNames(c("coding", "background", "coding", "background"),
                           stats$marker_names)
  fit <- sblr::stblr_csr_annot(
    stats = stats,
    ld_prefix = make_tiny_annotation_bayesc_chain_csr_prefix(stats$m),
    annotations = group,
    annotation_model = "group",
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
    convergence = "none",
    chain_seeds = c(21L, 22L)
  )

  expect_bayesc_chain_fit(fit, stats)
  expect_true(all(c(
    "group_pi", "group_vb_multiplier", "group_nincluded"
  ) %in% names(fit$chains[[1]])))
})
