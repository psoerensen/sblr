make_tiny_annotation_csr_prefix <- function(m = 4L) {
  prefix <- tempfile("tiny_annotation_csr_")
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

tiny_annotation_stats <- function(m = 4L) {
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

tiny_annotation_matrix <- function(m = 4L) {
  markers <- paste0("m", seq_len(m))
  A <- matrix(
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
  A
}

expect_annotation_fit_core <- function(fit, stats) {
  expect_true(all(c("dm", "bm", "input") %in% names(fit)))

  dm <- as.matrix(fit$dm)
  bm <- as.matrix(fit$bm)

  expect_equal(dim(dm), dim(bm))
  expect_equal(nrow(dm), stats$m)
  expect_equal(ncol(dm), length(stats$yy))
  expect_identical(rownames(dm), stats$marker_names)
  expect_identical(rownames(bm), stats$marker_names)
  expect_identical(colnames(dm), names(stats$yy))
  expect_identical(colnames(bm), names(stats$yy))
  expect_true(all(is.na(dm) | (dm >= -1e-8 & dm <= 1 + 1e-8)))
  expect_true(all(is.na(bm) | is.finite(bm)))

  expect_equal(length(colMeans(dm, na.rm = TRUE)), ncol(dm))
  expect_equal(length(apply(dm, 2, max, na.rm = TRUE)), ncol(dm))
}

expect_trace_matrix <- function(x, ntraits) {
  expect_true(is.matrix(x))
  expect_equal(ncol(x), ntraits)
  expect_true(all(is.na(x) | is.finite(x)))
}

expect_square_trait_matrix <- function(x, ntraits) {
  expect_true(is.matrix(x))
  expect_equal(dim(x), c(ntraits, ntraits))
  expect_true(all(is.na(x) | is.finite(x)))
}

test_that("annotation-aware exported wrappers and helpers are available", {
  exports <- getNamespaceExports("sblr")
  for (fn in c(
    "stblr_csr_prior_annot",
    "stblr_csr_learn_annot",
    "stblr_csr_group_annot",
    "stblr_csr_sbayesrc_generic",
    "make_sbayesrc_alpha_init",
    "sbayesrc_annotation_pi",
    "sbayesrc_annotation_gamma_mean",
    "sbayesrc_marker_pi",
    "sbayesrc_marker_gamma_mean",
    "mtsim_annotation",
    "summarize_annotation_signal"
  )) {
    expect_true(fn %in% exports, info = fn)
    expect_true(is.function(getExportedValue("sblr", fn)), info = fn)
  }
})

test_that("annotation backend design document records the backend contract", {
  path <- test_path("..", "..", "docs", "dev", "stblr_annotation_backend_design.md")
  expect_true(file.exists(path))
  txt <- readLines(path, warn = FALSE)
  for (backend in c(
    "csr_prior_bayesc",
    "csr_annot_bayesc",
    "csr_group_bayesc",
    "csr_sbayesrc"
  )) {
    expect_true(any(grepl(backend, txt, fixed = TRUE)), info = backend)
  }
})

test_that("SBayesRC diagnostic helpers preserve annotation and marker shapes", {
  gamma <- c(0, 0.01, 0.1, 1)
  A <- cbind(
    intercept = rep(1, 4),
    coding = c(1, 0, 1, 0),
    qtl = c(0, 1, 1, 0)
  )
  rownames(A) <- paste0("m", seq_len(nrow(A)))

  init <- sblr::make_sbayesrc_alpha_init(
    A = A,
    gamma = gamma,
    pi_init = 0.2
  )
  expect_equal(dim(init$alpha_init), c(ncol(A), length(gamma) - 1L))
  expect_true(all(is.finite(init$alpha_init)))
  expect_identical(rownames(init$alpha_init), colnames(A))

  ann_pi <- sblr::sbayesrc_annotation_pi(init$alpha_init, gamma)
  ann_gamma <- sblr::sbayesrc_annotation_gamma_mean(init$alpha_init, gamma)
  marker_pi <- sblr::sbayesrc_marker_pi(A, init$alpha_init, gamma)
  marker_gamma <- sblr::sbayesrc_marker_gamma_mean(A, init$alpha_init, gamma)

  expect_equal(dim(ann_pi), c(ncol(A), length(gamma)))
  expect_equal(dim(marker_pi), c(nrow(A), length(gamma)))
  expect_identical(rownames(ann_pi), colnames(A))
  expect_identical(rownames(marker_pi), rownames(A))
  expect_equal(unname(rowSums(ann_pi)), rep(1, nrow(ann_pi)), tolerance = 1e-12)
  expect_equal(unname(rowSums(marker_pi)), rep(1, nrow(marker_pi)), tolerance = 1e-12)
  expect_true(all(ann_pi >= 0 & ann_pi <= 1))
  expect_true(all(marker_pi >= 0 & marker_pi <= 1))
  expect_length(ann_gamma, ncol(A))
  expect_length(marker_gamma, nrow(A))
  expect_identical(names(ann_gamma), colnames(A))
  expect_identical(names(marker_gamma), rownames(A))
  expect_true(all(is.finite(ann_gamma)))
  expect_true(all(is.finite(marker_gamma)))
})

test_that("summarize_annotation_signal reports truth and fit summaries", {
  annot <- matrix(
    c(1, 0, 0, 1, 1, 1, 0, 0),
    nrow = 4,
    byrow = TRUE,
    dimnames = list(paste0("m", 1:4), c("coding", "qtl"))
  )
  sim <- list(
    annot = annot,
    causal_any = c(TRUE, FALSE, TRUE, FALSE)
  )
  fit <- list(
    dm = matrix(c(0.2, 0.1, 0.5, 0.0), ncol = 1,
                dimnames = list(rownames(annot), "trait1")),
    bm = matrix(c(0.3, -0.2, 0.4, 0.0), ncol = 1,
                dimnames = list(rownames(annot), "trait1"))
  )

  out_truth <- sblr::summarize_annotation_signal(sim)
  out_fit <- sblr::summarize_annotation_signal(sim, fit)

  expect_equal(nrow(out_truth), ncol(annot))
  expect_true(all(c("annotation", "size", "n_causal", "causal_rate_in_set") %in% names(out_truth)))
  expect_true(all(c("mean_dm_D1", "mean_abs_bm_D1") %in% names(out_fit)))
  expect_true(all(is.finite(out_fit$mean_dm_D1)))
  expect_true(all(is.finite(out_fit$mean_abs_bm_D1)))
})

test_that("csr_prior_bayesc wrapper returns standard marker summaries", {
  skip_if_not(
    exists("stblr_cpg_omp_csr_prior", mode = "function"),
    "native fixed-prior annotation CSR symbol is not loaded"
  )

  stats <- tiny_annotation_stats()
  A <- tiny_annotation_matrix()
  fit <- sblr::stblr_csr_prior_annot(
    stats = stats,
    ld_prefix = make_tiny_annotation_csr_prefix(stats$m),
    A = A,
    fixed_pi_marker = list(rep(0.35, stats$m)),
    fixed_vb_multiplier = list(c(1, 1.2, 0.8, 1)),
    use_pi_marker = TRUE,
    use_vb_multiplier = TRUE,
    pi_init = 0.35,
    pi_prior_mean = 0.35,
    pi_prior_strength = 2,
    updateB = FALSE,
    updateE = FALSE,
    updatePi = FALSE,
    nit = 2,
    nburn = 0,
    ncores = 1L,
    seed = 1L
  )

  expect_annotation_fit_core(fit, stats)
  expect_equal(fit$input$model, "prior")
  expect_identical(fit$input$annotation_names, colnames(A))
  expect_true(isTRUE(fit$input$use_pi_marker))
  expect_true(isTRUE(fit$input$use_vb_multiplier))
  expect_length(fit$input$pi_marker, 1)
  expect_length(fit$input$vb_multiplier, 1)
  expect_trace_matrix(fit$vbs, 1L)
  expect_trace_matrix(fit$vgs, 1L)
  expect_trace_matrix(fit$ves, 1L)
  expect_square_trait_matrix(fit$covb, 1L)
  # Target metadata for the next alignment task:
  # input$method = "bayesc", input$backend = "csr_prior_bayesc",
  # input$data_level = "summary", and input$annotation_model = "prior".
})

test_that("csr_annot_bayesc wrapper returns learned annotation fields", {
  skip_if_not(
    exists("stblr_cpg_omp_csr_annot", mode = "function"),
    "native learned-annotation CSR symbol is not loaded"
  )

  stats <- tiny_annotation_stats()
  A <- tiny_annotation_matrix()
  fit <- sblr::stblr_csr_learn_annot(
    stats = stats,
    ld_prefix = make_tiny_annotation_csr_prefix(stats$m),
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
    nit = 2,
    nburn = 0,
    ncores = 1L,
    seed = 2L
  )

  expect_annotation_fit_core(fit, stats)
  expect_equal(fit$input$model, "annot")
  expect_identical(fit$input$annotation_names, colnames(A))
  expect_true(is.matrix(fit$eta_pi))
  expect_true(is.matrix(fit$eta_vb))
  expect_equal(dim(fit$eta_pi), c(1L, ncol(A)))
  expect_equal(dim(fit$eta_vb), c(1L, ncol(A)))
  expect_identical(colnames(fit$eta_pi), colnames(A))
  expect_identical(colnames(fit$eta_vb), colnames(A))
  expect_true(all(is.finite(fit$eta_pi)))
  expect_true(all(is.finite(fit$eta_vb)))
  expect_trace_matrix(fit$vle, 1L)
  expect_trace_matrix(fit$vld, 1L)
  # Target aliases for the next alignment task: annotation_effects should
  # expose eta_pi/eta_vb while preserving these current native field names.
})

test_that("csr_group_bayesc wrapper returns group-level annotation fields", {
  skip_if_not(
    exists("stblr_cpg_omp_csr_group_annot", mode = "function"),
    "native group annotation CSR symbol is not loaded"
  )

  stats <- tiny_annotation_stats()
  group <- stats::setNames(c("coding", "background", "coding", "background"), stats$marker_names)
  fit <- sblr::stblr_csr_group_annot(
    stats = stats,
    ld_prefix = make_tiny_annotation_csr_prefix(stats$m),
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
    nit = 2,
    nburn = 0,
    ncores = 1L,
    seed = 3L
  )

  expect_annotation_fit_core(fit, stats)
  expect_equal(fit$input$model, "group")
  expect_identical(fit$input$group_names, c("coding", "background"))
  for (nm in c("group_pi", "group_vb_multiplier", "group_nincluded", "group_size")) {
    expect_true(is.matrix(fit[[nm]]), info = nm)
    expect_equal(dim(fit[[nm]]), c(1L, 2L), info = nm)
    expect_identical(colnames(fit[[nm]]), c("coding", "background"), info = nm)
    expect_true(all(is.finite(fit[[nm]])), info = nm)
  }
  expect_true(all(fit$group_pi >= 0 & fit$group_pi <= 1))
  expect_equal(unname(fit$group_size[1, ]), c(2, 2))
  # Target aliases for the next alignment task: annotation_pi and
  # annotation_variance should point to group summaries.
})

test_that("csr_sbayesrc wrapper returns component and annotation fields", {
  skip_if_not(
    exists("stblr_cpg_omp_csr_sbayesrc", mode = "function"),
    "native SBayesRC CSR symbol is not loaded"
  )

  stats <- tiny_annotation_stats()
  A <- tiny_annotation_matrix()
  gamma <- c(0, 0.1, 1)
  fit <- sblr::stblr_csr_sbayesrc_generic(
    stats = stats,
    ld_prefix = make_tiny_annotation_csr_prefix(stats$m),
    A = A,
    gamma = gamma,
    pi_init = 0.35,
    pi_prior_mean = 0.35,
    pi_prior_strength = 2,
    updateAlpha = FALSE,
    updateB = FALSE,
    updateE = FALSE,
    nit = 2,
    nburn = 0,
    ncores = 1L,
    seed = 4L
  )

  expect_annotation_fit_core(fit, stats)
  expect_equal(fit$input$model, "sbayesrc")
  expect_true("comp_prob" %in% names(fit))
  expect_length(fit$comp_prob, 1L)
  expect_equal(dim(fit$comp_prob$trait1), c(stats$m, length(gamma)))
  expect_equal(unname(rowSums(fit$comp_prob$trait1)), rep(1, stats$m), tolerance = 1e-8)
  expect_true(all(fit$comp_prob$trait1 >= 0 & fit$comp_prob$trait1 <= 1))
  expect_length(fit$alpha, 1L)
  expect_equal(dim(fit$alpha$trait1), c(ncol(fit$input$A), length(gamma) - 1L))
  expect_equal(dim(fit$sigmaSqAlpha), c(1L, length(gamma) - 1L))
  expect_equal(dim(fit$ncomp), c(1L, length(gamma)))
  expect_true(all(is.finite(fit$alpha$trait1)))
  expect_true(all(is.finite(fit$sigmaSqAlpha)))
  expect_true(all(is.finite(fit$ncomp)))
  # Target metadata for the next alignment task:
  # input$method = "sbayesrc" or a documented BayesR-like value,
  # input$backend = "csr_sbayesrc", input$data_level = "summary".
})
