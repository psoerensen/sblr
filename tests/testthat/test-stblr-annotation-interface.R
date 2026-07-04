make_tiny_annotation_interface_csr_prefix <- function(m = 4L) {
  prefix <- tempfile("tiny_annotation_interface_csr_")
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

tiny_annotation_interface_stats <- function(m = 4L) {
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

tiny_annotation_interface_matrix <- function(m = 4L) {
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

expect_annotation_interface_core <- function(fit, stats, backend, method, annotation_model) {
  expect_true(all(c("dm", "bm", "input") %in% names(fit)))
  expect_equal(dim(as.matrix(fit$dm)), c(stats$m, length(stats$yy)))
  expect_equal(dim(as.matrix(fit$bm)), c(stats$m, length(stats$yy)))
  expect_identical(rownames(fit$dm), stats$marker_names)
  expect_identical(colnames(fit$dm), names(stats$yy))
  expect_true(all(is.na(fit$dm) | (fit$dm >= -1e-8 & fit$dm <= 1 + 1e-8)))
  expect_true(all(is.na(fit$bm) | is.finite(fit$bm)))
  expect_length(colMeans(fit$dm, na.rm = TRUE), length(stats$yy))
  expect_length(apply(fit$dm, 2, max, na.rm = TRUE), length(stats$yy))

  expect_equal(fit$input$method, method)
  expect_equal(fit$input$model, method)
  expect_equal(fit$input$backend, backend)
  expect_equal(fit$input$data_level, "summary")
  expect_equal(fit$input$annotation_model, annotation_model)
  expect_true(isTRUE(fit$input$annotations))
  expect_equal(fit$input$nchains, 1L)
  expect_false(fit$input$keep_chains)
}

expect_annotation_raw_v1_top_level <- function(raw, backend) {
  expect_true(sblr:::.is_stblr_raw_v1(raw))
  expect_equal(raw$schema$class, "stblr_raw")
  expect_equal(as.integer(raw$schema$version), 1L)
  expect_equal(raw$meta$model, "bayesc")
  expect_equal(raw$meta$backend, backend)
  expect_true(all(c(
    "schema", "meta", "marker", "trace", "variance", "pi", "diagnostics",
    "chains", "prior", "group", "annotation", "component", "selection"
  ) %in% names(raw)))
}

make_annotation_raw_v1_common_args <- function() {
  stats <- tiny_annotation_interface_stats()
  list(
    stats = stats,
    wy = stats$wy,
    ww = stats$ww,
    yy = as.numeric(stats$yy),
    b_init = list(rep(0, stats$m)),
    d_init = list(rep(0, stats$m)),
    use_d_init = FALSE,
    r_init = list(rep(0, stats$m)),
    use_r_init = FALSE,
    rebuild_r_before_updateE = FALSE,
    ld_prefix = make_tiny_annotation_interface_csr_prefix(stats$m),
    B = matrix(0.2, 1, 1),
    E = matrix(0.8, 1, 1),
    ssb_prior = list(0.2),
    sse_prior = list(0.8),
    pi = c(0.65, 0.35),
    nub = 4,
    nue = 4,
    updateB = FALSE,
    updateE = FALSE,
    updatePi = FALSE,
    adjE = 0.9,
    n = 50L,
    nit = 2L,
    nburn = 0L,
    nthin = 1L,
    ncores = 1L,
    seed = 21L,
    nchains = 1L,
    keep_chains = FALSE,
    chain_seeds = NULL,
    updateLDswap = FALSE,
    ld_swap_prob = 0.05,
    ld_swap_r2 = 0.8,
    ld_swap_max_friends = 50L,
    ld_swap_moves = 1L
  )
}

test_that("stblr_csr_annot is exported and normalizes annotation model names", {
  expect_true("stblr_csr_annot" %in% getNamespaceExports("sblr"))

  expect_equal(sblr:::.stblr_match_annotation_model("prior"), "prior")
  expect_equal(sblr:::.stblr_match_annotation_model("fixed_prior"), "prior")
  expect_equal(sblr:::.stblr_match_annotation_model("fixed"), "prior")
  expect_equal(sblr:::.stblr_match_annotation_model("learned"), "learned")
  expect_equal(sblr:::.stblr_match_annotation_model("annot"), "learned")
  expect_equal(sblr:::.stblr_match_annotation_model("annotation"), "learned")
  expect_equal(sblr:::.stblr_match_annotation_model("group"), "group")
  expect_equal(sblr:::.stblr_match_annotation_model("groups"), "group")
  expect_equal(sblr:::.stblr_match_annotation_model("sbayesrc"), "sbayesrc")
  expect_equal(sblr:::.stblr_match_annotation_model("SBayesRC"), "sbayesrc")
})

test_that("marker-prior BayesC CSR backend returns raw v1 schema", {
  skip_if_not(
    exists("stblr_cpg_omp_csr_prior", mode = "function"),
    "native fixed-prior annotation CSR symbol is not loaded"
  )

  args <- make_annotation_raw_v1_common_args()
  raw <- do.call(stblr_cpg_omp_csr_prior, c(args[names(args) != "stats"], list(
    use_pi_marker = TRUE,
    pi_marker = list(rep(0.35, args$stats$m)),
    use_vb_multiplier = TRUE,
    vb_multiplier = list(c(1, 1.2, 0.8, 1)),
    pi_prior_a = 1,
    pi_prior_b = 1
  )))

  expect_annotation_raw_v1_top_level(raw, "csr_prior_bayesc")
  expect_true(all(c("marker_pi_mean", "marker_vb_multiplier_mean") %in% names(raw$prior)))
  expect_equal(dim(raw$prior$marker_pi_mean), c(args$stats$m, 1L))
  expect_null(raw$chains)

  fit <- sblr:::.format_stblr_raw_v1(
    raw,
    trait_names = args$stats$trait_names,
    variable_names = args$stats$marker_names
  )
  expect_true(all(c("bm", "dm", "vbs", "vgs", "ves", "vle", "vld", "pi", "pim", "pis") %in% names(fit)))
  expect_equal(dim(fit$bm), c(args$stats$m, 1L))
  expect_equal(dim(fit$pis), c(args$nit + args$nburn, 1L))
})

test_that("group-prior BayesC CSR backend returns raw v1 schema", {
  skip_if_not(
    exists("stblr_cpg_omp_csr_group_annot", mode = "function"),
    "native group annotation CSR symbol is not loaded"
  )

  args <- make_annotation_raw_v1_common_args()
  raw <- do.call(stblr_cpg_omp_csr_group_annot, c(args[names(args) != "stats"], list(
    group_index = c(0L, 1L, 0L, 1L),
    ngroup = 2L,
    group_pi_init = list(c(0.35, 0.25)),
    pi_group_prior_a = c(1, 1),
    pi_group_prior_b = c(1, 1),
    group_vb_multiplier_init = list(c(1.2, 0.8)),
    updateGroupVb = FALSE,
    nub_group = 4,
    ssb_group_prior = 1,
    normalize_group_vb = TRUE
  )))
  raw$group$group_names <- c("coding", "background")

  expect_annotation_raw_v1_top_level(raw, "csr_group_bayesc")
  expect_true(all(c("pi_mean", "vb_multiplier_mean", "n_included_mean", "size") %in% names(raw$group)))
  expect_equal(dim(raw$group$pi_mean), c(2L, 1L))
  expect_null(raw$chains)

  fit <- sblr:::.format_stblr_raw_v1(
    raw,
    trait_names = args$stats$trait_names,
    variable_names = args$stats$marker_names
  )
  expect_true(all(c("bm", "dm", "vbs", "vgs", "ves", "vle", "vld", "pi", "pim", "pis", "group_pi") %in% names(fit)))
  expect_equal(dim(fit$group_pi), c(1L, 2L))
})

test_that("learned-annotation BayesC CSR backend returns raw v1 schema", {
  skip_if_not(
    exists("stblr_cpg_omp_csr_annot", mode = "function"),
    "native learned-annotation CSR symbol is not loaded"
  )

  args <- make_annotation_raw_v1_common_args()
  A <- tiny_annotation_interface_matrix()
  raw <- do.call(stblr_cpg_omp_csr_annot, c(args[names(args) != "stats"], list(
    A = A,
    learn_pi_annot = TRUE,
    learn_vb_annot = TRUE,
    eta_pi_init = matrix(0, nrow = ncol(A), ncol = 1L),
    eta_vb_init = matrix(0, nrow = ncol(A), ncol = 1L),
    sigma_eta_pi = 1,
    sigma_eta_vb = 1,
    rw_sd_eta_pi = 0.01,
    rw_sd_eta_vb = 0.01,
    annot_update_every = 1L,
    pi_min = 1e-8,
    pi_max = 0.5,
    vb_multiplier_min = 1e-3,
    vb_multiplier_max = 1e3,
    pi_prior_a = 1,
    pi_prior_b = 1
  )))
  raw$annotation$annotation_names <- colnames(A)

  expect_annotation_raw_v1_top_level(raw, "csr_annot_bayesc")
  expect_true(all(c("eta_pi_mean", "eta_vb_mean") %in% names(raw$annotation)))
  expect_equal(dim(raw$annotation$eta_pi_mean), c(ncol(A), 1L))
  expect_null(raw$chains)

  fit <- sblr:::.format_stblr_raw_v1(
    raw,
    trait_names = args$stats$trait_names,
    variable_names = args$stats$marker_names
  )
  expect_true(all(c("bm", "dm", "vbs", "vgs", "ves", "vle", "vld", "pi", "pim", "pis", "eta_pi", "eta_vb") %in% names(fit)))
  expect_equal(dim(fit$eta_pi), c(1L, ncol(A)))
})

test_that("stblr_csr_annot dispatches fixed-prior BayesC annotations", {
  skip_if_not(
    exists("stblr_cpg_omp_csr_prior", mode = "function"),
    "native fixed-prior annotation CSR symbol is not loaded"
  )

  stats <- tiny_annotation_interface_stats()
  A <- tiny_annotation_interface_matrix()
  fit <- sblr::stblr_csr_annot(
    stats = stats,
    ld_prefix = make_tiny_annotation_interface_csr_prefix(stats$m),
    annotations = list(
      A = A,
      fixed_pi_marker = list(rep(0.35, stats$m)),
      fixed_vb_multiplier = list(c(1, 1.2, 0.8, 1))
    ),
    annotation_model = "fixed_prior",
    pi_init = 0.35,
    pi_prior_mean = 0.35,
    pi_prior_strength = 2,
    updateB = FALSE,
    updateE = FALSE,
    updatePi = FALSE,
    nit = 2,
    nburn = 0,
    ncores = 1L,
    seed = 11L
  )

  expect_annotation_interface_core(
    fit, stats, "csr_prior_bayesc", "bayesc", "prior"
  )
  expect_true(is.list(fit$annotation_prior))
  expect_true(is.data.frame(fit$annotation_summary))
  expect_true(all(c("covb", "vbs", "vgs", "ves") %in% names(fit)))
})

test_that("stblr_csr_annot dispatches learned BayesC annotations", {
  skip_if_not(
    exists("stblr_cpg_omp_csr_annot", mode = "function"),
    "native learned-annotation CSR symbol is not loaded"
  )

  stats <- tiny_annotation_interface_stats()
  A <- tiny_annotation_interface_matrix()
  fit <- sblr::stblr_csr_annot(
    stats = stats,
    ld_prefix = make_tiny_annotation_interface_csr_prefix(stats$m),
    annotations = A,
    annotation_model = "annot",
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
    seed = 12L
  )

  expect_annotation_interface_core(
    fit, stats, "csr_annot_bayesc", "bayesc", "learned"
  )
  expect_identical(fit$annotation_effects$pi, fit$eta_pi)
  expect_identical(fit$annotation_effects$variance, fit$eta_vb)
  expect_true(all(c("eta_pi", "eta_vb", "vle", "vld") %in% names(fit)))
})

test_that("stblr_csr_annot dispatches group BayesC annotations", {
  skip_if_not(
    exists("stblr_cpg_omp_csr_group_annot", mode = "function"),
    "native group annotation CSR symbol is not loaded"
  )

  stats <- tiny_annotation_interface_stats()
  group <- stats::setNames(
    c("coding", "background", "coding", "background"),
    stats$marker_names
  )
  fit <- sblr::stblr_csr_annot(
    stats = stats,
    ld_prefix = make_tiny_annotation_interface_csr_prefix(stats$m),
    annotations = group,
    annotation_model = "groups",
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
    seed = 13L
  )

  expect_annotation_interface_core(
    fit, stats, "csr_group_bayesc", "bayesc", "group"
  )
  expect_identical(fit$annotation_pi, fit$group_pi)
  expect_identical(fit$annotation_variance, fit$group_vb_multiplier)
  expect_true(all(c("group_pi", "group_vb_multiplier", "group_nincluded", "group_size") %in% names(fit)))
})

test_that("stblr_csr_annot dispatches SBayesRC annotations", {
  skip_if_not(
    exists("stblr_cpg_omp_csr_sbayesrc", mode = "function"),
    "native SBayesRC CSR symbol is not loaded"
  )

  stats <- tiny_annotation_interface_stats()
  A <- tiny_annotation_interface_matrix()
  gamma <- c(0, 0.1, 1)
  fit <- sblr::stblr_csr_annot(
    stats = stats,
    ld_prefix = make_tiny_annotation_interface_csr_prefix(stats$m),
    annotations = A,
    annotation_model = "SBayesRC",
    method = "bayesr",
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
    seed = 14L
  )

  expect_annotation_interface_core(
    fit, stats, "csr_sbayesrc", "sbayesrc", "sbayesrc"
  )
  expect_identical(fit$annotation_effects, fit$alpha)
  expect_identical(fit$annotation_variance, fit$sigmaSqAlpha)
  expect_length(fit$annotation_pi, 1L)
  expect_equal(unname(rowSums(fit$comp_prob$trait1)), rep(1, stats$m), tolerance = 1e-8)
  expect_true(all(c("comp_prob", "alpha", "sigmaSqAlpha", "ncomp") %in% names(fit)))
})
