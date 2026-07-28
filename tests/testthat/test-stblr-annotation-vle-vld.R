make_tiny_annotation_vle_vld_csr_prefix <- function(m = 4L) {
  prefix <- tempfile("tiny_annotation_vle_vld_csr_")
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

tiny_annotation_vle_vld_stats <- function(m = 4L) {
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

tiny_annotation_vle_vld_matrix <- function(m = 4L) {
  markers <- paste0("m", seq_len(m))
  matrix(
    c(1, 0, 0, 1, 1, 1, 0, 0),
    nrow = m,
    byrow = TRUE,
    dimnames = list(markers, c("coding", "qtl"))
  )
}

expect_vle_vld_trace_contract <- function(fit, stats) {
  for (nm in c("vle", "vld")) {
    expect_true(nm %in% names(fit), info = nm)
    expect_true(is.matrix(fit[[nm]]), info = nm)
    expect_equal(ncol(fit[[nm]]), length(stats$yy), info = nm)
    expect_identical(colnames(fit[[nm]]), names(stats$yy), info = nm)
    expect_identical(rownames(fit[[nm]]), paste0("Iter", seq_len(nrow(fit[[nm]]))), info = nm)
    expect_true(all(is.na(fit[[nm]]) | is.finite(fit[[nm]])), info = nm)
  }

  post <- sblr::summarise_posterior(
    fit,
    traces = c("vle", "vld"),
    derived = FALSE,
    include_diagnostics = FALSE
  )
  expect_setequal(post$parameter, c("vle", "vld"))

  chk <- sblr::check_stblr_consistency(fit, verbose = FALSE)
  expect_true(chk$ok)
  expect_true(all(c("trace.vle.finite", "trace.vld.finite") %in% chk$checks$check))
}

fit_tiny_annotation_vle_vld <- function(annotation_model,
                                        nchains = 1L,
                                        keep_chains = FALSE,
                                        updateLDswap = FALSE) {
  stats <- tiny_annotation_vle_vld_stats()
  A <- tiny_annotation_vle_vld_matrix()
  common <- list(
    stats = stats,
    ld_prefix = make_tiny_annotation_vle_vld_csr_prefix(stats$m),
    annotations = A,
    annotation_model = annotation_model,
    pi_init = 0.35,
    pi_prior_mean = 0.35,
    pi_prior_strength = 2,
    h2 = 0.3,
    updateB = FALSE,
    updateE = FALSE,
    nit = 3,
    nburn = 0,
    ncores = 1L,
    seed = 71L,
    nchains = nchains,
    keep_chains = keep_chains,
    chain_seeds = if (nchains > 1L) seq_len(nchains) + 100L else NULL,
    updateLDswap = updateLDswap,
    ld_swap_prob = 1,
    ld_swap_r2 = 0,
    ld_swap_moves = 1L
  )

  if (annotation_model == "fixed_marker") {
    common$annotations <- list(
      A = A,
      fixed_pi_marker = list(rep(0.35, stats$m)),
      fixed_vb_multiplier = list(c(1, 1.2, 0.8, 1))
    )
    common$updatePi <- FALSE
  } else if (annotation_model == "learned_logistic") {
    common$updatePi <- FALSE
    common$learn_pi_annot <- TRUE
    common$learn_vb_annot <- TRUE
    common$annot_update_every <- 1L
    common$rw_sd_eta_pi <- 0.01
    common$rw_sd_eta_vb <- 0.01
  } else if (annotation_model == "group") {
    common$annotations <- stats::setNames(
      c("coding", "background", "coding", "background"),
      stats$marker_names
    )
    common$group_names <- c("coding", "background")
    common$group_pi_init <- c(0.35, 0.25)
    common$group_vb_multiplier_init <- c(1.2, 0.8)
    common$updatePi <- FALSE
    common$updateGroupVb <- FALSE
  } else if (annotation_model == "annotation_probit_stick") {
    common$Glist <- list(
      rsidsLD = list(stats$marker_names),
      rsids = list(stats$marker_names),
      maf = list(rep(0.2, stats$m))
    )
    common$mixture_var <- c(0, 0.1, 1)
    common$updateAlpha <- FALSE
  }

  do.call(sblr::stblr_csr_annot, common)
}

test_that("annotation-aware CSR models expose vle and vld", {
  required <- c(
    fixed_marker = "stblr_cpg_omp_csr_prior",
    learned_logistic = "stblr_cpg_omp_csr_annot",
    group = "stblr_cpg_omp_csr_group_annot",
    annotation_probit_stick = "stblr_cpg_omp_csr_sbayesrc"
  )
  stats <- tiny_annotation_vle_vld_stats()
  backend_names <- c(
    fixed_marker = "prior",
    learned_logistic = "learned",
    group = "group",
    annotation_probit_stick = "sbayesrc"
  )

  for (model in names(required)) {
    skip_if_not(exists(required[[model]], mode = "function"))
    fit <- fit_tiny_annotation_vle_vld(model)
    expect_equal(fit$input$annotation_model, unname(backend_names[[model]]))
    expect_vle_vld_trace_contract(fit, stats)
  }
})

test_that("annotation-aware vle and vld survive chains and LD-swap", {
  required <- c(
    fixed_marker = "stblr_cpg_omp_csr_prior",
    learned_logistic = "stblr_cpg_omp_csr_annot",
    group = "stblr_cpg_omp_csr_group_annot",
    annotation_probit_stick = "stblr_cpg_omp_csr_sbayesrc"
  )
  stats <- tiny_annotation_vle_vld_stats()

  for (model in names(required)) {
    skip_if_not(exists(required[[model]], mode = "function"))
    fit <- fit_tiny_annotation_vle_vld(
      model,
      nchains = 2L,
      keep_chains = TRUE,
      updateLDswap = TRUE
    )
    expect_equal(fit$input$nchains, 2L)
    expect_true(isTRUE(fit$input$keep_chains))
    expect_true(isTRUE(fit$input$updateLDswap))
    expect_vle_vld_trace_contract(fit, stats)
  }
})
