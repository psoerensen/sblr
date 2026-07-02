make_selection_s_csr_prefix <- function() {
  prefix <- tempfile("tiny_selection_s_csr_")
  row_ptr <- c(0, 1, 1, 1)
  col_idx <- 1L
  values <- 0.4

  sblr:::.stblr_write_uint64_file(paste0(prefix, ".row_ptr.u64.bin"), row_ptr)
  sblr:::.stblr_write_uint32_file(
    paste0(prefix, ".col_idx.u32.0based.bin"),
    col_idx
  )
  writeBin(
    as.numeric(values),
    paste0(prefix, ".values.f32.bin"),
    size = 4,
    endian = "little"
  )

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

selection_s_csr_stats <- function() {
  markers <- paste0("m", 1:3)
  list(
    wy = list(trait1 = stats::setNames(c(8, 3, 0.4), markers)),
    ww = list(trait1 = stats::setNames(rep(100, 3), markers)),
    yy = stats::setNames(100, "trait1"),
    n = 100L,
    m = 3L,
    marker_names = markers,
    trait_names = "trait1"
  )
}

selection_s_csr_glist <- function() {
  list(
    rsidsLD = list(c("m1", "m2", "m3")),
    rsids = list(c("m3", "m1", "m2")),
    maf = list(c(0.40, 0.05, 0.20))
  )
}

make_selection_s_csr_fit <- function(selection_s = NULL, updateLDswap = FALSE,
                                     seed = 11) {
  stblr_csr(
    Glist = selection_s_csr_glist(),
    stats = selection_s_csr_stats(),
    ld_prefix = make_selection_s_csr_prefix(),
    selection_s = selection_s,
    updateLDswap = updateLDswap,
    ld_swap_prob = 1,
    ld_swap_r2 = 0.01,
    nit = 8,
    nburn = 2,
    nthin = 1,
    ncores = 1,
    seed = seed,
    pi_init = 0.5,
    pi_prior_mean = 0.5,
    pi_prior_strength = 2,
    updateB = FALSE,
    updateE = FALSE,
    updatePi = FALSE
  )
}

make_selection_s_csr_bayesr_fit <- function(selection_s = NULL,
                                            updateLDswap = FALSE,
                                            seed = 11) {
  stblr_csr(
    Glist = selection_s_csr_glist(),
    stats = selection_s_csr_stats(),
    ld_prefix = make_selection_s_csr_prefix(),
    method = "bayesR",
    selection_s = selection_s,
    updateLDswap = updateLDswap,
    ld_swap_prob = 1,
    ld_swap_r2 = 0.01,
    nit = 8,
    nburn = 2,
    nthin = 1,
    ncores = 1,
    seed = seed,
    mixture_var = c(0, 0.1, 1),
    pi = c(0.4, 0.3, 0.3),
    alpha = c(1, 1, 1),
    updateB = FALSE,
    updateE = FALSE,
    updatePi = FALSE
  )
}

expect_selection_s_bayesr_dm_matches_component0 <- function(fit,
                                                            tolerance = 1e-12) {
  for (trait in names(fit$comp_prob)) {
    expect_equal(
      unname(as.numeric(fit$dm[, trait])),
      unname(as.numeric(1 - fit$comp_prob[[trait]][, "component_0"])),
      tolerance = tolerance
    )
  }
}

test_that("fixed selection_s is optional and preserves default CSR BayesC behavior", {
  skip_if_not(
    exists("stblr_cpg_omp_csr", mode = "function"),
    "native BayesC CSR symbol is not loaded"
  )

  args <- list(
    Glist = selection_s_csr_glist(),
    stats = selection_s_csr_stats(),
    ld_prefix = make_selection_s_csr_prefix(),
    method = "bayesC",
    nit = 8,
    nburn = 2,
    nthin = 1,
    ncores = 1,
    seed = 101,
    pi_init = 0.5,
    pi_prior_mean = 0.5,
    pi_prior_strength = 2,
    updateB = FALSE,
    updateE = FALSE,
    updatePi = FALSE
  )
  fit_omitted <- do.call(stblr_csr, args)
  fit_null <- do.call(stblr_csr, c(args, list(selection_s = NULL)))

  expect_equal(fit_null$dm, fit_omitted$dm)
  expect_equal(fit_null$bm, fit_omitted$bm)
  expect_equal(fit_null$vbs, fit_omitted$vbs)
  expect_equal(fit_null$vgs, fit_omitted$vgs)
  expect_equal(fit_null$ves, fit_omitted$ves)
  expect_null(fit_null$input$selection_s)
  expect_false(fit_null$input$selection_s_fixed)
  expect_null(fit_null$input$selection_s_exponent)
  expect_equal(fit_null$input$selection_s_scale, "standardized_genotype_effect")
})

test_that("fixed selection_s CSR BayesC fits return finite outputs and metadata", {
  skip_if_not(
    exists("stblr_cpg_omp_csr", mode = "function"),
    "native BayesC CSR symbol is not loaded"
  )

  for (s in c(0, -0.25)) {
    fit <- make_selection_s_csr_fit(selection_s = s, seed = 110 + round(100 * s))
    expect_true(all(is.finite(fit$dm)))
    expect_true(all(is.finite(fit$bm)))
    expect_true(all(is.finite(fit$vle)))
    expect_true(all(is.finite(fit$vld)))
    expect_equal(fit$input$selection_s, s)
    expect_true(fit$input$selection_s_fixed)
    expect_equal(fit$input$selection_s_exponent, s + 1)
    expect_equal(fit$input$selection_s_scale, "standardized_genotype_effect")
  }
})

test_that("selection_s = -1 gives unit prior scale and matches default", {
  skip_if_not(
    exists("stblr_cpg_omp_csr", mode = "function"),
    "native BayesC CSR symbol is not loaded"
  )

  fit_default <- make_selection_s_csr_fit(seed = 121)
  fit_s_minus_one <- make_selection_s_csr_fit(selection_s = -1, seed = 121)

  expect_equal(fit_s_minus_one$dm, fit_default$dm)
  expect_equal(fit_s_minus_one$bm, fit_default$bm)
  expect_equal(fit_s_minus_one$vbs, fit_default$vbs)
  expect_equal(fit_s_minus_one$ves, fit_default$ves)
  expect_equal(fit_s_minus_one$input$selection_s, -1)
  expect_true(fit_s_minus_one$input$selection_s_fixed)
  expect_equal(fit_s_minus_one$input$selection_s_exponent, 0)
})

test_that("fixed selection_s is optional and preserves default CSR BayesR behavior", {
  skip_if_not(
    exists("stblr_cpg_omp_csr_bayesr", mode = "function"),
    "native BayesR CSR symbol is not loaded"
  )

  fit_omitted <- make_selection_s_csr_bayesr_fit(seed = 141)
  fit_null <- make_selection_s_csr_bayesr_fit(selection_s = NULL, seed = 141)

  expect_equal(fit_null$dm, fit_omitted$dm)
  expect_equal(fit_null$bm, fit_omitted$bm)
  expect_equal(fit_null$vbs, fit_omitted$vbs)
  expect_equal(fit_null$vle, fit_omitted$vle)
  expect_equal(fit_null$vld, fit_omitted$vld)
  expect_equal(fit_null$comp_prob, fit_omitted$comp_prob)
  expect_null(fit_null$input$selection_s)
  expect_false(fit_null$input$selection_s_fixed)
  expect_null(fit_null$input$selection_s_exponent)
  expect_equal(fit_null$input$selection_s_scale, "standardized_genotype_effect")
})

test_that("selection_s = -1 gives unit prior scale and matches default CSR BayesR", {
  skip_if_not(
    exists("stblr_cpg_omp_csr_bayesr", mode = "function"),
    "native BayesR CSR symbol is not loaded"
  )

  fit_default <- make_selection_s_csr_bayesr_fit(seed = 142)
  fit_s_minus_one <- make_selection_s_csr_bayesr_fit(selection_s = -1, seed = 142)

  expect_equal(fit_s_minus_one$dm, fit_default$dm)
  expect_equal(fit_s_minus_one$bm, fit_default$bm)
  expect_equal(fit_s_minus_one$vbs, fit_default$vbs)
  expect_equal(fit_s_minus_one$vle, fit_default$vle)
  expect_equal(fit_s_minus_one$vld, fit_default$vld)
  expect_equal(fit_s_minus_one$comp_prob, fit_default$comp_prob)
  expect_equal(fit_s_minus_one$input$selection_s, -1)
  expect_true(fit_s_minus_one$input$selection_s_fixed)
  expect_equal(fit_s_minus_one$input$selection_s_exponent, 0)
})

test_that("fixed selection_s CSR BayesR fits return finite outputs and component probabilities", {
  skip_if_not(
    exists("stblr_cpg_omp_csr_bayesr", mode = "function"),
    "native BayesR CSR symbol is not loaded"
  )

  for (s in c(0, -0.25)) {
    fit <- make_selection_s_csr_bayesr_fit(selection_s = s, seed = 150 + round(100 * s))
    expect_true(all(is.finite(fit$dm)))
    expect_true(all(is.finite(fit$bm)))
    expect_true(all(is.finite(fit$vle)))
    expect_true(all(is.finite(fit$vld)))
    expect_true(all(vapply(fit$comp_prob, function(cp) all(is.finite(cp)), logical(1))))
    expect_true(all(vapply(fit$comp_prob, function(cp) {
      all(abs(rowSums(cp) - 1) < 1e-8)
    }, logical(1))))
    expect_selection_s_bayesr_dm_matches_component0(fit, tolerance = 1e-8)
    expect_equal(fit$input$selection_s, s)
    expect_true(fit$input$selection_s_fixed)
    expect_equal(fit$input$selection_s_exponent, s + 1)
    expect_equal(fit$input$selection_s_scale, "standardized_genotype_effect")
  }
})

test_that("selection_s validates fixed-S inputs and unsupported CSR backends", {
  expect_error(
    make_selection_s_csr_fit(selection_s = c(0, 1)),
    "selection_s must be NULL or a finite numeric scalar"
  )
  expect_error(
    make_selection_s_csr_fit(selection_s = NA_real_),
    "selection_s must be NULL or a finite numeric scalar"
  )
  expect_error(
    make_selection_s_csr_fit(selection_s = NaN),
    "selection_s must be NULL or a finite numeric scalar"
  )
  expect_error(
    make_selection_s_csr_fit(selection_s = Inf),
    "selection_s must be NULL or a finite numeric scalar"
  )
  expect_error(
    make_selection_s_csr_fit(selection_s = "0"),
    "selection_s must be NULL or a finite numeric scalar"
  )
  expect_error(
    make_selection_s_csr_bayesr_fit(selection_s = c(0, 1)),
    "selection_s must be NULL or a finite numeric scalar"
  )
  expect_error(
    make_selection_s_csr_bayesr_fit(selection_s = NA_real_),
    "selection_s must be NULL or a finite numeric scalar"
  )
  expect_error(
    stblr_csr(
      Glist = selection_s_csr_glist(),
      stats = selection_s_csr_stats(),
      ld_prefix = make_selection_s_csr_prefix(),
      scheduled = TRUE,
      selection_s = 0,
      nit = 2,
      nburn = 0
    ),
    "selection_s is currently supported only for unscheduled CSR BayesC/BayesR"
  )
})

test_that("selection_s validates MAF alignment to CSR LD marker order", {
  bad_glist <- selection_s_csr_glist()
  bad_glist$rsids[[1]] <- c("m3", "m1", "missing")
  expect_error(
    stblr_csr(
      Glist = bad_glist,
      stats = selection_s_csr_stats(),
      ld_prefix = make_selection_s_csr_prefix(),
      selection_s = 0,
      nit = 2,
      nburn = 0
    ),
    "Could not align MAF to LD marker order for selection_s"
  )
})

test_that("selection_s works with CSR BayesC LD-swap and backend consistency", {
  skip_if_not(
    exists("stblr_cpg_omp_csr", mode = "function"),
    "native BayesC CSR symbol is not loaded"
  )
  if (!exists("check_stblr_backend_consistency", mode = "function")) {
    source_sblr_test_file("R/check-stblr-backend-consistency.R")
  }

  fit <- make_selection_s_csr_fit(selection_s = 0, updateLDswap = TRUE, seed = 131)

  expect_true(all(is.finite(fit$dm)))
  expect_true(all(is.finite(fit$bm)))
  expect_true(all(is.finite(fit$vle)))
  expect_true(all(is.finite(fit$vld)))
  expect_s3_class(fit$ld_swap, "data.frame")
  expect_true(all(c("attempted", "accepted", "acceptance_rate") %in%
                    names(fit$ld_swap)))

  chk <- check_stblr_backend_consistency(fit, require_ld_swap = TRUE, verbose = FALSE)
  expect_true(all(chk$checks$ok))
})

test_that("selection_s works with CSR BayesR LD-swap and backend consistency", {
  skip_if_not(
    exists("stblr_cpg_omp_csr_bayesr", mode = "function"),
    "native BayesR CSR symbol is not loaded"
  )
  if (!exists("check_stblr_backend_consistency", mode = "function")) {
    source_sblr_test_file("R/check-stblr-backend-consistency.R")
  }

  fit <- make_selection_s_csr_bayesr_fit(
    selection_s = 0,
    updateLDswap = TRUE,
    seed = 161
  )

  expect_true(all(is.finite(fit$dm)))
  expect_true(all(is.finite(fit$bm)))
  expect_true(all(is.finite(fit$vle)))
  expect_true(all(is.finite(fit$vld)))
  expect_s3_class(fit$ld_swap, "data.frame")
  expect_true(all(c("attempted", "accepted", "acceptance_rate") %in%
                    names(fit$ld_swap)))
  expect_true(all(vapply(fit$comp_prob, function(cp) {
    all(abs(rowSums(cp) - 1) < 1e-8)
  }, logical(1))))
  expect_selection_s_bayesr_dm_matches_component0(fit, tolerance = 1e-8)

  chk <- check_stblr_backend_consistency(fit, require_ld_swap = TRUE, verbose = FALSE)
  expect_true(all(chk$checks$ok))
})

make_selection_s_subset_fit <- function(n_markers = 20) {
  markers <- paste0("m", seq_len(n_markers))
  trait_names <- c("D1", "D2")
  pip_d1 <- seq(0.001, 0.2, length.out = n_markers)
  pip_d2 <- rev(pip_d1)
  effect_d1 <- seq(0.01, 0.08, length.out = n_markers)
  effect_d2 <- rev(effect_d1)
  dm <- cbind(D1 = pip_d1, D2 = pip_d2)
  bm <- cbind(D1 = effect_d1, D2 = effect_d2)
  rownames(dm) <- markers
  rownames(bm) <- markers
  list(
    dm = dm,
    bm = bm,
    maf = stats::setNames(seq(0.05, 0.45, length.out = n_markers), markers),
    markers = markers,
    trait_names = trait_names
  )
}

test_that("summarise_stblr_maf_architecture works for a small CSR fit", {
  fit <- make_selection_s_csr_fit()

  out <- summarise_stblr_maf_architecture(
    fit,
    maf = c(m1 = 0.05, m2 = 0.2, m3 = 0.4)
  )

  expect_s3_class(out, "data.frame")
  expect_equal(out$trait, "trait1")
  expect_equal(out$method, "posthoc_regression")
  expect_equal(out$response, "log_pip_weighted_bm2")
  expect_equal(out$n_markers, 3)
  expect_equal(out$n_effective_markers, 3)
  expect_true(is.na(out$selection_s_posthoc))
  expect_true(is.na(out$intercept))
})

test_that("summarise_stblr_maf_architecture handles h and unweighted response", {
  fit <- list(
    dm = matrix(seq(0.1, 0.9, length.out = 6), ncol = 1,
                dimnames = list(paste0("m", 1:6), "trait1")),
    bm = matrix(seq(0.01, 0.06, length.out = 6), ncol = 1,
                dimnames = list(paste0("m", 1:6), "trait1"))
  )
  fit_before <- fit

  out <- summarise_stblr_maf_architecture(
    fit,
    h = stats::setNames(seq(0.095, 0.48, length.out = 6), paste0("m", 1:6)),
    use_pip_weights = FALSE
  )

  expect_equal(out$method, "posthoc_regression")
  expect_equal(out$response, "log_bm2")
  expect_true(all(is.finite(out$selection_s_posthoc)))
  expect_identical(fit, fit_before)
})

test_that("summarise_stblr_maf_architecture supports all-marker and filtered subsets", {
  fit <- make_selection_s_subset_fit()

  out_all <- summarise_stblr_maf_architecture(fit, maf = fit$maf)
  expect_equal(out_all$n_markers, c(20L, 20L))
  expect_equal(out_all$marker_filter, c("all", "all"))
  expect_true(all(c("se", "p_value", "r2", "marker_filter") %in% names(out_all)))
  expect_true(all(is.finite(out_all$selection_s_posthoc)))

  out_min <- summarise_stblr_maf_architecture(fit, maf = fit$maf, min_pip = 0.01)
  expect_equal(out_min$marker_filter, rep("min_pip=0.01", 2))
  expect_true(all(out_min$n_markers < out_all$n_markers))
  expect_true(all(out_min$n_markers >= 5))

  out_top <- summarise_stblr_maf_architecture(fit, maf = fit$maf, top_n = 10)
  expect_equal(out_top$n_markers, c(10L, 10L))
  expect_equal(out_top$marker_filter, rep("top_n=10", 2))
})

test_that("summarise_stblr_maf_architecture supports character marker subsets", {
  fit <- make_selection_s_subset_fit()
  marker_subset <- paste0("m", 1:8)

  out <- summarise_stblr_maf_architecture(
    fit,
    maf = unname(fit$maf),
    markers = marker_subset
  )

  expect_equal(out$n_markers, c(8L, 8L))
  expect_equal(out$n_effective_markers, c(8L, 8L))
  expect_equal(out$marker_filter, rep("markers=character", 2))
  expect_true(all(is.finite(out$selection_s_posthoc)))
})

test_that("summarise_stblr_maf_architecture supports trait-specific marker subsets", {
  fit <- make_selection_s_subset_fit()
  causal_by_trait <- list(
    D1 = paste0("m", 1:6),
    D2 = paste0("m", 7:14)
  )

  out <- summarise_stblr_maf_architecture(
    fit,
    maf = fit$maf,
    markers = causal_by_trait
  )

  expect_equal(out$n_markers, c(6L, 8L))
  expect_equal(out$marker_filter, rep("markers=list", 2))
  expect_true(all(is.finite(out$selection_s_posthoc)))
})

test_that("summarise_stblr_maf_architecture combines marker and min_pip filters", {
  fit <- make_selection_s_subset_fit()
  marker_subset <- paste0("m", 1:12)

  out <- summarise_stblr_maf_architecture(
    fit,
    maf = fit$maf,
    markers = marker_subset,
    min_pip = 0.05
  )

  expected_d1 <- sum(fit$dm[marker_subset, "D1"] >= 0.05)
  expected_d2 <- sum(fit$dm[marker_subset, "D2"] >= 0.05)
  expect_equal(out$n_markers, c(expected_d1, expected_d2))
  expect_equal(out$marker_filter, rep("markers=character;min_pip=0.05", 2))
})

test_that("summarise_stblr_maf_architecture returns NA statistics with fewer than five markers", {
  fit <- make_selection_s_subset_fit()

  out <- summarise_stblr_maf_architecture(
    fit,
    maf = fit$maf,
    markers = paste0("m", 1:4)
  )

  expect_equal(out$n_markers, c(4L, 4L))
  expect_true(all(is.na(out$selection_s_posthoc)))
  expect_true(all(is.na(out$se)))
  expect_true(all(is.na(out$p_value)))
  expect_true(all(is.na(out$r2)))
})

test_that("summarise_stblr_maf_architecture validates marker filter inputs", {
  fit <- make_selection_s_subset_fit()

  expect_error(
    summarise_stblr_maf_architecture(fit, maf = fit$maf, min_pip = -0.1),
    "min_pip"
  )
  expect_error(
    summarise_stblr_maf_architecture(fit, maf = fit$maf, min_pip = 1.1),
    "min_pip"
  )
  expect_error(
    summarise_stblr_maf_architecture(fit, maf = fit$maf, top_n = 0),
    "top_n"
  )
  expect_error(
    summarise_stblr_maf_architecture(fit, maf = fit$maf, top_n = 1.5),
    "top_n"
  )
})

test_that("summarise_stblr_maf_architecture handles missing MAF or h clearly", {
  fit <- list(
    dm = matrix(c(0.1, 0.5, 0.9), ncol = 1),
    bm = matrix(c(0.01, -0.02, 0.03), ncol = 1)
  )

  expect_error(
    summarise_stblr_maf_architecture(fit),
    "Either maf or h must be supplied or recoverable from fit"
  )
  expect_error(
    summarise_stblr_maf_architecture(fit, maf = c(0.1, 0.2)),
    "length equal to the number of markers"
  )
  expect_error(
    summarise_stblr_maf_architecture(fit, maf = c(0.1, 0.2, 0.7)),
    "maf must contain values"
  )
})
