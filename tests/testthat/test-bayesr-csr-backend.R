.format_stblr_csr_bayesr_fit <- getFromNamespace(".format_stblr_csr_bayesr_fit", "sblr")
check_stblr_consistency <- getFromNamespace("check_stblr_consistency", "sblr")
extract_stblr_finemap_loci <- getFromNamespace("extract_stblr_finemap_loci", "sblr")

make_bayesr_csr_raw <- function(nchains = 1L) {
  comp_prob <- matrix(
    c(
      0.70, 0.20, 0.10,
      0.10, 0.30, 0.60,
      0.40, 0.40, 0.20
    ),
    nrow = 3,
    byrow = TRUE
  )
  dm <- 1 - comp_prob[, 1L]
  bm <- c(0.01, -0.02, 0.03)
  bm_chain <- if (nchains == 1L) {
    list(bm)
  } else {
    list(c(0.00, -0.03, 0.02), c(0.02, -0.01, 0.04))
  }
  dm_chain <- if (nchains == 1L) {
    list(dm)
  } else {
    list(c(0.20, 0.80, 0.50), c(0.40, 1.00, 0.70))
  }
  bm_mat <- do.call(cbind, bm_chain)
  dm_mat <- do.call(cbind, dm_chain)
  trace_len <- 4L

  list(
    bm = matrix(rowMeans(bm_mat), ncol = 1),
    dm = matrix(rowMeans(dm_mat), ncol = 1),
    wy = matrix(rep(0, 3), ncol = 1),
    r = matrix(rep(0, 3), ncol = 1),
    b = matrix(bm, ncol = 1),
    component = matrix(c(0, 2, 1), ncol = 1),
    vbs = matrix(seq_len(trace_len) / 10, ncol = 1),
    vgs = matrix(seq_len(trace_len) / 9, ncol = 1),
    ves = matrix(seq_len(trace_len) / 8, ncol = 1),
    covb = matrix(1),
    covg = matrix(1),
    cove = matrix(1),
    vb = matrix(1),
    vg = matrix(1),
    ve = matrix(1),
    pi = matrix(c(0.40, 0.30, 0.30), nrow = 1),
    pim = matrix(c(0.45, 0.25, 0.30), nrow = 1),
    vle = matrix(seq_len(trace_len) / 7, ncol = 1),
    vld = matrix(seq_len(trace_len) / 6, ncol = 1),
    bm_sd = matrix(if (nchains == 1L) rep(0, 3) else apply(bm_mat, 1, stats::sd), ncol = 1),
    bm_min = matrix(apply(bm_mat, 1, min), ncol = 1),
    bm_max = matrix(apply(bm_mat, 1, max), ncol = 1),
    dm_sd = matrix(if (nchains == 1L) rep(0, 3) else apply(dm_mat, 1, stats::sd), ncol = 1),
    dm_min = matrix(apply(dm_mat, 1, min), ncol = 1),
    dm_max = matrix(apply(dm_mat, 1, max), ncol = 1),
    comp_prob = list(comp_prob),
    dm_component_mean = matrix(c(0.4, 1.5, 0.8), ncol = 1),
    ncomp = matrix(colSums(comp_prob), nrow = 1),
    mixture_var = c(0, 0.01, 0.1)
  )
}

format_bayesr_csr_test_fit <- function(raw, nchains) {
  fit <- .format_stblr_csr_bayesr_fit(
    raw,
    nt = 1L,
    m = 3L,
    trait_names = "D1",
    variable_names = paste0("m", 1:3),
    n_components = 3L
  )
  fit$input <- list(
    method = "sbayesr",
    model = "sbayesr",
    backend = "csr_bayesr",
    data_level = "summary",
    scheduled = FALSE,
    keep_chains = FALSE,
    updateLDswap = FALSE,
    nchains = nchains
  )
  fit$chains <- sblr:::.blr_flatten_st_chains(
    fit$chains, "D1", "sbayesr", "csr")
  sblr:::.blr_finalize_fit(
    fit, "stblr", "sbayesr", "csr",
    data = list(marker_ids = paste0("m", 1:3), trait_names = "D1"),
    diagnostics = list(native = fit$diagnostics))
}

test_that("CSR BayesR formatter exposes standard and component fields", {
  raw <- make_bayesr_csr_raw(nchains = 1L)
  fit <- format_bayesr_csr_test_fit(raw, nchains = 1L)

  for (nm in c(
    "dm", "bm", "dm_chain_mean_sd", "dm_chain_mean_min",
    "dm_chain_mean_max", "bm_chain_mean_sd", "bm_chain_mean_min",
    "bm_chain_mean_max", "component_probabilities", "dm_component_mean"
  )) {
    expect_true(nm %in% names(fit))
  }

  expect_equal(as.numeric(fit$dm[, "D1"]), 1 - raw$comp_prob[[1]][, 1], tolerance = 1e-12)
  expect_true(all(fit$dm >= -1e-12 & fit$dm <= 1 + 1e-12))
  expect_identical(rownames(fit$component_probabilities$D1), paste0("m", 1:3))
  expect_identical(colnames(fit$component_probabilities$D1), paste0("component_", 0:2))
  expect_equal(unname(rowSums(fit$component_probabilities$D1)), rep(1, 3), tolerance = 1e-12)
  expect_true(all(fit$component_probabilities$D1 >= -1e-12 & fit$component_probabilities$D1 <= 1 + 1e-12))
  expect_equal(unname(fit$dm_component_mean[, "D1"]), c(0.4, 1.5, 0.8))
})

test_that("CSR BayesR formatter exposes single-chain summary convention", {
  fit <- format_bayesr_csr_test_fit(make_bayesr_csr_raw(nchains = 1L), nchains = 1L)

  for (nm in c("bm_chain_mean_sd", "bm_chain_mean_min", "bm_chain_mean_max",
               "dm_chain_mean_sd", "dm_chain_mean_min", "dm_chain_mean_max")) {
    expect_equal(dim(fit[[nm]]), dim(fit$dm))
    expect_identical(rownames(fit[[nm]]), rownames(fit$dm))
    expect_identical(colnames(fit[[nm]]), colnames(fit$dm))
  }
  expect_equal(fit$bm_chain_mean_sd, fit$bm * 0, tolerance = 1e-12)
  expect_equal(fit$dm_chain_mean_sd, fit$dm * 0, tolerance = 1e-12)
  expect_equal(fit$bm_chain_mean_min, fit$bm, tolerance = 1e-12)
  expect_equal(fit$bm_chain_mean_max, fit$bm, tolerance = 1e-12)
  expect_equal(fit$dm_chain_mean_min, fit$dm, tolerance = 1e-12)
  expect_equal(fit$dm_chain_mean_max, fit$dm, tolerance = 1e-12)
})

test_that("CSR BayesR formatter exposes finite multi-chain summaries", {
  fit <- format_bayesr_csr_test_fit(make_bayesr_csr_raw(nchains = 2L), nchains = 2L)

  expect_true(all(is.finite(fit$bm_chain_mean_sd)))
  expect_true(all(is.finite(fit$dm_chain_mean_sd)))
  expect_true(all(fit$bm_chain_mean_sd >= -1e-12))
  expect_true(all(fit$dm_chain_mean_sd >= -1e-12))
  expect_true(all(fit$bm_chain_mean_min <= fit$bm + 1e-12))
  expect_true(all(fit$bm <= fit$bm_chain_mean_max + 1e-12))
  expect_true(all(fit$dm_chain_mean_min <= fit$dm + 1e-12))
  expect_true(all(fit$dm <= fit$dm_chain_mean_max + 1e-12))

  chk <- check_stblr_consistency(
    fit,
    require_chain_summaries = TRUE,
    verbose = FALSE
  )
  expect_true(chk$ok)
})

test_that("CSR BayesR formatted fit is compatible with fine-mapping extractor", {
  fit <- format_bayesr_csr_test_fit(make_bayesr_csr_raw(nchains = 1L), nchains = 1L)
  glist <- list(
    rsids = list(paste0("m", 1:3)),
    chr = list(rep(1L, 3)),
    pos = list(c(100, 200, 300)),
    sparseLD = list(chr = 1L, cls = list(1:3), prefix = NULL)
  )

  fm <- extract_stblr_finemap_loci(
    fit = fit,
    Glist = glist,
    locus_sets = list(regionA = paste0("m", 1:3)),
    trait = "D1",
    credible_sets = FALSE
  )

  expect_equal(fm$markers$pip_sd, as.numeric(fit$dm_chain_mean_sd[, "D1"]))
  expect_equal(fm$markers$pip_min, as.numeric(fit$dm_chain_mean_min[, "D1"]))
  expect_equal(fm$markers$pip_max, as.numeric(fit$dm_chain_mean_max[, "D1"]))
  expect_equal(fm$markers$bm_sd, as.numeric(fit$bm_chain_mean_sd[, "D1"]))
  expect_equal(fm$markers$bm_min, as.numeric(fit$bm_chain_mean_min[, "D1"]))
  expect_equal(fm$markers$bm_max, as.numeric(fit$bm_chain_mean_max[, "D1"]))
})

test_that("public CSR BayesR API exists and rejects unsupported modes early", {
  expect_true(is.function(stblr_csr_bayesr))
  expect_error(
    stblr_csr(stats = list(), method = "sbayesr", scheduled = TRUE),
    "scheduled CSR BayesR"
  )
})

test_that("CSR BayesR LD-swap arguments are validated", {
  expect_error(
    stblr_csr(stats = list(), method = "sbayesr", updateLDswap = NA),
    "updateLDswap"
  )
  expect_error(
    stblr_csr(stats = list(), method = "sbayesr", ld_swap_prob = 2),
    "ld_swap_prob"
  )
  expect_error(
    stblr_csr(stats = list(), method = "sbayesr", ld_swap_r2 = -0.1),
    "ld_swap_r2"
  )
  expect_error(
    stblr_csr(stats = list(), method = "sbayesr", ld_swap_max_friends = 0),
    "ld_swap_max_friends"
  )
  expect_error(
    stblr_csr(stats = list(), method = "sbayesr", ld_swap_moves = -1),
    "ld_swap_moves"
  )
})

write_empty_csr_ld_fixture <- function(prefix, m) {
  writeLines(c(
    paste0("n_variants=", m),
    "nnz=0"
  ), paste0(prefix, ".meta.txt"))
  writeBin(rep(as.raw(0), 8L * (m + 1L)), paste0(prefix, ".row_ptr.u64.bin"))
  file.create(paste0(prefix, ".col_idx.u32.0based.bin"))
  file.create(paste0(prefix, ".values.f32.bin"))
}

write_high_ld_csr_ld_fixture <- function(prefix, m = 4L) {
  row_ptr <- c(0, 1, 1, 1, 1)
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
      paste0("n_variants=", m),
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
}

make_small_bayesr_csr_stats <- function() {
  m <- 4L
  nt <- 1L
  n <- 80L
  variable_names <- paste0("m", seq_len(m))
  yy <- stats::setNames(n - 1, "T1")
  ww <- stats::setNames(rep(n - 1, m), variable_names)
  wy <- stats::setNames(c(4.0, -3.0, 2.0, 1.0), variable_names)

  list(
    stats = list(
      yy = yy,
      ww = list(T1 = ww),
      wy = list(T1 = wy),
      n = n,
      m = m
    ),
    Glist = list(sparseLD = list(prefix = tempfile("bayesr-csr-ld-"))),
    nt = nt,
    m = m
  )
}

expect_bayesr_csr_conventions <- function(fit) {
  chk <- check_stblr_consistency(
    fit,
    require_chain_summaries = TRUE,
    verbose = FALSE
  )
  expect_true(chk$ok)

  expect_true(all(is.finite(fit$dm)))
  expect_true(all(fit$dm >= -1e-12 & fit$dm <= 1 + 1e-12))
  for (trait in names(fit$component_probabilities)) {
    cp <- fit$component_probabilities[[trait]]
    expect_equal(
      unname(as.numeric(fit$dm[rownames(cp), trait])),
      unname(as.numeric(1 - cp[, 1L])),
      tolerance = 1e-12
    )
    expect_equal(unname(rowSums(cp)), rep(1, nrow(cp)), tolerance = 1e-12)
    expect_true(all(cp >= -1e-12 & cp <= 1 + 1e-12))
  }
}

expect_bayesr_csr_chain_aggregation <- function(fit) {
  chk <- check_stblr_consistency(
    fit,
    require_chain_summaries = TRUE,
    require_chains = TRUE,
    verbose = FALSE
  )
  expect_true(chk$ok)

  for (trait in colnames(fit$dm)) {
    chains <- Filter(function(ch) identical(ch$trait_name, trait), fit$chains)
    dm_mat <- do.call(cbind, lapply(chains, function(ch) ch$dm))
    bm_mat <- do.call(cbind, lapply(chains, function(ch) ch$bm))
    expect_equal(unname(rowMeans(dm_mat)), as.numeric(fit$dm[, trait]), tolerance = 1e-8)
    expect_equal(unname(rowMeans(bm_mat)), as.numeric(fit$bm[, trait]), tolerance = 1e-8)

    cp_mean <- Reduce(`+`, lapply(chains, function(ch) ch$component_probabilities)) / length(chains)
    expect_equal(cp_mean, fit$component_probabilities[[trait]], tolerance = 1e-8)
    for (ch in chains) {
      expect_true(all(c(
        "dm", "bm", "component_probabilities", "dm_component_mean",
        "pi_final", "pi_mean",
        "vbs", "vgs", "ves", "updateE_diagnostics"
      ) %in% names(ch)))
    }
  }
}

expect_bayesr_ld_swap_diagnostics <- function(fit, require_chains = FALSE) {
  chk <- check_stblr_consistency(
    fit,
    require_chain_summaries = TRUE,
    require_chains = require_chains,
    require_ld_swap = TRUE,
    verbose = FALSE
  )
  expect_true(chk$ok)
  expect_true("ld_swap" %in% names(fit$diagnostics))
  expect_true(all(c("attempted", "accepted", "acceptance_rate") %in%
                    names(fit$diagnostics$ld_swap)))
  expect_true(all(fit$diagnostics$ld_swap$attempted >=
                    fit$diagnostics$ld_swap$accepted))
  expect_true(all(fit$diagnostics$ld_swap$accepted >= 0))
  expect_true(all(fit$diagnostics$ld_swap$acceptance_rate >= 0))
  expect_true(all(fit$diagnostics$ld_swap$acceptance_rate <= 1))
}

test_that("CSR SBayesR extended diagnostics use native task-private states", {
  fixture <- make_small_bayesr_csr_stats()
  write_empty_csr_ld_fixture(fixture$Glist$sparseLD$prefix, fixture$m)
  fit <- stblr_csr(
    stats = fixture$stats, Glist = fixture$Glist, method = "sbayesr",
    mixture_var = c(0, .1, 1), updateB = FALSE, updateE = FALSE,
    nit = 8L, nburn = 2L, nchains = 2L, ncores = 1L, seed = 901L,
    convergence = "extended", convergence_control = list(
      warn = FALSE, extended_groups = "probability",
      selected_markers = c("m3", "m1"),
      selected_marker_quantities = c("b", "d", "component"),
      keep_traces = TRUE))
  summary <- fit$convergence$summary
  expect_true(all(c("component_pi", "selected_b", "selected_d",
                    "selected_component") %in% summary$group))
  expect_identical(unique(summary$marker_id[summary$tier == 3L]),
                   c("m3", "m1"))
  probability <- fit$convergence_traces$values[, ,
    summary$quantity[summary$group == "component_pi"], drop = FALSE]
  expect_equal(unname(apply(probability, c(1, 2), sum)), matrix(1, 8, 2),
               tolerance = 1e-12)
  expect_true(all(fit$convergence_traces$values[, ,
    summary$quantity[summary$group == "selected_d"]] %in% c(0, 1)))
})

test_that("CSR SBayesR diagnostic capture is RNG-neutral and unthinned", {
  fixture <- make_small_bayesr_csr_stats()
  write_empty_csr_ld_fixture(fixture$Glist$sparseLD$prefix, fixture$m)
  common <- list(
    stats = fixture$stats, Glist = fixture$Glist, method = "sbayesr",
    mixture_var = c(0, .1, 1), updateB = FALSE, updateE = FALSE,
    nit = 8L, nburn = 2L, nchains = 2L, ncores = 1L, seed = 902L)
  none <- do.call(stblr_csr, c(common, list(convergence = "none")))
  run <- function(nthin, keep_chains, keep_traces) do.call(stblr_csr, c(
    common, list(nthin = nthin, keep_chains = keep_chains,
      convergence = "extended", convergence_control = list(
        warn = FALSE, extended_groups = "probability",
        selected_markers = c("m3", "m1"),
        selected_marker_quantities = c("b", "d", "component"),
        keep_traces = keep_traces))))
  retained <- run(1L, TRUE, TRUE)
  unretained <- run(1L, FALSE, FALSE)
  thinned <- run(2L, FALSE, TRUE)
  for (field in c("bm", "dm", "b_final", "d_final", "pi_final", "pi_mean",
                  "component_probabilities", "vbs", "vgs", "ves", "vle",
                  "vld"))
    expect_equal(retained[[field]], none[[field]], tolerance = 0, info = field)
  expect_equal(unretained$convergence$summary,
               retained$convergence$summary, tolerance = 0)
  expect_equal(thinned$convergence_traces$values,
               retained$convergence_traces$values, tolerance = 0)
  expect_null(unretained$convergence_traces)
  expect_null(unretained$chains)
})

test_that("binary fixed component probabilities are deduplicated and not updated", {
  fixture <- make_small_bayesr_csr_stats()
  write_empty_csr_ld_fixture(fixture$Glist$sparseLD$prefix, fixture$m)
  fit <- stblr_csr(
    stats = fixture$stats, Glist = fixture$Glist, method = "sbayesr",
    mixture_var = c(0, 1), pi = c(.7, .3), updatePi = FALSE,
    updateB = FALSE, updateE = FALSE, nit = 8L, nburn = 2L,
    nchains = 2L, ncores = 1L, seed = 903L,
    convergence = "extended", convergence_control = list(
      warn = FALSE, extended_groups = "probability", keep_traces = TRUE))
  probability <- fit$convergence$summary[
    fit$convergence$summary$group == "component_pi", , drop = FALSE]
  expect_equal(nrow(probability), 1L)
  expect_false(probability$updated)
  expect_true(all(probability$status == "not_updated"))
  expect_identical(probability$component_name, "component_1")
})

test_that("stblr_csr dispatches exact CSR BayesR through method argument", {
  skip_if_not(
    exists("stblr_cpg_omp_csr_bayesr", mode = "function"),
    "native CSR BayesR symbol is not loaded"
  )

  fixture <- make_small_bayesr_csr_stats()
  write_empty_csr_ld_fixture(fixture$Glist$sparseLD$prefix, fixture$m)

  common <- list(
    stats = fixture$stats,
    Glist = fixture$Glist,
    h2 = 0.3,
    adjE = 0.9,
    nit = 12,
    nburn = 4,
    ncores = 1,
    nchains = 1,
    seed = 10,
    updateE = FALSE
  )
  common$convergence <- "none"
  direct_common <- common
  direct_common$convergence <- NULL
  fit_direct <- do.call(stblr_csr_bayesr, direct_common)
  fit_bridge <- do.call(stblr_csr, c(common, list(method = "sbayesr")))

  expect_bayesr_csr_conventions(fit_bridge)
  expect_equal(fit_bridge$input$method, "sbayesr")
  expect_equal(fit_bridge$input$model, "sbayesr")
  expect_equal(fit_bridge$input$backend, "csr_bayesr")
  expect_equal(fit_bridge$input$data_level, "summary_statistics")
  expect_equal(fit_bridge$input$scheduled, FALSE)
  expect_equal(dim(fit_bridge$dm), dim(fit_direct$dm))
  expect_equal(dim(fit_bridge$bm), dim(fit_direct$bm))
  expect_equal(names(fit_bridge$component_probabilities),
               names(fit_direct$comp_prob))
  required_fields <- c(
    "bm", "dm", "wy", "r", "b", "d",
    "vbs", "vgs", "ves", "vle", "vld", "pi_trace",
    "cov_b_mean", "cov_g_mean", "cov_e_mean", "pi_final", "pi_mean",
    "input", "component_probabilities", "dm_component_mean"
  )
  missing_fields <- setdiff(required_fields, names(fit_bridge))
  expect_equal(missing_fields, character())
  expect_null(fit_bridge$chains)
  expect_null(fit_bridge$diagnostics$ld_swap_chains)

  keep_args <- common
  keep_args$nchains <- 2
  keep_args$keep_chains <- TRUE
  fit_keep <- do.call(stblr_csr, c(keep_args, list(method = "sbayesr")))
  expect_bayesr_csr_conventions(fit_keep)
  expect_bayesr_csr_chain_aggregation(fit_keep)
  expect_true("chains" %in% names(fit_keep))
})

test_that("stblr_csr method BayesR rejects unsupported high-level combinations", {
  expect_error(
    stblr_csr(stats = list(), method = "sbayesr", scheduled = TRUE),
    "scheduled CSR BayesR is not currently implemented"
  )
  expect_error(
    stblr_csr(stats = list(), method = "sbayesr", pi_init = 0.5),
    "BayesC-specific"
  )
  expect_error(
    stblr_csr(stats = list(), method = "sbayesr", pi_prior_a = 1),
    "BayesC-specific"
  )
})

test_that("CSR BayesR formatter preserves compact per-chain component summaries", {
  raw <- make_bayesr_csr_raw(nchains = 2L)
  cp1 <- matrix(
    c(0.80, 0.20, 0.00, 0.20, 0.30, 0.50, 0.50, 0.50, 0.00),
    nrow = 3,
    byrow = TRUE
  )
  cp2 <- matrix(
    c(0.60, 0.20, 0.20, 0.00, 0.30, 0.70, 0.30, 0.30, 0.40),
    nrow = 3,
    byrow = TRUE
  )
  raw$chains <- list(list(
    list(
      dm = 1 - cp1[, 1L],
      bm = c(0.00, -0.03, 0.02),
      comp_prob = cp1,
      dm_component_mean = c(0.2, 1.3, 0.5),
      final_pi = c(0.50, 0.25, 0.25),
      mean_pi = c(0.55, 0.20, 0.25),
      vbs = 1:4,
      vgs = 2:5,
      ves = 3:6,
      updateE_diagnostics = c(0, 0, 4, 1.1, 2, 1.2, 2, 0.1, 0.2)
    ),
    list(
      dm = 1 - cp2[, 1L],
      bm = c(0.02, -0.01, 0.04),
      comp_prob = cp2,
      dm_component_mean = c(0.6, 1.7, 1.1),
      final_pi = c(0.30, 0.35, 0.35),
      mean_pi = c(0.35, 0.30, 0.35),
      vbs = 2:5,
      vgs = 3:6,
      ves = 4:7,
      updateE_diagnostics = c(0, 1, 4, 1.0, 1, 1.1, 2, 0.2, 0.3)
    )
  ))

  fit <- format_bayesr_csr_test_fit(raw, nchains = 2L)
  fit$input$keep_chains <- TRUE

  expect_bayesr_csr_chain_aggregation(fit)
  expect_identical(
    colnames(fit$chains[[1L]]$updateE_diagnostics),
    c(
      "trait_index", "chain_index", "n_updateE", "min_sse",
      "min_sse_iter", "min_residual_scale", "max_nonzero_components",
      "max_abs_effect", "max_fitted_quadratic"
    )
  )
})

test_that("supported exact CSR BayesR public API supports strict updateE modes", {
  skip_if_not(
    exists("stblr_cpg_omp_csr_bayesr", mode = "function"),
    "native CSR BayesR symbol is not loaded"
  )

  fixture <- make_small_bayesr_csr_stats()
  write_empty_csr_ld_fixture(fixture$Glist$sparseLD$prefix, fixture$m)

  fit_noE <- stblr_csr(
    method = "sbayesr",
    stats = fixture$stats,
    Glist = fixture$Glist,
    h2 = 0.3,
    adjE = 0.9,
    nit = 12,
    nburn = 4,
    ncores = 1,
    nchains = 1,
    seed = 10,
    updateE = FALSE
  )
  expect_bayesr_csr_conventions(fit_noE)
  expect_null(fit_noE$diagnostics$ld_swap)
  expect_equal(sum(fit_noE$input$pi[-1L]), 0.001, tolerance = 1e-12)
  expect_equal(
    unname(fit_noE$input$alpha / sum(fit_noE$input$alpha)),
    unname(fit_noE$input$pi),
    tolerance = 1e-12
  )
  expect_equal(fit_noE$input$updateE_start, 0L)
  expect_equal(fit_noE$input$updateE_every, 1L)
  expect_equal(fit_noE$input$model, "sbayesr")
  expect_equal(fit_noE$input$backend, "csr_bayesr")
  expect_equal(fit_noE$input$data_level, "summary_statistics")
  expect_equal(fit_noE$input$scheduled, FALSE)
  expect_equal(fit_noE$input$keep_chains, FALSE)
  expect_equal(fit_noE$input$updateLDswap, FALSE)

  fit_E <- stblr_csr(
    method = "sbayesr",
    stats = fixture$stats,
    Glist = fixture$Glist,
    h2 = 0.3,
    adjE = 0.9,
    nit = 20,
    nburn = 5,
    ncores = 1,
    nchains = 1,
    seed = 10,
    updateE = TRUE
  )
  expect_bayesr_csr_conventions(fit_E)
  expect_equal(fit_E$input$updateE, TRUE)
  expect_equal(fit_E$input$updateE_start, 0L)
  expect_equal(fit_E$input$updateE_every, 1L)
  expect_true("updateE_diagnostics" %in% names(fit_E))
  expect_true("min_sse_iter" %in% colnames(fit_E$updateE_diagnostics))
  expect_equal(
    unname(fit_E$updateE_diagnostics[, "n_updateE"]),
    25,
    tolerance = 1e-12
  )
  expect_true(all(fit_E$updateE_diagnostics[, "min_residual_scale"] > 0))
  expect_true(all(fit_E$updateE_diagnostics[, "min_sse_iter"] >= 0))

  fit_E_two_chain <- stblr_csr(
    method = "sbayesr",
    stats = fixture$stats,
    Glist = fixture$Glist,
    h2 = 0.3,
    adjE = 0.9,
    nit = 20,
    nburn = 5,
    ncores = 1,
    nchains = 2,
    convergence = "none",
    seed = 10,
    updateE = TRUE
  )
  expect_bayesr_csr_conventions(fit_E_two_chain)
  expect_equal(nrow(fit_E_two_chain$updateE_diagnostics), 2L)

  fit_E_keep <- stblr_csr(
    method = "sbayesr",
    stats = fixture$stats,
    Glist = fixture$Glist,
    h2 = 0.3,
    adjE = 0.9,
    nit = 20,
    nburn = 5,
    ncores = 1,
    nchains = 2,
    keep_chains = TRUE,
    convergence = "none",
    seed = 10,
    updateE = TRUE
  )
  expect_bayesr_csr_conventions(fit_E_keep)
  expect_equal(fit_E_keep$input$keep_chains, TRUE)
  expect_bayesr_csr_chain_aggregation(fit_E_keep)

  fit_E_delayed <- stblr_csr(
    method = "sbayesr",
    stats = fixture$stats,
    Glist = fixture$Glist,
    h2 = 0.3,
    adjE = 0.9,
    nit = 20,
    nburn = 5,
    ncores = 1,
    nchains = 1,
    seed = 10,
    updateE = TRUE,
    updateE_start = 5,
    updateE_every = 2
  )
  expect_bayesr_csr_conventions(fit_E_delayed)
  expect_equal(fit_E_delayed$input$updateE_start, 5L)
  expect_equal(fit_E_delayed$input$updateE_every, 2L)
  expect_equal(
    unname(fit_E_delayed$updateE_diagnostics[, "n_updateE"]),
    10,
    tolerance = 1e-12
  )

  expect_error(
    stblr_csr(
      method = "sbayesr",
      stats = fixture$stats,
      Glist = fixture$Glist,
      h2 = 0.3,
      adjE = 0.9,
      nit = 12,
      nburn = 4,
      ncores = 1,
      nchains = 2,
      seed = 10,
      updateE = TRUE,
      updateE_every = 0
    ),
    "updateE_every must be a positive integer scalar"
  )

  expect_false(exists(
    ".stblr_csr_bayesr_experimental",
    envir = asNamespace("sblr"),
    inherits = FALSE
  ))
})

test_that("CSR BayesR LD-swap runs and returns diagnostics", {
  skip_if_not(
    exists("stblr_cpg_omp_csr_bayesr", mode = "function"),
    "native CSR BayesR symbol is not loaded"
  )

  fixture <- make_small_bayesr_csr_stats()
  write_high_ld_csr_ld_fixture(fixture$Glist$sparseLD$prefix, fixture$m)

  fit <- stblr_csr(
    method = "sbayesr",
    stats = fixture$stats,
    Glist = fixture$Glist,
    h2 = 0.3,
    adjE = 0.9,
    pi = c(0.5, 0.5, 0, 0),
    alpha = c(10, 10, 1, 1),
    nit = 12,
    nburn = 4,
    ncores = 1,
    nchains = 1,
    seed = 30,
    updateB = FALSE,
    updateE = FALSE,
    updatePi = FALSE,
    updateLDswap = TRUE,
    ld_swap_prob = 1,
    ld_swap_r2 = 0.8,
    ld_swap_max_friends = 2,
    ld_swap_moves = 2
  )

  expect_bayesr_csr_conventions(fit)
  expect_bayesr_ld_swap_diagnostics(fit)
  expect_true(isTRUE(fit$input$updateLDswap))
  expect_equal(fit$input$ld_swap_moves, 2L)

  fit_bridge <- stblr_csr(
    stats = fixture$stats,
    Glist = fixture$Glist,
    method = "sbayesr",
    h2 = 0.3,
    adjE = 0.9,
    pi = c(0.5, 0.5, 0, 0),
    alpha = c(10, 10, 1, 1),
    nit = 12,
    nburn = 4,
    ncores = 1,
    nchains = 1,
    seed = 30,
    updateB = FALSE,
    updateE = FALSE,
    updatePi = FALSE,
    updateLDswap = TRUE,
    ld_swap_prob = 1,
    ld_swap_r2 = 0.8,
    ld_swap_max_friends = 2,
    ld_swap_moves = 2
  )
  expect_bayesr_csr_conventions(fit_bridge)
  expect_bayesr_ld_swap_diagnostics(fit_bridge)
  expect_equal(fit_bridge$input$method, "sbayesr")
})

test_that("CSR BayesR LD-swap diagnostics aggregate across chains", {
  skip_if_not(
    exists("stblr_cpg_omp_csr_bayesr", mode = "function"),
    "native CSR BayesR symbol is not loaded"
  )

  fixture <- make_small_bayesr_csr_stats()
  write_high_ld_csr_ld_fixture(fixture$Glist$sparseLD$prefix, fixture$m)

  fit <- stblr_csr(
    method = "sbayesr",
    stats = fixture$stats,
    Glist = fixture$Glist,
    h2 = 0.3,
    adjE = 0.9,
    pi = c(0.5, 0.5, 0, 0),
    alpha = c(10, 10, 1, 1),
    nit = 12,
    nburn = 4,
    ncores = 1,
    nchains = 2,
    keep_chains = TRUE,
    convergence = "none",
    seed = 31,
    updateB = FALSE,
    updateE = FALSE,
    updatePi = FALSE,
    updateLDswap = TRUE,
    ld_swap_prob = 1,
    ld_swap_r2 = 0.8,
    ld_swap_max_friends = 2,
    ld_swap_moves = 2
  )

  expect_bayesr_csr_conventions(fit)
  expect_bayesr_csr_chain_aggregation(fit)
  expect_bayesr_ld_swap_diagnostics(fit, require_chains = TRUE)
  expect_true("ld_swap_chains" %in% names(fit$diagnostics))
  expect_equal(sum(fit$diagnostics$ld_swap_chains$T1$attempted),
               fit$diagnostics$ld_swap$attempted)
  expect_true(all(vapply(fit$chains, function(ch)
    "ld_swap" %in% names(ch), logical(1))))
})
