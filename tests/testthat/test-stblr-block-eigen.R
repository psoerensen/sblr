make_stblr_block_eigen_bed <- function(path, dosage) {
  dosage_to_code <- c(`0` = 3L, `1` = 2L, `2` = 0L)
  packed <- unlist(lapply(seq_len(nrow(dosage)), function(marker) {
    codes <- unname(dosage_to_code[as.character(dosage[marker, ])])
    codes <- c(codes, rep(0L, (-length(codes)) %% 4L))
    vapply(seq(1L, length(codes), by = 4L), function(i) {
      sum(codes[i:(i + 3L)] * c(1L, 4L, 16L, 64L))
    }, integer(1))
  }))
  writeBin(as.raw(c(0x6c, 0x1b, 0x01, packed)), path)
}

make_stblr_block_eigen_fixture <- function() {
  bed_file <- tempfile(fileext = ".bed")
  make_stblr_block_eigen_bed(
    bed_file,
    rbind(c(0, 1, 2, 0, 1, 2), c(2, 0, 1, 2, 0, 1))
  )
  Glist <- list(
    n = 6L,
    ids = paste0("id", seq_len(6L)),
    bedfiles = bed_file,
    rsids = list(c("rs1", "rs2")),
    rsidsLD = list(c("rs1", "rs2")),
    chr = list(c(1L, 1L)),
    pos = list(c(100, 200)),
    af = list(c(0.2, 0.3)),
    maf = list(c(0.2, 0.3))
  )
  y <- matrix(
    c(-1, 0, 1, -0.5, 0.5, 1.5),
    ncol = 1,
    dimnames = list(Glist$ids, "D1")
  )
  stats <- list(
    wy = list(D1 = stats::setNames(c(1, -0.5), c("rs1", "rs2"))),
    ww = list(D1 = stats::setNames(c(6, 6), c("rs1", "rs2"))),
    yy = stats::setNames(sum(y[, 1]^2), "D1"),
    n = 6L,
    m = 2L,
    bed_files = bed_file,
    cls = list(1:2),
    rows = seq_len(6L),
    af = list(c(0.2, 0.3)),
    marker_names = c("rs1", "rs2"),
    trait_names = "D1"
  )
  list(Glist = Glist, stats = stats)
}

make_stblr_block_eigen_csr_prefix <- function(m = 2L) {
  prefix <- tempfile("stblr_block_eigen_csr_")
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

skip_if_no_block_eigen_native <- function() {
  ok <- tryCatch({
    getNativeSymbolInfo("_sblr_stblr_cpg_omp_csr_block_eigen", PACKAGE = "sblr")
    TRUE
  }, error = function(e) FALSE)
  skip_if_not(ok, "native block-eigen BayesC symbol is not loaded")
}

skip_if_no_csr_native <- function() {
  ok <- tryCatch({
    getNativeSymbolInfo("_sblr_stblr_cpg_omp_csr", PACKAGE = "sblr")
    TRUE
  }, error = function(e) FALSE)
  skip_if_not(ok, "native CSR BayesC symbol is not loaded")
}

skip_if_no_block_eigen_bayesr_native <- function() {
  ok <- tryCatch({
    getNativeSymbolInfo(
      "_sblr_stblr_cpg_omp_csr_bayesr_block_eigen",
      PACKAGE = "sblr"
    )
    TRUE
  }, error = function(e) FALSE)
  skip_if_not(ok, "native block-eigen BayesR symbol is not loaded")
}

skip_if_no_csr_bayesr_native <- function() {
  ok <- tryCatch({
    getNativeSymbolInfo("_sblr_stblr_cpg_omp_csr_bayesr", PACKAGE = "sblr")
    TRUE
  }, error = function(e) FALSE)
  skip_if_not(ok, "native CSR BayesR symbol is not loaded")
}

skip_if_no_block_eigen_sbayesrc_native <- function() {
  ok <- tryCatch({
    getNativeSymbolInfo(
      "_sblr_stblr_cpg_omp_csr_sbayesrc_block_eigen",
      PACKAGE = "sblr"
    )
    TRUE
  }, error = function(e) FALSE)
  skip_if_not(ok, "native block-eigen SBayesRC symbol is not loaded")
}

skip_if_no_csr_sbayesrc_native <- function() {
  ok <- tryCatch({
    getNativeSymbolInfo("_sblr_stblr_cpg_omp_csr_sbayesrc", PACKAGE = "sblr")
    TRUE
  }, error = function(e) FALSE)
  skip_if_not(ok, "native CSR SBayesRC symbol is not loaded")
}

expect_stblr_block_eigen_diagnostics <- function(fit, ridge = FALSE) {
  diagnostics <- fit$input$eigen_diagnostics
  if (is.null(diagnostics)) return(invisible(NULL))
  blocks <- if (is.data.frame(diagnostics)) diagnostics else diagnostics$blocks
  expect_true(is.data.frame(blocks) || is.list(blocks))
  expect_equal(
    setdiff(c("start", "size", "n_kept", "mu_min", "shrink"), names(blocks)),
    character()
  )
  expect_true(all(blocks$size > 0))
  expect_true(all(blocks$n_kept > 0))
  expect_true(all(blocks$shrink >= 0 & blocks$shrink <= 1))
  if (ridge) {
    expect_equal(blocks$n_kept, blocks$size)
  } else {
    expect_true(all(blocks$n_kept <= blocks$size))
  }
  invisible(NULL)
}

test_that("internal block-eigen BayesC helper exists but is not exported", {
  expect_true(is.function(sblr:::.stblr_csr_bayesc_block_eigen))
  expect_false(".stblr_csr_bayesc_block_eigen" %in% getNamespaceExports("sblr"))
  expect_false("stblr_cpg_omp_csr_block_eigen" %in% getNamespaceExports("sblr"))
})

test_that("block-eigen block_start validation accepts 0- and 1-based input", {
  fixture <- make_stblr_block_eigen_fixture()
  bed_one <- sblr:::.stblr_csr_block_eigen_inputs(
    fixture$stats,
    fixture$Glist,
    block_start = 1L
  )
  bed_zero <- sblr:::.stblr_csr_block_eigen_inputs(
    fixture$stats,
    fixture$Glist,
    block_start = 0L
  )
  expect_equal(bed_one$block_start, 0L)
  expect_equal(bed_zero$block_start, 0L)
  for (bad_start in list(2L, integer(), c(1L, 1L))) {
    expect_error(
      sblr:::.stblr_csr_block_eigen_inputs(
        fixture$stats,
        fixture$Glist,
        block_start = bad_start
      ),
      "block_start"
    )
  }
})

test_that("internal block-eigen BayesC hard-truncate run returns CSR BayesC fields", {
  skip_if_no_block_eigen_native()
  fixture <- make_stblr_block_eigen_fixture()
  fit <- sblr:::.stblr_csr_bayesc_block_eigen(
    stats = fixture$stats,
    Glist = fixture$Glist,
    block_start = 1L,
    eigen_filter = "hard_truncate",
    eigen_tau = 0.01,
    pi_init = 0.5,
    pi_prior_mean = 0.5,
    pi_prior_strength = 2,
    updateB = FALSE,
    updateE = FALSE,
    updatePi = FALSE,
    nit = 5,
    nburn = 2,
    nchains = 1L,
    keep_chains = FALSE,
    ncores = 1L,
    seed = 1L
  )

  required_fields <- c(
    "bm", "dm", "wy", "r", "b", "d",
    "vbs", "vgs", "ves", "vle", "vld", "pis",
    "covb", "covg", "cove",
    "pi", "pim",
    "rb", "rg", "re",
    "input"
  )
  expect_equal(setdiff(required_fields, names(fit)), character())
  expect_null(fit$chains)
  expect_null(fit$ld_swap_chains)
  expect_equal(fit$input$ld_backend, "block_eigen")
  expect_equal(fit$input$eigen_filter, "hard_truncate")
  expect_true(all(is.finite(fit$bm)))
  expect_true(all(is.finite(fit$dm)))
  expect_true(all(is.finite(fit$vbs)))
  expect_true(all(is.finite(fit$ves)))
  expect_stblr_block_eigen_diagnostics(fit)
})

test_that("internal block-eigen BayesC ridge filters run on the tiny fixture", {
  skip_if_no_block_eigen_native()
  fixture <- make_stblr_block_eigen_fixture()
  common <- list(
    stats = fixture$stats, Glist = fixture$Glist, block_start = 1L,
    pi_init = 0.5, pi_prior_mean = 0.5, pi_prior_strength = 2,
    updateB = FALSE, updateE = FALSE, updatePi = FALSE,
    nit = 5, nburn = 2, nchains = 1L, keep_chains = FALSE,
    ncores = 1L, seed = 1L
  )
  fit_rg0 <- do.call(
    sblr:::.stblr_csr_bayesc_block_eigen,
    c(common, list(eigen_filter = "ridge_fixed", eigen_eta = 0))
  )
  expect_identical(fit_rg0$input$ld_backend, "block_eigen")
  expect_identical(fit_rg0$input$eigen_filter, "ridge_fixed")
  expect_equal(fit_rg0$input$eigen_eta, 0)
  expect_stblr_block_eigen_diagnostics(fit_rg0, ridge = TRUE)

  fit_lw <- do.call(
    sblr:::.stblr_csr_bayesc_block_eigen,
    c(common, list(eigen_filter = "ridge_lw"))
  )
  expect_identical(fit_lw$input$ld_backend, "block_eigen")
  expect_identical(fit_lw$input$eigen_filter, "ridge_lw")
  expect_stblr_block_eigen_diagnostics(fit_lw, ridge = TRUE)
})

test_that("block-eigen filter arguments are validated on the R side", {
  fixture <- make_stblr_block_eigen_fixture()
  common <- list(
    stats = fixture$stats, Glist = fixture$Glist, block_start = 1L,
    nit = 5, nburn = 2, seed = 1L
  )
  expect_error(
    do.call(sblr:::.stblr_csr_bayesc_block_eigen,
            c(common, list(eigen_filter = "bad_filter"))),
    "eigen_filter"
  )
  expect_error(
    do.call(sblr:::.stblr_csr_bayesc_block_eigen,
            c(common, list(eigen_tau = -1))),
    "eigen_tau"
  )
  expect_error(
    do.call(sblr:::.stblr_csr_bayesc_block_eigen,
            c(common, list(eigen_eta = -1))),
    "eigen_eta"
  )
})

test_that("block-eigen BayesC rejects LD-swap clearly", {
  fixture <- make_stblr_block_eigen_fixture()
  expect_error(
    sblr:::.stblr_csr_bayesc_block_eigen(
      stats = fixture$stats,
      Glist = fixture$Glist,
      block_start = 1L,
      updateLDswap = TRUE,
      nit = 2,
      nburn = 0
    ),
    "LD-swap is not yet supported with the experimental block-eigen operator"
  )
})

test_that("default stblr_csr BayesC remains on sparse CSR path", {
  skip_if_no_csr_native()
  fixture <- make_stblr_block_eigen_fixture()
  fit <- stblr_csr(
    stats = fixture$stats,
    ld_prefix = make_stblr_block_eigen_csr_prefix(m = fixture$stats$m),
    method = "sbayesc",
    pi_init = 0.5,
    pi_prior_mean = 0.5,
    pi_prior_strength = 2,
    updateB = FALSE,
    updateE = FALSE,
    updatePi = FALSE,
    nit = 2,
    nburn = 0,
    nchains = 1L,
    ncores = 1L,
    seed = 102L
  )
  expect_equal(fit$input$backend, "csr_bayesc")
  expect_null(fit$input$ld_backend)
})

test_that("internal block-eigen BayesR helper exists but is not exported", {
  expect_true(exists(
    ".stblr_csr_bayesr_block_eigen",
    envir = asNamespace("sblr"),
    inherits = FALSE
  ))
  expect_true(is.function(sblr:::.stblr_csr_bayesr_block_eigen))
  expect_false(".stblr_csr_bayesr_block_eigen" %in% getNamespaceExports("sblr"))
  expect_false(
    "stblr_cpg_omp_csr_bayesr_block_eigen" %in% getNamespaceExports("sblr")
  )
})

expect_stblr_block_eigen_bayesr_fit <- function(fit, eigen_filter) {
  required_fields <- c(
    "bm", "dm", "wy", "r", "b", "d",
    "vbs", "vgs", "ves", "vle", "vld", "pis",
    "covb", "covg", "cove", "pi", "pim", "input",
    "comp_prob", "dm_component_mean"
  )
  missing_fields <- setdiff(required_fields, names(fit))
  expect_equal(missing_fields, character())
  expect_null(fit$chains)
  expect_null(fit$ld_swap_chains)
  expect_identical(fit$input$ld_backend, "block_eigen")
  expect_identical(fit$input$eigen_filter, eigen_filter)
  expect_true(all(is.finite(fit$bm)))
  expect_true(all(is.finite(fit$dm)))
  expect_true(all(is.finite(fit$vbs)))
  expect_true(all(is.finite(fit$ves)))
  for (trait in names(fit$comp_prob)) {
    cp <- fit$comp_prob[[trait]]
    expect_equal(unname(rowSums(cp)), rep(1, nrow(cp)), tolerance = 1e-8)
    expect_true("component_0" %in% colnames(cp))
    expect_equal(
      unname(fit$dm[, trait]),
      unname(1 - cp[, "component_0"]),
      tolerance = 1e-8
    )
  }
  expect_stblr_block_eigen_diagnostics(
    fit,
    ridge = eigen_filter != "hard_truncate"
  )
}

test_that("internal block-eigen BayesR filter modes run on the tiny fixture", {
  skip_if_no_block_eigen_bayesr_native()
  fixture <- make_stblr_block_eigen_fixture()
  common <- list(
    stats = fixture$stats, Glist = fixture$Glist, block_start = 1L,
    updateB = FALSE, updateE = FALSE, updatePi = FALSE,
    nit = 5, nburn = 2, nchains = 1L, keep_chains = FALSE,
    ncores = 1L, seed = 1L
  )
  fits <- list(
    hard_truncate = do.call(
      sblr:::.stblr_csr_bayesr_block_eigen,
      c(common, list(eigen_filter = "hard_truncate", eigen_tau = 0.01))
    ),
    ridge_fixed = do.call(
      sblr:::.stblr_csr_bayesr_block_eigen,
      c(common, list(eigen_filter = "ridge_fixed", eigen_eta = 0))
    ),
    ridge_lw = do.call(
      sblr:::.stblr_csr_bayesr_block_eigen,
      c(common, list(eigen_filter = "ridge_lw"))
    )
  )
  for (filter in names(fits)) {
    expect_stblr_block_eigen_bayesr_fit(fits[[filter]], filter)
  }
  expect_equal(fits$ridge_fixed$input$eigen_eta, 0)
})

test_that("block-eigen BayesR validates filter arguments on the R side", {
  fixture <- make_stblr_block_eigen_fixture()
  common <- list(
    stats = fixture$stats,
    Glist = fixture$Glist,
    block_start = 1L,
    nit = 5,
    nburn = 2,
    seed = 1L
  )
  expect_error(
    do.call(
      sblr:::.stblr_csr_bayesr_block_eigen,
      c(common, list(eigen_filter = "bad_filter"))
    ),
    "eigen_filter"
  )
  expect_error(
    do.call(
      sblr:::.stblr_csr_bayesr_block_eigen,
      c(common, list(eigen_tau = -1))
    ),
    "eigen_tau"
  )
  expect_error(
    do.call(
      sblr:::.stblr_csr_bayesr_block_eigen,
      c(common, list(eigen_eta = -1))
    ),
    "eigen_eta"
  )
})

test_that("block-eigen BayesR rejects LD-swap clearly", {
  fixture <- make_stblr_block_eigen_fixture()
  expect_error(
    sblr:::.stblr_csr_bayesr_block_eigen(
      stats = fixture$stats,
      Glist = fixture$Glist,
      block_start = 1L,
      updateLDswap = TRUE,
      nit = 2,
      nburn = 0
    ),
    "LD-swap is not yet supported with the experimental block-eigen operator"
  )
})

test_that("default stblr_csr BayesR remains on sparse CSR path", {
  skip_if_no_csr_bayesr_native()
  fixture <- make_stblr_block_eigen_fixture()
  fit <- stblr_csr(
    stats = fixture$stats,
    ld_prefix = make_stblr_block_eigen_csr_prefix(m = fixture$stats$m),
    method = "sbayesr",
    updateB = FALSE,
    updateE = FALSE,
    updatePi = FALSE,
    nit = 2,
    nburn = 0,
    nchains = 1L,
    ncores = 1L,
    seed = 103L
  )
  expect_identical(fit$input$backend, "csr_bayesr")
  expect_false(identical(fit$input$ld_backend, "block_eigen"))
  expect_null(fit$input$ld_backend)
})

make_stblr_block_eigen_annotation <- function(fixture) {
  matrix(
    c(1, 1, 0, 1),
    nrow = fixture$stats$m,
    dimnames = list(fixture$stats$marker_names, c("intercept", "annot1"))
  )
}

test_that("internal block-eigen SBayesRC helper exists but is not exported", {
  expect_true(is.function(sblr:::.stblr_csr_sbayesrc_block_eigen))
  expect_false(".stblr_csr_sbayesrc_block_eigen" %in% getNamespaceExports("sblr"))
  expect_false(
    "stblr_cpg_omp_csr_sbayesrc_block_eigen" %in% getNamespaceExports("sblr")
  )
})

expect_stblr_block_eigen_sbayesrc_fit <- function(fit, eigen_filter, n_anno) {
  required_fields <- c(
    "bm", "dm", "wy", "r", "b", "d",
    "vbs", "vgs", "ves", "vle", "vld", "pis",
    "covb", "covg", "cove", "pi", "pim", "input",
    "comp_prob", "dm_component_mean", "alpha", "sigmaSqAlpha",
    "annotation_summary", "annotation_pi", "annotation_effects"
  )
  expect_equal(setdiff(required_fields, names(fit)), character())
  expect_null(fit$chains)
  expect_null(fit$ld_swap_chains)
  expect_identical(fit$input$ld_backend, "block_eigen")
  expect_identical(fit$input$eigen_filter, eigen_filter)
  expect_true(all(is.finite(fit$bm)))
  expect_true(all(is.finite(fit$dm)))
  expect_true(all(is.finite(fit$vbs)))
  expect_true(all(is.finite(fit$ves)))
  for (trait in names(fit$component_probabilities)) {
    cp <- fit$component_probabilities[[trait]]
    expect_identical(colnames(cp)[1L], "gamma_0.00")
    expect_equal(unname(rowSums(cp)), rep(1, nrow(cp)), tolerance = 1e-8)
    expect_equal(
      unname(fit$dm[, trait]),
      unname(1 - cp[, "gamma_0.00"]),
      tolerance = 1e-8
    )
    expect_equal(nrow(fit$alpha[[trait]]), n_anno)
    expect_equal(ncol(fit$alpha[[trait]]), ncol(cp) - 1L)
  }
  expect_equal(dim(fit$sigmaSqAlpha), c(length(fit$comp_prob), 3L))
  expect_stblr_block_eigen_diagnostics(
    fit,
    ridge = eigen_filter != "hard_truncate"
  )
}

test_that("internal block-eigen SBayesRC filter modes run on the tiny fixture", {
  skip_if_no_block_eigen_sbayesrc_native()
  fixture <- make_stblr_block_eigen_fixture()
  annotation <- make_stblr_block_eigen_annotation(fixture)
  common <- list(
    stats = fixture$stats, Glist = fixture$Glist,
    annotation = annotation, block_start = 1L,
    updateAlpha = FALSE, updateB = FALSE, updateE = FALSE,
    nit = 5, nburn = 2, ncores = 1L, seed = 1L,
    keep_chains = FALSE
  )
  fits <- list(
    hard_truncate = do.call(
      sblr:::.stblr_csr_sbayesrc_block_eigen,
      c(common, list(eigen_filter = "hard_truncate", eigen_tau = 0.01))
    ),
    ridge_fixed = do.call(
      sblr:::.stblr_csr_sbayesrc_block_eigen,
      c(common, list(eigen_filter = "ridge_fixed", eigen_eta = 0))
    ),
    ridge_lw = do.call(
      sblr:::.stblr_csr_sbayesrc_block_eigen,
      c(common, list(eigen_filter = "ridge_lw"))
    )
  )
  for (filter in names(fits)) {
    expect_stblr_block_eigen_sbayesrc_fit(fits[[filter]], filter, ncol(annotation))
  }
  expect_equal(fits$ridge_fixed$input$eigen_eta, 0)

  fit_zero <- do.call(
    sblr:::.stblr_csr_sbayesrc_block_eigen,
    c(
      common[names(common) != "block_start"],
      list(block_start = 0L, eigen_filter = "hard_truncate")
    )
  )
  expect_identical(fit_zero$input$eigen_blocks, 0L)
})

test_that("block-eigen SBayesRC validates block and filter arguments", {
  fixture <- make_stblr_block_eigen_fixture()
  common <- list(
    stats = fixture$stats,
    Glist = fixture$Glist,
    annotation = make_stblr_block_eigen_annotation(fixture),
    block_start = 1L
  )
  for (bad_start in list(integer(), 2L, c(1L, 1L))) {
    args <- common
    args$block_start <- bad_start
    expect_error(
      do.call(sblr:::.stblr_csr_sbayesrc_block_eigen, args),
      "block_start"
    )
  }
  expect_error(
    do.call(
      sblr:::.stblr_csr_sbayesrc_block_eigen,
      c(common, list(eigen_filter = "bad_filter"))
    ),
    "eigen_filter"
  )
  expect_error(
    do.call(
      sblr:::.stblr_csr_sbayesrc_block_eigen,
      c(common, list(eigen_tau = -1))
    ),
    "eigen_tau"
  )
  expect_error(
    do.call(
      sblr:::.stblr_csr_sbayesrc_block_eigen,
      c(common, list(eigen_eta = -1))
    ),
    "eigen_eta"
  )
})

test_that("sampled selection_s works with block-eigen SBayesRC", {
  skip_if_no_block_eigen_sbayesrc_native()
  fixture <- make_stblr_block_eigen_fixture()
  fit <- sblr:::.stblr_csr_sbayesrc_block_eigen(
    stats = fixture$stats,
    Glist = fixture$Glist,
    annotation = make_stblr_block_eigen_annotation(fixture),
    block_start = 1L,
    estimate_selection_s = TRUE,
    selection_s_init = 0,
    selection_s_proposal_sd = 0.25,
    updateAlpha = FALSE,
    updateB = FALSE,
    updateE = FALSE,
    nit = 5,
    nburn = 2,
    ncores = 1L,
    seed = 105L
  )
  expect_true("selection_s" %in% names(fit))
  expect_true("selection_s_trace" %in% names(fit))
  expect_true(all(is.finite(fit$selection_s)))
  expect_true(all(is.finite(fit$selection_s_trace)))
})

test_that("block-eigen SBayesRC rejects LD-swap clearly", {
  fixture <- make_stblr_block_eigen_fixture()
  expect_error(
    sblr:::.stblr_csr_sbayesrc_block_eigen(
      stats = fixture$stats,
      Glist = fixture$Glist,
      annotation = make_stblr_block_eigen_annotation(fixture),
      block_start = 1L,
      updateLDswap = TRUE,
      nit = 2,
      nburn = 0
    ),
    "LD-swap is not yet supported with the experimental block-eigen operator"
  )
})

test_that("public SBayesRC remains on the sparse CSR path", {
  skip_if_no_csr_sbayesrc_native()
  fixture <- make_stblr_block_eigen_fixture()
  fixture$Glist$sparseLD <- list(
    prefix = make_stblr_block_eigen_csr_prefix(fixture$stats$m)
  )
  fit <- stblr_csr_annot(
    stats = fixture$stats,
    Glist = fixture$Glist,
    annotations = make_stblr_block_eigen_annotation(fixture),
    annotation_model = "annotation_probit_stick",
    updateAlpha = FALSE,
    updateB = FALSE,
    updateE = FALSE,
    nit = 2,
    nburn = 0,
    ncores = 1L,
    seed = 104L
  )
  expect_identical(fit$input$backend, "csr_sbayesrc")
  expect_false(identical(fit$input$ld_backend, "block_eigen"))
})
