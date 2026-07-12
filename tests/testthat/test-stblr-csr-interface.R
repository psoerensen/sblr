make_stblr_csr_interface_prefix <- function(m = 3L) {
  prefix <- tempfile("stblr_csr_interface_ld_")
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

make_stblr_csr_interface_stats <- function() {
  markers <- paste0("m", 1:3)
  list(
    wy = list(trait1 = stats::setNames(c(2, -1, 0.5), markers)),
    ww = list(trait1 = stats::setNames(rep(50, 3), markers)),
    yy = stats::setNames(50, "trait1"),
    n = 50L,
    m = 3L,
    marker_names = markers,
    trait_names = "trait1"
  )
}

make_stblr_csr_interface_raw <- function(nit = 3L, nburn = 1L,
                                         nchains = 1L,
                                         keep_chains = FALSE,
                                         updateLDswap = FALSE) {
  stats <- make_stblr_csr_interface_stats()
  nt <- length(stats$yy)
  m <- stats$m
  n <- stats$n
  h2 <- 0.3
  pi_init <- 0.5
  vy <- as.numeric(stats$yy) / (n - 1)
  B <- diag((vy * h2) / (m * pi_init), nt, nt)
  E <- diag(vy * (1 - h2), nt, nt)
  ssb_prior <- diag(((4 - 2) / 4) * (vy * h2) / (m * pi_init), nt, nt)
  sse_prior <- diag(((4 - 2) / 4) * (vy * (1 - h2)), nt, nt)

  stblr_cpg_omp_csr(
    wy = stats$wy,
    ww = stats$ww,
    yy = stats$yy,
    b_init = list(rep(0, m)),
    d_init = list(rep(0, m)),
    use_d_init = FALSE,
    r_init = stats$wy,
    use_r_init = FALSE,
    rebuild_r_before_updateE = FALSE,
    ld_prefix = make_stblr_csr_interface_prefix(),
    B = B,
    E = E,
    ssb_prior = split(ssb_prior, rep(seq_len(nt), each = nt)),
    sse_prior = split(sse_prior, rep(seq_len(nt), each = nt)),
    pi = c(1 - pi_init, pi_init),
    nub = 4,
    nue = 4,
    updateB = FALSE,
    updateE = FALSE,
    updatePi = FALSE,
    adjE = 0.9,
    n = rep(as.integer(n), nt),
    nit = as.integer(nit),
    nburn = as.integer(nburn),
    nthin = 1L,
    pi_prior_a = 1,
    pi_prior_b = 1,
    ncores = 1L,
    seed = 101L,
    nchains = as.integer(nchains),
    keep_chains = keep_chains,
    chain_seeds = integer(),
    updateLDswap = updateLDswap,
    ld_swap_prob = 1,
    ld_swap_r2 = 0.01,
    ld_swap_max_friends = 10L,
    ld_swap_moves = 1L
  )
}

make_stblr_csr_interface_bayesr_raw <- function(nit = 3L, nburn = 1L,
                                                nchains = 1L,
                                                keep_chains = FALSE,
                                                updateLDswap = FALSE) {
  stats <- make_stblr_csr_interface_stats()
  nt <- length(stats$yy)
  m <- stats$m
  n <- stats$n
  h2 <- 0.3
  mixture_var <- c(0, 0.1, 1)
  pi <- c(0.4, 0.3, 0.3)
  alpha <- c(1, 1, 1)
  vy <- as.numeric(stats$yy) / (n - 1)
  pi_active <- sum(pi[-1L])
  B <- diag((vy * h2) / (m * pi_active), nt, nt)
  E <- diag(vy * (1 - h2), nt, nt)
  ssb_prior <- diag(((4 - 2) / 4) * (vy * h2) / (m * pi_active), nt, nt)
  sse_prior <- diag(((4 - 2) / 4) * (vy * (1 - h2)), nt, nt)

  stblr_cpg_omp_csr_bayesr(
    wy = stats$wy,
    ww = stats$ww,
    yy = stats$yy,
    b_init = list(rep(0, m)),
    comp_init = list(rep(0, m)),
    use_comp_init = FALSE,
    r_init = stats$wy,
    use_r_init = FALSE,
    rebuild_r_before_updateE = FALSE,
    ld_prefix = make_stblr_csr_interface_prefix(),
    B = B,
    E = E,
    ssb_prior = split(ssb_prior, rep(seq_len(nt), each = nt)),
    sse_prior = split(sse_prior, rep(seq_len(nt), each = nt)),
    pi = pi,
    mixture_var = mixture_var,
    alpha = alpha,
    nub = 4,
    nue = 4,
    updateB = FALSE,
    updateE = FALSE,
    updatePi = FALSE,
    adjE = 0.9,
    n = rep(as.integer(n), nt),
    nit = as.integer(nit),
    nburn = as.integer(nburn),
    nthin = 1L,
    ncores = 1L,
    seed = 102L,
    nchains = as.integer(nchains),
    keep_chains = keep_chains,
    chain_seeds = integer(),
    updateLDswap = updateLDswap,
    ld_swap_prob = 1,
    ld_swap_r2 = 0.01,
    ld_swap_max_friends = 10L,
    ld_swap_moves = 1L
  )
}

make_stblr_csr_interface_sbayesrc_raw <- function(nit = 3L, nburn = 1L,
                                                  nchains = 1L,
                                                  keep_chains = FALSE,
                                                  updateLDswap = FALSE,
                                                  estimate_selection_s = FALSE) {
  stats <- make_stblr_csr_interface_stats()
  nt <- length(stats$yy)
  m <- stats$m
  n <- stats$n
  h2 <- 0.3
  gamma <- c(0, 0.1, 1)
  pi_init <- 0.5
  A <- matrix(
    c(1, 0, 1, 1, 1, 0),
    nrow = m,
    ncol = 2,
    dimnames = list(stats$marker_names, c("intercept", "annot1"))
  )
  alpha <- sblr::make_sbayesrc_alpha_init(A, gamma = gamma, pi_init = pi_init)
  vy <- as.numeric(stats$yy) / (n - 1)
  B <- diag((vy * h2) / (m * pi_init), nt, nt)
  E <- diag(vy * (1 - h2), nt, nt)
  ssb_prior <- diag(((4 - 2) / 4) * (vy * h2) / (m * pi_init), nt, nt)
  sse_prior <- diag(((4 - 2) / 4) * (vy * (1 - h2)), nt, nt)

  stblr_cpg_omp_csr_sbayesrc(
    wy = stats$wy,
    ww = stats$ww,
    yy = stats$yy,
    b_init = list(rep(0, m)),
    comp_init = list(rep(0, m)),
    use_comp_init = FALSE,
    r_init = stats$wy,
    use_r_init = FALSE,
    rebuild_r_before_updateE = FALSE,
    ld_prefix = make_stblr_csr_interface_prefix(),
    B = B,
    E = E,
    ssb_prior = split(ssb_prior, rep(seq_len(nt), each = nt)),
    sse_prior = split(sse_prior, rep(seq_len(nt), each = nt)),
    A = A,
    gamma = gamma,
    alpha_init = alpha$alpha_init,
    sigmaSqAlpha_init = alpha$sigmaSqAlpha_init,
    intercept_flat = TRUE,
    sigmaSqAlpha_a = 2,
    sigmaSqAlpha_b = 2,
    pi_floor = 1e-12,
    nub = 4,
    nue = 4,
    updateAlpha = FALSE,
    updateB = FALSE,
    updateE = FALSE,
    alpha_update_every = 10L,
    adjE = 0.9,
    n = rep(as.integer(n), nt),
    nit = as.integer(nit),
    nburn = as.integer(nburn),
    nthin = 1L,
    ncores = 1L,
    seed = 103L,
    nchains = as.integer(nchains),
    keep_chains = keep_chains,
    chain_seeds = integer(),
    updateLDswap = updateLDswap,
    ld_swap_prob = 1,
    ld_swap_r2 = 0.01,
    ld_swap_max_friends = 10L,
    ld_swap_moves = 1L,
    selection_s_prior_scale = numeric(),
    estimate_selection_s = estimate_selection_s,
    selection_s_init = 0,
    selection_s_prior = c(-3, 2),
    selection_s_proposal_sd = 0.25,
    selection_s_log_h = if (estimate_selection_s) rep(log(0.3), m) else numeric()
  )
}

test_that("ordinary CSR BayesC native return uses raw schema v1", {
  skip_if_not(
    exists("stblr_cpg_omp_csr", mode = "function"),
    "native BayesC CSR symbol is not loaded"
  )

  raw <- make_stblr_csr_interface_raw()
  expect_equal(raw$schema$class, "stblr_raw")
  expect_equal(as.integer(raw$schema$version), 1L)
  expect_true(inherits(raw, "stblr_raw"))
  expect_true(all(c(
    "schema", "meta", "marker", "trace", "variance", "pi", "diagnostics",
    "chains", "prior", "group", "annotation", "component", "selection"
  ) %in% names(raw)))
})

test_that("CSR BayesR native return uses raw schema v1 component namespaces", {
  skip_if_not(
    exists("stblr_cpg_omp_csr_bayesr", mode = "function"),
    "native BayesR CSR symbol is not loaded"
  )

  raw <- make_stblr_csr_interface_bayesr_raw()
  expect_equal(raw$schema$class, "stblr_raw")
  expect_equal(as.integer(raw$schema$version), 1L)
  expect_true(inherits(raw, "stblr_raw"))
  expect_named(raw, c(
    "schema", "meta", "marker", "trace", "variance", "pi", "diagnostics",
    "chains", "prior", "group", "annotation", "component", "selection"
  ))
  expect_equal(raw$meta$model, "bayesr")
  expect_equal(raw$meta$backend, "csr_bayesr")
  expect_true(all(c("names", "mixture_var", "prob", "ncomp", "dm_component_mean") %in%
                    names(raw$component)))
  expect_equal(raw$component$names[1L], "component_0")
  expect_equal(dim(raw$component$prob[[1L]]), c(raw$meta$m, raw$meta$n_components))
  expect_equal(rowSums(raw$component$prob[[1L]]), rep(1, raw$meta$m), tolerance = 1e-8)
  expect_equal(
    as.numeric(raw$marker$dm[, 1L]),
    1 - raw$component$prob[[1L]][, 1L],
    tolerance = 1e-8
  )
  expect_equal(dim(raw$trace$pis), c(raw$meta$n_trace, raw$meta$nt))
})

test_that("CSR SBayesRC native return uses raw schema v1 component and annotation namespaces", {
  skip_if_not(
    exists("stblr_cpg_omp_csr_sbayesrc", mode = "function"),
    "native SBayesRC CSR symbol is not loaded"
  )

  raw <- make_stblr_csr_interface_sbayesrc_raw()
  expect_equal(raw$schema$class, "stblr_raw")
  expect_equal(as.integer(raw$schema$version), 1L)
  expect_true(inherits(raw, "stblr_raw"))
  expect_named(raw, c(
    "schema", "meta", "marker", "trace", "variance", "pi", "diagnostics",
    "chains", "prior", "group", "annotation", "component", "selection"
  ))
  expect_equal(raw$meta$model, "sbayesrc")
  expect_true(all(c("names", "gamma", "mixture_var", "prob", "ncomp", "dm_component_mean") %in%
                    names(raw$component)))
  expect_true(all(c("alpha_mean", "alpha_final", "sigmaSqAlpha_mean", "sigmaSqAlpha_final") %in%
                    names(raw$annotation)))
  expect_equal(raw$component$names[1L], "gamma_0.00")
  expect_equal(dim(raw$component$prob[[1L]]), c(raw$meta$m, raw$meta$n_components))
  expect_equal(rowSums(raw$component$prob[[1L]]), rep(1, raw$meta$m), tolerance = 1e-8)
  expect_equal(
    as.numeric(raw$marker$dm[, 1L]),
    1 - raw$component$prob[[1L]][, "gamma_0.00"],
    tolerance = 1e-8
  )
  expect_equal(dim(raw$trace$pis), c(raw$meta$n_trace, raw$meta$nt))
})

test_that("CSR SBayesRC raw formatter preserves formatted fields and no-chain behavior", {
  skip_if_not(
    exists("stblr_cpg_omp_csr_sbayesrc", mode = "function"),
    "native SBayesRC CSR symbol is not loaded"
  )

  raw <- make_stblr_csr_interface_sbayesrc_raw(keep_chains = FALSE)
  raw$annotation$annotation_names <- c("intercept", "annot1")
  fit <- sblr:::.as_stblr_fit(
    raw,
    trait_names = make_stblr_csr_interface_stats()$trait_names,
    variable_names = make_stblr_csr_interface_stats()$marker_names
  )

  expect_true(all(c(
    "bm", "dm", "vbs", "vgs", "ves", "vle", "vld", "pi", "pim", "pis",
    "comp_prob", "dm_component_mean", "alpha", "sigmaSqAlpha"
  ) %in% names(fit)))
  expect_equal(dim(fit$bm), c(3L, 1L))
  expect_equal(dim(fit$dm), c(3L, 1L))
  for (nm in c("vbs", "vgs", "ves", "vle", "vld", "pis")) {
    expect_equal(dim(fit[[nm]]), c(raw$meta$n_trace, raw$meta$nt))
  }
  expect_equal(names(fit$comp_prob), "trait1")
  expect_equal(unname(rowSums(fit$comp_prob$trait1)), rep(1, 3), tolerance = 1e-8)
  expect_equal(
    unname(as.numeric(fit$dm[, "trait1"])),
    unname(1 - fit$comp_prob$trait1[, "gamma_0.00"]),
    tolerance = 1e-8
  )
  expect_null(fit$chains)
  expect_null(fit$ld_swap_chains)
})

test_that("ordinary CSR BayesC raw formatting tolerates no-chain schema objects", {
  skip_if_not(
    exists("stblr_cpg_omp_csr", mode = "function"),
    "native BayesC CSR symbol is not loaded"
  )

  raw <- make_stblr_csr_interface_raw(keep_chains = FALSE)
  expect_false(isTRUE(raw$meta$keep_chains))

  fit_null <- sblr:::.as_stblr_fit(
    raw,
    trait_names = make_stblr_csr_interface_stats()$trait_names,
    variable_names = make_stblr_csr_interface_stats()$marker_names
  )
  expect_null(fit_null$chains)
  expect_null(fit_null$ld_swap_chains)

  raw_empty <- raw
  raw_empty$chains <- list()
  fit_empty <- sblr:::.as_stblr_fit(
    raw_empty,
    trait_names = make_stblr_csr_interface_stats()$trait_names,
    variable_names = make_stblr_csr_interface_stats()$marker_names
  )
  expect_null(fit_empty$chains)
  expect_null(fit_empty$ld_swap_chains)

  raw_null_trait <- raw
  raw_null_trait$chains <- list(NULL)
  fit_null_trait <- sblr:::.as_stblr_fit(
    raw_null_trait,
    trait_names = make_stblr_csr_interface_stats()$trait_names,
    variable_names = make_stblr_csr_interface_stats()$marker_names
  )
  expect_null(fit_null_trait$chains)
  expect_null(fit_null_trait$ld_swap_chains)
})

test_that("stblr_csr fits BayesC through public method interface", {
  skip_if_not(
    exists("stblr_cpg_omp_csr", mode = "function"),
    "native BayesC CSR symbol is not loaded"
  )

  fit <- stblr_csr(
    stats = make_stblr_csr_interface_stats(),
    ld_prefix = make_stblr_csr_interface_prefix(),
    method = "bayesC",
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
    seed = 101L
  )

  required_fields <- c(
    "bm", "dm", "wy", "r", "b", "d",
    "vbs", "vgs", "ves", "vle", "vld", "pis",
    "covb", "covg", "cove",
    "pi", "pim",
    "rb", "rg", "re",
    "input"
  )
  missing_fields <- setdiff(required_fields, names(fit))
  expect_equal(missing_fields, character())
  expect_equal(dim(fit$bm), c(3L, 1L))
  expect_equal(dim(fit$dm), c(3L, 1L))
  for (nm in c("vbs", "vgs", "ves", "vle", "vld", "pis")) {
    expect_equal(dim(fit[[nm]]), c(2L, 1L))
  }
  expect_equal(fit$input$method, "bayesc")
  expect_equal(fit$input$model, "bayesc")
  expect_equal(fit$input$backend, "csr_bayesc")
  expect_equal(fit$input$data_level, "summary")
  expect_equal(fit$input$scheduled, FALSE)
  expect_equal(fit$input$nchains, 1L)
  expect_null(fit$chains)
  expect_null(fit$ld_swap_chains)
})

test_that("stblr_csr BayesC no-chain formatting keeps LD-swap diagnostics", {
  skip_if_not(
    exists("stblr_cpg_omp_csr", mode = "function"),
    "native BayesC CSR symbol is not loaded"
  )

  fit <- stblr_csr(
    stats = make_stblr_csr_interface_stats(),
    ld_prefix = make_stblr_csr_interface_prefix(),
    method = "bayesC",
    pi_init = 0.5,
    pi_prior_mean = 0.5,
    pi_prior_strength = 2,
    updateB = FALSE,
    updateE = FALSE,
    updatePi = FALSE,
    nit = 2,
    nburn = 0,
    nchains = 1L,
    keep_chains = FALSE,
    ncores = 1L,
    seed = 103L,
    updateLDswap = TRUE,
    ld_swap_prob = 1,
    ld_swap_r2 = 0.01,
    ld_swap_max_friends = 10L,
    ld_swap_moves = 1L
  )

  expect_null(fit$chains)
  expect_null(fit$ld_swap_chains)
  expect_s3_class(fit$ld_swap, "data.frame")
  expect_true(all(c("attempted", "accepted", "acceptance_rate") %in%
                    names(fit$ld_swap)))
})

test_that("stblr_csr fits BayesR through public method interface", {
  skip_if_not(
    exists("stblr_cpg_omp_csr_bayesr", mode = "function"),
    "native BayesR CSR symbol is not loaded"
  )

  fit <- stblr_csr(
    stats = make_stblr_csr_interface_stats(),
    ld_prefix = make_stblr_csr_interface_prefix(),
    method = "bayesR",
    updateB = FALSE,
    updateE = FALSE,
    updatePi = FALSE,
    nit = 2,
    nburn = 0,
    nchains = 1L,
    ncores = 1L,
    seed = 102L
  )

  expect_true(all(c(
    "bm", "dm", "vbs", "vgs", "ves", "vle", "vld", "pi", "pim", "pis",
    "input", "comp_prob", "dm_component_mean"
  ) %in% names(fit)))
  expect_equal(dim(fit$bm), c(3L, 1L))
  expect_equal(dim(fit$dm), c(3L, 1L))
  for (nm in c("vbs", "vgs", "ves", "vle", "vld", "pis")) {
    expect_equal(dim(fit[[nm]]), c(2L, 1L))
  }
  expect_equal(fit$input$method, "bayesr")
  expect_equal(fit$input$model, "bayesr")
  expect_equal(fit$input$backend, "csr_bayesr")
  expect_equal(fit$input$data_level, "summary")
  expect_equal(fit$input$scheduled, FALSE)
  expect_equal(fit$input$nchains, 1L)
  expect_equal(
    unname(as.numeric(fit$dm[, "trait1"])),
    unname(1 - fit$comp_prob$trait1[, "component_0"]),
    tolerance = 1e-8
  )
  expect_equal(unname(rowSums(fit$comp_prob$trait1)), rep(1, 3), tolerance = 1e-8)
  expect_null(fit$chains)
  expect_null(fit$ld_swap_chains)
})

test_that("stblr_csr_annot routes SBayesRC through raw v1 formatter", {
  skip_if_not(
    exists("stblr_cpg_omp_csr_sbayesrc", mode = "function"),
    "native SBayesRC CSR symbol is not loaded"
  )

  stats <- make_stblr_csr_interface_stats()
  annotations <- matrix(
    c(1, 0, 1, 1, 1, 0),
    nrow = stats$m,
    ncol = 2,
    dimnames = list(stats$marker_names, c("intercept", "annot1"))
  )
  glist <- list(
    sparseLD = list(prefix = make_stblr_csr_interface_prefix()),
    rsidsLD = list(stats$marker_names),
    rsids = list(stats$marker_names),
    maf = list(c(0.1, 0.2, 0.3))
  )

  fit <- stblr_csr_annot(
    Glist = glist,
    stats = stats,
    annotations = annotations,
    annotation_model = "sbayesrc",
    gamma = c(0, 0.1, 1),
    pi_init = 0.5,
    updateAlpha = FALSE,
    updateB = FALSE,
    updateE = FALSE,
    nit = 3,
    nburn = 1,
    nthin = 1,
    ncores = 1,
    seed = 104,
    nchains = 1L,
    keep_chains = FALSE
  )

  required_fields <- c(
    "bm", "dm", "wy", "r", "b", "d",
    "vbs", "vgs", "ves", "vle", "vld", "pis",
    "covb", "covg", "cove", "pi", "pim", "input",
    "comp_prob", "dm_component_mean", "alpha", "sigmaSqAlpha",
    "annotation_summary", "annotation_pi", "annotation_effects"
  )
  missing_fields <- setdiff(required_fields, names(fit))
  expect_equal(missing_fields, character())
  expect_equal(fit$input$model, "sbayesrc")
  expect_equal(fit$input$backend, "csr_sbayesrc")
  expect_equal(dim(fit$bm), c(3L, 1L))
  expect_equal(dim(fit$dm), c(3L, 1L))
  expect_identical(colnames(fit$comp_prob$trait1)[1L], "gamma_0.00")
  expect_equal(unname(rowSums(fit$comp_prob$trait1)), rep(1, 3), tolerance = 1e-8)
  expect_equal(
    unname(as.numeric(fit$dm[, "trait1"])),
    unname(1 - fit$comp_prob$trait1[, "gamma_0.00"]),
    tolerance = 1e-8
  )
  expect_equal(dim(fit$alpha$trait1), c(ncol(fit$input$A), 2L))
  expect_equal(dim(fit$sigmaSqAlpha), c(1L, 2L))
  expect_null(fit$chains)
  expect_null(fit$ld_swap_chains)
})

test_that("CSR SBayesRC raw v1 keep_chains exposes compact chain summaries", {
  skip_if_not(
    exists("stblr_cpg_omp_csr_sbayesrc", mode = "function"),
    "native SBayesRC CSR symbol is not loaded"
  )

  raw <- make_stblr_csr_interface_sbayesrc_raw(
    nchains = 2L,
    keep_chains = TRUE,
    updateLDswap = TRUE
  )
  raw$annotation$annotation_names <- c("intercept", "annot1")
  fit <- sblr:::.as_stblr_fit(
    raw,
    trait_names = make_stblr_csr_interface_stats()$trait_names,
    variable_names = make_stblr_csr_interface_stats()$marker_names
  )

  expect_equal(length(fit$chains), 1L)
  expect_equal(length(fit$chains[[1L]]), 2L)
  expect_s3_class(fit$ld_swap_chains[[1L]], "data.frame")
  expect_true(all(vapply(fit$chains[[1L]], function(ch) {
    is.numeric(ch$bm) &&
      is.numeric(ch$dm) &&
      is.matrix(ch$comp_prob) &&
      is.matrix(ch$alpha) &&
      is.numeric(ch$sigmaSqAlpha) &&
      all(abs(rowSums(ch$comp_prob) - 1) < 1e-8)
  }, logical(1))))
})
