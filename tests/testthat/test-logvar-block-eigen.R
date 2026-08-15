make_logvar_block_gate_bed <- function(path, dosage) {
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

make_logvar_block_gate_fixture <- function() {
  bed_file <- tempfile(fileext = ".bed")
  make_logvar_block_gate_bed(
    bed_file,
    rbind(c(0, 1, 2, 0, 1, 2), c(2, 0, 1, 2, 0, 1))
  )
  glist <- list(
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
  y <- c(-1, 0, 1, -0.5, 0.5, 1.5)
  stats <- list(
    wy = list(D1 = stats::setNames(c(1, -0.5), c("rs1", "rs2"))),
    ww = list(D1 = stats::setNames(c(6, 6), c("rs1", "rs2"))),
    yy = stats::setNames(sum(y^2), "D1"),
    n = 6L,
    m = 2L,
    bed_files = bed_file,
    cls = list(1:2),
    rows = seq_len(6L),
    af = list(c(0.2, 0.3)),
    marker_names = c("rs1", "rs2"),
    trait_names = "D1"
  )
  list(Glist = glist, stats = stats)
}

logvar_block_gate_common <- function(fixture) {
  list(
    stats = fixture$stats,
    Glist = fixture$Glist,
    block_start = 1L,
    representation = "low_rank",
    eigen_prop = 0.999999,
    low_rank_residual_rebuild_every = 2L,
    nit = 10L,
    nburn = 4L,
    nthin = 1L,
    seed = 703L,
    nchains = 2L,
    ncores = 1L,
    chain_seeds = c(1703L, 2703L),
    keep_chains = TRUE,
    convergence = "extended"
  )
}

logvar_block_gate_trajectory <- function(fit, bayesr = FALSE) {
  out <- list(
    bm = fit$bm, dm = fit$dm, b = fit$b, d = fit$d,
    vbs = fit$vbs, vgs = fit$vgs, ves = fit$ves,
    vle = fit$vle, vld = fit$vld, pis = fit$pis,
    pi = fit$pi, pim = fit$pim, r = fit$r,
    chains = fit$chains,
    convergence = fit$convergence,
    convergence_traces = fit$convergence_traces,
    block_ve = fit$block_ve,
    block_ve_chains = fit$block_ve_chains
  )
  if (bayesr) {
    out <- c(out, list(
      comp_prob = fit$component_prob,
      dm_component_mean = fit$dm_component_mean,
      ncomp = fit$ncomp
    ))
  }
  out
}

logvar_block_gate_hash <- function(value) {
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path), add = TRUE)
  saveRDS(value, path, version = 2)
  unname(tools::md5sum(path))
}

test_that("ordinary retained block BayesC trajectory remains frozen", {
  fixture <- make_logvar_block_gate_fixture()
  fit <- blr_with_legacy_execution(function() {
    do.call(stblr_block_eigen, c(
      logvar_block_gate_common(fixture),
      list(
      method = "sbayesc",
      pi_init = 0.35,
      pi_prior_mean = 0.35,
      pi_prior_strength = 8,
      h2 = 0.4,
      convergence_control = list(
        warn = FALSE,
        extended_groups = "probability",
        selected_markers = 1:2,
        selected_marker_quantities = c("b", "d"),
        keep_traces = TRUE
      )
      )
    ))
  })
  expect_identical(
    logvar_block_gate_hash(logvar_block_gate_trajectory(fit)),
    "0e54c7ec18d4183c1351f657efcfeab2"
  )
})

test_that("ordinary retained block BayesR trajectory remains frozen", {
  fixture <- make_logvar_block_gate_fixture()
  fit <- blr_with_legacy_execution(function() {
    do.call(stblr_block_eigen, c(
      logvar_block_gate_common(fixture),
      list(
      method = "sbayesr",
      mixture_var = c(0, 0.05, 0.2),
      pi = c(0.6, 0.25, 0.15),
      alpha = c(3, 2, 1),
      h2 = 0.4,
      residual_policy = "gctb_block",
      block_ve_mode = "allMixVe",
      block_ve_keep_history = TRUE,
      convergence_control = list(
        warn = FALSE,
        extended_groups = "probability",
        selected_markers = 1:2,
        selected_marker_quantities = c("b", "d", "component"),
        keep_traces = TRUE
      )
      )
    ))
  })
  expect_identical(
    logvar_block_gate_hash(logvar_block_gate_trajectory(fit, bayesr = TRUE)),
    "aa329728ea6a3defd6d61f37ca1cea59"
  )
})

logvar_block_compare_fields <- function(fit) {
  list(
    bm = fit$bm, dm = fit$dm, b = fit$b, d = fit$d,
    vbs = fit$vbs, vgs = fit$vgs, ves = fit$ves,
    vle = fit$vle, vld = fit$vld, chains = fit$chains,
    component_prob = fit$component_prob,
    dm_component_mean = fit$dm_component_mean,
    ncomp = fit$ncomp, block_ve = fit$block_ve,
    block_ve_chains = fit$block_ve_chains
  )
}

test_that("zero-theta retained block BayesC-LV is an ordinary trajectory", {
  fixture <- blr_unified_fixture()
  on.exit(blr_unified_cleanup(fixture), add = TRUE)
  common <- list(
    stats = fixture$stats, Glist = fixture$Glist, block_start = 1L,
    method = "sbayesc", representation = "low_rank", eigen_prop = 0.999999,
    pi_init = 0.35, pi_prior_mean = 0.35, pi_prior_strength = 8,
    nit = 12L, nburn = 4L, seed = 811L, nchains = 2L, ncores = 1L,
    chain_seeds = c(1811L, 2811L), keep_chains = TRUE,
    convergence = "none"
  )
  ordinary <- blr_with_legacy_execution(function() {
    do.call(stblr_block_eigen, common)
  })
  lv <- do.call(stblr_block_eigen, c(common, list(
    annotations = matrix(c(0, 1, 1), ncol = 1,
      dimnames = list(fixture$stats$marker_names, "binary")),
    annotation_model = "log_variance", theta_init = 0,
    updateTheta = FALSE
  )))
  expect_identical(logvar_block_compare_fields(lv),
                   logvar_block_compare_fields(ordinary))
  expect_identical(drop(lv$theta), 0)
  expect_identical(unname(drop(lv$marker_prior_scale)), rep(1, 3))
})

test_that("zero-theta retained block BayesR-LV is an ordinary trajectory", {
  fixture <- blr_unified_fixture()
  on.exit(blr_unified_cleanup(fixture), add = TRUE)
  common <- list(
    stats = fixture$stats, Glist = fixture$Glist, block_start = 1L,
    method = "sbayesr", representation = "low_rank", eigen_prop = 0.999999,
    mixture_var = c(0, 0.05, 0.2), pi = c(0.6, 0.25, 0.15),
    alpha = c(3, 2, 1), nit = 12L, nburn = 4L, seed = 812L,
    nchains = 2L, ncores = 1L, chain_seeds = c(1812L, 2812L),
    keep_chains = TRUE, convergence = "none"
  )
  ordinary <- blr_with_legacy_execution(function() {
    do.call(stblr_block_eigen, common)
  })
  lv <- do.call(stblr_block_eigen, c(common, list(
    annotations = matrix(c(0, 1, 1), ncol = 1,
      dimnames = list(fixture$stats$marker_names, "binary")),
    annotation_model = "log_variance", theta_init = 0,
    updateTheta = FALSE
  )))
  expect_identical(logvar_block_compare_fields(lv),
                   logvar_block_compare_fields(ordinary))
  expect_identical(drop(lv$theta), 0)
  expect_identical(unname(drop(lv$marker_prior_scale)), rep(1, 3))
})

test_that("public retained block LV models expose theta and q diagnostics", {
  fixture <- blr_unified_fixture()
  on.exit(blr_unified_cleanup(fixture), add = TRUE)
  annotations <- cbind(
    binary = c(0, 1, 1),
    continuous = c(-1, 0.25, 2)
  )
  rownames(annotations) <- fixture$stats$marker_names
  common <- list(
    stats = fixture$stats, Glist = fixture$Glist, block_start = 1L,
    annotations = annotations, annotation_model = "log_variance",
    representation = "low_rank", eigen_prop = 0.9,
    nit = 24L, nburn = 8L, seed = 813L, nchains = 2L, ncores = 1L,
    chain_seeds = c(1813L, 2813L), keep_chains = TRUE,
    convergence = "core", convergence_control = list(warn = FALSE)
  )
  fits <- list(
    sbayesc = do.call(stblr_block_eigen, c(common, list(method = "sbayesc"))),
    sbayesr = do.call(stblr_block_eigen, c(common, list(
      method = "sbayesr", mixture_var = c(0, 0.05, 0.2),
      pi = c(0.6, 0.25, 0.15), alpha = c(3, 2, 1))))
  )
  for (method in names(fits)) {
    fit <- fits[[method]]
    expect_identical(fit$model, paste0(method, "_logvar"))
    expect_identical(fit$annotation_model, "log_variance")
    expect_identical(fit$input$ld_backend, "block_eigen")
    expect_equal(fit$input$theta_prior_sd, 0.7)
    expect_identical(fit$annotation_transform$type,
                     c("binary", "continuous"))
    expect_true(all(is.finite(fit$theta)))
    expect_true(all(is.finite(fit$marker_prior_scale)))
    expect_true(all(fit$marker_prior_scale > 0))
    expect_true(all(c("Rhat", "bulk_ESS", "tail_ESS", "MCSE") %in%
                      names(fit$theta_summary)))
    expect_gt(fit$diagnostics$logvar$theta_updates, 0)
    expect_true(all(is.finite(unlist(fit$diagnostics$logvar))))
    expect_true(!is.null(fit$theta_trace))
  }
})

test_that("effectively exact CSR and block LV paths are concordant", {
  fixture <- blr_unified_fixture()
  prefix <- blr_unified_exact_ld_prefix(fixture$dosage)
  on.exit(blr_unified_cleanup(fixture), add = TRUE)
  on.exit(blr_unified_cleanup_prefix(prefix), add = TRUE)
  annotations <- matrix(
    c(0, 1, 1), ncol = 1,
    dimnames = list(fixture$stats$marker_names, "binary")
  )
  common <- list(
    stats = fixture$stats, Glist = fixture$Glist,
    annotations = annotations, annotation_model = "log_variance",
    theta_init = 0, updateTheta = TRUE, nit = 60L, nburn = 20L,
    nchains = 2L, ncores = 1L, chain_seeds = c(1902L, 2902L),
    keep_chains = TRUE, convergence = "core",
    convergence_control = list(warn = FALSE), h2 = 0.4, adjE = 0
  )
  model_args <- list(
    sbayesc = list(
      pi_init = 0.35, pi_prior_mean = 0.35, pi_prior_strength = 8),
    sbayesr = list(
      mixture_var = c(0, 0.05, 0.2), pi = c(0.6, 0.25, 0.15),
      alpha = c(3, 2, 1))
  )
  for (method in names(model_args)) {
    csr <- do.call(stblr_csr_annot, c(
      common, list(method = method, ld_prefix = prefix), model_args[[method]]))
    block <- do.call(stblr_block_eigen, c(
      common,
      list(
        method = method, block_start = 1L,
        representation = "dense_reconstructed", eigen_tau = 0,
        residual_policy = "global_projected_legacy"
      ),
      model_args[[method]]
    ))
    expect_identical(block$theta, csr$theta, info = method)
    expect_identical(block$marker_prior_scale, csr$marker_prior_scale,
                     info = method)
    expect_identical(block$theta_chain_mean, csr$theta_chain_mean,
                     info = method)
    for (field in c("bm", "dm", "vbs", "vgs", "ves", "vle", "vld",
                    "pi_final", "pi_mean")) {
      expect_equal(block[[field]], csr[[field]], tolerance = 5e-7,
                   info = paste(method, field))
    }
    expect_true(all(is.finite(block$theta_summary$Rhat)))
    expect_true(all(block$theta_summary$bulk_ESS > 0))
    expect_true(all(block$theta_summary$tail_ESS > 0))
    if (identical(method, "sbayesr")) {
      expect_identical(block$component_prob, csr$component_prob)
      expect_identical(block$dm_component_mean, csr$dm_component_mean)
    }
  }
})

test_that("public block log-variance annotations fail before native dispatch", {
  fixture <- blr_unified_fixture()
  on.exit(blr_unified_cleanup(fixture), add = TRUE)
  intercept <- matrix(
    1, nrow = 3L, ncol = 1L,
    dimnames = list(fixture$stats$marker_names, "intercept")
  )
  expect_error(
    stblr_block_eigen(
      fixture$stats, fixture$Glist, 1L, method = "sbayesc",
      annotations = intercept, annotation_model = "log_variance",
      nit = 2L, nburn = 1L),
    "all-ones intercept"
  )
  rank_deficient <- cbind(a = c(0, 1, 0), duplicate = c(0, 1, 0))
  rownames(rank_deficient) <- fixture$stats$marker_names
  expect_error(
    stblr_block_eigen(
      fixture$stats, fixture$Glist, 1L, method = "sbayesr",
      annotations = rank_deficient, annotation_model = "log_variance",
      nit = 2L, nburn = 1L),
    "rank deficient|duplicate"
  )
})
