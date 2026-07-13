make_stblr_bed_interface_glist <- function(bed_file, n = 6L) {
  list(
    n = n,
    ids = paste0("id", seq_len(n)),
    bedfiles = bed_file,
    rsids = list(c("rs1", "rs2")),
    rsidsLD = list(c("rs1", "rs2")),
    chr = list(c(1L, 1L)),
    pos = list(c(100, 200)),
    af = list(c(0.2, 0.3))
  )
}

write_stblr_bed_interface_file <- function(path, dosage) {
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

make_stblr_bed_interface_fixture <- function() {
  bed_file <- tempfile(fileext = ".bed")
  write_stblr_bed_interface_file(
    bed_file,
    rbind(c(0, 1, 2, 0, 1, 2), c(2, 0, 1, 2, 0, 1))
  )
  list(
    bed_file = bed_file,
    Glist = make_stblr_bed_interface_glist(bed_file),
    y = matrix(
      c(-1, 0, 1, -0.5, 0.5, 1.5),
      ncol = 1,
      dimnames = list(NULL, "D1")
    )
  )
}

expect_stblr_bed_chain_fields <- function(fit) {
  for (nm in c("dm_sd", "dm_min", "dm_max", "bm_sd", "bm_min", "bm_max")) {
    expect_true(nm %in% names(fit))
    expect_equal(dim(fit[[nm]]), dim(fit$dm))
    expect_identical(rownames(fit[[nm]]), rownames(fit$dm))
    expect_identical(colnames(fit[[nm]]), colnames(fit$dm))
    expect_true(all(is.finite(fit[[nm]])))
  }
}

expect_stblr_bed_bayesr_convention <- function(fit) {
  expect_named(fit$comp_prob, colnames(fit$dm))
  for (trait in colnames(fit$dm)) {
    expect_true("component_0" %in% colnames(fit$comp_prob[[trait]]))
    expect_equal(
      unname(as.numeric(fit$dm[, trait])),
      unname(1 - fit$comp_prob[[trait]][, "component_0"]),
      tolerance = 1e-8
    )
  }
}

test_that("stblr_bed exists and is exported", {
  expect_true(is.function(stblr_bed))
  expect_true("stblr_bed" %in% getNamespaceExports("sblr"))
})

test_that("stblr_bed fits BayesC BED scheduled chains", {
  skip_if_not(
    exists("stblr_cpg_omp_bed_marker_scheduled_chains", mode = "function"),
    "native BayesC BED scheduled-chain symbol is not loaded"
  )
  fixture <- make_stblr_bed_interface_fixture()
  on.exit(unlink(fixture$bed_file), add = TRUE)

  fit <- stblr_bed(
    y = fixture$y,
    Glist = fixture$Glist,
    method = "bayesC",
    nit = 2L,
    nburn = 0L,
    full_sweep_every = 1L,
    seed = 10L,
    nchains = 1L,
    ncores = 1L
  )

  expect_true(all(c(
    "dm", "bm", "log_cpo", "mean_log_cpo", "final_pi", "mean_pi"
  ) %in% names(fit)))
  expect_stblr_bed_chain_fields(fit)
  expect_equal(fit$input$method, "bayesc")
  expect_equal(fit$input$model, "bayesc")
  expect_equal(fit$input$backend, "bed_bayesc")
  expect_equal(fit$input$data_level, "individual")
  expect_equal(fit$input$scheduled, TRUE)
  expect_equal(fit$input$nchains, 1L)
})

test_that("stblr_bed fits BayesR BED scheduled chains", {
  skip_if_not(
    exists("stblr_cpg_omp_bed_marker_scheduled_chains_bayesr", mode = "function"),
    "native BayesR BED scheduled-chain symbol is not loaded"
  )
  fixture <- make_stblr_bed_interface_fixture()
  on.exit(unlink(fixture$bed_file), add = TRUE)

  fit <- stblr_bed(
    y = fixture$y,
    Glist = fixture$Glist,
    method = "bayesR",
    nit = 2L,
    nburn = 0L,
    full_sweep_every = 1L,
    seed = 11L,
    nchains = 1L,
    ncores = 1L,
    updateE = FALSE
  )

  expect_true(all(c(
    "dm", "bm", "comp_prob", "dm_component_mean",
    "log_cpo", "mean_log_cpo", "final_pi", "mean_pi"
  ) %in% names(fit)))
  expect_stblr_bed_chain_fields(fit)
  expect_stblr_bed_bayesr_convention(fit)
  expect_equal(fit$input$method, "bayesr")
  expect_equal(fit$input$model, "bayesr")
  expect_equal(fit$input$backend, "bed_bayesr")
  expect_equal(fit$input$data_level, "individual")
  expect_equal(fit$input$scheduled, TRUE)
  expect_equal(fit$input$nchains, 1L)
})

test_that("stblr_bed normalizes method case", {
  fixture <- make_stblr_bed_interface_fixture()
  on.exit(unlink(fixture$bed_file), add = TRUE)
  common <- list(
    y = fixture$y,
    Glist = fixture$Glist,
    nit = 1L,
    nburn = 0L,
    full_sweep_every = 1L,
    nchains = 1L,
    ncores = 1L,
    updateE = FALSE
  )

  for (method in c("bayesC", "BayesC")) {
    fit <- do.call(stblr_bed, c(common, list(method = method, seed = 20L)))
    expect_equal(fit$input$method, "bayesc")
  }
  for (method in c("bayesR", "BayesR")) {
    fit <- do.call(stblr_bed, c(common, list(method = method, seed = 21L)))
    expect_equal(fit$input$method, "bayesr")
    expect_stblr_bed_bayesr_convention(fit)
  }
})

test_that("stblr_bed rejects invalid method and prior misuse clearly", {
  fixture <- make_stblr_bed_interface_fixture()
  on.exit(unlink(fixture$bed_file), add = TRUE)

  expect_error(
    stblr_bed(y = fixture$y, Glist = fixture$Glist, method = "bad"),
    "method must be one of"
  )
  expect_error(
    stblr_bed(y = fixture$y, Glist = fixture$Glist, method = "bayesr", pi_init = 0.5),
    "BayesC-specific"
  )
  expect_error(
    stblr_bed(y = fixture$y, Glist = fixture$Glist, method = "bayesc", pi = c(0.9, 0.1)),
    "BayesR-specific"
  )
  expect_error(
    stblr_bed(y = fixture$y, Glist = fixture$Glist, method = "bayesc", chain_seeds = 1L),
    "chain_seeds"
  )
})

expect_stblr_bed_bayesrc_convention <- function(fit) {
  required <- c(
    "bm", "dm", "b", "d", "vbs", "vgs", "ves", "vle", "vld", "pis",
    "covb", "covg", "cove", "pi", "pim", "comp_prob",
    "dm_component_mean", "ncomp", "alpha", "sigmaSqAlpha",
    "annotation_summary", "annotation_pi", "annotation_effects",
    "log_cpo", "mean_log_cpo", "chains", "diagnostics", "input"
  )
  expect_equal(setdiff(required, names(fit)), character())
  expect_named(fit$comp_prob, colnames(fit$dm))
  for (trait in colnames(fit$dm)) {
    cp <- fit$comp_prob[[trait]]
    expect_identical(colnames(cp)[1L], "gamma_0.00")
    expect_equal(unname(rowSums(cp)), rep(1, nrow(cp)), tolerance = 1e-12)
    expect_equal(unname(fit$dm[, trait]), unname(1 - cp[, 1L]), tolerance = 1e-12)
  }
  expect_equal(fit$vld, fit$vgs - fit$vle, tolerance = 1e-12)
}

test_that("stblr_bed publicly fits case-insensitive BayesRC", {
  skip_if_not(
    exists("stblr_cpg_omp_bed_marker_scheduled_chains_bayesrc", mode = "function"),
    "native BED BayesRC symbol is not loaded"
  )
  fixture <- make_stblr_bed_interface_fixture()
  on.exit(unlink(fixture$bed_file), add = TRUE)
  annotation <- data.frame(
    enriched = c(TRUE, FALSE),
    row.names = c("rs1", "rs2")
  )
  for (method in c("bayesrc", "BayesRC", "BAYESRC")) {
    fit <- stblr_bed(
      y = fixture$y, Glist = fixture$Glist, method = method,
      annotation = annotation, nit = 3L, nburn = 1L,
      updateAlpha = FALSE, updateE = FALSE, seed = 41L
    )
    expect_stblr_bed_bayesrc_convention(fit)
    expect_identical(fit$input$method, "bayesrc")
    expect_identical(fit$input$backend, "bed_bayesrc")
    expect_identical(fit$input$prior_type, "annotation_component")
    expect_true(fit$input$full_sweeps)
    expect_false(fit$input$adaptive_skipping)
    expect_false(fit$input$scheduled)
    expect_identical(fit$input$annotation_names, c("Intercept", "enriched"))
  }
  expect_identical(formals(stblr_bed)$method[[2L]], "bayesc")
})

test_that("public fixed-alpha BED BayesRC reduces exactly to public fixed-pi BayesR", {
  fixture <- make_stblr_bed_interface_fixture()
  on.exit(unlink(fixture$bed_file), add = TRUE)
  baseline_pi <- c(0.95, 0.03, 0.015, 0.005)
  common <- list(
    y = fixture$y, Glist = fixture$Glist,
    mixture_var = c(0, 0.01, 0.1, 1), pi = baseline_pi,
    nit = 5L, nburn = 2L, nthin = 1L, updateB = TRUE,
    updateE = TRUE, adjE = 0.9, rebuild_every = 1L,
    return_wy = TRUE, return_r = TRUE, nchains = 2L,
    ncores = 1L, seed = 17L
  )
  fit_rc <- do.call(stblr_bed, c(common, list(
    method = "bayesrc",
    annotation = matrix(1, 2L, 1L, dimnames = list(c("rs1", "rs2"), "Intercept")),
    updateAlpha = FALSE
  )))
  fit_r <- do.call(stblr_bed, c(common, list(
    method = "bayesr", alpha = rep(1, 4), updatePi = FALSE,
    full_sweep_every = 1L, null_skip_base = 1L, null_skip_max = 1L,
    candidate_threshold = 0, candidate_lifetime = 0L,
    skip_nulls_burnin_only = FALSE, progress_every = 0L
  )))
  expect_equal(fit_rc$bm, fit_r$bm, tolerance = 1e-12)
  expect_equal(fit_rc$dm, fit_r$dm, tolerance = 1e-12)
  expect_equal(
    unname(fit_rc$comp_prob[[1L]]), unname(fit_r$comp_prob[[1L]]),
    tolerance = 1e-12
  )
})

test_that("BED BayesRC aligns shuffled annotations and records unused rows", {
  fixture <- make_stblr_bed_interface_fixture()
  on.exit(unlink(fixture$bed_file), add = TRUE)
  shuffled <- data.frame(
    marker_id = c("extra", "rs2", "rs1"),
    score = c(99, 2, 1),
    enriched = c(FALSE, FALSE, TRUE)
  )
  fit <- stblr_bed(
    y = fixture$y, Glist = fixture$Glist, method = "bayesrc",
    annotation = shuffled, nit = 2L, nburn = 0L,
    updateAlpha = FALSE, updateE = FALSE, seed = 42L
  )
  expect_identical(rownames(fit$input$A), c("rs1", "rs2"))
  expect_equal(unname(fit$input$A[, "score"]), c(-sqrt(0.5), sqrt(0.5)))
  expect_identical(fit$input$annotation_alignment$unused_annotation_rows, 1L)

  aligned <- shuffled[c(3L, 2L), -1L, drop = FALSE]
  rownames(aligned) <- c("rs1", "rs2")
  fit_aligned <- stblr_bed(
    y = fixture$y, Glist = fixture$Glist, method = "bayesrc",
    annotation = aligned, nit = 2L, nburn = 0L,
    updateAlpha = FALSE, updateE = FALSE, seed = 42L
  )
  expect_identical(fit$bm, fit_aligned$bm)
  expect_identical(fit$dm, fit_aligned$dm)
  expect_identical(fit$comp_prob, fit_aligned$comp_prob)
})

test_that("BED BayesRC annotation alignment errors are informative", {
  fixture <- make_stblr_bed_interface_fixture()
  on.exit(unlink(fixture$bed_file), add = TRUE)
  run <- function(annotation, Glist = fixture$Glist) stblr_bed(
    y = fixture$y, Glist = Glist, method = "bayesrc",
    annotation = annotation, nit = 1L, nburn = 0L, updateE = FALSE
  )
  expect_error(run(data.frame(x = 1, row.names = "rs1")), "1 selected BED marker")
  expect_error(
    run(data.frame(marker_id = c("rs1", "rs1"), x = c(0, 1))),
    "Annotation marker IDs must be unique"
  )
  expect_error(
    run(data.frame(marker_id = c("rs1", NA), x = c(0, 1))),
    "non-missing"
  )
  duplicate_glist <- fixture$Glist
  duplicate_glist$rsids[[1L]] <- c("rs1", "rs1")
  duplicate_glist$rsidsLD[[1L]] <- c("rs1", "rs1")
  expect_error(
    run(matrix(c(0, 1), 2L, dimnames = list(NULL, "x")), duplicate_glist),
    "Selected BED marker IDs must be unique"
  )
})

test_that("BED BayesRC validates preprocessing and initialization", {
  fixture <- make_stblr_bed_interface_fixture()
  on.exit(unlink(fixture$bed_file), add = TRUE)
  run <- function(annotation, ...) stblr_bed(
    y = fixture$y, Glist = fixture$Glist, method = "bayesrc",
    annotation = annotation, nit = 1L, nburn = 0L,
    updateAlpha = FALSE, updateE = FALSE, ...
  )
  annotation <- data.frame(
    binary = c(TRUE, FALSE),
    category = factor(c("coding", "other")),
    row.names = c("rs1", "rs2")
  )
  fit <- run(annotation)
  expect_true(fit$input$annotation_preprocessing$intercept_added)
  expect_identical(colnames(fit$input$A)[1L], "Intercept")
  expect_true(all(fit$input$annot_alpha_init[-1L, ] == 0))
  expect_equal(unname(fit$input$annot_sigma_sq_alpha_init), rep(1, 3L))
  expect_equal(unname(fit$alpha[[1L]]), unname(fit$input$annot_alpha_init), tolerance = 1e-12)
  expect_error(run(data.frame(x = c(0, Inf), row.names = c("rs1", "rs2"))), "non-finite")
  expect_error(run(data.frame(x = c(2, 2), row.names = c("rs1", "rs2"))), "positive variance")
  bad_names <- matrix(c(0, 1, 1, 0), 2L, dimnames = list(c("rs1", "rs2"), c("x", "x")))
  expect_error(run(bad_names), "column names must be unique")
  expect_error(run(matrix(c(0, 1), 2L), annot_alpha_init = matrix(0, 3L, 3L)), "dimensions")
  expect_error(run(matrix(c(0, 1), 2L), annot_sigma_sq_alpha_init = c(1, 0, 1)), "positive finite")
})

test_that("public BED BayesRC supports traits, chains, CPO, and compact chains", {
  fixture <- make_stblr_bed_interface_fixture()
  on.exit(unlink(fixture$bed_file), add = TRUE)
  y <- cbind(D1 = fixture$y[, 1L], D2 = -0.5 * fixture$y[, 1L])
  fit <- stblr_bed(
    y = y, Glist = fixture$Glist, method = "bayesrc",
    annotation = data.frame(x = c(0, 1), row.names = c("rs1", "rs2")),
    nit = 3L, nburn = 1L, nchains = 2L, keep_chains = TRUE,
    updateAlpha = FALSE, seed = 43L
  )
  expect_stblr_bed_bayesrc_convention(fit)
  expect_equal(dim(fit$bm), c(2L, 2L))
  expect_named(fit$chains, c("D1", "D2"))
  expect_length(fit$chains$D1, 2L)
  expect_true(all(is.finite(fit$log_cpo)))
  expect_true(all(is.finite(fit$mean_log_cpo)))
})

test_that("BED BayesRC aligns annotations after chr and cls selection", {
  bed1 <- tempfile(fileext = ".bed")
  bed2 <- tempfile(fileext = ".bed")
  on.exit(unlink(c(bed1, bed2)), add = TRUE)
  write_stblr_bed_interface_file(
    bed1, rbind(c(0, 1, 2, 0, 1, 2), c(2, 0, 1, 2, 0, 1))
  )
  write_stblr_bed_interface_file(
    bed2, rbind(c(0, 2, 1, 0, 2, 1), c(1, 2, 0, 1, 2, 0))
  )
  glist <- list(
    n = 6L, ids = paste0("id", 1:6), bedfiles = c(bed1, bed2),
    rsids = list(c("c1m1", "c1m2"), c("c2m1", "c2m2")),
    rsidsLD = list(c("c1m1", "c1m2"), c("c2m1", "c2m2")),
    af = list(c(0.5, 0.5), c(0.5, 0.5))
  )
  y <- matrix(c(-1, 0, 1, -0.5, 0.5, 1.5), ncol = 1L,
              dimnames = list(NULL, "D1"))
  annotation <- data.frame(
    marker_id = c("c1m2", "c2m2", "unused", "c1m1", "c2m1"),
    enriched = c(0, 1, 0, 0, 1)
  )
  fit <- stblr_bed(
    y = y, Glist = glist, method = "bayesrc", annotation = annotation,
    chr = c(2L, 1L), cls = list(2L, 1L), nit = 2L, nburn = 0L,
    updateAlpha = FALSE, updateE = FALSE, seed = 44L
  )
  expect_identical(fit$input$selected_marker_ids, c("c2m2", "c1m1"))
  expect_identical(rownames(fit$input$A), c("c2m2", "c1m1"))
  expect_equal(unname(fit$input$A[, "enriched"]), c(1, 0))
  expect_identical(fit$input$annotation_alignment$unused_annotation_rows, 3L)
})

test_that("public BED BayesRC formatting preserves direct internal raw values", {
  fixture <- make_stblr_bed_interface_fixture()
  on.exit(unlink(fixture$bed_file), add = TRUE)
  annotation <- matrix(
    c(1, 1, 0, 1), 2L, 2L,
    dimnames = list(c("rs1", "rs2"), c("Intercept", "enriched"))
  )
  fit <- stblr_bed(
    y = fixture$y, Glist = fixture$Glist, method = "bayesrc",
    annotation = annotation, pi = c(0.95, 0.03, 0.015, 0.005),
    nit = 3L, nburn = 1L, updateAlpha = FALSE, updateE = FALSE,
    return_wy = TRUE, return_r = TRUE, seed = 45L
  )
  dat <- sblr:::.make_bed_marker_data(
    fixture$Glist, fixture$y, chr = 1L, cls = NULL, block_size = 1000L
  )
  raw <- sblr:::.stblr_bed_bayesrc_native(
    bed_files = dat$bed_files, n = dat$n_total, cls = dat$cls, y = dat$y,
    b_init = dat$b_init, sets = dat$sets, rows = dat$rows, af = dat$af,
    scale = TRUE, B = fit$input$B, E = fit$input$E,
    ssb_prior = split(fit$input$ssb_prior, rep(seq_len(dat$nt), each = dat$nt)),
    sse_prior = split(fit$input$sse_prior, rep(seq_len(dat$nt), each = dat$nt)),
    A = fit$input$A, gamma = fit$input$mixture_var,
    annot_alpha_init = fit$input$annot_alpha_init,
    annot_sigma_sq_alpha_init = fit$input$annot_sigma_sq_alpha_init,
    intercept_flat = fit$input$intercept_flat,
    sigmaSqAlpha_a = fit$input$sigmaSqAlpha_a,
    sigmaSqAlpha_b = fit$input$sigmaSqAlpha_b, pi_floor = fit$input$pi_floor,
    nub = fit$input$nub, nue = fit$input$nue, updateAlpha = FALSE,
    updateB = TRUE, updateE = FALSE,
    annot_alpha_update_every = fit$input$annot_alpha_update_every,
    adjE = fit$input$adjE, nit = fit$input$nit, nburn = fit$input$nburn,
    nthin = fit$input$nthin, rebuild_every = fit$input$rebuild_every,
    return_wy = TRUE, return_r = TRUE, read_block_size = fit$input$read_block_size,
    nchains = 1L, keep_chains = FALSE, ncores = 1L, seed = fit$input$seed
  )
  expect_equal(unname(fit$bm), unname(raw$marker$bm), tolerance = 1e-12)
  expect_equal(unname(fit$dm), unname(raw$marker$dm), tolerance = 1e-12)
  expect_equal(
    unname(fit$comp_prob[[1L]]), unname(raw$component$prob[[1L]]),
    tolerance = 1e-12
  )
})

test_that("BED BayesRC preserves explicit fixed annotation initialization", {
  fixture <- make_stblr_bed_interface_fixture()
  on.exit(unlink(fixture$bed_file), add = TRUE)
  A <- matrix(c(0, 1), 2L, 1L, dimnames = list(c("rs1", "rs2"), "enriched"))
  alpha_init <- matrix(seq_len(6L) / 10, 2L, 3L)
  sigma_init <- c(0.5, 1.5, 2.5)
  fit <- stblr_bed(
    y = fixture$y, Glist = fixture$Glist, method = "bayesrc", annotation = A,
    annot_alpha_init = alpha_init,
    annot_sigma_sq_alpha_init = sigma_init,
    updateAlpha = FALSE, updateE = FALSE, nit = 2L, nburn = 0L, seed = 46L
  )
  expect_equal(unname(fit$alpha[[1L]]), unname(alpha_init), tolerance = 1e-12)
  expect_equal(unname(fit$sigmaSqAlpha[1L, ]), sigma_init, tolerance = 1e-12)
})

test_that("BED BayesRC rejects deferred features and warns on scheduling controls", {
  fixture <- make_stblr_bed_interface_fixture()
  on.exit(unlink(fixture$bed_file), add = TRUE)
  common <- list(
    y = fixture$y, Glist = fixture$Glist, method = "bayesrc",
    annotation = matrix(c(0, 1), 2L), nit = 1L, nburn = 0L, updateE = FALSE
  )
  expect_error(do.call(stblr_bed, c(common, list(selection_s = 0))), "Unsupported argument")
  expect_error(do.call(stblr_bed, c(common, list(updateLDswap = TRUE))), "Unsupported argument")
  expect_error(do.call(stblr_bed, c(common, list(eigen_filter = "hard_truncate"))), "Unsupported argument")
  expect_warning(
    do.call(stblr_bed, c(common, list(full_sweep_every = 2L))),
    "always uses unscheduled full sweeps"
  )
})
