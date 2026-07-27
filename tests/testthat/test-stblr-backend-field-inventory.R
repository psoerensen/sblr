make_field_inventory_csr_prefix <- function(m = 4L) {
  prefix <- tempfile("stblr_field_inventory_csr_")
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

make_field_inventory_stats <- function(m = 4L) {
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

make_field_inventory_annotations <- function(m = 4L) {
  markers <- paste0("m", seq_len(m))
  matrix(
    c(1, 0, 0, 1, 1, 1, 0, 0),
    nrow = m,
    byrow = TRUE,
    dimnames = list(markers, c("coding", "qtl"))
  )
}

write_field_inventory_bed_file <- function(path, dosage) {
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

make_field_inventory_bed_fixture <- function() {
  bed_file <- tempfile(fileext = ".bed")
  write_field_inventory_bed_file(
    bed_file,
    rbind(c(0, 1, 2, 0, 1, 2), c(2, 0, 1, 2, 0, 1))
  )
  list(
    bed_file = bed_file,
    Glist = list(
      n = 6L,
      ids = paste0("id", seq_len(6L)),
      bedfiles = bed_file,
      rsids = list(c("rs1", "rs2")),
      rsidsLD = list(c("rs1", "rs2")),
      chr = list(c(1L, 1L)),
      pos = list(c(100, 200)),
      af = list(c(0.2, 0.3))
    ),
    y = matrix(
      c(-1, 0, 1, -0.5, 0.5, 1.5),
      ncol = 1,
      dimnames = list(NULL, "D1")
    )
  )
}

# Captures the raw stblr_raw object a backend call actually hands to
# .as_stblr_fit(), by tracing that internal formatter for the duration of
# `expr`. This lets the inventory test assert on the real native raw output
# (schema validity + .validate_stblr_raw() structural checks) in addition to
# the formatted fit, without altering any backend or formatter behaviour.
capture_stblr_raw <- function(expr) {
  ns <- asNamespace("sblr")
  assign(".stblr_last_raw_capture", NULL, envir = .GlobalEnv)
  suppressMessages(trace(
    ".as_stblr_fit",
    tracer = quote(assign(".stblr_last_raw_capture", raw, envir = .GlobalEnv)),
    where = ns, print = FALSE
  ))
  on.exit({
    suppressMessages(try(untrace(".as_stblr_fit", where = ns), silent = TRUE))
  }, add = TRUE)
  fit <- eval.parent(substitute(expr))
  list(fit = fit, raw = get(".stblr_last_raw_capture", envir = .GlobalEnv))
}

expect_field_inventory_raw <- function(raw, backend) {
  expect_true(sblr:::.is_stblr_raw(raw), info = backend)
  expect_true(sblr:::.validate_stblr_raw(raw, backend = backend), info = backend)
  expect_identical(raw$schema$class, "stblr_raw", info = backend)
  expect_identical(as.integer(raw$schema$version), 1L, info = backend)
}

expect_field_inventory_core <- function(fit, backend, data_level, annotations) {
  expect_true(all(c("dm", "bm", "input") %in% names(fit)))
  expect_true(is.matrix(fit$dm))
  expect_true(is.matrix(fit$bm))
  expect_equal(dim(fit$dm), dim(fit$bm))
  expect_false(is.null(rownames(fit$dm)))
  expect_false(is.null(colnames(fit$dm)))
  expect_true(all(is.na(fit$dm) | (fit$dm >= -1e-8 & fit$dm <= 1 + 1e-8)))
  expect_true(all(is.na(fit$bm) | is.finite(fit$bm)))
  expect_equal(fit$input$backend, backend)
  expect_equal(fit$input$data_level, data_level)
  expect_identical(fit$input$annotations, annotations)

  chk <- sblr::check_stblr_consistency(fit, verbose = FALSE)
  expect_true(chk$ok, info = paste(chk$checks$message[!chk$checks$ok],
                                   collapse = "; "))
}

expect_field_inventory_csr_variance <- function(fit) {
  for (nm in c("vbs", "vgs", "ves", "vle", "vld")) {
    expect_true(nm %in% names(fit), info = nm)
    expect_true(is.matrix(fit[[nm]]), info = nm)
    expect_identical(colnames(fit[[nm]]), colnames(fit$dm), info = nm)
    expect_true(all(is.na(fit[[nm]]) | is.finite(fit[[nm]])), info = nm)
  }
}

expect_field_inventory_chain_summaries <- function(fit) {
  for (nm in c("dm_chain_mean_sd", "dm_chain_mean_min", "dm_chain_mean_max",
               "bm_chain_mean_sd", "bm_chain_mean_min", "bm_chain_mean_max")) {
    expect_true(nm %in% names(fit), info = nm)
    expect_equal(dim(fit[[nm]]), dim(fit$dm), info = nm)
  }
}

expect_field_inventory_chains <- function(fit) {
  expect_true("chains" %in% names(fit))
  chain <- fit$chains[[1L]]
  expect_identical(names(chain$dm), rownames(fit$dm))
  expect_identical(names(chain$bm), rownames(fit$bm))
}

expect_field_inventory_bayesr <- function(fit) {
  expect_true(all(c("component_probabilities", "dm_component_mean") %in% names(fit)))
  expect_named(fit$component_probabilities, colnames(fit$dm))
  trait <- colnames(fit$dm)[1L]
  expect_equal(fit$dm[, trait], 1 - fit$component_probabilities[[trait]][, 1L], tolerance = 1e-8)
  expect_equal(dim(fit$dm_component_mean), dim(fit$dm))
}

test_that("CSR BayesC and BayesR expose common field inventory", {
  stats <- make_field_inventory_stats()

  skip_if_not(exists("stblr_cpg_omp_csr", mode = "function"))
  res_c <- capture_stblr_raw(sblr::stblr_csr(
    stats = stats,
    ld_prefix = make_field_inventory_csr_prefix(stats$m),
    method = "bayesc",
    pi_init = 0.35,
    pi_prior_mean = 0.35,
    pi_prior_strength = 2,
    updateB = FALSE,
    updateE = FALSE,
    updatePi = FALSE,
    nit = 3,
    nburn = 0,
    nchains = 2L,
    keep_chains = TRUE,
    chain_seeds = c(101L, 102L),
    ncores = 1L,
    seed = 101L,
    updateLDswap = TRUE,
    ld_swap_prob = 1,
    ld_swap_r2 = 0
  ))
  fit_c <- res_c$fit
  expect_field_inventory_raw(res_c$raw, "csr_bayesc")
  expect_field_inventory_core(fit_c, "csr_bayesc", "summary", FALSE)
  expect_field_inventory_csr_variance(fit_c)
  expect_field_inventory_chain_summaries(fit_c)
  expect_field_inventory_chains(fit_c)
  expect_true("ld_swap" %in% names(fit_c$diagnostics))

  skip_if_not(exists("stblr_cpg_omp_csr_bayesr", mode = "function"))
  res_r <- capture_stblr_raw(sblr::stblr_csr(
    stats = stats,
    ld_prefix = make_field_inventory_csr_prefix(stats$m),
    method = "bayesr",
    updateB = FALSE,
    updateE = FALSE,
    updatePi = FALSE,
    nit = 3,
    nburn = 0,
    nchains = 1L,
    ncores = 1L,
    seed = 102L
  ))
  fit_r <- res_r$fit
  expect_field_inventory_raw(res_r$raw, "csr_bayesr")
  expect_field_inventory_core(fit_r, "csr_bayesr", "summary", FALSE)
  expect_field_inventory_csr_variance(fit_r)
  expect_field_inventory_bayesr(fit_r)
})

test_that("annotation-aware CSR models expose aligned common and annotation fields", {
  stats <- make_field_inventory_stats()
  A <- make_field_inventory_annotations()

  skip_if_not(exists("stblr_cpg_omp_csr_prior", mode = "function"))
  res_prior <- capture_stblr_raw(sblr::stblr_csr_annot(
    stats = stats,
    ld_prefix = make_field_inventory_csr_prefix(stats$m),
    annotations = list(
      A = A,
      fixed_pi_marker = list(rep(0.35, stats$m)),
      fixed_vb_multiplier = list(rep(1, stats$m))
    ),
    annotation_model = "prior",
    pi_init = 0.35,
    pi_prior_mean = 0.35,
    pi_prior_strength = 2,
    updateB = FALSE,
    updateE = FALSE,
    updatePi = FALSE,
    nit = 3,
    nburn = 0,
    nchains = 2L,
    keep_chains = TRUE,
    chain_seeds = c(201L, 202L),
    convergence = "none",
    ncores = 1L,
    seed = 201L
  ))
  fit_prior <- res_prior$fit
  expect_field_inventory_raw(res_prior$raw, "csr_prior_bayesc")
  expect_field_inventory_core(fit_prior, "csr_prior_bayesc", "summary", TRUE)
  expect_field_inventory_csr_variance(fit_prior)
  expect_field_inventory_chain_summaries(fit_prior)
  expect_field_inventory_chains(fit_prior)
  expect_true("annotation_prior" %in% names(fit_prior))

  skip_if_not(exists("stblr_cpg_omp_csr_annot", mode = "function"))
  res_learned <- capture_stblr_raw(sblr::stblr_csr_annot(
    stats = stats,
    ld_prefix = make_field_inventory_csr_prefix(stats$m),
    annotations = A,
    annotation_model = "learned",
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
    seed = 202L
  ))
  fit_learned <- res_learned$fit
  expect_field_inventory_raw(res_learned$raw, "csr_annot_bayesc")
  expect_field_inventory_core(fit_learned, "csr_annot_bayesc", "summary", TRUE)
  expect_field_inventory_csr_variance(fit_learned)
  expect_true(all(c("annotation_effects", "eta_pi", "eta_vb") %in% names(fit_learned)))

  skip_if_not(exists("stblr_cpg_omp_csr_group_annot", mode = "function"))
  res_group <- capture_stblr_raw(sblr::stblr_csr_annot(
    stats = stats,
    ld_prefix = make_field_inventory_csr_prefix(stats$m),
    annotations = stats::setNames(
      c("coding", "background", "coding", "background"),
      stats$marker_names
    ),
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
    seed = 203L,
    updateLDswap = TRUE,
    ld_swap_prob = 1,
    ld_swap_r2 = 0
  ))
  fit_group <- res_group$fit
  expect_field_inventory_raw(res_group$raw, "csr_group_bayesc")
  expect_field_inventory_core(fit_group, "csr_group_bayesc", "summary", TRUE)
  expect_field_inventory_csr_variance(fit_group)
  expect_true("ld_swap" %in% names(fit_group$diagnostics))
  expect_true(all(c(
    "annotation_pi", "annotation_variance", "annotation_summary",
    "group_pi", "group_vb_multiplier", "group_nincluded", "group_size"
  ) %in% names(fit_group)))

  skip_if_not(exists("stblr_cpg_omp_csr_sbayesrc", mode = "function"))
  res_sbayesrc <- capture_stblr_raw(sblr::stblr_csr_annot(
    stats = stats,
    Glist = list(
      rsidsLD = list(stats$marker_names),
      rsids = list(stats$marker_names),
      maf = list(rep(0.2, stats$m))
    ),
    ld_prefix = make_field_inventory_csr_prefix(stats$m),
    annotations = A,
    annotation_model = "sbayesrc",
    method = "sbayesrc",
    gamma = c(0, 0.1, 1),
    pi_init = 0.35,
    pi_prior_mean = 0.35,
    pi_prior_strength = 2,
    updateAlpha = FALSE,
    updateB = FALSE,
    updateE = FALSE,
    nit = 3,
    nburn = 0,
    ncores = 1L,
    seed = 204L
  ))
  fit_sbayesrc <- res_sbayesrc$fit
  expect_field_inventory_raw(res_sbayesrc$raw, "csr_sbayesrc")
  expect_field_inventory_core(fit_sbayesrc, "csr_sbayesrc", "summary", TRUE)
  expect_field_inventory_csr_variance(fit_sbayesrc)
  expect_field_inventory_bayesr(fit_sbayesrc)
  expect_true(all(c(
    "annotation_pi", "annotation_summary", "annotation_effects",
    "annotation_variance", "alpha", "sigmaSqAlpha", "ncomp"
  ) %in% names(fit_sbayesrc)))
})

test_that("BED BayesC and BayesR expose supported field inventory", {
  fixture <- make_field_inventory_bed_fixture()
  on.exit(unlink(fixture$bed_file), add = TRUE)

  skip_if_not(
    exists("stblr_cpg_omp_bed_marker_scheduled_chains", mode = "function")
  )
  res_c <- capture_stblr_raw(sblr::stblr_bed(
    y = fixture$y,
    Glist = fixture$Glist,
    method = "bayesc",
    nit = 2L,
    nburn = 0L,
    full_sweep_every = 1L,
    seed = 301L,
    nchains = 1L,
    ncores = 1L
  ))
  fit_c <- res_c$fit
  expect_field_inventory_raw(res_c$raw, "bed_bayesc")
  expect_field_inventory_core(fit_c, "bed_bayesc", "individual", FALSE)
  expect_field_inventory_csr_variance(fit_c)
  expect_field_inventory_chain_summaries(fit_c)

  skip_if_not(
    exists("stblr_cpg_omp_bed_marker_scheduled_chains_bayesr", mode = "function")
  )
  res_r <- capture_stblr_raw(sblr::stblr_bed(
    y = fixture$y,
    Glist = fixture$Glist,
    method = "bayesr",
    nit = 2L,
    nburn = 0L,
    full_sweep_every = 1L,
    seed = 302L,
    nchains = 1L,
    ncores = 1L,
    updateE = FALSE
  ))
  fit_r <- res_r$fit
  expect_field_inventory_raw(res_r$raw, "bed_bayesr")
  expect_field_inventory_core(fit_r, "bed_bayesr", "individual", FALSE)
  expect_field_inventory_csr_variance(fit_r)
  expect_field_inventory_bayesr(fit_r)
  expect_field_inventory_chain_summaries(fit_r)
})

test_that("every inventoried backend rejects non-stblr_raw output instead of a silent legacy fallback", {
  inventoried_backends <- c(
    "csr_bayesc", "csr_bayesr", "csr_prior_bayesc", "csr_annot_bayesc",
    "csr_group_bayesc", "csr_sbayesrc", "bed_bayesc", "bed_bayesr"
  )
  legacy_positional_raw <- unname(list(1, 2, 3))
  expect_false(sblr:::.is_stblr_raw(legacy_positional_raw))

  for (backend in inventoried_backends) {
    expect_error(
      sblr:::.validate_stblr_raw(legacy_positional_raw, backend = backend),
      paste0("Unsupported backend output.*", backend),
      info = backend
    )
  }
})
