make_bed_marker_test_glist <- function(n = 6L) {
  list(
    n = n,
    ids = paste0("id", seq_len(n)),
    bedfiles = "unused.bed",
    rsids = list(c("rs1", "rs2")),
    rsidsLD = list(c("rs1", "rs2")),
    af = list(c(0.2, 0.3))
  )
}

make_bed_marker_test_data <- function(Glist, y, rows = NULL) {
  sblr:::.make_bed_marker_data(
    Glist = Glist,
    y = y,
    chr = 1L,
    cls = NULL,
    block_size = 1000L,
    rows = rows
  )
}

write_bed_marker_test_file <- function(path, dosage) {
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

test_that("BED marker data distinguishes total and selected sample counts", {
  Glist <- make_bed_marker_test_glist()
  y <- matrix(seq_len(Glist$n), ncol = 1)

  full <- make_bed_marker_test_data(Glist, y)

  expect_null(full$rows)
  expect_equal(full$n, Glist$n)
  expect_equal(full$n_total, Glist$n)
  expect_equal(full$n_used, Glist$n)

  idx <- c(2L, 4L, 6L)
  explicit <- make_bed_marker_test_data(Glist, y[idx, , drop = FALSE], idx)

  expect_identical(explicit$rows, idx)
  expect_equal(explicit$n, Glist$n)
  expect_equal(explicit$n_total, Glist$n)
  expect_equal(explicit$n_used, length(idx))
})

test_that("explicit and rowname-matched BED subsets resolve identically", {
  Glist <- make_bed_marker_test_glist()
  y <- matrix(seq_len(Glist$n), ncol = 1)
  idx <- c(2L, 4L, 6L)
  y_sub <- y[idx, , drop = FALSE]

  explicit <- make_bed_marker_test_data(Glist, y_sub, idx)
  rownames(y_sub) <- Glist$ids[idx]
  matched <- make_bed_marker_test_data(Glist, y_sub)

  expect_identical(explicit$rows, matched$rows)
  expect_equal(explicit$n, matched$n)
  expect_equal(explicit$n_used, matched$n_used)
})

test_that("BED marker subset rows and IDs are unambiguous", {
  Glist <- make_bed_marker_test_glist()
  y_sub <- matrix(1:3, ncol = 1)

  expect_error(
    make_bed_marker_test_data(Glist, y_sub),
    "nrow\\(y\\) != Glist\\$n, but rownames\\(y\\) are missing"
  )

  expect_error(
    make_bed_marker_test_data(Glist, y_sub, c(1L, 1L, 2L)),
    "rows must not contain duplicate"
  )

  rownames(y_sub) <- c("id1", "id1", "id2")
  expect_error(
    make_bed_marker_test_data(Glist, y_sub),
    "rownames\\(y\\) must not contain duplicates"
  )

  rownames(y_sub) <- c("id1", "id2", "id3")
  Glist$ids[2] <- Glist$ids[1]
  expect_error(
    make_bed_marker_test_data(Glist, y_sub),
    "Glist\\$ids must not contain duplicates"
  )

  Glist <- make_bed_marker_test_glist()
  rownames(y_sub) <- c("id1", "id2", "missing")
  expect_warning(
    missing <- make_bed_marker_test_data(Glist, y_sub),
    "phenotype IDs were not found in Glist\\$ids and will be dropped"
  )
  expect_identical(missing$rows, c(1L, 2L))
  expect_equal(missing$n_used, 2L)
})

test_that("BED marker fits report total and selected sample counts", {
  bed_file <- tempfile(fileext = ".bed")
  on.exit(unlink(bed_file), add = TRUE)
  write_bed_marker_test_file(
    bed_file,
    rbind(c(0, 1, 2, 0, 1, 2), c(2, 0, 1, 2, 0, 1))
  )

  Glist <- make_bed_marker_test_glist()
  Glist$bedfiles <- bed_file
  y <- matrix(c(-1, 0, 1, -0.5, 0.5, 1.5), ncol = 1)
  idx <- c(2L, 4L, 6L)
  y_sub <- y[idx, , drop = FALSE]
  fit_args <- list(
    Glist = Glist, chr = 1L, nit = 2L, nburn = 1L, ncores = 1L,
    seed = 10L, full_sweep_every = 1L, backend = "scheduled"
  )

  fit_explicit <- do.call(
    stblr_bed_marker,
    c(fit_args, list(y = y_sub, rows = idx))
  )
  rownames(y_sub) <- Glist$ids[idx]
  fit_matched <- do.call(
    stblr_bed_marker,
    c(fit_args, list(y = y_sub, rows = NULL))
  )

  expect_identical(fit_explicit$input$rows, fit_matched$input$rows)
  expect_equal(fit_explicit$input$n, Glist$n)
  expect_equal(fit_explicit$input$n_total, Glist$n)
  expect_equal(fit_explicit$input$n_used, nrow(y_sub))
})

make_bed_marker_chain_test_files <- function() {
  bed_file1 <- tempfile(fileext = ".bed")
  bed_file2 <- tempfile(fileext = ".bed")
  write_bed_marker_test_file(
    bed_file1,
    rbind(c(0, 1, 2, 0, 1, 2), c(2, 0, 1, 2, 0, 1))
  )
  write_bed_marker_test_file(
    bed_file2,
    rbind(c(1, 2, 0, 1, 2, 0))
  )
  list(bed_file1, bed_file2)
}

make_bed_marker_chain_test_glist <- function(files) {
  list(
    n = 6L,
    ids = paste0("id", seq_len(6L)),
    bedfiles = unlist(files, use.names = FALSE),
    rsids = list(c("rs1", "rs2"), "rs3"),
    rsidsLD = list(c("rs1", "rs2"), "rs3"),
    chr = list(c(1L, 1L), 2L),
    pos = list(c(100, 200), 300),
    af = list(c(0.2, 0.3), 0.4)
  )
}

expect_bed_chain_summary_shape <- function(fit) {
  for (nm in c("dm_sd", "dm_min", "dm_max", "bm_sd", "bm_min", "bm_max")) {
    expect_true(nm %in% names(fit))
    expect_equal(dim(fit[[nm]]), dim(fit$dm))
    expect_identical(rownames(fit[[nm]]), rownames(fit$dm))
    expect_identical(colnames(fit[[nm]]), colnames(fit$dm))
    expect_true(all(is.finite(fit[[nm]])))
  }
}

test_that("BED scheduled chains expose zero-width summaries for one chain", {
  files <- make_bed_marker_chain_test_files()
  on.exit(unlink(unlist(files, use.names = FALSE)), add = TRUE)

  fit <- stblr_bed_marker(
    Glist = make_bed_marker_chain_test_glist(files),
    y = matrix(c(-1, 0, 1, -0.5, 0.5, 1.5), ncol = 1,
               dimnames = list(NULL, "D1")),
    chr = 1:2,
    backend = "scheduled",
    nit = 2L,
    nburn = 0L,
    full_sweep_every = 1L,
    nchains = 1L,
    ncores = 1L,
    seed = 11L
  )

  expect_true(fit$input$use_chains_backend)
  expect_equal(fit$input$nchains, 1L)
  expect_bed_chain_summary_shape(fit)
  expect_equal(fit$dm_sd, fit$dm * 0, tolerance = 1e-12)
  expect_equal(fit$bm_sd, fit$bm * 0, tolerance = 1e-12)
  expect_equal(fit$dm_min, fit$dm, tolerance = 1e-12)
  expect_equal(fit$dm_max, fit$dm, tolerance = 1e-12)
  expect_equal(fit$bm_min, fit$bm, tolerance = 1e-12)
  expect_equal(fit$bm_max, fit$bm, tolerance = 1e-12)
})

test_that("BED scheduled chains expose finite summaries across chains", {
  bed_file <- tempfile(fileext = ".bed")
  on.exit(unlink(bed_file), add = TRUE)
  write_bed_marker_test_file(
    bed_file,
    rbind(c(0, 1, 2, 0, 1, 2), c(2, 0, 1, 2, 0, 1))
  )

  Glist <- make_bed_marker_test_glist()
  Glist$bedfiles <- bed_file
  Glist$chr <- list(c(1L, 1L))
  Glist$pos <- list(c(100, 200))

  fit <- stblr_bed_marker(
    Glist = Glist,
    y = matrix(c(-1, 0, 1, -0.5, 0.5, 1.5), ncol = 1,
               dimnames = list(NULL, "D1")),
    chr = 1L,
    backend = "scheduled",
    nit = 2L,
    nburn = 0L,
    full_sweep_every = 1L,
    nchains = 2L,
    ncores = 1L,
    seed = 12L
  )

  expect_true(fit$input$use_chains_backend)
  expect_equal(fit$input$nchains, 2L)
  expect_bed_chain_summary_shape(fit)
  expect_true(all(fit$dm_sd >= -1e-12))
  expect_true(all(fit$bm_sd >= -1e-12))
  expect_true(all(fit$dm_min <= fit$dm + 1e-12))
  expect_true(all(fit$dm <= fit$dm_max + 1e-12))
  expect_true(all(fit$bm_min <= fit$bm + 1e-12))
  expect_true(all(fit$bm <= fit$bm_max + 1e-12))

  fm <- extract_stblr_finemap_loci(
    fit = fit,
    Glist = Glist,
    locus_sets = list(regionA = rownames(fit$dm)),
    trait = "D1",
    credible_sets = FALSE
  )
  expect_equal(fm$markers$pip_sd, as.numeric(fit$dm_sd[, "D1"]))
  expect_equal(fm$markers$pip_min, as.numeric(fit$dm_min[, "D1"]))
  expect_equal(fm$markers$pip_max, as.numeric(fit$dm_max[, "D1"]))
  expect_equal(fm$markers$bm_sd, as.numeric(fit$bm_sd[, "D1"]))
  expect_equal(fm$markers$bm_min, as.numeric(fit$bm_min[, "D1"]))
  expect_equal(fm$markers$bm_max, as.numeric(fit$bm_max[, "D1"]))
})
