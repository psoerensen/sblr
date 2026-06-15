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
    chains = FALSE,
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
    "neither rows nor rownames\\(y\\) were supplied"
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
  expect_error(
    make_bed_marker_test_data(Glist, y_sub),
    "Some rownames\\(y\\) were not found"
  )
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
    seed = 10L, full_sweep_every = 1L, scheduled = TRUE
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
