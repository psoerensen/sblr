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
