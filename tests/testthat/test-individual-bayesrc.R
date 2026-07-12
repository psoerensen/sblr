test_that("individual BayesRC helpers are internal", {
  expect_true(exists(".stblr_bed_bayesrc_native", envir = asNamespace("sblr")))
  expect_false(".stblr_bed_bayesrc_native" %in% getNamespaceExports("sblr"))
  expect_true(exists(".bayesr_pi_to_probit_stick_intercepts", envir = asNamespace("sblr")))
  expect_false(".bayesr_pi_to_probit_stick_intercepts" %in% getNamespaceExports("sblr"))
})

test_that("individual BayesRC uses the shared BED utility header", {
  source_path <- testthat::test_path(
    "..", "..", "src", "stblr_cpg_omp_bed_marker_scheduled_chains_bayesrc.cpp"
  )
  source <- paste(readLines(source_path, warn = FALSE), collapse = "\n")
  expect_match(source, '#include "st_bed_bayesr_common.h"', fixed = TRUE)
  expect_false(grepl(
    'include "stblr_cpg_omp_bed_marker_scheduled_chains_bayesr.cpp"',
    source,
    fixed = TRUE
  ))
})

test_that("intercept conversion reproduces BayesR component probabilities", {
  target <- c(0.95, 0.03, 0.015, 0.005)
  intercept <- sblr:::.bayesr_pi_to_probit_stick_intercepts(target)
  stick <- stats::pnorm(intercept[1L, ])
  got <- numeric(length(target))
  remaining <- 1
  for (k in seq_along(stick)) {
    got[k] <- remaining * (1 - stick[k])
    remaining <- remaining * stick[k]
  }
  got[length(got)] <- remaining
  expect_equal(got, target, tolerance = 1e-12)
})

test_that("individual BayesRC native symbol is registered when compiled", {
  ok <- tryCatch({
    getNativeSymbolInfo(
      "_sblr_stblr_cpg_omp_bed_marker_scheduled_chains_bayesrc",
      PACKAGE = "sblr"
    )
    TRUE
  }, error = function(e) FALSE)
  skip_if_not(ok, "native individual BayesRC symbol is not loaded")
  expect_true(is.function(sblr:::.stblr_bed_bayesrc_native))
})

make_individual_bayesrc_fixture <- function() {
  bed <- tempfile(fileext = ".bed")
  dosage <- rbind(c(0, 1, 2, 0, 1, 2), c(2, 0, 1, 2, 0, 1))
  code <- c(`0` = 3L, `1` = 2L, `2` = 0L)
  packed <- unlist(lapply(seq_len(nrow(dosage)), function(j) {
    z <- unname(code[as.character(dosage[j, ])])
    z <- c(z, rep(0L, (-length(z)) %% 4L))
    vapply(seq(1L, length(z), 4L), function(i) {
      sum(z[i:(i + 3L)] * c(1L, 4L, 16L, 64L))
    }, integer(1))
  }))
  writeBin(as.raw(c(0x6c, 0x1b, 0x01, packed)), bed)
  list(
    bed = bed,
    y = matrix(c(-1, 0, 1, -0.5, 0.5, 1.5), ncol = 1L),
    gamma = c(0, 0.01, 0.1, 1),
    pi = c(0.95, 0.03, 0.015, 0.005)
  )
}

test_that("individual BayesRC native backend runs and reduces to fixed-pi BayesR", {
  bayesrc_ok <- tryCatch({
    getNativeSymbolInfo("_sblr_stblr_cpg_omp_bed_marker_scheduled_chains_bayesrc", PACKAGE = "sblr")
    TRUE
  }, error = function(e) FALSE)
  bayesr_ok <- tryCatch({
    getNativeSymbolInfo("_sblr_stblr_cpg_omp_bed_marker_scheduled_chains_bayesr", PACKAGE = "sblr")
    TRUE
  }, error = function(e) FALSE)
  skip_if_not(bayesrc_ok && bayesr_ok, "native BED BayesRC/BayesR symbols are not loaded")

  x <- make_individual_bayesrc_fixture()
  intercept <- sblr:::.bayesr_pi_to_probit_stick_intercepts(x$pi)
  common <- list(
    bed_files = x$bed, n = 6L, cls = list(1:2), y = x$y,
    b_init = list(c(0, 0)), sets = c(1L, 1L), rows = NULL,
    af = list(c(0.5, 0.5)), scale = TRUE,
    B = matrix(0.1, 1L, 1L), E = matrix(1, 1L, 1L),
    ssb_prior = list(0.05), sse_prior = list(0.5),
    nub = 4, nue = 4, updateB = TRUE, updateE = TRUE,
    adjE = 0.9, nit = 5L, nburn = 2L, nthin = 1L,
    rebuild_every = 1L, return_wy = TRUE, return_r = TRUE,
    read_block_size = 2L, nchains = 2L, ncores = 1L, seed = 17L
  )
  raw_rc <- do.call(sblr:::.stblr_bed_bayesrc_native, c(common, list(
    A = matrix(1, 2L, 1L), gamma = x$gamma,
    annot_alpha_init = intercept,
    annot_sigma_sq_alpha_init = rep(1, length(x$gamma) - 1L),
    updateAlpha = FALSE, annot_alpha_update_every = 10L
  )))
  expect_s3_class(raw_rc, "stblr_raw_v1")
  expect_identical(raw_rc$meta$model, "bayesrc")
  expect_identical(raw_rc$meta$backend, "bed_bayesrc")
  expect_true(isTRUE(raw_rc$diagnostics$full_sweeps))
  expect_false(isTRUE(raw_rc$diagnostics$adaptive_skipping))
  expect_null(raw_rc$chains)
  expect_identical(raw_rc$meta$nchains, 2L)
  cp <- raw_rc$component$prob[[1L]]
  expect_identical(colnames(cp), NULL)
  expect_true(all(is.finite(cp) & cp >= 0 & cp <= 1))
  expect_equal(rowSums(cp), rep(1, nrow(cp)), tolerance = 1e-12)
  expect_equal(raw_rc$marker$dm[, 1L], 1 - cp[, 1L], tolerance = 1e-12)
  expect_identical(raw_rc$component$names[1L], "gamma_0.00")
  expect_equal(raw_rc$pi$final[1L, ], x$pi, tolerance = 1e-12)
  expect_equal(raw_rc$annotation$alpha_mean[[1L]], intercept, tolerance = 1e-12)
  expect_equal(dim(raw_rc$annotation$sigmaSqAlpha_mean), c(3L, 1L))

  raw_r <- do.call(sblr:::stblr_cpg_omp_bed_marker_scheduled_chains_bayesr, c(common, list(
    pi = x$pi, c = x$gamma, alpha = rep(1, length(x$pi)),
    updatePi = FALSE, full_sweep_every = 1L, null_skip_base = 1L,
    null_skip_max = 1L, candidate_threshold = 0,
    candidate_lifetime = 0L, skip_nulls_burnin_only = FALSE,
    progress_every = 0L
  )))
  expect_equal(raw_rc$marker$bm, raw_r$marker$bm, tolerance = 1e-12)
  expect_equal(raw_rc$marker$dm, raw_r$marker$dm, tolerance = 1e-12)
  expect_equal(raw_rc$component$prob[[1L]], raw_r$component$prob[[1L]], tolerance = 1e-12)
})
