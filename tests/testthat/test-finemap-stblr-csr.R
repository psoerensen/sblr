if (!exists("make_credible_sets_from_ld", mode = "function")) {
  credible_sets_path <- if (file.exists("R/credible_sets.R")) {
    "R/credible_sets.R"
  } else {
    "../../R/credible_sets.R"
  }
  source(credible_sets_path)
}

if (!exists("finemap_stblr_csr", mode = "function")) {
  finemap_path <- if (file.exists("R/finemap-stblr-csr.R")) {
    "R/finemap-stblr-csr.R"
  } else {
    "../../R/finemap-stblr-csr.R"
  }
  source(finemap_path)
}

test_that("fine-mapping sets are cleaned and aligned to marker order", {
  marker_names <- paste0("m", 1:5)
  sets <- list(region = c("m5", "missing", "m2", "m2"), empty = "absent")

  expect_warning(
    cleaned <- .stblr_clean_finemap_sets(sets, marker_names),
    "Dropping empty"
  )

  expect_named(cleaned, "region")
  expect_equal(cleaned$region, c("m2", "m5"))
})

test_that("global parameters are extracted from traces after burn-in", {
  fit <- list(
    input = list(nburn = 2),
    ves = matrix(c(1, 2, 10, 20, 100, 200), ncol = 2, byrow = TRUE),
    vbs = matrix(c(3, 4, 30, 40, 300, 400), ncol = 2, byrow = TRUE),
    pis = matrix(c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6), ncol = 2, byrow = TRUE)
  )
  colnames(fit$ves) <- colnames(fit$vbs) <- colnames(fit$pis) <- c("t1", "t2")

  expect_equal(.stblr_get_global_parameter(fit, "ve", 2), 200)
  expect_equal(.stblr_get_global_parameter(fit, "vb", 1), 300)
  expect_equal(.stblr_get_global_parameter(fit, "pi", 2), 0.6)
  expect_equal(.stblr_get_global_parameter(fit, "pi", 2, override = 0.05), 0.05)
})

test_that("local marker tables aggregate PIPs and effects across runs", {
  fits <- list(
    list(
      dm = matrix(c(0.2, 0.7, 0.1), ncol = 1),
      bm = matrix(c(1, 2, 3), ncol = 1)
    ),
    list(
      dm = matrix(c(0.4, 0.5, 0.2), ncol = 1),
      bm = matrix(c(2, 4, 6), ncol = 1)
    )
  )
  markers <- c("m1", "m2", "m3")
  map <- data.frame(
    marker = markers,
    chr = c(1, 1, 1),
    pos = c(10, 20, 30),
    index = 1:3,
    stringsAsFactors = FALSE
  )

  out <- .stblr_aggregate_finemap_runs(fits, "locusA", "trait1", markers, map)

  expect_equal(out$markers$pip_mean, c(0.3, 0.6, 0.15))
  expect_equal(out$markers$bm_mean, c(1.5, 3, 4.5))
  expect_equal(out$summary$lead_marker, "m2")
  expect_equal(out$summary$total_pip, 1.05)
})

test_that("credible sets use aggregated local PIPs", {
  fits <- list(
    list(
      dm = matrix(c(0.1, 0.8, 0.1), ncol = 1),
      bm = matrix(0, nrow = 3, ncol = 1)
    ),
    list(
      dm = matrix(c(0.9, 0.2, 0.1), ncol = 1),
      bm = matrix(0, nrow = 3, ncol = 1)
    )
  )
  markers <- c("m1", "m2", "m3")
  agg <- .stblr_aggregate_finemap_runs(fits, "locusA", "trait1", markers)
  LD <- diag(3)
  rownames(LD) <- colnames(LD) <- markers

  cs <- make_credible_sets_from_ld(
    pip = stats::setNames(agg$markers$pip_mean, agg$markers$marker),
    LD = LD,
    coverage = 0.45,
    min_r2 = 0,
    pip_cutoff = 0.001
  )

  expect_equal(unname(cs$pip["m1"]), 0.5)
  expect_equal(unname(cs$pip["m2"]), 0.5)
  expect_true(cs$summary$lead_marker[1] %in% c("m1", "m2"))
})

test_that("credible-set thresholds do not remove markers before local fitting", {
  marker_names <- paste0("m", 1:4)
  sets <- list(locusA = marker_names)

  cleaned <- .stblr_clean_finemap_sets(sets, marker_names)

  expect_equal(cleaned$locusA, marker_names)

  low_pip <- c(m1 = 0.9, m2 = 0.0001, m3 = 0.0001, m4 = 0.0001)
  LD <- diag(4)
  rownames(LD) <- colnames(LD) <- marker_names
  cs <- make_credible_sets_from_ld(
    low_pip,
    LD,
    coverage = 0.8,
    min_r2 = 0.5,
    pip_cutoff = 0.001
  )

  expect_equal(length(cleaned$locusA), 4)
  expect_equal(unname(cs$pip[c("m2", "m3", "m4")]), c(0, 0, 0))
})
