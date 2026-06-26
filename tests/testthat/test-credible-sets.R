if (!exists("make_credible_sets_from_ld", mode = "function")) {
  credible_sets_path <- if (file.exists("R/credible_sets.R")) {
    "R/credible_sets.R"
  } else {
    "../../R/credible_sets.R"
  }
  source(credible_sets_path)
}

test_that("make_credible_sets_from_ld builds sets with both methods", {
  pip <- c(m1 = 0.60, m2 = 0.25, m3 = 0.08, m4 = 0.004, m5 = 0.0005)
  LD <- diag(5)
  LD[1, 2] <- LD[2, 1] <- sqrt(0.8)
  LD[1, 3] <- LD[3, 1] <- sqrt(0.6)
  LD[4, 5] <- LD[5, 4] <- sqrt(0.9)
  rownames(LD) <- colnames(LD) <- names(pip)

  out <- make_credible_sets_from_ld(
    pip = pip,
    LD = LD,
    coverage = 0.80,
    min_r2 = 0.5,
    pip_cutoff = 0.001,
    method = "pip"
  )

  expect_s3_class(out$summary, "data.frame")
  expect_gt(nrow(out$summary), 0)
  expect_gt(length(out$sets), 0)
  expect_equal(unname(out$pip["m5"]), 0)
  expect_true(out$summary$cs_pip[1] >= 0.80)
  expect_true(all(out$sets[[1]]$r2_to_lead >= 0.5))

  out_ld <- make_credible_sets_from_ld(
    pip = pip,
    LD = LD,
    coverage = 0.80,
    min_r2 = 0.5,
    method = "ld_pip"
  )
  expect_gt(nrow(out_ld$summary), 0)
  expect_equal(out_ld$summary$method[1], "ld_pip")
})

test_that("make_credible_sets_from_ld aligns by marker names", {
  pip <- c(m3 = 0.10, m1 = 0.60, m2 = 0.30)
  LD <- diag(3)
  rownames(LD) <- colnames(LD) <- c("m1", "m2", "m3")

  out <- make_credible_sets_from_ld(pip, LD, coverage = 0.50, min_r2 = 0)

  expect_equal(names(out$pip), c("m3", "m1", "m2"))
  expect_equal(unname(out$pip["m1"]), 0.60)
})

test_that("incomplete credible sets are not returned by default", {
  pip <- c(a = 0.1, b = 0.05, c = 0.02)
  LD <- diag(3)
  rownames(LD) <- colnames(LD) <- names(pip)

  cs <- make_credible_sets_from_ld(pip, LD, coverage = 0.95)

  expect_equal(nrow(cs$summary), 0)
})

test_that("incomplete credible sets can be returned when requested", {
  pip <- c(a = 0.1, b = 0.05, c = 0.02)
  LD <- diag(3)
  rownames(LD) <- colnames(LD) <- names(pip)

  cs <- make_credible_sets_from_ld(
    pip,
    LD,
    coverage = 0.95,
    allow_incomplete = TRUE
  )

  expect_true(nrow(cs$summary) > 0)
  expect_true(all(cs$summary$cs_pip < 0.95))
})

test_that("complete lead PIP suppresses extra incomplete singleton sets", {
  pip <- c(a = 1, b = 0.01, c = 0.005)
  LD <- diag(3)
  rownames(LD) <- colnames(LD) <- names(pip)

  cs <- make_credible_sets_from_ld(
    pip,
    LD,
    coverage = 0.95,
    pip_cutoff = 0.001
  )

  expect_equal(nrow(cs$summary), 1)
  expect_equal(cs$summary$lead_marker, "a")
  expect_equal(cs$summary$cs_pip, 1)
})

test_that("ld_neighborhood removal suppresses future LD-neighborhood sets", {
  pip <- c(a = 0.8, b = 0.2, c = 0.01)
  LD <- matrix(sqrt(0.6), nrow = 3, ncol = 3)
  diag(LD) <- 1
  rownames(LD) <- colnames(LD) <- names(pip)

  cs <- make_credible_sets_from_ld(
    pip,
    LD,
    coverage = 0.95,
    min_r2 = 0.5,
    pip_cutoff = 0.001,
    remove = "ld_neighborhood"
  )

  expect_equal(nrow(cs$summary), 1)
  expect_equal(cs$summary$lead_marker, "a")
  expect_equal(cs$sets[[1]]$marker, c("a", "b"))
  expect_equal(cs$summary$cs_pip, 1)
})

test_that("make_stblr_credible_sets works with predefined sets and dense LD", {
  fit <- list(dm = matrix(
    c(0.55, 0.20, 0.10, 0.03, 0.002),
    ncol = 1,
    dimnames = list(paste0("m", 1:5), "trait1")
  ))
  LD <- diag(5)
  LD[1, 2] <- LD[2, 1] <- sqrt(0.7)
  LD[2, 3] <- LD[3, 2] <- sqrt(0.6)
  rownames(LD) <- colnames(LD) <- paste0("m", 1:5)

  out <- make_stblr_credible_sets(
    fit = fit,
    LD = LD,
    sets = list(regionA = c("m1", "m2", "m3")),
    trait = "trait1",
    coverage = 0.70,
    min_r2 = 0.5
  )

  expect_s3_class(out$summary, "data.frame")
  expect_gt(nrow(out$summary), 0)
  expect_named(out$sets, "regionA")
  expect_s3_class(out$loci, "data.frame")
  expect_equal(out$trait, "trait1")
})

test_that("automatic loci use Glist marker maps", {
  fit <- list(dm = matrix(
    c(0.50, 0.20, 0.005, 0.40, 0.01),
    ncol = 1,
    dimnames = list(paste0("m", 1:5), "trait1")
  ))
  Glist <- list(
    rsids = list(c("m1", "m2", "m3"), c("m4", "m5")),
    pos = list(c(100, 200, 5000), c(1000, 1200)),
    sparseLD = list(chr = c(1L, 2L), cls = list(1:3, 1:2))
  )

  map <- .stblr_marker_map_from_Glist(Glist, fit = fit)
  loci <- .stblr_define_loci(
    pip = fit$dm[, 1],
    map = map,
    locus_pip_cutoff = 0.05,
    max_locus_distance = 500,
    max_loci = Inf
  )

  expect_s3_class(loci, "data.frame")
  expect_equal(nrow(loci), 2)
  expect_equal(loci$lead_marker, c("m1", "m4"))
})

test_that("Glist marker map uses biological chromosomes from one component", {
  fit <- list(dm = matrix(
    c(0.10, 0.20, 0.50, 0.05),
    ncol = 1,
    dimnames = list(c("1:100", "1:200", "2:100", "2:200"), "trait1")
  ))
  Glist <- list(
    rsids = list(c("1:100", "1:200", "2:100", "2:200")),
    chr = list(c(1, 1, 2, 2)),
    pos = list(c(100, 200, 100, 200)),
    sparseLD = list(chr = list(1), cls = list(1:4))
  )

  map <- .stblr_marker_map_from_Glist(Glist, fit = fit)
  loci <- .stblr_define_loci(
    pip = fit$dm[, 1],
    map = map,
    locus_pip_cutoff = 0.40,
    max_locus_distance = 500,
    max_loci = Inf
  )

  expect_equal(map$marker, c("1:100", "1:200", "2:100", "2:200"))
  expect_equal(map$chr, c(1, 1, 2, 2))
  expect_equal(map$pos, c(100, 200, 100, 200))
  expect_equal(map$index, 1:4)
  expect_equal(loci$chr, 2)
  expect_equal(loci$lead_marker, "2:100")
})

test_that("Glist marker map falls back to component index without chr", {
  Glist <- list(
    rsids = list(c("1:100", "1:200"), c("2:100", "2:200")),
    pos = list(c(100, 200), c(100, 200)),
    sparseLD = list(chr = list(1, 2), cls = list(1:2, 1:2))
  )

  map <- .stblr_marker_map_from_Glist(Glist)

  expect_equal(map$marker, c("1:100", "1:200", "2:100", "2:200"))
  expect_equal(map$chr, c(1, 1, 2, 2))
  expect_equal(map$pos, c(100, 200, 100, 200))
  expect_equal(map$index, 1:4)
})

test_that("sparse CSR extraction densifies selected markers", {
  csr <- list(
    indptr = c(0, 2, 3, 4, 4, 4),
    indices = c(2, 3, 3, 5),
    values = c(0.8, 0.4, 0.7, -0.6)
  )

  LD <- .extract_sparseLD_region_dense(
    csr,
    idx = c(1L, 2L, 3L),
    marker_names = c("m1", "m2", "m3")
  )

  expect_equal(dim(LD), c(3L, 3L))
  expect_equal(unname(diag(LD)), rep(1, 3))
  expect_equal(LD["m1", "m2"], 0.8)
  expect_equal(LD["m2", "m3"], 0.7)
  expect_equal(LD["m1", "m3"], 0.4)
  expect_equal(LD["m3", "m1"], 0.4)
})

test_that("sparse CSR extraction warns for large empty off-diagonal LD", {
  csr <- list(
    indptr = rep(0, 12),
    indices = integer(),
    values = numeric()
  )

  expect_warning(
    LD <- .extract_sparseLD_region_dense(
      csr,
      idx = 1:11,
      marker_names = paste0("m", 1:11)
    ),
    "LD extraction or sparseLD indexing may be wrong"
  )
  expect_equal(dim(LD), c(11L, 11L))
  expect_equal(unname(diag(LD)), rep(1, 11))
})
