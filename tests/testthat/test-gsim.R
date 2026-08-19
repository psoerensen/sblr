# Focused tests for the public gsim() contract.
# Run from the directory containing these files with:
#   testthat::test_file("test-gsim.R")

if (!exists("gsim", mode = "function")) {
  source("gsim_internal.R", local = FALSE)
  source("gsim.R", local = FALSE)
}

.gsim_without_multiplier_additions <- function(x) {
  x$marker_multipliers <- NULL
  x$settings$marker_multipliers <- NULL
  x
}

testthat::test_that("in-memory simulation is reproducible and internally exact", {
  set.seed(11)
  W <- matrix(rbinom(600 * 80, 2, 0.3), 600, 80)
  colnames(W) <- paste0("m", seq_len(ncol(W)))
  rownames(W) <- paste0("id", seq_len(nrow(W)))

  x <- gsim(
    W = W,
    architecture = "bayesr",
    n_causal = 12L,
    h2 = 0.5,
    seed = 123,
    return_genotypes = TRUE
  )
  y <- gsim(
    W = W,
    architecture = "bayesr",
    n_causal = 12L,
    h2 = 0.5,
    seed = 123,
    return_genotypes = TRUE
  )

  testthat::expect_equal(x$B, y$B)
  testthat::expect_equal(x$Y, y$Y)
  testthat::expect_equal(nrow(x$causal), 12L)
  testthat::expect_equal(
    unname(x$G),
    unname(x$W_causal %*% x$B_causal),
    tolerance = 1e-12
  )
  testthat::expect_lt(x$exactness$max_y_minus_g_plus_e, 1e-12)
  testthat::expect_lt(abs(x$vg_observed - 1), 1e-10)
})

testthat::test_that("mixed annotations produce valid SBayesRC probabilities", {
  set.seed(22)
  W <- matrix(rbinom(500 * 100, 2, 0.25), 500, 100)
  colnames(W) <- paste0("m", seq_len(ncol(W)))
  rownames(W) <- paste0("id", seq_len(nrow(W)))
  A <- cbind(
    binary = rep(c(0, 1), each = 50),
    continuous = seq(-2, 2, length.out = 100)
  )
  rownames(A) <- colnames(W)
  alpha <- matrix(
    c(1.2, 0, 0,
      0.4, 0, 0),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(colnames(A), paste0("stick", 1:3))
  )

  sim <- gsim(
    W = W,
    A = A,
    architecture = "bayesr",
    annotation_model = "sbayesrc",
    alpha = alpha,
    n_causal = 20L,
    seed = 456
  )

  P <- sim$marker_probabilities
  active <- 1 - P[, 1L]
  testthat::expect_equal(
    unname(rowSums(P)),
    rep(1, nrow(P)),
    tolerance = 1e-12
  )
  testthat::expect_true(all(P >= 0 & P <= 1))
  testthat::expect_gt(mean(active[A[, "binary"] == 1]),
                      mean(active[A[, "binary"] == 0]))
  testthat::expect_equal(sim$annotation_types,
                         c(binary = "binary", continuous = "continuous"))
})

testthat::test_that("unit marker multipliers are an exact RNG-neutral reduction", {
  set.seed(101)
  W <- matrix(rbinom(180 * 36, 2, 0.3), 180, 36)
  colnames(W) <- paste0("m", seq_len(ncol(W)))
  rownames(W) <- paste0("id", seq_len(nrow(W)))
  unit <- stats::setNames(rep(1, ncol(W)), colnames(W))
  common <- list(
    W = W, architecture = "bayesr", n_causal = 9L, h2 = 0.45,
    seed = 2468, return_genotypes = TRUE, compute_sumstats = TRUE
  )

  default <- do.call(gsim, common)
  rng_default <- .Random.seed
  supplied <- do.call(gsim, c(common, list(marker_multipliers = unit)))
  rng_supplied <- .Random.seed

  testthat::expect_identical(
    .gsim_without_multiplier_additions(default),
    .gsim_without_multiplier_additions(supplied)
  )
  testthat::expect_identical(rng_default, rng_supplied)
  testthat::expect_identical(default$marker_multipliers, unit)
  testthat::expect_identical(supplied$marker_multipliers, unit)
  testthat::expect_identical(default$settings$marker_multipliers$policy, "unit")
  testthat::expect_identical(
    supplied$settings$marker_multipliers$policy, "supplied"
  )
})

testthat::test_that("marker multipliers align once and reject invalid inputs", {
  set.seed(102)
  W <- matrix(rbinom(100 * 18, 2, 0.28), 100, 18)
  colnames(W) <- paste0("m", seq_len(ncol(W)))
  rownames(W) <- paste0("id", seq_len(nrow(W)))
  canonical <- stats::setNames(
    exp(seq(-0.5, 0.5, length.out = ncol(W))), colnames(W)
  )
  reordered <- canonical[rev(names(canonical))]
  common <- list(
    W = W, architecture = "bayesr", n_causal = 5L,
    seed = 975, scale_effects = FALSE
  )

  sim <- do.call(gsim, c(common, list(marker_multipliers = reordered)))
  testthat::expect_identical(sim$marker_multipliers, canonical)
  testthat::expect_identical(
    sim$settings$marker_multipliers$alignment, "canonical_marker_order"
  )

  duplicated <- unknown <- nonpositive <- nonfinite <- canonical
  names(duplicated)[1L] <- names(duplicated)[2L]
  names(unknown)[1L] <- "unknown_marker"
  nonpositive[1L] <- 0
  nonfinite[1L] <- Inf
  invalid <- list(
    unnamed = unname(canonical),
    missing = canonical[-1L],
    duplicated = duplicated,
    unknown = unknown,
    nonpositive = nonpositive,
    nonfinite = nonfinite
  )
  for (value in invalid) {
    testthat::expect_error(
      do.call(gsim, c(common, list(marker_multipliers = value))),
      "marker_multipliers"
    )
  }
})

testthat::test_that("marker multipliers change only active effect scales", {
  set.seed(103)
  W <- matrix(rbinom(160 * 30, 2, 0.31), 160, 30)
  colnames(W) <- paste0("m", seq_len(ncol(W)))
  rownames(W) <- paste0("id", seq_len(nrow(W)))
  q <- stats::setNames(exp(seq(-0.8, 0.8, length.out = ncol(W))), colnames(W))
  common <- list(
    W = W, architecture = "bayesr", nt = 2L, rg = 0.25,
    n_causal = 8L, seed = 8642, scale_effects = FALSE
  )

  unit <- do.call(gsim, common)
  rng_unit <- .Random.seed
  scaled <- do.call(gsim, c(common, list(marker_multipliers = q)))
  rng_scaled <- .Random.seed
  active <- unit$component != 1L
  expected <- sweep(
    unit$B[active, , drop = FALSE], 1L, sqrt(q[active]), "*"
  )

  testthat::expect_identical(unit$causal_rsids, scaled$causal_rsids)
  testthat::expect_identical(unit$component, scaled$component)
  testthat::expect_identical(
    unit$marker_probabilities, scaled$marker_probabilities
  )
  testthat::expect_equal(scaled$B[active, , drop = FALSE], expected,
                         tolerance = 1e-15)
  testthat::expect_true(all(scaled$B[!active, , drop = FALSE] == 0))
  testthat::expect_identical(rng_unit, rng_scaled)
  testthat::expect_identical(scaled$marker_multipliers, q)
  testthat::expect_equal(
    scaled$settings$marker_multipliers$geometric_mean,
    exp(mean(log(q))), tolerance = 1e-15
  )
  testthat::expect_identical(
    scaled$settings$marker_multipliers$n_markers, length(q)
  )
  testthat::expect_identical(scaled$settings$marker_multipliers$minimum, min(q))
  testthat::expect_identical(scaled$settings$marker_multipliers$maximum, max(q))
  testthat::expect_identical(scaled$settings$marker_multipliers$all_ones, FALSE)
})

testthat::test_that("membership and variance truth remain independent", {
  set.seed(104)
  W <- matrix(rbinom(140 * 24, 2, 0.27), 140, 24)
  colnames(W) <- paste0("m", seq_len(ncol(W)))
  rownames(W) <- paste0("id", seq_len(nrow(W)))
  A <- cbind(binary = rep(c(0, 1), each = 12),
             continuous = seq(-1, 1, length.out = 24))
  rownames(A) <- colnames(W)
  alpha <- matrix(
    c(0.8, 0, 0, 0.25, 0, 0), nrow = 2L, byrow = TRUE,
    dimnames = list(colnames(A), paste0("stick", 1:3))
  )
  q <- stats::setNames(exp(seq(-0.4, 0.4, length.out = 24)), colnames(W))
  common <- list(
    W = W, A = A, architecture = "bayesr",
    annotation_model = "sbayesrc", alpha = alpha,
    n_causal = 7L, seed = 1357, scale_effects = FALSE
  )

  inclusion_only <- do.call(gsim, common)
  combined <- do.call(gsim, c(common, list(marker_multipliers = q)))

  testthat::expect_identical(
    inclusion_only$marker_probabilities, combined$marker_probabilities
  )
  testthat::expect_identical(
    inclusion_only$continuation_probabilities,
    combined$continuation_probabilities
  )
  testthat::expect_identical(inclusion_only$component, combined$component)
  testthat::expect_identical(combined$marker_multipliers, q)
  testthat::expect_identical(
    combined$settings$annotation_model, "sbayesrc"
  )
  testthat::expect_identical(
    combined$settings$marker_multipliers$policy, "supplied"
  )
})

testthat::test_that("Glist phenotype generation requests causal markers only", {
  set.seed(33)
  W_store <- matrix(rbinom(120 * 60, 2, 0.3), 120, 60)
  colnames(W_store) <- paste0("m", seq_len(ncol(W_store)))
  rownames(W_store) <- paste0("id", seq_len(nrow(W_store)))
  Glist <- list(
    ids = rownames(W_store),
    rsids = list(`1` = colnames(W_store)[1:30],
                 `2` = colnames(W_store)[31:60])
  )
  q <- stats::setNames(
    seq(0.5, 1.5, length.out = ncol(W_store)), colnames(W_store)
  )

  requested <- new.env(parent = emptyenv())
  requested$rsids <- character(0)
  fake_getG <- function(Glist, rsids, ids, chr = NULL,
                        impute = TRUE, scale = FALSE) {
    requested$rsids <- c(requested$rsids, rsids)
    W_store[ids, rsids, drop = FALSE]
  }

  sim <- gsim(
    Glist = Glist,
    architecture = "bayesc",
    n_causal = 7L,
    h2 = 0.4,
    seed = 789,
    marker_multipliers = q[rev(names(q))],
    getG_fun = fake_getG,
    compute_sumstats = FALSE,
    return_genotypes = FALSE
  )

  testthat::expect_setequal(unique(requested$rsids), sim$causal_rsids)
  testthat::expect_equal(length(unique(requested$rsids)), 7L)
  testthat::expect_false(any(names(sim) == "W_causal"))
  testthat::expect_equal(dim(sim$Y), c(120L, 1L))
  testthat::expect_identical(sim$marker_multipliers, q)
})


testthat::test_that("gsim works with the qgg PLINK fixture", {
  testthat::skip_if_not_installed("qgg")

  bedfiles <- system.file(
    "extdata",
    paste0("sample_chr", 1:2, ".bed"),
    package = "qgg"
  )
  bimfiles <- sub("\\.bed$", ".bim", bedfiles)
  famfiles <- sub("\\.bed$", ".fam", bedfiles)

  Glist <- suppressMessages(
    qgg::gprep(
      study = "gsim-test",
      bedfiles = bedfiles,
      bimfiles = bimfiles,
      famfiles = famfiles
    )
  )

  sim <- gsim(
    Glist = Glist,
    architecture = "bayesr",
    n_causal = 10L,
    seed = 1,
    return_genotypes = TRUE
  )

  testthat::expect_equal(nrow(sim$Y), 489L)
  testthat::expect_equal(sim$settings$m, 2000L)
  testthat::expect_equal(sim$settings$n_causal, 10L)
  testthat::expect_equal(
    unname(sim$G),
    unname(sim$W_causal %*% sim$B_causal),
    tolerance = 1e-12
  )
})
