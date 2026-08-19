# Focused tests for the public gsim() contract.
# Run from the directory containing these files with:
#   testthat::test_file("test-gsim.R")

if (!exists("gsim", mode = "function")) {
  source("gsim_internal.R", local = FALSE)
  source("gsim.R", local = FALSE)
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
    getG_fun = fake_getG,
    compute_sumstats = FALSE,
    return_genotypes = FALSE
  )

  testthat::expect_setequal(unique(requested$rsids), sim$causal_rsids)
  testthat::expect_equal(length(unique(requested$rsids)), 7L)
  testthat::expect_false(any(names(sim) == "W_causal"))
  testthat::expect_equal(dim(sim$Y), c(120L, 1L))
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
