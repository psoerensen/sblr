make_stblr_low_rank_fixture <- function() {
  dosage <- rbind(
    c(0, 1, 2, 0, 1, 2, 1, 0),
    c(2, 0, 1, 2, 0, 1, 2, 1),
    c(0, 0, 1, 1, 2, 2, 1, 0),
    c(2, 2, 1, 1, 0, 0, 1, 2)
  )
  bed <- tempfile(fileext = ".bed")
  blr_block_fixture_write_bed(bed, dosage)
  af <- c(0.30, 0.35, 0.40, 0.25)
  marker <- paste0("rs", seq_len(4L))
  score <- c(1.2, -0.6, 0.5, 0.9)
  descriptor <- list(
    bed_files = bed, n_bed = 8L, cls = list(seq_len(4L)), rows = NULL,
    af = af, block_start = c(0L, 2L)
  )
  inspection <- sblr:::stblr_block_low_rank_contract_internal(
    descriptor$bed_files, descriptor$n_bed, descriptor$cls, descriptor$rows,
    descriptor$af, descriptor$block_start, matrix(score, nrow = 1L),
    c(0.1, -0.2, 0.05, 0.3), 0.999999
  )
  Glist <- list(
    n = 8L, ids = paste0("id", seq_len(8L)), bedfiles = bed,
    rsids = list(marker), rsidsLD = list(marker), chr = list(rep(1L, 4L)),
    pos = list(seq_len(4L) * 100), af = list(af), maf = list(af)
  )
  stats <- list(
    wy = list(T1 = stats::setNames(score, marker)),
    ww = list(T1 = stats::setNames(as.numeric(inspection$diagonal), marker)),
    yy = stats::setNames(20, "T1"), n = 8L, m = 4L,
    bed_files = bed, cls = list(seq_len(4L)), rows = seq_len(8L),
    af = list(af), marker_names = marker, trait_names = "T1"
  )
  list(descriptor = descriptor, inspection = inspection, Glist = Glist,
       stats = stats)
}

test_that("retained factors satisfy the reduced operator identities", {
  fixture <- make_stblr_low_rank_fixture()
  x <- fixture$inspection
  beta <- c(0.1, -0.2, 0.05, 0.3)
  for (block in seq_along(x$factor)) {
    Q <- x$factor[[block]]
    w <- as.numeric(x$transformed_score[[block]][1, ])
    index <- if (block == 1L) 1:2 else 3:4
    expect_equal(as.numeric(crossprod(Q, w)),
                 as.numeric(x$projected_score[1, index]), tolerance = 2e-6)
    expect_equal(colSums(Q * Q), as.numeric(x$diagonal[index]), tolerance = 2e-6)
    offset <- x$residual_offset[block] + seq_len(nrow(Q))
    expect_equal(as.numeric(x$residual[offset]),
                 as.numeric(w - Q %*% beta[index]), tolerance = 2e-6)
    expect_equal(as.numeric(crossprod(Q, x$residual[offset])),
                 as.numeric(x$marker_residual[index]), tolerance = 2e-6)
  }
  Qbeta_norm <- sum(vapply(seq_along(x$factor), function(block) {
    index <- if (block == 1L) 1:2 else 3:4
    sum((x$factor[[block]] %*% beta[index])^2)
  }, numeric(1)))
  expect_equal(x$quadratic_form, Qbeta_norm, tolerance = 2e-6)
  expect_equal(x$projected_score_dot,
               sum(beta * as.numeric(x$projected_score)), tolerance = 2e-6)
  expect_true(all(x$diagnostics$retained_mass_fraction > 0.999999))
})

test_that("full-rank retained factors match reconstructed dense blocks", {
  fixture <- make_stblr_low_rank_fixture()
  descriptor <- fixture$descriptor
  dense <- sblr:::stblr_block_eigen_contract_internal(
    descriptor$bed_files, descriptor$n_bed, descriptor$cls, descriptor$rows,
    descriptor$af, descriptor$block_start,
    matrix(c(1.2, -0.6, 0.5, 0.9), nrow = 1L), rep(0, 4L),
    "hard_truncate", 1e-12, 0, ""
  )
  low <- fixture$inspection
  for (block in seq_along(low$factor)) {
    Q <- low$factor[[block]]
    packed <- dense$packed_upper_triangle[[block]]
    oracle <- matrix(0, ncol(Q), ncol(Q)); cursor <- 1L
    for (i in seq_len(ncol(Q))) for (j in i:ncol(Q)) {
      oracle[i, j] <- oracle[j, i] <- packed[cursor]
      cursor <- cursor + 1L
    }
    expect_equal(crossprod(Q), oracle, tolerance = 3e-5)
  }
  expect_equal(as.numeric(low$projected_score),
               as.numeric(dense$transformed_wy), tolerance = 3e-5)
})

test_that("all scalar block-eigen models run with reduced residual state", {
  fixture <- make_stblr_low_rank_fixture()
  common <- list(
    stats = fixture$stats, Glist = fixture$Glist, block_start = c(1L, 3L),
    representation = "low_rank", eigen_prop = 0.995,
    updateB = FALSE, updateE = FALSE, nit = 4L, nburn = 1L,
    nchains = 1L, ncores = 1L, seed = 19L
  )
  bayesc <- do.call(stblr_block_eigen, c(common, list(
    method = "sbayesc", updatePi = FALSE, pi_init = 0.5,
    pi_prior_mean = 0.5, pi_prior_strength = 2
  )))
  bayesr <- do.call(stblr_block_eigen, c(common, list(
    method = "sbayesr", updatePi = FALSE
  )))
  annotation <- cbind(intercept = 1, informative = c(0, 0, 1, 1))
  sbayesrc <- do.call(stblr_block_eigen, c(common, list(
    method = "sbayesrc", annotation = annotation, updateAlpha = FALSE
  )))
  for (fit in list(bayesc, bayesr, sbayesrc)) {
    expect_identical(fit$input$operator_representation, "low_rank")
    expect_identical(fit$input$operator_contract, "block_low_rank_v1")
    expect_true(all(is.finite(fit$bm)))
    expect_true("r" %in% names(fit))
  }
  diagnostic <- bayesc$input$eigen_diagnostics
  expect_true(is.list(diagnostic))
  expect_true(all(c("blocks", "operator_contract", "operator_representation",
                    "operator_scale_contract", "eigen_policy", "eigen_prop",
                    "build") %in% names(diagnostic)))
  expect_true(all(c("block_start", "block_size", "positive_rank",
                    "retained_rank", "retained_mass_fraction",
                    "eigenvalue_tolerance") %in% names(diagnostic$blocks)))
  expect_gt(diagnostic$build$operator_storage_bytes, 0)
  expect_gt(diagnostic$build$chain_residual_storage_bytes, 0)
})

test_that("projected residual variance updates remain finite", {
  fixture <- make_stblr_low_rank_fixture()
  fit <- stblr_block_eigen(
    fixture$stats, fixture$Glist, c(1L, 3L), method = "sbayesc",
    representation = "low_rank", eigen_prop = 0.995,
    pi_init = 0.5, pi_prior_mean = 0.5, pi_prior_strength = 2,
    updateB = FALSE, updateE = TRUE, updatePi = FALSE,
    nit = 8L, nburn = 2L, seed = 71L, ncores = 1L
  )
  expect_true(all(is.finite(fit$ves)))
  expect_true(all(fit$ves > 0))
  expect_true(all(is.finite(fit$vgs)))
})

test_that("the public default resolves identically to explicit low rank", {
  fixture <- make_stblr_low_rank_fixture()
  args <- list(
    stats = fixture$stats, Glist = fixture$Glist, block_start = c(1L, 3L),
    method = "sbayesc", pi_init = 0.5, pi_prior_mean = 0.5,
    pi_prior_strength = 2, updateB = FALSE, updateE = FALSE,
    updatePi = FALSE, nit = 5L, nburn = 1L, seed = 808L, ncores = 1L
  )
  default <- do.call(stblr_block_eigen, args)
  explicit <- do.call(stblr_block_eigen, c(args, list(
    representation = "low_rank", eigen_policy = "cumulative_positive_mass",
    eigen_prop = 0.995
  )))
  for (field in c("bm", "dm", "vbs", "vgs", "ves", "vle", "vld"))
    expect_equal(default[[field]], explicit[[field]], tolerance = 0)
  expect_identical(default$input$operator_representation, "low_rank")
  expect_identical(explicit$input$operator_representation, "low_rank")
})

test_that("representation and eigen policy combinations fail before native execution", {
  expect_error(
    stblr_block_eigen(list(), NULL, 1L, representation = "low_rank",
                      eigen_policy = "absolute_threshold"),
    "Unsupported representation/eigen_policy"
  )
  expect_error(
    stblr_block_eigen(list(), NULL, 1L, representation = "low_rank",
                      eigen_prop = 1),
    "strictly between 0 and 1"
  )
})
