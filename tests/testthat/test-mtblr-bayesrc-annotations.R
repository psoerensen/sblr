test_that("MT BayesRC annotation alignment is exact and invariant", {
  x <- .mt_bayesrc_fixture(); on.exit(.mt_bayesrc_cleanup(x), add = TRUE)
  common <- .mt_bayesrc_common(); common$annotations <- x$annotations
  base <- do.call(mtblr_csr, c(list(stats = x$stats, ld_prefix = x$prefix,
    ld_metadata = x$ld_metadata), common))
  permuted <- common
  permuted$annotations <- x$annotations[c(3, 1, 2), , drop = FALSE]
  reordered <- do.call(mtblr_csr, c(list(stats = x$stats,
    ld_prefix = x$prefix, ld_metadata = x$ld_metadata), permuted))
  for (field in c("bm", "dm", "component_probabilities", "vgs", "vle", "vld"))
    expect_equal(reordered[[field]], base[[field]], tolerance = 0)
  bad <- common; rownames(bad$annotations)[1L] <- "missing"
  expect_error(do.call(mtblr_csr, c(list(stats = x$stats,
    ld_prefix = x$prefix, ld_metadata = x$ld_metadata), bad)),
    "missing annotation")
  bad <- common; rownames(bad$annotations)[2L] <- rownames(bad$annotations)[1L]
  expect_error(do.call(mtblr_csr, c(list(stats = x$stats,
    ld_prefix = x$prefix, ld_metadata = x$ld_metadata), bad)), "unique")
  bad <- common; colnames(bad$annotations) <- c("x", "x")
  expect_error(do.call(mtblr_csr, c(list(stats = x$stats,
    ld_prefix = x$prefix, ld_metadata = x$ld_metadata), bad)), "unique")
  bad <- common; rownames(bad$annotations) <- NULL
  expect_error(do.call(mtblr_csr, c(list(stats = x$stats,
    ld_prefix = x$prefix, ld_metadata = x$ld_metadata), bad)),
    "explicit marker IDs")
  bad <- common
  bad$annotations <- rbind(bad$annotations,
    extra = c(intercept = 1, coding = 0))
  expect_error(do.call(mtblr_csr, c(list(stats = x$stats,
    ld_prefix = x$prefix, ld_metadata = x$ld_metadata), bad)),
    "extra rows")
})

test_that("MT public method matrix requires annotations only for BayesRC", {
  x <- .mt_bayesrc_fixture(); on.exit(.mt_bayesrc_cleanup(x), add = TRUE)
  expect_error(mtblr_csr(x$stats, ld_prefix = x$prefix,
    ld_metadata = x$ld_metadata, method = "bayesrc"),
    "summary statistics")
  expect_error(mtblr_bed(x$fixture$y, x$fixture$Glist, method = "sbayesrc"),
    "summary statistics")
  expect_error(mtblr_csr(x$stats, ld_prefix = x$prefix,
    ld_metadata = x$ld_metadata, method = "sbayesrc"), "annotations")
  expect_error(mtblr_csr(x$stats, ld_prefix = x$prefix,
    ld_metadata = x$ld_metadata, method = "sbayesr",
    annotations = x$annotations), "Annotation controls")
})

test_that("MAF annotation overlap is explicit and advisory once", {
  x <- .mt_bayesrc_fixture(); on.exit(.mt_bayesrc_cleanup(x), add = TRUE)
  A <- cbind(intercept = 1, maf = c(.1, .2, .3))
  rownames(A) <- x$stats$marker_names
  args <- .mt_bayesrc_common(maf_effect_s = -1)
  args$annotations <- A
  expect_warning(fit <- do.call(mtblr_csr, c(list(stats = x$stats,
    ld_prefix = x$prefix, ld_metadata = x$ld_metadata), args)),
    "MAF-derived annotations")
  expect_true(fit$model_parameters$annotations$maf_annotation_overlap)
  expect_identical(fit$data$annotation_marker_alignment_status,
                   "exact_marker_id_match")
})

test_that("MT BayesRC annotation column order is coefficient aligned", {
  x <- .mt_bayesrc_fixture(); on.exit(.mt_bayesrc_cleanup(x), add = TRUE)
  A <- cbind(intercept = 1, coding = c(0, 1, 0), score = c(-1, 0, 1))
  rownames(A) <- x$stats$marker_names
  alpha <- matrix(c(-1, .2, -.1, .5, -.3, .25), 3L, 2L,
    dimnames = list(colnames(A), c("step_1", "step_2")))
  run <- function(annotation, alpha_init) {
    args <- .mt_bayesrc_common(alpha_init = alpha_init)
    args$annotations <- annotation
    do.call(mtblr_csr, c(list(stats = x$stats, ld_prefix = x$prefix,
      ld_metadata = x$ld_metadata), args))
  }
  base <- run(A, alpha)
  order <- c("intercept", "score", "coding")
  permuted <- run(A[, order, drop = FALSE], alpha[order, , drop = FALSE])
  for (field in c("bm", "dm", "component_probabilities", "vgs", "vle", "vld"))
    expect_equal(permuted[[field]], base[[field]], tolerance = 0)
  expect_equal(
    permuted$model_parameters$annotations$prior_component_probabilities,
    base$model_parameters$annotations$prior_component_probabilities,
    tolerance = 0)
})
