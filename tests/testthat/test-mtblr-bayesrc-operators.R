test_that("MT SBayesRC CSR and exact block eigen reduce numerically", {
  x <- .mt_bayesrc_fixture(); on.exit(.mt_bayesrc_cleanup(x), add = TRUE)
  common <- .mt_bayesrc_common()
  common$annotations <- x$annotations
  common$updateAlpha <- TRUE
  common$alpha_update_every <- 2L
  csr <- do.call(mtblr_csr, c(list(stats = x$stats, ld_prefix = x$prefix,
    ld_metadata = x$ld_metadata), common))
  block <- do.call(mtblr_block_eigen, c(list(stats = x$fixture$stats,
    Glist = x$fixture$Glist, block_start = 1L, eigen_filter = "hard_truncate",
    eigen_tau = 0), common))
  for (field in c("bm", "dm", "component_probabilities", "vbs", "vgs",
                  "ves", "vle", "vld"))
    expect_equal(block[[field]], csr[[field]], tolerance = 1e-7, info = field)
  for (field in c("annotation_coefficients_final",
                  "annotation_coefficients_mean", "annotation_variances_final",
                  "annotation_variances_mean", "pattern_pi_final",
                  "pattern_pi_mean", "prior_component_probabilities"))
    expect_equal(block$model_parameters$annotations[[field]],
                 csr$model_parameters$annotations[[field]], tolerance = 1e-7,
                 info = field)
})

test_that("MT BayesRC packed BED supports both residual policies", {
  x <- .mt_bayesrc_fixture(); on.exit(.mt_bayesrc_cleanup(x), add = TRUE)
  common <- .mt_bayesrc_common("bayesrc")
  common$annotations <- x$annotations
  common$updateE <- TRUE
  for (policy in c("diagonal", "full")) {
    fit <- do.call(mtblr_bed, c(list(y = x$fixture$y,
      Glist = x$fixture$Glist, residual_covariance = policy), common))
    expect_identical(fit$model, "bayesrc")
    expect_equal(unname(rowSums(fit$component_probabilities)), rep(1, 3),
                 tolerance = 1e-12)
    expect_equal(fit$vld, fit$vgs - fit$vle, tolerance = 1e-12)
    expect_null(fit$pi_final)
  }
})

test_that("MT BayesRC multichain results are worker and retention independent", {
  x <- .mt_bayesrc_fixture(); on.exit(.mt_bayesrc_cleanup(x), add = TRUE)
  common <- .mt_bayesrc_common()
  common$annotations <- x$annotations
  common$updateAlpha <- TRUE
  common$nchains <- 2L; common$chain_seeds <- c(101L, 202L)
  common$convergence <- "core"
  common$convergence_control <- list(warn = FALSE, keep_traces = TRUE)
  run <- function(cores, keep) {
    args <- common; args$ncores <- cores; args$keep_chains <- keep
    do.call(mtblr_csr, c(list(stats = x$stats, ld_prefix = x$prefix,
      ld_metadata = x$ld_metadata), args))
  }
  serial <- run(1L, FALSE); repeated <- run(1L, FALSE)
  retained <- run(1L, TRUE); parallel <- run(2L, FALSE)
  for (field in c("bm", "dm", "component_probabilities", "vbs", "vgs",
                  "ves", "vle", "vld", "convergence", "convergence_traces")) {
    expect_equal(repeated[[field]], serial[[field]], tolerance = 0)
    expect_equal(retained[[field]], serial[[field]], tolerance = 0)
    expect_equal(parallel[[field]], serial[[field]], tolerance = 0)
  }
  expect_null(serial$chains)
  expect_length(retained$chains, 2L)
})
