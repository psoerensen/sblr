test_that("MT BayesRC factorized annotation probabilities have one null state", {
  A <- cbind(intercept = 1, x = c(-1, 0, 1))
  alpha <- matrix(c(.2, -.1, .4, .3), 2L, 2L)
  theta <- sblr:::.mtblr_bayesrc_prior_probabilities(A, alpha)
  expect_equal(rowSums(theta), rep(1, 3), tolerance = 1e-14)
  omega <- c(.25, .75)
  joint <- cbind(theta[, 1L], theta[, 2L] * omega[1L],
                 theta[, 3L] * omega[1L], theta[, 2L] * omega[2L],
                 theta[, 3L] * omega[2L])
  expect_equal(rowSums(joint), rep(1, 3), tolerance = 1e-14)
  expect_equal(1 - theta[, 1L], rowSums(theta[, -1L, drop = FALSE]),
               tolerance = 1e-14)
})

test_that("fixed-alpha MT SBayesRC reduces exactly to MT SBayesR", {
  x <- .mt_bayesrc_fixture(); on.exit(.mt_bayesrc_cleanup(x), add = TRUE)
  init <- make_sbayesrc_alpha_init(
    x$annotations, gamma = c(0, .1, 1), pi_init = .3,
    active_comp_weights = c(.5, .5))
  common <- .mt_bayesrc_common(alpha_init = init$alpha_init)
  common$annotations <- x$annotations
  rc <- do.call(mtblr_csr, c(list(stats = x$stats, ld_prefix = x$prefix,
    ld_metadata = x$ld_metadata), common))
  r <- mtblr_csr(x$stats, ld_prefix = x$prefix,
    ld_metadata = x$ld_metadata, method = "sbayesr",
    mixture_var = c(0, .1, 1), models = matrix(c(0L, 1L), 2L, 1L),
    joint_pi = c(.7, .15, .15), vb = matrix(.1), ve = matrix(.5),
    updateB = FALSE, updateE = FALSE, updatePi = FALSE,
    nit = 8L, nburn = 2L, seed = 42L, convergence = "none")
  for (field in c("bm", "dm", "component_final",
                  "component_probabilities", "vgs", "vle", "vld"))
    expect_equal(rc[[field]], r[[field]], tolerance = 0, info = field)
  expect_null(rc$pi_final)
  expect_equal(rowSums(rc$model_parameters$annotations$
                         prior_component_probabilities), rep(1, 3),
               tolerance = 1e-14)
})

test_that("BayesRC update and selection-S controls are independent", {
  x <- .mt_bayesrc_fixture(); on.exit(.mt_bayesrc_cleanup(x), add = TRUE)
  call <- function(selection) {
    args <- .mt_bayesrc_common(maf_effect_s = selection)
    args$annotations <- x$annotations
    do.call(mtblr_csr, c(list(stats = x$stats, ld_prefix = x$prefix,
      ld_metadata = x$ld_metadata), args))
  }
  plain <- call(NULL); unit <- call(-1)
  for (field in c("bm", "dm", "component_probabilities", "vgs", "vle", "vld"))
    expect_equal(unit[[field]], plain[[field]], tolerance = 0, info = field)
  expect_identical(plain$input$effect_scale_policy, "component")
  expect_identical(unit$input$effect_scale_policy, "component_maf_s")
  expect_error(mtblr_csr(x$stats, ld_prefix = x$prefix,
    ld_metadata = x$ld_metadata, method = "sbayesrc",
    annotations = x$annotations, estimate_maf_effect_s = TRUE),
    "Sampled maf_effect_s")
})
