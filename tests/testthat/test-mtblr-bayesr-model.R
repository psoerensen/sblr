test_that("MT BayesR state specification is unique, ordered, and validated", {
  patterns <- sblr:::.mtblr_models(
    matrix(c(0L, 0L, 1L, 0L, 1L, 1L), ncol = 2L, byrow = TRUE),
    c(.7, .2, .1), .1, 2L)
  spec <- sblr:::.mtblr_bayesr_spec(
    "bayesr", patterns, NULL, 3L, c(0, .01, .1),
    component = c(0L, 1L, 2L))
  expect_identical(spec$joint_names, c(
    "null", "1_0__component_1", "1_0__component_2",
    "1_1__component_1", "1_1__component_2"))
  expect_identical(spec$joint_component, c(0L, 1L, 2L, 1L, 2L))
  expect_equal(spec$joint_multiplier, c(0, .01, .1, .01, .1))
  expect_equal(sum(spec$patterns$probabilities), 1)
  expect_error(sblr:::.mtblr_bayesr_spec(
    "bayesr", patterns, NULL, 3L, c(0, .1, .1)), "unique ascending")
  explicit_unit <- sblr:::.mtblr_bayesr_spec(
    "bayesr", patterns, rep(.25, 3L), 3L, c(0, 1), maf_effect_s = -1)
  expect_equal(explicit_unit$marker_scale, rep(1, 3L))
})

test_that("MT fixed maf_effect_s independently applies the canonical MAF-S scale", {
  patterns <- sblr:::.mtblr_models(matrix(c(0L, 1L), 2L, 1L),
                                    c(.8, .2), .2, 1L)
  frequency <- c(.1, .25, .4)
  scaled <- sblr:::.mtblr_bayesr_spec(
    "bayesr", patterns, frequency, 3L, c(0, 1), maf_effect_s = 0)
  expect_equal(scaled$marker_scale, 2 * frequency * (1 - frequency))
  unit <- sblr:::.mtblr_bayesr_spec(
    "bayesr", patterns, frequency, 3L, c(0, 1), maf_effect_s = -1)
  expect_equal(unit$marker_scale, rep(1, 3))
  expect_error(sblr:::.mtblr_bayesr_spec(
    "bayesr", patterns, frequency, 3L, c(0, 1),
    estimate_maf_effect_s = TRUE), "not implemented")
})

test_that("MT BayesR initialization enforces pattern-component consistency", {
  patterns <- matrix(c(0L, 0L, 1L, 0L, 1L, 1L), ncol = 2L,
                     byrow = TRUE)
  state <- matrix(c(0L, 0L, 1L, 0L), 2L, 2L, byrow = TRUE)
  beta <- matrix(c(0, 0, .2, 0), 2L, 2L, byrow = TRUE)
  good <- sblr:::.mtblr_bayesr_initialization(
    beta, beta, state, c(0L, 1L), patterns, 2L, 2L, "bayesr")
  expect_identical(good$component, c(0L, 1L))
  expect_error(sblr:::.mtblr_bayesr_initialization(
    beta, beta, state, c(1L, 1L), patterns, 2L, 2L, "bayesr"),
    "zero exactly")
})

test_that("shared MT BayesR state kernel matches an independent scalar oracle", {
  score <- 1.25
  diagonal <- 8
  B <- matrix(.2)
  E <- matrix(.5)
  gamma <- c(0, .1, 1)
  pi <- c(.7, .2, .1)
  native <- sblr:::mtblr_bayesr_marker_contract_internal(
    score, diagonal, B, E, list(0L, 1L, 1L), c(0L, 1L, 2L),
    gamma, c("null", "1__component_1", "1__component_2"),
    3L, pi, 1)
  prior_precision <- 1 / (B[1, 1] * gamma[-1L])
  precision <- prior_precision + diagonal / E[1, 1]
  rhs <- score / E[1, 1]
  weight <- c(pi[1L], pi[-1L] * sqrt(prior_precision / precision) *
                exp(.5 * rhs^2 / precision))
  probability <- weight / sum(weight)
  expect_equal(native$probability, probability, tolerance = 1e-14)
  expect_equal(vapply(native$mean[-1L], as.numeric, numeric(1)),
               rhs / precision, tolerance = 1e-14)
  expect_equal(vapply(native$covariance[-1L], as.numeric, numeric(1)),
               1 / precision, tolerance = 1e-14)
})

test_that("one unit positive component reduces to the BayesC marker kernel", {
  score <- 1.25; diagonal <- 8; B <- matrix(.2); E <- matrix(.5)
  probability <- c(.7, .3)
  bayesc <- sblr:::mtblr_bed_marker_contract_internal(
    score, diagonal, B, E, list(0L, 1L), probability)
  bayesr <- sblr:::mtblr_bayesr_marker_contract_internal(
    score, diagonal, B, E, list(0L, 1L), c(0L, 1L), c(0, 1),
    c("null", "1__component_1"), 2L, probability, 1)
  expect_equal(bayesr$probability, bayesc$probability, tolerance = 1e-14)
  expect_equal(as.numeric(bayesr$mean[[2L]]), bayesc$mean[[2L]],
               tolerance = 1e-14)
  expect_equal(bayesr$covariance[[2L]], bayesc$covariance[[2L]],
               tolerance = 1e-14)
})

test_that("MT BayesR public controls are aligned and failure-closed", {
  bayesr_controls <- c("mixture_var", "joint_pi", "joint_pi_prior", "component")
  bayesrc_controls <- c("annotations", "add_intercept",
    "standardize_annotations", "center_binary_annotations", "alpha_init",
    "sigmaSqAlpha_init", "intercept_flat", "sigmaSqAlpha_a",
    "sigmaSqAlpha_b", "pi_floor", "alpha_update_every", "updateAlpha")
  selection_controls <- c(
    "maf_effect_s", "effect_maf", "allow_reference_maf_for_maf_effect_s",
    "estimate_maf_effect_s", "maf_effect_s_init",
    "maf_effect_s_prior", "maf_effect_s_proposal_sd")
  controls <- c(bayesr_controls, bayesrc_controls, selection_controls)
  for (fun in list(mtblr_csr, mtblr_block_eigen, mtblr_bed)) {
    form <- names(formals(fun))
    start <- match("mixture_var", form)
    expect_identical(form[seq.int(start, length.out = length(controls))],
                     controls)
  }
  patterns <- sblr:::.mtblr_models(matrix(c(0L, 1L), 2L, 1L),
                                    c(.8, .2), .2, 1L)
  expect_error(sblr:::.mtblr_bayesr_spec(
    "bayesr", patterns, NULL, 2L, c(0, 1), joint_pi = c(.5, .3, .2)),
    "joint-state count")
  expect_error(sblr:::.mtblr_bayesr_spec(
    "bayesr", patterns, c(.2, NA), 2L, c(0, 1), maf_effect_s = 0),
    "allele frequencies")
  many <- list(matrix = rbind(matrix(0L, 1L, 13L),
    as.matrix(expand.grid(rep(list(0:1), 13L)))[-1L, , drop = FALSE]),
    probabilities = rep(1 / 8192, 8192),
    names = paste0("p", seq_len(8192)))
  expect_error(sblr:::.mtblr_bayesr_spec(
    "bayesr", many, NULL, 1L, c(0, 1)), "exceeds 4096")
})
