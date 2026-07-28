test_that("selected marker diagnostics preserve request order and meanings", {
  fixture <- blr_unified_fixture()
  on.exit(blr_unified_cleanup(fixture), add = TRUE)
  common <- list(
    y = fixture$y, Glist = fixture$Glist, method = "bayesr",
    mixture_var = c(0, .1, 1), models = matrix(c(0L, 1L), 2L, 1L),
    joint_pi = c(.7, .15, .15), joint_pi_prior = rep(1, 3),
    vb = matrix(.1), ve = matrix(.5), updateB = FALSE, updateE = FALSE,
    nit = 8L, nburn = 2L, nchains = 2L, ncores = 1L,
    residual_covariance = "diagonal")
  fit <- do.call(mtblr_bed, c(common, list(
    convergence = "extended", convergence_control = list(
      warn = FALSE, selected_markers = 2:1,
      selected_marker_quantities = c("b", "d", "component"),
      keep_traces = TRUE))))
  selected <- fit$convergence$summary[fit$convergence$summary$tier == 3L, ]
  expect_identical(selected$quantity,
    c("b[m2,T1]", "b[m1,T1]", "d[m2,T1]", "d[m1,T1]",
      "component[m2]", "component[m1]"))
  d_values <- fit$convergence_traces$values[, , selected$quantity[3:4], drop = FALSE]
  expect_true(all(d_values %in% c(0, 1)))
  expect_identical(unique(selected$marker_index), c(2L, 1L))
})

test_that("diagnostic capture is RNG neutral", {
  fixture <- blr_unified_fixture()
  on.exit(blr_unified_cleanup(fixture), add = TRUE)
  common <- list(
    y = fixture$y, Glist = fixture$Glist, method = "bayesr",
    mixture_var = c(0, .1, 1), models = matrix(c(0L, 1L), 2L, 1L),
    joint_pi = c(.7, .15, .15), joint_pi_prior = rep(1, 3),
    vb = matrix(.1), ve = matrix(.5), updateB = FALSE, updateE = FALSE,
    nit = 8L, nburn = 2L, nchains = 2L, ncores = 1L,
    residual_covariance = "diagonal")
  none <- do.call(mtblr_bed, c(common, list(convergence = "none")))
  extended <- do.call(mtblr_bed, c(common, list(
    convergence = "extended", convergence_control = list(
      warn = FALSE, selected_markers = 1L,
      selected_marker_quantities = c("b", "d", "component")))))
  for (field in c("bm", "dm", "b_final", "d_final", "component_final",
                  "component_probabilities", "vbs", "vgs", "ves", "vle",
                  "vld", "pi_final", "pi_mean"))
    expect_equal(extended[[field]], none[[field]], tolerance = 0, info = field)
})

test_that("selected marker validation rejects unsafe requests", {
  expect_error(.blr_convergence_controls(
    "extended", list(selected_markers = c(1L, 1L)), 2L), "unique")
  expect_error(.blr_convergence_controls(
    "extended", list(selected_markers = TRUE), 2L), "shortcuts")
  expect_error(.blr_resolve_selected_markers(c("m1", "missing"), c("m1", "m2")),
               "Unknown selected marker")
})
