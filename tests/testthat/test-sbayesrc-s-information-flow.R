test_that("Phase-4D conditional probabilities and stick algebra are exact", {
  probability <- rbind(
    c(0.40, 0.30, 0.20, 0.10),
    c(0.10, 0.20, 0.30, 0.40),
    c(0.25, 0.25, 0.25, 0.25)
  )
  prior <- rbind(
    c(0.25, 0.25, 0.25, 0.25),
    c(0.20, 0.30, 0.30, 0.20),
    c(0.40, 0.20, 0.20, 0.20)
  )
  annotation <- cbind(Intercept = 1, a = 1:3)
  out <- .st_bayesrc_information_summary(
    log(probability), prior, c(0L, 2L, 3L), annotation
  )

  expect_equal(out$rb_comp_prob, probability, tolerance = 1e-15)
  expect_equal(rowSums(out$rb_comp_prob), rep(1, 3L), tolerance = 1e-15)
  expect_true(all(out$rb_comp_prob >= 0))
  expect_equal(out$rb_dm, 1 - probability[, 1L], tolerance = 1e-15)
  expect_equal(
    as.numeric(out$hard_stick_trace),
    c(3, 2, 2 / 3, 2, 2, 1, 2, 1, 1 / 2),
    tolerance = 1e-15
  )
  expect_equal(
    as.numeric(out$soft_stick_trace),
    c(3, 2.25, 0.75, 2.25, 1.5, 2 / 3, 1.5, 0.75, 0.5),
    tolerance = 1e-15
  )
  hard <- matrix(out$hard_stick_trace, nrow = 1L)
  soft <- matrix(out$soft_stick_trace, nrow = 1L)
  expect_true(all(hard[, c(2L, 5L, 8L)] <= hard[, c(1L, 4L, 7L)]))
  expect_true(all(soft[, c(2L, 5L, 8L)] <= soft[, c(1L, 4L, 7L)]))
  expect_true(all(diff(soft[1L, c(1L, 4L, 7L)]) <= 0))

  hard_e1 <- c(1, 1, 1)
  hard_s1 <- c(0, 1, 1)
  expect_equal(as.numeric(out$hard_annotation_information[1L, 1:2]),
               as.numeric(colSums(annotation * hard_e1)), tolerance = 1e-15)
  expect_equal(as.numeric(out$hard_annotation_information[1L, 3:4]),
               as.numeric(colSums(annotation * hard_s1)), tolerance = 1e-15)
  expect_equal(as.numeric(out$hard_annotation_information[1L, 5:6]),
               as.numeric(colSums(annotation^2 * hard_e1)), tolerance = 1e-15)
  expect_true(all(is.finite(out$information_gain)))
  expect_true(all(out$information_gain >= 0))
})

test_that("Phase-4D diagnostics preserve the complete scientific RNG state", {
  fixture <- .sbs4b_fixture(36L, 20271150L)
  off <- .sbs4b_run(
    fixture, 20271151L, 100L, 20L,
    updateB = TRUE, updateE = FALSE,
    information_diagnostics = FALSE
  )
  on <- .sbs4b_run(
    fixture, 20271151L, 100L, 20L,
    updateB = TRUE, updateE = FALSE,
    information_diagnostics = TRUE
  )
  off_chain <- off$chains[[1L]][[1L]]
  on_chain <- on$chains[[1L]][[1L]]

  for (field in c("marker", "trace", "pi", "component", "annotation",
                  "convergence_trace", "selection")) {
    expect_identical(on_chain[[field]], off_chain[[field]], info = field)
  }
  for (field in c("marker", "trace", "variance", "pi", "annotation",
                  "component")) {
    expect_identical(on[[field]], off[[field]], info = field)
  }
  expect_false("information_flow" %in% names(off_chain))
  expect_true("information_flow" %in% names(on_chain))

  flow <- on_chain$information_flow
  expect_equal(rowSums(flow$rb_comp_prob), rep(1, 36L), tolerance = 1e-12)
  expect_true(all(is.finite(flow$rb_comp_prob)))
  expect_true(all(flow$rb_comp_prob >= 0 & flow$rb_comp_prob <= 1))
  expect_equal(flow$rb_dm, 1 - flow$rb_comp_prob[, 1L], tolerance = 1e-15)
  expect_equal(rowSums(on_chain$component$prob), rep(1, 36L),
               tolerance = 1e-12)
  expect_equal(on_chain$marker$dm,
               rowSums(on_chain$component$prob[, -1L, drop = FALSE]),
               tolerance = 1e-15)

  hard_stick <- flow$hard_stick_trace
  soft_stick <- flow$soft_stick_trace
  expect_true(all(hard_stick[, c(2L, 5L, 8L)] <=
                    hard_stick[, c(1L, 4L, 7L)]))
  expect_true(all(soft_stick[, c(2L, 5L, 8L)] <=
                    soft_stick[, c(1L, 4L, 7L)] + 1e-12))
  expect_true(all(soft_stick[, 1L] >= soft_stick[, 4L] - 1e-12))
  expect_true(all(soft_stick[, 4L] >= soft_stick[, 7L] - 1e-12))
  expect_true(all(is.finite(flow$information_gain)))
})

test_that("Phase-4D internal standard hierarchy diagnostics are coherent", {
  fixture <- .sbs4b_fixture(32L, 20271160L)
  fit <- .sbs4b_run(
    fixture, 20271161L, 80L, 20L,
    fixed_delta = rep(1L, 3L), update_pi_A = FALSE, update_tau2 = FALSE,
    selection_enabled = FALSE, information_diagnostics = TRUE
  )
  chain <- fit$chains[[1L]][[1L]]
  expect_identical(fit$meta$model, "sbayesrc")
  expect_true("information_flow" %in% names(chain))
  expect_equal(rowSums(chain$information_flow$rb_comp_prob), rep(1, 32L),
               tolerance = 1e-12)
  expect_null(chain$annotation$annotation_pip)
})
