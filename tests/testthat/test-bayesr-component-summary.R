summarise_components <- getFromNamespace("summarise_components", "sblr")

make_bayesr_component_summary_fit <- function(optional = TRUE) {
  markers <- paste0("m", 1:4)
  traits <- c("D1", "D2")
  cp1 <- matrix(
    c(
      0.90, 0.05, 0.05,
      0.40, 0.20, 0.40,
      0.10, 0.30, 0.60,
      0.70, 0.20, 0.10
    ),
    nrow = 4,
    byrow = TRUE,
    dimnames = list(markers, paste0("component_", 0:2))
  )
  cp2 <- matrix(
    c(
      0.80, 0.10, 0.10,
      0.60, 0.30, 0.10,
      0.20, 0.40, 0.40,
      0.95, 0.03, 0.02
    ),
    nrow = 4,
    byrow = TRUE,
    dimnames = list(markers, paste0("component_", 0:2))
  )
  dm <- cbind(D1 = 1 - cp1[, "component_0"], D2 = 1 - cp2[, "component_0"])
  rownames(dm) <- markers
  bm <- matrix(
    c(0.01, -0.02, 0.03, 0.04, -0.01, 0.02, -0.03, 0.05),
    nrow = 4,
    dimnames = list(markers, traits)
  )
  fit <- list(
    dm = dm,
    bm = bm,
    component_probabilities = list(D1 = cp1, D2 = cp2)
  )
  if (optional) {
    fit$dm_chain_mean_sd <- matrix(
      c(0.01, 0.02, 0.03, 0.04, 0.04, 0.03, 0.02, 0.01),
      nrow = 4,
      dimnames = dimnames(dm)
    )
    fit$bm_chain_mean_sd <- matrix(
      c(0.001, 0.002, 0.003, 0.004, 0.004, 0.003, 0.002, 0.001),
      nrow = 4,
      dimnames = dimnames(dm)
    )
    fit$dm_component_mean <- matrix(
      c(0.15, 1.00, 1.50, 0.40, 0.30, 0.50, 1.20, 0.07),
      nrow = 4,
      dimnames = dimnames(dm)
    )
  }
  class(fit) <- c("stblr_fit", "blr_fit", "list")
  fit
}

test_that("summarise_components returns basic PIP and component summaries", {
  fit <- make_bayesr_component_summary_fit(optional = FALSE)
  out <- summarise_components(
    fit,
    pip_thresholds = c(0.05, 0.5, 0.95)
  )

  expect_s3_class(out, "stblr_bayesr_component_summary")
  expect_equal(nrow(out), 2L)
  expect_equal(out$trait, c("D1", "D2"))
  expect_equal(out$n_markers, c(4L, 4L))
  expect_equal(out$mean_pip[out$trait == "D1"], mean(fit$dm[, "D1"]))
  expect_equal(out$sum_pip[out$trait == "D1"], sum(fit$dm[, "D1"]))
  expect_equal(out$max_pip[out$trait == "D1"], max(fit$dm[, "D1"]))
  expect_equal(out$n_pip_gt_0_05[out$trait == "D1"], sum(fit$dm[, "D1"] > 0.05))
  expect_equal(out$n_pip_gt_0_5[out$trait == "D1"], sum(fit$dm[, "D1"] > 0.5))
  expect_equal(out$n_pip_gt_0_95[out$trait == "D1"], sum(fit$dm[, "D1"] > 0.95))
  expect_equal(out$mean_component_0[out$trait == "D1"], mean(fit$component_probabilities$D1[, "component_0"]))
  expect_equal(out$mean_component_1[out$trait == "D1"], mean(fit$component_probabilities$D1[, "component_1"]))
  expect_equal(out$max_component_2[out$trait == "D1"], max(fit$component_probabilities$D1[, "component_2"]))
})

test_that("summarise_components includes chain stability fields", {
  fit <- make_bayesr_component_summary_fit(optional = TRUE)
  out <- summarise_components(fit)

  expect_equal(out$mean_pip_sd[out$trait == "D1"], mean(fit$dm_chain_mean_sd[, "D1"]))
  expect_equal(out$max_pip_sd[out$trait == "D1"], max(fit$dm_chain_mean_sd[, "D1"]))
  expect_equal(out$mean_effect_sd[out$trait == "D1"], mean(fit$bm_chain_mean_sd[, "D1"]))
  expect_equal(out$max_effect_sd[out$trait == "D1"], max(fit$bm_chain_mean_sd[, "D1"]))
  expect_equal(out$mean_abs_effect[out$trait == "D1"], mean(abs(fit$bm[, "D1"])))
  expect_equal(out$max_abs_effect[out$trait == "D1"], max(abs(fit$bm[, "D1"])))
  expect_equal(out$mean_component_index[out$trait == "D1"], mean(fit$dm_component_mean[, "D1"]))
  expect_equal(out$max_component_index[out$trait == "D1"], max(fit$dm_component_mean[, "D1"]))
})

test_that("summarise_components reports component convention differences", {
  fit <- make_bayesr_component_summary_fit(optional = FALSE)
  out <- summarise_components(fit)
  expect_equal(out$max_dm_component0_diff, c(0, 0), tolerance = 1e-12)

  fit$dm[1, "D1"] <- fit$dm[1, "D1"] + 0.05
  out_bad <- summarise_components(fit)
  expect_equal(out_bad$max_dm_component0_diff[out_bad$trait == "D1"], 0.05, tolerance = 1e-12)
})

test_that("summarise_components handles missing optional fields", {
  fit <- make_bayesr_component_summary_fit(optional = FALSE)
  fit$bm <- NULL
  out <- summarise_components(fit)

  expect_equal(nrow(out), 2L)
  expect_false("mean_pip_sd" %in% names(out))
  expect_false("mean_effect_sd" %in% names(out))
  expect_false("mean_component_index" %in% names(out))
  expect_false("mean_abs_effect" %in% names(out))
})

test_that("summarise_components validates required fields", {
  fit <- make_bayesr_component_summary_fit(optional = FALSE)

  fit_no_dm <- fit
  fit_no_dm$dm <- NULL
  expect_error(
    summarise_components(fit_no_dm),
    "pips must be present"
  )

  fit_no_comp <- fit
  fit_no_comp$component_probabilities <- NULL
  expect_error(
    summarise_components(fit_no_comp),
    "component probabilities must be present"
  )
})

test_that("summarise_components uses sensible missing trait names", {
  fit <- make_bayesr_component_summary_fit(optional = FALSE)
  colnames(fit$dm) <- NULL
  names(fit$component_probabilities) <- c("trait1", "trait2")

  out <- summarise_components(fit)

  expect_equal(out$trait, c("trait1", "trait2"))
})

test_that("summarise_components can return top unstable markers", {
  fit <- make_bayesr_component_summary_fit(optional = TRUE)
  out <- summarise_components(fit, top_unstable = 2L)

  expect_type(out, "list")
  expect_named(out, c("summary", "unstable"))
  expect_equal(nrow(out$unstable$D1), 2L)
  expect_equal(out$unstable$D1$marker[1], "m4")
})
