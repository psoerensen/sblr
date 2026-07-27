test_that("logical task seeds are independent of worker count", {
  one <- sblr:::.blr_chain_controls(10, 5, 1, 17, 3, 1, NULL, FALSE)
  many <- sblr:::.blr_chain_controls(10, 5, 1, 17, 3, 3, NULL, FALSE)
  expect_identical(sblr:::.blr_st_task_seeds(one, 2L),
                   sblr:::.blr_st_task_seeds(many, 2L))
  explicit <- sblr:::.blr_chain_controls(
    10, 5, 1, 17, 3, 2, c(41, 42, 43), FALSE)
  expect_identical(sblr:::.blr_st_task_seeds(explicit, 2L),
                   sblr:::blr_phase4_scalar_seeds_cpp(
                     17L, 2L, 3L, c(41L, 42L, 43L)))
})

.blr_unified_mt_comparable <- function(fit) {
  fields <- c(
    "bm", "dm", "b", "d", "vbs", "vgs", "ves", "vle", "vld",
    "pi_final", "pi_mean", "cov_b_mean", "cov_g_mean", "cov_e_mean",
    "cov_b_final", "cov_g_final", "cov_e_final",
    "bm_chain_mean_sd", "bm_chain_mean_min", "bm_chain_mean_max",
    "dm_chain_mean_sd", "dm_chain_mean_min", "dm_chain_mean_max",
    "convergence")
  fit[intersect(fields, names(fit))]
}

test_that("MT CSR native logical chains share one prepared operator", {
  case <- phase17j_public_case()
  on.exit(unlink(paste0(case$x$prefixes, c(
    ".row_ptr.u64.bin", ".col_idx.u32.0based.bin", ".values.f32.bin",
    ".meta.txt"))), add = TRUE)
  call <- function(ncores, keep_chains = FALSE, seeds = NULL) {
    phase17j_call(
      case, nchains = 4L, ncores = ncores, chain_seeds = seeds,
      keep_chains = keep_chains, convergence = "core",
      convergence_control = list(warn = FALSE, keep_traces = TRUE))
  }
  serial <- call(1L)
  parallel <- call(4L)
  retained <- call(2L, keep_chains = TRUE)
  expect_identical(.blr_unified_mt_comparable(parallel),
                   .blr_unified_mt_comparable(serial))
  expect_identical(.blr_unified_mt_comparable(retained),
                   .blr_unified_mt_comparable(serial))
  expect_identical(serial$diagnostics$multichain$operator_preparations, 1L)
  expect_identical(serial$diagnostics$multichain$implementation,
                   "native_static_chain_dispatch_shared_operator_preparation")
  expect_identical(serial$input$chain_seeds_resolved,
                   parallel$input$chain_seeds_resolved)
  explicit <- call(4L, seeds = c(31L, 41L, 59L, 26L))
  expect_identical(explicit$input$chain_seeds_resolved,
                   c(31, 41, 59, 26))
  expect_null(serial$chains)
  expect_length(retained$chains, 4L)
  expect_equal(serial$vld, serial$vgs - serial$vle, tolerance = 1e-12)
  expect_identical(unique(serial$convergence$summary$group),
                   c("vbs", "vgs", "ves", "vle", "vld"))
  expect_identical(dim(serial$convergence_traces$values),
                   c(8L, 4L, 5L * ncol(serial$bm)))
})

test_that("MT block-eigen native logical chains share one prepared operator", {
  case <- phase17m_public_case()
  on.exit(unlink(case$stats$bed_files), add = TRUE)
  call <- function(ncores, keep_chains = FALSE, seeds = NULL) {
    phase17m_call(
      case, nchains = 4L, ncores = ncores, chain_seeds = seeds,
      keep_chains = keep_chains, convergence = "core",
      convergence_control = list(warn = FALSE, keep_traces = TRUE))
  }
  serial <- call(1L)
  parallel <- call(4L)
  retained <- call(2L, keep_chains = TRUE)
  expect_identical(.blr_unified_mt_comparable(parallel),
                   .blr_unified_mt_comparable(serial))
  expect_identical(.blr_unified_mt_comparable(retained),
                   .blr_unified_mt_comparable(serial))
  expect_identical(serial$diagnostics$multichain$operator_preparations, 1L)
  expect_identical(serial$diagnostics$multichain$implementation,
                   "native_static_chain_dispatch_shared_operator_preparation")
  expect_identical(serial$input$chain_seeds_resolved,
                   parallel$input$chain_seeds_resolved)
  explicit <- call(4L, seeds = c(31L, 41L, 59L, 26L))
  expect_identical(explicit$input$chain_seeds_resolved,
                   c(31, 41, 59, 26))
  expect_null(serial$chains)
  expect_length(retained$chains, 4L)
  expect_equal(serial$vld, serial$vgs - serial$vle, tolerance = 1e-12)
  expect_identical(unique(serial$convergence$summary$group),
                   c("vbs", "vgs", "ves", "vle", "vld"))
  expect_identical(dim(serial$convergence_traces$values),
                   c(8L, 4L, 5L * ncol(serial$bm)))
})
