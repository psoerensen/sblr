test_that("Phase 17M public shared route equals the internal Phase 17L route", {
  case <- phase17m_public_case()
  fit <- phase17m_call(case, block_start = c(1L, 3L),
                       eigen_filter = "hard_truncate")
  expect_s3_class(fit, "mtblr_fit")
  expect_identical(fit$input$backend, "mt_block_eigen_bayesc")
  expect_identical(fit$input$summary_reference,
                   "same_bed_by_construction")
  expect_identical(fit$input$summary_ww_policy,
    "validated_by_construction_not_used_as_runtime_diagonal")
  expect_identical(fit$block_diagnostics$owner_count, 1L)
  expect_identical(fit$block_diagnostics$trait_owner, c(1L, 1L))
  expect_equal(fit$block_diagnostics$owners[[1L]]$blocks$start_1based,
               c(1L, 3L))
  expect_identical(fit$data$alignment$wy_transformation_status,
                   rep("projected_hard_truncate", 2L))

  mod <- sblr:::.mtblr_models(NULL, NULL, .001, 2L)
  vy <- case$stats$yy / (case$stats$n - 1)
  vg <- diag(vy * .5, 2L)
  ve <- diag(vy * .5, 2L)
  vb <- diag((vy * .5) / (4 * (1 - mod$probabilities[1L])), 2L)
  ssb <- ((4 - 2) / 4) * vg / (4 * (1 - mod$probabilities[1L]))
  sse <- ((4 - 2) / 4) * ve
  descriptor <- list(
    bed_files = normalizePath(case$stats$bed_files, winslash = "/"),
    n_bed = 8L, cls = list(1:4), rows = 1:8,
    af = unlist(case$stats$af), block_start = c(0L, 2L),
    eigen_filter = "hard_truncate", eigen_tau = .01, eigen_eta = .5)
  legacy <- sblr:::mtblr_block_eigen_internal(
    case$stats$wy, case$stats$yy, list(numeric(4), numeric(4)),
    list(descriptor), list(0:3), vb, ve,
    lapply(seq_len(2), function(i) ssb[, i]),
    lapply(seq_len(2), function(i) sse[, i]), mod$native,
    mod$probabilities, 4, 4, FALSE, FALSE, FALSE,
    case$stats$n, 8L, 3L, 2L, 17013L, 4L)
  expect_equal(unname(fit$bm), do.call(cbind, legacy[[1L]]), tolerance = 1e-12)
  expect_equal(unname(fit$dm), do.call(cbind, legacy[[2L]]), tolerance = 1e-12)
  expect_equal(unname(fit$wy), do.call(cbind, legacy[[3L]]), tolerance = 1e-12)
  expect_equal(unname(fit$r), do.call(cbind, legacy[[4L]]), tolerance = 1e-12)
  expect_equal(unname(fit$b), do.call(cbind, legacy[[5L]]), tolerance = 1e-12)
})

test_that("Phase 17M trait-specific filters retain consumed wy and diagnostics", {
  case <- phase17m_public_case(
    filters = c("hard_truncate", "ridge_fixed"),
    blocks = list(c(1L, 3L), c(1L, 2L, 4L)))
  fit <- phase17m_call(case, operator_sharing = "trait_specific",
                       eigen_tau = c(.9, .01))
  expect_identical(fit$block_diagnostics$owner_count, 2L)
  expect_identical(fit$block_diagnostics$trait_owner, 1:2)
  expect_identical(fit$operator_metadata$sharing_mode,
                   "trait_specific_operator")
  expect_identical(fit$data$alignment$wy_transformation_status,
                   c("projected_hard_truncate", "unchanged_ridge_fixed"))
  expect_equal(fit$wy[, 2L], case$stats$wy[[2L]], tolerance = 0)
  expect_false(isTRUE(all.equal(fit$wy[, 1L], case$stats$wy[[1L]])))
})

test_that("Phase 17M fails closed on provenance and public configuration", {
  case <- phase17m_public_case()
  external <- case$stats
  external$source <- "external"
  expect_error(phase17m_call(case, stats = external),
               "External GWAS/reference-panel")
  missing <- case$stats
  missing$bed_files <- NULL
  expect_error(phase17m_call(case, stats = missing),
               "construction provenance")
  bad_rows <- case$stats
  bad_rows$rows <- rev(bad_rows$rows)
  expect_error(phase17m_call(case, stats = bad_rows),
               "Analysis sample size|selected rows|row order|provenance")
  bad_af <- case$stats
  bad_af$af[[1L]][1L] <- .2
  expect_error(phase17m_call(case, stats = bad_af), "Allele frequencies")
  expect_error(phase17m_call(case, block_start = c(0L, 3L)),
               "begin at 1")
  expect_error(phase17m_call(case, block_start = c(1L, 1L)),
               "strictly ascending")
  expect_error(phase17m_call(case, eigen_filter = "bad"), "unsupported")
  expect_error(phase17m_call(case, sample_overlap = "modeled"),
               "not_modeled")
  expect_error(phase17m_call(case, summary_reference = "external"),
               "same_bed_by_construction")
})

test_that("Phase 17M named raw is schema version one and agrees with legacy", {
  case <- phase17m_public_case(nt = 1L, filters = "ridge_lw",
                              blocks = list(c(1L, 3L)))
  fit <- phase17m_call(case, block_start = c(1L, 3L),
                       eigen_filter = "ridge_lw")
  expect_identical(fit$raw_schema_version, 1L)
  expect_identical(fit$input$operator_sharing_mode,
                   "fully_shared_operator")
  expect_identical(fit$data$alignment$wy_transformation_status,
                   "unchanged_ridge_lw")
  expect_equal(fit$wy[, 1L], case$stats$wy[[1L]], tolerance = 0)
})

test_that("Phase 17M public blocks and sharing validation are strict", {
  case <- phase17m_public_case()
  expect_error(phase17m_call(case, block_start = c(2L, 3L)),
               "begin at 1")
  expect_error(phase17m_call(case, block_start = c(1L, 5L)),
               "in \\[1, m\\]")
  expect_error(phase17m_call(case, eigen_tau = -1), "nonnegative")
  expect_error(phase17m_call(case, eigen_eta = Inf), "finite")
  expect_error(phase17m_call(case, eigen_filter = c("ridge_fixed", "ridge_lw"),
                             operator_sharing = "shared"),
               "requires one common")
  replicated <- phase17m_call(case, operator_sharing = "trait_specific",
                              block_start = c(1L, 3L),
                              eigen_filter = "ridge_fixed")
  expect_identical(replicated$block_diagnostics$owner_count, 2L)
  expect_identical(replicated$operator_metadata$sharing_mode,
                   "trait_specific_shared_boundaries")
})

test_that("Phase 17M direct named raw agrees with the legacy adapter", {
  case <- phase17l_case(filters = c("hard_truncate", "ridge_fixed"),
                        shared = FALSE,
                        blocks = list(c(0L, 2L), c(0L, 1L, 3L)))
  raw <- do.call(sblr:::mtblr_block_eigen_raw_internal, case$block)
  legacy <- do.call(sblr:::mtblr_block_eigen_internal, case$block)
  expect_true(sblr:::.is_mtblr_raw(raw))
  expect_identical(raw$meta$backend, "mt_block_eigen_bayesc")
  expect_identical(raw$diagnostics$block_eigen$owner_count, 2L)
  expect_identical(raw$diagnostics$block_eigen$trait_owner, 1:2)
  expect_equal(raw$marker$bm, do.call(cbind, legacy[[1L]]), tolerance = 1e-12)
  expect_equal(raw$marker$dm, do.call(cbind, legacy[[2L]]), tolerance = 1e-12)
  expect_equal(raw$marker$wy, do.call(cbind, legacy[[3L]]), tolerance = 1e-12)
  expect_equal(raw$marker$r, do.call(cbind, legacy[[4L]]), tolerance = 1e-12)
  expect_equal(raw$marker$b, do.call(cbind, legacy[[5L]]), tolerance = 1e-12)
  expect_identical(raw$marker$state,
                   matrix(as.integer(unlist(legacy[[6L]])),
                          ncol = 2L))
  expect_silent(sblr:::.validate_mtblr_raw(raw))
})

test_that("Phase 17M reorder_stats preserves BED provenance", {
  case <- phase17m_public_case()
  permutation <- c(3L, 1L, 4L, 2L)
  reordered <- case$stats
  reordered$wy <- lapply(reordered$wy, `[`, permutation)
  reordered$ww <- lapply(reordered$ww, `[`, permutation)
  reordered$marker_names <- reordered$marker_names[permutation]
  reordered$marker_metadata <-
    reordered$marker_metadata[permutation, , drop = FALSE]
  reordered$cls[[1L]] <- reordered$cls[[1L]][permutation]
  reordered$af[[1L]] <- reordered$af[[1L]][permutation]
  expect_error(phase17m_call(case, stats = reordered),
               "strict marker policy")
  fit <- phase17m_call(case, stats = reordered,
                       marker_policy = "reorder_stats")
  expect_true(all(fit$alignment$per_trait$marker_order_status ==
                  "stats_reordered"))
  expect_identical(rownames(fit$wy), case$stats$marker_names)
})

test_that("Phase 17M accepts named single-trait stats and trait-specific n", {
  first <- phase17m_public_case(nt = 1L, rows = 1:8)
  second <- phase17m_public_case(nt = 1L, rows = 1:6)
  to_entry <- function(case, name) {
    x <- case$stats
    list(
      wy = x$wy[[1L]], ww = x$ww[[1L]], yy = x$yy[[1L]], n = x$n[[1L]],
      n_bed = x$n_bed, bed_files = x$bed_files, cls = x$cls, af = x$af,
      rows = x$rows, marker_names = x$marker_names,
      marker_metadata = x$marker_metadata, scale = x$scale,
      source = x$source, trait_names = name)
  }
  stats <- list(T1 = to_entry(first, "T1"), T2 = to_entry(second, "T2"))
  first$Glist$idsLD <- first$Glist$ids[1:8]
  second$Glist$idsLD <- second$Glist$ids[1:6]
  fit <- mtblr_block_eigen(
    stats, list(first$Glist, second$Glist),
    block_start = list(c(1L, 3L), c(1L, 2L, 4L)),
    operator_sharing = "trait_specific",
    eigen_filter = c("ridge_fixed", "ridge_lw"),
    eigen_eta = c(.5, 0), updateB = FALSE, updateE = FALSE,
    updatePi = FALSE, nit = 8L, nburn = 3L, nthin = 2L, seed = 17014L)
  expect_identical(fit$input$n, c(8L, 6L))
  expect_identical(fit$operator_metadata$selected_row_count, c(8L, 6L))
})
