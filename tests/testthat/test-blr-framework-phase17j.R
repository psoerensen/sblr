test_that("public MT CSR fit uses named raw and matches the internal route", {
  case <- phase17j_public_case()
  fit <- phase17j_call(case)
  legacy <- do.call(sblr:::mtblr_csr_internal, case$x$csr)
  expect_s3_class(fit, "mtblr_fit")
  expect_equal(unname(fit$bm), t(do.call(rbind, legacy[[1L]])), tolerance = 1e-12)
  expect_equal(unname(fit$dm), t(do.call(rbind, legacy[[2L]])), tolerance = 1e-12)
  expect_equal(unname(fit$b), t(do.call(rbind, legacy[[5L]])), tolerance = 1e-12)
  expect_equal(unname(fit$pi), legacy[[17L]][[1L]], tolerance = 1e-12)
  expect_equal(unname(fit$pim), legacy[[18L]][[1L]], tolerance = 1e-12)
  expect_identical(fit$input$backend, "mt_csr_bayesc")
  expect_identical(fit$input$sample_overlap, "not_modeled")
  expect_identical(fit$input$residual_covariance_policy, "diagonal")
  expect_false(any(grepl("row_ptr|column_index|offdiag", names(fit))))
})

test_that("public sharing modes resolve without changing the internal core", {
  shared <- phase17j_call(phase17j_public_case())
  expect_identical(shared$input$ld_sharing_mode, "fully_shared_operator")
  scaled_case <- phase17j_public_case(trait_specific = TRUE)
  scaled_case$x$csr$ld_prefixes <- scaled_case$x$csr$ld_prefixes[1L]
  scaled_case$metadata <- scaled_case$metadata[[1L]]
  scaled <- phase17j_call(scaled_case)
  expect_identical(scaled$input$ld_sharing_mode, "shared_correlation_reference")
  expect_length(scaled$input$ld_prefix, length(scaled_case$x$csr$wy))
  specific <- phase17j_call(phase17j_public_case(trait_specific = TRUE, independent = TRUE))
  expect_identical(specific$input$ld_sharing_mode, "trait_specific_reference")
})

test_that("marker policies and biological metadata fail closed", {
  case <- phase17j_public_case()
  bad <- case; bad$stats$marker_names[1L] <- bad$stats$marker_names[2L]
  expect_error(phase17j_call(bad), "unique")
  bad <- case; bad$stats$marker_names <- rev(bad$stats$marker_names)
  bad$stats$marker_metadata <- bad$stats$marker_metadata[rev(seq_len(nrow(bad$stats$marker_metadata))), ]
  bad$stats$wy <- lapply(bad$stats$wy, rev); bad$stats$ww <- lapply(bad$stats$ww, rev)
  expect_error(phase17j_call(bad), "strict marker")
  expect_s3_class(phase17j_call(bad, marker_policy = "reorder_stats"), "mtblr_fit")
  bad <- case; bad$metadata$marker_ids[1L] <- "extra"
  expect_error(phase17j_call(bad), "LD resources|marker sets|marker_metadata")
  bad <- case; bad$metadata$marker_metadata$effect_allele[1L] <- "C";
  bad$metadata$marker_metadata$other_allele[1L] <- "A"
  expect_error(phase17j_call(bad), "Swapped")
  bad <- case; bad$metadata$marker_metadata$effect_allele[1L] <- "T";
  bad$metadata$marker_metadata$other_allele[1L] <- "G"
  expect_error(phase17j_call(bad), "Strand")
  bad <- case; bad$metadata$marker_metadata <- data.frame(marker_id = bad$metadata$marker_ids)
  expect_error(phase17j_call(bad), "explicit effect_allele")
  bad <- case; bad$stats$scale <- "unknown"
  expect_error(phase17j_call(bad), "standardized_genotype")
})

test_that("overlap, residual covariance, patterns, and sets are validated", {
  case <- phase17j_public_case()
  expect_error(phase17j_call(case, sample_overlap = "known"), "not_modeled")
  bad <- case; bad$stats$yy <- matrix(c(50, 1, 1, 50), 2)
  expect_error(phase17j_call(bad), "off-diagonal yy")
  expect_error(phase17j_call(case, ve = matrix(c(.8, .1, .1, .8), 2)), "diagonal")
  expect_error(phase17j_call(case, sse_prior = matrix(c(.3, .1, .1, .3), 2)), "diagonal")
  expect_error(phase17j_call(case, sets = list(c(1L, 2L), c(2L, 3L, 4L))), "partition")
  expect_error(phase17j_call(case, models = matrix(c(0, 0, 2, 1), 2, 2, byrow=TRUE)), "binary")
})

test_that("explicit seed and stable metadata are reproducible", {
  case <- phase17j_public_case(trait_specific = TRUE)
  a <- phase17j_call(case); b <- phase17j_call(case)
  expect_equal(a, b, tolerance = 1e-12)
  expect_identical(a$input$seed, case$x$csr$seed)
  expect_identical(a$input$marker_intersection_policy, "error")
  expect_identical(a$alignment$per_trait$allele_status, rep("exact", 2L))
  expect_true(all(c("trait_id", "study_id", "ancestry", "population", "ld_reference",
                    "sample_size", "ld_prefix", "ld_sharing_mode") %in% names(a$trait_metadata)))
})

test_that("public route remains singular and research routes are absent", {
  r <- paste(readLines(blr_repo_path("R", "mtblr-csr.R"), warn=FALSE), collapse="\n")
  cpp <- paste(readLines(blr_repo_path("src", "mtblr.cpp"), warn=FALSE), collapse="\n")
  expect_match(r, "mtblr_csr_raw_internal", fixed=TRUE)
  expect_false(grepl("mtblr_eigen|mtblr_cpg_omp_csr", r))
  expect_identical(source_match_count("for ( int it =", phase17i_src("src/blr_mt_default_core_impl.h"), fixed=TRUE), 1L)
  expect_identical(source_match_count("mtblr_csr_raw_internal(", cpp, fixed=TRUE), 1L)
})
