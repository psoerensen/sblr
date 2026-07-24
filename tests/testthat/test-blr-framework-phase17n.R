test_that("Phase 17N independent BED decoder owns codes, padding, rows, and columns", {
  f <- phase17n_fixture()
  expect_identical(phase17n_bed_code(c(2, NA, 1, 0)), c(0L, 1L, 2L, 3L))
  expect_equal(file.info(f$paths)$size, 3 + 2 * c(4, 3))
  expect_equal(f$dosage[, 1], c(1, 1, 1, NA, 2))
  expect_equal(f$dosage[, 2], c(0, 2, 2, 1, 1))
  expect_equal(f$dosage[, 4], c(1, 2, 2, 1, NA))
  expect_identical(dim(f$dosage), c(5L, 5L))
})

test_that("Phase 17N PackedBedMatrix logical decoding matches the independent oracle", {
  f <- phase17n_fixture()
  native <- sblr:::test_read_bedfiles_to_dense_matrix(
    f$paths, f$n_bed, f$rows, f$cls
  )
  expect_equal(unname(native), f$dosage, tolerance = 0)
  packed <- sblr:::test_read_bedfiles_to_packed_matrix(
    f$paths, f$n_bed, f$rows, f$cls
  )
  expect_identical(packed$n_used, length(f$rows))
  expect_identical(packed$m, 5L)
  expect_equal(packed$nbytes_per_variant, 2)
  expect_equal(packed$stride, 64)
})

test_that("Phase 17N transformation and dense algebra oracles are exact", {
  f <- phase17n_fixture()
  raw <- phase17n_transform_genotypes(f$dosage, f$af, FALSE)
  expect_true(all(f$X[is.na(f$dosage)] == 0))
  for (j in seq_len(ncol(raw))) {
    expect_true(all(raw[is.na(f$dosage[, j]), j] == 2 * f$af[j]))
  }
  Y <- cbind(T1 = c(-2, -1, 0, 1, 2), T2 = c(1, -2, 2, -1, 0))
  effective <- rbind(c(.2, 0), c(0, -.1), c(.15, .05), c(0, .3), c(-.2, 0))
  R <- phase17n_rebuild_residual(Y, f$X, effective)
  expect_equal(crossprod(f$X), t(f$X) %*% f$X, tolerance = 1e-14)
  expect_equal(crossprod(f$X, Y), t(f$X) %*% Y, tolerance = 1e-14)
  expect_equal(R, Y - f$X %*% effective, tolerance = 1e-14)
  j <- 3L
  score <- phase17n_marker_score(f$X, R, j, effective[j, ])
  expect_equal(score, drop(crossprod(f$X[, j], R + tcrossprod(f$X[, j], effective[j, ]))),
               tolerance = 1e-14)
  delta <- c(.03, -.04)
  expect_equal(phase17n_update_residual(R, f$X[, j], delta),
               Y - f$X %*% effective - tcrossprod(f$X[, j], delta),
               tolerance = 1e-14)
  U <- phase17n_genetic_values(f$X, effective)
  expect_equal(unname(U), unname(Y - R), tolerance = 1e-14)
  expect_equal(phase17n_genetic_covariance(f$X, effective),
               crossprod(U) / nrow(U), tolerance = 1e-14)
  expect_equal(crossprod(R), t(R) %*% R, tolerance = 1e-14)
  expect_equal(diag(crossprod(R)), colSums(R^2), tolerance = 1e-14)
})

test_that("Phase 17N full-E conditional has explicit null and diagonal reductions", {
  B <- matrix(c(.8, .2, .2, .6), 2)
  E <- matrix(c(1.1, .25, .25, .9), 2)
  models <- rbind(null = c(0, 0), first = c(1, 0),
                  second = c(0, 1), both = c(1, 1))
  ans <- phase17n_marker_conditional(c(.7, -.4), 5.2, B, E, models,
                                     c(.7, .1, .1, .1))
  expect_equal(ans$models[[1]]$C, solve(B), tolerance = 1e-14)
  expect_equal(ans$models[[1]]$rhs, c(0, 0), tolerance = 0)
  expect_equal(ans$models[[1]]$mean, c(0, 0), tolerance = 1e-14)
  expect_equal(sum(ans$probabilities), 1, tolerance = 1e-14)

  Ed <- diag(diag(E))
  diagonal <- phase17n_marker_conditional(c(.7, -.4), 5.2, B, Ed, models,
                                          c(.7, .1, .1, .1))
  D <- diag(models["both", ])
  expect_equal(diagonal$models[[4]]$C,
               solve(B) + 5.2 * D %*% solve(Ed) %*% D, tolerance = 1e-14)
  expect_equal(diagonal$models[[4]]$rhs,
               drop(D %*% solve(Ed) %*% c(.7, -.4)), tolerance = 1e-14)
})

test_that("Phase 17N memory formulas distinguish aligned packed and dense storage", {
  z <- phase17n_memory_formula(101, 1000, 3, 8, 50)
  expect_identical(z$packed_owner, 64000)
  expect_identical(z$phenotype, 2424)
  expect_identical(z$sample_residual, z$phenotype)
  expect_identical(z$decoded_marker, 808)
  expect_gt(z$model_work, 0)
})

test_that("Phase 17N portable contract object makes Phase 17O unambiguous", {
  expect_identical(phase17n_contract$implementation_status, "audit_only_no_sampler")
  expect_identical(phase17n_contract$owner, "PackedBedMatrix")
  expect_identical(phase17n_contract$view, "BedPackedGenotypeView")
  expect_identical(phase17n_contract$residual_covariance,
                   "full_canonical_diagonal_reduction")
  expect_identical(phase17n_contract$raw_schema, "mtblr_raw_version_1")
  expect_identical(phase17n_contract$marker_wy, "X_transpose_Y")
  expect_identical(phase17n_contract$marker_r, "X_transpose_R_final")
  expect_identical(phase17n_contract$execution,
                   "serial_one_chain_fit_local_mt19937")
  expect_identical(phase17n_contract$missing_phenotypes, "unsupported")
  expect_identical(phase17n_contract$covariates,
                   "pre_adjusted_no_native_argument")
})

test_that("Phase 17N source architecture remains audit-only", {
  if (!exists("blr_repo_path", mode = "function")) {
    skip("source checkout helper is unavailable")
  }
  read <- function(path) paste(readLines(blr_repo_path(path), warn = FALSE), collapse = "\n")
  packed <- read("src/packed_bed.h")
  bayesc <- read("src/st_cpg_omp_individual_scheduled_chains.cpp")
  bayesr <- read("src/st_bed_bayesr_common.h")
  family <- read("src/blr_bed_family_types.h")
  public <- read("NAMESPACE")
  contract <- read("docs/dev/blr_mt_bed_internal_contract.md")
  expect_match(packed, "00.*dosage 2|case 0u: return 2.0", perl = TRUE)
  expect_match(bayesc, "map.val\\[1\\] = 0.0")
  expect_match(bayesc, "map.val\\[1\\] = 2.0 \\* p")
  expect_match(bayesc, "jbase \\+ 3 < n")
  expect_match(bayesc, "rows0\\[k\\]")
  expect_match(bayesr, "rows0\\[k\\]")
  expect_match(family, "const std::uint8_t\\* packed_markers")
  expect_match(family, "const PackedGenotype& storage")
  expect_false(grepl("export\\(mtblr_bed\\)", public))
  expect_false(exists("mtblr_bed", asNamespace("sblr"), inherits = FALSE))
  expect_match(contract, "Full-E marker conditional")
  expect_match(contract, "`mtblr_raw` version 1", fixed = TRUE)
  expect_match(contract, "complete finite")
  expect_match(contract, "pre-adjusted")
  expect_match(contract, "no individual-level MT\\s+sampler", ignore.case = TRUE)
})

test_that("current scalar BED interface remains trait-specific and rejects covariates", {
  expect_match(paste(deparse(body(stblr_bed)), collapse = "\n"),
               "covar is not currently supported")
  expect_false("mtblr_bed" %in% getNamespaceExports("sblr"))
})
