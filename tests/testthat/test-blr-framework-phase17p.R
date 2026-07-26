test_that("Phase 17P public and internal routes are numerically identical", {
  configurations <- list(
    list(nt = 2L, mode = "full", updates = FALSE),
    list(nt = 2L, mode = "full", updates = TRUE),
    list(nt = 2L, mode = "diagonal", updates = FALSE),
    list(nt = 2L, mode = "diagonal", updates = TRUE),
    list(nt = 3L, mode = "full", updates = TRUE),
    list(nt = 1L, mode = "diagonal", updates = TRUE)
  )
  for (spec in configurations) {
    case <- phase17p_case(nt = spec$nt)
    on.exit(phase17p_cleanup(case), add = TRUE)
    args <- phase17p_public_args(
      case, spec$mode, spec$updates, center = FALSE)
    phase17p_compare_public_internal(args)
  }

  case <- phase17p_case()
  on.exit(phase17p_cleanup(case), add = TRUE)
  m <- sum(lengths(case$fixture$cls))
  restrictive <- rbind(c(0L, 0L), c(1L, 1L))
  phase17p_compare_public_internal(phase17p_public_args(
    case, "full", TRUE, models = restrictive, pimodels = c(.7, .3)))
  phase17p_compare_public_internal(phase17p_public_args(
    case, "diagonal", TRUE, sets = list(c(1L, 3L, 5L), c(2L, 4L))))
  state <- beta <- b <- matrix(0, m, 2L)
  state[1L, ] <- 1L; beta[1L, ] <- b[1L, ] <- c(.04, -.03)
  phase17p_compare_public_internal(phase17p_public_args(
    case, "full", FALSE, beta = beta, b = b, state = state))
  beta[2L, ] <- c(.06, -.05)
  phase17p_compare_public_internal(phase17p_public_args(
    case, "full", FALSE, beta = beta, b = b, state = state))
  phase17p_compare_public_internal(phase17p_public_args(
    case, "full", FALSE, cls = case$fixture$cls))

  matched <- phase17p_case(matched_ids = TRUE)
  on.exit(phase17p_cleanup(matched), add = TRUE)
  phase17p_compare_public_internal(phase17p_public_args(
    matched, "full", FALSE))
  uncentered <- phase17p_case(uncentered = TRUE)
  on.exit(phase17p_cleanup(uncentered), add = TRUE)
  phase17p_compare_public_internal(phase17p_public_args(
    uncentered, "full", FALSE, center = TRUE))
})

test_that("Phase 17P reuses BED alignment and reports stable provenance", {
  case <- phase17p_case(matched_ids = TRUE)
  on.exit(phase17p_cleanup(case), add = TRUE)
  extra <- data.frame(trait_id = colnames(case$Y), cohort = c("A", "B"))
  fit <- do.call(mtblr_bed, phase17p_public_args(
    case, "full", FALSE, trait_metadata = extra))
  expect_equal(
    unname(fit$alignment[c("individual_policy", "row_selection_status",
      "sample_order_status", "phenotype_missingness_status",
      "marker_selection_status", "marker_order_status",
      "genotype_orientation_status", "genotype_scale_status")]),
    list("shared_individual_level", "matched_by_id",
      "phenotype_order_preserved", "complete", "default_rsidsLD",
      "selected_glist_order_preserved", "by_construction_same_glist",
      "standardized_genotype"))
  expect_identical(fit$alignment$selected_rows, case$fixture$rows)
  expect_identical(fit$trait_metadata$cohort, c("A", "B"))
  expect_named(fit$marker_metadata,
    c("marker_id", "chromosome_or_file", "bed_column", "allele_frequency"))
  expect_identical(fit$marker_metadata$marker_id,
    c("a4", "a2", "a1", "b3", "b1"))

  unmatched <- case
  unmatched$Y <- rbind(case$Y, missing = c(1, 2))
  expect_warning(
    fit2 <- do.call(mtblr_bed, phase17p_public_args(
      unmatched, "full", FALSE)), "will be dropped")
  expect_identical(fit2$alignment$unmatched_input_ids, "missing")

  all_rows <- phase17p_case(use_all_rows = TRUE)
  on.exit(phase17p_cleanup(all_rows), add = TRUE)
  fit3 <- do.call(mtblr_bed, phase17p_public_args(
    all_rows, "full", FALSE))
  expect_identical(fit3$alignment$row_selection_status, "all_glist_rows")
  expect_identical(fit3$input$set_source, "chromosome_or_file")
})

test_that("Phase 17P accepts documented phenotype forms and names", {
  case <- phase17p_case(nt = 1L, matched_ids = TRUE)
  on.exit(phase17p_cleanup(case), add = TRUE)
  y <- setNames(case$Y[, 1L], rownames(case$Y))
  fit <- mtblr_bed(y, case$Glist, center = FALSE,
    updateB = FALSE, updateE = FALSE, updatePi = FALSE,
    nit = 2L, nburn = 0L, memory_warning_gb = Inf)
  expect_identical(colnames(fit$bm), "T1")
  frame <- as.data.frame(cbind(A = case$Y[, 1L], B = -case$Y[, 1L]))
  rownames(frame) <- rownames(case$Y)
  fit2 <- mtblr_bed(frame, case$Glist, center = FALSE,
    updateB = FALSE, updateE = FALSE, updatePi = FALSE,
    nit = 2L, nburn = 0L, memory_warning_gb = Inf)
  expect_identical(colnames(fit2$bm), c("A", "B"))
  unnamed <- unname(case$Y)
  fit3 <- mtblr_bed(unnamed, case$Glist, rows = case$fixture$rows,
    center = FALSE, updateB = FALSE, updateE = FALSE, updatePi = FALSE,
    nit = 2L, nburn = 0L, memory_warning_gb = Inf)
  expect_identical(colnames(fit3$bm), "T1")
})

test_that("Phase 17P centers only after sample alignment", {
  case <- phase17p_case(matched_ids = TRUE, uncentered = TRUE)
  on.exit(phase17p_cleanup(case), add = TRUE)
  original_variance <- apply(case$Y, 2L, var)
  fit <- do.call(mtblr_bed, phase17p_public_args(
    case, "full", FALSE, center = TRUE))
  prep <- fit$phenotype_preprocessing
  expect_equal(unname(prep$mean_before), unname(colMeans(case$Y)),
               tolerance = 1e-15)
  expect_equal(unname(prep$mean_after), c(0, 0), tolerance = 1e-15)
  expect_equal(unname(prep$variance_after), unname(original_variance),
               tolerance = 1e-15)
  expect_true(prep$center_applied)
  expect_identical(prep$centering_status,
                   rep("centered_by_adapter", 2L))
  expect_error(do.call(mtblr_bed, phase17p_public_args(
    case, "full", FALSE, center = FALSE)), "already centered")
})

test_that("Phase 17P covariance modes and priors follow the public contract", {
  case <- phase17p_case()
  on.exit(phase17p_cleanup(case), add = TRUE)
  full <- matrix(c(1.1, .12, .12, .9), 2L)
  prior <- matrix(c(.4, .05, .05, .35), 2L)
  fit <- do.call(mtblr_bed, phase17p_public_args(
    case, "full", FALSE, ve = full, sse_prior = prior))
  expect_equal(fit$input$ve, full)
  expect_equal(fit$input$sse_prior, prior)
  expect_equal(unname(fit$ve), full)
  expect_error(do.call(mtblr_bed, phase17p_public_args(
    case, "diagonal", FALSE, ve = full)), "exactly diagonal")
  expect_error(do.call(mtblr_bed, phase17p_public_args(
    case, "diagonal", FALSE, sse_prior = prior)), "exactly diagonal")
  expect_error(do.call(mtblr_bed, phase17p_public_args(
    phase17p_case(nt = 5L), "full", FALSE)), "max\\(2, nt - 1\\)")
})

test_that("Phase 17P initialization accepts latent/effective state contracts", {
  case <- phase17p_case()
  on.exit(phase17p_cleanup(case), add = TRUE)
  m <- sum(lengths(case$fixture$cls))
  b <- matrix(0, m, 2L); b[1L, ] <- c(.02, -.03)
  expect_no_error(suppressWarnings(do.call(mtblr_bed, phase17p_public_args(
    case, "full", FALSE, b = b))))
  state <- matrix(0L, m, 2L); state[1L, ] <- 1L
  expect_no_error(suppressWarnings(do.call(mtblr_bed, phase17p_public_args(
    case, "full", FALSE, state = state))))
  beta <- matrix(0, m, 2L); beta[2L, ] <- c(.4, -.2)
  expect_silent(do.call(mtblr_bed, phase17p_public_args(
    case, "full", FALSE, beta = list(beta[, 1L], beta[, 2L]),
    b = list(rep(0, m), rep(0, m)),
    state = list(rep(0L, m), rep(0L, m)))))
  expect_error(do.call(mtblr_bed, phase17p_public_args(
    case, "full", FALSE, beta = matrix(0, m - 1L, 2L))), "m by nt")
  bad_state <- state; bad_state[1L, 1L] <- 2L
  expect_error(do.call(mtblr_bed, phase17p_public_args(
    case, "full", FALSE, state = bad_state)), "binary")
  bad_b <- matrix(0, m, 2L); bad_b[1L, 1L] <- .1
  expect_error(do.call(mtblr_bed, phase17p_public_args(
    case, "full", FALSE, b = bad_b, models = "restrictive")),
    "model pattern")
})

test_that("Phase 17P memory estimate is analytical and warning-only", {
  estimate <- sblr:::.mtblr_bed_memory_estimate(5L, 5L, 2L, 4L, 7L)
  expected <- c(320, 80, 80, 80, 80, 40, 40, 200, 40, 192, 320, 80,
                80, 336)
  expect_equal(unname(estimate$components_bytes), expected)
  expect_equal(estimate$estimated_total_bytes, sum(expected))
  expect_false(estimate$measured_rss)
  expect_false(estimate$measured_peak_rss)
  case <- phase17p_case()
  on.exit(phase17p_cleanup(case), add = TRUE)
  args <- phase17p_public_args(case, "full", FALSE)
  args$memory_warning_gb <- 1e-12
  expect_warning(do.call(mtblr_bed, args), "not measured peak RSS")
  args$memory_warning_gb <- Inf
  expect_no_warning(do.call(mtblr_bed, args))
  args$memory_warning_gb <- 0
  expect_error(do.call(mtblr_bed, args), "positive finite scalar or Inf")
})

test_that("Phase 17P rejects unsupported or malformed public input", {
  case <- phase17p_case()
  on.exit(phase17p_cleanup(case), add = TRUE)
  call <- function(...) do.call(mtblr_bed,
    phase17p_public_args(case, "full", FALSE, ...))
  expect_error(mtblr_bed(case$Y), "Glist")
  expect_error(mtblr_bed(case$Y, list(case$Glist, case$Glist)), "one Glist")
  expect_error(call(covar = matrix(1, nrow(case$Y), 1L)), "does not currently")
  expect_error(call(scale = FALSE), "scale = TRUE")
  expect_error(call(method = "bayesR"), "Only method")
  bad <- case; bad$Y[1L, 1L] <- NA_real_
  expect_error(do.call(mtblr_bed, phase17p_public_args(
    bad, "full", FALSE)), "complete finite")
  bad <- case; bad$Y[, 1L] <- 1
  expect_error(do.call(mtblr_bed, phase17p_public_args(
    bad, "full", FALSE, center = TRUE)), "positive finite variance")
  bad <- case; colnames(bad$Y) <- c("x", "x")
  expect_error(do.call(mtblr_bed, phase17p_public_args(
    bad, "full", FALSE)), "Trait names")
  bad <- case; bad$Glist$af[[1L]][4L] <- 0
  expect_error(do.call(mtblr_bed, phase17p_public_args(
    bad, "full", FALSE)), "strictly inside")
  bad <- case; bad$Glist$rsids[[1L]][4L] <- "a2"
  expect_error(do.call(mtblr_bed, phase17p_public_args(
    bad, "full", FALSE)), "not found|Marker IDs")
})

test_that("Phase 17P fit and raw metadata remain bounded", {
  case <- phase17p_case()
  on.exit(phase17p_cleanup(case), add = TRUE)
  fit <- do.call(mtblr_bed, phase17p_public_args(
    case, "full", TRUE, models = "restrictive"))
  expect_s3_class(fit, "mtblr_fit")
  expect_equal(fit$raw_schema_version, 1L)
  expect_equal(unname(fit$input[c("backend", "data_level", "residual_covariance",
    "genotype_scale", "phenotype_scaling", "cpo", "le_ld")]),
    list("mt_bed_bayesc", "individual", "full", "standardized_genotype",
         "not_performed", "unsupported", "unsupported"))
  expect_true(all(c("bed_diagnostics", "phenotype_preprocessing",
                    "memory_estimate") %in% names(fit)))
  expect_false(any(c("sample_residual", "genetic_values", "phenotype",
                     "packed_bytes", "raw_pointer") %in% names(fit)))
  expect_true(all(is.finite(fit$re)))
})

test_that("Phase 17P public architecture is singular and protected", {
  expect_true(is.function(mtblr_bed))
  expect_true("mtblr_bed" %in% getNamespaceExports("sblr"))
  expected_formals <- c(
    "y", "Glist", "covar", "chr", "cls", "rows", "scale", "center",
    "residual_covariance", "method", "trait_metadata", "sets",
    "block_size", "beta", "b", "state", "h2", "pi", "models",
    "pimodels", "vg", "vb", "ve", "ssb_prior", "sse_prior",
    "updateB", "updateE", "updatePi", "nub", "nue", "nit", "nburn",
    "nthin", "seed", "nchains", "ncores", "chain_seeds", "keep_chains",
    "memory_warning_gb", "verbose")
  expect_identical(names(formals(mtblr_bed)), expected_formals)
  root <- blr_repo_path()
  skip_if(is.null(root), "source architecture requires repository source")
  source <- readLines(file.path(root, "R", "mtblr-bed.R"), warn = FALSE)
  expect_equal(sum(grepl("mtblr_bed <- function", source, fixed = TRUE)), 1L)
  expect_equal(sum(grepl("raw <- mtblr_bed_chains_internal(", source,
                         fixed = TRUE)), 1L)
  expect_equal(sum(grepl(".as_mtblr_fit(", source, fixed = TRUE)), 1L)
  expect_true(any(grepl(".make_bed_marker_data(", source, fixed = TRUE)))
  expect_false(any(grepl("mtblr_eigen(", source, fixed = TRUE)))
  expect_false(any(grepl("mtblr_eigen(", source, fixed = TRUE)))
})
