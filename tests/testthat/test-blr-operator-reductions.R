test_that("operator reduction policy classifies exact and non-comparable cases", {
  expect_identical(sblr:::.blr_operator_reduction_policy(
    "stblr", "bayesc", "csr", "block_eigen", filtered = FALSE),
    "floating_point_equivalent")
  expect_identical(sblr:::.blr_operator_reduction_policy(
    "stblr", "bayesc", "packed_bed", "csr", filtered = FALSE),
    "comparable_when_likelihood_contracts_match")
  expect_identical(sblr:::.blr_operator_reduction_policy(
    "mtblr", "bayesc", "packed_bed", "csr", residual = "full"),
    "not_comparable_residual_covariance")
})

.blr_reduction_st_common <- function() list(
  updateB = FALSE, updateE = FALSE, updatePi = FALSE,
  nit = 8L, nburn = 2L, nthin = 1L, seed = 42L,
  convergence = "none")

test_that("ST CSR and unfiltered block eigen execute BayesC reduction", {
  fixture <- blr_unified_fixture()
  on.exit(blr_unified_cleanup(fixture), add = TRUE)
  prefix <- blr_unified_exact_ld_prefix(fixture$dosage)
  on.exit(blr_unified_cleanup_prefix(prefix), add = TRUE)
  common <- c(.blr_reduction_st_common(), list(
    method = "sbayesc", pi_init = .5, pi_prior_mean = .5,
    pi_prior_strength = 2))
  csr <- do.call(stblr_csr, c(list(stats = fixture$stats,
                                    ld_prefix = prefix), common))
  block <- do.call(stblr_block_eigen, c(list(
    stats = fixture$stats, Glist = fixture$Glist, block_start = 1L,
    eigen_filter = "hard_truncate", eigen_tau = 0), common))
  for (field in c("bm", "dm", "vbs", "vgs", "ves", "vle", "vld",
                  "pi_final", "pi_mean")) {
    expect_equal(block[[field]], csr[[field]], tolerance = 1e-7,
                 info = field)
  }
  diagnostic_common <- common
  diagnostic_common$pi_init <- .001
  diagnostic_common$pi_prior_mean <- .001
  diagnostic_common$pi_prior_strength <- 500000
  diagnostic_common$nchains <- 2L
  diagnostic_common$chain_seeds <- c(101L, 202L)
  diagnostic_common$convergence <- "core"
  diagnostic_common$convergence_control <- list(warn = FALSE)
  csr_diagnostic <- do.call(stblr_csr, c(list(
    stats = fixture$stats, ld_prefix = prefix), diagnostic_common))
  block_diagnostic <- do.call(stblr_block_eigen, c(list(
    stats = fixture$stats, Glist = fixture$Glist, block_start = 1L,
    eigen_filter = "hard_truncate", eigen_tau = 0), diagnostic_common))
  block_summary <- block_diagnostic$convergence$summary
  csr_summary <- csr_diagnostic$convergence$summary
  block_summary$diagnostic_key <- NULL
  csr_summary$diagnostic_key <- NULL
  expect_equal(block_summary, csr_summary, tolerance = 1e-12)
})

test_that("ST CSR and unfiltered block eigen execute BayesR reduction", {
  fixture <- blr_unified_fixture()
  on.exit(blr_unified_cleanup(fixture), add = TRUE)
  prefix <- blr_unified_exact_ld_prefix(fixture$dosage)
  on.exit(blr_unified_cleanup_prefix(prefix), add = TRUE)
  common <- c(.blr_reduction_st_common(), list(
    method = "sbayesr", mixture_var = c(0, .01, .1),
    pi = c(.7, .2, .1), alpha = c(1, 1, 1)))
  csr <- do.call(stblr_csr, c(list(stats = fixture$stats,
                                    ld_prefix = prefix), common))
  block <- do.call(stblr_block_eigen, c(list(
    stats = fixture$stats, Glist = fixture$Glist, block_start = 1L,
    eigen_filter = "hard_truncate", eigen_tau = 0), common))
  for (field in c("bm", "dm", "vbs", "vgs", "ves", "vle", "vld",
                  "pi_final", "pi_mean", "dm_component_mean")) {
    expect_lte(max(abs(block[[field]] - csr[[field]]), na.rm = TRUE), 1e-4)
  }
  expect_equal(block$component_probabilities,
               csr$component_probabilities, tolerance = 1e-8)
})

test_that("ST CSR and block eigen expose identical native BayesR diagnostics", {
  fixture <- blr_unified_fixture()
  on.exit(blr_unified_cleanup(fixture), add = TRUE)
  prefix <- blr_unified_exact_ld_prefix(fixture$dosage)
  on.exit(blr_unified_cleanup_prefix(prefix), add = TRUE)
  common <- list(
    method = "sbayesr", mixture_var = c(0, .01, .1),
    pi = c(.7, .2, .1), alpha = c(1, 1, 1),
    updateB = FALSE, updateE = FALSE, updatePi = TRUE,
    nit = 6L, nburn = 2L, nthin = 1L, seed = 421L,
    nchains = 2L, ncores = 1L, keep_chains = TRUE,
    convergence = "extended", convergence_control = list(
      warn = FALSE, extended_groups = "probability",
      selected_markers = c(3L, 1L),
      selected_marker_quantities = c("b", "d", "component"),
      keep_traces = TRUE))
  csr <- do.call(stblr_csr, c(list(stats = fixture$stats,
                                    ld_prefix = prefix), common))
  block <- do.call(stblr_block_eigen, c(list(
    stats = fixture$stats, Glist = fixture$Glist, block_start = 1L,
    eigen_filter = "hard_truncate", eigen_tau = 0), common))
  csr_desc <- csr$convergence_traces$quantities
  block_desc <- block$convergence_traces$quantities
  expect_identical(block_desc$quantity, csr_desc$quantity)
  expect_equal(block$convergence_traces$values,
               csr$convergence_traces$values, tolerance = 1e-7)
  expect_identical(unique(block_desc$marker_id[block_desc$tier == 3L]),
                   fixture$stats$marker_names[c(3L, 1L)])
  component <- block$convergence_traces$values[, ,
    which(block_desc$parameter_name == "selected_component"), drop = FALSE]
  expect_true(all(component %in% 0:2))
})

test_that("ST CSR and unfiltered block eigen execute SBayesRC reduction", {
  fixture <- blr_unified_fixture()
  on.exit(blr_unified_cleanup(fixture), add = TRUE)
  prefix <- blr_unified_exact_ld_prefix(fixture$dosage)
  on.exit(blr_unified_cleanup_prefix(prefix), add = TRUE)
  annotation <- matrix(
    c(1, 0, 1, 1, 1, 0), 3L, 2L,
    dimnames = list(fixture$stats$marker_names, c("intercept", "annot1")))
  common <- .blr_reduction_st_common()
  common$updatePi <- NULL
  common <- c(common, list(updateAlpha = FALSE, h2 = .5))
  csr <- do.call(stblr_csr_annot, c(list(
    stats = fixture$stats, Glist = fixture$Glist, ld_prefix = prefix,
    annotations = annotation,
    annotation_model = "sbayesrc"), common))
  block <- do.call(stblr_block_eigen, c(list(
    stats = fixture$stats, Glist = fixture$Glist, block_start = 1L,
    method = "sbayesrc", annotation = annotation,
    eigen_filter = "hard_truncate", eigen_tau = 0), common))
  for (field in c("bm", "dm", "vbs", "vgs", "ves", "vle", "vld",
                  "pi_final", "pi_mean", "dm_component_mean")) {
    expect_equal(block[[field]], csr[[field]], tolerance = 1e-12,
                 info = field)
  }
  expect_equal(block$component_probabilities,
               csr$component_probabilities, tolerance = 1e-12)
})

.blr_reduction_mt_inputs <- function(fixture, prefix) {
  stats <- fixture$stats
  marker_metadata <- stats$marker_metadata
  marker_metadata$effect_allele <- "A"
  marker_metadata$other_allele <- "C"
  stats$marker_metadata <- marker_metadata
  list(
    stats = stats,
    metadata = list(
      prefix = prefix, marker_ids = stats$marker_names,
      marker_metadata = marker_metadata, scale = "standardized_genotype",
      source = "make_summary_stats"))
}

test_that("MT CSR and unfiltered block eigen execute BayesC reduction", {
  fixture <- blr_unified_fixture()
  on.exit(blr_unified_cleanup(fixture), add = TRUE)
  prefix <- blr_unified_exact_ld_prefix(fixture$dosage)
  on.exit(blr_unified_cleanup_prefix(prefix), add = TRUE)
  inputs <- .blr_reduction_mt_inputs(fixture, prefix)
  common <- list(
    vb = matrix(.1), ve = matrix(.5), models = matrix(c(0L, 1L), 2L, 1L),
    pimodels = c(0, 1), updateB = FALSE, updateE = FALSE,
    updatePi = FALSE, nit = 8L, nburn = 2L, nthin = 1L, seed = 42L,
    convergence = "none")
  csr <- do.call(mtblr_csr, c(list(
    stats = inputs$stats, ld_prefix = prefix,
    ld_metadata = inputs$metadata), common))
  block <- do.call(mtblr_block_eigen, c(list(
    stats = fixture$stats, Glist = fixture$Glist, block_start = 1L,
    eigen_filter = "hard_truncate", eigen_tau = 0), common))
  expect_equal(block$dm, csr$dm, tolerance = 0)
  expect_equal(block$pi_mean, csr$pi_mean, tolerance = 0)
  expect_lte(max(abs(block$bm - csr$bm)), .02)
  for (field in c("vbs", "vgs", "ves", "vle", "vld")) {
    expect_lte(max(abs(block[[field]] - csr[[field]])), .02)
  }
  expect_equal(csr$vld, csr$vgs - csr$vle, tolerance = 1e-12)
  expect_equal(block$vld, block$vgs - block$vle, tolerance = 1e-12)
})

test_that("BED and summary BayesC execute a valid null-state reduction", {
  fixture <- blr_unified_fixture()
  on.exit(blr_unified_cleanup(fixture), add = TRUE)
  prefix <- blr_unified_exact_ld_prefix(fixture$dosage)
  on.exit(blr_unified_cleanup_prefix(prefix), add = TRUE)
  common <- list(
    method = "sbayesc", pi_init = .001, updateB = FALSE,
    updateE = FALSE, updatePi = FALSE, nit = 8L, nburn = 2L,
    nthin = 1L, seed = 42L, convergence = "none")
  csr <- do.call(stblr_csr, c(list(stats = fixture$stats,
                                    ld_prefix = prefix), common))
  bed_common <- common
  bed_common$method <- "bayesc"
  bed <- do.call(stblr_bed, c(list(y = fixture$y,
                                    Glist = fixture$Glist), bed_common))
  for (field in c("bm", "dm", "vbs", "vgs", "ves", "vle", "vld",
                  "pi_final")) {
    expect_equal(bed[[field]], csr[[field]], tolerance = 1e-12,
                 info = field)
  }
})

test_that("one-trait MT and ST BayesC execute matched null-state reduction", {
  fixture <- blr_unified_fixture()
  on.exit(blr_unified_cleanup(fixture), add = TRUE)
  prefix <- blr_unified_exact_ld_prefix(fixture$dosage)
  on.exit(blr_unified_cleanup_prefix(prefix), add = TRUE)
  st <- stblr_csr(
    fixture$stats, ld_prefix = prefix, method = "sbayesc", pi_init = .001,
    updateB = FALSE, updateE = FALSE, updatePi = FALSE,
    nit = 8L, nburn = 2L, seed = 42L, convergence = "none")
  inputs <- .blr_reduction_mt_inputs(fixture, prefix)
  mt <- mtblr_csr(
    inputs$stats, ld_prefix = prefix, ld_metadata = inputs$metadata,
    vb = st$cov_b_final, ve = st$cov_e_final,
    models = matrix(c(0L, 1L), 2L, 1L), pimodels = c(.999, .001),
    updateB = FALSE, updateE = FALSE, updatePi = FALSE,
    nit = 8L, nburn = 2L, seed = 42L, convergence = "none")
  for (field in c("bm", "dm", "vgs", "vle", "vld")) {
    expect_equal(mt[[field]], st[[field]], tolerance = 1e-12,
                 info = field)
  }
  expect_equal(unname(mt$pi_final), unname(st$pi_final), tolerance = 0)
  expect_false(isTRUE(all.equal(mt$vbs, st$vbs)))
  expect_false(isTRUE(all.equal(mt$ves, st$ves)))
})

test_that("S models reuse BayesC and BayesR kernels across exact operators", {
  fixture <- blr_unified_fixture()
  on.exit(blr_unified_cleanup(fixture), add = TRUE)
  prefix <- blr_unified_exact_ld_prefix(fixture$dosage)
  on.exit(blr_unified_cleanup_prefix(prefix), add = TRUE)
  cases <- list(
    list(method = "sbayesc", pi_init = .5, pi_prior_mean = .5,
         pi_prior_strength = 2, selection_s = 0),
    list(method = "sbayesr", mixture_var = c(0, .01, .1),
         pi = c(.7, .2, .1), alpha = c(1, 1, 1), selection_s = 0))
  for (case in cases) {
    common <- c(.blr_reduction_st_common(), case)
    csr <- do.call(stblr_csr, c(list(
      stats = fixture$stats, Glist = fixture$Glist, ld_prefix = prefix),
      common))
    block <- do.call(stblr_block_eigen, c(list(
      stats = fixture$stats, Glist = fixture$Glist, block_start = 1L,
      eigen_filter = "hard_truncate", eigen_tau = 0), common))
    expect_identical(csr$model, case$method)
    expect_identical(block$model, case$method)
    expect_true(all(csr$input$effect_scale == "maf_s" |
                    csr$input$effect_scale == "component_maf_s"))
    for (field in c("bm", "dm", "vbs", "vgs", "ves", "vle", "vld",
                    "pi_final", "pi_mean")) {
      expect_equal(block[[field]], csr[[field]], tolerance = 1e-4,
                   info = paste(case$method, field))
    }
  }
})
