.mt_bayesr_reduction_inputs <- function(fixture, prefix) {
  stats <- fixture$stats
  metadata <- stats$marker_metadata
  metadata$effect_allele <- "A"
  metadata$other_allele <- "C"
  stats$marker_metadata <- metadata
  list(stats = stats, metadata = list(
    prefix = prefix, marker_ids = stats$marker_names,
    marker_metadata = metadata, scale = "standardized_genotype",
    source = "make_summary_stats"))
}

.mt_bayesr_common <- function(method = "sbayesr", selection_s = NULL) {
  list(method = method, mixture_var = c(0, .1, 1),
       models = matrix(c(0L, 1L), 2L, 1L),
       joint_pi = c(.7, .15, .15), joint_pi_prior = rep(1, 3),
       selection_s = selection_s, vb = matrix(.1), ve = matrix(.5),
       updateB = FALSE, updateE = FALSE, updatePi = TRUE,
       nit = 8L, nburn = 2L, nthin = 1L, seed = 42L,
       convergence = "none")
}

test_that("MT BayesR executes in CSR and exact block eigen", {
  fixture <- blr_unified_fixture()
  on.exit(blr_unified_cleanup(fixture), add = TRUE)
  prefix <- blr_unified_exact_ld_prefix(fixture$dosage)
  on.exit(blr_unified_cleanup_prefix(prefix), add = TRUE)
  inputs <- .mt_bayesr_reduction_inputs(fixture, prefix)
  common <- .mt_bayesr_common()
  csr <- do.call(mtblr_csr, c(list(
    stats = inputs$stats, ld_prefix = prefix,
    ld_metadata = inputs$metadata), common))
  block <- do.call(mtblr_block_eigen, c(list(
    stats = fixture$stats, Glist = fixture$Glist, block_start = 1L,
    eigen_filter = "hard_truncate", eigen_tau = 0), common))
  expect_identical(csr$model, "sbayesr")
  expect_identical(block$model, "sbayesr")
  for (fit in list(csr, block)) {
    expect_identical(length(fit$component_final), nrow(fit$bm))
    expect_equal(unname(rowSums(fit$component_probabilities)), rep(1, nrow(fit$bm)),
                 tolerance = 1e-12)
    expect_equal(fit$vld, fit$vgs - fit$vle, tolerance = 1e-12)
    expect_identical(colnames(fit$pi_trace), fit$input$model_names)
    expect_equal(sum(fit$model_parameters$mixture$component_pi_mean), 1,
                 tolerance = 1e-12)
    expect_equal(sum(fit$model_parameters$mixture$pattern_pi_mean), 1,
                 tolerance = 1e-12)
    expect_identical(colnames(fit$component_probabilities),
                     fit$model_parameters$mixture$component_names)
  }
  for (field in c("bm", "dm", "vbs", "vgs", "ves", "vle", "vld",
                  "pi_final", "pi_mean", "component_probabilities"))
    expect_equal(block[[field]], csr[[field]], tolerance = 1e-7, info = field)
})

test_that("MT explicit selection_s=-1 uses the unit scale in every operator", {
  fixture <- blr_unified_fixture()
  on.exit(blr_unified_cleanup(fixture), add = TRUE)
  prefix <- blr_unified_exact_ld_prefix(fixture$dosage)
  on.exit(blr_unified_cleanup_prefix(prefix), add = TRUE)
  inputs <- .mt_bayesr_reduction_inputs(fixture, prefix)
  calls <- list(
    csr = list(fun = mtblr_csr, base = list(stats = inputs$stats,
      ld_prefix = prefix, ld_metadata = inputs$metadata)),
    block_eigen = list(fun = mtblr_block_eigen, base = list(
      stats = fixture$stats, Glist = fixture$Glist, block_start = 1L,
      eigen_filter = "hard_truncate", eigen_tau = 0)),
    packed_bed = list(fun = mtblr_bed, base = list(y = fixture$y,
      Glist = fixture$Glist, residual_covariance = "diagonal")))
  fields <- c("bm", "dm", "vbs", "vgs", "ves", "vle", "vld",
    "pi_final", "pi_mean", "component_final", "component_probabilities")
  for (entry in calls) {
    public_model <- if (identical(entry$fun, mtblr_bed)) "bayesr" else "sbayesr"
    bayesr <- do.call(entry$fun, c(entry$base,
      .mt_bayesr_common(public_model)))
    sbayesr <- do.call(entry$fun, c(entry$base,
      .mt_bayesr_common(public_model, -1)))
    for (field in fields)
      expect_equal(sbayesr[[field]], bayesr[[field]], tolerance = 0,
                   info = field)
  }
})

test_that("MT BayesR update controls and analytical mixture memory are explicit", {
  fixture <- blr_unified_fixture()
  on.exit(blr_unified_cleanup(fixture), add = TRUE)
  common <- .mt_bayesr_common("bayesr")
  common$updatePi <- FALSE
  common$updateB <- FALSE
  common$updateE <- FALSE
  fit <- do.call(mtblr_bed, c(list(y = fixture$y, Glist = fixture$Glist,
    residual_covariance = "diagonal"), common))
  expect_equal(unname(fit$pi_final), common$joint_pi)
  expect_equal(unname(fit$pi_mean), common$joint_pi)
  expect_gt(fit$memory_estimate$bayesr_components_bytes, 0)
  expect_true(all(c("bayesr_shared_state_descriptors_bytes",
    "bayesr_private_worker_state_bytes", "bayesr_chain_result_bytes",
    "bayesr_component_output_bytes") %in% names(fit$memory_estimate)))
})

test_that("one-trait ST and MT summary BayesR execute matched null reductions", {
  fixture <- blr_unified_fixture()
  on.exit(blr_unified_cleanup(fixture), add = TRUE)
  prefix <- blr_unified_exact_ld_prefix(fixture$dosage)
  on.exit(blr_unified_cleanup_prefix(prefix), add = TRUE)
  inputs <- .mt_bayesr_reduction_inputs(fixture, prefix)
  probability <- c(1 - 2e-12, 1e-12, 1e-12)
  for (selection in list(NULL, -1)) {
    method <- "sbayesr"
    s_control <- if (is.null(selection)) list() else list(selection_s = selection)
    st <- do.call(stblr_csr, c(list(stats = inputs$stats,
      Glist = fixture$Glist, ld_prefix = prefix, method = method,
      mixture_var = c(0, .1, 1),
      pi = probability, alpha = c(1, 1, 1), updateB = FALSE,
      updateE = FALSE, updatePi = FALSE, nit = 8L, nburn = 2L,
      seed = 42L, convergence = "none"), s_control))
    mt <- do.call(mtblr_csr, c(list(stats = inputs$stats,
      ld_prefix = prefix, ld_metadata = inputs$metadata, method = method,
      mixture_var = c(0, .1, 1), joint_pi = probability,
      models = matrix(c(0L, 1L), 2L, 1L), vb = st$cov_b_final,
      ve = st$cov_e_final, updateB = FALSE, updateE = FALSE,
      updatePi = FALSE, nit = 8L, nburn = 2L, seed = 42L,
      convergence = "none"), s_control))
    for (field in c("bm", "dm", "vgs", "vle", "vld"))
      expect_equal(mt[[field]], st[[field]], tolerance = 0,
                   info = paste(method, field))
    st_component <- if (is.list(st$component_probabilities))
      st$component_probabilities[[1L]] else st$component_probabilities
    expect_equal(unname(mt$component_probabilities), unname(st_component),
                 tolerance = 0)
    expect_equal(unname(mt$pi_final), as.numeric(st$pi_final), tolerance = 0)
  }
})

test_that("MT packed BED and summary BayesR execute a valid null-state reduction", {
  fixture <- blr_unified_fixture()
  on.exit(blr_unified_cleanup(fixture), add = TRUE)
  prefix <- blr_unified_exact_ld_prefix(fixture$dosage)
  on.exit(blr_unified_cleanup_prefix(prefix), add = TRUE)
  inputs <- .mt_bayesr_reduction_inputs(fixture, prefix)
  probability <- c(1 - 2e-12, 1e-12, 1e-12)
  summary_common <- list(method = "sbayesr", mixture_var = c(0, .1, 1),
    joint_pi = probability, models = matrix(c(0L, 1L), 2L, 1L),
    vb = matrix(.1), ve = matrix(.5), updateB = FALSE, updateE = FALSE,
    updatePi = FALSE, nit = 8L, nburn = 2L, seed = 42L,
    convergence = "none")
  csr <- do.call(mtblr_csr, c(list(stats = inputs$stats,
    ld_prefix = prefix, ld_metadata = inputs$metadata), summary_common))
  individual_common <- summary_common
  individual_common$method <- "bayesr"
  bed <- do.call(mtblr_bed, c(list(y = fixture$y, Glist = fixture$Glist,
    residual_covariance = "diagonal"), individual_common))
  for (field in c("bm", "dm", "vbs", "vgs", "ves", "vle", "vld",
                  "pi_final", "pi_mean", "component_probabilities"))
    expect_equal(bed[[field]], csr[[field]], tolerance = 0, info = field)
})

test_that("MT BayesR logical chains are reproducible and retention-independent", {
  fixture <- blr_unified_fixture()
  on.exit(blr_unified_cleanup(fixture), add = TRUE)
  prefix <- blr_unified_exact_ld_prefix(fixture$dosage)
  on.exit(blr_unified_cleanup_prefix(prefix), add = TRUE)
  inputs <- .mt_bayesr_reduction_inputs(fixture, prefix)
  common <- .mt_bayesr_common()
  common$nchains <- 2L; common$chain_seeds <- c(101L, 202L)
  common$ncores <- 1L; common$convergence <- "core"
  common$convergence_control <- list(warn = FALSE, keep_traces = TRUE)
  call <- function(keep, cores) {
    args <- c(list(stats = inputs$stats, ld_prefix = prefix,
                   ld_metadata = inputs$metadata), common)
    args$keep_chains <- keep; args$ncores <- cores
    do.call(mtblr_csr, args)
  }
  serial <- call(FALSE, 1L); repeated <- call(FALSE, 1L)
  retained <- call(TRUE, 1L); parallel <- call(FALSE, 2L)
  fields <- c("bm", "dm", "vbs", "vgs", "ves", "vle", "vld",
    "pi_final", "pi_mean", "component_final", "component_probabilities",
    "convergence", "convergence_traces")
  for (field in fields) {
    expect_equal(repeated[[field]], serial[[field]], tolerance = 0,
                 info = paste("repeat", field))
    expect_equal(retained[[field]], serial[[field]], tolerance = 0,
                 info = paste("retention", field))
    expect_equal(parallel[[field]], serial[[field]], tolerance = 0,
                 info = paste("workers", field))
  }
  expect_null(serial$chains)
  expect_identical(length(retained$chains), 2L)
})

test_that("MT summary BayesR selection_s=-1 reduces exactly to NULL", {
  fixture <- blr_unified_fixture()
  on.exit(blr_unified_cleanup(fixture), add = TRUE)
  prefix <- blr_unified_exact_ld_prefix(fixture$dosage)
  on.exit(blr_unified_cleanup_prefix(prefix), add = TRUE)
  inputs <- .mt_bayesr_reduction_inputs(fixture, prefix)
  bayesr <- do.call(mtblr_csr, c(list(
    stats = inputs$stats, ld_prefix = prefix,
    ld_metadata = inputs$metadata), .mt_bayesr_common()))
  sbayesr <- do.call(mtblr_csr, c(list(
    stats = inputs$stats, ld_prefix = prefix,
    ld_metadata = inputs$metadata), .mt_bayesr_common("sbayesr", -1)))
  for (field in c("bm", "dm", "vbs", "vgs", "ves", "vle", "vld",
                  "pi_final", "pi_mean", "component_final",
                  "component_probabilities"))
    expect_equal(sbayesr[[field]], bayesr[[field]], tolerance = 0, info = field)
})

test_that("MT BayesR multichain convergence and compact chains are coherent", {
  fixture <- blr_unified_fixture()
  on.exit(blr_unified_cleanup(fixture), add = TRUE)
  prefix <- blr_unified_exact_ld_prefix(fixture$dosage)
  on.exit(blr_unified_cleanup_prefix(prefix), add = TRUE)
  inputs <- .mt_bayesr_reduction_inputs(fixture, prefix)
  common <- .mt_bayesr_common()
  common$nchains <- 2L
  common$ncores <- 2L
  common$chain_seeds <- c(101L, 202L)
  common$keep_chains <- TRUE
  common$convergence <- "core"
  common$convergence_control <- list(warn = FALSE, keep_traces = TRUE)
  fit <- do.call(mtblr_csr, c(list(
    stats = inputs$stats, ld_prefix = prefix,
    ld_metadata = inputs$metadata), common))
  expect_identical(nrow(fit$convergence$summary), 5L)
  expect_identical(fit$convergence$summary$quantity,
                   c("vbs[T1]", "vgs[T1]", "ves[T1]",
                     "vle[T1]", "vld[T1]"))
  expect_identical(length(fit$chains), 2L)
  expect_true(all(vapply(fit$chains, function(chain)
    !is.null(chain$marker$component_probabilities), logical(1))))
  expect_identical(dim(fit$convergence_traces$values), c(8L, 2L, 5L))
})

test_that("MT BayesR packed BED returns coherent component output", {
  fixture <- blr_unified_fixture()
  on.exit(blr_unified_cleanup(fixture), add = TRUE)
  common <- .mt_bayesr_common("bayesr")
  common$vb <- matrix(.1)
  common$ve <- matrix(.5)
  fit <- do.call(mtblr_bed, c(list(y = fixture$y, Glist = fixture$Glist,
    residual_covariance = "diagonal"), common))
  expect_identical(fit$model, "bayesr")
  expect_equal(unname(rowSums(fit$component_probabilities)),
               rep(1, nrow(fit$bm)), tolerance = 1e-12)
  expect_equal(fit$vld, fit$vgs - fit$vle, tolerance = 1e-12)
  expect_identical(length(fit$component_final), nrow(fit$bm))
})
