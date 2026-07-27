test_that("Phase 17S exposes and validates the exact public controls", {
  expected <- c(
    "y", "Glist", "covar", "chr", "cls", "rows", "scale", "center",
    "residual_covariance", "method", "trait_metadata", "sets", "block_size",
    "beta", "b", "state", "h2", "pi", "models", "pimodels",
    "mixture_var", "joint_pi", "joint_pi_prior", "component",
    "selection_s", "selection_maf", "allow_reference_maf_for_selection_s",
    "estimate_selection_s", "selection_s_init",
    "selection_s_prior", "selection_s_proposal_sd", "vg", "vb",
    "ve", "ssb_prior", "sse_prior", "updateB", "updateE", "updatePi",
    "nub", "nue", "nit", "nburn", "nthin", "seed", "nchains", "ncores",
    "chain_seeds", "keep_chains", "convergence", "convergence_control",
    "memory_warning_gb", "verbose")
  expect_identical(names(formals(mtblr_bed)), expected)
  expect_identical(formals(mtblr_bed)$nchains, 1L)
  expect_identical(formals(mtblr_bed)$ncores, 1L)
  expect_null(formals(mtblr_bed)$chain_seeds)
  expect_identical(formals(mtblr_bed)$keep_chains, FALSE)

  bad_count <- list(0, -1, 1.5, Inf, NA_real_, c(1, 2), 2^31)
  for (value in bad_count) {
    expect_error(sblr:::.mtblr_bed_chain_controls(value, 1, NULL, FALSE),
                 "nchains")
    expect_error(sblr:::.mtblr_bed_chain_controls(1, value, NULL, FALSE),
                 "ncores")
  }
  expect_error(sblr:::.mtblr_bed_chain_controls(1, 1, NULL, NA),
               "keep_chains")
  expect_error(sblr:::.mtblr_bed_chain_controls(1, 1, NULL, 1),
               "keep_chains")
  expect_error(sblr:::.mtblr_bed_chain_controls(1, 1, NULL, c(TRUE, FALSE)),
               "keep_chains")
  invalid_seeds <- list(1:2, NA_real_, NaN, Inf, 1.5,
                        -2147483649, 2147483648)
  for (seeds in invalid_seeds) {
    expect_error(sblr:::.mtblr_bed_chain_controls(1, 1, seeds, FALSE),
                 "chain_seeds")
  }
  explicit <- sblr:::.mtblr_bed_chain_controls(
    3, 9, c(-2147483648, 7, 7), TRUE)
  expect_identical(explicit$ncores, 9L)
  expect_identical(explicit$native, c(-2147483648, 7, 7))
  expect_identical(explicit$requested, c(-2147483648, 7, 7))
  expect_identical(sblr:::.mtblr_bed_chain_controls(
    2, 1, NULL, FALSE)$native, integer())
})

test_that("Phase 17S default calls reduce exactly to the serial reference", {
  for (nt in 1:3) {
    for (mode in c("full", "diagonal")) {
      for (updates in c(FALSE, TRUE)) {
        case <- phase17p_case(nt = nt)
        on.exit(phase17p_cleanup(case), add = TRUE)
        args <- phase17s_public_args(case, residual_covariance = mode,
                                     updates = updates)
        fit <- suppressWarnings(do.call(mtblr_bed, args))
        serial <- do.call(sblr:::mtblr_bed_internal,
                          phase17p_native_args(args))
        expected <- phase17s_raw_numerics(phase17s_internal(args))
        expect_equal(phase17s_fit_numerics(fit), expected, tolerance = 1e-12)
        expect_equal(unname(fit$bm), serial$marker$bm, tolerance = 0)
        expect_equal(unname(fit$dm), serial$marker$dm, tolerance = 0)
        expect_identical(unname(fit$d), serial$marker$state)
        expect_true(all(fit$bm_chain_mean_sd == 0) && all(fit$dm_chain_mean_sd == 0))
        expect_identical(fit$bm_chain_mean_min, fit$bm)
        expect_identical(fit$bm_chain_mean_max, fit$bm)
        expect_identical(fit$dm_chain_mean_min, fit$dm)
        expect_identical(fit$dm_chain_mean_max, fit$dm)
        expect_null(fit$chains)
      }
    }
  }
})

test_that("Phase 17S public and internal multichain routes agree", {
  case2 <- phase17p_case(nt = 2L)
  case3 <- phase17p_case(nt = 3L)
  on.exit(phase17p_cleanup(case2), add = TRUE)
  on.exit(phase17p_cleanup(case3), add = TRUE)
  configurations <- list(
    phase17s_public_args(case2, 1L, 4L, residual_covariance = "full"),
    phase17s_public_args(case2, 2L, 1L, residual_covariance = "diagonal"),
    phase17s_public_args(case2, 2L, 2L, updates = TRUE,
                         keep_chains = TRUE),
    phase17s_public_args(case2, 4L, 1L,
                         chain_seeds = c(-1, 19, 19, 31)),
    phase17s_public_args(case2, 4L, 3L, models = "restrictive"),
    phase17s_public_args(case3, 2L, 2L, updates = TRUE,
                         sets = list(c(1L, 3L, 5L), c(2L, 4L))))
  for (args in configurations) phase17s_compare_public_internal(args)
})

test_that("Phase 17S exposes pooled, primary, trace, and compact semantics", {
  case <- phase17p_case(nt = 2L)
  on.exit(phase17p_cleanup(case), add = TRUE)
  args <- phase17s_public_args(case, 2L, 2L, updates = TRUE,
                               keep_chains = TRUE)
  fit <- do.call(mtblr_bed, args)
  expect_identical(names(fit$chains), c("chain1", "chain2"))
  bm <- simplify2array(lapply(fit$chains, function(x) x$marker$bm))
  dm <- simplify2array(lapply(fit$chains, function(x) x$marker$dm))
  expect_equal(fit$bm, apply(bm, 1:2, mean), tolerance = 1e-15)
  expect_equal(fit$dm, apply(dm, 1:2, mean), tolerance = 1e-15)
  expect_equal(fit$bm_chain_mean_sd, apply(bm, 1:2, sd), tolerance = 1e-15)
  expect_equal(fit$dm_chain_mean_sd, apply(dm, 1:2, sd), tolerance = 1e-15)
  expect_identical(fit$b, fit$chains$chain1$marker$b)
  expect_identical(fit$d, fit$chains$chain1$marker$state)
  expect_identical(fit$cov_b_final, fit$chains$chain1$variance$vb)
  expect_identical(fit$pi_final, fit$chains$chain1$pi$final)
  expect_identical(fit$vgs,
                   (fit$chains$chain1$trace$vgs +
                      fit$chains$chain2$trace$vgs) / 2)
  forbidden <- c("r", "wy", "order", "models", "sets", "phenotype",
                 "sample_residual", "genetic_values", "packed_genotype",
                 "marker_maps", "bed_files")
  expect_false(any(forbidden %in% unlist(lapply(fit$chains, names))))
  expect_identical(fit$input$posterior_summary_policy,
                   "pooled_retained_samples")
  expect_identical(fit$input$final_state_policy, "primary_chain")
  expect_identical(fit$input$trace_policy, "iterationwise_chain_mean")
})

test_that("Phase 17S memory scales by shared, workers, results, and retention", {
  one <- sblr:::.mtblr_bed_memory_estimate(5, 5, 2, 4, 7)
  many <- sblr:::.mtblr_bed_memory_estimate(5, 5, 2, 4, 7, 4, 2, FALSE)
  retained <- sblr:::.mtblr_bed_memory_estimate(5, 5, 2, 4, 7, 4, 9, TRUE)
  expect_equal(unname(one$components_bytes),
               c(320, 80, 80, 80, 80, 40, 40, 200, 40, 192, 320, 80,
                 80, 560))
  expect_identical(one$shared_bytes, many$shared_bytes)
  expect_identical(many$requested_worker_count, 2L)
  expect_identical(retained$requested_worker_count, 4L)
  expect_equal(many$estimated_chain_results_bytes,
               4 * many$result_bytes_per_chain)
  expect_identical(many$estimated_retained_output_bytes, 0)
  expect_equal(retained$estimated_retained_output_bytes,
               4 * retained$retained_chain_bytes_per_chain)
  expect_true(all(is.finite(unlist(one[c(
    "estimated_total_gib", "execution_estimated_total_gib")]))))

  case <- phase17p_case()
  on.exit(phase17p_cleanup(case), add = TRUE)
  args <- phase17s_public_args(case, 2L, 2L, keep_chains = TRUE)
  args$memory_warning_gb <- 1e-12
  expect_warning(do.call(mtblr_bed, args),
                 "nchains=2, ncores=2, requested workers=2, keep_chains=TRUE")
  expect_warning(do.call(mtblr_bed, args), "not measured RSS")
  args$memory_warning_gb <- Inf
  expect_no_warning(do.call(mtblr_bed, args))
})

test_that("Phase 17S is reproducible and preserves BLAS environment", {
  case <- phase17p_case()
  on.exit(phase17p_cleanup(case), add = TRUE)
  args <- phase17s_public_args(case, 2L, 1L, updates = TRUE,
                               chain_seeds = c(101L, 9277L))
  variables <- c("OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
                 "VECLIB_MAXIMUM_THREADS", "OMP_NUM_THREADS",
                 "OMP_THREAD_LIMIT")
  before <- Sys.getenv(variables, unset = NA_character_)
  first <- do.call(mtblr_bed, args)
  second <- do.call(mtblr_bed, args)
  expect_identical(phase17s_without_timing(first),
                   phase17s_without_timing(second))
  args$ncores <- 2L
  parallel <- do.call(mtblr_bed, args)
  expect_identical(phase17s_fit_numerics(first),
                   phase17s_fit_numerics(parallel))
  args$keep_chains <- TRUE
  retained <- do.call(mtblr_bed, args)
  expect_identical(phase17s_fit_numerics(parallel),
                   phase17s_fit_numerics(retained))
  expect_identical(Sys.getenv(variables, unset = NA_character_), before)
  expect_identical(first$input$blas_thread_environment, as.list(before))

  skip_if_not_installed("callr")
  root <- if (exists("blr_has_source_tree", mode = "function") &&
              blr_has_source_tree()) blr_source_root else NULL
  fresh_args <- phase17s_public_args(case, 2L, 1L, updates = TRUE,
                                     chain_seeds = c(101L, 9277L))
  fresh <- callr::r(function(args, root) {
    if (!is.null(root)) pkgload::load_all(root, compile = FALSE, quiet = TRUE)
    else library(sblr)
    do.call(mtblr_bed, args)
  }, list(args = fresh_args, root = root))
  expect_identical(phase17s_fit_numerics(first),
                   phase17s_fit_numerics(fresh))
})

test_that("Phase 17S changes only the public R activation surface", {
  root <- blr_repo_path()
  skip_if(is.null(root), "source checkout unavailable")
  public <- readLines(file.path(root, "R", "mtblr-bed.R"), warn = FALSE)
  expect_true(any(grepl("mtblr_bed_chains_internal", public,
                        fixed = TRUE)))
  expect_true(any(grepl("mtblr_bed_convergence_trace_internal", public,
                        fixed = TRUE)))
  expect_false(any(grepl("mtblr_bed_internal(", public, fixed = TRUE)))
  expect_equal(sum(grepl(".mtblr_bed_memory_estimate <-", public,
                         fixed = TRUE)), 1L)
  expect_equal(sum(readLines(file.path(root, "NAMESPACE"), warn = FALSE) ==
                     "export(mtblr_bed)"), 1L)
  expect_identical(names(formals(mtblr_csr)), c(
    "stats", "Glist", "ld_prefix", "ld_metadata", "trait_metadata",
    "marker_policy", "sample_overlap", "method", "n", "sets", "beta", "b",
    "state", "h2", "pi", "models", "pimodels", "mixture_var", "joint_pi",
    "joint_pi_prior", "component", "selection_s", "selection_maf",
    "allow_reference_maf_for_selection_s", "estimate_selection_s",
    "selection_s_init", "selection_s_prior", "selection_s_proposal_sd",
    "vg", "vb", "ve", "ssb_prior",
    "sse_prior", "updateB", "updateE", "updatePi", "nub", "nue", "nit",
    "nburn", "nthin", "seed", "nchains", "ncores", "chain_seeds",
    "keep_chains", "convergence", "convergence_control",
    "memory_warning_gb", "verbose"))
  expect_identical(names(formals(mtblr_block_eigen)), c(
    "stats", "Glist", "block_start", "operator_sharing", "eigen_filter",
    "eigen_tau", "eigen_eta", "summary_reference", "trait_metadata",
    "marker_policy", "sample_overlap", "method", "n", "sets", "beta", "b",
    "state", "h2", "pi", "models", "pimodels", "mixture_var", "joint_pi",
    "joint_pi_prior", "component", "selection_s", "selection_maf",
    "allow_reference_maf_for_selection_s", "estimate_selection_s",
    "selection_s_init", "selection_s_prior", "selection_s_proposal_sd",
    "vg", "vb", "ve", "ssb_prior",
    "sse_prior", "updateB", "updateE", "updatePi", "nub", "nue", "nit",
    "nburn", "nthin", "seed", "nchains", "ncores", "chain_seeds",
    "keep_chains", "convergence", "convergence_control",
    "memory_warning_gb", "verbose"))
  expect_false("sblr" %in% getNamespaceExports("sblr"))
  expect_identical(names(formals(stblr_bed)), c(
    "y", "Glist", "method", "...", "nit", "nburn", "nthin", "seed",
    "nchains", "ncores", "chain_seeds", "keep_chains", "convergence",
    "convergence_control", "memory_warning_gb", "verbose"))
})
