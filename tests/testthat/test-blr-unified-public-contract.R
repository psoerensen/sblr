test_that("unified public fitting routes and identifiers are canonical", {
  exports <- getNamespaceExports("sblr")
  expect_true(all(c(
    "stblr_csr", "stblr_block_eigen", "stblr_bed", "stblr_csr_annot",
    "mtblr_bed", "mtblr_csr", "mtblr_block_eigen",
    "extract_posterior", "extract_diagnostics") %in% exports))
  expect_false(any(c(
    "sblr", "stblr_bed_marker", "stblr_csr_bayesr",
    "stblr_csr_prior_annot", "stblr_csr_learn_annot",
    "stblr_csr_group_annot", "stblr_csr_sbayesrc_generic",
    "check_stblr_convergence") %in% exports))
  for (fun in list(stblr_csr, stblr_block_eigen, stblr_bed)) {
    args <- names(formals(fun))
    expect_true(all(c("nit", "nburn", "nthin", "seed", "nchains",
                      "ncores", "chain_seeds", "keep_chains",
                      "convergence", "convergence_control",
                      "memory_warning_gb", "verbose") %in% args))
  }
  mt_args <- names(formals(mtblr_bed))
  expect_true(all(c(
    "nit", "nburn", "nthin", "seed", "nchains", "ncores",
    "chain_seeds", "keep_chains", "convergence", "convergence_control",
    "keep_traces", "memory_limit_bytes") %in% mt_args))
  for (fun in list(mtblr_csr, mtblr_block_eigen)) {
    summary_args <- names(formals(fun))
    expect_true(all(c(
      "providers", "operator_resources", "global_marker_ids",
      "global_alleles", "trait_ids", "method", "nit", "nburn",
      "nthin", "seed", "nchains", "ncores", "chain_seeds",
      "keep_chains", "convergence", "convergence_control",
      "keep_traces", "memory_limit_bytes") %in% summary_args))
  }
  expect_identical(eval(formals(stblr_csr)$method),
                   c("sbayesc", "sbayesr"))
  expect_identical(eval(formals(stblr_bed)$method),
                   c("bayesc", "bayesr", "bayesrc"))
  expect_identical(eval(formals(stblr_block_eigen)$method),
                   c("sbayesc", "sbayesr", "sbayesrc"))
})

test_that("scientific models are compositional and capability-checked", {
  sbayesc <- sblr:::.blr_resolve_st_model(
    "sbayesc", list(), "sbayesc", operator = "csr")
  sbayesr <- sblr:::.blr_resolve_st_model(
    "sbayesr", list(maf_effect_s = -.25), "sbayesr",
    operator = "csr")
  expect_identical(sbayesc$kernel, "bayesc")
  expect_null(sbayesc$dots$maf_effect_s)
  expect_identical(sbayesc$effect_scale, "unit")
  expect_identical(sbayesr$kernel, "bayesr")
  expect_identical(sbayesr$dots$maf_effect_s, -.25)
  expect_identical(sbayesr$effect_scale, "component_maf_s")
  expect_error(sblr:::.blr_resolve_st_model(
    "bayesc", list(maf_effect_s = 0), "sbayesc", operator = "csr"),
    "summary statistics")
})

test_that("unsupported scientific model/operator combinations fail early", {
  fixture <- blr_unified_fixture()
  on.exit(blr_unified_cleanup(fixture), add = TRUE)
  expect_error(stblr_bed(
    fixture$y, fixture$Glist, method = "sbayesc",
    nit = 1L, nburn = 0L, convergence = "none"),
    "s.*prefix denotes summary statistics")
  expect_error(mtblr_csr(method = "bayesrc"), "method")
})

test_that("unified fit finalizer owns explicit names without aliases", {
  fit <- sblr:::.blr_finalize_fit(list(
    input = list(), pi = 1, pim = 2, pis = matrix(3),
    covb = matrix(1), vb = matrix(2), bm_sd = matrix(0)),
    "mtblr", "bayesc", "csr")
  expect_s3_class(fit, "blr_fit")
  expect_true(all(c("family", "model", "operator", "input", "data",
                    "diagnostics", "convergence", "convergence_traces",
                    "chains", "memory_estimate", "pi_final", "pi_mean",
                    "pi_trace", "cov_b_mean", "cov_b_final",
                    "bm_chain_mean_sd") %in% names(fit)))
  expect_false(any(c("pi", "pim", "pis", "covb", "vb", "bm_sd") %in%
                     names(fit)))
  expect_true(all(c("genotype_scale", "effect_scale", "phenotype_scale",
                    "ld_scale", "n_total", "n_used", "n_by_trait") %in%
                    names(fit$data)))
})
