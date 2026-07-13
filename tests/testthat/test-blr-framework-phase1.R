phase1_resolved_spec <- function(chain_seeds = c(101L, 202L)) {
  sblr:::.blr_resolved_spec(
    n_markers = 3L,
    n_traits = 2L,
    marker_ids = c("m3", "m1", "m2"),
    trait_ids = c("trait_b", "trait_a"),
    sample_size = c(120L, 100L),
    resource_id = "fixture/phase1_ld",
    nit = 20L,
    nburn = 5L,
    nthin = 2L,
    nchains = 2L,
    ncores = 2L,
    seed = 17L,
    chain_seeds = chain_seeds,
    keep_chain_summaries = TRUE
  )
}

test_that("Phase 1 R specification is explicit and preserves order", {
  spec <- phase1_resolved_spec()

  expect_identical(spec$schema, list(name = "blr_resolved_spec", version = 1L))
  expect_identical(spec$data$representation, "csr")
  expect_identical(spec$data$design, "independent_traits")
  expect_identical(spec$data$n_markers, 3L)
  expect_identical(spec$data$n_traits, 2L)
  expect_identical(spec$data$marker_ids, c("m3", "m1", "m2"))
  expect_identical(spec$data$trait_ids, c("trait_b", "trait_a"))
  expect_identical(spec$data$sample_size, c(120L, 100L))
  expect_identical(spec$data$scaling, "standardized_genotype")
  expect_identical(spec$model, list(
    kernel = "scalar", family = "bayesc", state = "binary",
    probability = "global_binary", scale = "unit",
    trait_covariance = "scalar_independent",
    residual_covariance = "scalar_independent"
  ))
  expect_identical(spec$mcmc$nit, 20L)
  expect_identical(spec$mcmc$nburn, 5L)
  expect_identical(spec$mcmc$nthin, 2L)
  expect_identical(spec$mcmc$nchains, 2L)
  expect_identical(spec$mcmc$ncores, 2L)
  expect_identical(spec$mcmc$seed, 17L)
  expect_identical(spec$mcmc$chain_seeds, c(101L, 202L))
  expect_true(all(unlist(spec$output)))
  expect_identical(spec$execution$operator, "csr")
  expect_identical(spec$execution$backend_reference, "stblr_cpg_omp_csr")
  expect_false(spec$execution$scheduled)
})

test_that("Phase 1 R validation reports field-specific errors", {
  mutate_spec <- function(path, value) {
    spec <- phase1_resolved_spec()
    spec[[path[[1L]]]][[path[[2L]]]] <- value
    spec
  }

  expect_error(sblr:::.blr_resolved_spec(
    0L, 2L, character(), c("a", "b"), c(10L, 10L), "x", 2L, 0L
  ), "data\\$n_markers")
  expect_error(sblr:::.blr_resolved_spec(
    2L, 0L, c("a", "b"), character(), integer(), "x", 2L, 0L
  ), "data\\$n_traits")
  expect_error(sblr:::.validate_blr_resolved_spec(
    mutate_spec(c("mcmc", "nburn"), -1L)
  ), "mcmc\\$nburn")
  expect_error(sblr:::.validate_blr_resolved_spec(
    mutate_spec(c("mcmc", "nthin"), 0L)
  ), "mcmc\\$nthin")
  expect_error(sblr:::.validate_blr_resolved_spec(
    mutate_spec(c("mcmc", "nchains"), 0L)
  ), "mcmc\\$nchains")
  expect_error(sblr:::.validate_blr_resolved_spec(
    mutate_spec(c("mcmc", "ncores"), 0L)
  ), "mcmc\\$ncores")
  expect_error(sblr:::.validate_blr_resolved_spec(
    mutate_spec(c("mcmc", "chain_seeds"), 1L)
  ), "mcmc\\$chain_seeds")
  expect_error(sblr:::.validate_blr_resolved_spec(
    mutate_spec(c("execution", "operator"), "bed")
  ), "execution\\$operator")
  expect_error(sblr:::.validate_blr_resolved_spec(
    mutate_spec(c("model", "family"), "bayesr")
  ), "model\\$family")
  expect_error(sblr:::.validate_blr_resolved_spec(
    mutate_spec(c("execution", "scheduled"), TRUE)
  ), "scheduled BayesC")
  expect_error(sblr:::.validate_blr_resolved_spec(
    mutate_spec(c("model", "trait_covariance"), "full")
  ), "model\\$trait_covariance")
  expect_error(sblr:::.validate_blr_resolved_spec(
    mutate_spec(c("model", "residual_covariance"), "full")
  ), "model\\$residual_covariance")
  expect_error(sblr:::.validate_blr_resolved_spec(
    mutate_spec(c("data", "representation"), "CSR")
  ), "data\\$representation")
  expect_error(sblr:::.validate_blr_resolved_spec(
    mutate_spec(c("data", "marker_ids"), c("m1", "m1", "m2"))
  ), "data\\$marker_ids")
  expect_error(sblr:::.validate_blr_resolved_spec(
    mutate_spec(c("data", "trait_ids"), "trait_a")
  ), "data\\$trait_ids")
})

test_that("Phase 1 C++ round trip preserves every supported value", {
  explicit <- phase1_resolved_spec()
  absent <- phase1_resolved_spec(chain_seeds = NULL)

  expect_identical(sblr:::.blr_validate_spec_cpp(explicit), explicit)
  expect_identical(sblr:::.blr_validate_spec_cpp(absent), absent)
  expect_identical(
    sblr:::.blr_validate_spec_cpp(explicit)$data$marker_ids,
    c("m3", "m1", "m2")
  )
  expect_identical(
    sblr:::.blr_validate_spec_cpp(explicit)$data$trait_ids,
    c("trait_b", "trait_a")
  )
})

test_that("Phase 1 core headers are binding-neutral", {
  headers <- testthat::test_path(
    "..", "..", "src",
    c("blr_spec.h", "blr_result.h", "blr_csr_contract.h")
  )
  source <- paste(vapply(headers, function(path) {
    paste(readLines(path, warn = FALSE), collapse = "\n")
  }, character(1)), collapse = "\n")
  forbidden <- c(
    "Rcpp", "RcppArmadillo", "SEXP", "RObject", "NumericVector",
    "NumericMatrix", "List", "Nullable", "pybind11", "Python.h", "NumPy",
    "R::rnorm", "R::rchisq", "R::pnorm", "R::qnorm", "arma::randn",
    "arma::randu", "Rcpp::stop", "Rcpp::Rcout"
  )
  for (token in forbidden) {
    expect_false(grepl(token, source, fixed = TRUE), info = token)
  }
  expect_match(source, "std::invalid_argument", fixed = TRUE)
})

phase1_result_dimensions <- function() {
  list(
    markers = 3L,
    traits = 2L,
    retained_samples = 5L,
    parameter_dimension = 4L,
    marker_effect_length = 6L,
    marker_pip_length = 6L,
    final_effect_length = 6L,
    final_state_length = 6L,
    trace_length = 20L,
    trait_covariance_length = 4L,
    residual_covariance_length = 4L,
    has_component_probability = FALSE,
    components = 0L,
    has_pattern_probability = FALSE,
    patterns = 0L
  )
}

test_that("Phase 1 typed result vocabulary enforces canonical dimensions", {
  dimensions <- phase1_result_dimensions()
  normalized <- sblr:::.blr_validate_result_dimensions_cpp(dimensions)
  expect_identical(normalized$marker_effect, c(3L, 2L))
  expect_identical(normalized$marker_pip, c(3L, 2L))
  expect_identical(normalized$trace, c(5L, 4L))
  expect_identical(normalized$trait_covariance, c(2L, 2L))
  expect_identical(normalized$residual_covariance, c(2L, 2L))
  expect_null(normalized$component_probability)
  expect_null(normalized$pattern_probability)

  component <- dimensions
  component$has_component_probability <- TRUE
  component$components <- 4L
  expect_identical(
    sblr:::.blr_validate_result_dimensions_cpp(component)$component_probability,
    c(3L, 4L, 2L)
  )
  pattern <- dimensions
  pattern$has_pattern_probability <- TRUE
  pattern$patterns <- 3L
  expect_identical(
    sblr:::.blr_validate_result_dimensions_cpp(pattern)$pattern_probability,
    c(3L, 3L)
  )

  bad <- dimensions
  bad$marker_effect_length <- 5L
  expect_error(
    sblr:::.blr_validate_result_dimensions_cpp(bad),
    "marker_effect_length"
  )
  bad <- dimensions
  bad$trace_length <- 19L
  expect_error(sblr:::.blr_validate_result_dimensions_cpp(bad), "trace_length")
  bad <- dimensions
  bad$trait_covariance_length <- 3L
  expect_error(
    sblr:::.blr_validate_result_dimensions_cpp(bad),
    "trait_covariance_length"
  )
  bad <- dimensions
  bad$has_component_probability <- TRUE
  expect_error(sblr:::.blr_validate_result_dimensions_cpp(bad), "component count")
})

test_that("Phase 1 CSR ownership is shared, read-only, and payload-free", {
  csr <- phase1_resolved_spec()$data$csr
  expect_true(csr$shared_read_only)
  expect_false(csr$per_chain_data)
  expect_true(csr$lifetime_exceeds_chains)
  expect_identical(csr$marker_count, 3L)
  expect_false(any(c("values", "indices", "row_ptr") %in% names(csr)))

  bad <- phase1_resolved_spec()
  bad$data$csr$resource_id <- ""
  expect_error(sblr:::.validate_blr_resolved_spec(bad), "resource_id")
  bad <- phase1_resolved_spec()
  bad$data$csr$marker_count <- 2L
  expect_error(sblr:::.validate_blr_resolved_spec(bad), "marker_count")
  bad <- phase1_resolved_spec()
  bad$data$csr$shared_read_only <- FALSE
  expect_error(sblr:::.validate_blr_resolved_spec(bad), "shared_read_only")
})

phase1_make_csr_prefix <- function() {
  prefix <- tempfile("blr_phase1_csr_")
  sblr:::.stblr_write_uint64_file(
    paste0(prefix, ".row_ptr.u64.bin"), c(0, 1, 1, 1)
  )
  sblr:::.stblr_write_uint32_file(
    paste0(prefix, ".col_idx.u32.0based.bin"), 1L
  )
  writeBin(0.4, paste0(prefix, ".values.f32.bin"),
           size = 4, endian = "little")
  writeLines(c(
    "format=sparse_ld_csr", "storage=streamed_upper_triangle",
    "n_bed=NA", "n_used=NA", "n_samples_used=NA", "n_variants=3",
    "nnz=1", "triangle=upper", "diagonal=implicit_1",
    paste0("row_ptr_file=", prefix, ".row_ptr.u64.bin"),
    paste0("col_idx_file=", prefix, ".col_idx.u32.0based.bin"),
    paste0("values_file=", prefix, ".values.f32.bin"),
    "row_ptr_type=uint64", "col_idx_type=uint32", "values_type=float32",
    "index_base=0", "value=r"
  ), paste0(prefix, ".meta.txt"))
  prefix
}

phase1_csr_stats <- function() {
  markers <- c("m1", "m2", "m3")
  list(
    wy = list(trait1 = stats::setNames(c(25, 15, 4), markers)),
    ww = list(trait1 = stats::setNames(rep(100, 3), markers)),
    yy = stats::setNames(100, "trait1"),
    n = 100L,
    m = 3L,
    marker_names = markers,
    trait_names = "trait1"
  )
}

phase1_csr_glist <- function(prefix) {
  list(
    rsidsLD = list(c("m1", "m2", "m3")),
    rsids = list(c("m3", "m1", "m2")),
    maf = list(c(0.40, 0.05, 0.20)),
    sparseLD = list(chr = 1L, cls = list(1:3), prefix = prefix)
  )
}

phase1_fit_csr <- function(ncores = 1L, selection_s = NULL,
                           chain_seeds = c(401L, 402L)) {
  prefix <- phase1_make_csr_prefix()
  stblr_csr(
    stats = phase1_csr_stats(),
    Glist = phase1_csr_glist(prefix),
    ld_prefix = prefix,
    pi_init = 0.5,
    pi_prior_mean = 0.5,
    pi_prior_strength = 2,
    updateB = FALSE,
    updateE = FALSE,
    updatePi = TRUE,
    nit = 8L,
    nburn = 2L,
    nthin = 1L,
    seed = 31L,
    nchains = 2L,
    keep_chains = TRUE,
    chain_seeds = chain_seeds,
    ncores = ncores,
    selection_s = selection_s,
    updateLDswap = FALSE,
    scheduled = FALSE
  )
}

phase1_stable_fit_fields <- function(fit) {
  fit[c(
    "bm", "dm", "wy", "r", "b", "d", "vbs", "vgs", "ves", "vle",
    "vld", "pis", "pi", "pim", "covb", "covg", "cove", "dm_sd",
    "dm_min", "dm_max", "bm_sd", "bm_min", "bm_max", "chains",
    "ld_swap", "ld_swap_chains"
  )]
}

test_that("existing unscheduled CSR BayesC remains the deterministic reference", {
  one <- phase1_fit_csr(ncores = 1L)
  repeated <- phase1_fit_csr(ncores = 1L)
  two_core <- phase1_fit_csr(ncores = 2L)
  fixed_unit <- phase1_fit_csr(ncores = 1L, selection_s = -1)

  expect_identical(phase1_stable_fit_fields(one),
                   phase1_stable_fit_fields(repeated))
  expect_identical(phase1_stable_fit_fields(one),
                   phase1_stable_fit_fields(two_core))
  expect_identical(phase1_stable_fit_fields(one),
                   phase1_stable_fit_fields(fixed_unit))
  expect_identical(rownames(one$bm), c("m1", "m2", "m3"))
  expect_identical(colnames(one$bm), "trait1")
  expect_identical(dim(one$bm), c(3L, 1L))
  expect_identical(dim(one$dm), c(3L, 1L))
  expect_equal(one$dm[, 1L], rowMeans(vapply(
    one$chains$trait1, `[[`, numeric(3), "dm"
  )), tolerance = 1e-12)
  expect_null(one$ld_swap)
  expect_null(one$ld_swap_chains)
  expect_identical(one$input$backend, "csr_bayesc")
  expect_false(one$input$scheduled)
  expect_false(any(c("blr_spec", "resolved_spec", "phase1") %in% names(one)))
})

test_that("explicit CSR chain seeds retain their current mapping", {
  first <- phase1_fit_csr(chain_seeds = c(501L, 502L))
  repeated <- phase1_fit_csr(chain_seeds = c(501L, 502L))
  changed <- phase1_fit_csr(chain_seeds = c(501L, 503L))

  expect_identical(phase1_stable_fit_fields(first),
                   phase1_stable_fit_fields(repeated))
  expect_false(identical(first$chains$trait1[[2L]]$bm,
                         changed$chains$trait1[[2L]]$bm))
})

test_that("the existing backend still emits stblr_raw_v1", {
  stats <- phase1_csr_stats()
  prefix <- phase1_make_csr_prefix()
  vy <- as.numeric(stats$yy) / (stats$n - 1)
  B <- diag((vy * 0.3) / (stats$m * 0.5), 1L, 1L)
  E <- diag(vy * 0.7, 1L, 1L)
  raw <- stblr_cpg_omp_csr(
    wy = stats$wy, ww = stats$ww, yy = stats$yy,
    b_init = list(rep(0, stats$m)), d_init = list(rep(0, stats$m)),
    use_d_init = FALSE, r_init = stats$wy, use_r_init = FALSE,
    rebuild_r_before_updateE = FALSE, ld_prefix = prefix,
    B = B, E = E,
    ssb_prior = list(((4 - 2) / 4) * B),
    sse_prior = list(((4 - 2) / 4) * E),
    pi = c(0.5, 0.5), nub = 4, nue = 4,
    updateB = FALSE, updateE = FALSE, updatePi = TRUE, adjE = 0.9,
    n = stats$n, nit = 4L, nburn = 1L, nthin = 1L,
    pi_prior_a = 1, pi_prior_b = 1, ncores = 1L, seed = 41L,
    nchains = 1L, keep_chains = FALSE, chain_seeds = integer(),
    updateLDswap = FALSE
  )

  expect_s3_class(raw, "stblr_raw_v1")
  expect_identical(raw$schema$class, "stblr_raw")
  expect_identical(raw$schema$version, 1L)
  expect_identical(dim(raw$marker$bm), c(3L, 1L))
  expect_identical(dim(raw$marker$dm), c(3L, 1L))
  expect_identical(raw$meta$backend, "csr_bayesc")
})
