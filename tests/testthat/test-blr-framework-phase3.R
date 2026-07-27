source(test_path("fixtures", "blr-phase2-reference.R"), local = TRUE)
source(test_path("fixtures", "blr-phase2-reference-hashes.R"), local = TRUE)

phase3_quiet <- function(expression) {
  invisible(capture.output(value <- suppressWarnings(force(expression))))
  value
}

phase3_native <- function(config) {
  prefix <- phase2_reference_write_csr()
  inputs <- phase2_reference_inputs(config$traits)
  phase3_quiet(phase2_reference_native(config, prefix, inputs))
}

phase3_fit <- function(config) {
  prefix <- phase2_reference_write_csr()
  inputs <- phase2_reference_inputs(config$traits)
  phase3_quiet(phase2_reference_formatted(config, prefix, inputs))
}

phase3_stable <- function(value) phase2_reference_normalize(value)

phase3_one_marker_raw <- function() {
  prefix <- tempfile("blr_phase3_one_marker_")
  sblr:::.stblr_write_uint64_file(
    paste0(prefix, ".row_ptr.u64.bin"), c(0, 0)
  )
  sblr:::.stblr_write_uint32_file(
    paste0(prefix, ".col_idx.u32.0based.bin"), integer()
  )
  writeBin(
    numeric(), paste0(prefix, ".values.f32.bin"),
    size = 4, endian = "little"
  )
  writeLines(c(
    "format=sparse_ld_csr", "storage=streamed_upper_triangle",
    "n_bed=NA", "n_used=NA", "n_samples_used=NA", "n_variants=1",
    "nnz=0", "triangle=upper", "diagonal=implicit_1",
    paste0("row_ptr_file=", prefix, ".row_ptr.u64.bin"),
    paste0("col_idx_file=", prefix, ".col_idx.u32.0based.bin"),
    paste0("values_file=", prefix, ".values.f32.bin"),
    "row_ptr_type=uint64", "col_idx_type=uint32", "values_type=float32",
    "index_base=0", "value=r"
  ), paste0(prefix, ".meta.txt"))

  phase3_quiet(stblr_cpg_omp_csr(
    wy = list(5), ww = list(100), yy = 100,
    b_init = list(0), d_init = list(0), use_d_init = FALSE,
    r_init = list(5), use_r_init = FALSE,
    rebuild_r_before_updateE = FALSE, ld_prefix = prefix,
    B = matrix(0.1, 1, 1), E = matrix(0.7, 1, 1),
    ssb_prior = list(0.05), sse_prior = list(0.35),
    pi = c(0.5, 0.5), nub = 4, nue = 4,
    updateB = FALSE, updateE = FALSE, updatePi = TRUE, adjE = 0.9,
    n = 100L, nit = 3L, nburn = 1L, nthin = 1L,
    pi_prior_a = 1, pi_prior_b = 1,
    ncores = 1L, seed = 31L, nchains = 1L, keep_chains = FALSE,
    chain_seeds = integer(), updateLDswap = FALSE
  ))
}

test_that("ordinary CSR BayesC has one canonical typed production route", {
  binding <- paste(readLines(
    blr_repo_path("src", "st_cpg_omp_csr.cpp"), warn = FALSE
  ), collapse = "\n")
  core <- paste(readLines(
    blr_repo_path("src", "blr_csr_bayesc_core_impl.h"),
    warn = FALSE
  ), collapse = "\n")
  public_r <- paste(readLines(
    blr_repo_path("R", "sparse_ld_bed_helper.R"), warn = FALSE
  ), collapse = "\n")

  expect_source_count("CsrBayesCResult run_csr_bayesc\\(", core, 1L, fixed = FALSE)
  expect_source_count("inline void sample_marker_scaled\\(", core, 1L, fixed = FALSE)
  expect_source_count("inline void sample_marker_unscaled\\(", core, 1L, fixed = FALSE)
  expect_source_count("run_csr_bayesc\\(input\\)", binding, 1L, fixed = FALSE)
  expect_match(
    binding,
    paste0(
      "if constexpr[\\s\\S]{0,300}CsrOperator[\\s\\S]{0,300}",
      "return stblr_csr_bayesc_run_canonical"
    ),
    perl = TRUE
  )
  expect_match(public_r, "do.call(stblr_cpg_omp_csr", fixed = TRUE)

  selectors <- c(
    "std::getenv(", "getenv(", "SBLR_CSR_BAYESC_USE_LEGACY",
    "legacy_csr_bayesc", "use_legacy_bayesc"
  )
  for (selector in selectors) {
    expect_false(grepl(selector, binding, fixed = TRUE), info = selector)
  }
})

test_that("the native CSR adapter and raw converter have single responsibilities", {
  source <- readLines(
    blr_repo_path("src", "st_cpg_omp_csr.cpp"), warn = FALSE
  )
  text <- paste(source, collapse = "\n")
  expect_length(
    grep("static Rcpp::List stblr_csr_bayesc_result_to_raw", source, fixed = TRUE),
    1L
  )
  expect_length(
    grep("static Rcpp::List stblr_csr_bayesc_run_canonical", source, fixed = TRUE),
    1L
  )

  first <- grep("static Rcpp::List stblr_csr_bayesc_run_canonical", source,
                fixed = TRUE)
  last <- grep("Main implementation: parallel single-trait BayesC", source,
               fixed = TRUE) - 2L
  adapter <- paste(source[seq.int(first, last)], collapse = "\n")
  expect_match(adapter, "CsrBayesCExecutionInput input", fixed = TRUE)
  expect_match(adapter, "run_csr_bayesc(input)", fixed = TRUE)
  expect_match(adapter, "stblr_csr_bayesc_result_to_raw(result, conversion)",
               fixed = TRUE)
  expect_false(grepl("sample_marker_", adapter, fixed = TRUE))
  expect_false(grepl("for (int it =", adapter, fixed = TRUE))
  expect_false(grepl("std::mt19937", adapter, fixed = TRUE))
  expect_source_count("stblr_csr_bayesc_result_to_raw\\(", text, 2L, fixed = FALSE)
})

test_that("the binding-neutral implementation header is single-TU and guarded", {
  src_dir <- blr_repo_path("src")
  files <- list.files(src_dir, pattern = "\\.(cpp|h)$", full.names = TRUE)
  include_pattern <- '#include "blr_csr_bayesc_core_impl.h"'
  includes <- vapply(files, function(file) {
    any(grepl(include_pattern, readLines(file, warn = FALSE), fixed = TRUE))
  }, logical(1))
  expect_identical(basename(files[includes]), "st_cpg_omp_csr.cpp")

  binding <- paste(readLines(
    file.path(src_dir, "st_cpg_omp_csr.cpp"), warn = FALSE
  ), collapse = "\n")
  core <- paste(readLines(
    file.path(src_dir, "blr_csr_bayesc_core_impl.h"), warn = FALSE
  ), collapse = "\n")
  expect_match(binding, "#define SBLR_CSR_BAYESC_CORE_IMPL_TRANSLATION_UNIT 1",
               fixed = TRUE)
  expect_match(binding, "#undef SBLR_CSR_BAYESC_CORE_IMPL_TRANSLATION_UNIT",
               fixed = TRUE)
  expect_match(core, "#ifndef SBLR_CSR_BAYESC_CORE_IMPL_TRANSLATION_UNIT",
               fixed = TRUE)
  expect_match(core, "#error", fixed = TRUE)
})

test_that("canonical contracts retain binding neutrality and borrowed ownership", {
  files <- c(
    "blr_spec.h", "blr_result.h", "blr_csr_contract.h",
    "blr_csr_bayesc_types.h", "blr_csr_bayesc_core_impl.h"
  )
  text <- paste(vapply(files, function(file) {
    paste(readLines(blr_repo_path("src", file), warn = FALSE),
          collapse = "\n")
  }, character(1)), collapse = "\n")
  forbidden <- c(
    "Rcpp", "RcppArmadillo", "SEXP", "RObject", "NumericVector",
    "NumericMatrix", "List", "Nullable", "Rcpp::stop", "Rcpp::Rcout",
    "R::rnorm", "R::rchisq", "R::pnorm", "R::qnorm", "arma::randn",
    "arma::randu", "pybind11", "Python.h"
  )
  for (token in forbidden) {
    expect_false(grepl(token, text, fixed = TRUE), info = token)
  }

  types <- readLines(
    blr_repo_path("src", "blr_csr_bayesc_types.h"), warn = FALSE
  )
  types_text <- paste(types, collapse = "\n")
  expect_match(types_text, "Borrowed immutable view", fixed = TRUE)
  expect_match(types_text, "must keep it alive", fixed = TRUE)
  shared_text <- paste(readLines(
    blr_repo_path("src", "blr_sparse_ld_csr.h"), warn = FALSE
  ), collapse = "\n")
  expect_match(shared_text, "const std::uint64_t* row_ptr", fixed = TRUE)
  expect_match(shared_text, "const int* column_index", fixed = TRUE)
  expect_match(shared_text, "const float* offdiag_xij", fixed = TRUE)

  first <- grep("struct CsrBayesCChainResult", types, fixed = TRUE)
  last <- grep("struct CsrBayesCResult", types, fixed = TRUE) - 1L
  chain_payload <- paste(types[seq.int(first, last)], collapse = "\n")
  expect_false(grepl("row_ptr =", chain_payload, fixed = TRUE))
  expect_false(grepl("column_index =", chain_payload, fixed = TRUE))
  expect_false(grepl("const float* offdiag_xij", chain_payload, fixed = TRUE))
})

test_that("all permanent pre-refactor raw numerical references remain exact", {
  blr_skip_if_no_source_tree()
  observed <- phase3_quiet(phase2_reference_hashes())
  expect_identical(names(observed), names(phase2_reference_expected_hashes))
  for (name in names(observed)) {
    if (!identical(name, "keep_chains"))
      expect_identical(observed[[name]][["raw"]],
        phase2_reference_expected_hashes[[name]][["raw"]], info = name)
  }
})

test_that("repeated, reversed-core, and intervening calls are exact", {
  one_config <- phase2_reference_configurations$one_trait_two_chains_one_core
  two_config <- phase2_reference_configurations$one_trait_two_chains_two_cores
  one_first <- phase3_stable(phase3_native(one_config))
  two_first <- phase3_stable(phase3_native(two_config))
  two_second <- phase3_stable(phase3_native(two_config))
  one_second <- phase3_stable(phase3_native(one_config))

  expect_identical(one_first, one_second)
  expect_identical(two_first, two_second)
  expect_identical(one_first, two_first)
  phase3_native(phase2_reference_configurations$multiple_traits)
  expect_identical(one_first, phase3_stable(phase3_native(one_config)))
})

test_that("explicit seeds, chain retention, and aggregation remain exact", {
  kept <- phase3_native(phase2_reference_configurations$keep_chains)
  dropped <- phase3_native(phase2_reference_configurations$explicit_chain_seeds)

  expect_named(kept$chains, "trait1")
  expect_named(kept$chains$trait1, c("chain1", "chain2"))
  chain_bm <- vapply(kept$chains$trait1, function(chain) chain$marker$bm,
                     numeric(3))
  chain_dm <- vapply(kept$chains$trait1, function(chain) chain$marker$dm,
                     numeric(3))
  expect_identical(as.numeric(rowMeans(chain_bm)), as.numeric(kept$marker$bm))
  expect_identical(as.numeric(rowMeans(chain_dm)), as.numeric(kept$marker$dm))

  kept_without_chains <- kept
  kept_without_chains["chains"] <- list(list())
  kept_without_chains$meta$keep_chains <- FALSE
  expect_identical(
    phase3_stable(kept_without_chains), phase3_stable(dropped)
  )
})

test_that("multiple traits, selection_s, LD-swap-off, and schemas are stable", {
  multi <- phase3_fit(phase2_reference_configurations$multiple_traits)
  fixed <- phase3_fit(phase2_reference_configurations$fixed_selection_s)
  ordinary <- phase3_fit(
    phase2_reference_configurations$one_trait_one_chain_one_core
  )

  expect_identical(colnames(multi$bm), c("trait1", "trait2"))
  expect_identical(rownames(multi$bm), c("m1", "m2", "m3"))
  expect_identical(dim(multi$bm), c(3L, 2L))
  expect_identical(dim(multi$dm), c(3L, 2L))
  expect_true(is.matrix(fixed$bm))
  expect_identical(dim(fixed$bm), c(3L, 1L))
  expect_true(isTRUE(fixed$input$selection_s_fixed))
  expect_null(ordinary$diagnostics$ld_swap)
  expect_null(ordinary$diagnostics$ld_swap_chains)
  expect_s3_class(ordinary, "stblr_fit")
  expect_s3_class(ordinary, "blr_fit")
  expect_false(any(grepl("blr_resolved_spec|typed", names(ordinary))))

  raw <- phase3_native(
    phase2_reference_configurations$one_trait_one_chain_one_core
  )
  expect_s3_class(raw, "stblr_raw_v1")
  expect_identical(raw$schema, list(class = "stblr_raw", version = 1L))
  expect_true(is.matrix(raw$marker$bm))
  expect_identical(dim(raw$marker$bm), c(3L, 1L))
  expect_null(raw$diagnostics$ld_swap)
  expect_null(raw$selection$trace)
})

test_that("one-marker and one-trait conversion preserves matrices and NULL", {
  raw <- phase3_one_marker_raw()
  fit <- sblr:::.as_stblr_fit(
    raw, trait_names = "trait1", variable_names = "m1"
  )

  for (field in c("bm", "dm", "wy", "r", "b", "state")) {
    expect_true(is.matrix(raw$marker[[field]]), info = field)
    expect_identical(dim(raw$marker[[field]]), c(1L, 1L), info = field)
  }
  expect_true(is.matrix(fit$bm))
  expect_true(is.matrix(fit$dm))
  expect_identical(dim(fit$bm), c(1L, 1L))
  expect_identical(dim(fit$dm), c(1L, 1L))
  expect_null(raw$diagnostics$ld_swap)
  expect_null(fit$ld_swap)
  expect_null(fit$ld_swap_chains)
})
