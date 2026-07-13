source(test_path("fixtures", "blr-phase2-reference.R"), local = TRUE)
source(test_path("fixtures", "blr-phase2-reference-hashes.R"), local = TRUE)

phase2_quiet <- function(expression) {
  suppressWarnings(force(expression))
}

phase2_native_object <- function(config) {
  prefix <- phase2_reference_write_csr()
  inputs <- phase2_reference_inputs(config$traits)
  phase2_quiet(phase2_reference_native(config, prefix, inputs))
}

phase2_stable <- function(x) phase2_reference_normalize(x)

test_that("Phase 2 exactly matches every frozen pre-refactor reference", {
  observed <- phase2_quiet(phase2_reference_hashes())
  expect_identical(names(observed), names(phase2_reference_expected_hashes))
  for (name in names(observed)) {
    expect_identical(observed[[name]], phase2_reference_expected_hashes[[name]],
                     info = name)
  }
})

test_that("fixed seeds are independent of cores and prior call order", {
  one <- phase2_reference_configurations$one_trait_two_chains_one_core
  two <- phase2_reference_configurations$one_trait_two_chains_two_cores

  one_first <- phase2_stable(phase2_native_object(one))
  two_first <- phase2_stable(phase2_native_object(two))
  two_second <- phase2_stable(phase2_native_object(two))
  one_second <- phase2_stable(phase2_native_object(one))

  expect_identical(one_first, one_second)
  expect_identical(two_first, two_second)
  expect_identical(one_first, two_first)

  phase2_native_object(phase2_reference_configurations$multiple_traits)
  expect_identical(one_first, phase2_stable(phase2_native_object(one)))
})

test_that("explicit chain seeds and optional chain retention preserve aggregates", {
  no_chains <- phase2_native_object(
    phase2_reference_configurations$explicit_chain_seeds
  )
  with_chains <- phase2_native_object(
    phase2_reference_configurations$keep_chains
  )

  expect_identical(no_chains$chains, list())
  expect_named(with_chains$chains, "trait1")
  expect_named(with_chains$chains$trait1, c("chain1", "chain2"))
  expect_named(with_chains$chains$trait1$chain1,
               c("marker", "trace", "pi", "selection", "diagnostics"))
  expect_named(with_chains$chains$trait1$chain2,
               c("marker", "trace", "pi", "selection", "diagnostics"))
  with_chains["chains"] <- list(list())
  with_chains$meta$keep_chains <- FALSE
  expect_identical(phase2_stable(no_chains), phase2_stable(with_chains))
})

test_that("multiple traits, selection_s, and disabled LD swap retain schemas", {
  multi <- phase2_native_object(phase2_reference_configurations$multiple_traits)
  fixed <- phase2_native_object(phase2_reference_configurations$fixed_selection_s)

  expect_s3_class(multi, "stblr_raw_v1")
  expect_identical(multi$schema$class, "stblr_raw")
  expect_identical(multi$schema$version, 1L)
  expect_identical(dim(multi$marker$bm), c(3L, 2L))
  expect_identical(dim(multi$marker$dm), c(3L, 2L))
  expect_identical(dim(multi$trace$vbs), c(10L, 2L))
  expect_null(multi$diagnostics$ld_swap)
  expect_true(fixed$selection$enabled)
  expect_true(fixed$selection$fixed)
  expect_null(fixed$selection$trace)
  expect_identical(dim(fixed$marker$bm), c(3L, 1L))
})

test_that("the migrated reusable files remain binding-neutral", {
  files <- c("blr_csr_bayesc_types.h", "blr_csr_bayesc_core.h")
  text <- paste(vapply(files, function(file) {
    paste(readLines(test_path("..", "..", "src", file), warn = FALSE),
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
})

test_that("CSR ownership and marker-loop allocation contracts are explicit", {
  types <- paste(readLines(
    test_path("..", "..", "src", "blr_csr_bayesc_types.h"), warn = FALSE
  ), collapse = "\n")
  core <- readLines(
    test_path("..", "..", "src", "blr_csr_bayesc_core.h"), warn = FALSE
  )

  expect_match(types, "Borrowed immutable view", fixed = TRUE)
  expect_match(types, "const std::uint64_t* row_ptr", fixed = TRUE)
  expect_match(types, "const int* column_index", fixed = TRUE)
  expect_match(types, "const float* values", fixed = TRUE)
  expect_match(types, "const arma::mat* wy", fixed = TRUE)
  expect_false(grepl("std::vector<float> values", types, fixed = TRUE))

  begin <- grep("PHASE2_MARKER_LOOP_BEGIN", core, fixed = TRUE)
  end <- grep("PHASE2_MARKER_LOOP_END", core, fixed = TRUE)
  expect_length(begin, 2L)
  expect_length(end, 2L)
  loops <- paste(unlist(Map(function(first, last) {
    core[seq.int(first, last)]
  }, begin, end), use.names = FALSE), collapse = "\n")
  allocation_tokens <- c(
    "new ", "malloc(", "calloc(", "realloc(", ".set_size(",
    ".resize(", ".push_back(", "std::vector<"
  )
  for (token in allocation_tokens) {
    expect_false(grepl(token, loops, fixed = TRUE), info = token)
  }
})
