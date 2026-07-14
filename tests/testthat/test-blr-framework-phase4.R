source(test_path("fixtures", "blr-phase2-reference.R"), local = TRUE)
source(test_path("fixtures", "blr-phase2-reference-hashes.R"), local = TRUE)

phase4_quiet <- function(expression) {
  invisible(capture.output(value <- suppressWarnings(force(expression))))
  value
}

phase4_native <- function(config) {
  prefix <- phase2_reference_write_csr()
  inputs <- phase2_reference_inputs(config$traits)
  phase4_quiet(phase2_reference_native(config, prefix, inputs))
}

phase4_stable <- function(value) phase2_reference_normalize(value)

test_that("scalar task construction preserves trait-major chain order", {
  expect_identical(
    blr_phase4_scalar_tasks_cpp(2L, 3L),
    data.frame(
      trait = c(0L, 0L, 0L, 1L, 1L, 1L),
      chain = c(0L, 1L, 2L, 0L, 1L, 2L),
      task = 0:5
    )
  )
  expect_identical(
    blr_phase4_scalar_tasks_cpp(1L, 1L),
    data.frame(trait = 0L, chain = 0L, task = 0L)
  )
  expect_error(blr_phase4_scalar_tasks_cpp(0L, 1L), "trait_count")
  expect_error(blr_phase4_scalar_tasks_cpp(1L, 0L), "chain_count")
  expect_error(blr_phase4_scalar_tasks_cpp(-1L, 1L), "non-negative")
})

test_that("scalar seed resolution exactly preserves current mappings", {
  expect_identical(
    blr_phase4_scalar_seeds_cpp(31L, 1L, 1L, integer()),
    1009210L
  )
  expect_identical(
    blr_phase4_scalar_seeds_cpp(31L, 1L, 2L, integer()),
    c(1009210L, 1018386L)
  )
  expect_identical(
    blr_phase4_scalar_seeds_cpp(31L, 2L, 2L, integer()),
    c(1009210L, 1018386L, 2009213L, 2018389L)
  )
  expect_identical(
    blr_phase4_scalar_seeds_cpp(31L, 2L, 2L, c(401L, 402L)),
    c(1000404L, 1000405L, 2000407L, 2000408L)
  )
  expect_identical(
    blr_phase4_scalar_seeds_cpp(1000000000L, 2L, 1L, integer()),
    c(1001009179L, 1002009182L)
  )
  expect_identical(
    blr_phase4_scalar_seeds_cpp(31L, 2L, 2L, integer()),
    blr_phase4_scalar_seeds_cpp(31L, 2L, 2L, integer())
  )
  expect_error(
    blr_phase4_scalar_seeds_cpp(31L, 1L, 2L, 401L),
    "must match chain_count"
  )
})

test_that("retained-iteration counting is the shared exact predicate", {
  expect_identical(
    blr_phase4_retained_iterations_cpp(10L, 2L, 3L),
    c(FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE)
  )
  expect_identical(
    which(blr_phase4_retained_iterations_cpp(10L, 2L, 3L)),
    c(3L, 6L, 9L)
  )
  expect_error(blr_phase4_retained_iterations_cpp(10L, 2L, 0L), "invalid")
})

test_that("shared scalar infrastructure is binding-neutral and allocation-safe", {
  header <- paste(readLines(
    test_path("..", "..", "src", "blr_scalar_execution.h"), warn = FALSE
  ), collapse = "\n")
  forbidden <- c(
    "Rcpp", "RcppArmadillo", "SEXP", "RObject", "NumericVector",
    "NumericMatrix", "List", "Nullable", "Rcpp::stop", "Rcpp::Rcout",
    "R::rnorm", "R::rchisq", "R::pnorm", "R::qnorm", "arma::randn",
    "arma::randu", "pybind11", "Python.h"
  )
  for (token in forbidden) {
    expect_false(grepl(token, header, fixed = TRUE), info = token)
  }
  expect_match(header, "tasks.reserve(task_count)", fixed = TRUE)
  expect_match(header, "RNG engine and every stateful normal", fixed = TRUE)
  expect_false(grepl("virtual ", header, fixed = TRUE))
  expect_false(grepl("std::function", header, fixed = TRUE))
})

test_that("comparative boundaries remain shared only where proven", {
  bayesc <- paste(readLines(
    test_path("..", "..", "src", "blr_csr_bayesc_core_impl.h"),
    warn = FALSE
  ), collapse = "\n")
  bayesr_path <- test_path("..", "..", "src", "st_cpg_omp_csr_bayesr.cpp")
  bayesr_core_path <- test_path("..", "..", "src", "blr_csr_bayesr_core_impl.h")
  bayesr <- paste(c(readLines(bayesr_path, warn = FALSE),
    readLines(bayesr_core_path, warn = FALSE)), collapse = "\n")

  expect_match(bayesc, "make_scalar_chain_tasks", fixed = TRUE)
  expect_match(bayesc, "resolve_scalar_chain_seed", fixed = TRUE)
  expect_match(bayesc, "scalar_iteration_is_retained", fixed = TRUE)
  expect_match(bayesr, "stblr_num_chain_tasks(nt, nchains)", fixed = TRUE)
  expect_match(bayesr, "stblr_task_trait(task, nchains)", fixed = TRUE)
  expect_match(bayesr, "stblr_seed_with_chain_base", fixed = TRUE)
  expect_false(grepl("blr_scalar_execution", bayesr, fixed = TRUE))
  expect_match(bayesc, "delta_log", fixed = TRUE)
  expect_match(bayesr, "sample_categorical_logprob_bayesr", fixed = TRUE)
  expect_equal(length(gregexpr("for (int it = 0; it < trace_len", bayesr, fixed = TRUE)[[1L]]), 1L)
  expect_match(bayesr, "run_csr_bayesr", fixed = TRUE)
})

test_that("all frozen BayesC raw and formatted references remain exact", {
  observed <- phase4_quiet(phase2_reference_hashes())
  expect_identical(names(observed), names(phase2_reference_expected_hashes))
  for (name in names(observed)) {
    expect_identical(
      observed[[name]], phase2_reference_expected_hashes[[name]], info = name
    )
  }
})

test_that("BayesC reproducibility remains input and seed only", {
  one <- phase2_reference_configurations$one_trait_two_chains_one_core
  two <- phase2_reference_configurations$one_trait_two_chains_two_cores
  one_first <- phase4_stable(phase4_native(one))
  two_first <- phase4_stable(phase4_native(two))
  two_second <- phase4_stable(phase4_native(two))
  one_second <- phase4_stable(phase4_native(one))

  expect_identical(one_first, one_second)
  expect_identical(two_first, two_second)
  expect_identical(one_first, two_first)
  phase4_native(phase2_reference_configurations$multiple_traits)
  expect_identical(one_first, phase4_stable(phase4_native(one)))

  explicit <- phase4_native(
    phase2_reference_configurations$explicit_chain_seeds
  )
  kept <- phase4_native(phase2_reference_configurations$keep_chains)
  expect_identical(phase4_stable(explicit$marker), phase4_stable(kept$marker))
})

test_that("ownership, routes, and public schemas remain unchanged", {
  types <- paste(readLines(
    test_path("..", "..", "src", "blr_csr_bayesc_types.h"), warn = FALSE
  ), collapse = "\n")
  public_r <- paste(readLines(
    test_path("..", "..", "R", "sparse_ld_bed_helper.R"), warn = FALSE
  ), collapse = "\n")
  expect_match(types, "Borrowed immutable view", fixed = TRUE)
  expect_match(types, "no chain result or chain state owns a CSR payload",
               fixed = TRUE)
  expect_match(public_r, "do.call(stblr_cpg_omp_csr", fixed = TRUE)
  expect_match(public_r, "stblr_cpg_omp_csr_bayesr", fixed = TRUE)

  raw <- phase4_native(
    phase2_reference_configurations$one_trait_one_chain_one_core
  )
  expect_s3_class(raw, "stblr_raw_v1")
  expect_identical(raw$schema, list(class = "stblr_raw", version = 1L))
  expect_true(is.matrix(raw$marker$bm))
  expect_true(is.matrix(raw$marker$dm))
  expect_null(raw$diagnostics$ld_swap)
})
