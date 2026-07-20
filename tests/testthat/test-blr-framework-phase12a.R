test_that("source architecture helpers distinguish absent and duplicate matches", {
  expect_identical(source_match_positions("z", "abc"), integer())
  expect_identical(source_match_count("z", "abc"), 0L)
  expect_identical(source_match_count("b", "abc"), 1L)
  expect_identical(source_match_count("a", "ababa"), 3L)
  expect_identical(source_match_count("[0-9]+", "ab12cd34", fixed = FALSE), 2L)
  expect_failure(expect_source_count("missing", "abc", 1L))
})

test_that("maintained developer Markdown has no unintended control bytes", {
  files <- list.files(test_path("..", "..", "docs", "dev"), "[.]md$",
                      full.names = TRUE)
  bad <- character()
  for (file in files) {
    bytes <- readBin(file, "raw", n = file.info(file)$size)
    offsets <- which(as.integer(bytes) < 32L & !as.integer(bytes) %in% c(9L, 10L, 13L))
    if (length(offsets)) {
      for (offset in offsets) {
        prior <- bytes[seq_len(max(0L, offset - 1L))]
        line <- 1L + sum(prior == as.raw(10L))
        last_lf <- tail(which(prior == as.raw(10L)), 1L)
        column <- if (length(last_lf)) offset - last_lf else offset
        bad <- c(bad, sprintf("%s: byte %d line %d column %d", basename(file),
                              as.integer(bytes[offset]), line, column))
      }
    }
  }
  expect_identical(bad, character())
  reduction <- paste(readLines(test_path("..", "..", "docs", "dev",
                                         "blr_reduction_test_matrix.md"), warn = FALSE),
                     collapse = "\n")
  expect_true(grepl("\\\\frac", reduction))
  expect_true(grepl("\\\\bmod", reduction))
  expect_true(grepl("\\\\theta", reduction))
})

test_that("current framework status remains synchronized after Phase 12A", {
  root <- test_path("..", "..")
  plan <- paste(readLines(file.path(root, "docs", "dev", "blr_framework_implementation_plan.md"),
                          warn = FALSE), collapse = "\n")
  matrix <- paste(readLines(file.path(root, "docs", "dev", "blr_model_capability_matrix.md"),
                            warn = FALSE), collapse = "\n")
  for (text in list(plan, matrix)) {
    expect_match(text, "Phase 16A experimental packed-BED BayesC routes explicitly disposed", fixed = TRUE)
    expect_match(text, "Canonical ordinary CSR", fixed = TRUE)
    expect_match(text, "Canonical scheduled", fixed = TRUE)
    expect_match(text, "historical", ignore.case = TRUE)
  }
})

phase12_scheduled_spec <- function(full) {
  list(
    execution = list(marker_count = 2L, trait_count = 1L, iterations = 4L,
      burnin = 1L, thinning = 1L, chains = 1L, cores = 1L, seed = 1L,
      chain_seeds = 2L, keep_chains = FALSE),
    rng_ownership = list(engine_owner = "chain", distribution_owner = "chain",
      lifetime = "one_chain_execution", worker_thread_owner = "none",
      fit_persistent_distribution_state = FALSE),
    sweep = list(full_sweep_every = full, iteration_zero_is_full = TRUE),
    skip = list(null_skip_base = 1L, null_skip_max = 1L, burnin_only = FALSE,
      growth_rule = "probability_adaptive"),
    candidate = list(threshold = 0, lifetime = 0L),
    neighbor = list(enabled = FALSE, difference_threshold = 0,
      maximum_neighbors = 0L, friend_marker_count = 2L,
      shared_read_only = TRUE, storage_outlives_execution = TRUE),
    state = list(scheduled_at = c(0L, 0L), last_updated = c(-1L, -1L),
      candidate = c(0L, 0L), in_candidate_list = c(0L, 0L),
      in_active_list = c(0L, 0L), last_interesting = c(-1L, -1L))
  )
}

test_that("full_sweep_every zero policy is consistent and behavior preserving", {
  expect_error(sblr:::blr_phase10a_validate_scheduled_execution_cpp(
    phase12_scheduled_spec(-1L)), "non-negative")
  for (value in c(0L, 1L, 99L)) {
    expect_true(sblr:::blr_phase10a_validate_scheduled_execution_cpp(
      phase12_scheduled_spec(value))$validated)
  }
  root <- test_path("..", "..")
  binding <- paste(readLines(file.path(root, "src", "st_cpg_omp_csr_scheduled.cpp"),
                             warn = FALSE), collapse = "\n")
  core <- paste(readLines(file.path(root, "src", "blr_csr_scheduled_bayesc_core_impl.h"),
                          warn = FALSE), collapse = "\n")
  expect_match(binding, "full_sweep_every < 0", fixed = TRUE)
  expect_match(core, "full_sweep_every <= 0", fixed = TRUE)
})

test_that("posterior pi aggregation is native and converter is presentation only", {
  root <- test_path("..", "..")
  types <- paste(readLines(file.path(root, "src", "blr_csr_scheduled_bayesc_types.h"),
                           warn = FALSE), collapse = "\n")
  core <- paste(readLines(file.path(root, "src", "blr_csr_scheduled_bayesc_core_impl.h"),
                          warn = FALSE), collapse = "\n")
  binding <- paste(readLines(file.path(root, "src", "st_cpg_omp_csr_scheduled.cpp"),
                             warn = FALSE), collapse = "\n")
  expect_source_count("arma::vec mean_pi;", types, 1L)
  expect_source_count("result.mean_pi=std::move(mean_pi);", core, 1L)
  expect_source_count("const double mean_pi=result.mean_pi(tu);", binding, 1L)
  expect_source_count("mean_pi+=result.pis", binding, 0L)
})

test_that("canonical scheduled cores have hardened diagnostic boundaries", {
  root <- test_path("..", "..")
  csr <- paste(readLines(file.path(root, "src", "blr_csr_scheduled_bayesc_core_impl.h"),
                         warn = FALSE), collapse = "\n")
  bed <- paste(readLines(file.path(root, "src", "blr_bed_scheduled_bayesc_core_impl.h"),
                         warn = FALSE), collapse = "\n")
  aggregate <- paste(readLines(file.path(root, "src", "blr_bed_scheduled_bayesc_aggregate_impl.h"),
                               warn = FALSE), collapse = "\n")
  expect_source_count("std::cout", csr, 0L)
  expect_source_count("Rcpp::Rcout", csr, 0L)
  expect_source_count("std::cout", aggregate, 0L)
  expect_source_count("Rcpp::Rcout", aggregate, 0L)
  expect_source_count("std::cout", bed, 1L)
  expect_match(bed, "progress_every > 0", fixed = TRUE)
})

test_that("canonical binding sources do not retain full commented snapshots", {
  root <- test_path("..", "..")
  for (file in c("st_cpg_omp_csr_scheduled.cpp",
                 "st_cpg_omp_individual_scheduled_chains.cpp")) {
    text <- paste(readLines(file.path(root, "src", file), warn = FALSE), collapse = "\n")
    expect_source_count("// // [[Rcpp::depends", text, 0L)
    commented_code <- strsplit(text, "\n", fixed = TRUE)[[1L]]
    code_like <- grepl("^\\s*//\\s+(for|if|while|struct|class|Rcpp::List)\\b", commented_code)
    expect_lt(max(rle(code_like)$lengths[rle(code_like)$values], 0L), 100L)
  }
})

test_that("fast and extended framework workflows are independently visible", {
  root <- test_path("..", "..")
  fast <- paste(readLines(file.path(root, ".github", "workflows", "blr-framework.yml"),
                          warn = FALSE), collapse = "\n")
  extended <- paste(readLines(file.path(root, ".github", "workflows",
                                      "blr-framework-extended.yml"), warn = FALSE),
                      collapse = "\n")
  for (needle in c("push:", "pull_request:", "workflow_dispatch:",
                   "Rcpp::compileAttributes", "devtools::test", "R CMD check"))
    expect_match(fast, needle, fixed = TRUE)
  expect_match(extended, "workflow_dispatch:", fixed = TRUE)
  expect_match(extended, "SBLR_RUN_EXTENDED_REPRODUCIBILITY", fixed = TRUE)
  expect_match(extended, "measure_peak_rss.R --smoke", fixed = TRUE)
})

test_that("peak RSS tool measures a child process when explicitly enabled", {
  tool <- test_path("..", "..", "tools", "benchmarks", "measure_peak_rss.R")
  expect_true(file.exists(tool))
  skip_if_not(identical(Sys.getenv("SBLR_RUN_PEAK_RSS"), "true"))
  skip_if_not_installed("processx")
  skip_if_not_installed("ps")
  env <- new.env(parent = globalenv())
  sys.source(tool, envir = env)
  result <- env$measure_peak_rss(file.path(R.home("bin"), "Rscript"),
    c("-e", "x<-raw(8*1024^2);Sys.sleep(0.25)"), 0.02)
  expect_identical(result$exit_status, 0L)
  expect_gt(result$peak_rss_bytes, 0)
  expect_gte(result$peak_rss_bytes, result$final_sampled_rss_bytes)
  expect_gt(result$sample_count, 0L)
  expect_true(nzchar(result$platform))
})
