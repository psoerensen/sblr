args <- commandArgs(trailingOnly = TRUE)
run_tests <- "--run" %in% args
args <- setdiff(args, "--run")
root <- normalizePath(if (length(args)) args[[1L]] else ".",
  winslash = "/", mustWork = TRUE)
test_dir <- file.path(root, "tests", "testthat")
files <- sort(list.files(test_dir, pattern = "^test.*\\.R$", full.names = TRUE))

count_pattern <- function(lines, pattern) sum(grepl(pattern, lines, perl = TRUE))
inventory <- do.call(rbind, lapply(files, function(path) {
  lines <- readLines(path, warn = FALSE)
  data.frame(
    file = basename(path),
    phase_or_subsystem = sub("^test-", "", sub("\\.R$", "", basename(path))),
    test_blocks = count_pattern(lines, "test_that\\s*\\("),
    expectation_calls = count_pattern(lines, "expect_[A-Za-z0-9_]+\\s*\\("),
    fixture_references = count_pattern(lines, "fixtures|readRDS"),
    fresh_process_calls = count_pattern(lines, "callr::r|fresh.process"),
    thread_environment_calls = count_pattern(lines, "OMP_NUM_THREADS|ncores.*2"),
    source_assertions = count_pattern(lines,
      "expect_source_|source_match_count|gregexpr|grepl\\("),
    source_hashes = count_pattern(lines, "md5sum.*src|region_md5"),
    fixture_hashes = count_pattern(lines, "md5sum.*fixture|fixture_md5"),
    schema_checks = count_pattern(lines, "schema|names\\(|dim\\("),
    scientific_checks = count_pattern(lines,
      "cov|pim|positive|finite|binary|residual|orientation"),
    protected_backend_checks = count_pattern(lines, "protected|unchanged"),
    current_unique_responsibility = sub("^test-", "",
      sub("\\.R$", "", basename(path))),
    duplication_candidates = "reviewed; none assigned after Phase 17F2",
    recommended_permanent_owner = basename(path),
    tier = if (grepl("extended", basename(path))) "extended" else "ordinary",
    stringsAsFactors = FALSE
  )
}))

inventory$actual_expectations <- NA_integer_
inventory$runtime_seconds <- NA_real_
if (run_tests) {
  pkgload::load_all(root, compile = FALSE, quiet = TRUE)
  results <- as.data.frame(testthat::test_dir(test_dir, reporter = "silent"))
  measured <- aggregate(cbind(passed, failed, warning, skipped, real) ~ file,
    results, sum)
  matched <- match(inventory$file, measured$file)
  inventory$actual_expectations <- ifelse(is.na(matched), 0L,
    measured$passed[matched] + measured$failed[matched] +
      measured$warning[matched])
  inventory$runtime_seconds <- ifelse(is.na(matched), 0, measured$real[matched])
}

out <- file.path(root, "tools", "audit", "blr_test_suite_inventory.csv")
write.csv(inventory, out, row.names = FALSE)
cat(sprintf("files=%d test_that_blocks=%d static_expectation_calls=%d\n",
  nrow(inventory), sum(inventory$test_blocks), sum(inventory$expectation_calls)))
cat(sprintf("wrote %s\n", out))
