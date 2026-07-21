root <- normalizePath(".", winslash = "/", mustWork = TRUE)
files <- list.files(file.path(root, "tests", "testthat"), "[.]R$", full.names = TRUE)
patterns <- list(
  direct_source_file_read = "readLines\\([^\n]*(src|R/|docs|[.]github|tools|examples)",
  production_R_source_fallback = "source_sblr_test_file|source\\([^\n]*R/",
  cross_test_file_parsing = "eval\\(parse|test-blr-framework-phase17i[.]R",
  source_or_document_hash = "md5sum",
  source_root_resolution = "blr_repo_path|blr_source_text",
  fixture_root_resolution = "blr_fixture_path|test_path\\(\"fixtures\""
)
inventory <- do.call(rbind, lapply(files, function(file) {
  lines <- readLines(file, warn = FALSE)
  do.call(rbind, lapply(names(patterns), function(category) {
    hit <- grep(patterns[[category]], lines)
    if (!length(hit)) return(NULL)
    data.frame(check_section="tests", test_file=basename(file), test_block=NA,
      failure_or_warning="static occurrence", category=category,
      path_or_field=trimws(lines[hit]), source_tree_required=grepl("source", category),
      portable_fixture=category == "fixture_root_resolution",
      installed_test_required=category == "fixture_root_resolution",
      proposed_resolution="classified by Phase 17J2 contract")
  }))
}))
print(inventory, row.names = FALSE)
