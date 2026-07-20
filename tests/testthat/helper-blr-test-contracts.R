blr_test_root <- normalizePath(file.path(testthat::test_path(), "..", ".."),
  winslash = "/", mustWork = TRUE)

blr_source_text <- function(path) {
  paste(readLines(file.path(blr_test_root, path), warn = FALSE), collapse = "\n")
}

blr_expect_fixture_md5 <- function(paths, expected) {
  testthat::expect_identical(
    unname(tools::md5sum(file.path(blr_test_root, paths))),
    unname(expected)
  )
}

blr_mt_public_source <- function() {
  adapter <- blr_source_text("src/mtblr.cpp")
  first <- regexpr("std::vector<std::vector<std::vector<double>>>  mtblr(",
    adapter, fixed = TRUE)[1]
  last <- regexpr("std::vector<std::vector<std::vector<double>>>  mtblr_eigen(",
    adapter, fixed = TRUE)[1]
  substr(adapter, first, last - 1L)
}
