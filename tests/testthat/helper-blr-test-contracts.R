blr_test_dir <- normalizePath(testthat::test_path(), winslash = "/",
                              mustWork = TRUE)

blr_fixture_path <- function(...) {
  testthat::test_path("fixtures", ...)
}

blr_is_source_root <- function(path) {
  if (is.null(path) || length(path) != 1L || is.na(path) || !nzchar(path) ||
      !dir.exists(path)) {
    return(FALSE)
  }
  required <- c("DESCRIPTION", "R", "src", file.path("tests", "testthat"))
  if (!all(file.exists(file.path(path, required)))) return(FALSE)
  desc <- tryCatch(read.dcf(file.path(path, "DESCRIPTION")), error = function(e) NULL)
  !is.null(desc) && identical(unname(desc[1L, "Package"]), "sblr")
}

blr_find_source_root <- function(explicit = Sys.getenv("SBLR_SOURCE_ROOT", ""),
                                 start = getwd(), test_dir = blr_test_dir) {
  candidates <- c(
    explicit,
    normalizePath(file.path(test_dir, "..", ".."), winslash = "/",
                  mustWork = FALSE),
    normalizePath(start, winslash = "/", mustWork = FALSE)
  )
  parent <- normalizePath(file.path(start, ".."), winslash = "/", mustWork = FALSE)
  candidates <- unique(c(candidates, parent,
                         normalizePath(file.path(parent, ".."), winslash = "/",
                                       mustWork = FALSE)))
  candidates <- candidates[nzchar(candidates)]
  found <- candidates[vapply(candidates, blr_is_source_root, logical(1))]
  if (length(found)) found[[1L]] else NA_character_
}

blr_source_root <- blr_find_source_root()

blr_has_source_tree <- function() {
  blr_is_source_root(blr_source_root)
}

blr_skip_if_no_source_tree <- function() {
  if (!blr_has_source_tree()) {
    testthat::skip(paste(
      "Source-tree architecture assertion:",
      "repository files are unavailable in installed-package check context."
    ))
  }
  invisible(TRUE)
}

blr_repo_path <- function(...) {
  blr_skip_if_no_source_tree()
  file.path(blr_source_root, ...)
}

blr_source_text <- function(path) {
  paste(readLines(blr_repo_path(path), warn = FALSE), collapse = "\n")
}

blr_expect_fixture_md5 <- function(paths, expected) {
  resolved <- vapply(paths, function(path) {
    parts <- strsplit(gsub("\\\\", "/", path), "/", fixed = TRUE)[[1L]]
    fixture_at <- match("fixtures", parts)
    if (!is.na(fixture_at)) parts <- parts[seq.int(fixture_at + 1L, length(parts))]
    do.call(blr_fixture_path, as.list(parts))
  }, character(1))
  testthat::expect_identical(unname(tools::md5sum(resolved)), unname(expected))
}

blr_mt_public_source <- function() {
  adapter <- blr_source_text("src/mtblr.cpp")
  first <- regexpr("std::vector<std::vector<std::vector<double>>>  mtblr(",
                   adapter, fixed = TRUE)[1]
  last <- regexpr("std::vector<std::vector<std::vector<double>>>  mtblr_eigen(",
                  adapter, fixed = TRUE)[1]
  substr(adapter, first, last - 1L)
}
