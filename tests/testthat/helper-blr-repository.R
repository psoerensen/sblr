# ---- consolidated from tests/testthat/helper-source-architecture.R ----
source_match_positions <- function(pattern, text, fixed = TRUE) {
  hits <- gregexpr(pattern, text, fixed = fixed, perl = !fixed)[[1L]]
  if (length(hits) == 1L && hits[[1L]] == -1L) integer() else hits
}

source_match_count <- function(pattern, text, fixed = TRUE) {
  length(source_match_positions(pattern, text, fixed = fixed))
}

expect_source_count <- function(pattern, text, expected, fixed = TRUE) {
  testthat::expect_identical(
    source_match_count(pattern, text, fixed = fixed),
    as.integer(expected)
  )
}

expect_source_forbidden <- function(text, patterns, fixed = TRUE) {
  for (pattern in patterns) {
    testthat::expect_identical(source_match_count(pattern, text, fixed), 0L)
  }
}

expect_source_hashes <- function(root, expected) {
  testthat::expect_identical(
    unname(tools::md5sum(file.path(root, names(expected)))),
    unname(expected)
  )
}

reference_first_difference <- function(expected, observed, path = "root") {
  if (!identical(typeof(expected), typeof(observed)))
    return(list(path=path,type=c(expected=typeof(expected),observed=typeof(observed))))
  if (!identical(dim(expected), dim(observed)))
    return(list(path=path,dimensions=list(expected=dim(expected),observed=dim(observed))))
  if (!identical(names(expected), names(observed)) ||
      !identical(dimnames(expected), dimnames(observed)) ||
      !identical(class(expected), class(observed)))
    return(list(path=path,names=list(expected=names(expected),observed=names(observed)),
      dimnames=list(expected=dimnames(expected),observed=dimnames(observed)),
      class=list(expected=class(expected),observed=class(observed))))
  if (is.list(expected)) {
    if (length(expected) != length(observed))
      return(list(path=path,length=c(expected=length(expected),observed=length(observed))))
    for (i in seq_along(expected)) {
      child <- reference_first_difference(expected[[i]], observed[[i]],
        paste0(path,"$",if(length(names(expected))) names(expected)[i] else i))
      if (!is.null(child)) return(child)
    }
    return(NULL)
  }
  if (!identical(expected, observed)) {
    at <- which(expected != observed | xor(is.na(expected),is.na(observed)))[1L]
    return(list(path=path,index=at,expected=expected[at],observed=observed[at]))
  }
  NULL
}

run_reference_fresh_process <- function(fun, args = list()) {
  testthat::skip_if_not_installed("callr")
  callr::r(fun, args)
}

# ---- consolidated from tests/testthat/helper-blr-test-contracts.R ----
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

