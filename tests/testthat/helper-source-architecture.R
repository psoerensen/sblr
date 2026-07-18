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
