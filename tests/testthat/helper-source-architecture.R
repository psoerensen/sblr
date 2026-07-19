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
