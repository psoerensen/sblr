#' Convert a Text LD Matrix to CSR Format
#'
#' Reads a text linkage disequilibrium matrix and returns a compressed sparse
#' row representation.
#'
#' @param filename Path to the text LD file.
#' @param mchr Chromosome index.
#' @param msize Number of markers.
#' @param r2 Minimum squared-correlation threshold.
#' @param onebased Return one-based column indices.
#' @return A list representing a CSR sparse matrix.
#' @export
readLD_to_CSR <- function(filename,
                          mchr,
                          msize,
                          r2 = 0.01,
                          onebased = FALSE) {
  stopifnot(is.character(filename), length(filename) == 1)
  stopifnot(file.exists(filename))
  stopifnot(is.numeric(mchr), mchr > 0)
  stopifnot(is.numeric(msize), msize >= 0)
  stopifnot(is.numeric(r2), r2 >= 0)

  res <- readLD_to_CSR_R(
    filename,
    as.integer(mchr),
    as.integer(msize),
    as.numeric(r2),
    as.logical(onebased)
  )

  class(res) <- "csr_matrix"
  res
}
