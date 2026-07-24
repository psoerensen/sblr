# Small same-BED public multivariate block-eigen workflow.
write_bed <- function(path, dosage) {
  code <- function(x) ifelse(is.na(x), 1L,
    c(`0` = 3L, `1` = 2L, `2` = 0L)[as.character(x)])
  packed <- unlist(lapply(seq_len(nrow(dosage)), function(marker) {
    values <- as.integer(code(dosage[marker, ]))
    values <- c(values, rep(0L, (-length(values)) %% 4L))
    vapply(seq(1L, length(values), by = 4L), function(i)
      sum(values[i:(i + 3L)] * c(1L, 4L, 16L, 64L)), integer(1))
  }))
  writeBin(as.raw(c(0x6c, 0x1b, 0x01, packed)), path)
}

dosage <- rbind(
  c(0, 1, 2, 0, 1, 2, 1, 0),
  c(2, 1, 0, 2, 1, 0, NA, 2),
  c(0, 0, 1, 1, 2, 2, 1, 0),
  c(2, 2, 1, 1, 0, 0, 1, 2)
)
bed <- tempfile(fileext = ".bed")
write_bed(bed, dosage)
ids <- paste0("id", seq_len(ncol(dosage)))
markers <- paste0("m", seq_len(nrow(dosage)))
Glist <- list(
  bedfiles = bed, n = ncol(dosage), ids = ids,
  rsids = list(markers), rsidsLD = list(markers),
  af = list(c(.25, .35, .40, .30))
)
y <- cbind(T1 = seq(-1, 1, length.out = 8),
           T2 = cos(seq(0, pi, length.out = 8)))
rownames(y) <- ids
stats <- make_summary_stats(Glist, y, nthreads = 1)

shared <- mtblr_block_eigen(
  stats, Glist, block_start = c(1L, 3L),
  eigen_filter = "hard_truncate",
  nit = 20, nburn = 5, seed = 17
)

trait_specific <- mtblr_block_eigen(
  stats, Glist, block_start = list(c(1L, 3L), c(1L, 2L, 4L)),
  operator_sharing = "trait_specific",
  eigen_filter = c("hard_truncate", "ridge_fixed"),
  eigen_eta = c(0, .5),
  nit = 20, nburn = 5, seed = 17
)

shared$block_diagnostics
trait_specific$block_diagnostics
trait_specific$alignment
