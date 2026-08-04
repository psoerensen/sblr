#!/usr/bin/env Rscript

# Read-only sparse spectral audit for the preserved Study 06 CSR operator.
# The matrix remains sparse; no genome-wide dense matrix is constructed.

root <- normalizePath(if (length(commandArgs(TRUE))) commandArgs(TRUE)[1L] else
  "../sblrbench", winslash = "/", mustWork = TRUE)
prefix <- file.path(root,
  "results/local/06_annotation_models/checkpoints/ld",
  "training_ld_train1400-test600-seed3101")
meta <- readLines(paste0(prefix, ".meta.txt"), warn = FALSE)
value <- function(key) sub(paste0("^", key, "="), "",
  meta[startsWith(meta, paste0(key, "="))])
m <- as.integer(value("n_variants"))
ptr_raw <- readBin(paste0(prefix, ".row_ptr.u64.bin"), "raw",
  n = 8L * (m + 1L))
ptr_bytes <- matrix(as.integer(ptr_raw), nrow = 8L)
# Study 06 has fewer than 2^32 entries, so the low four bytes are exact.
stopifnot(all(ptr_bytes[5:8, ] == 0L))
ptr <- ptr_bytes[1L, ] + 256 * ptr_bytes[2L, ] +
  65536 * ptr_bytes[3L, ] + 16777216 * ptr_bytes[4L, ]
col <- readBin(paste0(prefix, ".col_idx.u32.0based.bin"), "integer",
  n = as.integer(value("nnz")), size = 4L, signed = FALSE,
  endian = "little") + 1L
val <- readBin(paste0(prefix, ".values.f32.bin"), "numeric",
  n = length(col), size = 4L, endian = "little")
row <- rep.int(seq_len(m), diff(ptr))
upper <- Matrix::sparseMatrix(i = row, j = col, x = val,
  dims = c(m, m), giveCsparse = TRUE)
R <- upper + Matrix::t(upper) + Matrix::Diagonal(m)
rm(upper, row, col, val, ptr)

factorization <- tryCatch(
  Matrix::Cholesky(R, perm = TRUE, LDL = FALSE, super = FALSE),
  error = identity)
positive_definite <- !inherits(factorization, "error")
minimum_factor_diagonal <- if (positive_definite) {
  L <- Matrix::expand(factorization)$L
  min(Matrix::diag(L))
} else NA_real_
cat("Study 06 hard-sparse CSR spectral audit\n")
cat("matrix:", m, "x", m, "stored upper entries:", value("nnz"), "\n")
cat("sparse Cholesky positive definite:", positive_definite, "\n")
cat(sprintf("minimum factor diagonal=%.12g\n", minimum_factor_diagonal))
if (!positive_definite) cat("factorization error:", conditionMessage(factorization), "\n")
