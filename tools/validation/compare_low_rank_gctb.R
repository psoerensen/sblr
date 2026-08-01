#!/usr/bin/env Rscript

# Deterministic scale crosswalk to GCTB@cc7fa7d765c83a89c6375946cf77fe50ba1a317e.
# See docs/dev/stblr_low_rank_gctb_crosswalk.md for the source-level audit.

set.seed(4102)
Z <- scale(matrix(rnorm(20 * 7), 20, 7), center = TRUE, scale = FALSE)
A <- crossprod(Z); D <- diag(A)
C <- diag(1 / sqrt(D)) %*% A %*% diag(1 / sqrt(D))
eig <- eigen(C, symmetric = TRUE)
keep <- which(eig$values > 1e-10)
V <- eig$vectors[, keep, drop = FALSE]; lambda <- eig$values[keep]
s <- drop(crossprod(Z, sin(seq_len(nrow(Z)))))

Q_gctb <- diag(sqrt(lambda)) %*% t(V)
w_gctb <- drop(diag(1 / sqrt(lambda)) %*% t(V) %*% diag(1 / sqrt(D)) %*% s)
Q_sblr <- Q_gctb %*% diag(sqrt(D))
w_sblr <- w_gctb

stopifnot(
  max(abs(crossprod(Q_sblr) - A)) < 1e-8,
  max(abs(crossprod(Q_sblr, w_sblr) - s)) < 1e-8
)
beta <- seq(-0.1, 0.1, length.out = ncol(Z))
wcorr_gctb <- w_gctb - Q_gctb %*% (sqrt(D) * beta)
r_sblr <- w_sblr - Q_sblr %*% beta
stopifnot(max(abs(wcorr_gctb - r_sblr)) < 1e-11)
for (marker in seq_len(ncol(Z))) {
  rhs_gctb <- drop(crossprod(Q_gctb[, marker], wcorr_gctb) * sqrt(D[marker]) +
                     crossprod(Q_sblr[, marker]) * beta[marker])
  rhs_sblr <- drop(crossprod(Q_sblr[, marker], r_sblr) +
                     crossprod(Q_sblr[, marker]) * beta[marker])
  stopifnot(abs(rhs_gctb - rhs_sblr) < 1e-10)
}
cat("GCTB normalized-to-cross-product factor, score, residual, and conditional crosswalk passed\n")
