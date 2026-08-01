#!/usr/bin/env Rscript

# Deterministic algebraic crosswalk to qgg@bfac8b2c388afb7ae1c88019bcfef8588f81aedb.
# The pinned source remains external and read-only; Q is compared through its
# subspace and Gram operator, never by raw eigenvector signs.

set.seed(4101)
Z <- scale(matrix(rnorm(18 * 6), 18, 6), center = TRUE, scale = FALSE)
A <- crossprod(Z)
D <- diag(A)
C <- diag(1 / sqrt(D)) %*% A %*% diag(1 / sqrt(D))
eig <- eigen(C, symmetric = TRUE)
keep <- which(eig$values > 1e-10)
V <- eig$vectors[, keep, drop = FALSE]
lambda <- eig$values[keep]
Q_sblr <- diag(sqrt(lambda)) %*% t(V) %*% diag(sqrt(D))

# qgg's retained factor is the same eigen-by-marker orientation in normalized
# genotype/LD units. Restoring sblr cross-product columns gives this factor.
Q_qgg_cross_product <- diag(sqrt(lambda)) %*% t(V) %*% diag(sqrt(D))
s <- drop(crossprod(Z, seq_len(nrow(Z)) - mean(seq_len(nrow(Z)))))
w <- drop(diag(1 / sqrt(lambda)) %*% t(V) %*% diag(1 / sqrt(D)) %*% s)
beta <- seq(-0.15, 0.15, length.out = ncol(Z))
r <- w - Q_sblr %*% beta

stopifnot(
  max(abs(crossprod(Q_sblr) - crossprod(Q_qgg_cross_product))) < 1e-11,
  max(abs(crossprod(Q_sblr, w) - crossprod(Q_qgg_cross_product, w))) < 1e-11
)
for (marker in seq_len(ncol(Z))) {
  q <- Q_sblr[, marker]
  rhs_sblr <- drop(crossprod(q, r) + crossprod(q) * beta[marker])
  rhs_qgg <- drop(crossprod(Q_qgg_cross_product[, marker], r) +
                     crossprod(Q_qgg_cross_product[, marker]) * beta[marker])
  stopifnot(abs(rhs_sblr - rhs_qgg) < 1e-11)
  delta <- 0.01 * (-1)^marker
  r <- r - q * delta
  beta[marker] <- beta[marker] + delta
}
cat("qgg retained-subspace, Gram, transformed-score, conditional, and update comparison passed\n")

