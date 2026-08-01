tol <- 2e-12

select_rank <- function(values, prop, positive_tol = 1e-10) {
  stopifnot(length(prop) == 1L, is.finite(prop), prop > 0, prop < 1)
  positive <- sort(values[is.finite(values) & values > positive_tol],
                   decreasing = TRUE)
  if (!length(positive)) stop("no positive eigenvalues")
  cumulative <- cumsum(positive) / sum(positive)
  rank <- which(cumulative > prop)[1L]
  stopifnot(!is.na(rank), rank >= 1L)
  list(rank = rank, positive = positive, cumulative = cumulative)
}

make_operator <- function(A, score, prop, positive_tol = 1e-10) {
  stopifnot(is.matrix(A), nrow(A) == ncol(A), length(score) == nrow(A),
            all(is.finite(A)), all(is.finite(score)))
  d <- diag(A)
  stopifnot(all(d > 0))
  C <- diag(1 / sqrt(d)) %*% A %*% diag(1 / sqrt(d))
  eig <- eigen((C + t(C)) / 2, symmetric = TRUE)
  selected <- select_rank(eig$values, prop, positive_tol)
  keep <- which(eig$values > positive_tol)[seq_len(selected$rank)]
  V <- eig$vectors[, keep, drop = FALSE]
  lambda <- eig$values[keep]
  Q <- diag(sqrt(lambda), nrow = length(lambda)) %*% t(V) %*%
    diag(sqrt(d))
  w <- diag(1 / sqrt(lambda), nrow = length(lambda)) %*% t(V) %*%
    diag(1 / sqrt(d)) %*% score
  list(Q = Q, w = drop(w), diagonal = colSums(Q * Q),
       retained = lambda, positive = selected$positive,
       rank = selected$rank)
}

rebuild <- function(op, beta) drop(op$w - op$Q %*% beta)
corrected_rhs <- function(op, j, beta, residual) {
  drop(crossprod(op$Q[, j], residual)) + op$diagonal[j] * beta[j]
}
apply_difference <- function(op, j, difference, residual) {
  residual - op$Q[, j] * difference
}

# Known spectrum and exact strict boundary.
lambda <- c(5, 2, 1, 0.05, 0.005)
mass <- cumsum(lambda) / sum(lambda)
stopifnot(select_rank(lambda, 0.995)$rank == 4L,
          select_rank(lambda, mass[1L])$rank == 2L,
          select_rank(lambda, mass[3L])$rank == 4L,
          select_rank(lambda, mass[4L])$rank == 5L,
          select_rank(lambda, 1e-8)$rank == 1L,
          select_rank(lambda, 1 - 1e-12)$rank == 5L)

V <- qr.Q(qr(matrix(c(1, 2, 3, 4, 5,
                       2, -1, 0, 1, 3,
                       0, 1, 2, -2, 1,
                       3, 0, -1, 2, 1,
                       1, 1, 0, 1, -1), 5, 5)))
C <- V %*% diag(lambda) %*% t(V)
d <- c(7, 9, 11, 13, 15)
A <- diag(sqrt(d)) %*% C %*% diag(sqrt(d))
score <- c(2, -1, 0.5, 3, -2)

# Full positive rank via prop above the penultimate cumulative boundary.
full_prop <- (mass[4L] + 1) / 2
full <- make_operator(A, score, full_prop)
stopifnot(full$rank == 5L,
          max(abs(crossprod(full$Q) - A)) < tol,
          max(abs(drop(crossprod(full$Q, full$w)) - score)) < tol,
          max(abs(full$diagonal - diag(A))) < tol)

beta <- c(0.1, -0.2, 0.3, 0, 0.15)
r <- rebuild(full, beta)
j <- 3L
u_direct <- score[j] - sum(A[j, -j] * beta[-j])
stopifnot(abs(corrected_rhs(full, j, beta, r) - u_direct) < tol)
beta_new <- beta
beta_new[j] <- -0.05
r_updated <- apply_difference(full, j, beta_new[j] - beta[j], r)
stopifnot(max(abs(r_updated - rebuild(full, beta_new))) < tol)

yy <- 40
sse_a <- yy - 2 * sum(beta * score) + drop(crossprod(beta, A %*% beta))
sse_q <- yy - sum(full$w^2) + sum(r^2)
vg_a <- drop(crossprod(beta, A %*% beta))
vg_q <- sum(drop(full$Q %*% beta)^2)
stopifnot(abs(sse_a - sse_q) < tol, abs(vg_a - vg_q) < tol)

# Truncated dense oracle constructed from the exact retained factor.
truncated <- make_operator(A, score, 0.90)
A_tilde <- crossprod(truncated$Q)
s_tilde <- drop(crossprod(truncated$Q, truncated$w))
r <- rebuild(truncated, beta)
stopifnot(max(abs(truncated$diagonal - diag(A_tilde))) < tol,
          abs(corrected_rhs(truncated, j, beta, r) -
              (s_tilde[j] - sum(A_tilde[j, -j] * beta[-j]))) < tol,
          abs(sum(drop(truncated$Q %*% beta)^2) -
              drop(crossprod(beta, A_tilde %*% beta))) < tol,
          abs(sum(beta * s_tilde) -
              drop(crossprod(truncated$Q %*% beta, truncated$w))) < tol)

updates <- list(c(2L, 0.4), c(5L, -0.3), c(1L, 0), c(4L, 0.25))
r_seq <- r
b_seq <- beta
for (update in updates) {
  idx <- update[1L]
  value <- update[2L]
  r_seq <- apply_difference(truncated, idx, value - b_seq[idx], r_seq)
  b_seq[idx] <- value
  stopifnot(max(abs(r_seq - rebuild(truncated, b_seq))) < tol)
}

sse_dense <- yy - 2 * sum(b_seq * s_tilde) +
  drop(crossprod(b_seq, A_tilde %*% b_seq))
sse_reduced <- yy - sum(truncated$w^2) + sum(r_seq^2)
stopifnot(abs(sse_dense - sse_reduced) < tol)

# Two blocks are independent and quadratic forms add.
op1 <- make_operator(A[1:2, 1:2], score[1:2], 0.8)
op2 <- make_operator(A[3:5, 3:5], score[3:5], 0.8)
b1 <- beta[1:2]
b2 <- beta[3:5]
qsum <- sum(drop(op1$Q %*% b1)^2) + sum(drop(op2$Q %*% b2)^2)
block_A <- matrix(0, 5, 5)
block_A[1:2, 1:2] <- crossprod(op1$Q)
block_A[3:5, 3:5] <- crossprod(op2$Q)
stopifnot(abs(qsum - drop(crossprod(beta, block_A %*% beta))) < tol)
r2_before <- rebuild(op2, b2)
r1 <- rebuild(op1, b1)
r1 <- apply_difference(op1, 1L, 0.2, r1)
stopifnot(identical(r2_before, rebuild(op2, b2)), length(r1) == nrow(op1$Q))

# A within-block permutation covaries Q'Q, Q'w, diagonals, and conditionals.
perm <- c(3, 1, 5, 2, 4)
permuted <- make_operator(A[perm, perm], score[perm], 0.90)
stopifnot(max(abs(crossprod(permuted$Q) - A_tilde[perm, perm])) < tol,
          max(abs(drop(crossprod(permuted$Q, permuted$w)) -
                  s_tilde[perm])) < tol)

cat("STAGE A COMPLETE — LOW-RANK MATHEMATICS AND OPERATOR CONTRACT VERIFIED\n")
