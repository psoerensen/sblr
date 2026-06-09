estimate_vg_by_set <- function(fit, sets, n) {
 b <- fit$bm
 q <- fit$wy - fit$r   # approximately X'X b

 out <- matrix(NA_real_, nrow = length(sets), ncol = ncol(b))
 rownames(out) <- names(sets)
 colnames(out) <- colnames(b)

 for (s in seq_along(sets)) {
  idx <- sets[[s]]
  for (t in seq_len(ncol(b))) {
   out[s, t] <- sum(b[idx, t] * q[idx, t]) / n
  }
 }

 out
}

estimate_pip_enrichment <- function(dm, sets) {
 out <- matrix(NA_real_, nrow = length(sets), ncol = ncol(dm))
 rownames(out) <- names(sets)
 colnames(out) <- colnames(dm)

 all_idx <- seq_len(nrow(dm))

 for (s in seq_along(sets)) {
  idx <- sets[[s]]
  bg <- setdiff(all_idx, idx)

  for (t in seq_len(ncol(dm))) {
   mean_in <- mean(dm[idx, t])
   mean_bg <- mean(dm[bg, t])
   out[s, t] <- mean_in / mean_bg
  }
 }

 out
}

estimate_weighted_vg_by_set <- function(fit, sets, n) {
 b <- fit$bm
 dm <- fit$dm
 q <- fit$wy - fit$r

 out <- matrix(NA_real_, nrow = length(sets), ncol = ncol(b))
 rownames(out) <- names(sets)
 colnames(out) <- colnames(b)

 for (s in seq_along(sets)) {
  idx <- sets[[s]]
  for (t in seq_len(ncol(b))) {
   out[s, t] <- sum(dm[idx, t] * b[idx, t] * q[idx, t]) / n
  }
 }

 out
}


estimate_G_by_set <- function(fit, sets, n) {
 b <- fit$bm
 Q <- fit$wy - fit$r
 nt <- ncol(b)

 out <- vector("list", length(sets))
 names(out) <- names(sets)

 for (s in seq_along(sets)) {
  idx <- sets[[s]]
  G <- matrix(0, nt, nt)
  rownames(G) <- colnames(G) <- colnames(b)

  for (t1 in seq_len(nt)) {
   for (t2 in t1:nt) {
    g12 <- sum(b[idx, t1] * Q[idx, t2]) / n
    g21 <- sum(b[idx, t2] * Q[idx, t1]) / n
    gij <- 0.5 * (g12 + g21)

    G[t1, t2] <- gij
    G[t2, t1] <- gij
   }
  }

  out[[s]] <- G
 }

 out
}


estimate_annotation_vb <- function(fit, sets, n, min_eff = 1, shrink = 0.5) {
 b  <- fit$bm
 dm <- fit$dm
 q  <- fit$wy - fit$r

 nt <- ncol(b)
 K <- length(sets)

 vb_ann <- matrix(NA_real_, K, nt)
 vg_ann <- matrix(NA_real_, K, nt)
 meff   <- matrix(NA_real_, K, nt)

 rownames(vb_ann) <- rownames(vg_ann) <- rownames(meff) <- names(sets)
 colnames(vb_ann) <- colnames(vg_ann) <- colnames(meff) <- colnames(b)

 for (k in seq_along(sets)) {
  idx <- sets[[k]]

  for (t in seq_len(nt)) {
   vg <- sum(b[idx, t] * q[idx, t]) / n
   neff <- sum(dm[idx, t])

   vg_ann[k, t] <- max(vg, 0)
   meff[k, t] <- neff
   vb_ann[k, t] <- max(vg_ann[k, t] / max(neff, min_eff), 1e-12)
  }
 }

 # shrink annotation-specific scales toward global scale
 vb_global <- colMeans(vb_ann, na.rm = TRUE)

 for (t in seq_len(nt)) {
  vb_ann[, t] <- shrink * vb_ann[, t] + (1 - shrink) * vb_global[t]
 }

 list(vb = vb_ann, vg = vg_ann, meff = meff)
}


estimate_annotation_pi <- function(dm, sets, alpha = 1) {
 nt <- ncol(dm)
 K <- length(sets)

 pi_ann <- matrix(NA_real_, K, nt)
 rownames(pi_ann) <- names(sets)
 colnames(pi_ann) <- colnames(dm)

 for (k in seq_along(sets)) {
  idx <- sets[[k]]

  for (t in seq_len(nt)) {
   pi_ann[k, t] <- (alpha + sum(dm[idx, t])) / (2 * alpha + length(idx))
  }
 }

 pmin(pmax(pi_ann, 1e-8), 1 - 1e-8)
}

derive_marker_vb_from_annotation_model <- function(global_vb, annot, gamma, cap = c(0.1, 10)) {
 # annot: m x K annotation matrix
 # gamma: K x nt annotation coefficients on log-scale
 # global_vb: length nt

 eta <- annot %*% gamma
 vb_marker <- matrix(NA_real_, nrow(annot), ncol(gamma))

 for (t in seq_len(ncol(gamma))) {
  mult <- exp(eta[, t])
  mult <- pmin(pmax(mult, cap[1]), cap[2])
  vb_marker[, t] <- global_vb[t] * mult
 }

 colnames(vb_marker) <- colnames(gamma)
 vb_marker
}
