#' Simulate Multi-Trait Phenotypes
#'
#' Simulates multi-trait phenotypes and marker effects from a supplied genotype
#' matrix or a qgg genotype list.
#'
#' @param Glist Optional qgg genotype list.
#' @param chr,rsids,ids Optional chromosome, marker, and individual selections.
#' @param causal_rsids Markers eligible to be sampled as causal.
#' @param W Optional genotype matrix.
#' @param n,m Simulated sample and marker counts when `W` is omitted.
#' @param nt Number of traits.
#' @param n_shared,n_specific Numbers of shared and trait-specific causal
#'   markers.
#' @param h2 Trait heritabilities.
#' @param rg,re Genetic-effect and residual correlation specifications.
#' @param effect_sd Standard deviation used to simulate marker effects.
#' @param maf_min,maf_max Minor-allele-frequency range when simulating `W`.
#' @param standardize_W Standardize marker columns before simulation.
#' @param seed Optional simulation seed.
#' @param exact_shared_cor Force the sampled shared effects to have exactly
#'   correlation `rg`.
#' @return A list containing phenotypes, marker effects, genetic values,
#'   residuals, and simulation metadata.
#' @export
mtsim <- function(
  Glist = NULL,
  chr = NULL,
  rsids = NULL,
  ids = NULL,
  causal_rsids = NULL,
  W = NULL,
  n = 1000,
  m = 1000,
  nt = 2,
  n_shared = 10,
  n_specific = 10,
  h2 = 0.5,
  rg = NULL,
  re = 0.0,
  effect_sd = 1,
  maf_min = 0.05,
  maf_max = 0.5,
  standardize_W = TRUE,
  seed = NULL,
  exact_shared_cor = FALSE
) {

 if (!is.null(seed)) {
  set.seed(seed)
 }

 if (length(h2) == 1) {
  h2 <- rep(h2, nt)
 }

 if (length(h2) != nt) {
  stop("h2 must be either a scalar or a vector of length nt.")
 }

 if (is.null(rg)) {
  rg <- diag(nt)
 }

 if (length(rg) == 1 && nt == 2) {
  rg_mat <- matrix(rg, nt, nt)
  diag(rg_mat) <- 1
  rg <- rg_mat
 }

 if (!is.matrix(rg) || !all(dim(rg) == c(nt, nt))) {
  stop("rg must be NULL, a scalar for nt = 2, or an nt x nt correlation matrix.")
 }

 if (any(abs(rg) > 1)) {
  stop("All genetic/effect correlations in rg must be between -1 and 1.")
 }

 if (!isTRUE(all.equal(diag(rg), rep(1, nt)))) {
  stop("The diagonal of rg must be 1.")
 }

 if (any(eigen(rg, symmetric = TRUE)$values <= 1e-8)) {
  stop("rg must be positive definite.")
 }

 # ------------------------------------------------------------
 # Load genotype matrix from Glist
 # ------------------------------------------------------------

 if (!is.null(Glist)) {

  if (is.null(rsids)) {
   if (!is.null(chr) && !is.null(Glist$rsidsLD)) {
    rsids <- Glist$rsidsLD[[chr]]
   } else {
    rsids <- Glist$rsids
   }
  }

  if (is.null(ids)) {
   ids <- Glist$ids
  }

  if (is.null(chr)) {
   W <- qgg::getG(Glist = Glist, rsids = rsids, ids = ids)
  } else {
   W <- qgg::getG(Glist = Glist, chr = chr, rsids = rsids, ids = ids)
  }
 }

 # ------------------------------------------------------------
 # Simulate W if not supplied
 # ------------------------------------------------------------

 if (is.null(W)) {

  maf <- runif(m, maf_min, maf_max)

  W <- sapply(maf, function(p) {
   rbinom(n, size = 2, prob = p)
  })

  W <- as.matrix(W)

  colnames(W) <- paste0("m", seq_len(m))
  rownames(W) <- paste0("id", seq_len(n))
 }

 n <- nrow(W)
 m <- ncol(W)

 if (is.null(colnames(W))) {
  colnames(W) <- paste0("m", seq_len(m))
 }

 if (is.null(rownames(W))) {
  if (!is.null(ids) && length(ids) == n) {
   rownames(W) <- ids
  } else {
   rownames(W) <- paste0("id", seq_len(n))
  }
 }

 # ------------------------------------------------------------
 # Define causal marker pool
 # ------------------------------------------------------------

 if (is.null(causal_rsids)) {
  causal_rsids <- colnames(W)
 }

 causal_pool <- which(colnames(W) %in% causal_rsids)

 if (length(causal_pool) == 0) {
  stop("None of the supplied causal_rsids were found in colnames(W).")
 }

 total_needed <- n_shared + nt * n_specific

 if (total_needed > length(causal_pool)) {
  stop(
   "Not enough eligible causal markers. Requested ",
   total_needed,
   " causal markers, but only ",
   length(causal_pool),
   " markers are available in causal_rsids."
  )
 }

 # ------------------------------------------------------------
 # Standardize marker matrix
 # ------------------------------------------------------------

 if (standardize_W) {
  W <- scale(W)
 }

 # ------------------------------------------------------------
 # Select causal markers
 # ------------------------------------------------------------

 causal_all <- sample(causal_pool, total_needed)

 shared_idx <- if (n_shared > 0) {
  causal_all[seq_len(n_shared)]
 } else {
  integer(0)
 }

 specific_idx <- vector("list", nt)
 names(specific_idx) <- paste0("D", seq_len(nt))

 start <- n_shared + 1

 for (j in seq_len(nt)) {
  if (n_specific > 0) {
   specific_idx[[j]] <- causal_all[start:(start + n_specific - 1)]
   start <- start + n_specific
  } else {
   specific_idx[[j]] <- integer(0)
  }
 }

 # ------------------------------------------------------------
 # Simulate marker effects directly
 # ------------------------------------------------------------

 B <- matrix(0, nrow = m, ncol = nt)
 colnames(B) <- paste0("D", seq_len(nt))
 rownames(B) <- colnames(W)

 # Shared causal markers: multivariate effects with correlation rg
 if (n_shared > 0) {

  if (nt == 1) {
   B[shared_idx, 1] <- rnorm(n_shared, mean = 0, sd = effect_sd)
  } else {

   if (exact_shared_cor) {

    if (n_shared <= nt) {
     stop("exact_shared_cor = TRUE requires n_shared > nt.")
    }

    B_shared <- matrix(rnorm(n_shared * nt), nrow = n_shared, ncol = nt)
    B_shared <- force_correlation(B_shared, target_cor = rg)
    B_shared <- B_shared * effect_sd

   } else {

    Sigma_b <- effect_sd^2 * rg
    Lb <- chol(Sigma_b)

    B_shared <- matrix(rnorm(n_shared * nt), nrow = n_shared, ncol = nt) %*% Lb
   }

   B[shared_idx, ] <- B_shared
  }
 }

 # Trait-specific causal markers: effects only on one trait
 if (n_specific > 0) {
  for (j in seq_len(nt)) {
   B[specific_idx[[j]], j] <- rnorm(n_specific, mean = 0, sd = effect_sd)
  }
 }

 # ------------------------------------------------------------
 # Genetic values
 # ------------------------------------------------------------

 G_raw <- W %*% B
 colnames(G_raw) <- paste0("D", seq_len(nt))
 rownames(G_raw) <- rownames(W)

 var_g_raw <- apply(G_raw, 2, var)

 if (any(var_g_raw <= 0 | !is.finite(var_g_raw))) {
  stop("At least one trait has zero or non-finite genetic variance.")
 }

 # Scale marker effects so that Var(G_j) = 1 for all traits.
 # This preserves G = W %*% B.
 scale_b <- 1 / sqrt(var_g_raw)

 B <- sweep(B, 2, scale_b, "*")

 G <- W %*% B
 colnames(G) <- paste0("D", seq_len(nt))
 rownames(G) <- rownames(W)

 var_g <- apply(G, 2, var)

 # ------------------------------------------------------------
 # Residual covariance calibrated to target h2
 # ------------------------------------------------------------

 var_e <- var_g * (1 - h2) / h2

 if (length(re) == 1) {
  Sigma_e_cor <- matrix(re, nt, nt)
  diag(Sigma_e_cor) <- 1
 } else {
  Sigma_e_cor <- re
 }

 if (!is.matrix(Sigma_e_cor) || !all(dim(Sigma_e_cor) == c(nt, nt))) {
  stop("re must be a scalar or an nt x nt residual correlation matrix.")
 }

 if (any(abs(Sigma_e_cor) > 1)) {
  stop("All residual correlations in re must be between -1 and 1.")
 }

 if (!isTRUE(all.equal(diag(Sigma_e_cor), rep(1, nt)))) {
  stop("The diagonal of re must be 1.")
 }

 if (any(eigen(Sigma_e_cor, symmetric = TRUE)$values <= 1e-8)) {
  stop("Residual correlation matrix re must be positive definite.")
 }

 sd_e <- sqrt(var_e)

 # Explicit dimensions prevent a scalar value from being treated as a matrix size.
 D_e <- diag(as.numeric(sd_e), nrow = nt, ncol = nt)

 Sigma_e <- D_e %*% Sigma_e_cor %*% D_e

 Le <- chol(Sigma_e)

 E <- matrix(rnorm(n * nt), nrow = n, ncol = nt) %*% Le
 colnames(E) <- paste0("D", seq_len(nt))
 rownames(E) <- rownames(W)

 # ------------------------------------------------------------
 # Phenotypes
 # ------------------------------------------------------------

 Y <- G + E
 colnames(Y) <- paste0("D", seq_len(nt))
 rownames(Y) <- rownames(W)

 causal <- list(
  shared = colnames(W)[shared_idx],
  specific = lapply(specific_idx, function(idx) colnames(W)[idx]),
  all = colnames(W)[sort(c(shared_idx, unlist(specific_idx)))]
 )

 B_shared_cor <- if (nt == 1 || n_shared == 0) {
  NA
 } else {
  cor(B[shared_idx, , drop = FALSE])
 }

 B_all_cor <- if (nt == 1) {
  1
 } else {
  cor(B)
 }

 list(
  y = if (nt == 1) Y[, 1] else Y,
  W = W,
  B = B,
  G = G,
  E = E,
  h2_target = h2,
  h2_observed = apply(G, 2, var) / apply(Y, 2, var),
  rg_target = rg,
  rg_observed = if (nt == 1) 1 else cor(G),
  rb_shared_observed = B_shared_cor,
  rb_all_observed = B_all_cor,
  re_target = Sigma_e_cor,
  re_observed = if (nt == 1) 1 else cor(E),
  Sigma_e = Sigma_e,
  causal = causal,
  shared_idx = shared_idx,
  specific_idx = specific_idx,
  rsids = colnames(W),
  ids = rownames(W),
  causal_rsids = causal_rsids
 )
}

force_correlation <- function(X, target_cor) {

 X <- scale(X, center = TRUE, scale = TRUE)

 current_cor <- cor(X)

 eig_current <- eigen(current_cor, symmetric = TRUE)
 eig_target <- eigen(target_cor, symmetric = TRUE)

 if (any(eig_current$values <= 1e-8)) {
  stop("Current correlation matrix is not positive definite.")
 }

 if (any(eig_target$values <= 1e-8)) {
  stop("Target correlation matrix is not positive definite.")
 }

 current_inv_sqrt <- eig_current$vectors %*%
  diag(1 / sqrt(eig_current$values)) %*%
  t(eig_current$vectors)

 target_sqrt <- eig_target$vectors %*%
  diag(sqrt(eig_target$values)) %*%
  t(eig_target$vectors)

 X_new <- X %*% current_inv_sqrt %*% target_sqrt

 X_new
}


# =============================================================================
# Multi-trait simulation with annotation-enriched causal markers
# =============================================================================
#
# Adds support for testing stblr_cpg_omp_csr_prior():
#   - marker annotations / sets, including overlapping sets
#   - enriched causal sampling from selected annotations
#   - marker-specific pi_marker matrix
#   - marker-specific vb_multiplier matrix
#
# Returned objects useful for the new C++ function:
#   sim$annot                 m x K 0/1 annotation matrix
#   sim$sets                  list of marker indices per annotation
#   sim$annotation_id         optional primary annotation per marker
#   sim$pi_marker             list format: length nt, each length m
#   sim$vb_multiplier         list format: length nt, each length m
#   sim$pi_marker_mat         m x nt matrix
#   sim$vb_multiplier_mat     m x nt matrix
#
# Typical call:
#   sim <- mtsim_annotation(
#     Glist = Glist,
#     chr = chr,
#     nt = 3,
#     n_shared = 30,
#     n_specific = 30,
#     n_annotations = 5,
#     enriched_annotations = c(1, 2),
#     annotation_enrichment = 10,
#     base_pi = 0.001,
#     enriched_pi_multiplier = 10,
#     enriched_vb_multiplier = 4,
#     seed = 10
#   )
#
# Then:
#   fit <- stblr_cpg_omp_csr_prior(
#     ...,
#     use_pi_marker = TRUE,
#     pi_marker = sim$pi_marker,
#     use_vb_multiplier = TRUE,
#     vb_multiplier = sim$vb_multiplier,
#     ...
#   )
#
# =============================================================================

make_overlapping_annotations <- function(
  m,
  n_annotations = 5,
  annotation_prob = 0.05,
  force_at_least_one = FALSE,
  seed = NULL
) {
 if (!is.null(seed)) set.seed(seed)

 if (length(annotation_prob) == 1) {
  annotation_prob <- rep(annotation_prob, n_annotations)
 }

 if (length(annotation_prob) != n_annotations) {
  stop("annotation_prob must be scalar or length n_annotations.")
 }

 annot <- matrix(0L, nrow = m, ncol = n_annotations)
 colnames(annot) <- paste0("A", seq_len(n_annotations))

 for (k in seq_len(n_annotations)) {
  annot[, k] <- rbinom(m, size = 1, prob = annotation_prob[k])
 }

 if (force_at_least_one) {
  none <- which(rowSums(annot) == 0)
  if (length(none) > 0) {
   kk <- sample(seq_len(n_annotations), length(none), replace = TRUE)
   annot[cbind(none, kk)] <- 1L
  }
 }

 sets <- lapply(seq_len(n_annotations), function(k) which(annot[, k] == 1L))
 names(sets) <- colnames(annot)

 primary <- rep(NA_integer_, m)
 for (i in seq_len(m)) {
  idx <- which(annot[i, ] == 1L)
  if (length(idx) > 0) primary[i] <- idx[1]
 }

 list(
  annot = annot,
  sets = sets,
  annotation_id = primary
 )
}

make_marker_priors_from_annotations <- function(
  annot,
  nt,
  base_pi = 0.001,
  enriched_annotations = NULL,
  enriched_traits = NULL,
  pi_multiplier = 1,
  vb_multiplier = 1,
  pi_cap = c(1e-8, 0.05),
  vb_cap = c(0.25, 4),
  center_log_pi = TRUE,
  center_log_vb = TRUE
) {
 annot <- as.matrix(annot)
 m <- nrow(annot)
 K <- ncol(annot)

 if (is.null(enriched_annotations)) {
  enriched_annotations <- integer(0)
 }

 if (is.null(enriched_traits)) {
  enriched_traits <- seq_len(nt)
 }

 if (length(base_pi) == 1) base_pi <- rep(base_pi, nt)
 if (length(base_pi) != nt) stop("base_pi must be scalar or length nt.")

 if (length(pi_multiplier) == 1) {
  pi_multiplier <- rep(pi_multiplier, length(enriched_annotations))
 }
 if (length(vb_multiplier) == 1) {
  vb_multiplier <- rep(vb_multiplier, length(enriched_annotations))
 }

 if (length(enriched_annotations) != length(pi_multiplier)) {
  stop("pi_multiplier must be scalar or same length as enriched_annotations.")
 }
 if (length(enriched_annotations) != length(vb_multiplier)) {
  stop("vb_multiplier must be scalar or same length as enriched_annotations.")
 }

 eta_pi <- matrix(0, nrow = K, ncol = nt)
 eta_vb <- matrix(0, nrow = K, ncol = nt)

 if (length(enriched_annotations) > 0) {
  for (a in seq_along(enriched_annotations)) {
   k <- enriched_annotations[a]
   if (k < 1 || k > K) stop("enriched_annotations contains invalid index.")

   eta_pi[k, enriched_traits] <- log(pi_multiplier[a])
   eta_vb[k, enriched_traits] <- log(vb_multiplier[a])
  }
 }

 lin_pi <- annot %*% eta_pi
 lin_vb <- annot %*% eta_vb

 if (center_log_pi) {
  lin_pi <- sweep(lin_pi, 2, colMeans(lin_pi), "-")
 }

 if (center_log_vb) {
  lin_vb <- sweep(lin_vb, 2, colMeans(lin_vb), "-")
 }

 logit_base <- qlogis(pmin(pmax(base_pi, pi_cap[1]), pi_cap[2]))

 pi_marker_mat <- matrix(NA_real_, m, nt)
 vb_multiplier_mat <- matrix(NA_real_, m, nt)

 for (t in seq_len(nt)) {
  pi_marker_mat[, t] <- plogis(logit_base[t] + lin_pi[, t])
  vb_multiplier_mat[, t] <- exp(lin_vb[, t])
 }

 pi_marker_mat <- pmin(pmax(pi_marker_mat, pi_cap[1]), pi_cap[2])
 vb_multiplier_mat <- pmin(pmax(vb_multiplier_mat, vb_cap[1]), vb_cap[2])

 colnames(pi_marker_mat) <- paste0("D", seq_len(nt))
 colnames(vb_multiplier_mat) <- paste0("D", seq_len(nt))

 pi_marker <- lapply(seq_len(nt), function(t) pi_marker_mat[, t])
 vb_multiplier_list <- lapply(seq_len(nt), function(t) vb_multiplier_mat[, t])

 names(pi_marker) <- colnames(pi_marker_mat)
 names(vb_multiplier_list) <- colnames(vb_multiplier_mat)

 list(
  pi_marker = pi_marker,
  vb_multiplier = vb_multiplier_list,
  pi_marker_mat = pi_marker_mat,
  vb_multiplier_mat = vb_multiplier_mat,
  eta_pi = eta_pi,
  eta_vb = eta_vb
 )
}

sample_causal_pool_weighted <- function(
  eligible,
  n_select,
  annot,
  enriched_annotations = NULL,
  annotation_enrichment = 1
) {
 if (n_select == 0) return(integer(0))

 if (length(eligible) < n_select) {
  stop("Not enough eligible markers to sample causal markers.")
 }

 w <- rep(1, length(eligible))

 if (!is.null(enriched_annotations) && length(enriched_annotations) > 0) {
  enriched_hit <- rowSums(annot[eligible, enriched_annotations, drop = FALSE]) > 0
  w[enriched_hit] <- w[enriched_hit] * annotation_enrichment
 }

 w <- w / sum(w)
 sample(eligible, n_select, replace = FALSE, prob = w)
}

mtsim_annotation <- function(
  Glist = NULL,
  chr = NULL,
  rsids = NULL,
  ids = NULL,
  causal_rsids = NULL,
  W = NULL,
  n = 1000,
  m = 1000,
  nt = 2,
  n_shared = 10,
  n_specific = 10,
  h2 = 0.5,
  rg = NULL,
  re = 0.0,
  effect_sd = 1,
  maf_min = 0.05,
  maf_max = 0.5,
  standardize_W = TRUE,
  seed = NULL,
  exact_shared_cor = FALSE,

  # Annotation simulation
  annot = NULL,
  sets = NULL,
  n_annotations = 5,
  annotation_prob = 0.05,
  force_at_least_one_annotation = FALSE,

  # Causal enrichment
  enriched_annotations = NULL,
  annotation_enrichment = 1,

  # Prior construction for stblr_cpg_omp_csr_prior()
  base_pi = 0.001,
  enriched_traits = NULL,
  enriched_pi_multiplier = 1,
  enriched_vb_multiplier = 1,
  pi_cap = c(1e-8, 0.05),
  vb_cap = c(0.25, 4),
  center_log_pi = TRUE,
  center_log_vb = TRUE
) {

 if (!is.null(seed)) {
  set.seed(seed)
 }

 if (length(h2) == 1) {
  h2 <- rep(h2, nt)
 }

 if (length(h2) != nt) {
  stop("h2 must be either a scalar or a vector of length nt.")
 }

 if (is.null(rg)) {
  rg <- diag(nt)
 }

 if (length(rg) == 1 && nt == 2) {
  rg_mat <- matrix(rg, nt, nt)
  diag(rg_mat) <- 1
  rg <- rg_mat
 }

 if (!is.matrix(rg) || !all(dim(rg) == c(nt, nt))) {
  stop("rg must be NULL, a scalar for nt = 2, or an nt x nt correlation matrix.")
 }

 if (any(abs(rg) > 1)) {
  stop("All genetic/effect correlations in rg must be between -1 and 1.")
 }

 if (!isTRUE(all.equal(diag(rg), rep(1, nt)))) {
  stop("The diagonal of rg must be 1.")
 }

 if (any(eigen(rg, symmetric = TRUE)$values <= 1e-8)) {
  stop("rg must be positive definite.")
 }

 # ------------------------------------------------------------
 # Load genotype matrix from Glist
 # ------------------------------------------------------------

 if (!is.null(Glist)) {
  if (is.null(rsids)) {
   if (!is.null(chr) && !is.null(Glist$rsidsLD)) {
    rsids <- Glist$rsidsLD[[chr]]
   } else {
    rsids <- Glist$rsids
   }
  }

  if (is.null(ids)) {
   ids <- Glist$ids
  }

  if (is.null(chr)) {
   W <- qgg::getG(Glist = Glist, rsids = rsids, ids = ids)
  } else {
   W <- qgg::getG(Glist = Glist, chr = chr, rsids = rsids, ids = ids)
  }
 }

 # ------------------------------------------------------------
 # Simulate W if not supplied
 # ------------------------------------------------------------

 if (is.null(W)) {
  maf <- runif(m, maf_min, maf_max)

  W <- sapply(maf, function(p) {
   rbinom(n, size = 2, prob = p)
  })

  W <- as.matrix(W)

  colnames(W) <- paste0("m", seq_len(m))
  rownames(W) <- paste0("id", seq_len(n))
 }

 n <- nrow(W)
 m <- ncol(W)

 if (is.null(colnames(W))) {
  colnames(W) <- paste0("m", seq_len(m))
 }

 if (is.null(rownames(W))) {
  if (!is.null(ids) && length(ids) == n) {
   rownames(W) <- ids
  } else {
   rownames(W) <- paste0("id", seq_len(n))
  }
 }

 # ------------------------------------------------------------
 # Build annotations if not supplied
 # ------------------------------------------------------------

 if (is.null(annot)) {
  ann <- make_overlapping_annotations(
   m = m,
   n_annotations = n_annotations,
   annotation_prob = annotation_prob,
   force_at_least_one = force_at_least_one_annotation,
   seed = NULL
  )

  annot <- ann$annot
  sets <- ann$sets
  annotation_id <- ann$annotation_id
 } else {
  annot <- as.matrix(annot)

  if (nrow(annot) != m) {
   stop("annot must have one row per marker.")
  }

  if (is.null(colnames(annot))) {
   colnames(annot) <- paste0("A", seq_len(ncol(annot)))
  }

  if (is.null(sets)) {
   sets <- lapply(seq_len(ncol(annot)), function(k) which(annot[, k] != 0))
   names(sets) <- colnames(annot)
  }

  annotation_id <- rep(NA_integer_, m)
  for (i in seq_len(m)) {
   idx <- which(annot[i, ] != 0)
   if (length(idx) > 0) annotation_id[i] <- idx[1]
  }
 }

 # ------------------------------------------------------------
 # Define causal marker pool
 # ------------------------------------------------------------

 if (is.null(causal_rsids)) {
  causal_rsids <- colnames(W)
 }

 causal_pool <- which(colnames(W) %in% causal_rsids)

 if (length(causal_pool) == 0) {
  stop("None of the supplied causal_rsids were found in colnames(W).")
 }

 total_needed <- n_shared + nt * n_specific

 if (total_needed > length(causal_pool)) {
  stop(
   "Not enough eligible causal markers. Requested ",
   total_needed,
   " causal markers, but only ",
   length(causal_pool),
   " markers are available in causal_rsids."
  )
 }

 # ------------------------------------------------------------
 # Standardize marker matrix
 # ------------------------------------------------------------

 if (standardize_W) {
  W <- scale(W)
 }

 # ------------------------------------------------------------
 # Select causal markers with optional annotation enrichment
 # ------------------------------------------------------------

 causal_all <- sample_causal_pool_weighted(
  eligible = causal_pool,
  n_select = total_needed,
  annot = annot,
  enriched_annotations = enriched_annotations,
  annotation_enrichment = annotation_enrichment
 )

 shared_idx <- if (n_shared > 0) {
  causal_all[seq_len(n_shared)]
 } else {
  integer(0)
 }

 specific_idx <- vector("list", nt)
 names(specific_idx) <- paste0("D", seq_len(nt))

 start <- n_shared + 1

 for (j in seq_len(nt)) {
  if (n_specific > 0) {
   specific_idx[[j]] <- causal_all[start:(start + n_specific - 1)]
   start <- start + n_specific
  } else {
   specific_idx[[j]] <- integer(0)
  }
 }

 # ------------------------------------------------------------
 # Simulate marker effects directly
 # ------------------------------------------------------------

 B <- matrix(0, nrow = m, ncol = nt)
 colnames(B) <- paste0("D", seq_len(nt))
 rownames(B) <- colnames(W)

 # Shared causal markers: multivariate effects with correlation rg
 if (n_shared > 0) {
  if (nt == 1) {
   B[shared_idx, 1] <- rnorm(n_shared, mean = 0, sd = effect_sd)
  } else {
   if (exact_shared_cor) {
    if (n_shared <= nt) {
     stop("exact_shared_cor = TRUE requires n_shared > nt.")
    }

    B_shared <- matrix(rnorm(n_shared * nt), nrow = n_shared, ncol = nt)
    B_shared <- force_correlation(B_shared, target_cor = rg)
    B_shared <- B_shared * effect_sd
   } else {
    Sigma_b <- effect_sd^2 * rg
    Lb <- chol(Sigma_b)
    B_shared <- matrix(rnorm(n_shared * nt), nrow = n_shared, ncol = nt) %*% Lb
   }

   B[shared_idx, ] <- B_shared
  }
 }

 # Trait-specific causal markers: effects only on one trait
 if (n_specific > 0) {
  for (j in seq_len(nt)) {
   B[specific_idx[[j]], j] <- rnorm(n_specific, mean = 0, sd = effect_sd)
  }
 }

 # ------------------------------------------------------------
 # Genetic values
 # ------------------------------------------------------------

 G_raw <- W %*% B
 colnames(G_raw) <- paste0("D", seq_len(nt))
 rownames(G_raw) <- rownames(W)

 var_g_raw <- apply(G_raw, 2, var)

 if (any(var_g_raw <= 0 | !is.finite(var_g_raw))) {
  stop("At least one trait has zero or non-finite genetic variance.")
 }

 # Scale marker effects so that Var(G_j) = 1 for all traits.
 scale_b <- 1 / sqrt(var_g_raw)

 B <- sweep(B, 2, scale_b, "*")

 G <- W %*% B
 colnames(G) <- paste0("D", seq_len(nt))
 rownames(G) <- rownames(W)

 var_g <- apply(G, 2, var)

 # ------------------------------------------------------------
 # Residual covariance calibrated to target h2
 # ------------------------------------------------------------

 var_e <- var_g * (1 - h2) / h2

 if (length(re) == 1) {
  Sigma_e_cor <- matrix(re, nt, nt)
  diag(Sigma_e_cor) <- 1
 } else {
  Sigma_e_cor <- re
 }

 if (!is.matrix(Sigma_e_cor) || !all(dim(Sigma_e_cor) == c(nt, nt))) {
  stop("re must be a scalar or an nt x nt residual correlation matrix.")
 }

 if (any(abs(Sigma_e_cor) > 1)) {
  stop("All residual correlations in re must be between -1 and 1.")
 }

 if (!isTRUE(all.equal(diag(Sigma_e_cor), rep(1, nt)))) {
  stop("The diagonal of re must be 1.")
 }

 if (any(eigen(Sigma_e_cor, symmetric = TRUE)$values <= 1e-8)) {
  stop("Residual correlation matrix re must be positive definite.")
 }

 Sigma_e <- diag(sqrt(var_e)) %*% Sigma_e_cor %*% diag(sqrt(var_e))
 
 Le <- chol(Sigma_e)

 E <- matrix(rnorm(n * nt), nrow = n, ncol = nt) %*% Le
 colnames(E) <- paste0("D", seq_len(nt))
 rownames(E) <- rownames(W)

 # ------------------------------------------------------------
 # Phenotypes
 # ------------------------------------------------------------

 Y <- G + E
 colnames(Y) <- paste0("D", seq_len(nt))
 rownames(Y) <- rownames(W)

 causal <- list(
  shared = colnames(W)[shared_idx],
  specific = lapply(specific_idx, function(idx) colnames(W)[idx]),
  all = colnames(W)[sort(c(shared_idx, unlist(specific_idx)))]
 )

 B_shared_cor <- if (nt == 1 || n_shared == 0) {
  NA
 } else {
  cor(B[shared_idx, , drop = FALSE])
 }

 B_all_cor <- if (nt == 1) {
  1
 } else {
  cor(B)
 }

 # ------------------------------------------------------------
 # Build marker priors from annotations
 # ------------------------------------------------------------

 pri <- make_marker_priors_from_annotations(
  annot = annot,
  nt = nt,
  base_pi = base_pi,
  enriched_annotations = enriched_annotations,
  enriched_traits = enriched_traits,
  pi_multiplier = enriched_pi_multiplier,
  vb_multiplier = enriched_vb_multiplier,
  pi_cap = pi_cap,
  vb_cap = vb_cap,
  center_log_pi = center_log_pi,
  center_log_vb = center_log_vb
 )

 causal_any <- rep(FALSE, m)
 causal_any[sort(c(shared_idx, unlist(specific_idx)))] <- TRUE

 annotation_causal_summary <- data.frame(
  annotation = colnames(annot),
  size = colSums(annot != 0),
  n_causal = colSums(annot[causal_any, , drop = FALSE] != 0),
  causal_rate = colSums(annot[causal_any, , drop = FALSE] != 0) / pmax(colSums(annot != 0), 1)
 )

 list(
  y = if (nt == 1) Y[, 1] else Y,
  W = W,
  B = B,
  G = G,
  E = E,
  h2_target = h2,
  h2_observed = apply(G, 2, var) / apply(Y, 2, var),
  rg_target = rg,
  rg_observed = if (nt == 1) 1 else cor(G),
  rb_shared_observed = B_shared_cor,
  rb_all_observed = B_all_cor,
  re_target = Sigma_e_cor,
  re_observed = if (nt == 1) 1 else cor(E),
  Sigma_e = Sigma_e,
  causal = causal,
  causal_any = causal_any,
  shared_idx = shared_idx,
  specific_idx = specific_idx,
  rsids = colnames(W),
  ids = rownames(W),
  causal_rsids = causal_rsids,

  annot = annot,
  sets = sets,
  annotation_id = annotation_id,
  enriched_annotations = enriched_annotations,
  annotation_enrichment = annotation_enrichment,
  annotation_causal_summary = annotation_causal_summary,

  pi_marker = pri$pi_marker,
  vb_multiplier = pri$vb_multiplier,
  pi_marker_mat = pri$pi_marker_mat,
  vb_multiplier_mat = pri$vb_multiplier_mat,
  eta_pi = pri$eta_pi,
  eta_vb = pri$eta_vb
 )
}

# =============================================================================
# Diagnostic helper: compare causal enrichment in annotations
# =============================================================================

summarize_annotation_signal <- function(sim, fit = NULL) {
 annot <- sim$annot
 causal_any <- sim$causal_any

 out <- data.frame(
  annotation = colnames(annot),
  size = colSums(annot != 0),
  n_causal = colSums(annot[causal_any, , drop = FALSE] != 0),
  causal_rate_in_set = colSums(annot[causal_any, , drop = FALSE] != 0) / pmax(colSums(annot != 0), 1)
 )

 if (!is.null(fit)) {
  dm <- as.matrix(fit$dm)
  bm <- as.matrix(fit$bm)

  for (t in seq_len(ncol(dm))) {
   out[[paste0("mean_dm_D", t)]] <- colMeans(annot * dm[, t]) / pmax(colMeans(annot), 1e-12)
   out[[paste0("mean_abs_bm_D", t)]] <- colMeans(annot * abs(bm[, t])) / pmax(colMeans(annot), 1e-12)
  }
 }

 out
}
