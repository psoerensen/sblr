# Development notebook containing annotation-aware wrapper prototypes,
# simulations, and benchmarks. It is not intended to run during package checks.
#
# Set SBLR_EXAMPLE_DATA_DIR before running this notebook, or edit the fallback
# path below to point to a local directory containing the example data.
data_dir <- Sys.getenv("SBLR_EXAMPLE_DATA_DIR", unset = "path/to/example/data")
#
# =============================================================================
# Generic R interface for annotation-aware CSR STBLR samplers
# =============================================================================
#
# Supported models:
#   model = "prior"    : fixed marker-specific pi_marker and/or vb_multiplier
#   model = "annot"    : learned continuous/binary annotation effects on pi/vb
#   model = "sbayesrc" : SBayesRC-style mixture components with annotations
#   model = "group"    : group-level annotation priors for pi and variance
#
# Expected C++ functions:
#   stblr_cpg_omp_csr_prior()
#   stblr_cpg_omp_csr_annot()
#   stblr_cpg_omp_csr_sbayesrc()
#   stblr_cpg_omp_csr_group_annot()
#
# Expected generic CSR formatter:
#   format_stblr_fit()
#
# Expected SBayesRC formatter/helpers from stblr_csr_sbayesrc_wrappers_R:
#   format_sbayesrc_csr_fit()
#   make_sbayesrc_component_prior()
#   make_sbayesrc_alpha_init()
#
# =============================================================================

.stblr_match_annotation_model <- function(model) {
 model <- match.arg(model, c("prior", "annot", "sbayesrc", "group"))
 model
}

.stblr_get_nt_m_names <- function(stats, n = NULL, m = NULL) {
 nt <- length(stats$yy)

 if (is.null(n)) {
  if (!is.null(stats$n)) n <- stats$n
  else stop("n must be supplied or available as stats$n.")
 }

 if (is.null(m)) {
  if (!is.null(stats$m)) m <- stats$m
  else m <- length(stats$ww[[1]])
 }

 trait_names <- names(stats$yy)
 if (is.null(trait_names)) trait_names <- paste0("T", seq_len(nt))

 variable_names <- names(stats$ww[[1]])
 if (is.null(variable_names)) variable_names <- paste0("V", seq_len(m))

 list(
  nt = nt,
  n = n,
  m = m,
  trait_names = trait_names,
  variable_names = variable_names
 )
}

.stblr_validate_stats <- function(stats, nt, m) {
 req <- c("wy", "ww", "yy")
 miss <- setdiff(req, names(stats))
 if (length(miss) > 0) {
  stop("stats is missing required element(s): ", paste(miss, collapse = ", "))
 }

 if (length(stats$wy) != nt || length(stats$ww) != nt || length(stats$yy) != nt) {
  stop("stats$wy, stats$ww and stats$yy must all have length nt.")
 }

 if (any(lengths(stats$wy) != m) || any(lengths(stats$ww) != m)) {
  stop("Each stats$wy[[t]] and stats$ww[[t]] must have length m.")
 }

 invisible(TRUE)
}

.stblr_resolve_architecture <- function(
  pi_marker = 0.001,
  pi_init = NULL,
  pi_vb_init = NULL,
  pi_prior_mean = NULL,
  pi_prior_strength = NULL,
  pi_prior_a = NULL,
  pi_prior_b = NULL
) {
 if (is.null(pi_init)) pi_init <- pi_marker
 if (is.null(pi_vb_init)) pi_vb_init <- pi_init
 if (is.null(pi_prior_mean)) pi_prior_mean <- pi_init

 for (nm in c("pi_init", "pi_vb_init", "pi_prior_mean")) {
  val <- get(nm)
  if (!is.numeric(val) || length(val) != 1 || !is.finite(val) || val <= 0 || val >= 1) {
   stop(nm, " must be a finite scalar in (0, 1).")
  }
 }

 if (is.null(pi_prior_a) || is.null(pi_prior_b)) {
  if (is.null(pi_prior_strength)) {
   pi_prior_strength <- 2
  }

  if (!is.numeric(pi_prior_strength) || length(pi_prior_strength) != 1 ||
      !is.finite(pi_prior_strength) || pi_prior_strength <= 0) {
   stop("pi_prior_strength must be positive when pi_prior_a/pi_prior_b are not supplied.")
  }

  pi_prior_a <- pi_prior_mean * pi_prior_strength
  pi_prior_b <- (1 - pi_prior_mean) * pi_prior_strength
 }

 if (!is.numeric(pi_prior_a) || length(pi_prior_a) != 1 ||
     !is.finite(pi_prior_a) || pi_prior_a <= 0) {
  stop("pi_prior_a must be a positive finite scalar.")
 }

 if (!is.numeric(pi_prior_b) || length(pi_prior_b) != 1 ||
     !is.finite(pi_prior_b) || pi_prior_b <= 0) {
  stop("pi_prior_b must be a positive finite scalar.")
 }

 list(
  pi_marker = pi_marker,
  pi_init = pi_init,
  pi_vb_init = pi_vb_init,
  pi_prior_mean = pi_prior_mean,
  pi_prior_strength = pi_prior_strength,
  pi_prior_a = pi_prior_a,
  pi_prior_b = pi_prior_b,
  pi = c(1 - pi_init, pi_init)
 )
}

.stblr_make_csr_variance_priors <- function(
  stats,
  n,
  m,
  nt,
  h2 = 0.5,
  nub = 4,
  nue = 4,
  pi_vb_init = 0.001,
  pi_prior_mean = 0.001,
  trait_names = NULL,
  B = NULL,
  E = NULL,
  ssb_prior = NULL,
  sse_prior = NULL
) {
 if (is.null(trait_names)) trait_names <- names(stats$yy)
 if (is.null(trait_names)) trait_names <- paste0("T", seq_len(nt))

 vy <- as.numeric(stats$yy) / (n - 1)

 if (is.null(B)) {
  B <- diag((vy * h2) / (m * pi_vb_init), nrow = nt, ncol = nt)
 }

 if (is.null(E)) {
  E <- diag(vy * (1 - h2), nrow = nt, ncol = nt)
 }

 if (is.null(ssb_prior)) {
  ssb_prior <- diag(
   ((nub - 2) / nub) * (vy * h2) / (m * pi_prior_mean),
   nrow = nt,
   ncol = nt
  )
 }

 if (is.null(sse_prior)) {
  sse_prior <- diag(
   ((nue - 2) / nue) * (vy * (1 - h2)),
   nrow = nt,
   ncol = nt
  )
 }

 if (!all(dim(B) == c(nt, nt))) stop("B must be nt x nt.")
 if (!all(dim(E) == c(nt, nt))) stop("E must be nt x nt.")
 if (!all(dim(ssb_prior) == c(nt, nt))) stop("ssb_prior must be nt x nt.")
 if (!all(dim(sse_prior) == c(nt, nt))) stop("sse_prior must be nt x nt.")

 rownames(B) <- colnames(B) <- trait_names
 rownames(E) <- colnames(E) <- trait_names
 rownames(ssb_prior) <- colnames(ssb_prior) <- trait_names
 rownames(sse_prior) <- colnames(sse_prior) <- trait_names

 list(
  vy = vy,
  B = B,
  E = E,
  ssb_prior = ssb_prior,
  sse_prior = sse_prior,
  ssb_prior_list = split(ssb_prior, rep(seq_len(ncol(ssb_prior)), each = nrow(ssb_prior))),
  sse_prior_list = split(sse_prior, rep(seq_len(ncol(sse_prior)), each = nrow(sse_prior)))
 )
}

.stblr_init_marker_state <- function(nt, m, b_init = NULL, d_init = NULL) {
 if (is.null(b_init)) b_init <- lapply(seq_len(nt), function(i) rep(0, m))
 if (is.null(d_init)) d_init <- lapply(seq_len(nt), function(i) rep(0, m))

 if (length(b_init) != nt || any(lengths(b_init) != m)) {
  stop("b_init must be a list of length nt, each element length m.")
 }

 if (length(d_init) != nt || any(lengths(d_init) != m)) {
  stop("d_init must be a list of length nt, each element length m.")
 }

 list(b_init = b_init, d_init = d_init)
}

.stblr_init_r_state <- function(stats, nt, m, use_r_init = FALSE, r_init = NULL) {
 if (is.null(r_init)) r_init <- stats$wy

 if (use_r_init) {
  if (length(r_init) != nt || any(lengths(r_init) != m)) {
   stop("r_init must be a list of length nt, each element length m, when use_r_init = TRUE.")
  }
 }

 r_init
}

.stblr_prepare_annotation_matrix <- function(
  A = NULL,
  m,
  variable_names = NULL,
  add_intercept = FALSE,
  standardize = TRUE,
  center_binary = FALSE,
  intercept_name = "Intercept"
) {
 if (is.null(A)) {
  if (!add_intercept) {
   A <- matrix(numeric(0), nrow = m, ncol = 0)
   rownames(A) <- variable_names
   return(A)
  }

  A <- matrix(1, nrow = m, ncol = 1)
  colnames(A) <- intercept_name
  rownames(A) <- variable_names
  return(A)
 }

 A <- as.matrix(A)
 storage.mode(A) <- "double"

 if (!is.null(rownames(A)) && !is.null(variable_names) &&
     all(variable_names %in% rownames(A))) {
  A <- A[variable_names, , drop = FALSE]
 }

 if (nrow(A) != m) {
  stop("A must have one row per marker, or rownames(A) must contain all variable_names.")
 }

 if (any(!is.finite(A))) stop("A contains non-finite values.")

 if (is.null(colnames(A))) colnames(A) <- paste0("Anno", seq_len(ncol(A)))

 has_intercept <- ncol(A) >= 1 && all(abs(A[, 1] - 1) < 1e-12)

 if (standardize && ncol(A) > 0) {
  for (j in seq_len(ncol(A))) {
   is_intercept <- all(abs(A[, j] - 1) < 1e-12)
   is_binary <- all(A[, j] %in% c(0, 1))

   if (is_intercept) next
   if (is_binary && !center_binary) next

   s <- stats::sd(A[, j])
   if (is.finite(s) && s > 0) {
    A[, j] <- (A[, j] - mean(A[, j])) / s
   }
  }
 }

 if (add_intercept && !has_intercept) {
  A <- cbind(Intercept = 1, A)
 }

 if (!is.null(variable_names)) rownames(A) <- variable_names

 A
}

.stblr_make_prior_from_annotations <- function(
  A,
  nt,
  pi_base,
  beta_pi = NULL,
  beta_vb = NULL,
  pi_min = 1e-8,
  pi_max = 0.5,
  vb_multiplier_min = 1e-3,
  vb_multiplier_max = 1e3
) {
 m <- nrow(A)
 K <- ncol(A)

 if (is.null(beta_pi)) beta_pi <- rep(0, K)
 if (is.null(beta_vb)) beta_vb <- rep(0, K)

 beta_pi <- as.matrix(beta_pi)
 beta_vb <- as.matrix(beta_vb)

 if (nrow(beta_pi) == 1 && K > 1) beta_pi <- t(beta_pi)
 if (nrow(beta_vb) == 1 && K > 1) beta_vb <- t(beta_vb)

 if (ncol(beta_pi) == 1 && nt > 1) beta_pi <- beta_pi[, rep(1, nt), drop = FALSE]
 if (ncol(beta_vb) == 1 && nt > 1) beta_vb <- beta_vb[, rep(1, nt), drop = FALSE]

 if (!all(dim(beta_pi) == c(K, nt))) stop("beta_pi must be K x nt or length K.")
 if (!all(dim(beta_vb) == c(K, nt))) stop("beta_vb must be K x nt or length K.")

 pi_marker <- vector("list", nt)
 vb_multiplier <- vector("list", nt)

 for (t in seq_len(nt)) {
  lp_pi <- as.numeric(A %*% beta_pi[, t])
  lp_pi <- lp_pi - mean(lp_pi)

  p <- stats::plogis(stats::qlogis(pi_base) + lp_pi)
  p <- pmin(pmax(p, pi_min), pi_max)

  lp_vb <- as.numeric(A %*% beta_vb[, t])
  lp_vb <- lp_vb - mean(lp_vb)

  mult <- exp(lp_vb)
  mult <- pmin(pmax(mult, vb_multiplier_min), vb_multiplier_max)

  pi_marker[[t]] <- p
  vb_multiplier[[t]] <- mult
 }

 list(
  pi_marker = pi_marker,
  vb_multiplier = vb_multiplier
 )
}

format_csr_annot_fit <- function(
  fit,
  nt,
  m,
  annotation_names,
  trait_names = NULL,
  variable_names = NULL
) {
 if (length(fit) < 22) {
  stop("format_csr_annot_fit() expects the updated 22-slot annotation CSR return object.")
 }

 names(fit)[1:22] <- c(
  "bm", "dm", "wy", "r", "b", "d", "o",
  "vbs", "vgs", "ves",
  "covb", "covg", "cove",
  "vb", "vg", "ve",
  "pi", "pim", "eta_pi", "eta_vb",
  "vle", "vld"
 )

 if (is.null(trait_names)) trait_names <- paste0("T", seq_len(nt))
 if (is.null(variable_names)) variable_names <- paste0("V", seq_len(m))

 K <- length(annotation_names)
 if (K == 0) annotation_names <- "none"

 for (i in 1:7) {
  fit[[i]] <- as.matrix(as.data.frame(fit[[i]]))
  rownames(fit[[i]]) <- variable_names
  colnames(fit[[i]]) <- trait_names
 }

 for (i in 8:10) {
  fit[[i]] <- as.matrix(as.data.frame(fit[[i]]))
  rownames(fit[[i]]) <- paste0("Iter", seq_len(nrow(fit[[i]])))
  colnames(fit[[i]]) <- trait_names
 }

 for (i in 11:16) {
  fit[[i]] <- matrix(unlist(fit[[i]]), ncol = nt, byrow = TRUE)
  rownames(fit[[i]]) <- trait_names
  colnames(fit[[i]]) <- trait_names
 }

 fit$pi <- matrix(unlist(fit$pi), ncol = 2, byrow = TRUE)
 rownames(fit$pi) <- trait_names
 colnames(fit$pi) <- c("pi0", "pi1")

 fit$pim <- matrix(unlist(fit$pim), ncol = 2, byrow = TRUE)
 rownames(fit$pim) <- trait_names
 colnames(fit$pim) <- c("pi0", "pi1")

 eta_pi <- matrix(unlist(fit$eta_pi), nrow = nt, byrow = TRUE)
 eta_vb <- matrix(unlist(fit$eta_vb), nrow = nt, byrow = TRUE)
 rownames(eta_pi) <- rownames(eta_vb) <- trait_names
 colnames(eta_pi) <- colnames(eta_vb) <- annotation_names

 vle <- as.matrix(as.data.frame(fit$vle))
 vld <- as.matrix(as.data.frame(fit$vld))
 rownames(vle) <- paste0("Iter", seq_len(nrow(vle)))
 rownames(vld) <- paste0("Iter", seq_len(nrow(vld)))
 colnames(vle) <- trait_names
 colnames(vld) <- trait_names

 out <- fit[1:18]
 out$eta_pi <- eta_pi
 out$eta_vb <- eta_vb
 out$vle <- vle
 out$vld <- vld

 if (sum(diag(out$covb)) > 0) out$rb <- cov2cor(out$covb)
 if (sum(diag(out$covg)) > 0) out$rg <- cov2cor(out$covg)
 if (sum(diag(out$cove)) > 0) out$re <- cov2cor(out$cove)

 out
}

format_csr_prior_fit <- function(
  fit,
  nt,
  m,
  trait_names = NULL,
  variable_names = NULL
) {
 format_stblr_fit(
  fit = fit,
  nt = nt,
  m = m,
  trait_names = trait_names,
  variable_names = variable_names
 )
}

.stblr_prepare_group_index <- function(
  group = NULL,
  m,
  variable_names = NULL,
  group_names = NULL
) {
 if (is.null(group)) {
  group <- rep("all", m)
 }

 if (is.data.frame(group)) {
  if (ncol(group) != 1) {
   stop("group data frame must have exactly one column.")
  }
  group <- group[[1]]
 }

 if (!is.null(names(group)) && !is.null(variable_names) &&
     all(variable_names %in% names(group))) {
  group <- group[variable_names]
 }

 if (length(group) != m) {
  stop("group must have length m, or names(group) must contain all variable_names.")
 }

 if (any(is.na(group))) {
  stop("group contains NA values.")
 }

 if (is.null(group_names)) {
  if (is.factor(group)) {
   group_names <- levels(group)
  } else {
   group_names <- unique(as.character(group))
  }
 }

 group_factor <- factor(as.character(group), levels = as.character(group_names))
 if (any(is.na(group_factor))) {
  stop("Some group values are not present in group_names.")
 }

 group_index <- as.integer(group_factor) - 1L
 ngroup <- length(levels(group_factor))
 group_names <- levels(group_factor)
 group_size <- as.integer(tabulate(group_index + 1L, nbins = ngroup))

 if (any(group_size <= 0)) {
  stop("Every group must contain at least one marker.")
 }

 list(
  group = group_factor,
  group_index0 = group_index,
  ngroup = ngroup,
  group_names = group_names,
  group_size = group_size
 )
}

.stblr_expand_group_trait_matrix <- function(
  x,
  ngroup,
  nt,
  default,
  name,
  group_names = NULL,
  trait_names = NULL
) {
 if (is.null(x)) {
  out <- matrix(default, nrow = nt, ncol = ngroup)
 } else if (is.list(x) && !is.data.frame(x)) {
  if (length(x) != nt) stop(name, " list must have length nt.")
  out <- do.call(rbind, lapply(x, as.numeric))
 } else {
  x <- as.matrix(x)
  storage.mode(x) <- "double"

  if (length(x) == ngroup && nrow(x) == ngroup && ncol(x) == 1) {
   out <- matrix(as.numeric(x), nrow = nt, ncol = ngroup, byrow = TRUE)
  } else if (length(x) == ngroup && nrow(x) == 1 && ncol(x) == ngroup) {
   out <- matrix(as.numeric(x), nrow = nt, ncol = ngroup, byrow = TRUE)
  } else if (all(dim(x) == c(nt, ngroup))) {
   out <- x
  } else if (all(dim(x) == c(ngroup, nt))) {
   out <- t(x)
  } else {
   stop(name, " must be length ngroup, nt x ngroup, ngroup x nt, or a list of length nt.")
  }
 }

 if (!all(dim(out) == c(nt, ngroup))) {
  stop(name, " could not be expanded to nt x ngroup.")
 }

 if (any(!is.finite(out))) {
  stop(name, " contains non-finite values.")
 }

 rownames(out) <- trait_names %||% paste0("T", seq_len(nt))
 colnames(out) <- group_names %||% paste0("G", seq_len(ngroup))
 out
}

.stblr_make_group_priors <- function(
  ngroup,
  nt,
  pi_init,
  group_pi_init = NULL,
  group_vb_multiplier_init = NULL,
  pi_group_prior_mean = NULL,
  pi_group_prior_strength = NULL,
  pi_group_prior_a = NULL,
  pi_group_prior_b = NULL,
  group_names = NULL,
  trait_names = NULL
) {
 if (is.null(pi_group_prior_mean)) pi_group_prior_mean <- pi_init
 if (length(pi_group_prior_mean) == 1) {
  pi_group_prior_mean <- rep(pi_group_prior_mean, ngroup)
 }
 if (length(pi_group_prior_mean) != ngroup ||
     any(!is.finite(pi_group_prior_mean)) ||
     any(pi_group_prior_mean <= 0 | pi_group_prior_mean >= 1)) {
  stop("pi_group_prior_mean must be a scalar or length-ngroup vector in (0, 1).")
 }

 if (is.null(pi_group_prior_strength)) pi_group_prior_strength <- 2
 if (length(pi_group_prior_strength) == 1) {
  pi_group_prior_strength <- rep(pi_group_prior_strength, ngroup)
 }
 if (length(pi_group_prior_strength) != ngroup ||
     any(!is.finite(pi_group_prior_strength)) ||
     any(pi_group_prior_strength <= 0)) {
  stop("pi_group_prior_strength must be positive scalar or length-ngroup vector.")
 }

 if (is.null(pi_group_prior_a) || is.null(pi_group_prior_b)) {
  pi_group_prior_a <- pi_group_prior_mean * pi_group_prior_strength
  pi_group_prior_b <- (1 - pi_group_prior_mean) * pi_group_prior_strength
 }

 pi_group_prior_a <- as.numeric(pi_group_prior_a)
 pi_group_prior_b <- as.numeric(pi_group_prior_b)

 if (length(pi_group_prior_a) == 1) pi_group_prior_a <- rep(pi_group_prior_a, ngroup)
 if (length(pi_group_prior_b) == 1) pi_group_prior_b <- rep(pi_group_prior_b, ngroup)

 if (length(pi_group_prior_a) != ngroup || length(pi_group_prior_b) != ngroup ||
     any(!is.finite(pi_group_prior_a)) || any(!is.finite(pi_group_prior_b)) ||
     any(pi_group_prior_a <= 0) || any(pi_group_prior_b <= 0)) {
  stop("pi_group_prior_a and pi_group_prior_b must be positive scalars or length-ngroup vectors.")
 }

 group_pi_init_mat <- .stblr_expand_group_trait_matrix(
  group_pi_init,
  ngroup = ngroup,
  nt = nt,
  default = pi_init,
  name = "group_pi_init",
  group_names = group_names,
  trait_names = trait_names
 )

 if (any(group_pi_init_mat <= 0 | group_pi_init_mat >= 1)) {
  stop("group_pi_init values must be in (0, 1).")
 }

 group_vb_multiplier_init_mat <- .stblr_expand_group_trait_matrix(
  group_vb_multiplier_init,
  ngroup = ngroup,
  nt = nt,
  default = 1,
  name = "group_vb_multiplier_init",
  group_names = group_names,
  trait_names = trait_names
 )

 if (any(group_vb_multiplier_init_mat <= 0)) {
  stop("group_vb_multiplier_init values must be positive.")
 }

 list(
  group_pi_init = split(group_pi_init_mat, row(group_pi_init_mat)),
  group_vb_multiplier_init = split(group_vb_multiplier_init_mat, row(group_vb_multiplier_init_mat)),
  group_pi_init_matrix = group_pi_init_mat,
  group_vb_multiplier_init_matrix = group_vb_multiplier_init_mat,
  pi_group_prior_a = pi_group_prior_a,
  pi_group_prior_b = pi_group_prior_b,
  pi_group_prior_mean = pi_group_prior_mean,
  pi_group_prior_strength = pi_group_prior_strength
 )
}

format_csr_group_annot_fit <- function(
  fit,
  nt,
  m,
  ngroup,
  group_names = NULL,
  trait_names = NULL,
  variable_names = NULL
) {
 if (length(fit) < 26) {
  stop("format_csr_group_annot_fit() expects the updated 26-slot group-annotation CSR return object.")
 }

 names(fit)[1:26] <- c(
  "bm", "dm", "wy", "r", "b", "d", "o",
  "vbs", "vgs", "ves",
  "covb", "covg", "cove",
  "vb", "vg", "ve",
  "pi", "pim", "diagnostics", "pimarker",
  "vle", "vld",
  "group_pi", "group_vb_multiplier", "group_nincluded", "group_size"
 )

 if (is.null(trait_names)) trait_names <- paste0("T", seq_len(nt))
 if (is.null(variable_names)) variable_names <- paste0("V", seq_len(m))
 if (is.null(group_names)) group_names <- paste0("G", seq_len(ngroup))

 for (i in 1:7) {
  fit[[i]] <- as.matrix(as.data.frame(fit[[i]]))
  rownames(fit[[i]]) <- variable_names
  colnames(fit[[i]]) <- trait_names
 }

 for (i in 8:10) {
  fit[[i]] <- as.matrix(as.data.frame(fit[[i]]))
  rownames(fit[[i]]) <- paste0("Iter", seq_len(nrow(fit[[i]])))
  colnames(fit[[i]]) <- trait_names
 }

 for (i in 11:16) {
  fit[[i]] <- matrix(unlist(fit[[i]]), ncol = nt, byrow = TRUE)
  rownames(fit[[i]]) <- trait_names
  colnames(fit[[i]]) <- trait_names
 }

 fit$pi <- matrix(unlist(fit$pi), ncol = 2, byrow = TRUE)
 rownames(fit$pi) <- trait_names
 colnames(fit$pi) <- c("pi0", "pi1")

 fit$pim <- matrix(unlist(fit$pim), ncol = 2, byrow = TRUE)
 rownames(fit$pim) <- trait_names
 colnames(fit$pim) <- c("pi0", "pi1")

 fit$diagnostics <- matrix(unlist(fit$diagnostics), ncol = 4, byrow = TRUE)
 rownames(fit$diagnostics) <- trait_names
 colnames(fit$diagnostics) <- c("log_cpo", "mean_log_cpo", "seconds_mean", "seconds_max")

 fit$pimarker <- matrix(unlist(fit$pimarker), ncol = 2, byrow = TRUE)
 rownames(fit$pimarker) <- trait_names
 colnames(fit$pimarker) <- c("nsamples", "n")

 vle <- as.matrix(as.data.frame(fit$vle))
 vld <- as.matrix(as.data.frame(fit$vld))
 rownames(vle) <- paste0("Iter", seq_len(nrow(vle)))
 rownames(vld) <- paste0("Iter", seq_len(nrow(vld)))
 colnames(vle) <- trait_names
 colnames(vld) <- trait_names

 group_pi <- matrix(unlist(fit$group_pi), nrow = nt, ncol = ngroup, byrow = TRUE)
 group_vb_multiplier <- matrix(unlist(fit$group_vb_multiplier), nrow = nt, ncol = ngroup, byrow = TRUE)
 group_nincluded <- matrix(unlist(fit$group_nincluded), nrow = nt, ncol = ngroup, byrow = TRUE)
 group_size <- matrix(unlist(fit$group_size), nrow = nt, ncol = ngroup, byrow = TRUE)

 rownames(group_pi) <- rownames(group_vb_multiplier) <- rownames(group_nincluded) <- rownames(group_size) <- trait_names
 colnames(group_pi) <- colnames(group_vb_multiplier) <- colnames(group_nincluded) <- colnames(group_size) <- group_names

 out <- fit[1:20]
 out$vle <- vle
 out$vld <- vld
 out$group_pi <- group_pi
 out$group_vb_multiplier <- group_vb_multiplier
 out$group_nincluded <- group_nincluded
 out$group_size <- group_size
 out$log_cpo <- out$diagnostics[, "log_cpo"]
 out$mean_log_cpo <- out$diagnostics[, "mean_log_cpo"]

 if (sum(diag(out$covb)) > 0) out$rb <- cov2cor(out$covb)
 if (sum(diag(out$covg)) > 0) out$rg <- cov2cor(out$covg)
 if (sum(diag(out$cove)) > 0) out$re <- cov2cor(out$cove)

 out
}

`%||%` <- function(x, y) {
 if (is.null(x)) y else x
}

# =============================================================================
# Main generic interface
# =============================================================================

stblr_csr_annotation <- function(
  stats,
  ld_prefix,
  model = c("prior", "annot", "sbayesrc", "group"),
  A = NULL,
  n = NULL,
  m = NULL,

  # Shared architecture controls
  pi_marker = 0.001,
  pi_init = NULL,
  pi_vb_init = NULL,
  pi_prior_mean = NULL,
  pi_prior_strength = NULL,
  pi_prior_a = NULL,
  pi_prior_b = NULL,

  h2 = 0.5,
  nub = 4,
  nue = 4,
  B = NULL,
  E = NULL,
  ssb_prior = NULL,
  sse_prior = NULL,

  updateB = TRUE,
  updateE = TRUE,
  updatePi = TRUE,
  adjE = 0.9,

  nit = 1000,
  nburn = 100,
  nthin = 1,
  ncores = 3,
  seed = 10,

  b_init = NULL,
  d_init = NULL,
  use_d_init = FALSE,
  r_init = NULL,
  use_r_init = FALSE,
  rebuild_r_before_updateE = FALSE,

  # Annotation preprocessing
  add_intercept = FALSE,
  standardize_annotations = TRUE,
  center_binary_annotations = FALSE,

  # model = "prior" options
  use_pi_marker = FALSE,
  use_vb_multiplier = FALSE,
  fixed_pi_marker = NULL,
  fixed_vb_multiplier = NULL,
  beta_pi = NULL,
  beta_vb = NULL,

  # model = "annot" options
  learn_pi_annot = TRUE,
  learn_vb_annot = FALSE,
  eta_pi_init = NULL,
  eta_vb_init = NULL,
  sigma_eta_pi = 1,
  sigma_eta_vb = 1,
  rw_sd_eta_pi = 0.05,
  rw_sd_eta_vb = 0.05,
  annot_update_every = 10,
  pi_min = 1e-8,
  pi_max = 0.5,
  vb_multiplier_min = 1e-3,
  vb_multiplier_max = 1e3,

  # model = "group" options
  group = NULL,
  group_names = NULL,
  group_pi_init = NULL,
  group_vb_multiplier_init = NULL,
  pi_group_prior_mean = NULL,
  pi_group_prior_strength = NULL,
  pi_group_prior_a = NULL,
  pi_group_prior_b = NULL,
  updateGroupVb = FALSE,
  nub_group = 4,
  ssb_group_prior = 1,
  normalize_group_vb = TRUE,

  # model = "sbayesrc" options
  gamma = c(0, 0.01, 0.1, 1),
  active_comp_weights = NULL,
  alpha_init = NULL,
  sigmaSqAlpha_init = NULL,
  intercept_flat = TRUE,
  sigmaSqAlpha_a = 2,
  sigmaSqAlpha_b = 2,
  pi_floor = 1e-12,
  updateAlpha = TRUE,
  alpha_update_every = 10,
  comp_init = NULL,
  use_comp_init = FALSE
) {
 model <- .stblr_match_annotation_model(model)

 dims <- .stblr_get_nt_m_names(stats, n = n, m = m)
 nt <- dims$nt
 n <- dims$n
 m <- dims$m
 trait_names <- dims$trait_names
 variable_names <- dims$variable_names

 .stblr_validate_stats(stats, nt = nt, m = m)

 arch <- .stblr_resolve_architecture(
  pi_marker = pi_marker,
  pi_init = pi_init,
  pi_vb_init = pi_vb_init,
  pi_prior_mean = pi_prior_mean,
  pi_prior_strength = pi_prior_strength,
  pi_prior_a = pi_prior_a,
  pi_prior_b = pi_prior_b
 )

 pri <- .stblr_make_csr_variance_priors(
  stats = stats,
  n = n,
  m = m,
  nt = nt,
  h2 = h2,
  nub = nub,
  nue = nue,
  pi_vb_init = arch$pi_vb_init,
  pi_prior_mean = arch$pi_prior_mean,
  trait_names = trait_names,
  B = B,
  E = E,
  ssb_prior = ssb_prior,
  sse_prior = sse_prior
 )

 state <- .stblr_init_marker_state(nt = nt, m = m, b_init = b_init, d_init = d_init)
 r_init <- .stblr_init_r_state(stats, nt = nt, m = m, use_r_init = use_r_init, r_init = r_init)

 if (model %in% c("prior", "annot", "sbayesrc")) {
  A <- .stblr_prepare_annotation_matrix(
   A = A,
   m = m,
   variable_names = variable_names,
   add_intercept = add_intercept || model == "sbayesrc",
   standardize = standardize_annotations,
   center_binary = center_binary_annotations
  )
 }

 annotation_names <- colnames(A)

 if (model == "prior") {
  if (!is.null(fixed_pi_marker)) {
   use_pi_marker <- TRUE
   pi_marker_list <- fixed_pi_marker
  } else if (use_pi_marker || !is.null(beta_pi)) {
   ann_prior <- .stblr_make_prior_from_annotations(
    A = A,
    nt = nt,
    pi_base = arch$pi_init,
    beta_pi = beta_pi,
    beta_vb = beta_vb,
    pi_min = pi_min,
    pi_max = pi_max,
    vb_multiplier_min = vb_multiplier_min,
    vb_multiplier_max = vb_multiplier_max
   )
   pi_marker_list <- ann_prior$pi_marker
  } else {
   pi_marker_list <- lapply(seq_len(nt), function(i) rep(arch$pi_init, m))
  }

  if (!is.null(fixed_vb_multiplier)) {
   use_vb_multiplier <- TRUE
   vb_multiplier_list <- fixed_vb_multiplier
  } else if (use_vb_multiplier || !is.null(beta_vb)) {
   if (!exists("ann_prior", inherits = FALSE)) {
    ann_prior <- .stblr_make_prior_from_annotations(
     A = A,
     nt = nt,
     pi_base = arch$pi_init,
     beta_pi = beta_pi,
     beta_vb = beta_vb,
     pi_min = pi_min,
     pi_max = pi_max,
     vb_multiplier_min = vb_multiplier_min,
     vb_multiplier_max = vb_multiplier_max
    )
   }
   vb_multiplier_list <- ann_prior$vb_multiplier
  } else {
   vb_multiplier_list <- lapply(seq_len(nt), function(i) rep(1, m))
  }

  raw_fit <- stblr_cpg_omp_csr_prior(
   wy = stats$wy,
   ww = stats$ww,
   yy = stats$yy,
   b_init = state$b_init,
   d_init = state$d_init,
   use_d_init = use_d_init,
   r_init = r_init,
   use_r_init = use_r_init,
   rebuild_r_before_updateE = rebuild_r_before_updateE,
   ld_prefix = ld_prefix,
   B = pri$B,
   E = pri$E,
   ssb_prior = pri$ssb_prior_list,
   sse_prior = pri$sse_prior_list,
   pi = arch$pi,
   use_pi_marker = use_pi_marker,
   pi_marker = pi_marker_list,
   use_vb_multiplier = use_vb_multiplier,
   vb_multiplier = vb_multiplier_list,
   nub = nub,
   nue = nue,
   updateB = updateB,
   updateE = updateE,
   updatePi = updatePi,
   adjE = adjE,
   n = rep(as.integer(n), nt),
   nit = as.integer(nit),
   nburn = as.integer(nburn),
   nthin = as.integer(nthin),
   pi_prior_a = arch$pi_prior_a,
   pi_prior_b = arch$pi_prior_b,
   ncores = as.integer(ncores),
   seed = as.integer(seed)
  )

  fit <- format_csr_prior_fit(
   fit = raw_fit,
   nt = nt,
   m = m,
   trait_names = trait_names,
   variable_names = variable_names
  )

  fit$input <- list(
   model = model,
   n = n,
   m = m,
   nt = nt,
   A = A,
   annotation_names = annotation_names,
   use_pi_marker = use_pi_marker,
   use_vb_multiplier = use_vb_multiplier,
   pi_marker = if (use_pi_marker) pi_marker_list else NULL,
   vb_multiplier = if (use_vb_multiplier) vb_multiplier_list else NULL
  )
 }

 if (model == "annot") {
  K <- ncol(A)

  if (K == 0) {
   stop("model = 'annot' requires at least one annotation column in A.")
  }

  if (is.null(eta_pi_init)) eta_pi_init <- matrix(0, nrow = K, ncol = nt)
  if (is.null(eta_vb_init)) eta_vb_init <- matrix(0, nrow = K, ncol = nt)

  eta_pi_init <- as.matrix(eta_pi_init)
  eta_vb_init <- as.matrix(eta_vb_init)
  storage.mode(eta_pi_init) <- "double"
  storage.mode(eta_vb_init) <- "double"

  if (!all(dim(eta_pi_init) == c(K, nt))) stop("eta_pi_init must be ncol(A) x nt.")
  if (!all(dim(eta_vb_init) == c(K, nt))) stop("eta_vb_init must be ncol(A) x nt.")

  raw_fit <- stblr_cpg_omp_csr_annot(
   wy = stats$wy,
   ww = stats$ww,
   yy = stats$yy,
   b_init = state$b_init,
   d_init = state$d_init,
   use_d_init = use_d_init,
   r_init = r_init,
   use_r_init = use_r_init,
   rebuild_r_before_updateE = rebuild_r_before_updateE,
   ld_prefix = ld_prefix,
   B = pri$B,
   E = pri$E,
   ssb_prior = pri$ssb_prior_list,
   sse_prior = pri$sse_prior_list,
   pi = arch$pi,
   A = A,
   learn_pi_annot = learn_pi_annot,
   learn_vb_annot = learn_vb_annot,
   eta_pi_init = eta_pi_init,
   eta_vb_init = eta_vb_init,
   sigma_eta_pi = sigma_eta_pi,
   sigma_eta_vb = sigma_eta_vb,
   rw_sd_eta_pi = rw_sd_eta_pi,
   rw_sd_eta_vb = rw_sd_eta_vb,
   annot_update_every = as.integer(annot_update_every),
   pi_min = pi_min,
   pi_max = pi_max,
   vb_multiplier_min = vb_multiplier_min,
   vb_multiplier_max = vb_multiplier_max,
   nub = nub,
   nue = nue,
   updateB = updateB,
   updateE = updateE,
   updatePi = updatePi,
   adjE = adjE,
   n = rep(as.integer(n), nt),
   nit = as.integer(nit),
   nburn = as.integer(nburn),
   nthin = as.integer(nthin),
   pi_prior_a = arch$pi_prior_a,
   pi_prior_b = arch$pi_prior_b,
   ncores = as.integer(ncores),
   seed = as.integer(seed)
  )

  fit <- format_csr_annot_fit(
   fit = raw_fit,
   nt = nt,
   m = m,
   annotation_names = annotation_names,
   trait_names = trait_names,
   variable_names = variable_names
  )

  fit$input <- list(
   model = model,
   n = n,
   m = m,
   nt = nt,
   A = A,
   annotation_names = annotation_names,
   learn_pi_annot = learn_pi_annot,
   learn_vb_annot = learn_vb_annot,
   eta_pi_init = eta_pi_init,
   eta_vb_init = eta_vb_init,
   sigma_eta_pi = sigma_eta_pi,
   sigma_eta_vb = sigma_eta_vb,
   rw_sd_eta_pi = rw_sd_eta_pi,
   rw_sd_eta_vb = rw_sd_eta_vb,
   annot_update_every = annot_update_every
  )
 }

 if (model == "group") {
  group_info <- .stblr_prepare_group_index(
   group = group,
   m = m,
   variable_names = variable_names,
   group_names = group_names
  )

  group_priors <- .stblr_make_group_priors(
   ngroup = group_info$ngroup,
   nt = nt,
   pi_init = arch$pi_init,
   group_pi_init = group_pi_init,
   group_vb_multiplier_init = group_vb_multiplier_init,
   pi_group_prior_mean = pi_group_prior_mean,
   pi_group_prior_strength = pi_group_prior_strength,
   pi_group_prior_a = pi_group_prior_a,
   pi_group_prior_b = pi_group_prior_b,
   group_names = group_info$group_names,
   trait_names = trait_names
  )

  raw_fit <- stblr_cpg_omp_csr_group_annot(
   wy = stats$wy,
   ww = stats$ww,
   yy = stats$yy,
   b_init = state$b_init,
   d_init = state$d_init,
   use_d_init = use_d_init,
   r_init = r_init,
   use_r_init = use_r_init,
   rebuild_r_before_updateE = rebuild_r_before_updateE,
   ld_prefix = ld_prefix,
   B = pri$B,
   E = pri$E,
   ssb_prior = pri$ssb_prior_list,
   sse_prior = pri$sse_prior_list,
   pi = arch$pi,
   group_index = as.integer(group_info$group_index0),
   ngroup = as.integer(group_info$ngroup),
   group_pi_init = group_priors$group_pi_init,
   pi_group_prior_a = group_priors$pi_group_prior_a,
   pi_group_prior_b = group_priors$pi_group_prior_b,
   group_vb_multiplier_init = group_priors$group_vb_multiplier_init,
   updateGroupVb = updateGroupVb,
   nub_group = nub_group,
   ssb_group_prior = ssb_group_prior,
   normalize_group_vb = normalize_group_vb,
   nub = nub,
   nue = nue,
   updateB = updateB,
   updateE = updateE,
   updatePi = updatePi,
   adjE = adjE,
   n = rep(as.integer(n), nt),
   nit = as.integer(nit),
   nburn = as.integer(nburn),
   nthin = as.integer(nthin),
   ncores = as.integer(ncores),
   seed = as.integer(seed)
  )

  fit <- format_csr_group_annot_fit(
   fit = raw_fit,
   nt = nt,
   m = m,
   ngroup = group_info$ngroup,
   group_names = group_info$group_names,
   trait_names = trait_names,
   variable_names = variable_names
  )

  fit$input <- list(
   model = model,
   n = n,
   m = m,
   nt = nt,
   group = group_info$group,
   group_index0 = group_info$group_index0,
   ngroup = group_info$ngroup,
   group_names = group_info$group_names,
   group_size = group_info$group_size,
   group_pi_init = group_priors$group_pi_init_matrix,
   group_vb_multiplier_init = group_priors$group_vb_multiplier_init_matrix,
   pi_group_prior_mean = group_priors$pi_group_prior_mean,
   pi_group_prior_strength = group_priors$pi_group_prior_strength,
   pi_group_prior_a = group_priors$pi_group_prior_a,
   pi_group_prior_b = group_priors$pi_group_prior_b,
   updateGroupVb = updateGroupVb,
   nub_group = nub_group,
   ssb_group_prior = ssb_group_prior,
   normalize_group_vb = normalize_group_vb
  )
 }

 if (model == "sbayesrc") {
  gamma <- as.numeric(gamma)
  Kgamma <- length(gamma)
  if (Kgamma < 2) stop("gamma must have at least two elements.")
  if (!isTRUE(all.equal(gamma[1], 0))) stop("gamma[1] must be 0.")

  if (is.null(comp_init)) {
   comp_init <- lapply(seq_len(nt), function(i) rep(0, m))
  }
  if (length(comp_init) != nt || any(lengths(comp_init) != m)) {
   stop("comp_init must be a list of length nt, each element length m.")
  }

  alpha <- make_sbayesrc_alpha_init(
   A = A,
   gamma = gamma,
   pi_init = arch$pi_init,
   active_comp_weights = active_comp_weights,
   alpha_init = alpha_init,
   sigmaSqAlpha_init = sigmaSqAlpha_init
  )

  raw_fit <- stblr_cpg_omp_csr_sbayesrc(
   wy = stats$wy,
   ww = stats$ww,
   yy = stats$yy,
   b_init = state$b_init,
   comp_init = comp_init,
   use_comp_init = use_comp_init,
   r_init = r_init,
   use_r_init = use_r_init,
   rebuild_r_before_updateE = rebuild_r_before_updateE,
   ld_prefix = ld_prefix,
   B = pri$B,
   E = pri$E,
   ssb_prior = pri$ssb_prior_list,
   sse_prior = pri$sse_prior_list,
   A = A,
   gamma = gamma,
   alpha_init = alpha$alpha_init,
   sigmaSqAlpha_init = alpha$sigmaSqAlpha_init,
   intercept_flat = intercept_flat,
   sigmaSqAlpha_a = sigmaSqAlpha_a,
   sigmaSqAlpha_b = sigmaSqAlpha_b,
   pi_floor = pi_floor,
   nub = nub,
   nue = nue,
   updateAlpha = updateAlpha,
   updateB = updateB,
   updateE = updateE,
   alpha_update_every = as.integer(alpha_update_every),
   adjE = adjE,
   n = rep(as.integer(n), nt),
   nit = as.integer(nit),
   nburn = as.integer(nburn),
   nthin = as.integer(nthin),
   ncores = as.integer(ncores),
   seed = as.integer(seed)
  )

  fit <- format_sbayesrc_csr_fit(
   fit = raw_fit,
   nt = nt,
   m = m,
   gamma = gamma,
   n_anno = ncol(A),
   trait_names = trait_names,
   variable_names = variable_names,
   annotation_names = annotation_names
  )

  fit$input <- list(
   model = model,
   n = n,
   m = m,
   nt = nt,
   A = A,
   annotation_names = annotation_names,
   gamma = gamma,
   component_names = colnames(fit$ncomp),
   active_comp_weights = alpha$active_comp_weights,
   component_prob_init = alpha$component_prob_init,
   step_prob_init = alpha$step_prob_init,
   alpha_init = alpha$alpha_init,
   sigmaSqAlpha_init = alpha$sigmaSqAlpha_init,
   intercept_flat = intercept_flat,
   sigmaSqAlpha_a = sigmaSqAlpha_a,
   sigmaSqAlpha_b = sigmaSqAlpha_b,
   pi_floor = pi_floor,
   updateAlpha = updateAlpha,
   alpha_update_every = alpha_update_every
  )
 }

 shared_input <- list(
  pi_marker = arch$pi_marker,
  pi_init = arch$pi_init,
  pi_vb_init = arch$pi_vb_init,
  pi_prior_mean = arch$pi_prior_mean,
  pi_prior_strength = arch$pi_prior_strength,
  pi_prior_a = arch$pi_prior_a,
  pi_prior_b = arch$pi_prior_b,
  h2 = h2,
  nub = nub,
  nue = nue,
  vy = pri$vy,
  B = pri$B,
  E = pri$E,
  ssb_prior = pri$ssb_prior,
  sse_prior = pri$sse_prior,
  updateB = updateB,
  updateE = updateE,
  updatePi = updatePi,
  adjE = adjE,
  nit = nit,
  nburn = nburn,
  nthin = nthin,
  ncores = ncores,
  seed = seed,
  use_d_init = use_d_init,
  use_r_init = use_r_init,
  rebuild_r_before_updateE = rebuild_r_before_updateE,
  ld_prefix = ld_prefix
 )

 fit$input <- c(fit$input, shared_input)
 fit
}

# Convenience wrappers ---------------------------------------------------------

stblr_csr_prior_annot <- function(...) {
 stblr_csr_annotation(..., model = "prior")
}

stblr_csr_learn_annot <- function(...) {
 stblr_csr_annotation(..., model = "annot")
}

stblr_csr_sbayesrc_generic <- function(...) {
 stblr_csr_annotation(..., model = "sbayesrc")
}

stblr_csr_group_annot <- function(...) {
 stblr_csr_annotation(..., model = "group")
}

# # =============================================================================
# # Generic R interface for annotation-aware CSR STBLR samplers
# # =============================================================================
# #
# # Supported models:
# #   model = "prior"    : fixed marker-specific pi_marker and/or vb_multiplier
# #   model = "annot"    : learned continuous/binary annotation effects on pi/vb
# #   model = "sbayesrc" : SBayesRC-style mixture components with annotations
# #
# # Expected C++ functions:
# #   stblr_cpg_omp_csr_prior()
# #   stblr_cpg_omp_csr_annot()
# #   stblr_cpg_omp_csr_sbayesrc()
# #
# # Expected generic CSR formatter:
# #   format_stblr_fit()
# #
# # Expected SBayesRC formatter/helpers from stblr_csr_sbayesrc_wrappers_R:
# #   format_sbayesrc_csr_fit()
# #   make_sbayesrc_component_prior()
# #   make_sbayesrc_alpha_init()
# #
# # =============================================================================
#
# .stblr_match_annotation_model <- function(model) {
#  model <- match.arg(model, c("prior", "annot", "sbayesrc"))
#  model
# }
#
# .stblr_get_nt_m_names <- function(stats, n = NULL, m = NULL) {
#  nt <- length(stats$yy)
#
#  if (is.null(n)) {
#   if (!is.null(stats$n)) n <- stats$n
#   else stop("n must be supplied or available as stats$n.")
#  }
#
#  if (is.null(m)) {
#   if (!is.null(stats$m)) m <- stats$m
#   else m <- length(stats$ww[[1]])
#  }
#
#  trait_names <- names(stats$yy)
#  if (is.null(trait_names)) trait_names <- paste0("T", seq_len(nt))
#
#  variable_names <- names(stats$ww[[1]])
#  if (is.null(variable_names)) variable_names <- paste0("V", seq_len(m))
#
#  list(
#   nt = nt,
#   n = n,
#   m = m,
#   trait_names = trait_names,
#   variable_names = variable_names
#  )
# }
#
# .stblr_validate_stats <- function(stats, nt, m) {
#  req <- c("wy", "ww", "yy")
#  miss <- setdiff(req, names(stats))
#  if (length(miss) > 0) {
#   stop("stats is missing required element(s): ", paste(miss, collapse = ", "))
#  }
#
#  if (length(stats$wy) != nt || length(stats$ww) != nt || length(stats$yy) != nt) {
#   stop("stats$wy, stats$ww and stats$yy must all have length nt.")
#  }
#
#  if (any(lengths(stats$wy) != m) || any(lengths(stats$ww) != m)) {
#   stop("Each stats$wy[[t]] and stats$ww[[t]] must have length m.")
#  }
#
#  invisible(TRUE)
# }
#
# .stblr_resolve_architecture <- function(
#   pi_marker = 0.001,
#   pi_init = NULL,
#   pi_vb_init = NULL,
#   pi_prior_mean = NULL,
#   pi_prior_strength = NULL,
#   pi_prior_a = NULL,
#   pi_prior_b = NULL
# ) {
#  if (is.null(pi_init)) pi_init <- pi_marker
#  if (is.null(pi_vb_init)) pi_vb_init <- pi_init
#  if (is.null(pi_prior_mean)) pi_prior_mean <- pi_init
#
#  for (nm in c("pi_init", "pi_vb_init", "pi_prior_mean")) {
#   val <- get(nm)
#   if (!is.numeric(val) || length(val) != 1 || !is.finite(val) || val <= 0 || val >= 1) {
#    stop(nm, " must be a finite scalar in (0, 1).")
#   }
#  }
#
#  if (is.null(pi_prior_a) || is.null(pi_prior_b)) {
#   if (is.null(pi_prior_strength)) {
#    pi_prior_strength <- 2
#   }
#
#   if (!is.numeric(pi_prior_strength) || length(pi_prior_strength) != 1 ||
#       !is.finite(pi_prior_strength) || pi_prior_strength <= 0) {
#    stop("pi_prior_strength must be positive when pi_prior_a/pi_prior_b are not supplied.")
#   }
#
#   pi_prior_a <- pi_prior_mean * pi_prior_strength
#   pi_prior_b <- (1 - pi_prior_mean) * pi_prior_strength
#  }
#
#  if (!is.numeric(pi_prior_a) || length(pi_prior_a) != 1 ||
#      !is.finite(pi_prior_a) || pi_prior_a <= 0) {
#   stop("pi_prior_a must be a positive finite scalar.")
#  }
#
#  if (!is.numeric(pi_prior_b) || length(pi_prior_b) != 1 ||
#      !is.finite(pi_prior_b) || pi_prior_b <= 0) {
#   stop("pi_prior_b must be a positive finite scalar.")
#  }
#
#  list(
#   pi_marker = pi_marker,
#   pi_init = pi_init,
#   pi_vb_init = pi_vb_init,
#   pi_prior_mean = pi_prior_mean,
#   pi_prior_strength = pi_prior_strength,
#   pi_prior_a = pi_prior_a,
#   pi_prior_b = pi_prior_b,
#   pi = c(1 - pi_init, pi_init)
#  )
# }
#
# .stblr_make_csr_variance_priors <- function(
#   stats,
#   n,
#   m,
#   nt,
#   h2 = 0.5,
#   nub = 4,
#   nue = 4,
#   pi_vb_init = 0.001,
#   pi_prior_mean = 0.001,
#   trait_names = NULL,
#   B = NULL,
#   E = NULL,
#   ssb_prior = NULL,
#   sse_prior = NULL
# ) {
#  if (is.null(trait_names)) trait_names <- names(stats$yy)
#  if (is.null(trait_names)) trait_names <- paste0("T", seq_len(nt))
#
#  vy <- as.numeric(stats$yy) / (n - 1)
#
#  if (is.null(B)) {
#   B <- diag((vy * h2) / (m * pi_vb_init), nrow = nt, ncol = nt)
#  }
#
#  if (is.null(E)) {
#   E <- diag(vy * (1 - h2), nrow = nt, ncol = nt)
#  }
#
#  if (is.null(ssb_prior)) {
#   ssb_prior <- diag(
#    ((nub - 2) / nub) * (vy * h2) / (m * pi_prior_mean),
#    nrow = nt,
#    ncol = nt
#   )
#  }
#
#  if (is.null(sse_prior)) {
#   sse_prior <- diag(
#    ((nue - 2) / nue) * (vy * (1 - h2)),
#    nrow = nt,
#    ncol = nt
#   )
#  }
#
#  if (!all(dim(B) == c(nt, nt))) stop("B must be nt x nt.")
#  if (!all(dim(E) == c(nt, nt))) stop("E must be nt x nt.")
#  if (!all(dim(ssb_prior) == c(nt, nt))) stop("ssb_prior must be nt x nt.")
#  if (!all(dim(sse_prior) == c(nt, nt))) stop("sse_prior must be nt x nt.")
#
#  rownames(B) <- colnames(B) <- trait_names
#  rownames(E) <- colnames(E) <- trait_names
#  rownames(ssb_prior) <- colnames(ssb_prior) <- trait_names
#  rownames(sse_prior) <- colnames(sse_prior) <- trait_names
#
#  list(
#   vy = vy,
#   B = B,
#   E = E,
#   ssb_prior = ssb_prior,
#   sse_prior = sse_prior,
#   ssb_prior_list = split(ssb_prior, rep(seq_len(ncol(ssb_prior)), each = nrow(ssb_prior))),
#   sse_prior_list = split(sse_prior, rep(seq_len(ncol(sse_prior)), each = nrow(sse_prior)))
#  )
# }
#
# .stblr_init_marker_state <- function(nt, m, b_init = NULL, d_init = NULL) {
#  if (is.null(b_init)) b_init <- lapply(seq_len(nt), function(i) rep(0, m))
#  if (is.null(d_init)) d_init <- lapply(seq_len(nt), function(i) rep(0, m))
#
#  if (length(b_init) != nt || any(lengths(b_init) != m)) {
#   stop("b_init must be a list of length nt, each element length m.")
#  }
#
#  if (length(d_init) != nt || any(lengths(d_init) != m)) {
#   stop("d_init must be a list of length nt, each element length m.")
#  }
#
#  list(b_init = b_init, d_init = d_init)
# }
#
# .stblr_init_r_state <- function(stats, nt, m, use_r_init = FALSE, r_init = NULL) {
#  if (is.null(r_init)) r_init <- stats$wy
#
#  if (use_r_init) {
#   if (length(r_init) != nt || any(lengths(r_init) != m)) {
#    stop("r_init must be a list of length nt, each element length m, when use_r_init = TRUE.")
#   }
#  }
#
#  r_init
# }
#
# .stblr_prepare_annotation_matrix <- function(
#   A = NULL,
#   m,
#   variable_names = NULL,
#   add_intercept = FALSE,
#   standardize = TRUE,
#   center_binary = FALSE,
#   intercept_name = "Intercept"
# ) {
#  if (is.null(A)) {
#   if (!add_intercept) {
#    A <- matrix(numeric(0), nrow = m, ncol = 0)
#    rownames(A) <- variable_names
#    return(A)
#   }
#
#   A <- matrix(1, nrow = m, ncol = 1)
#   colnames(A) <- intercept_name
#   rownames(A) <- variable_names
#   return(A)
#  }
#
#  A <- as.matrix(A)
#  storage.mode(A) <- "double"
#
#  if (!is.null(rownames(A)) && !is.null(variable_names) &&
#      all(variable_names %in% rownames(A))) {
#   A <- A[variable_names, , drop = FALSE]
#  }
#
#  if (nrow(A) != m) {
#   stop("A must have one row per marker, or rownames(A) must contain all variable_names.")
#  }
#
#  if (any(!is.finite(A))) stop("A contains non-finite values.")
#
#  if (is.null(colnames(A))) colnames(A) <- paste0("Anno", seq_len(ncol(A)))
#
#  has_intercept <- ncol(A) >= 1 && all(abs(A[, 1] - 1) < 1e-12)
#
#  if (standardize && ncol(A) > 0) {
#   for (j in seq_len(ncol(A))) {
#    is_intercept <- all(abs(A[, j] - 1) < 1e-12)
#    is_binary <- all(A[, j] %in% c(0, 1))
#
#    if (is_intercept) next
#    if (is_binary && !center_binary) next
#
#    s <- stats::sd(A[, j])
#    if (is.finite(s) && s > 0) {
#     A[, j] <- (A[, j] - mean(A[, j])) / s
#    }
#   }
#  }
#
#  if (add_intercept && !has_intercept) {
#   A <- cbind(Intercept = 1, A)
#  }
#
#  if (!is.null(variable_names)) rownames(A) <- variable_names
#
#  A
# }
#
# .stblr_make_prior_from_annotations <- function(
#   A,
#   nt,
#   pi_base,
#   beta_pi = NULL,
#   beta_vb = NULL,
#   pi_min = 1e-8,
#   pi_max = 0.5,
#   vb_multiplier_min = 1e-3,
#   vb_multiplier_max = 1e3
# ) {
#  m <- nrow(A)
#  K <- ncol(A)
#
#  if (is.null(beta_pi)) beta_pi <- rep(0, K)
#  if (is.null(beta_vb)) beta_vb <- rep(0, K)
#
#  beta_pi <- as.matrix(beta_pi)
#  beta_vb <- as.matrix(beta_vb)
#
#  if (nrow(beta_pi) == 1 && K > 1) beta_pi <- t(beta_pi)
#  if (nrow(beta_vb) == 1 && K > 1) beta_vb <- t(beta_vb)
#
#  if (ncol(beta_pi) == 1 && nt > 1) beta_pi <- beta_pi[, rep(1, nt), drop = FALSE]
#  if (ncol(beta_vb) == 1 && nt > 1) beta_vb <- beta_vb[, rep(1, nt), drop = FALSE]
#
#  if (!all(dim(beta_pi) == c(K, nt))) stop("beta_pi must be K x nt or length K.")
#  if (!all(dim(beta_vb) == c(K, nt))) stop("beta_vb must be K x nt or length K.")
#
#  pi_marker <- vector("list", nt)
#  vb_multiplier <- vector("list", nt)
#
#  for (t in seq_len(nt)) {
#   lp_pi <- as.numeric(A %*% beta_pi[, t])
#   lp_pi <- lp_pi - mean(lp_pi)
#
#   p <- stats::plogis(stats::qlogis(pi_base) + lp_pi)
#   p <- pmin(pmax(p, pi_min), pi_max)
#
#   lp_vb <- as.numeric(A %*% beta_vb[, t])
#   lp_vb <- lp_vb - mean(lp_vb)
#
#   mult <- exp(lp_vb)
#   mult <- pmin(pmax(mult, vb_multiplier_min), vb_multiplier_max)
#
#   pi_marker[[t]] <- p
#   vb_multiplier[[t]] <- mult
#  }
#
#  list(
#   pi_marker = pi_marker,
#   vb_multiplier = vb_multiplier
#  )
# }
#
# format_csr_annot_fit <- function(
#   fit,
#   nt,
#   m,
#   annotation_names,
#   trait_names = NULL,
#   variable_names = NULL
# ) {
#  if (length(fit) < 22) {
#   stop("format_csr_annot_fit() expects the updated 22-slot annotation CSR return object.")
#  }
#
#  names(fit)[1:22] <- c(
#   "bm", "dm", "wy", "r", "b", "d", "o",
#   "vbs", "vgs", "ves",
#   "covb", "covg", "cove",
#   "vb", "vg", "ve",
#   "pi", "pim", "eta_pi", "eta_vb",
#   "vle", "vld"
#  )
#
#  if (is.null(trait_names)) trait_names <- paste0("T", seq_len(nt))
#  if (is.null(variable_names)) variable_names <- paste0("V", seq_len(m))
#
#  K <- length(annotation_names)
#  if (K == 0) annotation_names <- "none"
#
#  for (i in 1:7) {
#   fit[[i]] <- as.matrix(as.data.frame(fit[[i]]))
#   rownames(fit[[i]]) <- variable_names
#   colnames(fit[[i]]) <- trait_names
#  }
#
#  for (i in 8:10) {
#   fit[[i]] <- as.matrix(as.data.frame(fit[[i]]))
#   rownames(fit[[i]]) <- paste0("Iter", seq_len(nrow(fit[[i]])))
#   colnames(fit[[i]]) <- trait_names
#  }
#
#  for (i in 11:16) {
#   fit[[i]] <- matrix(unlist(fit[[i]]), ncol = nt, byrow = TRUE)
#   rownames(fit[[i]]) <- trait_names
#   colnames(fit[[i]]) <- trait_names
#  }
#
#  fit$pi <- matrix(unlist(fit$pi), ncol = 2, byrow = TRUE)
#  rownames(fit$pi) <- trait_names
#  colnames(fit$pi) <- c("pi0", "pi1")
#
#  fit$pim <- matrix(unlist(fit$pim), ncol = 2, byrow = TRUE)
#  rownames(fit$pim) <- trait_names
#  colnames(fit$pim) <- c("pi0", "pi1")
#
#  eta_pi <- matrix(unlist(fit$eta_pi), nrow = nt, byrow = TRUE)
#  eta_vb <- matrix(unlist(fit$eta_vb), nrow = nt, byrow = TRUE)
#  rownames(eta_pi) <- rownames(eta_vb) <- trait_names
#  colnames(eta_pi) <- colnames(eta_vb) <- annotation_names
#
#  vle <- as.matrix(as.data.frame(fit$vle))
#  vld <- as.matrix(as.data.frame(fit$vld))
#  rownames(vle) <- paste0("Iter", seq_len(nrow(vle)))
#  rownames(vld) <- paste0("Iter", seq_len(nrow(vld)))
#  colnames(vle) <- trait_names
#  colnames(vld) <- trait_names
#
#  out <- fit[1:18]
#  out$eta_pi <- eta_pi
#  out$eta_vb <- eta_vb
#  out$vle <- vle
#  out$vld <- vld
#
#  if (sum(diag(out$covb)) > 0) out$rb <- cov2cor(out$covb)
#  if (sum(diag(out$covg)) > 0) out$rg <- cov2cor(out$covg)
#  if (sum(diag(out$cove)) > 0) out$re <- cov2cor(out$cove)
#
#  out
# }
#
# format_csr_prior_fit <- function(
#   fit,
#   nt,
#   m,
#   trait_names = NULL,
#   variable_names = NULL
# ) {
#  format_stblr_fit(
#   fit = fit,
#   nt = nt,
#   m = m,
#   trait_names = trait_names,
#   variable_names = variable_names
#  )
# }
#
# # =============================================================================
# # Main generic interface
# # =============================================================================
#
# stblr_csr_annotation <- function(
#   stats,
#   ld_prefix,
#   model = c("prior", "annot", "sbayesrc"),
#   A = NULL,
#   n = NULL,
#   m = NULL,
#
#   # Shared architecture controls
#   pi_marker = 0.001,
#   pi_init = NULL,
#   pi_vb_init = NULL,
#   pi_prior_mean = NULL,
#   pi_prior_strength = NULL,
#   pi_prior_a = NULL,
#   pi_prior_b = NULL,
#
#   h2 = 0.5,
#   nub = 4,
#   nue = 4,
#   B = NULL,
#   E = NULL,
#   ssb_prior = NULL,
#   sse_prior = NULL,
#
#   updateB = TRUE,
#   updateE = TRUE,
#   updatePi = TRUE,
#   adjE = 0.9,
#
#   nit = 1000,
#   nburn = 100,
#   nthin = 1,
#   ncores = 3,
#   seed = 10,
#
#   b_init = NULL,
#   d_init = NULL,
#   use_d_init = FALSE,
#   r_init = NULL,
#   use_r_init = FALSE,
#   rebuild_r_before_updateE = FALSE,
#
#   # Annotation preprocessing
#   add_intercept = FALSE,
#   standardize_annotations = TRUE,
#   center_binary_annotations = FALSE,
#
#   # model = "prior" options
#   use_pi_marker = FALSE,
#   use_vb_multiplier = FALSE,
#   fixed_pi_marker = NULL,
#   fixed_vb_multiplier = NULL,
#   beta_pi = NULL,
#   beta_vb = NULL,
#
#   # model = "annot" options
#   learn_pi_annot = TRUE,
#   learn_vb_annot = FALSE,
#   eta_pi_init = NULL,
#   eta_vb_init = NULL,
#   sigma_eta_pi = 1,
#   sigma_eta_vb = 1,
#   rw_sd_eta_pi = 0.05,
#   rw_sd_eta_vb = 0.05,
#   annot_update_every = 10,
#   pi_min = 1e-8,
#   pi_max = 0.5,
#   vb_multiplier_min = 1e-3,
#   vb_multiplier_max = 1e3,
#
#   # model = "sbayesrc" options
#   gamma = c(0, 0.01, 0.1, 1),
#   active_comp_weights = NULL,
#   alpha_init = NULL,
#   sigmaSqAlpha_init = NULL,
#   intercept_flat = TRUE,
#   sigmaSqAlpha_a = 2,
#   sigmaSqAlpha_b = 2,
#   pi_floor = 1e-12,
#   updateAlpha = TRUE,
#   alpha_update_every = 10,
#   comp_init = NULL,
#   use_comp_init = FALSE
# ) {
#  model <- .stblr_match_annotation_model(model)
#
#  dims <- .stblr_get_nt_m_names(stats, n = n, m = m)
#  nt <- dims$nt
#  n <- dims$n
#  m <- dims$m
#  trait_names <- dims$trait_names
#  variable_names <- dims$variable_names
#
#  .stblr_validate_stats(stats, nt = nt, m = m)
#
#  arch <- .stblr_resolve_architecture(
#   pi_marker = pi_marker,
#   pi_init = pi_init,
#   pi_vb_init = pi_vb_init,
#   pi_prior_mean = pi_prior_mean,
#   pi_prior_strength = pi_prior_strength,
#   pi_prior_a = pi_prior_a,
#   pi_prior_b = pi_prior_b
#  )
#
#  pri <- .stblr_make_csr_variance_priors(
#   stats = stats,
#   n = n,
#   m = m,
#   nt = nt,
#   h2 = h2,
#   nub = nub,
#   nue = nue,
#   pi_vb_init = arch$pi_vb_init,
#   pi_prior_mean = arch$pi_prior_mean,
#   trait_names = trait_names,
#   B = B,
#   E = E,
#   ssb_prior = ssb_prior,
#   sse_prior = sse_prior
#  )
#
#  state <- .stblr_init_marker_state(nt = nt, m = m, b_init = b_init, d_init = d_init)
#  r_init <- .stblr_init_r_state(stats, nt = nt, m = m, use_r_init = use_r_init, r_init = r_init)
#
#  if (model %in% c("prior", "annot", "sbayesrc")) {
#   A <- .stblr_prepare_annotation_matrix(
#    A = A,
#    m = m,
#    variable_names = variable_names,
#    add_intercept = add_intercept || model == "sbayesrc",
#    standardize = standardize_annotations,
#    center_binary = center_binary_annotations
#   )
#  }
#
#  annotation_names <- colnames(A)
#
#  if (model == "prior") {
#   if (!is.null(fixed_pi_marker)) {
#    use_pi_marker <- TRUE
#    pi_marker_list <- fixed_pi_marker
#   } else if (use_pi_marker || !is.null(beta_pi)) {
#    ann_prior <- .stblr_make_prior_from_annotations(
#     A = A,
#     nt = nt,
#     pi_base = arch$pi_init,
#     beta_pi = beta_pi,
#     beta_vb = beta_vb,
#     pi_min = pi_min,
#     pi_max = pi_max,
#     vb_multiplier_min = vb_multiplier_min,
#     vb_multiplier_max = vb_multiplier_max
#    )
#    pi_marker_list <- ann_prior$pi_marker
#   } else {
#    pi_marker_list <- lapply(seq_len(nt), function(i) rep(arch$pi_init, m))
#   }
#
#   if (!is.null(fixed_vb_multiplier)) {
#    use_vb_multiplier <- TRUE
#    vb_multiplier_list <- fixed_vb_multiplier
#   } else if (use_vb_multiplier || !is.null(beta_vb)) {
#    if (!exists("ann_prior", inherits = FALSE)) {
#     ann_prior <- .stblr_make_prior_from_annotations(
#      A = A,
#      nt = nt,
#      pi_base = arch$pi_init,
#      beta_pi = beta_pi,
#      beta_vb = beta_vb,
#      pi_min = pi_min,
#      pi_max = pi_max,
#      vb_multiplier_min = vb_multiplier_min,
#      vb_multiplier_max = vb_multiplier_max
#     )
#    }
#    vb_multiplier_list <- ann_prior$vb_multiplier
#   } else {
#    vb_multiplier_list <- lapply(seq_len(nt), function(i) rep(1, m))
#   }
#
#   raw_fit <- stblr_cpg_omp_csr_prior(
#    wy = stats$wy,
#    ww = stats$ww,
#    yy = stats$yy,
#    b_init = state$b_init,
#    d_init = state$d_init,
#    use_d_init = use_d_init,
#    r_init = r_init,
#    use_r_init = use_r_init,
#    rebuild_r_before_updateE = rebuild_r_before_updateE,
#    ld_prefix = ld_prefix,
#    B = pri$B,
#    E = pri$E,
#    ssb_prior = pri$ssb_prior_list,
#    sse_prior = pri$sse_prior_list,
#    pi = arch$pi,
#    use_pi_marker = use_pi_marker,
#    pi_marker = pi_marker_list,
#    use_vb_multiplier = use_vb_multiplier,
#    vb_multiplier = vb_multiplier_list,
#    nub = nub,
#    nue = nue,
#    updateB = updateB,
#    updateE = updateE,
#    updatePi = updatePi,
#    adjE = adjE,
#    n = rep(as.integer(n), nt),
#    nit = as.integer(nit),
#    nburn = as.integer(nburn),
#    nthin = as.integer(nthin),
#    pi_prior_a = arch$pi_prior_a,
#    pi_prior_b = arch$pi_prior_b,
#    ncores = as.integer(ncores),
#    seed = as.integer(seed)
#   )
#
#   fit <- format_csr_prior_fit(
#    fit = raw_fit,
#    nt = nt,
#    m = m,
#    trait_names = trait_names,
#    variable_names = variable_names
#   )
#
#   fit$input <- list(
#    model = model,
#    n = n,
#    m = m,
#    nt = nt,
#    A = A,
#    annotation_names = annotation_names,
#    use_pi_marker = use_pi_marker,
#    use_vb_multiplier = use_vb_multiplier,
#    pi_marker = if (use_pi_marker) pi_marker_list else NULL,
#    vb_multiplier = if (use_vb_multiplier) vb_multiplier_list else NULL
#   )
#  }
#
#  if (model == "annot") {
#   K <- ncol(A)
#
#   if (K == 0) {
#    stop("model = 'annot' requires at least one annotation column in A.")
#   }
#
#   if (is.null(eta_pi_init)) eta_pi_init <- matrix(0, nrow = K, ncol = nt)
#   if (is.null(eta_vb_init)) eta_vb_init <- matrix(0, nrow = K, ncol = nt)
#
#   eta_pi_init <- as.matrix(eta_pi_init)
#   eta_vb_init <- as.matrix(eta_vb_init)
#   storage.mode(eta_pi_init) <- "double"
#   storage.mode(eta_vb_init) <- "double"
#
#   if (!all(dim(eta_pi_init) == c(K, nt))) stop("eta_pi_init must be ncol(A) x nt.")
#   if (!all(dim(eta_vb_init) == c(K, nt))) stop("eta_vb_init must be ncol(A) x nt.")
#
#   raw_fit <- stblr_cpg_omp_csr_annot(
#    wy = stats$wy,
#    ww = stats$ww,
#    yy = stats$yy,
#    b_init = state$b_init,
#    d_init = state$d_init,
#    use_d_init = use_d_init,
#    r_init = r_init,
#    use_r_init = use_r_init,
#    rebuild_r_before_updateE = rebuild_r_before_updateE,
#    ld_prefix = ld_prefix,
#    B = pri$B,
#    E = pri$E,
#    ssb_prior = pri$ssb_prior_list,
#    sse_prior = pri$sse_prior_list,
#    pi = arch$pi,
#    A = A,
#    learn_pi_annot = learn_pi_annot,
#    learn_vb_annot = learn_vb_annot,
#    eta_pi_init = eta_pi_init,
#    eta_vb_init = eta_vb_init,
#    sigma_eta_pi = sigma_eta_pi,
#    sigma_eta_vb = sigma_eta_vb,
#    rw_sd_eta_pi = rw_sd_eta_pi,
#    rw_sd_eta_vb = rw_sd_eta_vb,
#    annot_update_every = as.integer(annot_update_every),
#    pi_min = pi_min,
#    pi_max = pi_max,
#    vb_multiplier_min = vb_multiplier_min,
#    vb_multiplier_max = vb_multiplier_max,
#    nub = nub,
#    nue = nue,
#    updateB = updateB,
#    updateE = updateE,
#    updatePi = updatePi,
#    adjE = adjE,
#    n = rep(as.integer(n), nt),
#    nit = as.integer(nit),
#    nburn = as.integer(nburn),
#    nthin = as.integer(nthin),
#    pi_prior_a = arch$pi_prior_a,
#    pi_prior_b = arch$pi_prior_b,
#    ncores = as.integer(ncores),
#    seed = as.integer(seed)
#   )
#
#   fit <- format_csr_annot_fit(
#    fit = raw_fit,
#    nt = nt,
#    m = m,
#    annotation_names = annotation_names,
#    trait_names = trait_names,
#    variable_names = variable_names
#   )
#
#   fit$input <- list(
#    model = model,
#    n = n,
#    m = m,
#    nt = nt,
#    A = A,
#    annotation_names = annotation_names,
#    learn_pi_annot = learn_pi_annot,
#    learn_vb_annot = learn_vb_annot,
#    eta_pi_init = eta_pi_init,
#    eta_vb_init = eta_vb_init,
#    sigma_eta_pi = sigma_eta_pi,
#    sigma_eta_vb = sigma_eta_vb,
#    rw_sd_eta_pi = rw_sd_eta_pi,
#    rw_sd_eta_vb = rw_sd_eta_vb,
#    annot_update_every = annot_update_every
#   )
#  }
#
#  if (model == "sbayesrc") {
#   gamma <- as.numeric(gamma)
#   Kgamma <- length(gamma)
#   if (Kgamma < 2) stop("gamma must have at least two elements.")
#   if (!isTRUE(all.equal(gamma[1], 0))) stop("gamma[1] must be 0.")
#
#   if (is.null(comp_init)) {
#    comp_init <- lapply(seq_len(nt), function(i) rep(0, m))
#   }
#   if (length(comp_init) != nt || any(lengths(comp_init) != m)) {
#    stop("comp_init must be a list of length nt, each element length m.")
#   }
#
#   alpha <- make_sbayesrc_alpha_init(
#    A = A,
#    gamma = gamma,
#    pi_init = arch$pi_init,
#    active_comp_weights = active_comp_weights,
#    alpha_init = alpha_init,
#    sigmaSqAlpha_init = sigmaSqAlpha_init
#   )
#
#   raw_fit <- stblr_cpg_omp_csr_sbayesrc(
#    wy = stats$wy,
#    ww = stats$ww,
#    yy = stats$yy,
#    b_init = state$b_init,
#    comp_init = comp_init,
#    use_comp_init = use_comp_init,
#    r_init = r_init,
#    use_r_init = use_r_init,
#    rebuild_r_before_updateE = rebuild_r_before_updateE,
#    ld_prefix = ld_prefix,
#    B = pri$B,
#    E = pri$E,
#    ssb_prior = pri$ssb_prior_list,
#    sse_prior = pri$sse_prior_list,
#    A = A,
#    gamma = gamma,
#    alpha_init = alpha$alpha_init,
#    sigmaSqAlpha_init = alpha$sigmaSqAlpha_init,
#    intercept_flat = intercept_flat,
#    sigmaSqAlpha_a = sigmaSqAlpha_a,
#    sigmaSqAlpha_b = sigmaSqAlpha_b,
#    pi_floor = pi_floor,
#    nub = nub,
#    nue = nue,
#    updateAlpha = updateAlpha,
#    updateB = updateB,
#    updateE = updateE,
#    alpha_update_every = as.integer(alpha_update_every),
#    adjE = adjE,
#    n = rep(as.integer(n), nt),
#    nit = as.integer(nit),
#    nburn = as.integer(nburn),
#    nthin = as.integer(nthin),
#    ncores = as.integer(ncores),
#    seed = as.integer(seed)
#   )
#
#   fit <- format_sbayesrc_csr_fit(
#    fit = raw_fit,
#    nt = nt,
#    m = m,
#    gamma = gamma,
#    n_anno = ncol(A),
#    trait_names = trait_names,
#    variable_names = variable_names,
#    annotation_names = annotation_names
#   )
#
#   fit$input <- list(
#    model = model,
#    n = n,
#    m = m,
#    nt = nt,
#    A = A,
#    annotation_names = annotation_names,
#    gamma = gamma,
#    component_names = colnames(fit$ncomp),
#    active_comp_weights = alpha$active_comp_weights,
#    component_prob_init = alpha$component_prob_init,
#    step_prob_init = alpha$step_prob_init,
#    alpha_init = alpha$alpha_init,
#    sigmaSqAlpha_init = alpha$sigmaSqAlpha_init,
#    intercept_flat = intercept_flat,
#    sigmaSqAlpha_a = sigmaSqAlpha_a,
#    sigmaSqAlpha_b = sigmaSqAlpha_b,
#    pi_floor = pi_floor,
#    updateAlpha = updateAlpha,
#    alpha_update_every = alpha_update_every
#   )
#  }
#
#  shared_input <- list(
#   pi_marker = arch$pi_marker,
#   pi_init = arch$pi_init,
#   pi_vb_init = arch$pi_vb_init,
#   pi_prior_mean = arch$pi_prior_mean,
#   pi_prior_strength = arch$pi_prior_strength,
#   pi_prior_a = arch$pi_prior_a,
#   pi_prior_b = arch$pi_prior_b,
#   h2 = h2,
#   nub = nub,
#   nue = nue,
#   vy = pri$vy,
#   B = pri$B,
#   E = pri$E,
#   ssb_prior = pri$ssb_prior,
#   sse_prior = pri$sse_prior,
#   updateB = updateB,
#   updateE = updateE,
#   updatePi = updatePi,
#   adjE = adjE,
#   nit = nit,
#   nburn = nburn,
#   nthin = nthin,
#   ncores = ncores,
#   seed = seed,
#   use_d_init = use_d_init,
#   use_r_init = use_r_init,
#   rebuild_r_before_updateE = rebuild_r_before_updateE,
#   ld_prefix = ld_prefix
#  )
#
#  fit$input <- c(fit$input, shared_input)
#  fit
# }
#
# # Convenience wrappers ---------------------------------------------------------
#
# stblr_csr_prior_annot <- function(...) {
#  stblr_csr_annotation(..., model = "prior")
# }
#
# stblr_csr_learn_annot <- function(...) {
#  stblr_csr_annotation(..., model = "annot")
# }
#
# stblr_csr_sbayesrc_generic <- function(...) {
#  stblr_csr_annotation(..., model = "sbayesrc")
# }

fit_prior <- stblr_csr_annotation(
 stats = stats,
 ld_prefix = ld_prefix,
 model = "prior",
 A = A,
 n = Glist$n,
 beta_pi = c(0.5, -0.2),
 beta_vb = c(0.0, 0.3),
 use_pi_marker = TRUE,
 use_vb_multiplier = TRUE
)

fit_annot <- stblr_csr_annotation(
 stats = stats,
 ld_prefix = ld_prefix,
 model = "annot",
 A = A,
 n = Glist$n,
 learn_pi_annot = TRUE,
 learn_vb_annot = TRUE,
 rw_sd_eta_pi = 0.05,
 rw_sd_eta_vb = 0.05
)

fit_rc <- stblr_csr_annotation(
 stats = stats,
 ld_prefix = ld_prefix,
 model = "sbayesrc",
 A = A,
 n = Glist$n,
 gamma = c(0, 0.01, 0.1, 1)
)


fit_group <- stblr_csr_annotation(
 stats = stats,
 ld_prefix = ld_prefix,
 model = "group",
 group = marker_groups,
 n = Glist$n,

 pi_init = 0.001,
 pi_vb_init = 0.001,
 pi_prior_mean = 0.001,

 updatePi = TRUE,
 updateGroupVb = TRUE,

 nit = 1000,
 nburn = 100,
 nthin = 1,
 ncores = 3,
 seed = 10
)

fit_group <- stblr_csr_group_annot(
 stats = stats,
 ld_prefix = ld_prefix,
 group = marker_groups,
 n = Glist$n,
 updatePi = TRUE,
 updateGroupVb = TRUE
)










format_stblr_fit <- function(
  fit,
  nt,
  m,
  trait_names = NULL,
  variable_names = NULL,
  annotation_names = NULL
) {
 names(fit) <- c(
  "bm", "dm", "wy", "r", "b", "d", "o",
  "vbs", "vgs", "ves",
  "covb", "covg", "cove",
  "vb", "vg", "ve",
  "pi", "pim",
  "eta_pi", "eta_vb"
 )

 if (is.null(trait_names)) {
  trait_names <- paste0("T", seq_len(nt))
 }

 if (is.null(variable_names)) {
  variable_names <- paste0("V", seq_len(m))
 }

 for (i in 1:7) {
  fit[[i]] <- as.matrix(as.data.frame(fit[[i]]))
  rownames(fit[[i]]) <- variable_names
  colnames(fit[[i]]) <- trait_names
 }

 for (i in 8:10) {
  fit[[i]] <- as.matrix(as.data.frame(fit[[i]]))
  rownames(fit[[i]]) <- paste0("Iter", seq_len(nrow(fit[[i]])))
  colnames(fit[[i]]) <- trait_names
 }

 for (i in 11:16) {
  fit[[i]] <- matrix(unlist(fit[[i]]), ncol = nt, byrow = TRUE)
  colnames(fit[[i]]) <- rownames(fit[[i]]) <- trait_names
 }

 fit[[17]] <- do.call(rbind, fit[[17]])
 fit[[18]] <- do.call(rbind, fit[[18]])

 rownames(fit[[17]]) <- trait_names
 rownames(fit[[18]]) <- trait_names
 colnames(fit[[17]]) <- c("pi0", "pi1")
 colnames(fit[[18]]) <- c("pi0", "pi1")

 # Annotation posterior means: currently returned as nt x K lists.
 fit[[19]] <- do.call(rbind, fit[[19]])
 fit[[20]] <- do.call(rbind, fit[[20]])

 rownames(fit[[19]]) <- trait_names
 rownames(fit[[20]]) <- trait_names

 if (is.null(annotation_names)) {
  annotation_names <- paste0("A", seq_len(ncol(fit[[19]])))
 }

 colnames(fit[[19]]) <- annotation_names
 colnames(fit[[20]]) <- annotation_names

 if (sum(diag(fit$covb)) > 0) fit$rb <- cov2cor(fit$covb)
 if (sum(diag(fit$covg)) > 0) fit$rg <- cov2cor(fit$covg)
 if (sum(diag(fit$cove)) > 0) fit$re <- cov2cor(fit$cove)

 fit
}

format_stblr_fit <- function(
  fit,
  nt,
  m,
  trait_names = NULL,
  variable_names = NULL,
  annotation_names = NULL,
  gamma = NULL
) {
 names(fit) <- c(
  "bm", "dm", "wy", "r", "b", "d", "o",
  "vbs", "vgs", "ves",
  "covb", "covg", "cove",
  "vb", "vg", "ve",
  "pi", "pim", "pitrait", "pimarker"
 )

 if (is.null(trait_names)) {
  trait_names <- paste0("T", seq_len(nt))
 }

 if (is.null(variable_names)) {
  variable_names <- paste0("V", seq_len(m))
 }

 for (i in 1:7) {
  fit[[i]] <- as.matrix(as.data.frame(fit[[i]]))
  rownames(fit[[i]]) <- variable_names
  colnames(fit[[i]]) <- trait_names
 }

 for (i in 8:10) {
  fit[[i]] <- as.matrix(as.data.frame(fit[[i]]))
  rownames(fit[[i]]) <- paste0("Iter", seq_len(nrow(fit[[i]])))
  colnames(fit[[i]]) <- trait_names
 }

 for (i in 11:16) {
  fit[[i]] <- matrix(unlist(fit[[i]]), ncol = nt, byrow = TRUE)
  colnames(fit[[i]]) <- rownames(fit[[i]]) <- trait_names
 }

 fit[[17]] <- do.call(rbind, fit[[17]])
 fit[[18]] <- do.call(rbind, fit[[18]])

 rownames(fit[[17]]) <- trait_names
 rownames(fit[[18]]) <- trait_names
 colnames(fit[[17]]) <- c("pi0", "pi1")
 colnames(fit[[18]]) <- c("pi0", "pi1")

 # Keep base 18-slot output first
 out <- fit[1:18]

 if (sum(diag(out$covb)) > 0) out$rb <- cov2cor(out$covb)
 if (sum(diag(out$covg)) > 0) out$rg <- cov2cor(out$covg)
 if (sum(diag(out$cove)) > 0) out$re <- cov2cor(out$cove)

 # Optional SBayesRC annotation output
 if (length(fit) >= 20 && !is.null(annotation_names) && !is.null(gamma)) {
  ndist <- length(gamma)
  n_alpha <- ndist - 1
  K <- length(annotation_names)

  alpha_raw <- do.call(rbind, fit[[19]])
  sigma_raw <- do.call(rbind, fit[[20]])

  if (ncol(alpha_raw) == K * n_alpha) {
   alpha_names <- as.vector(outer(
    annotation_names,
    paste0("p", 2:ndist),
    paste,
    sep = ":"
   ))

   colnames(alpha_raw) <- alpha_names
   rownames(alpha_raw) <- trait_names

   out$alpha_flat <- alpha_raw

   out$alpha <- vector("list", nt)
   names(out$alpha) <- trait_names

   for (t in seq_len(nt)) {
    out$alpha[[t]] <- matrix(
     alpha_raw[t, ],
     nrow = K,
     ncol = n_alpha,
     byrow = FALSE
    )
    rownames(out$alpha[[t]]) <- annotation_names
    colnames(out$alpha[[t]]) <- paste0("p", 2:ndist)
   }
  }

  if (ncol(sigma_raw) == n_alpha) {
   colnames(sigma_raw) <- paste0("p", 2:ndist)
   rownames(sigma_raw) <- trait_names
   out$sigmaSqAlpha <- sigma_raw
  }
 }

 out
}

recovery_summary <- function(
  fit,
  Btrue,
  Ks = c(40, 60, 100, 150, 250),
  trait_names = NULL
) {
 if (is.null(fit$dm)) stop("fit$dm is missing.")
 if (is.null(fit$bm)) stop("fit$bm is missing.")

 if (!is.matrix(Btrue)) Btrue <- as.matrix(Btrue)
 if (!is.matrix(fit$dm)) fit$dm <- as.matrix(fit$dm)
 if (!is.matrix(fit$bm)) fit$bm <- as.matrix(fit$bm)

 if (!all(dim(fit$dm) == dim(Btrue))) {
  stop("fit$dm and Btrue must have the same dimensions.")
 }

 if (!all(dim(fit$bm) == dim(Btrue))) {
  stop("fit$bm and Btrue must have the same dimensions.")
 }

 m <- nrow(Btrue)
 nt <- ncol(Btrue)

 if (is.null(trait_names)) {
  trait_names <- paste0("D", seq_len(nt))
 }

 if (length(trait_names) != nt) {
  stop("trait_names must have length equal to ncol(Btrue).")
 }

 Ks <- sort(unique(as.integer(Ks)))
 Ks <- Ks[Ks > 0]

 if (length(Ks) == 0) {
  stop("Ks must contain at least one positive integer.")
 }

 out <- do.call(rbind, lapply(seq_len(nt), function(t) {
  causal <- Btrue[, t] != 0
  causal[is.na(causal)] <- FALSE

  n_causal <- sum(causal)

  dm_score <- fit$dm[, t]
  bm_score <- abs(fit$bm[, t])

  dm_score[is.na(dm_score)] <- -Inf
  bm_score[is.na(bm_score)] <- -Inf

  rank_dm <- order(dm_score, decreasing = TRUE)
  rank_bm <- order(bm_score, decreasing = TRUE)

  do.call(rbind, lapply(Ks, function(K) {
   K_eff <- min(K, m)

   top_dm <- rank_dm[seq_len(K_eff)]
   top_bm <- rank_bm[seq_len(K_eff)]

   recovered_dm <- sum(causal[top_dm])
   recovered_bm <- sum(causal[top_bm])

   data.frame(
    trait = t,
    trait_name = trait_names[t],
    K = K,
    K_eff = K_eff,
    n_causal = n_causal,

    recovered_dm = recovered_dm,
    precision_dm = recovered_dm / K_eff,
    recall_dm = ifelse(n_causal > 0, recovered_dm / n_causal, NA_real_),

    recovered_bm = recovered_bm,
    precision_bm = recovered_bm / K_eff,
    recall_bm = ifelse(n_causal > 0, recovered_bm / n_causal, NA_real_)
   )
  }))
 }))

 rownames(out) <- NULL
 out
}

derive_marker_priors_from_fit <- function(
  fit,
  anno,
  base_pi = NULL,
  use_traits = NULL,
  pi_floor = 1e-8,
  pi_cap = 0.20,
  vb_floor = 0.25,
  vb_cap = 4,
  shrink = 0.5,
  normalize = TRUE,
  return_list = TRUE
) {
 if (is.null(fit$dm)) stop("fit$dm is missing.")
 if (is.null(fit$bm)) stop("fit$bm is missing.")

 dm <- as.matrix(fit$dm)
 bm <- as.matrix(fit$bm)
 anno <- as.matrix(anno)

 m <- nrow(dm)
 nt <- ncol(dm)

 if (nrow(anno) != m) {
  stop("anno must have the same number of rows as fit$dm.")
 }

 if (is.null(use_traits)) {
  use_traits <- seq_len(nt)
 }

 if (is.null(base_pi)) {
  base_pi <- pmax(colMeans(dm, na.rm = TRUE), pi_floor)
 }

 if (length(base_pi) == 1) {
  base_pi <- rep(base_pi, nt)
 }

 if (length(base_pi) != nt) {
  stop("base_pi must be scalar or length ncol(fit$dm).")
 }

 A <- ncol(anno)

 if (is.null(colnames(anno))) {
  annotation_names <- paste0("A", seq_len(A))
 } else {
  annotation_names <- colnames(anno)
 }

 anno_size <- colSums(anno != 0)

 pi_marker_mat <- matrix(NA_real_, m, nt)
 vb_multiplier_mat <- matrix(NA_real_, m, nt)

 colnames(pi_marker_mat) <- colnames(dm)
 colnames(vb_multiplier_mat) <- colnames(dm)

 rownames(pi_marker_mat) <- rownames(dm)
 rownames(vb_multiplier_mat) <- rownames(dm)

 annotation_summary <- vector("list", nt)

 for (t in seq_len(nt)) {
  marker_signal <- dm[, t]
  marker_effect <- abs(bm[, t])

  marker_signal[is.na(marker_signal)] <- 0
  marker_effect[is.na(marker_effect)] <- 0

  anno_dm <- rep(NA_real_, A)
  anno_bm <- rep(NA_real_, A)

  for (a in seq_len(A)) {
   idx <- anno[, a] != 0

   if (sum(idx) > 0) {
    anno_dm[a] <- mean(marker_signal[idx], na.rm = TRUE)
    anno_bm[a] <- mean(marker_effect[idx], na.rm = TRUE)
   }
  }

  global_dm <- mean(marker_signal, na.rm = TRUE)
  global_bm <- mean(marker_effect, na.rm = TRUE)

  if (!is.finite(global_dm) || global_dm <= 0) {
   global_dm <- pi_floor
  }

  if (!is.finite(global_bm) || global_bm <= 0) {
   positive_bm <- marker_effect[marker_effect > 0]
   global_bm <- ifelse(length(positive_bm) > 0, mean(positive_bm), 1)
  }

  anno_pi_enrichment <- anno_dm / global_dm
  anno_vb_enrichment <- anno_bm / global_bm

  anno_pi_enrichment[!is.finite(anno_pi_enrichment)] <- 1
  anno_vb_enrichment[!is.finite(anno_vb_enrichment)] <- 1

  anno_pi_enrichment <- 1 + shrink * (anno_pi_enrichment - 1)
  anno_vb_enrichment <- 1 + shrink * (anno_vb_enrichment - 1)

  # Log-scale combination is safer for overlapping annotations.
  log_pi_weight <- as.numeric(anno %*% log(pmax(anno_pi_enrichment, 1e-8)))
  log_vb_weight <- as.numeric(anno %*% log(pmax(anno_vb_enrichment, 1e-8)))

  pi_weight <- exp(log_pi_weight)
  vb_weight <- exp(log_vb_weight)

  pi_weight[!is.finite(pi_weight) | pi_weight <= 0] <- 1
  vb_weight[!is.finite(vb_weight) | vb_weight <= 0] <- 1

  if (normalize) {
   pi_weight <- pi_weight / mean(pi_weight, na.rm = TRUE)
   vb_weight <- vb_weight / mean(vb_weight, na.rm = TRUE)
  }

  pi_marker_mat[, t] <- pmin(
   pmax(base_pi[t] * pi_weight, pi_floor),
   pi_cap
  )

  vb_multiplier_mat[, t] <- pmin(
   pmax(vb_weight, vb_floor),
   vb_cap
  )

  annotation_summary[[t]] <- data.frame(
   trait = t,
   trait_name = if (!is.null(colnames(dm))) colnames(dm)[t] else paste0("T", t),
   annotation = annotation_names,
   size = as.integer(anno_size),
   mean_dm = anno_dm,
   mean_abs_bm = anno_bm,
   pi_enrichment = anno_pi_enrichment,
   vb_enrichment = anno_vb_enrichment
  )
 }

 annotation_summary <- do.call(rbind, annotation_summary)
 rownames(annotation_summary) <- NULL

 if (return_list) {
  pi_marker <- lapply(seq_len(nt), function(t) pi_marker_mat[, t])
  vb_multiplier <- lapply(seq_len(nt), function(t) vb_multiplier_mat[, t])

  names(pi_marker) <- colnames(dm)
  names(vb_multiplier) <- colnames(dm)
 } else {
  pi_marker <- pi_marker_mat
  vb_multiplier <- vb_multiplier_mat
 }

 list(
  pi_marker = pi_marker,
  vb_multiplier = vb_multiplier,
  pi_marker_mat = pi_marker_mat,
  vb_multiplier_mat = vb_multiplier_mat,
  annotation_summary = annotation_summary
 )
}

library(qgg)
library(sblr)

Glist <- readRDS(file = file.path(data_dir, "Glist_sparseLD_1k.RDS"))

chr <- 1
rsids <- Glist$rsidsLD[[chr]]
h2 <- c(0.4, 0.5, 0.3)
rg <- matrix(
 c(
  1.0, 0.7, 0.3,
  0.7, 1.0, 0.5,
  0.3, 0.5, 1.0
 ),
 nrow = 3,
 byrow = TRUE
)

sim <- mtsim(
 Glist = Glist, chr=chr, rsids = rsids,
 nt = 3, n_shared = 30, n_specific = 10,
 h2 = h2, rg = rg, re = 0,
 seed = 1
)

sim <- mtsim_annotation(
 Glist = Glist,
 chr = chr,
 nt = 3,
 n_shared = 30,
 n_specific = 30,
 h2 = 0.5,
 rg = matrix(c(
  1, 0.2, 0.2,
  0.2, 1, 0.2,
  0.2, 0.2, 1
 ), 3, 3),

 n_annotations = 5,
 annotation_prob = c(0.03, 0.05, 0.10, 0.15, 0.20),

 enriched_annotations = c(1, 2),
 annotation_enrichment = 10,

 base_pi = 0.001,
 enriched_pi_multiplier = c(10, 5),
 enriched_vb_multiplier = c(4, 2),

 seed = 10
)

stat <- glma(y = scale(sim$y[,1]), rsids=Glist$rsidsLD[[1]], Glist = Glist)
system.time(fitC <- gbayes(stat = stat, Glist = Glist, method = "bayesC", nit = 1000))


# Compute sumstats
y <- scale(sim$y)
cls <- match(Glist$rsidsLD[[chr]], Glist$rsids[[chr]])
chr <- 1
system.time(stats <- bed_xtx_xty(
 bed_file = Glist$bedfiles[chr],
 n = Glist$n,
 cls = cls,
 af = Glist$af[[chr]][cls],
 y = y,
 scale = TRUE,
 nthreads = 4
))

# Compute sparse ld
cls <- match(Glist$rsidsLD[[1]],(Glist$rsids[[1]]))
system.time(out <- sparseLD_stream_CSR(
 bed_files = Glist$bedfiles[1],
 n = Glist$n,
 cls = list(cls),
 out_prefix = file.path(
  data_dir,
  "ld_test"
 ),
 rows = NULL,
 af = list(Glist$af[[1]][cls]),
 pos_bp = list(Glist$pos[[1]][cls]),
 max_distance_bp = 0,            # disables bp-distance filtering
 max_distance_variants = 1000,   # local LD window; 0 disables this filter
 r2_threshold = 0.0001,
 block_size = 1024,
 nthreads = 1
))


ld <- sparseLD_read_CSR(file.path(
 data_dir,
 "ld_test"), one_based = TRUE)


m <- length(cls)
n <- Glist$n
nt <- length(stats$yy)

b <- lapply(seq_len(nt), function(x) rep(0, m))
d <- lapply(seq_len(nt), function(x) rep(0, m))

nub <- nue <- 4

vy <- stats$yy / (n - 1)
h2 <- 0.5
pi_marker_global <- 0.001

vb <- diag((vy * h2) / (m * pi_marker_global))
ve <- diag(vy * (1 - h2))

ssb_prior <- diag(((nub - 2) / nub) * (vy * h2) / (m * pi_marker_global))
sse_prior <- diag(((nue - 2) / nue) * (vy * (1 - h2)))

ssb_prior_list <- split(
 ssb_prior,
 rep(seq_len(ncol(ssb_prior)), each = nrow(ssb_prior))
)

sse_prior_list <- split(
 sse_prior,
 rep(seq_len(ncol(sse_prior)), each = nrow(sse_prior))
)

stopifnot(length(sim$pi_marker) == nt)
stopifnot(all(lengths(sim$pi_marker) == m))
stopifnot(length(sim$vb_multiplier) == nt)
stopifnot(all(lengths(sim$vb_multiplier) == m))

Sys.setenv(
 MKL_NUM_THREADS = "1",
 MKL_DYNAMIC = "FALSE",
 OMP_NUM_THREADS = "3",
 OMP_DYNAMIC = "FALSE"
)

system.time(fitSTA <- stblr_cpg_omp_csr_prior(
 wy = stats$wy,
 ww = stats$ww,
 yy = stats$yy,
 b_init = b,
 d_init = d,
 use_d_init = FALSE,
 r_init = stats$wy,
 use_r_init = FALSE,
 rebuild_r_before_updateE = FALSE,

 ld_prefix = file.path(
  data_dir,
  "ld_test"
 ),

 B = vb,
 E = ve,
 ssb_prior = ssb_prior_list,
 sse_prior = sse_prior_list,

 pi = c(1 - pi_marker_global, pi_marker_global),

 use_pi_marker = TRUE,
 pi_marker = sim$pi_marker,

 use_vb_multiplier = TRUE,
 vb_multiplier = sim$vb_multiplier,

 nub = nub,
 nue = nue,
 updateB = TRUE,
 updateE = TRUE,
 updatePi = FALSE,

 adjE = 0.9,
 n = rep(as.integer(n), nt),

 nit = 1000,
 nburn = 100,
 nthin = 1,

 ncores = 3,
 seed = 10
))

system.time(fitST <- stblr_cpg_omp_csr(
 wy = stats$wy,
 ww = stats$ww,
 yy = stats$yy,
 b_init = b,
 d_init = d,
 use_d_init = FALSE,
 r_init = stats$wy,
 use_r_init = FALSE,
 rebuild_r_before_updateE = FALSE,
 #ld_row_ptr = ld$row_ptr,
 #ld_col_idx = ld$col_idx,
 #ld_values = ld$values,
 #ld_col_idx_one_based = TRUE,
 ld_prefix = file.path(
  data_dir,
  "ld_test"
 ),
 B = vb,
 E = ve,
 ssb_prior = split(ssb_prior, rep(1:ncol(ssb_prior), each = nrow(ssb_prior))),
 sse_prior = split(sse_prior, rep(1:ncol(sse_prior), each = nrow(sse_prior))),
 pi = c(1 - pi, pi),
 nub = 4,
 nue = 4,
 updateB = TRUE,
 updateE = TRUE,
 updatePi = TRUE,
 adjE = 0.9,
 n = rep(as.integer(n), nt),
 nit = 1000,
 nburn = 100,
 nthin = 1,
 ncores=3,
 seed = 10
))


fitST <- format_stblr_fit(fitST, nt=3, m=m)
fitSTA <- format_stblr_fit(fitSTA, nt=3, m=m)

summarize_annotation_signal(sim, fitST)
summarize_annotation_signal(sim, fitSTA)

recST <- recovery_summary(
 fit = fitST,
 Btrue = sim$B,
 Ks = c(40, 60, 100, 150, 250),
 trait_names = paste0("D", 1:3)
)

recSTA <- recovery_summary(
 fit = fitSTA,
 Btrue = sim$B,
 Ks = c(40, 60, 100, 150, 250),
 trait_names = paste0("D", 1:3)
)

recST$model <- "ST-BLR"
recSTA$model <- "ST-BLR annotation prior"

rec <- rbind(recST, recSTA)
rec


rec_summary <- rec %>%
 group_by(model, K) %>%
 summarise(
  mean_recovered_dm = mean(recovered_dm),
  mean_precision_dm = mean(precision_dm),
  mean_recall_dm = mean(recall_dm),

  mean_recovered_bm = mean(recovered_bm),
  mean_precision_bm = mean(precision_bm),
  mean_recall_bm = mean(recall_bm),

  .groups = "drop"
 )

rec_summary


rec_summary_delta <- rec_summary %>%
 tidyr::pivot_wider(
  names_from = model,
  values_from = c(
   mean_recovered_dm,
   mean_precision_dm,
   mean_recall_dm,
   mean_recovered_bm,
   mean_precision_bm,
   mean_recall_bm
  )
 ) %>%
 mutate(
  delta_mean_recovered_dm =
   `mean_recovered_dm_ST-BLR annotation prior` -
   `mean_recovered_dm_ST-BLR`,

  delta_mean_recall_dm =
   `mean_recall_dm_ST-BLR annotation prior` -
   `mean_recall_dm_ST-BLR`,

  delta_mean_recovered_bm =
   `mean_recovered_bm_ST-BLR annotation prior` -
   `mean_recovered_bm_ST-BLR`,

  delta_mean_recall_bm =
   `mean_recall_bm_ST-BLR annotation prior` -
   `mean_recall_bm_ST-BLR`
 )

rec_summary_delta

rank_similarity <- lapply(seq_len(ncol(sim$B)), function(t) {
 cor(
  fitST$dm[, t],
  abs(fitST$bm[, t]),
  method = "spearman"
 )
})

unlist(rank_similarity)

rank_similarity <- lapply(seq_len(ncol(sim$B)), function(t) {
 cor(
  fitSTA$dm[, t],
  abs(fitSTA$bm[, t]),
  method = "spearman"
 )
})

unlist(rank_similarity)



pri <- derive_marker_priors_from_fit(
 fit = fitST,
 anno = sim$annot,
 base_pi = rep(pi_marker_global, nt),
 shrink = 0.5,
 pi_cap = 0.05,
 vb_floor = 0.25,
 vb_cap = 4,
 return_list = TRUE
)

system.time(fitSTA_EB <- stblr_cpg_omp_csr_prior(
 wy = stats$wy,
 ww = stats$ww,
 yy = stats$yy,
 b_init = b,
 d_init = d,
 use_d_init = FALSE,
 r_init = stats$wy,
 use_r_init = FALSE,
 rebuild_r_before_updateE = FALSE,

 ld_prefix = file.path(
  data_dir,
  "ld_test"
 ),

 B = vb,
 E = ve,
 ssb_prior = ssb_prior_list,
 sse_prior = sse_prior_list,

 pi = c(1 - pi_marker_global, pi_marker_global),

 use_pi_marker = TRUE,
 pi_marker = pri$pi_marker,

 use_vb_multiplier = TRUE,
 vb_multiplier = pri$vb_multiplier,

 nub = nub,
 nue = nue,
 updateB = TRUE,
 updateE = TRUE,
 updatePi = FALSE,

 adjE = 0.9,
 n = rep(as.integer(n), nt),

 nit = 1000,
 nburn = 100,
 nthin = 1,

 ncores = 3,
 seed = 10
))


fitSTA_EB <- format_stblr_fit(fitSTA_EB, nt=3, m=m)

summarize_annotation_signal(sim, fitSTA_EB)

recSTA_EB <- recovery_summary(
 fit = fitSTA_EB,
 Btrue = sim$B,
 Ks = c(40, 60, 100, 150, 250),
 trait_names = paste0("D", 1:3)
)




K  <- ncol(sim$annot)
nt <- length(stats$yy)
m  <- length(stats$wy[[1]])

eta_pi_init <- matrix(0, nrow = K, ncol = nt)
eta_vb_init <- matrix(0, nrow = K, ncol = nt)

rownames(eta_pi_init) <- colnames(sim$annot)
rownames(eta_vb_init) <- colnames(sim$annot)
colnames(eta_pi_init) <- paste0("D", seq_len(nt))
colnames(eta_vb_init) <- paste0("D", seq_len(nt))

stopifnot(nrow(sim$annot) == m)
stopifnot(ncol(eta_pi_init) == nt)
stopifnot(ncol(eta_vb_init) == nt)


system.time(fitSTA_annot <- stblr_cpg_omp_csr_annot(
 wy = stats$wy,
 ww = stats$ww,
 yy = stats$yy,
 b_init = b,
 d_init = d,
 use_d_init = FALSE,
 r_init = stats$wy,
 use_r_init = FALSE,
 rebuild_r_before_updateE = FALSE,

 ld_prefix = file.path(
  data_dir,
  "ld_test"
 ),

 B = vb,
 E = ve,
 ssb_prior = ssb_prior_list,
 sse_prior = sse_prior_list,

 pi = c(1 - pi_marker_global, pi_marker_global),

 A = sim$annot,

 learn_pi_annot = TRUE,
 learn_vb_annot = FALSE,

 eta_pi_init = eta_pi_init,
 eta_vb_init = eta_vb_init,

 sigma_eta_pi = 1.0,
 sigma_eta_vb = 0.5,

 rw_sd_eta_pi = 0.2,
 rw_sd_eta_vb = 0.2,

 annot_update_every = 10,

 pi_min = 1e-8,
 pi_max = 0.05,

 vb_multiplier_min = 0.25,
 vb_multiplier_max = 4.0,

 nub = nub,
 nue = nue,
 updateB = TRUE,
 updateE = TRUE,
 updatePi = FALSE,

 adjE = 0.9,
 n = rep(as.integer(n), nt),

 nit = 1000,
 nburn = 100,
 nthin = 1,

 ncores = 3,
 seed = 10
))

fitSTA_annot <- format_stblr_fit(
 fitSTA_annot,
 nt = nt,
 m = m,
 annotation_names = colnames(sim$annot)
)

colSums(fitSTA_annot$dm)
colMeans(fitSTA_annot$vbs)
colMeans(fitSTA_annot$vgs)
colMeans(fitSTA_annot$ves)

recSTA_annot <- recovery_summary(
 fit = fitSTA_annot,
 Btrue = sim$B,
 Ks = c(40, 60, 100, 150, 250),
 trait_names = paste0("D", seq_len(nt))
)

recSTA_annot


# If eta outputs are available from the C++ result:
fitSTA_annot$eta_pi
fitSTA_annot$eta_vb

# Annotation sizes and causal enrichment
colSums(sim$annot)

causal_any <- rowSums(sim$B != 0) > 0
apply(sim$annot[causal_any, , drop = FALSE], 2, mean)
apply(sim$annot[!causal_any, , drop = FALSE], 2, mean)


annot_enrichment <-
 apply(sim$annot[causal_any, , drop = FALSE], 2, mean) /
 apply(sim$annot[!causal_any, , drop = FALSE], 2, mean)

round(annot_enrichment, 2)

annotation_diagnostics <- function(fit, annot, Btrue) {
 causal_any <- rowSums(Btrue != 0) > 0

 true_enrichment <-
  apply(annot[causal_any, , drop = FALSE], 2, mean) /
  apply(annot[!causal_any, , drop = FALSE], 2, mean)

 out <- list(
  true_enrichment = true_enrichment,
  eta_pi = fit$eta_pi,
  eta_vb = fit$eta_vb,
  mean_eta_pi = colMeans(fit$eta_pi),
  mean_eta_vb = colMeans(fit$eta_vb)
 )

 out
}

diag_annot <- annotation_diagnostics(
 fit = fitSTA_annot,
 annot = sim$annot,
 Btrue = sim$B
)

round(diag_annot$true_enrichment, 2)
round(diag_annot$eta_pi, 3)
round(diag_annot$eta_vb, 3)
round(diag_annot$mean_eta_pi, 3)


make_pi_from_eta_R <- function(A, eta, base_pi, pi_min = 1e-8, pi_max = 0.05) {
 eta_marker <- A %*% eta
 eta_marker <- sweep(eta_marker, 2, colMeans(eta_marker), "-")
 pi_marker <- plogis(sweep(eta_marker, 2, qlogis(base_pi), "+"))
 pmin(pmax(pi_marker, pi_min), pi_max)
}

pi_learned <- make_pi_from_eta_R(
 A = sim$annot,
 eta = t(fitSTA_annot$eta_pi),
 base_pi = pi_marker_global,
 pi_min = 1e-8,
 pi_max = 0.05
)

apply(pi_learned, 2, summary)

for (a in colnames(sim$annot)) {
 idx <- sim$annot[, a] != 0
 cat("\n", a, "\n")
 print(rbind(
  in_set  = colMeans(pi_learned[idx, , drop = FALSE]),
  out_set = colMeans(pi_learned[!idx, , drop = FALSE])
 ))
}


trait_annotation_enrichment <- function(annot, Btrue) {
 nt <- ncol(Btrue)
 out <- matrix(NA_real_, nrow = ncol(annot), ncol = nt)
 rownames(out) <- colnames(annot)
 colnames(out) <- colnames(Btrue)

 for (t in seq_len(nt)) {
  causal_t <- Btrue[, t] != 0

  out[, t] <-
   apply(annot[causal_t, , drop = FALSE], 2, mean) /
   apply(annot[!causal_t, , drop = FALSE], 2, mean)
 }

 out
}

round(trait_annotation_enrichment(sim$annot, sim$B), 2)

make_vb_multiplier_from_eta_R <- function(A, eta, mult_min = 0.25, mult_max = 4) {
 A <- as.matrix(A)
 eta <- as.matrix(eta)

 if (ncol(A) != nrow(eta)) {
  stop("ncol(A) must equal nrow(eta). eta should be K x nt.")
 }

 z <- A %*% eta

 # Center annotation contribution trait-wise, matching the C++ implementation
 z <- sweep(z, 2, colMeans(z, na.rm = TRUE), "-")

 z <- pmin(pmax(z, log(mult_min)), log(mult_max))
 mult <- exp(z)

 pmin(pmax(mult, mult_min), mult_max)
}

vb_mult_learned <- make_vb_multiplier_from_eta_R(
 A = sim$annot,
 eta = t(fitSTA_annot$eta_vb),
 mult_min = 0.25,
 mult_max = 4
)

for (a in colnames(sim$annot)) {
 idx <- sim$annot[, a] != 0
 cat("\n", a, "\n")
 print(rbind(
  in_set  = colMeans(vb_mult_learned[idx, , drop = FALSE]),
  out_set = colMeans(vb_mult_learned[!idx, , drop = FALSE]),
  ratio   = colMeans(vb_mult_learned[idx, , drop = FALSE]) /
   colMeans(vb_mult_learned[!idx, , drop = FALSE])
 ))
}



# =============================================================================
# R-side helper snippets
# =============================================================================
# The helper definitions below are executable R code.
 # Convert overlapping binary annotations to one primary class for this comparison.
 # This is only one possible encoding. For SBayesRC-like comparison, you need
 # one annotation_class per marker.

 make_primary_annotation_class <- function(A, priority = seq_len(ncol(A))) {
  A <- as.matrix(A)
  cls <- integer(nrow(A))

  # class 0 = background / no selected annotation
  # classes 1..K = first active annotation according to priority
  for (i in seq_len(nrow(A))) {
   hit <- which(A[i, priority] != 0)
   if (length(hit) > 0) {
    cls[i] <- hit[1]
   } else {
    cls[i] <- 0L
   }
  }

  cls
 }

format_stblr_fit_sbayesrc <- function(
  fit,
  nt,
  m,
  n_classes,
  mixture_var,
  trait_names = NULL,
  variable_names = NULL,
  class_names = NULL
) {
 names(fit) <- c(
  "bm", "dm", "wy", "r", "b", "d", "o",
  "vbs", "vgs", "ves",
  "covb", "covg", "cove",
  "vb", "vg", "ve",
  "pi", "pim", "pi_class_mean", "pi_class_final"
 )

 if (is.null(trait_names)) trait_names <- paste0("T", seq_len(nt))
 if (is.null(variable_names)) variable_names <- paste0("V", seq_len(m))
 if (is.null(class_names)) class_names <- paste0("C", seq_len(n_classes) - 1L)

 Kmix <- length(mixture_var)
 mix_names <- paste0("M", seq_len(Kmix) - 1L, "_v", mixture_var)

 for (i in 1:7) {
  fit[[i]] <- as.matrix(as.data.frame(fit[[i]]))
  rownames(fit[[i]]) <- variable_names
  colnames(fit[[i]]) <- trait_names
 }

 for (i in 8:10) {
  fit[[i]] <- as.matrix(as.data.frame(fit[[i]]))
  rownames(fit[[i]]) <- paste0("Iter", seq_len(nrow(fit[[i]])))
  colnames(fit[[i]]) <- trait_names
 }

 for (i in 11:16) {
  fit[[i]] <- matrix(unlist(fit[[i]]), ncol = nt, byrow = TRUE)
  colnames(fit[[i]]) <- rownames(fit[[i]]) <- trait_names
 }

 fit[[17]] <- do.call(rbind, fit[[17]])
 fit[[18]] <- do.call(rbind, fit[[18]])
 rownames(fit[[17]]) <- trait_names
 rownames(fit[[18]]) <- trait_names
 colnames(fit[[17]]) <- c("pi0", "pi_nonnull")
 colnames(fit[[18]]) <- c("pi0", "pi_nonnull")

 raw_mean <- fit[[19]]
 raw_final <- fit[[20]]

 pi_class_mean <- vector("list", nt)
 pi_class_final <- vector("list", nt)
 names(pi_class_mean) <- names(pi_class_final) <- trait_names

 for (t in seq_len(nt)) {
  pi_class_mean[[t]] <- matrix(unlist(raw_mean[[t]]), nrow = n_classes, ncol = Kmix, byrow = TRUE)
  pi_class_final[[t]] <- matrix(unlist(raw_final[[t]]), nrow = n_classes, ncol = Kmix, byrow = TRUE)

  rownames(pi_class_mean[[t]]) <- rownames(pi_class_final[[t]]) <- class_names
  colnames(pi_class_mean[[t]]) <- colnames(pi_class_final[[t]]) <- mix_names
 }

 fit <- fit[1:18]
 fit$pi_class_mean <- pi_class_mean
 fit$pi_class_final <- pi_class_final

 if (sum(diag(fit$covb)) > 0) fit$rb <- cov2cor(fit$covb)
 if (sum(diag(fit$covg)) > 0) fit$rg <- cov2cor(fit$covg)
 if (sum(diag(fit$cove)) > 0) fit$re <- cov2cor(fit$cove)

 fit
}

# Example usage:
#
# annotation_class <- make_primary_annotation_class(sim$annot)
# n_classes <- max(annotation_class) + 1L
#
# mixture_var <- c(0, 0.01, 0.1, 1.0)
# Kmix <- length(mixture_var)
#
# pi0 <- c(0.999, 0.0007, 0.0002, 0.0001)
# pi0 <- pi0 / sum(pi0)
#
# pi_class_init <- matrix(rep(pi0, each = n_classes), nrow = n_classes, byrow = FALSE)
# alpha_class <- matrix(1, nrow = n_classes, ncol = Kmix)
#
# gamma_init <- lapply(seq_len(nt), function(t) rep(0, m))
#
# fitSRC <- stblr_cpg_omp_csr_sbayesrc(
#   wy = stats$wy,
#   ww = stats$ww,
#   yy = stats$yy,
#   b_init = b,
#   gamma_init = gamma_init,
#   use_gamma_init = FALSE,
#   r_init = stats$wy,
#   use_r_init = FALSE,
#   rebuild_r_before_updateE = FALSE,
#   ld_prefix = file.path(data_dir, "ld_test"),
#   B = vb,
#   E = ve,
#   ssb_prior = ssb_prior_list,
#   sse_prior = sse_prior_list,
#   annotation_class = annotation_class,
#   n_classes = n_classes,
#   mixture_var = mixture_var,
#   pi_class_init = pi_class_init,
#   alpha_class = alpha_class,
#   updatePiClass = TRUE,
#   nub = nub,
#   nue = nue,
#   updateB = TRUE,
#   updateE = TRUE,
#   adjE = 0.9,
#   n = rep(as.integer(n), nt),
#   nit = 1000,
#   nburn = 100,
#   nthin = 1,
#   ncores = 3,
#   seed = 10
# )
#
# fitSRC <- format_stblr_fit_sbayesrc(
#   fitSRC,
#   nt = nt,
#   m = m,
#   n_classes = n_classes,
#   mixture_var = mixture_var,
#   trait_names = paste0("T", seq_len(nt)),
#   variable_names = paste0("V", seq_len(m)),
#   class_names = c("background", colnames(sim$annot))
# )
#
# fitSRC$pi_class_mean


A <- cbind(intercept = 1, sim$annot)
A <- as.matrix(A)

K  <- ncol(A)
nt <- length(stats$yy)
m  <- nrow(A)

gamma  <- c(0, 0.01, 0.1, 1)
ndist  <- length(gamma)

startPi <- c(0.99, 0.006, 0.003, 0.001)

# Stick-breaking conditional probabilities:
# pi1 = 1 - p1
# pi2 = p1 * (1 - p2)
# pi3 = p1 * p2 * (1 - p3)
# pi4 = p1 * p2 * p3
p1 <- 1 - startPi[1]
p2 <- 1 - startPi[2] / p1
p3 <- startPi[4] / (p1 * p2)

alpha_init <- matrix(0, nrow = K, ncol = ndist - 1)
rownames(alpha_init) <- colnames(A)
colnames(alpha_init) <- paste0("p", 2:ndist)

alpha_init[1, ] <- qnorm(c(p1, p2, p3))
p <- pnorm(alpha_init[1, ])

pi_check <- c(
 1 - p[1],
 p[1] * (1 - p[2]),
 p[1] * p[2] * (1 - p[3]),
 p[1] * p[2] * p[3]
)

round(pi_check, 6)
round(startPi, 6)
sum(pi_check)

sigmaSqAlpha_init <- rep(1, ndist - 1)

b    <- lapply(seq_len(nt), function(t) rep(0, m))
comp <- lapply(seq_len(nt), function(t) rep(0, m))

nub <- 4
nue <- 4

n  <- Glist$n
vy <- stats$yy / (n - 1)
h2 <- 0.5

mean_gamma_active <- sum(gamma * startPi)

vb <- diag((vy * h2) / (m * mean_gamma_active))
ve <- diag(vy * (1 - h2))

ssb_prior <- diag(((nub - 2) / nub) * (vy * h2) / (m * mean_gamma_active))
sse_prior <- diag(((nue - 2) / nue) * (vy * (1 - h2)))

ssb_prior_list <- split(
 ssb_prior,
 rep(seq_len(ncol(ssb_prior)), each = nrow(ssb_prior))
)

sse_prior_list <- split(
 sse_prior,
 rep(seq_len(ncol(sse_prior)), each = nrow(sse_prior))
)

Sys.setenv(
 MKL_NUM_THREADS = "1",
 MKL_DYNAMIC = "FALSE",
 OMP_NUM_THREADS = "3",
 OMP_DYNAMIC = "FALSE"
)

fitSRC <- stblr_cpg_omp_csr_sbayesrc(
 wy = stats$wy,
 ww = stats$ww,
 yy = stats$yy,

 b_init = b,
 comp_init = comp,
 use_comp_init = FALSE,

 r_init = stats$wy,
 use_r_init = FALSE,
 rebuild_r_before_updateE = FALSE,

 ld_prefix = file.path(
  data_dir,
  "ld_test"
 ),

 B = vb,
 E = ve,
 ssb_prior = ssb_prior_list,
 sse_prior = sse_prior_list,

 A = A,
 gamma = c(0, 0.001, 0.01, 0.1),

 alpha_init = alpha_init,
 sigmaSqAlpha_init = sigmaSqAlpha_init,

 intercept_flat = TRUE,
 sigmaSqAlpha_a = 2,
 sigmaSqAlpha_b = 2,

 pi_floor = 1e-12,

 nub = nub,
 nue = nue,

 updateAlpha = TRUE,
 updateB = TRUE,
 updateE = TRUE,

 alpha_update_every = 1,

 adjE = 0.9,
 n = rep(as.integer(n), nt),

 nit = 100,
 nburn = 10,
 nthin = 1,

 ncores = 3,
 seed = 10
)

fitSRC <- format_stblr_fit(
 fitSRC,
 nt = nt,
 m = m,
 trait_names = paste0("T", seq_len(nt)),
 variable_names = rownames(sim$B),
 annotation_names = colnames(A),
 gamma = c(0, 0.001, 0.01, 0.1)
)
colSums(fitSRC$dm)
colMeans(fitSRC$vbs)
colMeans(fitSRC$vgs)
colMeans(fitSRC$ves)

recSRC <- recovery_summary(
 fit = fitSRC,
 Btrue = sim$B,
 Ks = c(40, 60, 100, 150, 250),
 trait_names = paste0("D", seq_len(nt))
)

recSRC


fitSRC$alpha$T1
fitSRC$alpha$T2
fitSRC$alpha$T3

fitSRC$sigmaSqAlpha
colSums(fitSRC$dm)
colMeans(fitSRC$vbs)
colMeans(fitSRC$vgs)
colMeans(fitSRC$ves)


sbayesrc_annotation_pi <- function(alpha, gamma = c(0, 0.001, 0.01, 0.1)) {
 ndist <- length(gamma)

 p <- pnorm(alpha)

 pi <- matrix(
  NA_real_,
  nrow = nrow(alpha),
  ncol = ndist,
  dimnames = list(rownames(alpha), paste0("Pi", seq_len(ndist)))
 )

 pi[, 1] <- 1 - p[, 1]
 pi[, 2] <- p[, 1] * (1 - p[, 2])
 pi[, 3] <- p[, 1] * p[, 2] * (1 - p[, 3])
 pi[, 4] <- p[, 1] * p[, 2] * p[, 3]

 pi
}

sbayesrc_annotation_gamma_mean <- function(alpha, gamma = c(0, 0.001, 0.01, 0.1)) {
 pi <- sbayesrc_annotation_pi(alpha, gamma)
 as.numeric(pi %*% gamma)
}

round(sbayesrc_annotation_pi(fitSRC$alpha$T1), 4)
round(sbayesrc_annotation_pi(fitSRC$alpha$T2), 4)
round(sbayesrc_annotation_pi(fitSRC$alpha$T3), 4)

data.frame(
 annotation = rownames(fitSRC$alpha$T1),
 T1 = sbayesrc_annotation_gamma_mean(fitSRC$alpha$T1),
 T2 = sbayesrc_annotation_gamma_mean(fitSRC$alpha$T2),
 T3 = sbayesrc_annotation_gamma_mean(fitSRC$alpha$T3)
)

sbayesrc_marker_pi <- function(A, alpha, gamma = c(0, 0.001, 0.01, 0.1)) {
 eta <- A %*% alpha
 p <- pnorm(eta)

 pi <- matrix(NA_real_, nrow = nrow(A), ncol = length(gamma))
 colnames(pi) <- paste0("Pi", seq_along(gamma))

 pi[, 1] <- 1 - p[, 1]
 pi[, 2] <- p[, 1] * (1 - p[, 2])
 pi[, 3] <- p[, 1] * p[, 2] * (1 - p[, 3])
 pi[, 4] <- p[, 1] * p[, 2] * p[, 3]

 pi
}

sbayesrc_marker_gamma_mean <- function(A, alpha, gamma = c(0, 0.001, 0.01, 0.1)) {
 pi <- sbayesrc_marker_pi(A, alpha, gamma)
 as.numeric(pi %*% gamma)
}

marker_gmean_T1 <- sbayesrc_marker_gamma_mean(A, fitSRC$alpha$T1, gamma)
marker_gmean_T2 <- sbayesrc_marker_gamma_mean(A, fitSRC$alpha$T2, gamma)
marker_gmean_T3 <- sbayesrc_marker_gamma_mean(A, fitSRC$alpha$T3, gamma)

summary(marker_gmean_T1)
summary(marker_gmean_T2)
summary(marker_gmean_T3)

causal_any <- rowSums(sim$B != 0) > 0

rbind(
 causal = c(
  mean(marker_gmean_T1[causal_any]),
  mean(marker_gmean_T2[causal_any]),
  mean(marker_gmean_T3[causal_any])
 ),
 noncausal = c(
  mean(marker_gmean_T1[!causal_any]),
  mean(marker_gmean_T2[!causal_any]),
  mean(marker_gmean_T3[!causal_any])
 )
)

prior_architecture_diagnostics <- function(A, fit, Btrue, gamma) {
 causal_any <- rowSums(Btrue != 0) > 0
 traits <- names(fit$alpha)

 out <- do.call(rbind, lapply(traits, function(tr) {
  gmean <- sbayesrc_marker_gamma_mean(A, fit$alpha[[tr]], gamma)

  data.frame(
   trait = tr,
   mean_causal = mean(gmean[causal_any]),
   mean_noncausal = mean(gmean[!causal_any]),
   enrichment = mean(gmean[causal_any]) / mean(gmean[!causal_any]),
   median_causal = median(gmean[causal_any]),
   median_noncausal = median(gmean[!causal_any]),
   max = max(gmean)
  )
 }))

 rownames(out) <- NULL
 out
}

prior_architecture_diagnostics(
 A = A,
 fit = fitSRC,
 Btrue = sim$B,
 gamma = gamma
)

gmean <- data.frame(
 marker = rownames(sim$B),
 causal_any = rowSums(sim$B != 0) > 0,
 T1_prior = sbayesrc_marker_gamma_mean(A, fitSRC$alpha$T1, gamma),
 T2_prior = sbayesrc_marker_gamma_mean(A, fitSRC$alpha$T2, gamma),
 T3_prior = sbayesrc_marker_gamma_mean(A, fitSRC$alpha$T3, gamma),
 T1_pip = fitSRC$dm[, "T1"],
 T2_pip = fitSRC$dm[, "T2"],
 T3_pip = fitSRC$dm[, "T3"]
)

cor(gmean$T1_prior, gmean$T1_pip, method = "spearman")
cor(gmean$T2_prior, gmean$T2_pip, method = "spearman")
cor(gmean$T3_prior, gmean$T3_pip, method = "spearman")

plot(gmean$T1_prior, gmean$T1_pip,
     pch = 16, cex = 0.5,
     xlab = "Annotation-derived prior E(gamma)",
     ylab = "Posterior inclusion probability",
     main = "T1: SBayesRC prior vs posterior evidence")

points(gmean$T1_prior[gmean$causal_any],
       gmean$T1_pip[gmean$causal_any],
       pch = 16, col = "red")
