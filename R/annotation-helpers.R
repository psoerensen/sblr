`%||%` <- function(x, y) {
 if (is.null(x)) y else x
}

.stblr_match_annotation_model <- function(model) {
 if (length(model) != 1L || is.na(model)) {
  stop("annotation_model must be a single string.")
 }

 key <- tolower(gsub("-", "_", model))
 aliases <- c(
  prior = "prior",
  fixed_prior = "prior",
  fixed = "prior",
  learned = "learned",
  annot = "learned",
  annotation = "learned",
  group = "group",
  groups = "group",
  sbayesrc = "sbayesrc"
 )

 out <- unname(aliases[[key]])
 if (is.null(out)) {
  stop(
   "annotation_model must be one of prior, learned, group, or sbayesrc.",
   call. = FALSE
  )
 }
 out
}

.stblr_validate_annotation_chain_args <- function(nchains, keep_chains, chain_seeds) {
 if (!is.numeric(nchains) || length(nchains) != 1L ||
     !is.finite(nchains) || nchains < 1 || nchains != floor(nchains)) {
  stop("nchains must be a positive integer scalar.", call. = FALSE)
 }
 nchains <- as.integer(nchains)
 if (!is.logical(keep_chains) || length(keep_chains) != 1L ||
     is.na(keep_chains)) {
  stop("keep_chains must be TRUE or FALSE.", call. = FALSE)
 }
 if (!is.null(chain_seeds)) {
  if (!is.numeric(chain_seeds) || length(chain_seeds) != nchains ||
      anyNA(chain_seeds) || any(!is.finite(chain_seeds)) ||
      any(chain_seeds != floor(chain_seeds))) {
   stop("chain_seeds must be NULL or an integer/numeric vector of length nchains.",
        call. = FALSE)
  }
  chain_seeds <- as.integer(chain_seeds)
 } else {
  chain_seeds <- integer()
 }
 list(nchains = nchains, keep_chains = keep_chains, chain_seeds = chain_seeds)
}

.stblr_stop_bayesc_annotation_ld_swap <- function(updateLDswap) {
 if (isTRUE(updateLDswap)) {
  stop(
   "LD-swap/MH is currently implemented only for annotation_model = 'prior', 'group', and 'learned' among annotation-aware CSR models.",
   call. = FALSE
  )
 }
 invisible(TRUE)
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

 if (length(stats$wy) != nt || length(stats$ww) != nt ||
     length(stats$yy) != nt) {
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
  if (!is.numeric(val) || length(val) != 1 || !is.finite(val) ||
      val <= 0 || val >= 1) {
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
  ssb_prior_list = split(
   ssb_prior,
   rep(seq_len(ncol(ssb_prior)), each = nrow(ssb_prior))
  ),
  sse_prior_list = split(
   sse_prior,
   rep(seq_len(ncol(sse_prior)), each = nrow(sse_prior))
  )
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

 if (ncol(beta_pi) == 1 && nt > 1) {
  beta_pi <- beta_pi[, rep(1, nt), drop = FALSE]
 }
 if (ncol(beta_vb) == 1 && nt > 1) {
  beta_vb <- beta_vb[, rep(1, nt), drop = FALSE]
 }

 if (!all(dim(beta_pi) == c(K, nt))) {
  stop("beta_pi must be K x nt or length K.")
 }
 if (!all(dim(beta_vb) == c(K, nt))) {
  stop("beta_vb must be K x nt or length K.")
 }

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

 if (length(pi_group_prior_a) == 1) {
  pi_group_prior_a <- rep(pi_group_prior_a, ngroup)
 }
 if (length(pi_group_prior_b) == 1) {
  pi_group_prior_b <- rep(pi_group_prior_b, ngroup)
 }

 if (length(pi_group_prior_a) != ngroup ||
     length(pi_group_prior_b) != ngroup ||
     any(!is.finite(pi_group_prior_a)) ||
     any(!is.finite(pi_group_prior_b)) ||
     any(pi_group_prior_a <= 0) ||
     any(pi_group_prior_b <= 0)) {
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
  group_vb_multiplier_init = split(
   group_vb_multiplier_init_mat,
   row(group_vb_multiplier_init_mat)
  ),
  group_pi_init_matrix = group_pi_init_mat,
  group_vb_multiplier_init_matrix = group_vb_multiplier_init_mat,
  pi_group_prior_a = pi_group_prior_a,
  pi_group_prior_b = pi_group_prior_b,
  pi_group_prior_mean = pi_group_prior_mean,
  pi_group_prior_strength = pi_group_prior_strength
 )
}
