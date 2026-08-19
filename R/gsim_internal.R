# Internal helpers for gsim().
#
# These functions deliberately use base R only. qgg is needed only when gsim()
# is called with a Glist and no custom getG_fun is supplied.

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

.gsim_stop <- function(...) {
  stop(..., call. = FALSE)
}

.gsim_check_scalar <- function(x, name, lower = -Inf, upper = Inf,
                               integer = FALSE) {
  if (length(x) != 1L || !is.finite(x) || x < lower || x > upper ||
      (integer && x != as.integer(x))) {
    .gsim_stop(name, " must be a finite ",
               if (integer) "integer " else "",
               "scalar in [", lower, ", ", upper, "].")
  }
  invisible(TRUE)
}

.gsim_as_named_vector <- function(x, ids, name, allow_null = TRUE,
                                  numeric_only = TRUE) {
  if (is.null(x)) {
    if (allow_null) return(NULL)
    .gsim_stop(name, " must not be NULL.")
  }

  x <- unlist(x, use.names = TRUE)
  if (numeric_only && !is.numeric(x)) .gsim_stop(name, " must be numeric.")

  if (!is.null(names(x)) && all(nzchar(names(x)))) {
    if (anyDuplicated(names(x))) .gsim_stop(name, " has duplicated names.")
    pos <- match(ids, names(x))
    if (anyNA(pos)) .gsim_stop(name, " does not contain every requested marker.")
    x <- x[pos]
  } else if (length(x) != length(ids)) {
    .gsim_stop(name, " must have length ", length(ids),
               " or be named by marker ID.")
  }

  names(x) <- ids
  x
}

.gsim_marker_ids_from_glist <- function(Glist, rsids = NULL) {
  if (!is.null(rsids)) {
    marker_ids <- as.character(unlist(rsids, use.names = FALSE))
  } else if (!is.null(Glist$rsidsLD)) {
    marker_ids <- as.character(unlist(Glist$rsidsLD, use.names = FALSE))
  } else if (!is.null(Glist$rsids)) {
    marker_ids <- as.character(unlist(Glist$rsids, use.names = FALSE))
  } else {
    .gsim_stop("Glist must contain rsidsLD or rsids, or rsids must be supplied.")
  }

  if (!length(marker_ids) || anyNA(marker_ids) || any(!nzchar(marker_ids))) {
    .gsim_stop("Marker IDs must be non-empty and non-missing.")
  }
  if (anyDuplicated(marker_ids)) .gsim_stop("Marker IDs must be unique.")
  marker_ids
}

.gsim_chr_map_from_glist <- function(Glist, marker_ids) {
  source <- Glist$rsidsLD %||% Glist$rsids
  if (is.null(source) || !is.list(source)) {
    return(setNames(rep(NA_character_, length(marker_ids)), marker_ids))
  }

  chr_names <- names(source)
  if (is.null(chr_names) || any(!nzchar(chr_names))) {
    chr_names <- as.character(seq_along(source))
  }

  out <- rep(NA_character_, length(marker_ids))
  names(out) <- marker_ids
  for (i in seq_along(source)) {
    hit <- intersect(marker_ids, as.character(source[[i]]))
    out[hit] <- chr_names[i]
  }
  out
}

.gsim_prepare_W <- function(W, ids = NULL) {
  W <- as.matrix(W)
  storage.mode(W) <- "double"
  if (!nrow(W) || !ncol(W)) .gsim_stop("W must have at least one row and column.")

  if (is.null(colnames(W))) colnames(W) <- paste0("m", seq_len(ncol(W)))
  if (is.null(rownames(W))) rownames(W) <- paste0("id", seq_len(nrow(W)))
  if (anyDuplicated(colnames(W))) .gsim_stop("W marker names must be unique.")
  if (anyDuplicated(rownames(W))) .gsim_stop("W sample names must be unique.")

  if (!is.null(ids)) {
    if (is.numeric(ids)) {
      if (any(ids < 1 | ids > nrow(W))) .gsim_stop("Numeric ids are out of range.")
      W <- W[ids, , drop = FALSE]
    } else {
      pos <- match(as.character(ids), rownames(W))
      if (anyNA(pos)) .gsim_stop("Some ids were not found in rownames(W).")
      W <- W[pos, , drop = FALSE]
    }
  }
  W
}

.gsim_impute_and_standardize <- function(W, standardize = TRUE) {
  W <- as.matrix(W)
  storage.mode(W) <- "double"

  if (anyNA(W)) {
    means <- colMeans(W, na.rm = TRUE)
    if (any(!is.finite(means))) .gsim_stop("A genotype column contains only missing values.")
    missing <- which(is.na(W), arr.ind = TRUE)
    W[missing] <- means[missing[, 2L]]
  }

  if (!standardize) return(W)

  means <- colMeans(W)
  sds <- apply(W, 2L, stats::sd)
  if (any(!is.finite(sds) | sds <= 0)) {
    bad <- colnames(W)[which(!is.finite(sds) | sds <= 0)]
    .gsim_stop("Causal genotype columns have zero/non-finite variance: ",
               paste(head(bad, 10L), collapse = ", "))
  }
  sweep(sweep(W, 2L, means, "-"), 2L, sds, "/")
}

.gsim_maf_from_W <- function(W) {
  centers <- attr(W, "scaled:center")
  if (!is.null(centers) && length(centers) == ncol(W) &&
      all(is.finite(centers)) && all(centers >= 0 & centers <= 2)) {
    p <- centers / 2
  } else if (all(W[is.finite(W)] >= 0 & W[is.finite(W)] <= 2)) {
    p <- colMeans(W, na.rm = TRUE) / 2
  } else {
    out <- rep(NA_real_, ncol(W))
    names(out) <- colnames(W)
    return(out)
  }
  p <- pmin(pmax(p, 0), 1)
  pmin(p, 1 - p)
}

.gsim_prepare_annotations <- function(A, marker_ids, n_annotations = 0L,
                                      annotation_types = NULL,
                                      annotation_prob = 0.1) {
  m <- length(marker_ids)

  if (is.null(A) && n_annotations > 0L) {
    n_annotations <- as.integer(n_annotations)
    types <- annotation_types %||% rep(c("binary", "continuous"),
                                        length.out = n_annotations)
    types <- rep(types, length.out = n_annotations)
    A <- matrix(0, m, n_annotations)
    for (j in seq_len(n_annotations)) {
      if (types[j] == "binary") {
        prob <- rep(annotation_prob, length.out = n_annotations)[j]
        A[, j] <- stats::rbinom(m, 1L, prob)
      } else if (types[j] == "continuous") {
        A[, j] <- stats::rnorm(m)
      } else {
        .gsim_stop("annotation_types must contain only 'binary' or 'continuous'.")
      }
    }
    colnames(A) <- paste0("A", seq_len(n_annotations))
    rownames(A) <- marker_ids
  }

  if (is.null(A)) {
    return(list(A = NULL, A_centered = NULL, types = character(0)))
  }

  if (is.data.frame(A) && "rsid" %in% names(A)) {
    rn <- as.character(A$rsid)
    A$rsid <- NULL
    A <- as.matrix(A)
    rownames(A) <- rn
  } else {
    A <- as.matrix(A)
  }

  if (!is.null(rownames(A)) &&
      identical(rownames(A), as.character(seq_len(nrow(A)))) &&
      !identical(rownames(A), marker_ids)) {
    rownames(A) <- NULL
  }

  storage.mode(A) <- "double"
  if (is.null(colnames(A))) colnames(A) <- paste0("A", seq_len(ncol(A)))
  if (anyDuplicated(colnames(A))) .gsim_stop("Annotation names must be unique.")

  if (!is.null(rownames(A))) {
    if (anyDuplicated(rownames(A))) .gsim_stop("Annotation marker IDs must be unique.")
    pos <- match(marker_ids, rownames(A))
    if (anyNA(pos)) .gsim_stop("A does not contain every requested marker ID.")
    A <- A[pos, , drop = FALSE]
  } else if (nrow(A) != m) {
    .gsim_stop("A must have one row per marker or row names containing marker IDs.")
  }

  if (any(!is.finite(A))) .gsim_stop("A contains non-finite values.")
  rownames(A) <- marker_ids

  types <- vapply(seq_len(ncol(A)), function(j) {
    x <- unique(A[, j])
    if (all(x %in% c(0, 1))) "binary" else "continuous"
  }, character(1))
  names(types) <- colnames(A)

  centered <- sweep(A, 2L, colMeans(A), "-")
  list(A = A, A_centered = centered, types = types)
}

.gsim_architecture_defaults <- function(architecture, m, n_causal, pi,
                                        mixture_variances) {
  default <- switch(
    architecture,
    bayesc = list(pi = c(0.99, 0.01), v = c(0, 1)),
    bayesr = list(pi = c(0.95, 0.04, 0.009, 0.001),
                  v = c(0, 0.01, 0.1, 1)),
    major_polygenic = list(pi = c(0.20, 0.795, 0.005),
                           v = c(0, 0.01, 1)),
    maf_dependent = list(pi = c(0.95, 0.04, 0.009, 0.001),
                         v = c(0, 0.01, 0.1, 1)),
    clustered = list(pi = c(0.95, 0.04, 0.009, 0.001),
                     v = c(0, 0.01, 0.1, 1)),
    .gsim_stop("No defaults exist for architecture '", architecture, "'.")
  )

  pi_was_null <- is.null(pi)
  pi <- pi %||% default$pi
  mixture_variances <- mixture_variances %||% default$v

  if (!is.numeric(pi) || length(pi) < 2L || any(!is.finite(pi)) || any(pi < 0) ||
      sum(pi) <= 0) {
    .gsim_stop("pi must contain non-negative finite values with positive sum.")
  }
  pi <- pi / sum(pi)

  if (!is.numeric(mixture_variances) || length(mixture_variances) != length(pi) ||
      any(!is.finite(mixture_variances)) || any(mixture_variances < 0)) {
    .gsim_stop("mixture_variances must be non-negative, finite, and match pi.")
  }
  if (mixture_variances[1L] != 0 || any(mixture_variances[-1L] <= 0)) {
    .gsim_stop("The first mixture variance must be zero and all others positive.")
  }

  if (!is.null(n_causal)) {
    .gsim_check_scalar(n_causal, "n_causal", 1, m, integer = TRUE)
    if (pi_was_null) {
      active <- n_causal / m
      rel <- pi[-1L] / sum(pi[-1L])
      pi <- c(1 - active, active * rel)
    }
  }

  list(pi = pi, mixture_variances = mixture_variances)
}

.gsim_prepare_alpha <- function(alpha, A, n_sticks) {
  p <- if (is.null(A)) 0L else ncol(A)
  if (is.null(alpha)) {
    out <- matrix(0, p, n_sticks)
    if (p) rownames(out) <- colnames(A)
    colnames(out) <- paste0("stick", seq_len(n_sticks))
    return(out)
  }

  if (is.null(A)) .gsim_stop("alpha requires an annotation matrix A.")
  if (is.vector(alpha) && n_sticks == 1L) alpha <- matrix(alpha, ncol = 1L)
  alpha <- as.matrix(alpha)
  storage.mode(alpha) <- "double"

  if (!is.null(rownames(alpha))) {
    pos <- match(colnames(A), rownames(alpha))
    if (anyNA(pos)) .gsim_stop("alpha row names do not match colnames(A).")
    alpha <- alpha[pos, , drop = FALSE]
  }
  if (!all(dim(alpha) == c(p, n_sticks))) {
    .gsim_stop("alpha must have dimensions ncol(A) x (length(pi) - 1).")
  }
  if (any(!is.finite(alpha))) .gsim_stop("alpha contains non-finite values.")
  rownames(alpha) <- colnames(A)
  colnames(alpha) <- paste0("stick", seq_len(n_sticks))
  alpha
}

.gsim_probability_surface <- function(pi, mixture_variances, A_centered, alpha) {
  K <- length(pi)
  # Active components are visited from largest to smallest variance; null is last.
  stick_order <- c(order(mixture_variances[-1L], decreasing = TRUE) + 1L, 1L)
  pi_ordered <- pi[stick_order]
  n <- if (is.null(A_centered)) 1L else nrow(A_centered)

  q0 <- vapply(seq_len(K - 1L), function(k) {
    pi_ordered[k] / sum(pi_ordered[k:K])
  }, numeric(1))

  if (is.null(A_centered)) {
    eta <- matrix(rep(stats::qlogis(q0), each = n), n, K - 1L)
  } else {
    eta <- matrix(rep(stats::qlogis(q0), each = n), n, K - 1L) +
      A_centered %*% alpha
  }
  continuation <- stats::plogis(eta)

  ordered <- matrix(0, n, K)
  remaining <- rep(1, n)
  for (k in seq_len(K - 1L)) {
    ordered[, k] <- remaining * continuation[, k]
    remaining <- remaining * (1 - continuation[, k])
  }
  ordered[, K] <- remaining

  probability <- matrix(0, n, K)
  probability[, stick_order] <- ordered
  probability <- probability / rowSums(probability)
  colnames(probability) <- paste0("component", seq_len(K))
  colnames(continuation) <- paste0("stick", seq_len(K - 1L))

  list(probability = probability,
       continuation = continuation,
       stick_order = stick_order,
       baseline_continuation = q0)
}

.gsim_probability_to_continuation <- function(probability, stick_order) {
  ordered <- probability[, stick_order, drop = FALSE]
  K <- ncol(ordered)
  out <- matrix(0, nrow(ordered), K - 1L)
  for (k in seq_len(K - 1L)) {
    out[, k] <- ordered[, k] / pmax(rowSums(ordered[, k:K, drop = FALSE]),
                                     .Machine$double.eps)
  }
  colnames(out) <- paste0("stick", seq_len(K - 1L))
  rownames(out) <- rownames(probability)
  out
}

.gsim_apply_cluster_weights <- function(probability, block_id,
                                        n_hot_blocks = NULL,
                                        cluster_enrichment = 20) {
  if (is.null(block_id)) .gsim_stop("clustered architecture requires block_id.")
  blocks <- unique(block_id[!is.na(block_id)])
  if (!length(blocks)) .gsim_stop("block_id contains no usable blocks.")
  n_hot_blocks <- n_hot_blocks %||% max(1L, round(sqrt(length(blocks))))
  n_hot_blocks <- min(as.integer(n_hot_blocks), length(blocks))
  hot <- sample(blocks, n_hot_blocks)

  p_active <- 1 - probability[, 1L]
  odds <- p_active / pmax(1 - p_active, .Machine$double.eps)
  odds[block_id %in% hot] <- odds[block_id %in% hot] * cluster_enrichment
  p_new <- odds / (1 + odds)

  active_conditional <- probability[, -1L, drop = FALSE] /
    pmax(p_active, .Machine$double.eps)
  probability[, 1L] <- 1 - p_new
  probability[, -1L] <- active_conditional * p_new
  probability <- probability / rowSums(probability)
  list(probability = probability, hot_blocks = hot)
}

.gsim_sample_one <- function(probability) {
  u <- stats::runif(1L)
  which(u <= cumsum(probability))[1L]
}

.gsim_draw_components <- function(probability, n_causal = NULL) {
  m <- nrow(probability)
  K <- ncol(probability)

  if (is.null(n_causal)) {
    u <- stats::runif(m)
    cs <- t(apply(probability, 1L, cumsum))
    component <- 1L + rowSums(u > cs)
    return(pmin(component, K))
  }

  active_weight <- 1 - probability[, 1L]
  if (sum(active_weight) <= 0) .gsim_stop("All markers have zero active probability.")
  causal <- sample(seq_len(m), n_causal, replace = FALSE, prob = active_weight)
  component <- rep(1L, m)
  for (j in causal) {
    p <- probability[j, -1L]
    component[j] <- 1L + .gsim_sample_one(p / sum(p))
  }
  component
}

.gsim_validate_corr <- function(x, nt, name) {
  if (nt == 1L) return(matrix(1, 1L, 1L))
  if (is.null(x)) return(diag(nt))
  if (length(x) == 1L) {
    out <- matrix(x, nt, nt)
    diag(out) <- 1
    x <- out
  }
  x <- as.matrix(x)
  if (!all(dim(x) == c(nt, nt)) || any(!is.finite(x)) ||
      max(abs(x - t(x))) > 1e-10 || any(abs(x) > 1 + 1e-12) ||
      max(abs(diag(x) - 1)) > 1e-10) {
    .gsim_stop(name, " must be a symmetric nt x nt correlation matrix.")
  }
  if (min(eigen(x, symmetric = TRUE, only.values = TRUE)$values) <= 1e-8) {
    .gsim_stop(name, " must be positive definite.")
  }
  x
}

.gsim_draw_beta <- function(component, mixture_variances, nt, rg,
                            maf = NULL, maf_exponent = -0.5) {
  m <- length(component)
  B <- matrix(0, m, nt)
  L <- chol(rg)

  for (k in seq_along(mixture_variances)[-1L]) {
    idx <- which(component == k)
    if (!length(idx)) next
    z <- matrix(stats::rnorm(length(idx) * nt), length(idx), nt) %*% L
    B[idx, ] <- z * sqrt(mixture_variances[k])
  }

  if (!is.null(maf)) {
    active <- which(component != 1L)
    if (any(!is.finite(maf[active]) | maf[active] <= 0 | maf[active] > 0.5)) {
      .gsim_stop("maf must lie in (0, 0.5] for active markers.")
    }
    scale <- (2 * maf[active] * (1 - maf[active]))^(maf_exponent / 2)
    scale <- scale / sqrt(mean(scale^2))
    B[active, ] <- B[active, , drop = FALSE] * scale
  }
  B
}

.gsim_align_fixed_beta <- function(beta, marker_ids, nt) {
  if (is.vector(beta)) beta <- matrix(beta, ncol = 1L,
                                      dimnames = list(names(beta), NULL))
  beta <- as.matrix(beta)
  storage.mode(beta) <- "double"
  if (ncol(beta) != nt) .gsim_stop("beta must have nt columns.")

  if (!is.null(rownames(beta))) {
    if (anyDuplicated(rownames(beta))) .gsim_stop("beta marker names must be unique.")
    B <- matrix(0, length(marker_ids), nt,
                dimnames = list(marker_ids, colnames(beta)))
    pos <- match(rownames(beta), marker_ids)
    if (anyNA(pos)) .gsim_stop("Some beta marker names are absent from the marker catalog.")
    B[pos, ] <- beta
  } else {
    if (nrow(beta) != length(marker_ids)) {
      .gsim_stop("Unnamed beta must have one row per marker.")
    }
    B <- beta
    rownames(B) <- marker_ids
  }
  if (any(!is.finite(B))) .gsim_stop("beta contains non-finite values.")
  B
}

.gsim_resolve_getG <- function(getG_fun) {
  if (!is.null(getG_fun)) return(getG_fun)
  if (!requireNamespace("qgg", quietly = TRUE)) {
    .gsim_stop("qgg is required for Glist simulation, or supply getG_fun.")
  }
  qgg::getG
}

.gsim_call_getG <- function(getG_fun, Glist, marker_ids, sample_ids,
                            chr = NULL) {
  args <- list(Glist = Glist, rsids = marker_ids, ids = sample_ids,
               impute = TRUE, scale = FALSE)
  if (!is.null(chr) && !is.na(chr)) args$chr <- suppressWarnings(type.convert(chr, as.is = TRUE))
  W <- do.call(getG_fun, args)
  if (is.null(dim(W))) W <- matrix(W, ncol = 1L)
  W <- as.matrix(W)

  if (is.null(colnames(W))) {
    if (ncol(W) != length(marker_ids)) .gsim_stop("getG returned unexpected columns.")
    colnames(W) <- marker_ids
  }
  pos <- match(marker_ids, colnames(W))
  if (anyNA(pos)) .gsim_stop("getG did not return every requested causal marker.")
  W <- W[, pos, drop = FALSE]

  if (is.null(rownames(W))) {
    if (nrow(W) != length(sample_ids)) .gsim_stop("getG returned unexpected rows.")
    rownames(W) <- sample_ids
  }
  rpos <- match(sample_ids, rownames(W))
  if (anyNA(rpos)) .gsim_stop("getG did not return every requested sample.")
  W[rpos, , drop = FALSE]
}

.gsim_load_glist_markers <- function(Glist, marker_ids, sample_ids, chr_map,
                                     getG_fun, chunk_size = 5000L) {
  if (!length(marker_ids)) .gsim_stop("No causal markers were selected.")
  chunk_size <- as.integer(chunk_size)
  if (chunk_size < 1L) .gsim_stop("chunk_size must be positive.")

  groups <- split(marker_ids, chr_map[marker_ids], drop = TRUE)
  if (!length(groups)) groups <- list(`NA` = marker_ids)
  chunks <- list()
  z <- 0L

  for (chr_name in names(groups)) {
    ids_chr <- groups[[chr_name]]
    starts <- seq.int(1L, length(ids_chr), by = chunk_size)
    for (start in starts) {
      z <- z + 1L
      take <- ids_chr[start:min(start + chunk_size - 1L, length(ids_chr))]
      chunks[[z]] <- .gsim_call_getG(
        getG_fun, Glist, take, sample_ids,
        chr = if (identical(chr_name, "NA")) NULL else chr_name
      )
    }
  }

  W <- do.call(cbind, chunks)
  W[, match(marker_ids, colnames(W)), drop = FALSE]
}

.gsim_compute_sumstats <- function(W, Y, marker_ids = colnames(W)) {
  W <- as.matrix(W)
  Y <- as.matrix(Y)
  W <- sweep(W, 2L, colMeans(W), "-")
  Y <- sweep(Y, 2L, colMeans(Y), "-")
  ssx <- colSums(W^2)
  if (any(ssx <= 0)) .gsim_stop("Cannot compute statistics for invariant markers.")

  out <- vector("list", ncol(Y))
  for (t in seq_len(ncol(Y))) {
    xy <- as.numeric(crossprod(W, Y[, t]))
    bhat <- xy / ssx
    sse <- pmax(sum(Y[, t]^2) - bhat^2 * ssx, 0)
    sigma2 <- sse / max(nrow(W) - 2L, 1L)
    se <- sqrt(sigma2 / ssx)
    out[[t]] <- data.frame(
      rsid = marker_ids,
      trait = colnames(Y)[t] %||% paste0("D", t),
      beta = bhat,
      se = se,
      z = bhat / se,
      n = nrow(W),
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, out)
}

.gsim_stream_sumstats <- function(Glist, marker_ids, sample_ids, chr_map,
                                  getG_fun, Y, chunk_size, standardize) {
  groups <- split(marker_ids, chr_map[marker_ids], drop = TRUE)
  if (!length(groups)) groups <- list(`NA` = marker_ids)
  out <- list()
  z <- 0L

  for (chr_name in names(groups)) {
    ids_chr <- groups[[chr_name]]
    starts <- seq.int(1L, length(ids_chr), by = chunk_size)
    for (start in starts) {
      take <- ids_chr[start:min(start + chunk_size - 1L, length(ids_chr))]
      W <- .gsim_call_getG(
        getG_fun, Glist, take, sample_ids,
        chr = if (identical(chr_name, "NA")) NULL else chr_name
      )
      W <- .gsim_impute_and_standardize(W, standardize = standardize)
      z <- z + 1L
      out[[z]] <- .gsim_compute_sumstats(W, Y, take)
    }
  }
  do.call(rbind, out)
}

.gsim_exactness <- function(G, E, Y, h2_target, probability = NULL) {
  list(
    max_y_minus_g_plus_e = max(abs(Y - (G + E))),
    max_probability_sum_error = if (is.null(probability)) NA_real_ else
      max(abs(rowSums(probability) - 1)),
    h2_observed = apply(G, 2L, stats::var) / apply(Y, 2L, stats::var),
    h2_target = h2_target
  )
}
