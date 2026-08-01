.stblr_calibration_matrix <- function(x, nt, m, name, positive = FALSE) {
 if (is.null(x)) x <- 1
 if (is.list(x)) {
  if (length(x) != nt || any(lengths(x) != m)) {
   stop(name, " must have one length-m vector per trait.", call. = FALSE)
  }
  x <- do.call(rbind, lapply(x, as.numeric))
 } else if (is.matrix(x)) {
  if (identical(dim(x), c(m, nt))) x <- t(x)
  if (!identical(dim(x), c(nt, m))) {
   stop(name, " must be nt by m (or m by nt).", call. = FALSE)
  }
 } else {
  x <- as.numeric(x)
  if (length(x) == 1L) x <- matrix(x, nt, m)
  else if (length(x) == m) x <- matrix(rep(x, each = nt), nt, m)
  else if (length(x) == nt) x <- matrix(rep(x, m), nt, m)
  else stop(name, " has incompatible length.", call. = FALSE)
 }
 storage.mode(x) <- "double"
 if (any(!is.finite(x)) || any(x < 0) || (positive && any(x <= 0))) {
  stop(name, " must contain ", if (positive) "positive " else "non-negative ",
       "finite values.", call. = FALSE)
 }
 x
}

.stblr_scalar_prior_calibration <- function(
  vy, h2, nub, nue, expected_multiplier_initial,
  expected_multiplier_prior = expected_multiplier_initial,
  marker_scale = 1, variance_multiplier = 1, trait_names = NULL,
  B = NULL, E = NULL, ssb_prior = NULL, sse_prior = NULL,
  component_probability_source = "global",
  annotation_probability_policy = "not_applicable"
) {
 vy <- as.numeric(vy)
 nt <- length(vy)
 h2 <- rep(as.numeric(h2), length.out = nt)
 if (any(!is.finite(vy)) || any(vy <= 0) || any(!is.finite(h2)) ||
     any(h2 <= 0 | h2 >= 1)) {
  stop("vy must be positive and h2 must be in (0, 1).", call. = FALSE)
 }
 if (is.null(trait_names)) trait_names <- names(vy)
 if (is.null(trait_names)) trait_names <- paste0("T", seq_len(nt))
 candidates <- list(expected_multiplier_initial, expected_multiplier_prior,
                    marker_scale, variance_multiplier)
 lengths_seen <- unlist(lapply(candidates, function(x) {
  if (is.list(x)) lengths(x) else if (is.matrix(x)) dim(x) else length(x)
 }))
 m <- max(lengths_seen[lengths_seen > nt], 1L)
 if (any(vapply(candidates, is.list, logical(1)))) {
  m <- max(vapply(candidates[vapply(candidates, is.list, logical(1))],
                  function(x) max(lengths(x)), numeric(1)))
 } else if (any(vapply(candidates, is.matrix, logical(1)))) {
  dims <- lapply(candidates[vapply(candidates, is.matrix, logical(1))], dim)
  m <- max(vapply(dims, max, numeric(1)))
 }
 initial <- .stblr_calibration_matrix(
  expected_multiplier_initial, nt, m, "expected_multiplier_initial")
 prior <- .stblr_calibration_matrix(
  expected_multiplier_prior, nt, m, "expected_multiplier_prior")
 q <- .stblr_calibration_matrix(marker_scale, nt, m, "marker_scale", TRUE)
 v <- .stblr_calibration_matrix(
  variance_multiplier, nt, m, "variance_multiplier", TRUE)
 weight_initial <- rowSums(initial * q * v)
 weight_prior <- rowSums(prior * q * v)
 if (any(!is.finite(weight_initial)) || any(weight_initial <= 0) ||
     any(!is.finite(weight_prior)) || any(weight_prior <= 0)) {
  stop("Resolved scalar prior weights must be positive and finite.",
       call. = FALSE)
 }
 target <- vy * h2
 B_default <- diag(target / weight_initial, nt)
 E_default <- diag(vy * (1 - h2), nt)
 ssb_default <- diag(((nub - 2) / nub) * target / weight_prior, nt)
 sse_default <- diag(((nue - 2) / nue) * vy * (1 - h2), nt)
 validate <- function(x, default, name) {
  if (is.null(x)) x <- default
  x <- as.matrix(x)
  if (!identical(dim(x), c(nt, nt)) || any(!is.finite(x))) {
   stop(name, " must be a finite nt by nt matrix.", call. = FALSE)
  }
  rownames(x) <- colnames(x) <- trait_names
  x
 }
 B <- validate(B, B_default, "B")
 E <- validate(E, E_default, "E")
 ssb_prior <- validate(ssb_prior, ssb_default, "ssb_prior")
 sse_prior <- validate(sse_prior, sse_default, "sse_prior")
 marker_scale_sum <- rowSums(q)
 variance_multiplier_weight <- rowSums(initial * q * v) /
  pmax(rowSums(initial * q), .Machine$double.eps)
 metadata <- list(
  prior_calibration_policy = "resolved_expected_genetic_variance",
  prior_calibration_version = 1L,
  prior_weight_initial = stats::setNames(weight_initial, trait_names),
  prior_weight_prior_mean = stats::setNames(weight_prior, trait_names),
  mixture_weight_initial = stats::setNames(rowSums(initial), trait_names),
  mixture_weight_prior_mean = stats::setNames(rowSums(prior), trait_names),
  marker_scale_sum = stats::setNames(marker_scale_sum, trait_names),
  variance_multiplier_weight = stats::setNames(
   variance_multiplier_weight, trait_names),
  component_probability_source = component_probability_source,
  annotation_probability_policy = annotation_probability_policy,
  h2_prior_interpretation =
   "requested initial expected genetic variance under resolved prior weights"
 )
 list(vy = vy, B = B, E = E, ssb_prior = ssb_prior,
      sse_prior = sse_prior,
      ssb_prior_list = split(ssb_prior, rep(seq_len(nt), each = nt)),
      sse_prior_list = split(sse_prior, rep(seq_len(nt), each = nt)),
      calibration = metadata)
}

.stblr_calibration_maf_scale <- function(info, estimate, init, m) {
 if (isTRUE(info$fixed)) return(as.numeric(info$prior_scale))
 if (isTRUE(estimate)) {
  if (length(info$log_h) != m) {
   stop("Sampled MAF-S calibration requires aligned marker heterozygosity.",
        call. = FALSE)
  }
  return(exp((as.numeric(init) + 1) * as.numeric(info$log_h)))
 }
 rep(1, m)
}

.mtblr_prior_weight_matrix <- function(patterns, probability, gamma,
                                       marker_scale, name) {
 patterns <- as.matrix(patterns)
 nt <- ncol(patterns)
 ns <- nrow(patterns)
 m <- length(marker_scale)
 if (is.vector(probability)) probability <- matrix(
  rep(as.numeric(probability), each = m), m, ns)
 probability <- as.matrix(probability)
 if (!identical(dim(probability), c(m, ns)) ||
     any(!is.finite(probability)) || any(probability < 0) ||
     any(abs(rowSums(probability) - 1) > 1e-8)) {
  stop(name, " probabilities must be an m by state row-stochastic matrix.",
       call. = FALSE)
 }
 gamma <- as.numeric(gamma)
 if (length(gamma) != ns || any(!is.finite(gamma)) || any(gamma < 0)) {
  stop(name, " component multipliers must match states.", call. = FALSE)
 }
 marker_scale <- as.numeric(marker_scale)
 if (length(marker_scale) != m || any(!is.finite(marker_scale)) ||
     any(marker_scale <= 0)) stop("Invalid MT marker scale.", call. = FALSE)
 W <- matrix(0, nt, nt)
 for (s in seq_len(ns)) {
  state_weight <- sum(probability[, s] * gamma[s] * marker_scale)
  W <- W + state_weight * tcrossprod(patterns[s, ])
 }
 W
}

.mtblr_prior_calibration <- function(
  vy, h2, nub, nue, patterns, probability_initial,
  probability_prior = probability_initial, gamma = NULL,
  marker_scale = NULL, vg = NULL, vb = NULL, ve = NULL,
  ssb_prior = NULL, sse_prior = NULL,
  component_probability_source = "joint_global",
  annotation_probability_policy = "not_applicable"
) {
 vy <- as.numeric(vy); nt <- length(vy)
 h2 <- rep(as.numeric(h2), length.out = nt)
 patterns <- as.matrix(patterns)
 if (ncol(patterns) != nt) stop("MT patterns must have nt columns.", call. = FALSE)
 ns <- nrow(patterns)
 if (is.null(gamma)) gamma <- rep(1, ns)
 if (is.null(marker_scale)) {
  m <- if (is.matrix(probability_initial)) nrow(probability_initial) else 1L
  marker_scale <- rep(1, m)
 }
 m <- length(marker_scale)
 expand_probability <- function(p) {
  if (is.vector(p)) matrix(rep(as.numeric(p), each = m), m, ns) else as.matrix(p)
 }
 p0 <- expand_probability(probability_initial)
 pp <- expand_probability(probability_prior)
 W0 <- .mtblr_prior_weight_matrix(patterns, p0, gamma, marker_scale, "Initial")
 Wp <- .mtblr_prior_weight_matrix(patterns, pp, gamma, marker_scale, "Prior-mean")
 if (any(diag(W0) <= 0) || any(diag(Wp) <= 0)) {
  stop("Every MT trait must have positive resolved prior weight.", call. = FALSE)
 }
 target <- vy * h2
 vg_default <- diag(target, nt)
 ve_default <- diag(vy * (1 - h2), nt)
 if (is.null(vb) && !is.null(vg)) {
  vg_matrix <- as.matrix(vg)
  if (any(vg_matrix[row(vg_matrix) != col(vg_matrix)] != 0)) {
   stop("Automatic calibration of a full vg is not guaranteed positive semidefinite; supply an explicit positive-definite vb.",
        call. = FALSE)
  }
 }
 vb_default <- diag(target / diag(W0), nt)
 ssb_default <- diag(((nub - 2) / nub) * target / diag(Wp), nt)
 vg <- .mtblr_cov(vg, vg_default, "vg", nt)
 vb <- .mtblr_cov(vb, vb_default, "vb", nt)
 ve <- .mtblr_cov(ve, ve_default, "ve", nt)
 ssb_prior <- .mtblr_cov(ssb_prior, ssb_default, "ssb_prior", nt)
 sse_prior <- .mtblr_cov(
  sse_prior, ((nue - 2) / nue) * ve, "sse_prior", nt)
 metadata <- list(
  prior_calibration_policy = "joint_state_expected_covariance_diagonal_default",
  prior_calibration_version = 1L,
  prior_weight_initial = W0,
  prior_weight_prior_mean = Wp,
  marker_scale_sum = sum(marker_scale),
  component_probability_source = component_probability_source,
  annotation_probability_policy = annotation_probability_policy,
  h2_prior_interpretation =
   "requested initial expected trait genetic variances under resolved joint-state weights"
 )
 list(vg = vg, vb = vb, ve = ve, ssb_prior = ssb_prior,
      sse_prior = sse_prior, calibration = metadata)
}
