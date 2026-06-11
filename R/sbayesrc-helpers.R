#' Initialize SBayesRC-Style Annotation Mixture Probabilities
#'
#' Creates deterministic initial values for the probit stick-breaking
#' annotation model used by the SBayesRC-style CSR sampler. The requested
#' overall active probability is distributed across the non-null mixture
#' components, then converted to conditional stick-breaking probabilities.
#'
#' When `A` contains an all-ones column, that column receives probit
#' coefficients encoding the baseline stick-breaking probabilities and all
#' other annotation coefficients are initialized to zero. When `A` has no
#' all-ones column, all coefficients are initialized to zero. Boundary
#' stick-breaking probabilities caused by zero component weights are clamped
#' to finite values before applying [stats::qnorm()].
#'
#' @param A Numeric marker annotation matrix with markers in rows and
#'   annotations in columns.
#' @param gamma Numeric mixture variance multipliers. The first value must be
#'   zero and all remaining values must be positive.
#' @param pi_init Initial probability that a marker belongs to any active
#'   component.
#' @param active_comp_weights Optional non-negative weights for distributing
#'   `pi_init` across the active components. The weights are normalized to sum
#'   to one.
#' @param alpha_init Optional annotation coefficient matrix with `ncol(A)` rows
#'   and `length(gamma) - 1` columns.
#' @param sigmaSqAlpha_init Optional positive initial annotation-effect
#'   variances, one per stick-breaking step.
#'
#' @return A list containing the initialized coefficients, variances, component
#'   probabilities, stick-breaking probabilities, mixture multipliers, and
#'   annotation and component names.
#' @export
make_sbayesrc_alpha_init <- function(
  A,
  gamma = c(0, 0.01, 0.1, 1),
  pi_init = 0.001,
  active_comp_weights = NULL,
  alpha_init = NULL,
  sigmaSqAlpha_init = NULL
) {
 if (!is.numeric(gamma) || length(gamma) < 2 ||
     any(!is.finite(gamma))) {
  stop("gamma must be a finite numeric vector with at least two elements.")
 }
 gamma <- as.numeric(gamma)
 if (gamma[1] != 0) stop("gamma[1] must be 0.")
 if (any(gamma[-1] <= 0)) stop("All active gamma values must be positive.")

 if (is.null(dim(A)) || length(dim(A)) != 2) {
  stop("A must be a matrix-like object with two dimensions.")
 }
 A <- as.matrix(A)
 if (!is.numeric(A)) stop("A must be numeric.")
 storage.mode(A) <- "double"
 if (nrow(A) < 1 || ncol(A) < 1) {
  stop("A must have at least one row and one column.")
 }
 if (any(!is.finite(A))) stop("A must contain only finite values.")

 if (!is.numeric(pi_init) || length(pi_init) != 1 ||
     !is.finite(pi_init) || pi_init <= 0 || pi_init >= 1) {
  stop("pi_init must be a finite numeric scalar in (0, 1).")
 }

 K <- ncol(A)
 C <- length(gamma) - 1L

 annotation_names <- colnames(A)
 if (is.null(annotation_names)) {
  annotation_names <- paste0("A", seq_len(K))
 }
 component_names <- c(
  "null",
  paste0("gamma_", format(gamma[-1], trim = TRUE, scientific = FALSE))
 )
 step_names <- paste0("step_", seq_len(C))

 if (is.null(active_comp_weights)) {
  active_comp_weights <- rep(1 / C, C)
 } else {
  if (!is.numeric(active_comp_weights) ||
      length(active_comp_weights) != C ||
      any(!is.finite(active_comp_weights)) ||
      any(active_comp_weights < 0) ||
      sum(active_comp_weights) <= 0) {
   stop(
    "active_comp_weights must be a finite, non-negative numeric vector ",
    "of length length(gamma) - 1 with positive sum."
   )
  }
  active_comp_weights <- as.numeric(active_comp_weights)
  active_comp_weights <- active_comp_weights / sum(active_comp_weights)
 }
 names(active_comp_weights) <- component_names[-1]

 component_prob_init <- c(
  1 - pi_init,
  pi_init * active_comp_weights
 )
 names(component_prob_init) <- component_names

 remaining_prob <- rev(cumsum(rev(component_prob_init)))
 step_prob_init <- numeric(C)
 for (j in seq_len(C)) {
  if (remaining_prob[j] > 0) {
   step_prob_init[j] <- remaining_prob[j + 1L] / remaining_prob[j]
  } else {
   step_prob_init[j] <- 0.5
  }
 }

 prob_epsilon <- .Machine$double.eps
 step_prob_init <- pmin(pmax(step_prob_init, prob_epsilon), 1 - prob_epsilon)
 names(step_prob_init) <- step_names

 if (is.null(alpha_init)) {
  alpha_init <- matrix(
   0,
   nrow = K,
   ncol = C,
   dimnames = list(annotation_names, step_names)
  )
  intercept_columns <- which(
   vapply(
    seq_len(K),
    function(j) all(abs(A[, j] - 1) < 1e-12),
    logical(1)
   )
  )
  if (length(intercept_columns) > 0) {
   alpha_init[intercept_columns[1], ] <- stats::qnorm(step_prob_init)
  }
 } else {
  if (!is.numeric(alpha_init) || is.null(dim(alpha_init)) ||
      length(dim(alpha_init)) != 2 ||
      !all(dim(alpha_init) == c(K, C)) ||
      any(!is.finite(alpha_init))) {
   stop(
    "alpha_init must be a finite numeric matrix with dimensions ",
    "ncol(A) x (length(gamma) - 1)."
   )
  }
  alpha_init <- as.matrix(alpha_init)
  storage.mode(alpha_init) <- "double"
  if (is.null(rownames(alpha_init))) rownames(alpha_init) <- annotation_names
  if (is.null(colnames(alpha_init))) colnames(alpha_init) <- step_names
 }

 if (is.null(sigmaSqAlpha_init)) {
  sigmaSqAlpha_init <- rep(1, C)
 } else {
  if (!is.numeric(sigmaSqAlpha_init) ||
      length(sigmaSqAlpha_init) != C ||
      any(!is.finite(sigmaSqAlpha_init)) ||
      any(sigmaSqAlpha_init <= 0)) {
   stop(
    "sigmaSqAlpha_init must be a positive finite numeric vector of length ",
    "length(gamma) - 1."
   )
  }
  sigmaSqAlpha_init <- as.numeric(sigmaSqAlpha_init)
 }
 names(sigmaSqAlpha_init) <- step_names

 list(
  alpha_init = alpha_init,
  sigmaSqAlpha_init = sigmaSqAlpha_init,
  active_comp_weights = active_comp_weights,
  component_prob_init = component_prob_init,
  step_prob_init = step_prob_init,
  gamma = gamma,
  annotation_names = annotation_names,
  component_names = component_names
 )
}

.sbayesrc_validate_gamma <- function(gamma) {
 if (!is.numeric(gamma) || length(gamma) < 2 ||
     any(!is.finite(gamma))) {
  stop("gamma must be a finite numeric vector with at least two elements.")
 }
 gamma <- as.numeric(gamma)
 if (gamma[1] != 0) stop("gamma[1] must be 0.")
 if (any(gamma[-1] <= 0)) stop("All active gamma values must be positive.")
 gamma
}

.sbayesrc_validate_alpha <- function(alpha, gamma) {
 if (!is.matrix(alpha) || !is.numeric(alpha) ||
     any(!is.finite(alpha))) {
  stop("alpha must be a finite numeric matrix.")
 }
 if (ncol(alpha) != length(gamma) - 1L) {
  stop("ncol(alpha) must equal length(gamma) - 1.")
 }
 alpha
}

.sbayesrc_stick_breaking_pi <- function(eta, gamma) {
 p <- stats::pnorm(eta)
 n <- nrow(p)
 ncomp <- length(gamma)
 component_names <- paste0(
  "gamma_", format(gamma, trim = TRUE, scientific = FALSE)
 )
 out <- matrix(
  0,
  nrow = n,
  ncol = ncomp,
  dimnames = list(rownames(eta), component_names)
 )

 remaining <- rep(1, n)
 for (j in seq_len(ncomp - 1L)) {
  out[, j] <- remaining * (1 - p[, j])
  remaining <- remaining * p[, j]
 }
 out[, ncomp] <- remaining
 out
}

#' Convert SBayesRC-Style Annotation Coefficients to Component Probabilities
#'
#' Applies the generalized probit stick-breaking transform to annotation
#' coefficient rows.
#'
#' @param alpha Numeric annotation coefficient matrix with one row per
#'   annotation and `length(gamma) - 1` columns.
#' @param gamma Numeric mixture variance multipliers. The first value must be
#'   zero and all remaining values must be positive.
#'
#' @return An annotation by mixture-component probability matrix.
#' @export
sbayesrc_annotation_pi <- function(alpha, gamma = c(0, 0.01, 0.1, 1)) {
 gamma <- .sbayesrc_validate_gamma(gamma)
 alpha <- .sbayesrc_validate_alpha(alpha, gamma)
 .sbayesrc_stick_breaking_pi(alpha, gamma)
}

#' Calculate Expected SBayesRC-Style Gamma by Annotation
#'
#' @inheritParams sbayesrc_annotation_pi
#'
#' @return A named numeric vector containing expected gamma per annotation.
#' @export
sbayesrc_annotation_gamma_mean <- function(
  alpha,
  gamma = c(0, 0.01, 0.1, 1)
) {
 gamma <- .sbayesrc_validate_gamma(gamma)
 pi <- sbayesrc_annotation_pi(alpha, gamma)
 out <- as.numeric(pi %*% gamma)
 names(out) <- rownames(pi)
 out
}

#' Convert SBayesRC-Style Marker Annotations to Component Probabilities
#'
#' Combines marker annotations with SBayesRC-style annotation coefficients,
#' then applies the generalized probit stick-breaking transform.
#'
#' @param A Numeric marker annotation matrix with markers in rows and
#'   annotations in columns.
#' @inheritParams sbayesrc_annotation_pi
#'
#' @return A marker by mixture-component probability matrix.
#' @export
sbayesrc_marker_pi <- function(A, alpha, gamma = c(0, 0.01, 0.1, 1)) {
 gamma <- .sbayesrc_validate_gamma(gamma)
 alpha <- .sbayesrc_validate_alpha(alpha, gamma)
 if (!is.matrix(A) || !is.numeric(A) || any(!is.finite(A))) {
  stop("A must be a finite numeric matrix.")
 }
 if (ncol(A) != nrow(alpha)) {
  stop("ncol(A) must equal nrow(alpha).")
 }

 eta <- A %*% alpha
 rownames(eta) <- rownames(A)
 .sbayesrc_stick_breaking_pi(eta, gamma)
}

#' Calculate Expected SBayesRC-Style Gamma by Marker
#'
#' @inheritParams sbayesrc_marker_pi
#'
#' @return A named numeric vector containing expected gamma per marker.
#' @export
sbayesrc_marker_gamma_mean <- function(
  A,
  alpha,
  gamma = c(0, 0.01, 0.1, 1)
) {
 gamma <- .sbayesrc_validate_gamma(gamma)
 pi <- sbayesrc_marker_pi(A, alpha, gamma)
 out <- as.numeric(pi %*% gamma)
 names(out) <- rownames(pi)
 out
}

format_sbayesrc_csr_fit <- function(
  fit,
  nt,
  m,
  gamma,
  n_anno,
  trait_names = NULL,
  variable_names = NULL,
  annotation_names = NULL
) {
 validate_positive_integer <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1 || !is.finite(x) ||
      x <= 0 || x != as.integer(x)) {
   stop(name, " must be a positive integer.")
  }
  as.integer(x)
 }

 validate_trait_list <- function(x, lengths_expected, name) {
  if (!is.list(x) || length(x) != nt ||
      any(lengths(x) != lengths_expected)) {
   stop(
    name, " must be a list of length nt with element length ",
    lengths_expected, "."
   )
  }
 }

 nt <- validate_positive_integer(nt, "nt")
 m <- validate_positive_integer(m, "m")
 n_anno <- validate_positive_integer(n_anno, "n_anno")

 if (!is.list(fit) || length(fit) != 24) {
  stop("format_sbayesrc_csr_fit() expects the 24-slot SBayesRC CSR return object.")
 }

 if (!is.numeric(gamma) || length(gamma) < 2 ||
     any(!is.finite(gamma))) {
  stop("gamma must be a finite numeric vector with at least two elements.")
 }
 gamma <- as.numeric(gamma)
 if (gamma[1] != 0) stop("gamma[1] must be 0.")
 if (any(gamma[-1] <= 0)) stop("All active gamma values must be positive.")

 if (is.null(trait_names)) trait_names <- paste0("T", seq_len(nt))
 if (is.null(variable_names)) variable_names <- paste0("V", seq_len(m))
 if (is.null(annotation_names)) annotation_names <- paste0("A", seq_len(n_anno))

 if (length(trait_names) != nt) stop("trait_names must have length nt.")
 if (length(variable_names) != m) stop("variable_names must have length m.")
 if (length(annotation_names) != n_anno) {
  stop("annotation_names must have length n_anno.")
 }

 Kgamma <- length(gamma)
 nstep <- Kgamma - 1L
 component_names <- paste0(
  "gamma_", format(gamma, trim = TRUE, scientific = FALSE)
 )
 step_names <- paste0("step_", seq_len(nstep))
 alpha_names <- as.vector(outer(annotation_names, step_names, paste, sep = ":"))

 for (i in seq_len(7)) {
  validate_trait_list(fit[[i]], m, paste0("fit[[", i, "]]"))
 }
 for (i in 8:10) {
  if (!is.list(fit[[i]]) || length(fit[[i]]) != nt ||
      length(unique(lengths(fit[[i]]))) != 1) {
   stop("fit[[", i, "]] must contain equally sized trace vectors for each trait.")
  }
 }
 for (i in 11:16) {
  validate_trait_list(fit[[i]], nt, paste0("fit[[", i, "]]"))
 }
 for (i in 17:18) {
  validate_trait_list(fit[[i]], 2, paste0("fit[[", i, "]]"))
 }
 validate_trait_list(fit[[19]], n_anno * nstep, "fit[[19]]")
 validate_trait_list(fit[[20]], nstep, "fit[[20]]")
 for (i in 21:22) {
  if (!is.list(fit[[i]]) || length(fit[[i]]) != nt ||
      length(unique(lengths(fit[[i]]))) != 1) {
   stop("fit[[", i, "]] must contain equally sized trace vectors for each trait.")
  }
 }
 validate_trait_list(fit[[23]], m * Kgamma, "fit[[23]]")
 validate_trait_list(fit[[24]], Kgamma, "fit[[24]]")

 names(fit) <- c(
  "bm", "dm", "wy", "r", "b", "component", "marker_index",
  "vbs", "vgs", "ves",
  "covb", "covg", "cove",
  "vb", "vg", "ve",
  "pi", "pim", "alpha_raw", "sigmaSqAlpha_raw",
  "vle", "vld", "comp_prob_raw", "ncomp_raw"
 )

 for (i in seq_len(7)) {
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
  rownames(fit[[i]]) <- colnames(fit[[i]]) <- trait_names
 }

 for (i in 17:18) {
  fit[[i]] <- matrix(unlist(fit[[i]]), ncol = 2, byrow = TRUE)
  rownames(fit[[i]]) <- trait_names
  colnames(fit[[i]]) <- c("pi0", "pi_active")
 }

 alpha_flat <- matrix(
  unlist(fit$alpha_raw),
  nrow = nt,
  ncol = n_anno * nstep,
  byrow = TRUE,
  dimnames = list(trait_names, alpha_names)
 )
 alpha <- lapply(seq_len(nt), function(t) {
  matrix(
   alpha_flat[t, ],
   nrow = n_anno,
   ncol = nstep,
   dimnames = list(annotation_names, step_names)
  )
 })
 names(alpha) <- trait_names

 sigmaSqAlpha <- matrix(
  unlist(fit$sigmaSqAlpha_raw),
  nrow = nt,
  ncol = nstep,
  byrow = TRUE,
  dimnames = list(trait_names, step_names)
 )

 for (i in 21:22) {
  fit[[i]] <- as.matrix(as.data.frame(fit[[i]]))
  rownames(fit[[i]]) <- paste0("Iter", seq_len(nrow(fit[[i]])))
  colnames(fit[[i]]) <- trait_names
 }

 comp_prob <- lapply(seq_len(nt), function(t) {
  matrix(
   fit$comp_prob_raw[[t]],
   nrow = m,
   ncol = Kgamma,
   dimnames = list(variable_names, component_names)
  )
 })
 names(comp_prob) <- trait_names

 ncomp <- matrix(
  unlist(fit$ncomp_raw),
  nrow = nt,
  ncol = Kgamma,
  byrow = TRUE,
  dimnames = list(trait_names, component_names)
 )

 out <- fit[1:18]
 out$alpha_flat <- alpha_flat
 out$alpha <- alpha
 out$sigmaSqAlpha <- sigmaSqAlpha
 out$vle <- fit$vle
 out$vld <- fit$vld
 out$comp_prob <- comp_prob
 out$ncomp <- ncomp

 valid_covariance <- function(x) {
  is.matrix(x) && all(dim(x) == c(nt, nt)) &&
   all(is.finite(x)) && all(diag(x) > 0)
 }
 if (valid_covariance(out$covb)) out$rb <- stats::cov2cor(out$covb)
 if (valid_covariance(out$covg)) out$rg <- stats::cov2cor(out$covg)
 if (valid_covariance(out$cove)) out$re <- stats::cov2cor(out$cove)

 out
}
