#' Initialize SBayesRC Annotation Mixture Probabilities
#'
#' Creates deterministic initial values for the probit stick-breaking
#' annotation model used by the SBayesRC CSR sampler. The requested overall
#' active probability is distributed across the non-null mixture components,
#' then converted to conditional stick-breaking probabilities.
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
