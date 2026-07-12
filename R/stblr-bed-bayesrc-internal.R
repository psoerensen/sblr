.bayesr_pi_to_probit_stick_intercepts <- function(pi, pi_floor = 1e-12) {
 if (!is.numeric(pi) || length(pi) < 2L || anyNA(pi) ||
     any(!is.finite(pi)) || any(pi <= 0)) {
  stop("pi must contain at least two positive finite component probabilities.")
 }
 if (!is.numeric(pi_floor) || length(pi_floor) != 1L ||
     !is.finite(pi_floor) || pi_floor <= 0 || pi_floor >= 0.5) {
  stop("pi_floor must be a finite scalar in (0, 0.5).")
 }
 pi <- as.numeric(pi / sum(pi))
 remaining <- 1
 stick <- numeric(length(pi) - 1L)
 for (k in seq_along(stick)) {
  stick[k] <- 1 - pi[k] / remaining
  remaining <- remaining * stick[k]
 }
 stick <- pmin(pmax(stick, pi_floor), 1 - pi_floor)
 matrix(stats::qnorm(stick), nrow = 1L,
        dimnames = list("intercept", paste0("step_", seq_along(stick))))
}

.stblr_bed_bayesrc_native <- function(...) {
 args <- list(...)
 if (is.null(args$A)) stop("A is required and must be aligned to selected BED markers.")
 raw <- do.call(stblr_cpg_omp_bed_marker_scheduled_chains_bayesrc, args)
 raw$annotation$annotation_names <- colnames(args$A) %||%
  paste0("A", seq_len(ncol(args$A)))
 raw
}
