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

.stblr_align_bed_bayesrc_annotations <- function(
  annotation,
  selected_marker_ids,
  add_intercept = TRUE,
  standardize_annotations = TRUE,
  center_binary_annotations = FALSE
) {
 if (is.null(annotation)) {
  stop("annotation is required for method = 'bayesrc'.")
 }
 selected_marker_ids <- as.character(selected_marker_ids)
 if (anyNA(selected_marker_ids) || any(!nzchar(selected_marker_ids))) {
  stop("Selected BED marker IDs must be non-missing and non-empty for BayesRC alignment.")
 }
 if (anyDuplicated(selected_marker_ids)) {
  dup <- unique(selected_marker_ids[duplicated(selected_marker_ids)])
  stop(
   "Selected BED marker IDs must be unique. First duplicates: ",
   paste(utils::head(dup, 10L), collapse = ", ")
  )
 }

 x <- annotation
 marker_ids <- NULL
 id_source <- "row_names"
 if (is.data.frame(x)) {
  id_names <- names(x)[tolower(names(x)) %in%
   c("marker_id", "marker", "rsid", "rsids")]
  if (length(id_names) > 1L) {
   stop("annotation contains multiple marker-ID columns: ", paste(id_names, collapse = ", "))
  }
  if (length(id_names) == 1L) {
   marker_ids <- as.character(x[[id_names]])
   x[[id_names]] <- NULL
   id_source <- paste0("column:", id_names)
  }
 }
 if (is.null(marker_ids)) {
  rn <- rownames(x)
  default_rn <- !is.null(rn) && identical(rn, as.character(seq_len(nrow(x))))
  if (!is.null(rn) && !default_rn) marker_ids <- as.character(rn)
 }

 aligned_by_id <- !is.null(marker_ids)
 unused_rows <- 0L
 if (aligned_by_id) {
  if (anyNA(marker_ids) || any(!nzchar(marker_ids))) {
   stop("Annotation marker IDs must be non-missing and non-empty.")
  }
  if (anyDuplicated(marker_ids)) {
   dup <- unique(marker_ids[duplicated(marker_ids)])
   stop(
    "Annotation marker IDs must be unique. First duplicates: ",
    paste(utils::head(dup, 10L), collapse = ", ")
   )
  }
  idx <- match(selected_marker_ids, marker_ids)
  if (anyNA(idx)) {
   missing_ids <- selected_marker_ids[is.na(idx)]
   stop(
    length(missing_ids), " selected BED marker(s) are missing annotation rows. ",
    "First missing marker IDs: ",
    paste(utils::head(missing_ids, 10L), collapse = ", ")
   )
  }
  unused_rows <- length(marker_ids) - length(unique(idx))
  x <- x[idx, , drop = FALSE]
 } else if (nrow(x) != length(selected_marker_ids)) {
  stop(
   "Annotation rows have no marker IDs; nrow(annotation) must equal the ",
   "final selected BED marker count (", length(selected_marker_ids), ")."
  )
 }

 if (is.data.frame(x)) {
  if (ncol(x) < 1L) stop("annotation must contain at least one annotation column.")
  pieces <- lapply(names(x), function(nm) {
   value <- x[[nm]]
   if (is.factor(value)) {
    mm <- tryCatch(
     stats::model.matrix(~ value - 1, na.action = stats::na.pass),
     error = function(e) stop(
      "Could not expand annotation factor '", nm, "': ", conditionMessage(e)
     )
    )
    colnames(mm) <- paste0(nm, levels(value))
    return(mm)
   }
   if (!is.numeric(value) && !is.integer(value) && !is.logical(value)) {
    stop("Annotation data-frame column '", nm, "' must be numeric, logical, or a factor.")
   }
   matrix(as.numeric(value), ncol = 1L, dimnames = list(NULL, nm))
  })
  x <- do.call(cbind, pieces)
 } else {
  if (is.null(dim(x)) || length(dim(x)) != 2L) {
   stop("annotation must be a numeric matrix or data frame.")
  }
  if (!is.numeric(x) && !is.logical(x) && !is.integer(x)) {
   stop("A matrix annotation must be numeric, integer, or logical.")
  }
  x <- as.matrix(x)
  storage.mode(x) <- "double"
 }
 if (ncol(x) < 1L) stop("annotation must contain at least one annotation column.")
 if (is.null(colnames(x))) colnames(x) <- paste0("Anno", seq_len(ncol(x)))
 if (anyNA(colnames(x)) || any(!nzchar(colnames(x)))) {
  stop("Annotation column names must be non-missing and non-empty.")
 }
 if (anyDuplicated(colnames(x))) {
  stop("Annotation column names must be unique.")
 }
 if (any(!is.finite(x))) {
  stop("annotation contains non-finite values after preprocessing.")
 }

 intercept_cols <- which(vapply(
  seq_len(ncol(x)), function(j) all(abs(x[, j] - 1) < 1e-12), logical(1)
 ))
 if (length(intercept_cols) > 1L) {
  stop("annotation contains more than one all-ones intercept column.")
 }
 if (length(intercept_cols) == 1L && intercept_cols != 1L) {
  x <- x[, c(intercept_cols, setdiff(seq_len(ncol(x)), intercept_cols)), drop = FALSE]
  intercept_cols <- 1L
 }
 non_intercept <- setdiff(seq_len(ncol(x)), intercept_cols)
 zero_variance <- non_intercept[vapply(
  non_intercept, function(j) {
   s <- stats::sd(x[, j])
   !is.finite(s) || s == 0
  }, logical(1)
 )]
 if (length(zero_variance)) {
  stop(
   "Non-intercept annotation columns must have positive variance. Invalid columns: ",
   paste(colnames(x)[zero_variance], collapse = ", ")
  )
 }

 intercept_added <- isTRUE(add_intercept) && length(intercept_cols) == 0L
 x <- .stblr_prepare_annotation_matrix(
  A = x,
  m = length(selected_marker_ids),
  variable_names = selected_marker_ids,
  add_intercept = add_intercept,
  standardize = standardize_annotations,
  center_binary = center_binary_annotations
 )
 if (nrow(x) != length(selected_marker_ids) ||
     !identical(rownames(x), selected_marker_ids)) {
  stop("Internal BayesRC annotation alignment failed.")
 }

 list(
  A = x,
  alignment = list(
   matched_by_id = aligned_by_id,
   marker_id_source = if (aligned_by_id) id_source else "prealigned_row_order",
   selected_marker_count = length(selected_marker_ids),
   unused_annotation_rows = as.integer(unused_rows)
  ),
  preprocessing = list(
   add_intercept = isTRUE(add_intercept),
   intercept_added = intercept_added,
   intercept_name = if (ncol(x) && all(abs(x[, 1L] - 1) < 1e-12)) colnames(x)[1L] else NULL,
   standardize_annotations = isTRUE(standardize_annotations),
   center_binary_annotations = isTRUE(center_binary_annotations),
   factor_expansion = is.data.frame(annotation)
  )
 )
}

.stblr_initialize_bed_bayesrc_prior <- function(
  A,
  mixture_var,
  pi,
  annot_alpha_init = NULL,
  annot_sigma_sq_alpha_init = NULL,
  pi_floor = 1e-12
) {
 if (!is.numeric(mixture_var) || length(mixture_var) < 2L ||
     any(!is.finite(mixture_var)) || mixture_var[1L] != 0 ||
     any(mixture_var[-1L] <= 0)) {
  stop("mixture_var must start with 0 and have positive finite non-null components.")
 }
 if (!is.numeric(pi) || length(pi) != length(mixture_var) ||
     any(!is.finite(pi)) || any(pi <= 0) || sum(pi) <= 0) {
  stop("For BayesRC, pi must be a positive finite vector matching mixture_var.")
 }
 pi <- as.numeric(pi / sum(pi))
 P <- ncol(A)
 K_minus_one <- length(mixture_var) - 1L
 step_names <- paste0("step_", seq_len(K_minus_one))
 if (is.null(annot_alpha_init)) {
  intercept <- which(vapply(
   seq_len(P), function(j) all(abs(A[, j] - 1) < 1e-12), logical(1)
  ))
  if (!length(intercept)) {
   stop("An intercept column or explicit annot_alpha_init is required for BayesRC initialization.")
  }
  annot_alpha_init <- matrix(
   0, P, K_minus_one,
   dimnames = list(colnames(A), step_names)
  )
  annot_alpha_init[intercept[1L], ] <-
   drop(.bayesr_pi_to_probit_stick_intercepts(pi, pi_floor = pi_floor))
 } else {
  if (!is.numeric(annot_alpha_init) || is.null(dim(annot_alpha_init)) ||
      length(dim(annot_alpha_init)) != 2L ||
      !identical(dim(annot_alpha_init), c(P, K_minus_one)) ||
      any(!is.finite(annot_alpha_init))) {
   stop(
    "annot_alpha_init must be a finite numeric matrix with dimensions ",
    P, " x ", K_minus_one, "."
   )
  }
  annot_alpha_init <- as.matrix(annot_alpha_init)
  rownames(annot_alpha_init) <- colnames(A)
  if (is.null(colnames(annot_alpha_init))) colnames(annot_alpha_init) <- step_names
 }
 if (is.null(annot_sigma_sq_alpha_init)) {
  annot_sigma_sq_alpha_init <- rep(1, K_minus_one)
 } else if (!is.numeric(annot_sigma_sq_alpha_init) ||
            length(annot_sigma_sq_alpha_init) != K_minus_one ||
            any(!is.finite(annot_sigma_sq_alpha_init)) ||
            any(annot_sigma_sq_alpha_init <= 0)) {
  stop(
   "annot_sigma_sq_alpha_init must be a positive finite numeric vector of length ",
   K_minus_one, "."
  )
 }
 annot_sigma_sq_alpha_init <- stats::setNames(
  as.numeric(annot_sigma_sq_alpha_init), step_names
 )
 list(
  pi = pi,
  annot_alpha_init = annot_alpha_init,
  annot_sigma_sq_alpha_init = annot_sigma_sq_alpha_init
 )
}
