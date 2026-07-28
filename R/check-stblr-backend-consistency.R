#' Check ST-BLR Backend Output Consistency
#'
#' Performs lightweight structural checks on a fitted ST-BLR object. The checker
#' is intended to verify backend output conventions such as marker summary
#' dimensions, optional chain-stability summaries, compact per-chain summaries,
#' and LD-swap diagnostics. It is not a convergence diagnostic and does not
#' inspect MCMC mixing.
#'
#' @param fit Fitted ST-BLR object.
#' @param require_chain_summaries Logical or `NULL`. If `NULL`, chain summary
#'   fields are required when `fit$input$nchains > 1` and optional otherwise.
#'   If `TRUE`, all six chain summary fields are required. If `FALSE`, they are
#'   optional but validated when present.
#' @param require_chains Logical; require compact per-chain summaries in
#'   `fit$chains`.
#' @param require_ld_swap Logical; require LD-swap diagnostics in
#'   `fit$diagnostics$ld_swap`.
#' @param verbose Logical; retained for API symmetry. Printing is handled by
#'   [print.stblr_backend_check()].
#'
#' @return A list of class `"stblr_backend_check"` with:
#' \describe{
#'   \item{ok}{Logical overall pass/fail indicator.}
#'   \item{checks}{Data frame with `check`, `ok`, and `message` columns.}
#'   \item{metadata}{List of detected backend metadata.}
#' }
#'
#' @examples
#' markers <- paste0("m", 1:3)
#' fit <- list(
#'   dm = matrix(c(0.1, 0.2, 0.3), ncol = 1,
#'               dimnames = list(markers, "trait1")),
#'   bm = matrix(c(0.01, -0.02, 0.03), ncol = 1,
#'               dimnames = list(markers, "trait1")),
#'   input = list(nchains = 1)
#' )
#' check_stblr_consistency(fit)
#'
#' @export
check_stblr_consistency <- function(
    fit,
    require_chain_summaries = NULL,
    require_chains = FALSE,
    require_ld_swap = FALSE,
    verbose = TRUE
) {
  if (!is.list(fit)) {
    stop("fit must be a fitted ST-BLR list-like object.")
  }
  if (!is.null(require_chain_summaries) &&
      (!is.logical(require_chain_summaries) ||
       length(require_chain_summaries) != 1L ||
       is.na(require_chain_summaries))) {
    stop("require_chain_summaries must be TRUE, FALSE, or NULL.")
  }
  if (!is.logical(require_chains) || length(require_chains) != 1L ||
      is.na(require_chains)) {
    stop("require_chains must be TRUE or FALSE.")
  }
  if (!is.logical(require_ld_swap) || length(require_ld_swap) != 1L ||
      is.na(require_ld_swap)) {
    stop("require_ld_swap must be TRUE or FALSE.")
  }
  if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
    stop("verbose must be TRUE or FALSE.")
  }

  checks <- list()
  add_check <- function(check, ok, message) {
    checks[[length(checks) + 1L]] <<- data.frame(
      check = check,
      ok = isTRUE(ok),
      message = as.character(message),
      stringsAsFactors = FALSE
    )
  }

  dm <- .stblr_backend_as_matrix(fit$dm)
  bm <- .stblr_backend_as_matrix(fit$bm)
  input <- fit$input

  add_check("required.dm", !is.null(dm), "fit$dm is present and matrix-like")
  add_check("required.bm", !is.null(bm), "fit$bm is present and matrix-like")
  add_check("required.input", is.list(input), "fit$input is present")

  dims_ok <- !is.null(dm) && !is.null(bm) && identical(dim(dm), dim(bm))
  add_check("dimensions.dm_bm", dims_ok, "dim(fit$dm) matches dim(fit$bm)")

  if (!is.null(dm) && !is.null(bm) && dims_ok) {
    add_check(
      "dimnames.rows.dm_bm",
      .stblr_backend_dimnames_match(rownames(dm), rownames(bm)),
      "rownames(fit$dm) match rownames(fit$bm) when present"
    )
    add_check(
      "dimnames.cols.dm_bm",
      .stblr_backend_dimnames_match(colnames(dm), colnames(bm)),
      "colnames(fit$dm) match colnames(fit$bm) when present"
    )
    add_check(
      "values.dm.range",
      all(is.na(dm) | (dm >= -1e-8 & dm <= 1 + 1e-8)),
      "fit$dm values are within [0, 1] allowing tolerance"
    )
  }

  .stblr_backend_check_optional_trace(
    fit = fit,
    name = "vle",
    dm = dm,
    add_check = add_check
  )
  .stblr_backend_check_optional_trace(
    fit = fit,
    name = "vld",
    dm = dm,
    add_check = add_check
  )
  .stblr_backend_check_optional_trace(
    fit = fit,
    name = "pi_trace",
    dm = dm,
    add_check = add_check
  )

  nchains_missing <- !is.list(input) || is.null(input$nchains)
  nchains <- if (nchains_missing) {
    1L
  } else {
    suppressWarnings(as.integer(input$nchains[1L]))
  }
  nchains_ok <- length(nchains) == 1L && !is.na(nchains) && nchains >= 1L
  add_check(
    "metadata.nchains",
    nchains_ok,
    if (nchains_missing) {
      "fit$input$nchains is missing; assuming 1"
    } else {
      "fit$input$nchains is a positive integer"
    }
  )
  if (!nchains_ok) nchains <- NA_integer_

  has_chain_summaries <- all(vapply(
    c("dm_chain_mean_sd", "dm_chain_mean_min", "dm_chain_mean_max", "bm_chain_mean_sd", "bm_chain_mean_min", "bm_chain_mean_max"),
    function(nm) !is.null(fit[[nm]]),
    logical(1)
  ))
  summaries_required <- if (is.null(require_chain_summaries)) {
    isTRUE(nchains_ok && nchains > 1L)
  } else {
    isTRUE(require_chain_summaries)
  }
  .stblr_backend_check_chain_summaries(
    fit = fit,
    dm = dm,
    bm = bm,
    nchains = nchains,
    summaries_required = summaries_required,
    add_check = add_check
  )

  .stblr_backend_check_chains(
    fit = fit,
    dm = dm,
    bm = bm,
    nchains = nchains,
    require_chains = require_chains,
    add_check = add_check
  )

  .stblr_backend_check_ld_swap(
    fit = fit,
    require_ld_swap = require_ld_swap,
    add_check = add_check
  )

  .stblr_backend_check_component_probabilities(
    fit = fit,
    dm = dm,
    add_check = add_check
  )

  .stblr_backend_check_maf_effect_s(
    fit = fit,
    dm = dm,
    add_check = add_check
  )

  .stblr_backend_check_backend_specific(
    fit = fit,
    input = input,
    add_check = add_check
  )

  checks_df <- if (length(checks)) {
    do.call(rbind, checks)
  } else {
    data.frame(check = character(), ok = logical(), message = character(),
               stringsAsFactors = FALSE)
  }
  rownames(checks_df) <- NULL

  metadata <- list(
    nmarkers = if (!is.null(dm)) nrow(dm) else NA_integer_,
    ntraits = if (!is.null(dm)) ncol(dm) else NA_integer_,
    nchains = nchains,
    scheduled = if (is.list(input) && !is.null(input$scheduled)) input$scheduled else NA,
    backend = if (is.list(input) && !is.null(input$backend)) input$backend else NA,
    keep_chains = if (is.list(input) && !is.null(input$keep_chains)) input$keep_chains else NA,
    updateLDswap = if (is.list(input) && !is.null(input$updateLDswap)) input$updateLDswap else NA,
    has_chain_summaries = has_chain_summaries,
    has_chains = !is.null(fit$chains),
    has_ld_swap = is.list(fit$diagnostics) &&
      (!is.null(fit$diagnostics$ld_swap) ||
       (is.list(fit$diagnostics$native) &&
        !is.null(fit$diagnostics$native$ld_swap))),
    nchains_missing = nchains_missing
  )

  out <- list(
    ok = all(checks_df$ok),
    checks = checks_df,
    metadata = metadata
  )
  class(out) <- c("stblr_backend_check", class(out))
  out
}

#' Print an ST-BLR Backend Consistency Check
#'
#' @param x Object returned by [check_stblr_consistency()].
#' @param ... Unused.
#'
#' @return Invisibly returns `x`.
#'
#' @export
print.stblr_backend_check <- function(x, ...) {
  cat("ST-BLR backend consistency check: ",
      if (isTRUE(x$ok)) "OK" else "FAILED", "\n", sep = "")

  meta <- x$metadata
  cat(
    "  nmarkers=", .stblr_backend_meta_value(meta$nmarkers),
    ", ntraits=", .stblr_backend_meta_value(meta$ntraits),
    ", nchains=", .stblr_backend_meta_value(meta$nchains),
    ", scheduled=", .stblr_backend_meta_value(meta$scheduled),
    ", backend=", .stblr_backend_meta_value(meta$backend),
    ", keep_chains=", .stblr_backend_meta_value(meta$keep_chains),
    ", updateLDswap=", .stblr_backend_meta_value(meta$updateLDswap),
    "\n",
    sep = ""
  )

  failed <- x$checks[!x$checks$ok, , drop = FALSE]
  if (nrow(failed) == 0L) {
    cat("  All checks passed.\n")
  } else {
    cat("  Failed checks:\n")
    for (i in seq_len(nrow(failed))) {
      cat("  - ", failed$check[i], ": ", failed$message[i], "\n", sep = "")
    }
  }
  invisible(x)
}

.stblr_backend_as_matrix <- function(x) {
  if (is.null(x)) return(NULL)
  if (is.matrix(x)) return(x)
  if (is.data.frame(x)) return(as.matrix(x))
  NULL
}

.stblr_backend_dimnames_match <- function(a, b) {
  if (is.null(a) || is.null(b)) return(TRUE)
  identical(a, b)
}

.stblr_backend_meta_value <- function(x) {
  if (length(x) == 0L || all(is.na(x))) return("NA")
  paste(as.character(x), collapse = ",")
}

.stblr_backend_matrix_values_ok <- function(x) {
  all(is.na(x) | is.finite(x))
}

.stblr_backend_check_optional_trace <- function(fit, name, dm, add_check) {
  if (is.null(fit[[name]])) return(invisible(NULL))

  x <- .stblr_backend_as_matrix(fit[[name]])
  add_check(
    paste0("trace.", name, ".matrix"),
    !is.null(x),
    paste0("fit$", name, " is matrix-like")
  )
  if (is.null(x)) return(invisible(NULL))

  add_check(
    paste0("trace.", name, ".nonempty"),
    nrow(x) > 0L && ncol(x) > 0L,
    paste0("fit$", name, " has at least one iteration and one trait")
  )
  if (!is.null(dm)) {
    add_check(
      paste0("trace.", name, ".ntraits"),
      ncol(x) == ncol(dm),
      paste0("ncol(fit$", name, ") matches number of traits")
    )
    add_check(
      paste0("trace.", name, ".colnames"),
      .stblr_backend_dimnames_match(colnames(x), colnames(dm)),
      paste0("colnames(fit$", name, ") match trait names when present")
    )
  }
  add_check(
    paste0("trace.", name, ".finite"),
    .stblr_backend_matrix_values_ok(x),
    paste0("fit$", name, " contains only finite values or NA")
  )
  invisible(x)
}

.stblr_backend_check_summary_matrix <- function(
    fit, name, target, add_check
) {
  x <- .stblr_backend_as_matrix(fit[[name]])
  add_check(
    paste0("chain_summary.", name, ".matrix"),
    !is.null(x),
    paste0("fit$", name, " is matrix-like")
  )
  if (is.null(x) || is.null(target)) return(NULL)

  add_check(
    paste0("chain_summary.", name, ".dim"),
    identical(dim(x), dim(target)),
    paste0("dim(fit$", name, ") matches target dimensions")
  )
  add_check(
    paste0("chain_summary.", name, ".rownames"),
    .stblr_backend_dimnames_match(rownames(x), rownames(target)),
    paste0("rownames(fit$", name, ") match target rownames when present")
  )
  add_check(
    paste0("chain_summary.", name, ".colnames"),
    .stblr_backend_dimnames_match(colnames(x), colnames(target)),
    paste0("colnames(fit$", name, ") match target colnames when present")
  )
  add_check(
    paste0("chain_summary.", name, ".finite"),
    .stblr_backend_matrix_values_ok(x),
    paste0("fit$", name, " contains only finite values or NA")
  )
  x
}

.stblr_backend_check_chain_summaries <- function(
    fit, dm, bm, nchains, summaries_required, add_check
) {
  fields <- c("dm_chain_mean_sd", "dm_chain_mean_min", "dm_chain_mean_max", "bm_chain_mean_sd", "bm_chain_mean_min", "bm_chain_mean_max")
  present <- vapply(fields, function(nm) !is.null(fit[[nm]]), logical(1))

  add_check(
    "chain_summaries.required",
    !summaries_required || all(present),
    if (summaries_required) {
      "all chain summary fields are present"
    } else {
      "chain summary fields are optional"
    }
  )

  dm_chain_mean_sd <- dm_chain_mean_min <- dm_chain_mean_max <- bm_chain_mean_sd <- bm_chain_mean_min <- bm_chain_mean_max <- NULL
  if (present[["dm_chain_mean_sd"]]) dm_chain_mean_sd <- .stblr_backend_check_summary_matrix(fit, "dm_chain_mean_sd", dm, add_check)
  if (present[["dm_chain_mean_min"]]) dm_chain_mean_min <- .stblr_backend_check_summary_matrix(fit, "dm_chain_mean_min", dm, add_check)
  if (present[["dm_chain_mean_max"]]) dm_chain_mean_max <- .stblr_backend_check_summary_matrix(fit, "dm_chain_mean_max", dm, add_check)
  if (present[["bm_chain_mean_sd"]]) bm_chain_mean_sd <- .stblr_backend_check_summary_matrix(fit, "bm_chain_mean_sd", bm, add_check)
  if (present[["bm_chain_mean_min"]]) bm_chain_mean_min <- .stblr_backend_check_summary_matrix(fit, "bm_chain_mean_min", bm, add_check)
  if (present[["bm_chain_mean_max"]]) bm_chain_mean_max <- .stblr_backend_check_summary_matrix(fit, "bm_chain_mean_max", bm, add_check)

  tol <- 1e-8
  if (!is.null(dm_chain_mean_sd)) {
    add_check("chain_summary.dm_chain_mean_sd.nonnegative",
              all(is.na(dm_chain_mean_sd) | dm_chain_mean_sd >= -tol),
              "fit$dm_chain_mean_sd is non-negative")
  }
  if (!is.null(bm_chain_mean_sd)) {
    add_check("chain_summary.bm_chain_mean_sd.nonnegative",
              all(is.na(bm_chain_mean_sd) | bm_chain_mean_sd >= -tol),
              "fit$bm_chain_mean_sd is non-negative")
  }
  if (!is.null(dm) && !is.null(dm_chain_mean_min) && !is.null(dm_chain_mean_max) &&
      identical(dim(dm), dim(dm_chain_mean_min)) && identical(dim(dm), dim(dm_chain_mean_max))) {
    add_check("chain_summary.dm.bounds",
              all(is.na(dm) | is.na(dm_chain_mean_min) | dm_chain_mean_min <= dm + tol) &&
                all(is.na(dm) | is.na(dm_chain_mean_max) | dm <= dm_chain_mean_max + tol),
              "fit$dm lies between fit$dm_chain_mean_min and fit$dm_chain_mean_max")
    add_check("chain_summary.dm_chain_mean_minmax.range",
              all(is.na(dm_chain_mean_min) | (dm_chain_mean_min >= -tol & dm_chain_mean_min <= 1 + tol)) &&
                all(is.na(dm_chain_mean_max) | (dm_chain_mean_max >= -tol & dm_chain_mean_max <= 1 + tol)),
              "fit$dm_chain_mean_min and fit$dm_chain_mean_max are within [0, 1] allowing tolerance")
  }
  if (!is.null(bm) && !is.null(bm_chain_mean_min) && !is.null(bm_chain_mean_max) &&
      identical(dim(bm), dim(bm_chain_mean_min)) && identical(dim(bm), dim(bm_chain_mean_max))) {
    add_check("chain_summary.bm.bounds",
              all(is.na(bm) | is.na(bm_chain_mean_min) | bm_chain_mean_min <= bm + tol) &&
                all(is.na(bm) | is.na(bm_chain_mean_max) | bm <= bm_chain_mean_max + tol),
              "fit$bm lies between fit$bm_chain_mean_min and fit$bm_chain_mean_max")
  }
  if (identical(as.integer(nchains), 1L) &&
      all(present) && !is.null(dm) && !is.null(bm)) {
    add_check("chain_summary.single.dm_chain_mean_sd_zero",
              all(abs(dm_chain_mean_sd) <= tol, na.rm = TRUE),
              "single-chain fit$dm_chain_mean_sd is approximately zero")
    add_check("chain_summary.single.bm_chain_mean_sd_zero",
              all(abs(bm_chain_mean_sd) <= tol, na.rm = TRUE),
              "single-chain fit$bm_chain_mean_sd is approximately zero")
    add_check("chain_summary.single.dm_equal",
              isTRUE(all.equal(dm_chain_mean_min, dm, tolerance = tol, check.attributes = FALSE)) &&
                isTRUE(all.equal(dm_chain_mean_max, dm, tolerance = tol, check.attributes = FALSE)),
              "single-chain fit$dm_chain_mean_min and fit$dm_chain_mean_max equal fit$dm")
    add_check("chain_summary.single.bm_equal",
              isTRUE(all.equal(bm_chain_mean_min, bm, tolerance = tol, check.attributes = FALSE)) &&
                isTRUE(all.equal(bm_chain_mean_max, bm, tolerance = tol, check.attributes = FALSE)),
              "single-chain fit$bm_chain_mean_min and fit$bm_chain_mean_max equal fit$bm")
  }
}

.stblr_backend_check_chains <- function(
    fit, dm, bm, nchains, require_chains, add_check
) {
  has_chains <- !is.null(fit$chains)
  add_check(
    "chains.required",
    !require_chains || has_chains,
    if (require_chains) "fit$chains is present" else "fit$chains is optional"
  )
  if (!has_chains || is.null(dm) || is.null(bm)) return(invisible(NULL))

  trait_names <- colnames(dm) %||% paste0("trait", seq_len(ncol(dm)))
  chains <- fit$chains
  expected_tasks <- as.integer(nchains) * ncol(dm)
  add_check("chains.task_count", length(chains) == expected_tasks,
            "fit$chains has one flat record per trait-by-chain task")

  for (tt in seq_len(ncol(dm))) {
    cname <- trait_names[tt]
    trait_chains <- Filter(function(x) {
      is.list(x) && identical(as.integer(x$trait_index), as.integer(tt))
    }, chains)
    add_check(
      paste0("chains.", cname, ".chain_count"),
      is.list(trait_chains) && length(trait_chains) == nchains,
      paste0("fit$chains has nchains records for ", cname)
    )
    if (!is.list(trait_chains)) next

    dm_chain <- list()
    bm_chain <- list()
    component_probabilities_chain <- list()
    for (cc in seq_along(trait_chains)) {
      ch <- trait_chains[[cc]]
      if (!is.list(ch)) {
        add_check(paste0("chains.", cname, ".chain", cc, ".list"),
                  FALSE, "chain entry is list-like")
        next
      }
      marker <- if (is.list(ch$marker)) ch$marker else list()
      chain_dm <- ch$dm %||% marker$dm
      chain_bm <- ch$bm %||% marker$bm
      chain_component_probabilities <-
        ch$component_probabilities %||% marker$component_probabilities
      if (!is.null(chain_dm)) {
        dm_chain[[length(dm_chain) + 1L]] <- as.numeric(chain_dm)
        add_check(paste0("chains.", cname, ".chain", cc, ".dm_length"),
                  length(chain_dm) == nrow(dm),
                  "chain dm length matches number of markers")
      }
      if (!is.null(chain_bm)) {
        bm_chain[[length(bm_chain) + 1L]] <- as.numeric(chain_bm)
        add_check(paste0("chains.", cname, ".chain", cc, ".bm_length"),
                  length(chain_bm) == nrow(bm),
                  "chain bm length matches number of markers")
      }
      if (!is.null(chain_component_probabilities)) {
        cp <- as.matrix(chain_component_probabilities)
        component_probabilities_chain[[length(component_probabilities_chain) + 1L]] <- cp
        add_check(paste0("chains.", cname, ".chain", cc, ".component_probabilities_rows"),
                  nrow(cp) == nrow(dm),
                  "chain component_probabilities row count matches number of markers")
      }
    }
    if (length(dm_chain) == nchains &&
        all(vapply(dm_chain, length, integer(1)) == nrow(dm))) {
      dm_mat <- do.call(cbind, dm_chain)
      add_check(paste0("chains.", cname, ".dm_mean"),
                isTRUE(all.equal(rowMeans(dm_mat), dm[, tt], tolerance = 1e-8,
                                 check.attributes = FALSE)),
                "mean chain dm equals fit$dm")
    }
    if (length(bm_chain) == nchains &&
        all(vapply(bm_chain, length, integer(1)) == nrow(bm))) {
      bm_mat <- do.call(cbind, bm_chain)
      add_check(paste0("chains.", cname, ".bm_mean"),
                isTRUE(all.equal(rowMeans(bm_mat), bm[, tt], tolerance = 1e-8,
                                 check.attributes = FALSE)),
                "mean chain bm equals fit$bm")
    }
    if (!is.null(fit$component_probabilities) && !is.null(fit$component_probabilities[[cname]]) &&
        length(component_probabilities_chain) == nchains) {
      target_cp <- as.matrix(fit$component_probabilities[[cname]])
      same_dim <- all(vapply(
        component_probabilities_chain,
        function(x) identical(dim(x), dim(target_cp)),
        logical(1)
      ))
      add_check(paste0("chains.", cname, ".component_probabilities_dim"),
                same_dim,
                "chain component_probabilities dimensions match fit$component_probabilities")
      if (same_dim) {
        cp_mean <- Reduce(`+`, component_probabilities_chain) / length(component_probabilities_chain)
        add_check(paste0("chains.", cname, ".component_probabilities_mean"),
                  isTRUE(all.equal(cp_mean, target_cp, tolerance = 1e-8,
                                   check.attributes = FALSE)),
                  "mean chain component_probabilities equals fit$component_probabilities")
      }
    }
  }
  invisible(NULL)
}

.stblr_backend_check_ld_swap <- function(fit, require_ld_swap, add_check) {
  diagnostics <- if (is.list(fit$diagnostics)) fit$diagnostics else list()
  native <- if (is.list(diagnostics$native)) diagnostics$native else list()
  ld_swap <- diagnostics$ld_swap %||% native$ld_swap
  ld_swap_chains <- diagnostics$ld_swap_chains %||% native$ld_swap_chains
  has_ld_swap <- !is.null(ld_swap)
  add_check(
    "ld_swap.required",
    !require_ld_swap || has_ld_swap,
    if (require_ld_swap) "fit$diagnostics$ld_swap is present" else "fit$diagnostics$ld_swap is optional"
  )
  if (has_ld_swap) {
    .stblr_backend_validate_ld_swap_table(ld_swap, "ld_swap", add_check)
  }
  if (!is.null(ld_swap_chains)) {
    if (is.list(ld_swap_chains) && !is.data.frame(ld_swap_chains)) {
      for (nm in names(ld_swap_chains)) {
        .stblr_backend_validate_ld_swap_table(
          ld_swap_chains[[nm]],
          paste0("ld_swap_chains.", nm),
          add_check
        )
      }
    } else {
      .stblr_backend_validate_ld_swap_table(
        ld_swap_chains,
        "ld_swap_chains",
        add_check
      )
    }
  }
}

.stblr_backend_validate_ld_swap_table <- function(x, label, add_check) {
  tab <- as.data.frame(x)
  required <- c("attempted", "accepted", "acceptance_rate")
  has_cols <- all(required %in% names(tab))
  add_check(paste0(label, ".columns"), has_cols,
            paste0(label, " has attempted, accepted, and acceptance_rate"))
  if (!has_cols) return(invisible(NULL))

  attempted <- suppressWarnings(as.numeric(tab$attempted))
  accepted <- suppressWarnings(as.numeric(tab$accepted))
  rate <- suppressWarnings(as.numeric(tab$acceptance_rate))
  add_check(paste0(label, ".counts"),
            all(is.finite(attempted) & is.finite(accepted) &
                  attempted >= -1e-8 & accepted >= -1e-8 &
                  accepted <= attempted + 1e-8),
            paste0(label, " has attempted >= accepted >= 0"))
  add_check(paste0(label, ".acceptance_rate"),
            all(is.finite(rate) & rate >= -1e-8 & rate <= 1 + 1e-8),
            paste0(label, " acceptance_rate is finite and within [0, 1]"))
  positive <- attempted > 0
  if (any(positive)) {
    add_check(paste0(label, ".acceptance_rate_ratio"),
              isTRUE(all.equal(rate[positive], accepted[positive] / attempted[positive],
                               tolerance = 1e-8, check.attributes = FALSE)),
              paste0(label, " acceptance_rate matches accepted / attempted"))
  }
  invisible(NULL)
}

.stblr_backend_check_component_probabilities <- function(fit, dm, add_check) {
  cp <- fit$component_probabilities
  if (is.null(cp)) return(invisible(NULL))

  if (!is.list(cp) || is.null(names(cp)) || any(!nzchar(names(cp)))) {
    add_check(
      "component_probabilities.named_list", FALSE,
      "fit$component_probabilities is a named list of marker-by-component matrices by trait"
    )
    return(invisible(NULL))
  }
  add_check(
    "component_probabilities.named_list", TRUE,
    "fit$component_probabilities is a named list of marker-by-component matrices by trait"
  )

  tol <- 1e-6
  for (trait in names(cp)) {
    x <- .stblr_backend_as_matrix(cp[[trait]])
    label <- paste0("component_probabilities.", trait)
    add_check(
      paste0(label, ".matrix"), !is.null(x),
      paste0("fit$component_probabilities[['", trait, "']] is matrix-like")
    )
    if (is.null(x)) next

    if (!is.null(dm)) {
      add_check(
        paste0(label, ".nrow"), nrow(x) == nrow(dm),
        paste0("fit$component_probabilities[['", trait, "']] has nrow(fit$dm) rows")
      )
    }
    add_check(
      paste0(label, ".range"),
      all(is.na(x) | (x >= -tol & x <= 1 + tol)),
      paste0("fit$component_probabilities[['", trait, "']] values are within [0, 1]")
    )
    row_sums <- rowSums(x)
    add_check(
      paste0(label, ".rowsums"),
      all(is.na(row_sums) | abs(row_sums - 1) <= 1e-4),
      paste0("fit$component_probabilities[['", trait, "']] rows sum to 1")
    )
  }
  invisible(NULL)
}

.stblr_backend_check_maf_effect_s <- function(fit, dm, add_check) {
  has_any <- !is.null(fit$maf_effect_s) || !is.null(fit$maf_effect_s_trace)
  if (!has_any) return(invisible(NULL))

  if (!is.null(fit$maf_effect_s_trace)) {
    .stblr_backend_check_optional_trace(
      fit = fit, name = "maf_effect_s_trace", dm = dm, add_check = add_check
    )
  }
  if (!is.null(fit$maf_effect_s)) {
    x <- fit$maf_effect_s
    add_check(
      "maf_effect_s.finite", all(is.finite(x)),
      "fit$maf_effect_s contains only finite values"
    )
    if (!is.null(dm)) {
      add_check(
        "maf_effect_s.length", length(x) == ncol(dm),
        "length(fit$maf_effect_s) matches number of traits"
      )
    }
  }
  invisible(NULL)
}

.stblr_backend_check_backend_specific <- function(fit, input, add_check) {
  backend <- if (is.list(input) && !is.null(input$backend)) {
    input$backend
  } else {
    NA_character_
  }
  if (length(backend) != 1L || is.na(backend)) return(invisible(NULL))

  expect_fields <- switch(
    backend,
    csr_prior_bayesc = "prior",
    csr_group_bayesc = c("group", "group_pi", "group_vb_multiplier"),
    csr_annot_bayesc = c("eta_pi", "eta_vb"),
    csr_sbayesrc = c("alpha", "sigmaSqAlpha", "component_probabilities"),
    csr_bayesr = "component_probabilities",
    bed_bayesr = "component_probabilities",
    bed_bayesrc = c("alpha", "sigmaSqAlpha", "component_probabilities"),
    character()
  )
  for (nm in expect_fields) {
    add_check(
      paste0("backend.", backend, ".", nm),
      !is.null(fit[[nm]]),
      paste0("fit$", nm, " is present for backend '", backend, "'")
    )
  }
  invisible(NULL)
}
