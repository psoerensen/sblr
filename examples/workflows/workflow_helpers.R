## Shared helpers for sblr example workflows.
##
## These helpers intentionally use base R only. They are meant for example scripts,
## not for package internals.

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

workflow_data_dir <- function() {
  data_dir <- Sys.getenv("SBLR_EXAMPLE_DATA_DIR")
  if (!nzchar(data_dir)) {
    stop(
      "Set SBLR_EXAMPLE_DATA_DIR to a directory containing example data, ",
      "or create Glist/y objects before running this workflow."
    )
  }
  data_dir
}

workflow_should_run_heavy <- function() {
  identical(tolower(Sys.getenv("SBLR_RUN_HEAVY_EXAMPLES")), "true")
}

read_example_glist <- function(filename = "Glist_sparseLD_1k.RDS") {
  path <- file.path(workflow_data_dir(), filename)
  if (!file.exists(path)) {
    stop("Example Glist file not found: ", path)
  }
  readRDS(path)
}

check_fit <- function(fit, require_ld_swap = FALSE, require_chain_summaries = FALSE) {
  chk <- check_stblr_consistency(
    fit,
    require_ld_swap = require_ld_swap,
    require_chain_summaries = require_chain_summaries
  )
  print(chk)
  stopifnot(isTRUE(chk$ok))
  invisible(chk)
}

compact_input <- function(fit) {
  x <- fit$input
  x[c(
    "prior_kernel", "probability_policy", "effect_scale_policy",
    "selection_s", "estimate_selection_s", "updateLDswap", "ld_swap_prob", "ld_swap_r2",
    "ld_swap_moves", "nchains", "keep_chains", "n", "m", "nt", "nit",
    "nburn", "nthin", "seed"
  )]
}

fit_metadata_table <- function(fits) {
  do.call(
    rbind,
    lapply(names(fits), function(model_name) {
      fit <- fits[[model_name]]
      data.frame(
        model_name = model_name,
        family = fit$family %||% NA_character_,
        model = fit$model %||% NA_character_,
        operator = fit$operator %||% NA_character_,
        data_level = fit$data$data_level %||% NA_character_,
        prior_kernel = fit$input$prior_kernel %||% NA_character_,
        probability_policy = fit$input$probability_policy %||% NA_character_,
        updateLDswap = fit$input$updateLDswap %||% FALSE,
        nchains = fit$input$nchains %||% 1L,
        keep_chains = fit$input$keep_chains %||% FALSE,
        n_markers = nrow(fit$dm),
        n_traits = ncol(fit$dm),
        stringsAsFactors = FALSE
      )
    })
  )
}

fit_field_inventory <- function(fits) {
  fields <- c(
    "dm", "bm", "vbs", "vgs", "ves", "vle", "vld", "pi_trace",
    "pi_final", "pi_mean", "chains", "bm_chain_mean_sd",
    "bm_chain_mean_min", "bm_chain_mean_max", "dm_chain_mean_sd",
    "dm_chain_mean_min", "dm_chain_mean_max", "component_probabilities",
    "model_parameters", "diagnostics", "convergence", "convergence_traces",
    "chains", "memory_estimate"
  )

  do.call(
    rbind,
    lapply(names(fits), function(model_name) {
      fit <- fits[[model_name]]
      out <- data.frame(
        model = model_name,
        operator = fit$operator %||% NA_character_,
        data_level = fit$data$data_level %||% NA_character_,
        probability_policy = fit$input$probability_policy %||% NA_character_,
        stringsAsFactors = FALSE
      )
      for (field in fields) {
        out[[paste0("has_", field)]] <- field %in% names(fit) && !is.null(fit[[field]])
      }
      out
    })
  )
}

posterior_signal_summary <- function(fit, model_name) {
  dm <- as.matrix(fit$dm)
  bm <- as.matrix(fit$bm)
  trait_names <- colnames(dm) %||% paste0("T", seq_len(ncol(dm)))

  data.frame(
    model = model_name,
    trait = trait_names,
    n_markers = nrow(dm),
    mean_pip = colMeans(dm, na.rm = TRUE),
    sum_pip = colSums(dm, na.rm = TRUE),
    max_pip = apply(dm, 2, max, na.rm = TRUE),
    n_pip_gt_0_001 = colSums(dm > 0.001, na.rm = TRUE),
    n_pip_gt_0_01 = colSums(dm > 0.01, na.rm = TRUE),
    n_pip_gt_0_05 = colSums(dm > 0.05, na.rm = TRUE),
    n_pip_gt_0_5 = colSums(dm > 0.5, na.rm = TRUE),
    mean_abs_bm = colMeans(abs(bm), na.rm = TRUE),
    max_abs_bm = apply(abs(bm), 2, max, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

summarise_fit_list <- function(fits) {
  do.call(
    rbind,
    lapply(names(fits), function(model_name) {
      posterior_signal_summary(fits[[model_name]], model_name)
    })
  )
}

ld_swap_summary_table <- function(fits) {
  out <- lapply(names(fits), function(model_name) {
    x <- fits[[model_name]]$ld_swap
    if (is.null(x)) {
      return(NULL)
    }
    x <- as.data.frame(x)
    x$model <- model_name
    x$trait <- rownames(x) %||% seq_len(nrow(x))
    x[, c("model", "trait", setdiff(names(x), c("model", "trait"))), drop = FALSE]
  })
  out <- Filter(Negate(is.null), out)
  if (!length(out)) {
    return(data.frame())
  }
  do.call(rbind, out)
}

top_markers <- function(fit, model_name = fit$model %||% "model", trait = 1, n = 20) {
  dm <- as.matrix(fit$dm)
  bm <- as.matrix(fit$bm)
  trait_index <- if (is.character(trait)) match(trait, colnames(dm)) else trait
  if (length(trait_index) != 1L || is.na(trait_index)) {
    stop("trait must be a valid column name or column index.")
  }

  marker <- rownames(dm) %||% paste0("V", seq_len(nrow(dm)))
  trait_name <- colnames(dm)[trait_index] %||% paste0("T", trait_index)

  out <- data.frame(
    model = model_name,
    trait = trait_name,
    marker = marker,
    pip = dm[, trait_index],
    bm = bm[, trait_index],
    abs_bm = abs(bm[, trait_index]),
    stringsAsFactors = FALSE
  )
  out[order(-out$pip, -out$abs_bm), ][seq_len(min(n, nrow(out))), ]
}

marker_aligned_maf <- function(Glist, chr = 1L) {
  idx <- match(Glist$rsidsLD[[chr]], Glist$rsids[[chr]])
  if (anyNA(idx)) {
    stop("Could not align all LD markers to Glist$rsids for chromosome ", chr, ".")
  }
  maf <- Glist$maf[[chr]][idx]
  names(maf) <- Glist$rsidsLD[[chr]]
  maf
}

run_credible_sets <- function(fit, Glist, trait = 1, coverage = 0.95) {
  make_credible_sets(
    fit = fit,
    Glist = Glist,
    trait = trait,
    coverage = coverage,
    min_r2 = 0.5,
    pip_cutoff = 0.001,
    locus_pip_cutoff = 0.01,
    max_locus_distance = 1e6,
    method = "pip"
  )
}
