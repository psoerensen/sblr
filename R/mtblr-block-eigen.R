.mtblr_block_eigen_reference_error <- function(detail = NULL) {
  message <- paste(
    "mtblr_block_eigen() currently requires sufficient statistics produced",
    "from the same selected BED genotypes used to construct the block operator.",
    "External GWAS/reference-panel projection is not yet supported."
  )
  if (!is.null(detail)) message <- paste0(message, " ", detail)
  stop(message, call. = FALSE)
}

.mtblr_block_eigen_paths <- function(paths, field) {
  if (is.null(paths) || !length(paths) || anyNA(paths) || any(!nzchar(paths))) {
    .mtblr_block_eigen_reference_error(paste0(field, " is missing."))
  }
  normalizePath(as.character(paths), winslash = "/", mustWork = TRUE)
}

.mtblr_block_eigen_rows <- function(rows, n_bed) {
  if (is.null(rows)) return(seq_len(n_bed))
  rows <- as.integer(rows)
  if (!length(rows) || anyNA(rows) || any(rows < 1L | rows > n_bed) ||
      anyDuplicated(rows)) {
    .mtblr_block_eigen_reference_error("Selected BED rows are invalid.")
  }
  rows
}

.mtblr_block_eigen_equal <- function(x, y, label) {
  if (!identical(x, y)) {
    .mtblr_block_eigen_reference_error(paste0(label, " differs."))
  }
}

.mtblr_block_eigen_reference <- function(Glist, provenance) {
  required <- c("bed_files", "n_bed", "cls", "af", "marker_ids",
                "marker_metadata", "source", "scale", "sample_size")
  if (!is.list(provenance) ||
      any(vapply(required, function(x) is.null(provenance[[x]]), logical(1)))) {
    .mtblr_block_eigen_reference_error("BED construction provenance is incomplete.")
  }
  if (!identical(provenance$source, "make_summary_stats")) {
    .mtblr_block_eigen_reference_error("stats source is not make_summary_stats.")
  }
  if (!identical(.mtblr_normalize_scale(provenance$scale),
                 "standardized_genotype")) {
    .mtblr_block_eigen_reference_error("The standardized-genotype scale is required.")
  }
  if (!is.list(Glist) || is.null(Glist$bedfiles) || is.null(Glist$n) ||
      is.null(Glist$rsids) || is.null(Glist$af)) {
    .mtblr_block_eigen_reference_error("Glist lacks BED, sample, marker, or frequency data.")
  }
  n_bed <- as.integer(Glist$n)
  if (length(n_bed) != 1L || is.na(n_bed) || n_bed < 1L) {
    .mtblr_block_eigen_reference_error("Glist$n is invalid.")
  }
  .mtblr_block_eigen_equal(as.integer(provenance$n_bed), n_bed,
                           "Original BED sample count")
  stats_paths <- .mtblr_block_eigen_paths(provenance$bed_files,
                                          "Stats BED files")
  all_paths <- as.character(Glist$bedfiles)
  available <- !is.na(all_paths) & nzchar(all_paths)
  normalized_all <- rep(NA_character_, length(all_paths))
  normalized_all[available] <- .mtblr_block_eigen_paths(
    all_paths[available], "Glist BED files")
  chr <- match(stats_paths, normalized_all)
  if (anyNA(chr) || anyDuplicated(chr)) {
    .mtblr_block_eigen_reference_error("BED files or BED-file order differ.")
  }
  cls <- provenance$cls
  if (!is.list(cls)) cls <- list(cls)
  cls <- lapply(cls, as.integer)
  if (length(cls) != length(chr) || any(lengths(cls) == 0L) ||
      anyNA(unlist(cls, use.names = FALSE)) ||
      any(unlist(cls, use.names = FALSE) < 1L)) {
    .mtblr_block_eigen_reference_error("Selected BED marker columns are invalid.")
  }
  af <- provenance$af
  if (!is.list(af)) af <- list(af)
  af <- lapply(af, as.numeric)
  if (!identical(lengths(af), lengths(cls)) ||
      any(!is.finite(unlist(af, use.names = FALSE)))) {
    .mtblr_block_eigen_reference_error("Allele-frequency provenance is invalid.")
  }
  reference_cls <- Map(function(cc, columns) {
    if (!is.null(Glist$rsidsLD) && cc <= length(Glist$rsidsLD) &&
        !is.null(Glist$rsidsLD[[cc]])) {
      resolved <- match(Glist$rsidsLD[[cc]], Glist$rsids[[cc]])
      if (anyNA(resolved)) {
        .mtblr_block_eigen_reference_error(
          "Glist$rsidsLD cannot be resolved to BED marker columns.")
      }
      as.integer(resolved)
    } else columns
  }, chr, cls)
  if (!identical(sort(unlist(cls, use.names = FALSE)),
                 sort(unlist(reference_cls, use.names = FALSE)))) {
    .mtblr_block_eigen_reference_error("Selected BED marker columns differ.")
  }
  glist_ids <- unlist(Map(function(cc, columns) {
    if (cc > length(Glist$rsids) || any(columns > length(Glist$rsids[[cc]]))) {
      .mtblr_block_eigen_reference_error("Selected BED marker columns differ.")
    }
    as.character(Glist$rsids[[cc]][columns])
  }, chr, reference_cls), use.names = FALSE)
  glist_af <- unlist(Map(function(cc, columns) {
    if (cc > length(Glist$af) || any(columns > length(Glist$af[[cc]]))) {
      .mtblr_block_eigen_reference_error("Selected BED marker columns differ.")
    }
    as.numeric(Glist$af[[cc]][columns])
  }, chr, reference_cls), use.names = FALSE)
  marker_ids <- as.character(provenance$marker_ids)
  if (length(marker_ids) != length(glist_ids) ||
      !setequal(marker_ids, glist_ids)) {
    .mtblr_block_eigen_reference_error("Marker IDs differ.")
  }
  reorder <- match(glist_ids, marker_ids)
  stats_cls <- unlist(cls, use.names = FALSE)
  .mtblr_block_eigen_equal(as.integer(stats_cls[reorder]),
                           as.integer(unlist(reference_cls, use.names = FALSE)),
                           "Selected BED marker columns")
  .mtblr_block_eigen_equal(unname(unlist(af, use.names = FALSE)[reorder]),
                           unname(glist_af), "Allele frequencies")
  rows <- .mtblr_block_eigen_rows(provenance$rows, n_bed)
  reference_rows <- if (!is.null(Glist$sparseLD) &&
                        "rows" %in% names(Glist$sparseLD)) {
    .mtblr_block_eigen_rows(Glist$sparseLD$rows, n_bed)
  } else if (!is.null(Glist$idsLD)) {
    resolved <- match(Glist$idsLD, Glist$ids)
    if (anyNA(resolved)) {
      .mtblr_block_eigen_reference_error(
        "Glist$idsLD cannot be resolved to BED rows.")
    }
    as.integer(resolved)
  } else {
    seq_len(n_bed)
  }
  .mtblr_block_eigen_equal(rows, reference_rows, "Selected BED rows or row order")
  selected_count <- length(rows)
  .mtblr_block_eigen_equal(as.integer(provenance$sample_size),
                           as.integer(selected_count),
                           "Analysis sample size")
  marker_metadata <- .mtblr_marker_metadata(marker_ids,
                                             provenance$marker_metadata)
  marker_metadata <- marker_metadata[reorder, , drop = FALSE]
  expected_metadata <- data.frame(
    marker_id = glist_ids,
    chromosome_or_file = rep(chr, lengths(reference_cls)),
    bed_column = unlist(reference_cls, use.names = FALSE),
    allele_frequency = glist_af,
    stringsAsFactors = FALSE
  )
  for (field in intersect(names(expected_metadata), names(marker_metadata))) {
    .mtblr_block_eigen_equal(unname(marker_metadata[[field]]),
                             unname(expected_metadata[[field]]),
                             paste("Marker metadata", field))
  }
  sparse <- Glist$sparseLD
  if (!is.null(sparse)) {
    if (!is.null(sparse$bed_files)) {
      .mtblr_block_eigen_equal(
        .mtblr_block_eigen_paths(sparse$bed_files, "Glist sparseLD BED files"),
        stats_paths, "Glist sparseLD BED files")
    }
    if (!is.null(sparse$cls)) {
      sparse_cls <- unname(lapply(sparse$cls, as.integer))
      .mtblr_block_eigen_equal(sparse_cls, unname(reference_cls),
                               "Glist sparseLD marker columns")
    }
    if (!is.null(sparse$af)) {
      .mtblr_block_eigen_equal(
        unname(unlist(sparse$af, use.names = FALSE)),
        unname(unlist(af, use.names = FALSE)),
        "Glist sparseLD allele frequencies")
    }
  }
  list(
    bed_files = stats_paths,
    n_bed = n_bed,
    cls = unname(reference_cls),
    rows = if (is.null(provenance$rows)) NULL else as.integer(provenance$rows),
    af = unname(unlist(af, use.names = FALSE)),
    marker_ids = glist_ids,
    marker_metadata = marker_metadata,
    scale = "standardized_genotype",
    source = "same_bed_by_construction",
    by_construction = TRUE,
    selected_row_count = selected_count,
    provenance_key = list(
      bed_files = stats_paths, n_bed = n_bed, cls = unname(reference_cls),
      rows = rows,
      af = unname(unlist(af, use.names = FALSE)), marker_ids = glist_ids,
      marker_metadata = marker_metadata, scale = "standardized_genotype")
  )
}

.mtblr_block_eigen_specs <- function(block_start, m, nt) {
  values <- if (is.list(block_start)) block_start else list(block_start)
  if (!(length(values) %in% c(1L, nt))) {
    stop("block_start must contain one specification or one per trait.",
         call. = FALSE)
  }
  values <- lapply(values, function(x) {
    if (!is.numeric(x) || !length(x) || any(!is.finite(x)) ||
        any(x != as.integer(x))) {
      stop("Every block_start must be a nonempty integer vector.",
           call. = FALSE)
    }
    x <- as.integer(x)
    if (x[1L] != 1L || any(x < 1L | x > m) ||
        any(diff(x) <= 0L)) {
      stop("Public block_start must begin at 1 and be strictly ascending in [1, m].",
           call. = FALSE)
    }
    x
  })
  list(values = if (length(values) == 1L) rep(values, nt) else values,
       input_count = length(values))
}

.mtblr_block_eigen_vector <- function(x, nt, name, choices = NULL) {
  if (!(length(x) %in% c(1L, nt))) {
    stop(name, " must have length one or nt.", call. = FALSE)
  }
  if (!is.null(choices)) {
    x <- as.character(x)
    if (anyNA(x) || any(!x %in% choices)) {
      stop(name, " contains an unsupported value.", call. = FALSE)
    }
  } else {
    x <- as.numeric(x)
    if (any(!is.finite(x)) || any(x < 0)) {
      stop(name, " must be finite and nonnegative.", call. = FALSE)
    }
  }
  list(values = rep(x, length.out = nt), input_count = length(x))
}

.mtblr_block_eigen_same_provenance <- function(references) {
  all(vapply(references[-1L], function(x)
    identical(x$provenance_key, references[[1L]]$provenance_key),
    logical(1)))
}

#' Joint multivariate BayesC, BayesR, and SBayesR with block-eigen operators
#'
#' Fits the aligned joint multivariate mixture model using canonical
#' block-filtered operators built from PLINK BED data. The function currently
#' requires sufficient statistics produced from the same selected BED
#' genotypes used to construct the operator. External GWAS/reference-panel
#' projection is not supported.
#'
#' Hard truncation projects `wy`; fixed and Ledoit-Wolf ridge filters leave it
#' unchanged. Returned `wy` is exactly the value consumed by the sampler.
#' Runtime storage contains reconstructed dense blocks, not low-rank factors.
#' Blocks are mandatory and `block_start` uses one-based R indexing. Sample
#' overlap is not modeled and residual covariance is diagonal.
#'
#' @inheritParams mtblr_csr
#' @param Glist One genotype list or one genotype list per trait.
#' @param block_start One mandatory one-based block-start vector, a one-entry
#'   list, or a list with one vector per trait.
#' @param operator_sharing Operator ownership policy: automatic, shared, or
#'   trait-specific.
#' @param eigen_filter One filter name or one per trait: `"hard_truncate"`,
#'   `"ridge_fixed"`, or `"ridge_lw"`.
#' @param eigen_tau Nonnegative hard-truncation threshold, scalar or per trait.
#' @param eigen_eta Nonnegative fixed-ridge parameter, scalar or per trait.
#' @param summary_reference Must be `"same_bed_by_construction"`.
#' @return An object of class `mtblr_fit` with block diagnostics and operator
#'   metadata.
#' @export
mtblr_block_eigen <- function(
  stats, Glist, block_start,
  operator_sharing = c("auto", "shared", "trait_specific"),
  eigen_filter = "hard_truncate", eigen_tau = 0.01, eigen_eta = 0,
  summary_reference = "same_bed_by_construction", trait_metadata = NULL,
  marker_policy = c("strict", "reorder_stats"),
  sample_overlap = "not_modeled", method = "sbayesc", n = NULL, sets = NULL,
  beta = NULL, b = NULL, state = NULL, h2 = 0.5, pi = 0.001,
  models = NULL, pimodels = NULL, mixture_var = NULL, joint_pi = NULL,
  joint_pi_prior = NULL, component = NULL, selection_s = NULL,
  selection_maf = NULL, allow_reference_maf_for_selection_s = FALSE,
  estimate_selection_s = FALSE, selection_s_init = NULL,
  selection_s_prior = NULL, selection_s_proposal_sd = NULL,
  vg = NULL, vb = NULL, ve = NULL, ssb_prior = NULL, sse_prior = NULL,
  updateB = TRUE, updateE = TRUE, updatePi = TRUE, nub = 4, nue = 4,
  nit = 1000, nburn = 500, nthin = 1, seed = 1,
  nchains = 1L, ncores = 1L, chain_seeds = NULL,
  keep_chains = FALSE, convergence = c("auto", "none", "core"),
  convergence_control = NULL, memory_warning_gb = 8, verbose = FALSE
) {
  semantics <- .mtblr_resolve_public_method(method, "block_eigen")
  if (missing(Glist) || is.null(Glist)) {
    stop("Glist is mandatory for mtblr_block_eigen().", call. = FALSE)
  }
  if (missing(block_start) || is.null(block_start)) {
    stop("block_start is mandatory for mtblr_block_eigen().", call. = FALSE)
  }
  operator_sharing <- match.arg(operator_sharing)
  marker_policy <- match.arg(marker_policy)
  chain <- .blr_chain_controls(
    nit, nburn, nthin, seed, nchains, ncores, chain_seeds, keep_chains)
  conv <- .blr_convergence_controls(
    convergence, convergence_control, chain$nchains)
  if (!identical(summary_reference, "same_bed_by_construction")) {
    .mtblr_block_eigen_reference_error(
      "summary_reference must be 'same_bed_by_construction'.")
  }
  if (!identical(sample_overlap, "not_modeled")) {
    stop("sample_overlap must be exactly 'not_modeled'.", call. = FALSE)
  }
  st <- .mtblr_normalize_stats(stats)
  if (!is.null(n) &&
      !identical(as.integer(rep(n, length.out = st$nt)), st$n)) {
    stop("n conflicts with stats sample sizes.", call. = FALSE)
  }
  glists <- if (is.list(Glist) && !is.null(Glist$bedfiles)) list(Glist) else Glist
  if (!is.list(glists) || !(length(glists) %in% c(1L, st$nt))) {
    stop("Glist must be one genotype list or one genotype list per trait.",
         call. = FALSE)
  }
  glist_input_count <- length(glists)
  glists <- if (length(glists) == 1L) rep(glists, st$nt) else glists
  references <- Map(.mtblr_block_eigen_reference, glists,
                    st$genotype_provenance)
  ld <- list(descriptors = lapply(references, function(x) list(
    marker_ids = x$marker_ids, marker_metadata = x$marker_metadata,
    scale = x$scale, source = x$source, by_construction = TRUE)))
  aligned <- .mtblr_align(st, ld, marker_policy)
  st <- aligned$stats
  blocks <- .mtblr_block_eigen_specs(block_start, st$m, st$nt)
  filters <- .mtblr_block_eigen_vector(
    eigen_filter, st$nt, "eigen_filter",
    c("hard_truncate", "ridge_fixed", "ridge_lw"))
  taus <- .mtblr_block_eigen_vector(eigen_tau, st$nt, "eigen_tau")
  etas <- .mtblr_block_eigen_vector(eigen_eta, st$nt, "eigen_eta")
  common_provenance <- .mtblr_block_eigen_same_provenance(references)
  unambiguously_shared <- glist_input_count == 1L &&
    blocks$input_count == 1L && filters$input_count == 1L &&
    taus$input_count == 1L && etas$input_count == 1L && common_provenance
  if (operator_sharing == "shared" && !unambiguously_shared) {
    stop("operator_sharing = 'shared' requires one common Glist, block, filter, tau, eta, and identical genotype provenance.",
         call. = FALSE)
  }
  resolved_policy <- if (operator_sharing == "auto") {
    if (unambiguously_shared) "shared" else "trait_specific"
  } else operator_sharing
  boundary_shared <- all(vapply(blocks$values[-1L], identical, logical(1),
                                blocks$values[[1L]]))
  sharing_mode <- if (resolved_policy == "shared") {
    "fully_shared_operator"
  } else if (boundary_shared) {
    "trait_specific_shared_boundaries"
  } else {
    "trait_specific_operator"
  }
  descriptor_for_trait <- lapply(seq_len(st$nt), function(t) list(
    bed_files = references[[t]]$bed_files,
    n_bed = references[[t]]$n_bed,
    cls = references[[t]]$cls,
    rows = references[[t]]$rows,
    af = references[[t]]$af,
    block_start = blocks$values[[t]] - 1L,
    eigen_filter = filters$values[t],
    eigen_tau = taus$values[t],
    eigen_eta = etas$values[t]))
  operator_descriptors <- if (resolved_policy == "shared") {
    descriptor_for_trait[1L]
  } else descriptor_for_trait
  owner_count <- length(operator_descriptors)
  trait_owner <- if (owner_count == 1L) rep(1L, st$nt) else seq_len(st$nt)
  pattern_spec <- .mtblr_models(models, pimodels, pi, st$nt)
  maf_info <- .mtblr_resolve_selection_maf(
    selection_maf, !is.null(selection_s) || isTRUE(estimate_selection_s),
    st$m, summary_marker_metadata = aligned$marker_metadata,
    reference_marker_metadata = references[[1L]]$marker_metadata,
    allow_reference_maf_for_selection_s =
      allow_reference_maf_for_selection_s)
  mixture <- .mtblr_bayesr_spec(
    semantics$prior_kernel, pattern_spec, maf_info$values,
    st$m, mixture_var, joint_pi, joint_pi_prior, component, selection_s,
    estimate_selection_s, selection_s_init, selection_s_prior,
    selection_s_proposal_sd)
  mod <- mixture$patterns
  set_spec <- .mtblr_sets(sets, st$m)
  h2 <- rep(h2, length.out = st$nt)
  if (any(!is.finite(h2)) || any(h2 <= 0 | h2 >= 1)) {
    stop("h2 must be in (0, 1).", call. = FALSE)
  }
  vy <- st$yy / (st$n - 1)
  vg0 <- diag(vy * h2, st$nt)
  ve0 <- diag(vy * (1 - h2), st$nt)
  vb0 <- diag((vy * h2) /
    (st$m * max(1e-12, 1 - mod$probabilities[1L])), st$nt)
  vg <- .mtblr_cov(vg, vg0, "vg", st$nt)
  vb <- .mtblr_cov(vb, vb0, "vb", st$nt)
  ve <- .mtblr_cov(ve, ve0, "ve", st$nt, TRUE)
  ssb_prior <- .mtblr_cov(
    ssb_prior, ((nub - 2) / nub) * vg /
      (st$m * max(1e-12, 1 - mod$probabilities[1L])),
    "ssb_prior", st$nt)
  sse_prior <- .mtblr_cov(
    sse_prior, ((nue - 2) / nue) * ve,
    "sse_prior", st$nt, TRUE)
  initialization <- .mtblr_bayesr_initialization(
    beta,b,state,mixture$component_init,pattern_spec$matrix,st$m,st$nt,
    semantics$prior_kernel)
  b <- initialization$b
  trait_metadata <- if (is.null(trait_metadata)) {
    data.frame(trait_id = st$trait_names)
  } else as.data.frame(trait_metadata, stringsAsFactors = FALSE)
  if (is.null(trait_metadata$trait_id)) {
    trait_metadata$trait_id <- st$trait_names
  }
  if (nrow(trait_metadata) != st$nt ||
      !identical(as.character(trait_metadata$trait_id), st$trait_names) ||
      anyDuplicated(trait_metadata$trait_id)) {
    stop("trait_metadata trait_id must uniquely match trait order.",
         call. = FALSE)
  }
  for (x in c("study_id", "ancestry", "population", "ld_reference")) {
    if (is.null(trait_metadata[[x]])) trait_metadata[[x]] <- NA_character_
  }
  trait_metadata$sample_size <- st$n
  trait_metadata$operator_sharing_mode <- sharing_mode
  memory <- .blr_memory_estimate(
    "mtblr", "block_eigen", st$m, st$nt, chain$nchains, chain$ncores,
    chain$nit, chain$nit + chain$nburn, chain$keep_chains,
    convergence_quantities = if (conv$compute || conv$keep_traces) 5L * st$nt else 0L,
    keep_traces = conv$keep_traces,
    operator_bytes = 8 * st$m * st$m * owner_count)
  memory <- .mtblr_bayesr_memory(
    memory,method,st$m,chain$nchains,chain$ncores,chain$nit+chain$nburn,
    nrow(mod$matrix),mixture$component_count)
  .blr_memory_warning(memory, memory_warning_gb, conv$mode,
                      conv$compute || conv$keep_traces, conv$keep_traces)
  native_execution <- mtblr_block_eigen_chains_raw_internal(
    st$wy, st$yy, lapply(seq_len(st$nt), function(t) b[, t]),
    operator_descriptors, set_spec$native, vb, ve,
    lapply(seq_len(st$nt), function(i) ssb_prior[, i]),
    lapply(seq_len(st$nt), function(i) sse_prior[, i]),
    mod$native, mod$probabilities, nub, nue, updateB, updateE, updatePi,
    st$n, chain$nit, chain$nburn, chain$nthin,
    chain$seed, mixture$method_code, chain$nchains, chain$ncores,
    chain$chain_seeds_native,mixture$joint_component,
    mixture$joint_multiplier,mixture$joint_names,mixture$component_count,
    mixture$marker_scale,initialization$component,mixture$pi_prior,
    lapply(seq_len(st$nt), function(t) initialization$beta[,t]),
    lapply(seq_len(st$nt), function(t) initialization$state[,t]))
  execution <- .mtblr_summary_multichain(
    native_execution, chain, conv, st$trait_names, method, "block_eigen",
    updateB, updateE, mixture$model_parameters)
  raw <- execution$raw
  raw$diagnostics$block_eigen$sharing_mode <- sharing_mode
  raw$model$names <- mod$names
  raw$pi$names <- mod$names
  wy_status <- c(
    hard_truncate = "projected_hard_truncate",
    ridge_fixed = "unchanged_ridge_fixed",
    ridge_lw = "unchanged_ridge_lw")[filters$values]
  aligned$report$bed_provenance_status <- "exact"
  aligned$report$row_provenance_status <- "exact"
  aligned$report$column_provenance_status <- "exact"
  aligned$report$allele_frequency_provenance_status <- "exact"
  aligned$report$wy_transformation_status <- unname(wy_status)
  raw$data <- list(
    marker_metadata = aligned$marker_metadata,
    trait_metadata = trait_metadata,
    sample_size = st$n,
    reference_sample_size = vapply(references, `[[`, integer(1), "n_bed"),
    selected_row_count = vapply(references, `[[`, integer(1),
                                "selected_row_count"),
    operator_sharing_mode = sharing_mode,
    data_level = "summary_statistics",
    selection_maf_source = maf_info$selection_maf_source,
    selection_maf_population = maf_info$selection_maf_population,
    selection_maf_alignment_status = maf_info$selection_maf_alignment_status,
    selection_maf_fallback_used = maf_info$selection_maf_fallback_used,
    trait_owner = trait_owner,
    scale = "standardized_genotype",
    summary_reference = summary_reference,
    bed_reference_source = lapply(references, `[[`, "bed_files"),
    block_start = blocks$values,
    eigen_filter = filters$values,
    eigen_tau = taus$values,
    eigen_eta = etas$values,
    mu_floor = 0.01,
    summary_ww_policy =
      "validated_by_construction_not_used_as_runtime_diagonal")
  raw$alignment <- list(
    marker_policy = marker_policy,
    intersection_policy = "error",
    per_trait = aligned$report,
    orientation_status = aligned$report$allele_status,
    bed_provenance_status = rep("exact", st$nt),
    row_provenance_status = rep("exact", st$nt),
    column_provenance_status = rep("exact", st$nt),
    allele_frequency_provenance_status = rep("exact", st$nt),
    wy_transformation_status = unname(wy_status))
  input <- list(
    method = method, model = method,
    backend = paste0("mt_block_eigen_", semantics$prior_kernel),
    prior_kernel = semantics$prior_kernel,
    data_level = "summary_statistics",
    effect_scale_policy = if (!is.null(selection_s))
      if (semantics$prior_kernel == "bayesr") "component_maf_s" else "maf_s"
      else if (semantics$prior_kernel == "bayesr") "component" else "unit",
    model_semantics_version = 2L,
    model_semantics = "s_prefix_means_summary_statistics",
    summary_reference = summary_reference,
    sample_overlap = "not_modeled",
    phenotype_crossproduct_policy = "marginal_yy_only",
    residual_covariance_policy = "diagonal",
    summary_ww_policy =
      "validated_by_construction_not_used_as_runtime_diagonal",
    m = st$m, nt = st$nt, n = st$n, h2 = h2, nub = nub, nue = nue,
    nit = chain$nit, nburn = chain$nburn, nthin = chain$nthin,
    seed = chain$seed, nchains = chain$nchains, ncores = chain$ncores,
    chain_seeds_requested = chain$chain_seeds_requested,
    chain_seeds_resolved = execution$seeds,
    keep_chains = chain$keep_chains, convergence = conv$mode,
    convergence_control = conv[c("warn", "rhat_threshold",
      "ess_per_chain_threshold", "mcse_mean_over_sd_threshold",
      "keep_traces")], memory_warning_gb = memory_warning_gb,
    updateB = updateB, updateE = updateE, updatePi = updatePi,
    models = mod$matrix, model_names = mod$names,
    pimodels = mod$probabilities, mixture_var = mixture$mixture_var,
    joint_pi_prior = mixture$pi_prior, selection_s = mixture$selection_s,
    estimate_selection_s = estimate_selection_s,
    selection_maf_source = maf_info$selection_maf_source,
    selection_maf_population = maf_info$selection_maf_population,
    selection_maf_alignment_status = maf_info$selection_maf_alignment_status,
    selection_maf_fallback_used = maf_info$selection_maf_fallback_used,
    sets = set_spec$public,
    operator_sharing_mode = sharing_mode, trait_owner = trait_owner,
    block_start = blocks$values, eigen_filter = filters$values,
    eigen_tau = taus$values, eigen_eta = etas$values, mu_floor = 0.01,
    reference_sample_size = raw$data$reference_sample_size,
    selected_row_count = raw$data$selected_row_count,
    marker_policy = marker_policy, marker_intersection_policy = "error",
    scale = "standardized_genotype", alignment = raw$alignment,
    trait_metadata = trait_metadata)
  if (isTRUE(verbose)) {
    print(input[c("backend", "m", "nt", "operator_sharing_mode", "seed",
                  "nchains", "ncores")])
  }
  fit <- .as_mtblr_fit(raw, aligned$marker_ids, st$trait_names,
                       aligned$marker_metadata, trait_metadata,
                       raw$alignment, input)
  diagnostics <- raw$diagnostics$block_eigen
  diagnostics$owners <- lapply(diagnostics$owners, function(owner) {
    owner$blocks$start_1based <- owner$blocks$start + 1L
    owner
  })
  fit$block_diagnostics <- diagnostics
  fit$operator_metadata <- list(
    summary_reference = summary_reference,
    sharing_mode = sharing_mode,
    owner_count = owner_count,
    trait_owner = trait_owner,
    block_boundary_sharing = boundary_shared,
    numerical_sharing = resolved_policy == "shared",
    filter_sharing = length(unique(filters$values)) == 1L,
    reference_sharing = common_provenance,
    bed_reference_source = raw$data$bed_reference_source,
    reference_sample_size = raw$data$reference_sample_size,
    selected_row_count = raw$data$selected_row_count)
  fit$chains <- execution$chains
  fit$convergence <- execution$convergence
  fit["convergence_traces"] <- list(execution$convergence_traces)
  fit$input$memory_estimate <- memory
  fit <- .mtblr_bayesr_format_fit(fit, mixture$model_parameters)
  messages <- .blr_convergence_warning_messages(
    fit$convergence, if (conv$mode == "core") "core" else "auto",
    "mtblr", "block_eigen")
  if (isTRUE(conv$warn) && conv$mode != "none" &&
      !(conv$mode == "auto" && chain$nchains == 1L) && length(messages)) {
    warning(messages[[1L]], call. = FALSE)
  }
  .blr_finalize_fit(
    fit, "mtblr", method, "block_eigen", data = raw$data,
    diagnostics = c(raw$diagnostics, list(
      block_eigen_public = diagnostics,
      operator_metadata = fit$operator_metadata)),
    memory_estimate = memory)
}
