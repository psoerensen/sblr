.mtblr_marker_metadata <- function(marker_ids, metadata = NULL) {
  marker_ids <- as.character(marker_ids)
  if (!length(marker_ids) || anyNA(marker_ids) || any(!nzchar(marker_ids)) || anyDuplicated(marker_ids)) {
    stop("Marker IDs must be unique, nonempty, and non-missing.", call. = FALSE)
  }
  if (is.null(metadata)) metadata <- data.frame(marker_id = marker_ids)
  metadata <- as.data.frame(metadata, stringsAsFactors = FALSE)
  if (is.null(metadata$marker_id)) metadata$marker_id <- marker_ids
  if (nrow(metadata) != length(marker_ids) || !identical(as.character(metadata$marker_id), marker_ids)) {
    stop("marker_metadata must have one row per marker in marker order.", call. = FALSE)
  }
  metadata
}

.mtblr_stats_provenance <- function(x, marker_ids, marker_metadata, source,
                                    scale, sample_size) {
  fields <- c("bed_files", "n_bed", "cls", "rows", "af")
  out <- setNames(lapply(fields, function(field) x[[field]]), fields)
  out$marker_ids <- marker_ids
  out$marker_metadata <- marker_metadata
  out$source <- unname(source)
  out$scale <- unname(scale)
  out$sample_size <- sample_size
  out
}

.mtblr_normalize_stats <- function(stats) {
  is_multi <- is.list(stats) && all(c("wy", "ww", "yy") %in% names(stats))
  if (is_multi) {
    wy <- stats$wy; ww <- stats$ww
    if (!is.list(wy) || !is.list(ww) || length(wy) != length(ww) || !length(wy)) {
      stop("A multi-trait stats object must contain wy and ww lists of equal positive length.", call. = FALSE)
    }
    nt <- length(wy)
    trait_names <- stats$trait_names %||% names(wy) %||% names(ww) %||% names(stats$yy)
    if (is.null(trait_names)) trait_names <- paste0("T", seq_len(nt))
    marker_by_trait <- lapply(seq_len(nt), function(t) {
      stats$marker_names %||% names(wy[[t]]) %||% names(ww[[t]])
    })
    marker_metadata <- replicate(nt, stats$marker_metadata, simplify = FALSE)
    source <- rep(stats$source %||% "external", nt)
    scale <- rep(stats$scale %||% NA_character_, nt)
    n <- stats$n
    yy <- stats$yy
    provenance_source <- replicate(nt, stats, simplify = FALSE)
  } else if (is.list(stats) && length(stats) && !is.null(names(stats)) && all(nzchar(names(stats)))) {
    entries <- stats; nt <- length(entries); trait_names <- names(entries)
    valid <- vapply(entries, function(x) is.list(x) && all(c("wy", "ww", "yy", "n") %in% names(x)), logical(1))
    if (!all(valid)) stop("A named list of single-trait stats must contain wy, ww, yy, and n in every entry.", call. = FALSE)
    one <- function(x, field) {
      value <- x[[field]]
      if (is.list(value)) {
        if (length(value) != 1L) stop("Single-trait wy and ww lists must have length one.", call. = FALSE)
        value <- value[[1L]]
      }
      as.numeric(value)
    }
    wy <- lapply(entries, one, field = "wy")
    ww <- lapply(entries, one, field = "ww")
    marker_by_trait <- lapply(seq_len(nt), function(t) entries[[t]]$marker_names %||%
      names(if (is.list(entries[[t]]$wy)) entries[[t]]$wy[[1L]] else entries[[t]]$wy))
    marker_metadata <- lapply(entries, `[[`, "marker_metadata")
    source <- vapply(entries, function(x) x$source %||% "external", character(1))
    scale <- vapply(entries, function(x) x$scale %||% NA_character_, character(1))
    yy <- vapply(entries, function(x) as.numeric(x$yy)[1L], numeric(1))
    n <- vapply(entries, function(x) as.integer(x$n)[1L], integer(1))
    provenance_source <- entries
  } else {
    stop("stats must be one multi-trait stats object or a named list of single-trait stats objects.", call. = FALSE)
  }
  if (length(trait_names) != nt || anyNA(trait_names) || any(!nzchar(trait_names)) || anyDuplicated(trait_names)) {
    stop("Trait names must be unique, nonempty, and non-missing.", call. = FALSE)
  }
  if (length(n) == 1L) n <- rep(n, nt)
  if (length(n) != nt || any(!is.finite(n)) || any(n <= 1) || any(n != as.integer(n))) stop("n must contain positive integer sample sizes greater than one.", call. = FALSE)
  if (is.matrix(yy)) {
    if (!all(dim(yy) == nt)) stop("yy matrix must be nt by nt.", call. = FALSE)
    if (any(yy[row(yy) != col(yy)] != 0)) stop("Nonzero off-diagonal yy is unsupported: sample-overlap and cross-trait SSY likelihoods are not modeled.", call. = FALSE)
    yy <- diag(yy)
  }
  if (length(yy) != nt || any(!is.finite(yy))) stop("yy must be a finite vector of length nt.", call. = FALSE)
  m <- length(wy[[1L]])
  if (!m || any(lengths(wy) != m) || any(lengths(ww) != m)) stop("Every wy and ww vector must share one positive marker count.", call. = FALSE)
  if (any(!is.finite(unlist(wy))) || any(!is.finite(unlist(ww))) || any(unlist(ww) <= 0)) stop("wy must be finite and ww must be finite and positive.", call. = FALSE)
  marker_by_trait <- lapply(marker_by_trait, function(x) {
    if (is.null(x)) stop("Every stats trait requires marker IDs.", call. = FALSE)
    .mtblr_marker_metadata(x)$marker_id
  })
  marker_metadata <- Map(.mtblr_marker_metadata, marker_by_trait, marker_metadata)
  genotype_provenance <- lapply(seq_len(nt), function(t)
    .mtblr_stats_provenance(provenance_source[[t]], marker_by_trait[[t]],
      marker_metadata[[t]], source[t], scale[t], as.integer(n[t])))
  names(wy) <- names(ww) <- trait_names
  list(wy = lapply(wy, as.numeric), ww = lapply(ww, as.numeric), yy = as.numeric(yy),
       n = as.integer(n), m = m, nt = nt, trait_names = trait_names,
       marker_ids = marker_by_trait, marker_metadata = marker_metadata,
       source = source, scale = scale,
       genotype_provenance = genotype_provenance)
}

.mtblr_glist_descriptor <- function(Glist) {
  sparse <- Glist$sparseLD
  if (is.null(sparse$prefix) || !nzchar(sparse$prefix)) stop("Glist$sparseLD$prefix is required.", call. = FALSE)
  ids <- sparse$marker_names %||% unlist(Glist$rsidsLD, use.names = FALSE)
  if (is.null(ids) && !is.null(sparse$chr) && !is.null(sparse$cls)) {
    ids <- unlist(Map(function(cc, cl) Glist$rsids[[cc]][cl], sparse$chr, sparse$cls), use.names = FALSE)
  }
  ids <- .mtblr_marker_metadata(ids)$marker_id
  metadata <- sparse$marker_metadata %||% data.frame(marker_id = ids)
  list(prefix = sparse$prefix, marker_ids = ids,
       marker_metadata = .mtblr_marker_metadata(ids, metadata),
       scale = sparse$scale %||% "standardized_genotype",
       source = sparse$source %||% "make_sparse_ld",
       reference_id = sparse$reference_id %||% NA_character_,
       ancestry = sparse$ancestry %||% NA_character_, population = sparse$population %||% NA_character_,
       by_construction = identical(sparse$source %||% "make_sparse_ld", "make_sparse_ld"))
}

.mtblr_resolve_ld <- function(Glist, ld_prefix, ld_metadata, nt) {
  glists <- if (is.null(Glist)) NULL else if (is.list(Glist) && !is.null(Glist$sparseLD)) list(Glist) else Glist
  descriptors <- if (!is.null(glists)) lapply(glists, .mtblr_glist_descriptor) else ld_metadata
  if (!is.null(descriptors) && !is.list(descriptors[[1L]])) descriptors <- list(descriptors)
  prefixes <- ld_prefix
  if (is.null(prefixes) && !is.null(descriptors)) prefixes <- vapply(descriptors, `[[`, character(1), "prefix")
  if (is.null(prefixes) || !(length(prefixes) %in% c(1L, nt))) stop("ld_prefix must have length one or nt, or resolve from Glist.", call. = FALSE)
  prefixes <- as.character(prefixes)
  if (anyNA(prefixes) || any(!nzchar(prefixes))) stop("LD prefixes must be nonempty.", call. = FALSE)
  if (is.null(descriptors)) stop("Explicit ld_prefix without Glist requires ld_metadata.", call. = FALSE)
  if (!(length(descriptors) %in% c(1L, nt))) stop("LD metadata must contain one descriptor or one per trait.", call. = FALSE)
  descriptors <- if (length(descriptors) == 1L) rep(descriptors, nt) else descriptors
  expected <- if (length(prefixes) == 1L) rep(prefixes, nt) else prefixes
  supplied <- vapply(descriptors, function(x) as.character(x$prefix %||% ""), character(1))
  if (any(nzchar(supplied) & supplied != expected)) stop("Glist/ld_metadata prefixes conflict with ld_prefix.", call. = FALSE)
  for (i in seq_len(nt)) {
    descriptors[[i]]$prefix <- expected[i]
    descriptors[[i]]$marker_ids <- .mtblr_marker_metadata(descriptors[[i]]$marker_ids)$marker_id
    descriptors[[i]]$marker_metadata <- .mtblr_marker_metadata(descriptors[[i]]$marker_ids, descriptors[[i]]$marker_metadata)
  }
  list(prefixes = prefixes, descriptors = descriptors)
}

.mtblr_normalize_scale <- function(x) {
  x <- unname(x)
  if (isTRUE(x)) return("standardized_genotype")
  if (identical(x, "standardized_genotype")) return(x)
  NA_character_
}

.mtblr_align <- function(stats, ld, marker_policy) {
  canonical <- ld$descriptors[[1L]]$marker_ids
  if (any(vapply(ld$descriptors, function(x) !identical(x$marker_ids, canonical), logical(1)))) stop("All LD resources must already use identical marker IDs and order.", call. = FALSE)
  report <- vector("list", stats$nt)
  for (t in seq_len(stats$nt)) {
    ids <- stats$marker_ids[[t]]
    if (marker_policy == "strict" && !identical(ids, canonical)) stop("strict marker policy requires identical stats and LD marker order.", call. = FALSE)
    if (marker_policy == "reorder_stats" && !identical(ids, canonical)) {
      if (length(ids) != length(canonical) || !setequal(ids, canonical)) stop("Stats and LD marker sets differ; automatic intersection is not supported.", call. = FALSE)
      idx <- match(canonical, ids)
      stats$wy[[t]] <- stats$wy[[t]][idx]; stats$ww[[t]] <- stats$ww[[t]][idx]
      stats$marker_metadata[[t]] <- stats$marker_metadata[[t]][idx, , drop = FALSE]
      stats$marker_ids[[t]] <- canonical
    }
    sm <- stats$marker_metadata[[t]]; lm <- ld$descriptors[[t]]$marker_metadata
    explicit <- all(c("effect_allele", "other_allele") %in% names(sm)) && all(c("effect_allele", "other_allele") %in% names(lm))
    allele_status <- "unresolved"
    if (explicit) {
      se <- toupper(sm$effect_allele); so <- toupper(sm$other_allele)
      le <- toupper(lm$effect_allele); lo <- toupper(lm$other_allele)
      if (any(se == lo & so == le)) stop("Swapped effect/other alleles are not accepted; harmonize upstream.", call. = FALSE)
      comp <- chartr("ACGT", "TGCA", se) == le & chartr("ACGT", "TGCA", so) == lo
      if (any(comp & !(se == le & so == lo))) stop("Strand-complement-only allele matches are not accepted.", call. = FALSE)
      if (any(se != le | so != lo)) stop("Stats and LD allele orientation differs.", call. = FALSE)
      allele_status <- "exact"
    } else if (identical(unname(stats$source[t]), "make_summary_stats") && isTRUE(ld$descriptors[[t]]$by_construction)) {
      provenance_fields <- c("chromosome_or_file", "bed_column")
      if (!all(provenance_fields %in% names(sm)) || !all(provenance_fields %in% names(lm)) ||
          any(vapply(provenance_fields, function(field)
            !identical(sm[[field]], lm[[field]]), logical(1)))) {
        stop("BED by-construction orientation requires identical marker IDs, chromosome/file, and BED columns.", call. = FALSE)
      }
      allele_status <- "by_construction_same_glist"
    } else stop("External stats/LD require explicit effect_allele and other_allele metadata.", call. = FALSE)
    for (metadata in list(sm, lm)) {
      if ("allele_frequency" %in% names(metadata) && any(!is.finite(metadata$allele_frequency))) {
        stop("Allele frequencies must be finite when supplied.", call. = FALSE)
      }
    }
    ss <- .mtblr_normalize_scale(stats$scale[t]); ls <- .mtblr_normalize_scale(ld$descriptors[[t]]$scale)
    if (is.na(ss) || is.na(ls) || ss != ls) stop("Stats and LD must use scale = 'standardized_genotype'.", call. = FALSE)
    report[[t]] <- data.frame(trait_id = stats$trait_names[t], n_markers = stats$m,
      marker_order_status = if (identical(ids, canonical)) "identical" else "stats_reordered",
      allele_status = allele_status, scale_status = "standardized_genotype",
      stats_source = stats$source[t], ld_source = ld$descriptors[[t]]$source %||% "external")
  }
  list(stats = stats, marker_ids = canonical, marker_metadata = stats$marker_metadata[[1L]], report = do.call(rbind, report))
}

.mtblr_models <- function(models, pimodels, pi, nt) {
  if (is.null(models)) models <- as.matrix(expand.grid(rep(list(0:1), nt)))
  if (is.character(models)) {
    if (!identical(models, "restrictive")) stop("models character value must be 'restrictive'.", call. = FALSE)
    models <- rbind(rep(0L, nt), rep(1L, nt))
  }
  models <- as.matrix(models)
  if (ncol(models) != nt || any(!models %in% 0:1) || anyDuplicated(as.data.frame(models))) stop("models must contain unique binary patterns of length nt.", call. = FALSE)
  if (!any(rowSums(models) == 0L)) stop("models must include the null pattern.", call. = FALSE)
  if (nrow(models) > 4096L) stop("Automatic model-pattern expansion is impractically large; supply models explicitly.", call. = FALSE)
  names <- apply(models, 1L, paste, collapse = "_")
  probs <- pimodels
  if (is.null(probs)) probs <- if (length(pi) == nrow(models)) pi else c(1 - pi[1L], rep(pi[1L] / (nrow(models) - 1L), nrow(models) - 1L))
  if (length(probs) != nrow(models) || any(!is.finite(probs)) || any(probs < 0) || sum(probs) <= 0) stop("pimodels must be finite, nonnegative, and match models.", call. = FALSE)
  list(matrix = models, native = lapply(seq_len(nrow(models)), function(i) as.integer(models[i, ])), probabilities = as.numeric(probs / sum(probs)), names = names)
}

.mtblr_sets <- function(sets, m) {
  if (is.null(sets)) sets <- list(seq_len(m))
  if (!is.list(sets) || !length(sets) || any(lengths(sets) == 0L)) stop("sets must be a nonempty list of nonempty marker-index vectors.", call. = FALSE)
  sets <- lapply(sets, function(x) {
    if (!is.numeric(x) || any(!is.finite(x)) || any(x != as.integer(x)) || any(x < 1L | x > m) || anyDuplicated(x)) stop("sets must contain unique integer indices in [1, m].", call. = FALSE)
    as.integer(x)
  })
  flat <- unlist(sets)
  if (anyDuplicated(flat) || !identical(sort(flat), seq_len(m))) stop("sets must form a disjoint complete partition of 1:m.", call. = FALSE)
  list(public = sets, native = lapply(sets, function(x) x - 1L))
}

.mtblr_cov <- function(x, default, name, nt, diagonal = FALSE) {
  if (is.null(x)) x <- default
  x <- as.matrix(x)
  if (!all(dim(x) == nt) || any(!is.finite(x)) ||
      !isTRUE(all.equal(unname(x), unname(t(x)), tolerance = 0))) {
    stop(name, " must be a finite symmetric nt by nt matrix.", call. = FALSE)
  }
  if (diagonal && any(x[row(x) != col(x)] != 0)) stop(name, " must be diagonal for the public MT CSR route.", call. = FALSE)
  if (any(diag(x) <= 0) || min(eigen(x, symmetric = TRUE, only.values = TRUE)$values) <= 0) stop(name, " must be positive definite.", call. = FALSE)
  x
}

.is_mtblr_raw <- function(raw) is.list(raw) && identical(raw$schema$class, "mtblr_raw") && identical(as.integer(raw$schema$version), 1L)
.validate_mtblr_raw <- function(raw) {
  if (!.is_mtblr_raw(raw)) stop("Expected mtblr_raw schema version 1.", call. = FALSE)
  required <- c("meta", "marker", "trace", "variance", "pi", "model", "diagnostics", "data", "alignment")
  if (!all(required %in% names(raw))) stop("mtblr_raw is missing required namespaces.", call. = FALSE)
  m <- raw$meta$m; nt <- raw$meta$nt; ntr <- raw$meta$n_trace; nm <- raw$meta$nmodels
  for (x in c("bm", "dm", "wy", "r", "b", "state")) if (!identical(dim(raw$marker[[x]]), c(m, nt))) stop("Invalid mtblr_raw marker dimensions.", call. = FALSE)
  for (x in c("vbs", "vgs", "ves", "vle", "vld")) if (!identical(dim(raw$trace[[x]]), c(ntr, nt))) stop("Invalid mtblr_raw trace dimensions.", call. = FALSE)
  for (x in c("covb", "covg", "cove", "vb", "vg", "ve")) if (!identical(dim(raw$variance[[x]]), c(nt, nt))) stop("Invalid mtblr_raw covariance dimensions.", call. = FALSE)
  is_bayesrc <- raw$meta$model %in% c("bayesrc", "sbayesrc")
  valid_pi <- if (is_bayesrc) {
    is.null(raw$pi$final) && is.null(raw$pi$mean) && is.null(raw$pi$trace)
  } else length(raw$pi$final) == nm && length(raw$pi$mean) == nm
  if (!valid_pi || !identical(dim(raw$model$patterns), c(nm, nt)))
    stop("Invalid mtblr_raw probability/model dimensions.", call. = FALSE)
  has_components <- !is.null(raw$marker$component_final) ||
    !is.null(raw$marker$component_probabilities)
  if (raw$meta$model %in% c("bayesr", "sbayesr", "bayesrc", "sbayesrc") && !has_components)
    stop("BayesR mtblr_raw objects require component fields.", call. = FALSE)
  if (raw$meta$model %in% c("bayesc", "sbayesc") && has_components)
    stop("BayesC mtblr_raw objects must not contain component fields.", call. = FALSE)
  if (!is.null(raw$meta$model_semantics_version) &&
      (!identical(raw$meta$model_semantics_version, 2L) ||
       !identical(raw$meta$model_semantics,
                  "s_prefix_means_summary_statistics"))) {
    stop("Invalid model-semantics metadata in mtblr_raw.", call. = FALSE)
  }
  if (has_components) {
    if (is.null(raw$marker$component_final) ||
        length(raw$marker$component_final) != m ||
        !identical(nrow(raw$marker$component_probabilities), m) ||
        any(!is.finite(raw$marker$component_probabilities)) ||
        any(abs(rowSums(raw$marker$component_probabilities) - 1) > 1e-10) ||
        (!is_bayesrc && is.null(raw$pi$trace)) || is.null(raw$model$mixture) ||
        length(raw$model$mixture$joint_state_names) != nm ||
        length(raw$model$mixture$joint_component_index) != nm ||
        (!is_bayesrc && !identical(dim(raw$pi$trace), c(ntr, nm)))) {
      stop("Invalid mtblr_raw BayesR component fields.", call. = FALSE)
    }
  }
  if (is_bayesrc) {
    ann <- raw$annotations
    component_count <- ncol(raw$marker$component_probabilities)
    if (!is.list(ann) || !identical(ann$policy, "annotation_probit_stick") ||
        !identical(dim(ann$annotation_coefficients_final),
                   c(length(ann$metadata$processed_annotation_names),
                     component_count - 1L)) ||
        !identical(dim(ann$annotation_coefficients_mean),
                   dim(ann$annotation_coefficients_final)) ||
        length(ann$annotation_variances_final) != component_count - 1L ||
        length(ann$annotation_variances_mean) != component_count - 1L ||
        !identical(nrow(ann$prior_component_probabilities), m) ||
        !identical(ncol(ann$prior_component_probabilities), component_count) ||
        any(!is.finite(ann$prior_component_probabilities)) ||
        any(abs(rowSums(ann$prior_component_probabilities) - 1) > 1e-10) ||
        any(!is.finite(ann$pattern_pi_final)) ||
        abs(sum(ann$pattern_pi_final) - 1) > 1e-10) {
      stop("Invalid mtblr_raw BayesRC annotation fields.", call. = FALSE)
    }
  } else if (!is.null(raw$annotations)) {
    stop("BayesC/BayesR mtblr_raw objects must not contain BayesRC annotations.",
         call. = FALSE)
  }
  if (any(!raw$model$patterns %in% 0:1) || length(raw$marker$order) != m) stop("Invalid mtblr_raw patterns or marker order.", call. = FALSE)
  if (raw$meta$backend %in% c("mt_block_eigen_bayesc", "mt_block_eigen_bayesr",
                              "mt_block_eigen_bayesrc")) {
    block <- raw$diagnostics$block_eigen
    if (!is.list(block) || length(block$owner_count) != 1L ||
        !is.finite(block$owner_count) || block$owner_count < 1L) {
      stop("Invalid mtblr_raw block-eigen owner count.", call. = FALSE)
    }
    owner_count <- as.integer(block$owner_count)
    if (length(block$trait_owner) != nt ||
        any(!is.finite(block$trait_owner)) ||
        any(block$trait_owner < 1L | block$trait_owner > owner_count) ||
        !is.list(block$owners) || length(block$owners) != owner_count) {
      stop("Invalid mtblr_raw block-eigen owner mapping.", call. = FALSE)
    }
    required_blocks <- c("start", "size", "n_kept", "mu_min", "shrink")
    for (owner in block$owners) {
      blocks <- owner$blocks
      if (!is.data.frame(blocks) || !all(required_blocks %in% names(blocks)) ||
          !nrow(blocks) || any(!is.finite(blocks$start)) ||
          any(!is.finite(blocks$size)) || any(blocks$size <= 0) ||
          any(!is.finite(blocks$n_kept)) || any(blocks$n_kept <= 0) ||
          any(blocks$n_kept > blocks$size) ||
          any(!is.finite(blocks$shrink)) ||
          any(blocks$shrink < 0 | blocks$shrink > 1) ||
          !identical(as.integer(blocks$start), cumsum(c(0L, head(as.integer(blocks$size), -1L)))) ||
          sum(blocks$size) != m) {
        stop("Invalid mtblr_raw block-eigen block diagnostics.", call. = FALSE)
      }
    }
  }
  if (raw$meta$backend %in% c("mt_bed_bayesc", "mt_bed_bayesr",
                              "mt_bed_bayesrc")) {
    bed <- raw$diagnostics$mt_bed
    nchains <- if (is.null(raw$meta$nchains)) 1L else as.integer(raw$meta$nchains)
    expected_updates <- if (isTRUE(raw$diagnostics$cove > 0)) {
      (raw$meta$nit + raw$meta$nburn) * nchains
    } else {
      0
    }
    if (!identical(raw$meta$data_level, "individual") ||
        !is.list(bed) ||
        !bed$residual_covariance %in% c("full", "diagonal") ||
        length(bed$sample_count) != 1L || !is.finite(bed$sample_count) ||
        bed$sample_count <= 1 ||
        !identical(as.integer(bed$marker_count), as.integer(m)) ||
        !identical(as.integer(bed$trait_count), as.integer(nt)) ||
        !identical(bed$owner, "PackedBedMatrix") ||
        !identical(bed$view, "BedPackedGenotypeView") ||
        !identical(bed$genotype_scale, "standardized_genotype") ||
        !identical(bed$marker_workspace, "double") ||
        length(bed$marker_cholesky_jitter_attempts) != 1L ||
        !is.finite(bed$marker_cholesky_jitter_attempts) ||
        bed$marker_cholesky_jitter_attempts < 0 ||
        length(bed$marker_cholesky_max_increment) != 1L ||
        !is.finite(bed$marker_cholesky_max_increment) ||
        bed$marker_cholesky_max_increment < 0 ||
        length(bed$full_e_updates) != 1L ||
        !is.finite(bed$full_e_updates) || bed$full_e_updates < 0 ||
        length(bed$diagonal_e_updates) != 1L ||
        !is.finite(bed$diagonal_e_updates) ||
        bed$diagonal_e_updates < 0 ||
        (bed$residual_covariance == "full" &&
         (bed$full_e_updates != expected_updates ||
          bed$diagonal_e_updates != 0)) ||
        (bed$residual_covariance == "diagonal" &&
         (bed$diagonal_e_updates != expected_updates ||
          bed$full_e_updates != 0)) ||
        !identical(bed$sample_residual_returned, FALSE) ||
        !identical(bed$genetic_values_returned, FALSE) ||
        !identical(bed$cpo, "unsupported") ||
        !identical(bed$le_ld, "trait_diagonal_decomposition")) {
      stop("Invalid mtblr_raw individual-level MT BED diagnostics.",
           call. = FALSE)
    }
    if (!is.null(raw$meta$nchains)) {
      stability <- c("bm_sd", "bm_min", "bm_max", "dm_sd", "dm_min", "dm_max")
      valid_stability <- all(stability %in% names(raw$marker)) &&
        all(vapply(raw$marker[stability], function(x) {
          identical(dim(x), c(m, nt)) && all(is.finite(x))
        }, logical(1))) &&
        all(raw$marker$bm_sd >= 0) && all(raw$marker$dm_sd >= 0) &&
        all(raw$marker$bm_min <= raw$marker$bm_max) &&
        all(raw$marker$dm_min <= raw$marker$dm_max)
      vector_fields <- c(
        "chain_seeds", "chain_seconds",
        "chain_marker_cholesky_jitter_attempts",
        "chain_marker_cholesky_max_increment",
        "chain_full_e_updates", "chain_diagonal_e_updates"
      )
      valid_vectors <- all(vapply(vector_fields, function(x) {
        length(bed[[x]]) == nchains && all(is.finite(bed[[x]])) &&
          all(bed[[x]] >= 0)
      }, logical(1)))
      valid_policies <- identical(as.integer(bed$primary_chain), 1L) &&
        identical(bed$final_state_policy, "primary_chain") &&
        identical(bed$posterior_summary_policy, "pooled_retained_samples") &&
        identical(bed$trace_policy, "iterationwise_chain_mean")
      valid_counts <- isTRUE(all.equal(
        bed$marker_cholesky_jitter_attempts,
        sum(bed$chain_marker_cholesky_jitter_attempts), tolerance = 0)) &&
        isTRUE(all.equal(
          bed$marker_cholesky_max_increment,
          max(bed$chain_marker_cholesky_max_increment), tolerance = 0)) &&
        isTRUE(all.equal(bed$full_e_updates,
                         sum(bed$chain_full_e_updates), tolerance = 0)) &&
        isTRUE(all.equal(bed$diagonal_e_updates,
                         sum(bed$chain_diagonal_e_updates), tolerance = 0)) &&
        (bed$residual_covariance == "full" &&
           all(bed$chain_full_e_updates == expected_updates / nchains) &&
           all(bed$chain_diagonal_e_updates == 0) ||
         bed$residual_covariance == "diagonal" &&
           all(bed$chain_diagonal_e_updates == expected_updates / nchains) &&
           all(bed$chain_full_e_updates == 0))
      valid_runtime <- length(raw$meta$keep_chains) == 1L &&
        is.logical(raw$meta$keep_chains) && !is.na(raw$meta$keep_chains) &&
        nchains > 0L && length(bed$requested_cores) == 1L &&
        is.finite(bed$requested_cores) && bed$requested_cores > 0 &&
        length(bed$used_workers) == 1L && is.finite(bed$used_workers) &&
        bed$used_workers > 0 && bed$used_workers <= nchains &&
        length(bed$openmp_available) == 1L &&
        is.logical(bed$openmp_available) && !is.na(bed$openmp_available) &&
        all(bed$chain_seeds <= 2^32 - 1) &&
        all(is.finite(c(bed$seconds_mean, bed$seconds_max,
                        bed$dispatch_seconds))) &&
        all(c(bed$seconds_mean, bed$seconds_max,
              bed$dispatch_seconds) >= 0)
      valid_chains <- if (!isTRUE(raw$meta$keep_chains)) {
        "chains" %in% names(raw) && is.null(raw$chains)
      } else {
        is.list(raw$chains) && length(raw$chains) == nchains &&
          identical(names(raw$chains), paste0("chain", seq_len(nchains))) &&
          all(vapply(seq_len(nchains), function(i) {
            chain <- raw$chains[[i]]
            component_fields <- if (has_components)
              c("component_final", "component_probabilities") else character()
            pi_fields <- if (has_components) c("final", "mean", "trace") else
              c("final", "mean")
            chain_names <- c("chain", "seed", "marker", "trace", "variance",
                             "pi", "diagnostics",
                             if (is_bayesrc) "model_parameters")
            is.list(chain) && identical(names(chain), c(
              chain_names)) &&
              identical(as.integer(chain$chain), i) &&
              identical(as.numeric(chain$seed), as.numeric(bed$chain_seeds[i])) &&
              identical(names(chain$marker),
                        c("bm", "dm", "b", "state", component_fields)) &&
              identical(names(chain$trace), c("vbs", "vgs", "ves", "vle", "vld")) &&
              identical(names(chain$variance),
                        c("covb", "covg", "cove", "vb", "vg", "ve")) &&
              identical(names(chain$pi), pi_fields) &&
              identical(dim(chain$marker$bm), c(m, nt)) &&
              identical(dim(chain$marker$dm), c(m, nt)) &&
              identical(dim(chain$marker$b), c(m, nt)) &&
              identical(dim(chain$marker$state), c(m, nt)) &&
              (!has_components ||
                 (length(chain$marker$component_final) == m &&
                  identical(nrow(chain$marker$component_probabilities), m) &&
                  (is_bayesrc || identical(dim(chain$pi$trace), c(ntr, nm))))) &&
              identical(dim(chain$trace$vbs), c(ntr, nt)) &&
              identical(dim(chain$trace$vgs), c(ntr, nt)) &&
              identical(dim(chain$trace$ves), c(ntr, nt)) &&
              all(vapply(chain$variance, function(x) {
                identical(dim(x), c(nt, nt))
              }, logical(1))) &&
              (is_bayesrc ||
                 (length(chain$pi$final) == nm && length(chain$pi$mean) == nm))
          }, logical(1)))
      }
      if (!valid_stability || !valid_vectors || !valid_policies ||
          !valid_counts || !valid_runtime || !valid_chains) {
        stop("Invalid mtblr_raw multichain MT BED extension.", call. = FALSE)
      }
    }
  }
  if (!is.null(raw$diagnostics$convergence)) {
    raw$diagnostics$convergence <-
      .blr_validate_convergence_result(raw$diagnostics$convergence)
  }
  raw
}

.as_mtblr_fit <- function(raw, marker_ids, trait_names, marker_metadata, trait_metadata, alignment, input) {
  raw <- .validate_mtblr_raw(raw)
  dimnames_marker <- list(marker_ids, trait_names); dimnames_trace <- list(paste0("Iter", seq_len(raw$meta$n_trace)), trait_names)
  marker <- lapply(raw$marker[c("bm", "dm", "wy", "r", "b", "state")], function(x) { dimnames(x) <- dimnames_marker; x })
  trace <- lapply(raw$trace, function(x) { dimnames(x) <- dimnames_trace; x })
  variance <- lapply(raw$variance, function(x) { dimnames(x) <- list(trait_names, trait_names); x })
  model_names <- raw$model$names
  if (!is.null(raw$pi$final)) names(raw$pi$final) <- model_names
  if (!is.null(raw$pi$mean)) names(raw$pi$mean) <- model_names
  if (!is.null(raw$pi$trace)) colnames(raw$pi$trace) <- model_names
  o <- matrix(raw$marker$order, nrow = length(marker_ids), ncol = length(trait_names), dimnames = dimnames_marker)
  fit <- c(marker[c("bm", "dm", "wy", "r", "b")], list(d = marker$state, o = o), trace, variance,
           list(pi = raw$pi$final, pim = raw$pi$mean, marker_order = raw$marker$order,
                model_patterns = raw$model$patterns, marker_metadata = marker_metadata,
                trait_metadata = trait_metadata, alignment = alignment, input = input,
                raw_schema_version = 1L))
  if (!is.null(raw$marker$component_final)) {
    fit$component_final <- as.integer(raw$marker$component_final)
    names(fit$component_final) <- marker_ids
    fit$component_probabilities <- raw$marker$component_probabilities
    rownames(fit$component_probabilities) <- marker_ids
    fit$pi_trace <- raw$pi$trace
  }
  fit["rb"] <- list(if (sum(diag(fit$covb)) > 0) cov2cor(fit$covb) else NULL)
  fit["rg"] <- list(if (sum(diag(fit$covg)) > 0) cov2cor(fit$covg) else NULL)
  fit["re"] <- list(if (sum(diag(fit$cove)) > 0) cov2cor(fit$cove) else NULL)
  stability <- c("bm_sd", "bm_min", "bm_max", "dm_sd", "dm_min", "dm_max")
  if (all(stability %in% names(raw$marker))) {
    for (field in stability) {
      value <- raw$marker[[field]]
      dimnames(value) <- dimnames_marker
      fit[[field]] <- value
    }
  }
  if (!is.null(raw$meta$nchains)) {
    bed <- raw$diagnostics$mt_bed
    fit$nchains <- as.integer(raw$meta$nchains)
    fit$chain_seeds <- bed$chain_seeds
    fit$chain_diagnostics <- bed[c(
      "requested_cores", "used_workers", "openmp_available",
      "chain_seconds", "seconds_mean", "seconds_max", "dispatch_seconds",
      "primary_chain", "final_state_policy", "posterior_summary_policy",
      "trace_policy", "chain_marker_cholesky_jitter_attempts",
      "chain_marker_cholesky_max_increment", "chain_full_e_updates",
      "chain_diagonal_e_updates"
    )]
    if (is.null(raw$chains)) {
      fit["chains"] <- list(NULL)
    } else {
      formatted_chains <- lapply(raw$chains, function(chain) {
        for (field in intersect(c("bm", "dm", "b", "state"),
                                names(chain$marker)))
          dimnames(chain$marker[[field]]) <- dimnames_marker
        if (!is.null(chain$marker$component_final))
          names(chain$marker$component_final) <- marker_ids
        if (!is.null(chain$marker$component_probabilities))
          rownames(chain$marker$component_probabilities) <- marker_ids
        for (field in names(chain$trace))
          dimnames(chain$trace[[field]]) <- dimnames_trace
        for (field in names(chain$variance))
          dimnames(chain$variance[[field]]) <- list(trait_names, trait_names)
        if (!is.null(chain$pi$final)) names(chain$pi$final) <- model_names
        if (!is.null(chain$pi$mean)) names(chain$pi$mean) <- model_names
        if (!is.null(chain$pi$trace)) colnames(chain$pi$trace) <- model_names
        chain
      })
      fit$chains <- formatted_chains
    }
  }
  if (!is.null(raw$diagnostics$convergence)) {
    fit$convergence <- raw$diagnostics$convergence
  }
  class(fit) <- c("mtblr_fit", "list")
  fit
}

#' Joint multivariate BayesC, BayesR, and SBayesR with sparse LD
#'
#' Fits aligned joint multivariate mixture models from standardized
#' summary statistics and one shared or one-per-trait disk-backed CSR LD
#' resource. Marker identity, order, scale, and available allele orientation
#' metadata are validated before native execution. The function uses marginal
#' `yy` only; sample overlap is not modeled and residual covariance is diagonal.
#'
#' @param stats A multi-trait sufficient-statistics object or named list of
#'   single-trait objects containing `wy`, `ww`, `yy`, `n`, and marker IDs.
#' @param Glist One Glist or a list of one Glist per trait.
#' @param ld_prefix One CSR prefix or one prefix per trait.
#' @param ld_metadata LD descriptors required when no Glist is supplied.
#' @param trait_metadata Optional data frame with trait/study provenance.
#' @param marker_policy Either strict order matching or R-side stats reordering.
#' @param sample_overlap Must be `"not_modeled"`.
#' @param method Exactly `"sbayesc"`, `"sbayesr"`, or `"sbayesrc"`; the `s` prefix denotes
#'   summary-statistics data and does not activate MAF scaling.
#' @param n Optional sample sizes overriding identical values only.
#' @param sets Optional disjoint complete marker partition (1-based).
#' @param beta,b,state Optional BayesR latent effects, effective effects, and
#'   binary trait-pattern states. BayesC accepts `b` only.
#' @param h2 Requested initial expected genetic-variance fraction, scalar or
#'   one value per trait, under resolved joint trait-pattern, component,
#'   annotation, and fixed MAF-S weights.
#' @param pi Initial non-null probability or complete pattern probabilities.
#' @param models,pimodels Pattern matrix and probabilities.
#' @param mixture_var Fixed BayesR component-variance multipliers: one leading
#'   zero followed by unique ascending positive values.
#' @param joint_pi Optional initial probability vector over the deterministic
#'   joint pattern-by-component states.
#' @param joint_pi_prior Optional positive Dirichlet prior over joint states.
#' @param component Optional zero-based component initialization per marker.
#' @param annotations Required marker-by-annotation numeric matrix or data
#'   frame for `"sbayesrc"`. Explicit unique marker IDs are required and are
#'   matched to the final marker order.
#' @param add_intercept Add one intercept when none is supplied.
#' @param standardize_annotations Standardize eligible non-intercept columns.
#' @param center_binary_annotations Center and scale binary annotations when
#'   standardization is enabled.
#' @param alpha_init Optional processed-annotation-by-stick coefficient matrix.
#' @param sigmaSqAlpha_init Optional positive variance initialization per stick.
#' @param intercept_flat Use a flat prior for the first intercept coefficient.
#' @param sigmaSqAlpha_a,sigmaSqAlpha_b Positive annotation-variance prior
#'   hyperparameters.
#' @param pi_floor Probability floor used by probit stick-breaking.
#' @param alpha_update_every Positive iteration interval for coefficient updates.
#' @param updateAlpha Update annotation coefficients and their variances.
#' @param maf_effect_s Optional fixed scalar MAF-S exponent, independent of the
#'   summary-statistics `sbayesr` model name.
#' @param effect_maf Optional allele frequencies aligned to the final marker
#'   order for the independent `maf_effect_s` scale policy.
#' @param allow_reference_maf_for_maf_effect_s Allow explicit fallback to
#'   reference-panel MAF when GWAS-summary MAF is unavailable.
#' @param estimate_maf_effect_s Logical; sampled MT S is currently unsupported.
#' @param maf_effect_s_init,maf_effect_s_prior,maf_effect_s_proposal_sd Reserved
#'   sampled-S controls, rejected while `estimate_maf_effect_s` is unsupported.
#' @param vg,vb,ve Initial covariance matrices.
#' @param ssb_prior,sse_prior Covariance prior scale matrices.
#' @param updateB,updateE,updatePi Update controls.
#' @param nub,nue Prior degrees of freedom.
#' @param nit,nburn,nthin MCMC controls.
#' @param seed Explicit native RNG seed.
#' @param nchains Number of complete joint-MT logical chains.
#' @param ncores Requested logical-chain workers.
#' @param chain_seeds Optional signed integer seed per chain.
#' @param keep_chains Retain compact logical-chain records.
#' @param convergence Convergence mode: `"auto"`, `"none"`, `"core"`, or
#'   `"extended"`. Automatic mode remains core-only.
#' @param convergence_control Optional uniquely named convergence-control list.
#'   Extended mode can request covariance, probability, annotation, or
#'   selection-S groups plus explicitly selected marker traces, subject to a
#'   pre-execution hard memory guard.
#' @param memory_warning_gb Analytical memory-warning threshold in GiB.
#' @param verbose Print resolved execution metadata.
#' @return An object of class `mtblr_fit`.
#' @export
mtblr_csr <- function(stats, Glist = NULL, ld_prefix = NULL, ld_metadata = NULL,
  trait_metadata = NULL, marker_policy = c("strict", "reorder_stats"),
  sample_overlap = "not_modeled", method = "sbayesc", n = NULL, sets = NULL,
  beta = NULL, b = NULL, state = NULL, h2 = 0.5, pi = 0.001,
  models = NULL, pimodels = NULL, mixture_var = NULL, joint_pi = NULL,
  joint_pi_prior = NULL, component = NULL,
  annotations = NULL, add_intercept = TRUE,
  standardize_annotations = TRUE, center_binary_annotations = FALSE,
  alpha_init = NULL, sigmaSqAlpha_init = NULL, intercept_flat = TRUE,
  sigmaSqAlpha_a = 2, sigmaSqAlpha_b = 2, pi_floor = 1e-12,
  alpha_update_every = 1L, updateAlpha = TRUE, maf_effect_s = NULL,
  effect_maf = NULL, allow_reference_maf_for_maf_effect_s = FALSE,
  estimate_maf_effect_s = FALSE, maf_effect_s_init = NULL,
  maf_effect_s_prior = NULL, maf_effect_s_proposal_sd = NULL,
  vg = NULL, vb = NULL, ve = NULL, ssb_prior = NULL, sse_prior = NULL,
  updateB = TRUE, updateE = TRUE, updatePi = TRUE, nub = 4, nue = 4,
  nit = 1000, nburn = 500, nthin = 1, seed = 1,
  nchains = 1L, ncores = 1L, chain_seeds = NULL,
  keep_chains = FALSE, convergence = c("auto", "none", "core", "extended"),
  convergence_control = NULL, memory_warning_gb = 8, verbose = FALSE) {
  marker_policy <- match.arg(marker_policy)
  chain <- .blr_chain_controls(
    nit, nburn, nthin, seed, nchains, ncores, chain_seeds, keep_chains)
  conv <- .blr_convergence_controls(
    convergence, convergence_control, chain$nchains)
  if (!identical(sample_overlap, "not_modeled")) stop("sample_overlap must be exactly 'not_modeled'.", call. = FALSE)
  semantics <- .mtblr_resolve_public_method(method, "csr")
  st <- .mtblr_normalize_stats(stats); if (!is.null(n) && !identical(as.integer(rep(n, length.out = st$nt)), st$n)) stop("n conflicts with stats sample sizes.", call. = FALSE)
  ld <- .mtblr_resolve_ld(Glist, ld_prefix, ld_metadata, st$nt)
  aligned <- .mtblr_align(st, ld, marker_policy); st <- aligned$stats
  shared <- length(ld$prefixes) == 1L
  sharing <- if (shared && all(vapply(st$ww[-1L], identical, logical(1), st$ww[[1L]]))) "fully_shared_operator" else if (shared) "shared_correlation_reference" else "trait_specific_reference"
  native_prefix <- if (sharing == "shared_correlation_reference") rep(ld$prefixes, st$nt) else ld$prefixes
  pattern_spec <- .mtblr_models(models, pimodels, pi, st$nt)
  maf_info <- .mtblr_resolve_effect_maf(
    effect_maf, !is.null(maf_effect_s) || isTRUE(estimate_maf_effect_s),
    st$m, summary_marker_metadata = aligned$marker_metadata,
    reference_marker_metadata = ld$descriptors[[1L]]$marker_metadata,
    allow_reference_maf_for_maf_effect_s =
      allow_reference_maf_for_maf_effect_s)
  mixture <- .mtblr_bayesr_spec(
    semantics$prior_kernel, pattern_spec, maf_info$values,
    st$m, mixture_var, joint_pi, joint_pi_prior, component, maf_effect_s,
    estimate_maf_effect_s, maf_effect_s_init, maf_effect_s_prior,
    maf_effect_s_proposal_sd)
  if (identical(semantics$prior_kernel, "bayesrc")) {
    if (!is.null(joint_pi) || !is.null(joint_pi_prior))
      stop("joint_pi and joint_pi_prior are not BayesRC controls; use pimodels for conditional pattern initialization.", call. = FALSE)
    mixture$method_code <- 6L
  }
  bayesrc <- .mtblr_bayesrc_controls(
    semantics$prior_kernel, annotations, aligned$marker_ids, pattern_spec,
    mixture, add_intercept, standardize_annotations,
    center_binary_annotations, alpha_init, sigmaSqAlpha_init,
    intercept_flat, sigmaSqAlpha_a, sigmaSqAlpha_b, pi_floor,
    alpha_update_every, updateAlpha)
  extended_plan <- .blr_mtblr_extended_plan(
    conv, aligned$marker_ids, st$trait_names, method, mixture, bayesrc,
    updateB, updateE, updatePi, "diagonal", chain$nchains, chain$nit)
  conv$selected_markers_resolved <- extended_plan$selected
  mod <- mixture$patterns
  set_spec <- .mtblr_sets(sets, st$m)
  h2 <- rep(h2, length.out = st$nt); if (any(!is.finite(h2)) || any(h2 <= 0 | h2 >= 1)) stop("h2 must be in (0, 1).", call. = FALSE)
  vy <- st$yy / (st$n - 1)
  calibration_inputs <- .mtblr_calibration_inputs(mixture, bayesrc, st$m)
  calibrated <- .mtblr_prior_calibration(
   vy, h2, nub, nue, calibration_inputs$patterns,
   calibration_inputs$initial, calibration_inputs$prior,
   calibration_inputs$gamma, calibration_inputs$marker_scale,
   vg, vb, ve, ssb_prior, sse_prior,
   calibration_inputs$component_probability_source,
   calibration_inputs$annotation_probability_policy)
  vg <- calibrated$vg; vb <- calibrated$vb; ve <- calibrated$ve
  ssb_prior <- calibrated$ssb_prior; sse_prior <- calibrated$sse_prior
  if (any(ve[row(ve) != col(ve)] != 0) ||
      any(sse_prior[row(sse_prior) != col(sse_prior)] != 0)) {
   stop("ve and sse_prior must be diagonal for the public MT CSR route.",
        call. = FALSE)
  }
  initialization <- .mtblr_bayesr_initialization(
    beta,b,state,mixture$component_init,pattern_spec$matrix,st$m,st$nt,
    semantics$prior_kernel)
  b <- initialization$b
  trait_metadata <- if (is.null(trait_metadata)) data.frame(trait_id = st$trait_names) else as.data.frame(trait_metadata, stringsAsFactors = FALSE)
  if (is.null(trait_metadata$trait_id)) trait_metadata$trait_id <- st$trait_names
  if (nrow(trait_metadata) != st$nt || !identical(as.character(trait_metadata$trait_id), st$trait_names) || anyDuplicated(trait_metadata$trait_id)) stop("trait_metadata trait_id must uniquely match trait order.", call. = FALSE)
  for (x in c("study_id", "ancestry", "population", "ld_reference")) if (is.null(trait_metadata[[x]])) trait_metadata[[x]] <- NA_character_
  trait_metadata$sample_size <- st$n; trait_metadata$ld_prefix <- rep(native_prefix, length.out = st$nt); trait_metadata$ld_sharing_mode <- sharing
  operator_bytes <- sum(file.info(unlist(lapply(unique(native_prefix), function(x)
    paste0(x, c(".row_ptr.u64.bin", ".col_idx.u32.0based.bin",
                ".values.f32.bin")))))$size, na.rm = TRUE)
  memory <- .blr_memory_estimate(
    "mtblr", "csr", st$m, st$nt, chain$nchains, chain$ncores,
    chain$nit, chain$nit + chain$nburn, chain$keep_chains,
    convergence_quantities = if (conv$compute || conv$keep_traces) 5L * st$nt else 0L,
    keep_traces = conv$keep_traces, operator_bytes = operator_bytes)
  memory <- .mtblr_bayesr_memory(
    memory,method,st$m,chain$nchains,chain$ncores,chain$nit+chain$nburn,
    nrow(mod$matrix),mixture$component_count)
  memory <- .mtblr_bayesrc_memory(
    memory,bayesrc,st$m,chain$nchains,chain$ncores,
    chain$nit+chain$nburn)
  memory <- .blr_add_extended_memory(memory, extended_plan)
  .blr_memory_warning(memory, memory_warning_gb, conv$mode,
                      conv$compute || conv$keep_traces, conv$keep_traces)
  native_execution <- mtblr_csr_chains_raw_internal(
    st$wy, st$ww, st$yy,
    lapply(seq_len(st$nt), function(t) b[, t]), native_prefix,
    set_spec$native, vb, ve,
    lapply(seq_len(st$nt), function(i) ssb_prior[, i]),
    lapply(seq_len(st$nt), function(i) sse_prior[, i]), mod$native,
    mod$probabilities, nub, nue, updateB, updateE, updatePi, st$n,
    chain$nit, chain$nburn, chain$nthin, chain$seed, mixture$method_code,
    chain$nchains, chain$ncores, chain$chain_seeds_native,
    mixture$joint_component,mixture$joint_multiplier,mixture$joint_names,
    mixture$component_count,mixture$marker_scale,initialization$component,
    mixture$pi_prior,
    lapply(seq_len(st$nt), function(t) initialization$beta[,t]),
    lapply(seq_len(st$nt), function(t) initialization$state[,t]),
    bayesrc$annotations,bayesrc$alpha_init,bayesrc$sigma_alpha_init,
    bayesrc$pattern_pi_init,bayesrc$pattern_pi_prior,bayesrc$updateAlpha,
    bayesrc$intercept_flat,bayesrc$sigma_alpha_a,bayesrc$sigma_alpha_b,
    bayesrc$pi_floor,bayesrc$alpha_update_every,
    extended_plan$native$convergence_covariance,
    extended_plan$native$convergence_probability,
    extended_plan$native$convergence_annotations,
    extended_plan$native$convergence_full_probability,
    extended_plan$native$convergence_markers,
    extended_plan$native$convergence_b,
    extended_plan$native$convergence_d,
    extended_plan$native$convergence_component)
  if (!is.null(bayesrc$model_parameters))
    native_execution$raws <- lapply(native_execution$raws,
      .mtblr_bayesrc_enrich_raw, bayesrc = bayesrc, method = method,
      updatePi = updatePi)
  execution <- .mtblr_summary_multichain(
    native_execution, chain, conv, st$trait_names, method, "csr",
    updateB, updateE, mixture$model_parameters, extended_plan)
  raw <- execution$raw
  raw$model$names <- mod$names; raw$pi$names <- mod$names
  raw$data <- list(marker_metadata = aligned$marker_metadata, trait_metadata = trait_metadata,
    sample_size = st$n, ld_prefix = native_prefix, ld_sharing_mode = sharing,
    scale = "standardized_genotype", data_level = "summary_statistics",
    effect_maf_source = maf_info$effect_maf_source,
    effect_maf_population = maf_info$effect_maf_population,
    effect_maf_alignment_status = maf_info$effect_maf_alignment_status,
    effect_maf_fallback_used = maf_info$effect_maf_fallback_used,
    annotation_source = bayesrc$metadata$annotation_source %||% "not_applicable",
    annotation_marker_alignment_status =
      bayesrc$metadata$annotation_marker_alignment_status %||% "not_applicable")
  raw$alignment <- list(marker_policy = marker_policy, intersection_policy = "error", per_trait = aligned$report,
    orientation_status = aligned$report$allele_status)
  input <- list(method = method, model = method,
    backend = paste0("mt_csr_", semantics$prior_kernel),
    prior_kernel = semantics$prior_kernel,
    data_level = "summary_statistics",
    effect_scale_policy = if (!is.null(maf_effect_s))
      if (semantics$prior_kernel %in% c("bayesr", "bayesrc")) "component_maf_s" else "maf_s"
      else if (semantics$prior_kernel %in% c("bayesr", "bayesrc")) "component" else "unit",
    model_semantics_version = 2L,
    model_semantics = "s_prefix_means_summary_statistics",
    m = st$m, nt = st$nt, n = st$n, h2 = h2, nub = nub, nue = nue,
    nit = chain$nit, nburn = chain$nburn, nthin = chain$nthin,
    seed = chain$seed, nchains = chain$nchains, ncores = chain$ncores,
    chain_seeds_requested = chain$chain_seeds_requested,
    chain_seeds_resolved = execution$seeds,
    keep_chains = chain$keep_chains, convergence = conv$mode,
    convergence_control = conv[c("warn", "rhat_threshold",
      "ess_per_chain_threshold", "mcse_mean_over_sd_threshold",
      "keep_traces", "extended_groups_requested",
      "extended_groups_resolved", "selected_markers",
      "selected_markers_resolved",
      "selected_marker_quantities", "full_probability_states",
      "max_trace_gb", "allow_large_traces")], memory_warning_gb = memory_warning_gb,
    updateB = updateB, updateE = updateE, updatePi = updatePi,
    updateAlpha = bayesrc$updateAlpha,
    annotation_policy = if (is.null(bayesrc$model_parameters)) "global" else
      "annotation_probit_stick",
    models = mod$matrix, model_names = mod$names, pimodels = mod$probabilities,
    mixture_var = mixture$mixture_var, joint_pi_prior = mixture$pi_prior,
    maf_effect_s = mixture$maf_effect_s,
    estimate_maf_effect_s = estimate_maf_effect_s,
    effect_maf_source = maf_info$effect_maf_source,
    effect_maf_population = maf_info$effect_maf_population,
    effect_maf_alignment_status = maf_info$effect_maf_alignment_status,
    effect_maf_fallback_used = maf_info$effect_maf_fallback_used,
    sets = set_spec$public,
    ld_prefix = native_prefix, ld_sharing_mode = sharing, ld_reference = trait_metadata$ld_reference,
    marker_policy = marker_policy, marker_intersection_policy = "error", scale = "standardized_genotype",
    sample_overlap = "not_modeled", phenotype_crossproduct_policy = "marginal_yy_only",
    residual_covariance_policy = "diagonal", trait_metadata = trait_metadata, alignment = raw$alignment)
  input <- c(input, calibrated$calibration)
  if (isTRUE(verbose)) print(input[c("backend", "m", "nt", "ld_sharing_mode", "seed", "nchains", "ncores")])
  fit <- .as_mtblr_fit(raw, aligned$marker_ids, st$trait_names,
                       aligned$marker_metadata, trait_metadata,
                       raw$alignment, input)
  fit$chains <- execution$chains
  fit$convergence <- execution$convergence
  fit["convergence_traces"] <- list(execution$convergence_traces)
  fit$input$memory_estimate <- memory
  fit <- .mtblr_bayesr_format_fit(fit, mixture$model_parameters)
  fit <- .mtblr_bayesrc_format_fit(fit, raw$annotations, bayesrc)
  if (isTRUE(bayesrc$maf_annotation_overlap) && !is.null(maf_effect_s))
    warning("MAF-derived annotations and maf_effect_s are both active; MAF may influence component probabilities and effect-size variance.", call. = FALSE)
  messages <- .blr_convergence_warning_messages(
    fit$convergence, if (conv$mode == "core") "core" else "auto",
    "mtblr", "csr")
  if (isTRUE(conv$warn) && conv$mode != "none" &&
      !(conv$mode == "auto" && chain$nchains == 1L) && length(messages)) {
    warning(messages[[1L]], call. = FALSE)
  }
  .blr_finalize_fit(
    fit,
    "mtblr", method, "csr", data = raw$data,
    diagnostics = raw$diagnostics, memory_estimate = memory)
}
