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
  names(wy) <- names(ww) <- trait_names
  list(wy = lapply(wy, as.numeric), ww = lapply(ww, as.numeric), yy = as.numeric(yy),
       n = as.integer(n), m = m, nt = nt, trait_names = trait_names,
       marker_ids = marker_by_trait, marker_metadata = marker_metadata,
       source = source, scale = scale)
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
    } else if (identical(stats$source[t], "make_summary_stats") && isTRUE(ld$descriptors[[t]]$by_construction)) {
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
  for (x in c("vbs", "vgs", "ves")) if (!identical(dim(raw$trace[[x]]), c(ntr, nt))) stop("Invalid mtblr_raw trace dimensions.", call. = FALSE)
  for (x in c("covb", "covg", "cove", "vb", "vg", "ve")) if (!identical(dim(raw$variance[[x]]), c(nt, nt))) stop("Invalid mtblr_raw covariance dimensions.", call. = FALSE)
  if (length(raw$pi$final) != nm || length(raw$pi$mean) != nm || !identical(dim(raw$model$patterns), c(nm, nt))) stop("Invalid mtblr_raw probability/model dimensions.", call. = FALSE)
  if (any(!raw$model$patterns %in% 0:1) || length(raw$marker$order) != m) stop("Invalid mtblr_raw patterns or marker order.", call. = FALSE)
  raw
}

.as_mtblr_fit <- function(raw, marker_ids, trait_names, marker_metadata, trait_metadata, alignment, input) {
  raw <- .validate_mtblr_raw(raw)
  dimnames_marker <- list(marker_ids, trait_names); dimnames_trace <- list(paste0("Iter", seq_len(raw$meta$n_trace)), trait_names)
  marker <- lapply(raw$marker[c("bm", "dm", "wy", "r", "b", "state")], function(x) { dimnames(x) <- dimnames_marker; x })
  trace <- lapply(raw$trace, function(x) { dimnames(x) <- dimnames_trace; x })
  variance <- lapply(raw$variance, function(x) { dimnames(x) <- list(trait_names, trait_names); x })
  model_names <- raw$model$names; names(raw$pi$final) <- names(raw$pi$mean) <- model_names
  o <- matrix(raw$marker$order, nrow = length(marker_ids), ncol = length(trait_names), dimnames = dimnames_marker)
  fit <- c(marker[c("bm", "dm", "wy", "r", "b")], list(d = marker$state, o = o), trace, variance,
           list(pi = raw$pi$final, pim = raw$pi$mean, marker_order = raw$marker$order,
                model_patterns = raw$model$patterns, marker_metadata = marker_metadata,
                trait_metadata = trait_metadata, alignment = alignment, input = input,
                raw_schema_version = 1L))
  fit["rb"] <- list(if (sum(diag(fit$covb)) > 0) cov2cor(fit$covb) else NULL)
  fit["rg"] <- list(if (sum(diag(fit$covg)) > 0) cov2cor(fit$covg) else NULL)
  fit["re"] <- list(if (sum(diag(fit$cove)) > 0) cov2cor(fit$cove) else NULL)
  class(fit) <- c("mtblr_fit", "list")
  fit
}

#' Trait-specific multivariate BayesC with sparse LD
#'
#' Fits the corrected serial multivariate BayesC model from standardized
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
#' @param method Must be `"bayesC"`.
#' @param n Optional sample sizes overriding identical values only.
#' @param sets Optional disjoint complete marker partition (1-based).
#' @param b Optional initial marker-by-trait effects.
#' @param h2 Heritability scalar or one value per trait.
#' @param pi Initial non-null probability or complete pattern probabilities.
#' @param models,pimodels Pattern matrix and probabilities.
#' @param vg,vb,ve Initial covariance matrices.
#' @param ssb_prior,sse_prior Covariance prior scale matrices.
#' @param updateB,updateE,updatePi Update controls.
#' @param nub,nue Prior degrees of freedom.
#' @param nit,nburn,nthin MCMC controls.
#' @param seed Explicit native RNG seed.
#' @param verbose Print resolved execution metadata.
#' @return An object of class `mtblr_fit`.
#' @export
mtblr_csr <- function(stats, Glist = NULL, ld_prefix = NULL, ld_metadata = NULL,
  trait_metadata = NULL, marker_policy = c("strict", "reorder_stats"),
  sample_overlap = "not_modeled", method = "bayesC", n = NULL, sets = NULL,
  b = NULL, h2 = 0.5, pi = 0.001, models = NULL, pimodels = NULL,
  vg = NULL, vb = NULL, ve = NULL, ssb_prior = NULL, sse_prior = NULL,
  updateB = TRUE, updateE = TRUE, updatePi = TRUE, nub = 4, nue = 4,
  nit = 1000, nburn = 500, nthin = 1, seed = 1, verbose = FALSE) {
  marker_policy <- match.arg(marker_policy)
  if (!identical(sample_overlap, "not_modeled")) stop("sample_overlap must be exactly 'not_modeled'.", call. = FALSE)
  if (!identical(method, "bayesC")) stop("Only method = 'bayesC' is supported.", call. = FALSE)
  st <- .mtblr_normalize_stats(stats); if (!is.null(n) && !identical(as.integer(rep(n, length.out = st$nt)), st$n)) stop("n conflicts with stats sample sizes.", call. = FALSE)
  ld <- .mtblr_resolve_ld(Glist, ld_prefix, ld_metadata, st$nt)
  aligned <- .mtblr_align(st, ld, marker_policy); st <- aligned$stats
  shared <- length(ld$prefixes) == 1L
  sharing <- if (shared && all(vapply(st$ww[-1L], identical, logical(1), st$ww[[1L]]))) "fully_shared_operator" else if (shared) "shared_correlation_reference" else "trait_specific_reference"
  native_prefix <- if (sharing == "shared_correlation_reference") rep(ld$prefixes, st$nt) else ld$prefixes
  mod <- .mtblr_models(models, pimodels, pi, st$nt); set_spec <- .mtblr_sets(sets, st$m)
  h2 <- rep(h2, length.out = st$nt); if (any(!is.finite(h2)) || any(h2 <= 0 | h2 >= 1)) stop("h2 must be in (0, 1).", call. = FALSE)
  vy <- st$yy / (st$n - 1); vg0 <- diag(vy * h2, st$nt); ve0 <- diag(vy * (1 - h2), st$nt)
  vb0 <- diag((vy * h2) / (st$m * max(1e-12, 1 - mod$probabilities[1L])), st$nt)
  vg <- .mtblr_cov(vg, vg0, "vg", st$nt); vb <- .mtblr_cov(vb, vb0, "vb", st$nt)
  ve <- .mtblr_cov(ve, ve0, "ve", st$nt, TRUE)
  ssb_prior <- .mtblr_cov(ssb_prior, ((nub - 2) / nub) * vg / (st$m * max(1e-12, 1 - mod$probabilities[1L])), "ssb_prior", st$nt)
  sse_prior <- .mtblr_cov(sse_prior, ((nue - 2) / nue) * ve, "sse_prior", st$nt, TRUE)
  if (is.null(b)) b <- matrix(0, st$m, st$nt)
  if (is.list(b)) b <- do.call(cbind, b); b <- as.matrix(b)
  if (!identical(dim(b), c(st$m, st$nt)) || any(!is.finite(b))) stop("b must be a finite m by nt matrix or trait list.", call. = FALSE)
  for (x in c("nit", "nburn", "nthin", "seed")) { value <- get(x); if (length(value) != 1L || !is.finite(value) || value != as.integer(value) || (x != "nburn" && value <= 0) || (x == "nburn" && value < 0)) stop(x, " must be an integer-compatible scalar in its valid range.", call. = FALSE) }
  trait_metadata <- if (is.null(trait_metadata)) data.frame(trait_id = st$trait_names) else as.data.frame(trait_metadata, stringsAsFactors = FALSE)
  if (is.null(trait_metadata$trait_id)) trait_metadata$trait_id <- st$trait_names
  if (nrow(trait_metadata) != st$nt || !identical(as.character(trait_metadata$trait_id), st$trait_names) || anyDuplicated(trait_metadata$trait_id)) stop("trait_metadata trait_id must uniquely match trait order.", call. = FALSE)
  for (x in c("study_id", "ancestry", "population", "ld_reference")) if (is.null(trait_metadata[[x]])) trait_metadata[[x]] <- NA_character_
  trait_metadata$sample_size <- st$n; trait_metadata$ld_prefix <- rep(native_prefix, length.out = st$nt); trait_metadata$ld_sharing_mode <- sharing
  raw <- mtblr_csr_raw_internal(st$wy, st$ww, st$yy,
    lapply(seq_len(st$nt), function(t) b[, t]), native_prefix, set_spec$native,
    vb, ve, lapply(seq_len(st$nt), function(i) ssb_prior[, i]),
    lapply(seq_len(st$nt), function(i) sse_prior[, i]), mod$native, mod$probabilities,
    nub, nue, updateB, updateE, updatePi, st$n, as.integer(nit), as.integer(nburn),
    as.integer(nthin), as.integer(seed), 4L)
  raw$model$names <- mod$names; raw$pi$names <- mod$names
  raw$data <- list(marker_metadata = aligned$marker_metadata, trait_metadata = trait_metadata,
    sample_size = st$n, ld_prefix = native_prefix, ld_sharing_mode = sharing, scale = "standardized_genotype")
  raw$alignment <- list(marker_policy = marker_policy, intersection_policy = "error", per_trait = aligned$report,
    orientation_status = aligned$report$allele_status)
  input <- list(method = "bayesc", model = "bayesc", backend = "mt_csr_bayesc", data_level = "summary",
    m = st$m, nt = st$nt, n = st$n, h2 = h2, nub = nub, nue = nue, nit = nit, nburn = nburn,
    nthin = nthin, seed = seed, updateB = updateB, updateE = updateE, updatePi = updatePi,
    models = mod$matrix, model_names = mod$names, pimodels = mod$probabilities, sets = set_spec$public,
    ld_prefix = native_prefix, ld_sharing_mode = sharing, ld_reference = trait_metadata$ld_reference,
    marker_policy = marker_policy, marker_intersection_policy = "error", scale = "standardized_genotype",
    sample_overlap = "not_modeled", phenotype_crossproduct_policy = "marginal_yy_only",
    residual_covariance_policy = "diagonal", trait_metadata = trait_metadata, alignment = raw$alignment)
  if (isTRUE(verbose)) print(input[c("backend", "m", "nt", "ld_sharing_mode", "seed")])
  .as_mtblr_fit(raw, aligned$marker_ids, st$trait_names, aligned$marker_metadata, trait_metadata, raw$alignment, input)
}
