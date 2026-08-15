.blr_phase2_operator_types <- c(
  "bed", "dense", "dense_cross_product", "csr", "block_eigen",
  "full_rank_block_eigen", "retained_rank_block_eigen"
)

.blr_phase2_likelihood_regimes <- c(
  "common_sample", "independent_summary", "overlap_aware"
)

.blr_new_global_marker_map <- function(marker_ids, alleles) {
  marker_ids <- .blr_ids(marker_ids, "global marker IDs")
  .blr_validate_allele_table(
    alleles, marker_ids, "global marker alleles", require_coding = TRUE)
  structure(list(
    marker_ids = marker_ids,
    alleles = alleles,
    global_index = stats::setNames(seq_along(marker_ids), marker_ids)
  ), class = c("blr_global_marker_map_v1", "list"))
}

.blr_validate_global_marker_map <- function(x) {
  .blr_exact_fields(x, "global marker map",
                    c("marker_ids", "alleles", "global_index"))
  marker_ids <- .blr_ids(x$marker_ids, "global marker map$marker_ids")
  .blr_validate_allele_table(
    x$alleles, marker_ids, "global marker map$alleles",
    require_coding = TRUE)
  expected <- stats::setNames(seq_along(marker_ids), marker_ids)
  if (!identical(x$global_index, expected)) {
    stop("global marker map$global_index must preserve declared marker order.",
         call. = FALSE)
  }
  invisible(TRUE)
}

.blr_new_operator_storage_ref <- function(kind, payload) {
  .blr_character_scalar(kind, "operator storage kind")
  owner <- new.env(parent = emptyenv())
  owner$kind <- kind
  owner$payload <- payload
  class(owner) <- c("blr_operator_storage_ref_v1", "environment")
  lockEnvironment(owner, bindings = TRUE)
  owner
}

.blr_operator_storage_payload <- function(storage) {
  if (!inherits(storage, "blr_operator_storage_ref_v1")) {
    stop("operator storage is not a Phase 2 immutable storage reference.",
         call. = FALSE)
  }
  storage$payload
}

.blr_operator_alleles <- function(marker_ids, effect, other) {
  data.frame(
    marker_id = marker_ids,
    effect = rep_len(effect, length(marker_ids)),
    other = rep_len(other, length(marker_ids)),
    stringsAsFactors = FALSE
  )
}

.blr_new_operator_resource <- function(
    resource_id, operator_type, marker_ids, alleles,
    genotype_coding, centering, standardization, operator_scale,
    storage, block_eigen = NULL, approximation, provenance = list()) {
  resource <- list(
    resource_id = resource_id,
    operator_type = operator_type,
    marker_ids = marker_ids,
    alleles = alleles,
    genotype_coding = genotype_coding,
    centering = centering,
    standardization = standardization,
    operator_scale = operator_scale,
    storage = storage,
    block_eigen = block_eigen,
    approximation = approximation,
    provenance = provenance
  )
  .blr_validate_operator_resource(resource)
  class(resource) <- c("blr_operator_resource_v1", "list")
  resource
}

.blr_validate_dense_storage <- function(storage, marker_ids) {
  matrix <- .blr_operator_storage_payload(storage)
  m <- length(marker_ids)
  if (!is.matrix(matrix) || !is.numeric(matrix) ||
      !identical(dim(matrix), c(m, m)) || any(!is.finite(matrix)) ||
      !isTRUE(all.equal(matrix, t(matrix), tolerance = 1e-12)) ||
      any(diag(matrix) <= 0)) {
    stop("dense cross-product storage must be finite, symmetric, have positive diagonals, and match marker dimensions.",
         call. = FALSE)
  }
  eigenvalues <- eigen(matrix, symmetric = TRUE, only.values = TRUE)$values
  tolerance <- 1e-10 * max(1, max(abs(eigenvalues)))
  if (min(eigenvalues) < -tolerance) {
    stop("dense cross-product storage must be positive semidefinite.",
         call. = FALSE)
  }
  if (!is.null(dimnames(matrix)) &&
      !identical(dimnames(matrix), list(marker_ids, marker_ids))) {
    stop("dense cross-product dimnames must follow resource marker order.",
         call. = FALSE)
  }
  invisible(TRUE)
}

.blr_validate_csr_storage <- function(storage, marker_ids) {
  csr <- .blr_operator_storage_payload(storage)
  .blr_exact_fields(csr, "CSR storage", c(
    "row_ptr", "column_index", "values", "diagonal", "complete"))
  m <- length(marker_ids)
  valid <- is.numeric(csr$row_ptr) && length(csr$row_ptr) == m + 1L &&
    !anyNA(csr$row_ptr) && all(is.finite(csr$row_ptr)) &&
    all(csr$row_ptr == floor(csr$row_ptr)) && csr$row_ptr[[1L]] == 0 &&
    all(diff(csr$row_ptr) >= 0) &&
    tail(csr$row_ptr, 1L) == length(csr$values) &&
    is.numeric(csr$column_index) &&
    length(csr$column_index) == length(csr$values) &&
    !anyNA(csr$column_index) && all(is.finite(csr$column_index)) &&
    all(csr$column_index == floor(csr$column_index)) &&
    all(csr$column_index >= 1L & csr$column_index <= m) &&
    is.numeric(csr$values) && !anyNA(csr$values) &&
    all(is.finite(csr$values)) &&
    is.numeric(csr$diagonal) && length(csr$diagonal) == m &&
    !anyNA(csr$diagonal) && all(is.finite(csr$diagonal)) &&
    all(csr$diagonal > 0) && is.logical(csr$complete) &&
    length(csr$complete) == 1L && !is.na(csr$complete)
  if (!valid) {
    stop("CSR storage has invalid dimensions, indices, values, diagonal, or completeness metadata.",
         call. = FALSE)
  }
  reconstructed <- matrix(0, m, m)
  for (row in seq_len(m)) {
    first <- csr$row_ptr[[row]] + 1L
    last <- csr$row_ptr[[row + 1L]]
    if (first <= last) {
      columns <- csr$column_index[first:last]
      if (anyDuplicated(columns)) {
        stop("CSR storage cannot repeat a column within one row.",
             call. = FALSE)
      }
      reconstructed[row, columns] <- csr$values[first:last]
    }
  }
  if (!isTRUE(all.equal(reconstructed, t(reconstructed), tolerance = 1e-12)) ||
      !isTRUE(all.equal(unname(diag(reconstructed)),
                       unname(csr$diagonal),
                       tolerance = 1e-12))) {
    stop("CSR storage must be symmetric and its explicit diagonal must agree with stored entries.",
         call. = FALSE)
  }
  invisible(TRUE)
}

.blr_validate_bed_storage <- function(storage, marker_ids) {
  bed <- .blr_operator_storage_payload(storage)
  .blr_exact_fields(bed, "packed-BED storage", c(
    "bed_files", "source_sample_count", "sample_ids", "selected_rows",
    "selected_columns"))
  valid_files <- is.character(bed$bed_files) && length(bed$bed_files) > 0L &&
    !anyNA(bed$bed_files) && all(nzchar(bed$bed_files))
  valid_n <- is.numeric(bed$source_sample_count) &&
    length(bed$source_sample_count) == 1L &&
    !is.na(bed$source_sample_count) && is.finite(bed$source_sample_count) &&
    bed$source_sample_count == floor(bed$source_sample_count) &&
    bed$source_sample_count > 1L
  valid_samples <- is.character(bed$sample_ids) &&
    length(bed$sample_ids) > 1L && !anyNA(bed$sample_ids) &&
    all(nzchar(bed$sample_ids)) && !anyDuplicated(bed$sample_ids)
  valid_rows <- is.numeric(bed$selected_rows) &&
    length(bed$selected_rows) == length(bed$sample_ids) &&
    !anyNA(bed$selected_rows) && all(is.finite(bed$selected_rows)) &&
    all(bed$selected_rows == floor(bed$selected_rows)) &&
    all(bed$selected_rows >= 1L &
          bed$selected_rows <= bed$source_sample_count) &&
    !anyDuplicated(bed$selected_rows)
  valid_columns <- is.list(bed$selected_columns) &&
    length(bed$selected_columns) == length(bed$bed_files) &&
    sum(lengths(bed$selected_columns)) == length(marker_ids) &&
    all(vapply(bed$selected_columns, function(columns) {
      is.numeric(columns) && length(columns) > 0L && !anyNA(columns) &&
        all(is.finite(columns)) && all(columns == floor(columns)) &&
        all(columns >= 1L) && !anyDuplicated(columns)
    }, logical(1)))
  if (!valid_files || !valid_n || !valid_samples || !valid_rows ||
      !valid_columns) {
    stop("packed-BED storage has invalid files, selected samples, or marker columns.",
         call. = FALSE)
  }
  invisible(TRUE)
}

.blr_validate_block_eigen_storage <- function(resource) {
  descriptor <- resource$block_eigen
  .blr_exact_fields(descriptor, "block-eigen descriptor", c("blocks"))
  if (!is.list(descriptor$blocks) || !length(descriptor$blocks)) {
    stop("block-eigen descriptor must contain one or more blocks.",
         call. = FALSE)
  }
  m <- length(resource$marker_ids)
  seen <- integer()
  retained <- logical(length(descriptor$blocks))
  block_ids <- character(length(descriptor$blocks))
  for (index in seq_along(descriptor$blocks)) {
    block <- descriptor$blocks[[index]]
    .blr_exact_fields(block, "block-eigen block", c(
      "block_id", "marker_indices", "eigenvectors", "eigenvalues",
      "retained_rank"))
    .blr_character_scalar(block$block_id, "block-eigen block$block_id")
    block_ids[[index]] <- block$block_id
    ids <- block$marker_indices
    size <- length(ids)
    valid_ids <- is.numeric(ids) && size > 0L && !anyNA(ids) &&
      all(is.finite(ids)) && all(ids == floor(ids)) &&
      all(ids >= 1L & ids <= m) && !anyDuplicated(ids)
    rank <- block$retained_rank
    valid_rank <- is.numeric(rank) && length(rank) == 1L && !is.na(rank) &&
      is.finite(rank) && rank == floor(rank) && rank >= 1L && rank <= size
    valid_values <- is.numeric(block$eigenvalues) && valid_rank &&
      length(block$eigenvalues) == rank && !anyNA(block$eigenvalues) &&
      all(is.finite(block$eigenvalues)) && all(block$eigenvalues >= 0)
    valid_vectors <- is.matrix(block$eigenvectors) &&
      is.numeric(block$eigenvectors) && valid_rank &&
      identical(dim(block$eigenvectors), c(size, as.integer(rank))) &&
      all(is.finite(block$eigenvectors))
    if (!valid_ids || !valid_values || !valid_vectors) {
      stop("block-eigen block has invalid membership, eigenvectors, eigenvalues, or retained rank.",
           call. = FALSE)
    }
    if (!isTRUE(all.equal(
        crossprod(block$eigenvectors), diag(as.integer(rank)),
        tolerance = 1e-10))) {
      stop("block-eigen eigenvectors must have orthonormal retained columns.",
           call. = FALSE)
    }
    seen <- c(seen, as.integer(ids))
    retained[[index]] <- rank < size
  }
  if (!identical(sort(seen), seq_len(m)) || anyDuplicated(seen)) {
    stop("block-eigen blocks must partition all resource markers exactly once.",
         call. = FALSE)
  }
  if (anyDuplicated(block_ids)) {
    stop("block-eigen block IDs must be unique.", call. = FALSE)
  }
  if (identical(resource$operator_type, "full_rank_block_eigen") &&
      any(retained)) {
    stop("a full-rank block-eigen resource cannot retain fewer than all block components.",
         call. = FALSE)
  }
  if (identical(resource$operator_type, "retained_rank_block_eigen") &&
      !any(retained)) {
    stop("a retained-rank block-eigen resource must truncate at least one block.",
         call. = FALSE)
  }
  invisible(TRUE)
}

.blr_validate_operator_resource <- function(resource) {
  .blr_exact_fields(resource, "operator resource", c(
    "resource_id", "operator_type", "marker_ids", "alleles",
    "genotype_coding", "centering", "standardization", "operator_scale",
    "storage", "block_eigen", "approximation", "provenance"))
  .blr_character_scalar(resource$resource_id, "operator resource$resource_id")
  .blr_character_scalar(resource$operator_type,
                        "operator resource$operator_type",
                        .blr_phase2_operator_types)
  marker_ids <- .blr_ids(resource$marker_ids,
                         "operator resource$marker_ids")
  .blr_validate_allele_table(resource$alleles, marker_ids,
                             "operator resource$alleles")
  for (field in c("genotype_coding", "centering", "standardization",
                  "operator_scale", "approximation")) {
    .blr_character_scalar(resource[[field]],
                          paste0("operator resource$", field))
  }
  .blr_exact_names(resource$provenance, "operator resource$provenance")
  type <- resource$operator_type
  if (inherits(resource$storage, "blr_operator_storage_ref_v1")) {
    storage_kind <- resource$storage$kind
    if (type %in% c("dense", "dense_cross_product")) {
      if (!identical(storage_kind, "dense_cross_product")) {
        stop("dense resources require dense_cross_product storage.",
             call. = FALSE)
      }
      .blr_validate_dense_storage(resource$storage, marker_ids)
    } else if (identical(type, "csr")) {
      if (!identical(storage_kind, "csr_cross_product")) {
        stop("CSR resources require csr_cross_product storage.", call. = FALSE)
      }
      .blr_validate_csr_storage(resource$storage, marker_ids)
    } else if (type %in% c("full_rank_block_eigen",
                           "retained_rank_block_eigen")) {
      if (!identical(storage_kind, "block_eigen")) {
        stop("block-eigen resources require block_eigen storage.",
             call. = FALSE)
      }
      .blr_validate_block_eigen_storage(resource)
    } else if (identical(type, "bed")) {
      if (!identical(storage_kind, "packed_bed")) {
        stop("BED resources require packed_bed storage.", call. = FALSE)
      }
      .blr_validate_bed_storage(resource$storage, marker_ids)
    }
  } else if (!is.list(resource$storage)) {
    stop("operator resource$storage must be an immutable reference or a validated view descriptor.",
         call. = FALSE)
  }
  if (type %in% c("full_rank_block_eigen", "retained_rank_block_eigen") &&
      is.null(resource$block_eigen)) {
    stop("block-eigen resources require block metadata.", call. = FALSE)
  }
  invisible(TRUE)
}

.blr_new_likelihood_provider <- function(
    provider_id, trait_ids, operator_resource_id, local_to_global,
    sufficient_statistics, sample_size, likelihood_regime,
    residual_contract, population, effect_scale,
    overlap_group = NULL, provenance = list()) {
  provider <- list(
    provider_id = provider_id,
    trait_ids = trait_ids,
    operator_resource_id = operator_resource_id,
    local_to_global = local_to_global,
    sufficient_statistics = sufficient_statistics,
    sample_size = sample_size,
    likelihood_regime = likelihood_regime,
    residual_contract = residual_contract,
    population = population,
    effect_scale = effect_scale,
    overlap_group = overlap_group,
    provenance = provenance
  )
  .blr_validate_likelihood_provider(provider)
  class(provider) <- c("blr_likelihood_provider_v1", "list")
  provider
}

.blr_validate_likelihood_provider <- function(provider) {
  .blr_exact_fields(provider, "likelihood provider", c(
    "provider_id", "trait_ids", "operator_resource_id", "local_to_global",
    "sufficient_statistics", "sample_size", "likelihood_regime",
    "residual_contract", "population", "effect_scale", "overlap_group",
    "provenance"))
  .blr_character_scalar(provider$provider_id, "provider$provider_id")
  .blr_ids(provider$trait_ids, "provider$trait_ids")
  .blr_character_scalar(provider$operator_resource_id,
                        "provider$operator_resource_id")
  .blr_character_scalar(provider$likelihood_regime,
                        "provider$likelihood_regime",
                        .blr_phase2_likelihood_regimes)
  .blr_character_scalar(provider$residual_contract,
                        "provider$residual_contract")
  if (!is.numeric(provider$sample_size) ||
      length(provider$sample_size) != length(provider$trait_ids) ||
      anyNA(provider$sample_size) || any(!is.finite(provider$sample_size)) ||
      any(provider$sample_size <= 0) ||
      !identical(names(provider$sample_size), provider$trait_ids)) {
    stop("provider$sample_size must be a finite positive trait-named vector.",
         call. = FALSE)
  }
  if (!is.null(provider$population)) {
    .blr_character_scalar(provider$population, "provider$population")
  }
  .blr_character_scalar(provider$effect_scale, "provider$effect_scale")
  if (!is.null(provider$overlap_group)) {
    .blr_character_scalar(provider$overlap_group, "provider$overlap_group")
  }
  .blr_exact_names(provider$sufficient_statistics,
                   "provider$sufficient_statistics")
  .blr_exact_names(provider$provenance, "provider$provenance")
  invisible(TRUE)
}

.blr_validate_provider_statistics <- function(provider, resource) {
  statistics <- provider$sufficient_statistics
  if (identical(statistics$source %||% NULL, "legacy_wrapper")) {
    return(invisible(TRUE))
  }
  marker_ids <- resource$marker_ids
  trait_ids <- provider$trait_ids
  if (!is.null(statistics$score)) {
    score <- statistics$score
    if (!is.matrix(score) || !is.numeric(score) ||
        !identical(dim(score), c(length(marker_ids), length(trait_ids))) ||
        any(!is.finite(score)) ||
        !identical(dimnames(score),
                   list(marker = marker_ids, trait = trait_ids))) {
      stop("provider score must be finite and follow resource marker and trait order.",
           call. = FALSE)
    }
  }
  if (!is.null(statistics$phenotype)) {
    phenotype <- statistics$phenotype
    if (!is.matrix(phenotype) || !is.numeric(phenotype) ||
        ncol(phenotype) != length(trait_ids) || any(!is.finite(phenotype)) ||
        !identical(colnames(phenotype), trait_ids) ||
        any(provider$sample_size != nrow(phenotype))) {
      stop("provider phenotype must be finite and follow declared sample and trait order.",
           call. = FALSE)
    }
    if (identical(resource$operator_type, "bed") &&
        inherits(resource$storage, "blr_operator_storage_ref_v1") &&
        !identical(rownames(phenotype),
                   .blr_operator_storage_payload(resource$storage)$sample_ids)) {
      stop("provider phenotype sample order must match packed-BED selected sample order.",
           call. = FALSE)
    }
  }
  if (is.null(statistics$score) && is.null(statistics$phenotype)) {
    stop("provider sufficient statistics require an aligned score or phenotype payload.",
         call. = FALSE)
  }
  invisible(TRUE)
}

.blr_validate_phase2_provider_alignment <- function(
    provider, resource, global_map) {
  .blr_validate_likelihood_provider(provider)
  .blr_validate_operator_resource(resource)
  .blr_validate_global_marker_map(global_map)
  map <- provider$local_to_global
  marker_ids <- resource$marker_ids
  if (!is.numeric(map) || length(map) != length(marker_ids) || anyNA(map) ||
      any(!is.finite(map)) || any(map != floor(map)) ||
      any(map < 1L | map > length(global_map$marker_ids)) ||
      anyDuplicated(map) || !identical(names(map), marker_ids) ||
      !identical(unname(global_map$marker_ids[map]), unname(marker_ids))) {
    stop("provider local_to_global must be a unique resource-ordered map into the global marker universe.",
         call. = FALSE)
  }
  aligned <- global_map$alleles[map, , drop = FALSE]
  for (field in c("effect", "other")) {
    if (!identical(as.character(resource$alleles[[field]]),
                   as.character(aligned[[field]]))) {
      stop("provider/resource alleles disagree with the global marker map.",
           call. = FALSE)
    }
  }
  known_coding <- unique(as.character(aligned$coding))
  if (length(known_coding) != 1L ||
      !identical(known_coding, resource$genotype_coding)) {
    stop("provider/resource genotype coding disagrees with the global marker map.",
         call. = FALSE)
  }
  invisible(TRUE)
}

.blr_new_provider_collection <- function(global_map, resources, providers,
                                         likelihood_regime, analysis_mode) {
  .blr_validate_global_marker_map(global_map)
  .blr_character_scalar(likelihood_regime, "likelihood_regime",
                        .blr_phase2_likelihood_regimes)
  .blr_character_scalar(analysis_mode, "analysis_mode", c(
    "single_trait", "independent_traits", "joint_multitrait"))
  .blr_exact_names(resources, "operator resources")
  .blr_exact_names(providers, "likelihood providers")
  if (!length(resources) || !length(providers)) {
    stop("provider collection requires resources and providers.", call. = FALSE)
  }
  resource_ids <- unname(vapply(resources, function(resource) {
    .blr_validate_operator_resource(resource)
    resource$resource_id
  }, character(1)))
  provider_ids <- unname(vapply(providers, function(provider) {
    .blr_validate_likelihood_provider(provider)
    provider$provider_id
  }, character(1)))
  if (anyDuplicated(resource_ids) || !identical(names(resources), resource_ids)) {
    stop("operator resources must be named by unique resource IDs.",
         call. = FALSE)
  }
  if (anyDuplicated(provider_ids) || !identical(names(providers), provider_ids)) {
    stop("likelihood providers must be named by unique provider IDs.",
         call. = FALSE)
  }
  for (provider in providers) {
    if (!provider$operator_resource_id %in% resource_ids) {
      stop("likelihood provider references an unknown operator resource.",
           call. = FALSE)
    }
    .blr_validate_phase2_provider_alignment(
      provider, resources[[provider$operator_resource_id]], global_map)
    .blr_validate_provider_statistics(
      provider, resources[[provider$operator_resource_id]])
    if (!identical(provider$likelihood_regime, likelihood_regime)) {
      stop("provider likelihood regimes must match the collection regime.",
           call. = FALSE)
    }
  }
  for (trait_id in unique(unlist(
      lapply(providers, `[[`, "trait_ids"), use.names = FALSE))) {
    owners <- providers[vapply(
      providers, function(provider) trait_id %in% provider$trait_ids,
      logical(1))]
    if (length(owners) > 1L) {
      effect_scales <- vapply(owners, `[[`, character(1), "effect_scale")
      if (length(unique(effect_scales)) != 1L) {
        stop("providers for one shared effect must declare compatible effect scales.",
             call. = FALSE)
      }
    }
  }
  if (identical(likelihood_regime, "independent_summary") &&
      any(lengths(lapply(providers, `[[`, "trait_ids")) != 1L)) {
    stop("independent_summary providers must each own one trait.",
         call. = FALSE)
  }
  if (identical(analysis_mode, "joint_multitrait") &&
      identical(likelihood_regime, "common_sample") &&
      !any(vapply(providers, function(x) length(x$trait_ids) > 1L,
                  logical(1)))) {
    stop("common-sample joint_multitrait data require one non-factorized multi-trait provider.",
         call. = FALSE)
  }
  structure(list(
    global_marker_map = global_map,
    operator_resources = resources,
    providers = providers,
    likelihood_regime = likelihood_regime,
    analysis_mode = analysis_mode
  ), class = c("blr_provider_collection_v1", "list"))
}

.blr_provider_presence <- function(collection) {
  provider_ids <- names(collection$providers)
  marker_ids <- collection$global_marker_map$marker_ids
  out <- matrix(FALSE, length(provider_ids), length(marker_ids),
                dimnames = list(provider = provider_ids, marker = marker_ids))
  for (id in provider_ids) {
    out[id, collection$providers[[id]]$local_to_global] <- TRUE
  }
  out
}

.blr_provider_resource <- function(collection, provider_id) {
  .blr_character_scalar(provider_id, "provider_id")
  if (!provider_id %in% names(collection$providers)) {
    stop("unknown provider_id.", call. = FALSE)
  }
  resource_id <- collection$providers[[provider_id]]$operator_resource_id
  collection$operator_resources[[resource_id]]
}

.blr_dense_resource <- function(resource_id, cross_product, marker_ids,
                                alleles, genotype_coding = "effect_allele_count",
                                provenance = list()) {
  .blr_new_operator_resource(
    resource_id, "dense_cross_product", marker_ids, alleles,
    genotype_coding, "declared", "declared", "cross_product",
    .blr_new_operator_storage_ref("dense_cross_product", cross_product),
    block_eigen = NULL, approximation = "exact_declared_operator",
    provenance = provenance)
}

.blr_csr_from_dense_resource <- function(
    resource_id, cross_product, marker_ids, alleles,
    zero_tolerance = 0, genotype_coding = "effect_allele_count",
    provenance = list()) {
  if (!is.numeric(zero_tolerance) || length(zero_tolerance) != 1L ||
      is.na(zero_tolerance) || !is.finite(zero_tolerance) ||
      zero_tolerance < 0) {
    stop("zero_tolerance must be one finite nonnegative value.",
         call. = FALSE)
  }
  dense_ref <- .blr_new_operator_storage_ref(
    "dense_cross_product", cross_product)
  .blr_validate_dense_storage(dense_ref, marker_ids)
  m <- length(marker_ids)
  row_ptr <- numeric(m + 1L)
  column <- integer()
  values <- numeric()
  for (row in seq_len(m)) {
    keep <- which(abs(cross_product[row, ]) > zero_tolerance)
    column <- c(column, keep)
    values <- c(values, cross_product[row, keep])
    row_ptr[row + 1L] <- length(values)
  }
  complete <- zero_tolerance == 0
  payload <- list(
    row_ptr = row_ptr, column_index = column, values = values,
    diagonal = diag(cross_product), complete = complete)
  .blr_new_operator_resource(
    resource_id, "csr", marker_ids, alleles, genotype_coding,
    "declared", "declared", "cross_product",
    .blr_new_operator_storage_ref("csr_cross_product", payload),
    block_eigen = NULL,
    approximation = if (complete) "exact_declared_operator" else
      "thresholded_cross_product_approximation",
    provenance = provenance)
}

.blr_block_eigen_resource <- function(
    resource_id, cross_product, marker_ids, alleles,
    blocks = list(seq_along(marker_ids)), retained_ranks = NULL,
    genotype_coding = "effect_allele_count", provenance = list()) {
  dense_ref <- .blr_new_operator_storage_ref(
    "dense_cross_product", cross_product)
  .blr_validate_dense_storage(dense_ref, marker_ids)
  if (!is.list(blocks) || !length(blocks)) {
    stop("blocks must be a nonempty list of marker-index vectors.",
         call. = FALSE)
  }
  if (is.null(retained_ranks)) retained_ranks <- lengths(blocks)
  if (!is.numeric(retained_ranks) || length(retained_ranks) != length(blocks)) {
    stop("retained_ranks must have one value per block.", call. = FALSE)
  }
  block_data <- lapply(seq_along(blocks), function(index) {
    ids <- as.integer(blocks[[index]])
    rank <- as.integer(retained_ranks[[index]])
    if (!length(ids) || anyNA(ids) || anyDuplicated(ids) ||
        any(ids < 1L | ids > length(marker_ids)) || rank < 1L ||
        rank > length(ids)) {
      stop("blocks or retained_ranks are invalid.", call. = FALSE)
    }
    decomposition <- eigen(
      cross_product[ids, ids, drop = FALSE], symmetric = TRUE)
    keep <- seq_len(rank)
    values <- pmax(decomposition$values[keep], 0)
    list(
      block_id = paste0("block", index), marker_indices = ids,
      eigenvectors = decomposition$vectors[, keep, drop = FALSE],
      eigenvalues = values, retained_rank = rank)
  })
  truncated <- any(retained_ranks < lengths(blocks))
  type <- if (truncated) "retained_rank_block_eigen" else
    "full_rank_block_eigen"
  .blr_new_operator_resource(
    resource_id, type, marker_ids, alleles, genotype_coding,
    "declared", "declared", "cross_product",
    .blr_new_operator_storage_ref("block_eigen", NULL),
    block_eigen = list(blocks = block_data),
    approximation = if (truncated) "retained_rank_approximation" else
      "exact_declared_block_diagonal_operator",
    provenance = provenance)
}

.blr_operator_apply <- function(resource, vector) {
  .blr_validate_operator_resource(resource)
  if (!is.numeric(vector) || length(vector) != length(resource$marker_ids) ||
      anyNA(vector) || any(!is.finite(vector))) {
    stop("operator vector must be finite and follow resource marker order.",
         call. = FALSE)
  }
  type <- resource$operator_type
  if (type %in% c("dense", "dense_cross_product")) {
    return(as.numeric(.blr_operator_storage_payload(resource$storage) %*%
                        vector))
  }
  if (identical(type, "csr")) {
    csr <- .blr_operator_storage_payload(resource$storage)
    out <- numeric(length(vector))
    for (row in seq_along(vector)) {
      first <- csr$row_ptr[[row]] + 1L
      last <- csr$row_ptr[[row + 1L]]
      if (first <= last) {
        out[[row]] <- sum(csr$values[first:last] *
                            vector[csr$column_index[first:last]])
      }
    }
    return(out)
  }
  if (type %in% c("full_rank_block_eigen",
                  "retained_rank_block_eigen")) {
    out <- numeric(length(vector))
    for (block in resource$block_eigen$blocks) {
      ids <- block$marker_indices
      coordinates <- crossprod(block$eigenvectors, vector[ids])
      out[ids] <- out[ids] + as.numeric(
        block$eigenvectors %*% (block$eigenvalues * coordinates))
    }
    return(out)
  }
  stop("This operator resource does not expose a cross-product vector action.",
       call. = FALSE)
}

.blr_operator_matrix <- function(resource) {
  m <- length(resource$marker_ids)
  vapply(seq_len(m), function(column) {
    .blr_operator_apply(resource, diag(m)[, column])
  }, numeric(m))
}

.blr_provider_global_score <- function(provider, resource, global_map) {
  .blr_validate_phase2_provider_alignment(provider, resource, global_map)
  score <- provider$sufficient_statistics$score
  expected <- c(length(resource$marker_ids), length(provider$trait_ids))
  if (!is.matrix(score) || !is.numeric(score) ||
      !identical(dim(score), expected) || any(!is.finite(score)) ||
      !identical(dimnames(score), list(
        marker = resource$marker_ids, trait = provider$trait_ids))) {
    stop("summary score must be finite and follow provider marker and trait order.",
         call. = FALSE)
  }
  out <- matrix(0, length(global_map$marker_ids), length(provider$trait_ids),
                dimnames = list(marker = global_map$marker_ids,
                                trait = provider$trait_ids))
  out[provider$local_to_global, ] <- score
  out
}

.blr_provider_residual_score <- function(collection, provider_id,
                                         global_effects) {
  .blr_character_scalar(provider_id, "provider_id")
  if (!provider_id %in% names(collection$providers)) {
    stop("unknown provider_id.", call. = FALSE)
  }
  provider <- collection$providers[[provider_id]]
  resource <- .blr_provider_resource(collection, provider_id)
  marker_ids <- collection$global_marker_map$marker_ids
  expected <- c(length(marker_ids), length(provider$trait_ids))
  if (!is.matrix(global_effects) || !is.numeric(global_effects) ||
      !identical(dim(global_effects), expected) ||
      any(!is.finite(global_effects)) ||
      !identical(dimnames(global_effects), list(
        marker = marker_ids, trait = provider$trait_ids))) {
    stop("global_effects must be finite and follow global marker and provider trait order.",
         call. = FALSE)
  }
  score <- provider$sufficient_statistics$score
  local_effects <- global_effects[provider$local_to_global, , drop = FALSE]
  fitted <- vapply(seq_along(provider$trait_ids), function(trait) {
    .blr_operator_apply(resource, local_effects[, trait])
  }, numeric(length(resource$marker_ids)))
  dimnames(fitted) <- list(marker = resource$marker_ids,
                           trait = provider$trait_ids)
  if (!is.matrix(score) || !identical(dim(score), dim(fitted)) ||
      !identical(dimnames(score), dimnames(fitted)) ||
      any(!is.finite(score))) {
    stop("provider score must be finite and follow resource marker and trait order.",
         call. = FALSE)
  }
  score - fitted
}
