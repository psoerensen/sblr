# Standalone summary-provider and likelihood-operator reference implementation.
#
# This research code deliberately does not depend on production sblr classes.
# Operators are on the genotype cross-product scale.  Provider likelihoods are
# evaluated separately and then added only after metadata compatibility and
# independence checks.

mt_assert_psd <- function(x, name = deparse(substitute(x)), tolerance = 1e-10) {
  if (!is.matrix(x) || nrow(x) != ncol(x) || any(!is.finite(x))) {
    stop(name, " must be a finite square matrix.", call. = FALSE)
  }
  if (!isTRUE(all.equal(x, t(x), tolerance = tolerance))) {
    stop(name, " must be symmetric.", call. = FALSE)
  }
  values <- eigen((x + t(x)) / 2, symmetric = TRUE, only.values = TRUE)$values
  if (min(values) < -tolerance * max(1, max(abs(values)))) {
    stop(name, " must be positive semidefinite.", call. = FALSE)
  }
  invisible(values)
}

mt_dense_operator <- function(cross_product) {
  mt_assert_psd(cross_product, "cross_product")
  structure(
    list(kind = "dense", dimension = nrow(cross_product),
         cross_product = (cross_product + t(cross_product)) / 2),
    class = "mt_summary_operator"
  )
}

mt_csr_operator <- function(cross_product, zero_tolerance = 0) {
  mt_assert_psd(cross_product, "cross_product")
  if (length(zero_tolerance) != 1L || !is.finite(zero_tolerance) ||
      zero_tolerance < 0) {
    stop("zero_tolerance must be one finite nonnegative number.",
         call. = FALSE)
  }
  sparse_matrix <- cross_product
  sparse_matrix[abs(sparse_matrix) <= zero_tolerance] <- 0
  mt_assert_psd(sparse_matrix, "thresholded CSR cross_product")
  p <- nrow(cross_product)
  row_ptr <- integer(p + 1L)
  column <- integer(0)
  value <- numeric(0)
  for (i in seq_len(p)) {
    keep <- which(sparse_matrix[i, ] != 0)
    column <- c(column, keep)
    value <- c(value, sparse_matrix[i, keep])
    row_ptr[i + 1L] <- length(column)
  }
  structure(
    list(kind = "csr", dimension = p, row_ptr = row_ptr,
         column = column, value = value),
    class = "mt_summary_operator"
  )
}

mt_eigen_operator <- function(cross_product, blocks = list(seq_len(nrow(cross_product))),
                              retained_rank = NULL, tolerance = 1e-10) {
  mt_assert_psd(cross_product, "cross_product", tolerance)
  p <- nrow(cross_product)
  if (!is.list(blocks) || !length(blocks) || any(lengths(blocks) == 0L)) {
    stop("blocks must be a nonempty list of nonempty local-index vectors.",
         call. = FALSE)
  }
  blocks <- lapply(blocks, function(index) {
    if (!is.numeric(index) || any(!is.finite(index)) ||
        any(index != as.integer(index))) {
      stop("Every block index must be an integer.", call. = FALSE)
    }
    as.integer(index)
  })
  flat <- unlist(blocks, use.names = FALSE)
  if (!identical(sort(flat), seq_len(p)) || anyDuplicated(flat)) {
    stop("blocks must form a disjoint complete partition of local markers.",
         call. = FALSE)
  }
  if (is.null(retained_rank)) retained_rank <- lengths(blocks)
  if (length(retained_rank) == 1L) retained_rank <- rep(retained_rank, length(blocks))
  if (length(retained_rank) != length(blocks) ||
      any(!is.finite(retained_rank)) ||
      any(retained_rank != as.integer(retained_rank)) ||
      any(retained_rank < 0L | retained_rank > lengths(blocks))) {
    stop("retained_rank must give a valid integer rank for every block.",
         call. = FALSE)
  }

  representation <- Map(function(index, rank) {
    block <- (cross_product[index, index, drop = FALSE] +
                t(cross_product[index, index, drop = FALSE])) / 2
    decomposition <- eigen(block, symmetric = TRUE)
    if (min(decomposition$values) <
        -tolerance * max(1, max(abs(decomposition$values)))) {
      stop("A block cross-product is not positive semidefinite.", call. = FALSE)
    }
    values <- pmax(decomposition$values, 0)
    keep <- if (rank == 0L) integer(0) else seq_len(rank)
    list(index = index,
         Q = decomposition$vectors[, keep, drop = FALSE],
         lambda = values[keep],
         full_dimension = length(index), retained_rank = rank)
  }, blocks, as.integer(retained_rank))

  structure(
    list(kind = "block_eigen", dimension = p, blocks = representation,
         exact_for_supplied_matrix =
           all(as.integer(retained_rank) == lengths(blocks)) &&
           max(abs(cross_product -
                     mt_block_diagonal_projection(cross_product, blocks))) <= tolerance),
    class = "mt_summary_operator"
  )
}

mt_block_diagonal_projection <- function(x, blocks) {
  out <- matrix(0, nrow(x), ncol(x))
  for (index in blocks) out[index, index] <- x[index, index, drop = FALSE]
  out
}

mt_operator_apply <- function(operator, vector) {
  if (!inherits(operator, "mt_summary_operator")) {
    stop("operator must inherit from mt_summary_operator.", call. = FALSE)
  }
  if (!is.numeric(vector) || length(vector) != operator$dimension ||
      any(!is.finite(vector))) {
    stop("vector must be finite and match the operator dimension.",
         call. = FALSE)
  }
  if (operator$kind == "dense") {
    return(as.numeric(operator$cross_product %*% vector))
  }
  if (operator$kind == "csr") {
    out <- numeric(operator$dimension)
    for (i in seq_len(operator$dimension)) {
      first <- operator$row_ptr[i] + 1L
      last <- operator$row_ptr[i + 1L]
      if (first <= last) {
        out[i] <- sum(operator$value[first:last] *
                        vector[operator$column[first:last]])
      }
    }
    return(out)
  }
  if (operator$kind == "block_eigen") {
    out <- numeric(operator$dimension)
    for (block in operator$blocks) {
      if (block$retained_rank > 0L) {
        local <- vector[block$index]
        out[block$index] <- as.numeric(
          block$Q %*% (block$lambda * as.numeric(crossprod(block$Q, local)))
        )
      }
    }
    return(out)
  }
  stop("Unknown summary-operator kind: ", operator$kind, call. = FALSE)
}

mt_operator_matrix <- function(operator) {
  basis <- diag(operator$dimension)
  vapply(seq_len(operator$dimension), function(j) {
    mt_operator_apply(operator, basis[, j])
  }, numeric(operator$dimension))
}

mt_operator_quadratic <- function(operator, vector) {
  sum(vector * mt_operator_apply(operator, vector))
}

mt_summary_provider <- function(id, trait, score, sample_size, markers,
                                global_markers, operator, effect_allele,
                                global_effect_allele,
                                genotype_coding = "centered_dosage",
                                genotype_scale = "cross_product",
                                phenotype_scale = "common_effect_scale",
                                population = "declared",
                                error_group = NA_character_) {
  if (length(id) != 1L || is.na(id) || !nzchar(id)) {
    stop("id must be one nonempty string.", call. = FALSE)
  }
  if (length(trait) != 1L || !is.finite(trait) || trait != as.integer(trait) ||
      trait < 1L) {
    stop("trait must be one positive integer.", call. = FALSE)
  }
  if (length(sample_size) != 1L || !is.finite(sample_size) || sample_size <= 0) {
    stop("sample_size must be one finite positive number.", call. = FALSE)
  }
  if (!inherits(operator, "mt_summary_operator")) {
    stop("operator must inherit from mt_summary_operator.", call. = FALSE)
  }
  if (!is.character(markers) || anyNA(markers) || any(!nzchar(markers)) ||
      anyDuplicated(markers)) {
    stop("markers must contain unique nonempty marker identifiers.",
         call. = FALSE)
  }
  if (!is.character(global_markers) || anyNA(global_markers) ||
      any(!nzchar(global_markers)) || anyDuplicated(global_markers)) {
    stop("global_markers must contain unique nonempty marker identifiers.",
         call. = FALSE)
  }
  map <- match(markers, global_markers)
  if (anyNA(map)) stop("Every provider marker must occur in global_markers.",
                       call. = FALSE)
  if (!is.numeric(score) || length(score) != length(markers) ||
      any(!is.finite(score)) || operator$dimension != length(markers)) {
    stop("score and operator dimensions must match provider markers.",
         call. = FALSE)
  }
  if (length(effect_allele) != length(markers) || anyNA(effect_allele) ||
      length(global_effect_allele) != length(global_markers) ||
      anyNA(global_effect_allele) ||
      !identical(as.character(effect_allele),
                 as.character(global_effect_allele[map]))) {
    stop("Provider alleles must be aligned to the global effect allele.",
         call. = FALSE)
  }
  metadata <- c(genotype_coding = genotype_coding,
                genotype_scale = genotype_scale,
                phenotype_scale = phenotype_scale)
  if (anyNA(metadata) || any(!nzchar(metadata))) {
    stop("Coding and scale metadata must be declared.", call. = FALSE)
  }
  if (!identical(genotype_scale, "cross_product")) {
    stop("This prototype requires operators on the cross-product scale.",
         call. = FALSE)
  }
  structure(
    list(id = id, trait = as.integer(trait), score = as.numeric(score),
         sample_size = sample_size, markers = markers,
         global_markers = global_markers, map = map, operator = operator,
         effect_allele = as.character(effect_allele),
         global_effect_allele = as.character(global_effect_allele),
         genotype_coding = genotype_coding,
         genotype_scale = genotype_scale,
         phenotype_scale = phenotype_scale, population = population,
         error_group = error_group),
    class = "mt_summary_provider"
  )
}

mt_validate_provider_collection <- function(providers, require_independent = TRUE) {
  if (!is.list(providers) || !length(providers) ||
      !all(vapply(providers, inherits, logical(1), "mt_summary_provider"))) {
    stop("providers must be a nonempty list of mt_summary_provider objects.",
         call. = FALSE)
  }
  ids <- vapply(providers, `[[`, character(1), "id")
  if (anyDuplicated(ids)) stop("Provider identifiers must be unique.",
                               call. = FALSE)
  reference_markers <- providers[[1]]$global_markers
  if (!all(vapply(providers, function(x) identical(x$global_markers,
                                                   reference_markers),
                  logical(1)))) {
    stop("Providers must use the same ordered global marker universe.",
         call. = FALSE)
  }
  reference_alleles <- providers[[1]]$global_effect_allele
  if (!all(vapply(providers, function(x) identical(x$global_effect_allele,
                                                   reference_alleles),
                  logical(1)))) {
    stop("Providers must use the same ordered global effect alleles.",
         call. = FALSE)
  }
  for (trait in unique(vapply(providers, `[[`, integer(1), "trait"))) {
    same_trait <- providers[vapply(providers, function(x) x$trait == trait,
                                   logical(1))]
    fields <- c("genotype_coding", "genotype_scale", "phenotype_scale")
    for (field in fields) {
      values <- vapply(same_trait, `[[`, character(1), field)
      if (length(unique(values)) != 1L) {
        stop("Providers for one shared effect have incompatible ", field, ".",
             call. = FALSE)
      }
    }
  }
  if (isTRUE(require_independent)) {
    groups <- vapply(providers, function(x) {
      if (length(x$error_group) != 1L || is.na(x$error_group) ||
          !nzchar(x$error_group)) paste0("independent:", x$id) else x$error_group
    }, character(1))
    duplicated_group <- duplicated(groups) | duplicated(groups, fromLast = TRUE)
    if (any(duplicated_group)) {
      stop("Providers sharing an error_group require an overlap-aware likelihood.",
           call. = FALSE)
    }
  }
  invisible(TRUE)
}

mt_provider_residual_score <- function(provider, effect_global) {
  if (!inherits(provider, "mt_summary_provider")) {
    stop("provider must be an mt_summary_provider.", call. = FALSE)
  }
  if (!is.numeric(effect_global) ||
      length(effect_global) != length(provider$global_markers) ||
      any(!is.finite(effect_global))) {
    stop("effect_global must be finite and match the global marker universe.",
         call. = FALSE)
  }
  local <- effect_global[provider$map]
  provider$score - mt_operator_apply(provider$operator, local)
}

mt_provider_loglik <- function(providers, effects) {
  mt_validate_provider_collection(providers, require_independent = TRUE)
  if (!is.matrix(effects) || nrow(effects) != length(providers[[1]]$global_markers) ||
      ncol(effects) < max(vapply(providers, `[[`, integer(1), "trait")) ||
      any(!is.finite(effects))) {
    stop("effects must be a finite global-marker by trait matrix.",
         call. = FALSE)
  }
  contributions <- vapply(providers, function(provider) {
    local <- effects[provider$map, provider$trait]
    sum(provider$score * local) -
      0.5 * mt_operator_quadratic(provider$operator, local)
  }, numeric(1))
  structure(sum(contributions), provider_contributions = contributions)
}

mt_combine_same_trait_providers <- function(providers, id = "combined") {
  mt_validate_provider_collection(providers, require_independent = TRUE)
  traits <- vapply(providers, `[[`, integer(1), "trait")
  if (length(unique(traits)) != 1L) {
    stop("Only providers for one shared trait effect can be combined.",
         call. = FALSE)
  }
  maps <- lapply(providers, `[[`, "map")
  if (!all(vapply(maps, identical, logical(1), maps[[1]]))) {
    stop("Explicit operator combination requires identical local marker maps.",
         call. = FALSE)
  }
  matrices <- lapply(providers, function(x) mt_operator_matrix(x$operator))
  first <- providers[[1]]
  mt_summary_provider(
    id = id, trait = first$trait,
    score = Reduce(`+`, lapply(providers, `[[`, "score")),
    sample_size = sum(vapply(providers, `[[`, numeric(1), "sample_size")),
    markers = first$markers, global_markers = first$global_markers,
    operator = mt_dense_operator(Reduce(`+`, matrices)),
    effect_allele = first$effect_allele,
    global_effect_allele = first$global_effect_allele,
    genotype_coding = first$genotype_coding,
    genotype_scale = first$genotype_scale,
    phenotype_scale = first$phenotype_scale,
    population = "combined-independent-providers",
    error_group = NA_character_
  )
}
