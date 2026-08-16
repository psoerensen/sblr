# Shared public resolution for corrected independent-summary Cheng MT-BayesC-Pi.

.mtblr_summary_descriptor <- function(x, what, required, optional = character()) {
  if (!is.list(x) || is.null(names(x)) || anyNA(names(x)) ||
      any(!nzchar(names(x))) || anyDuplicated(names(x))) {
    stop(what, " must be a uniquely named list.", call. = FALSE)
  }
  missing <- setdiff(required, names(x))
  if (length(missing)) {
    stop(what, " is missing: ", paste(missing, collapse = ", "), ".",
         call. = FALSE)
  }
  unknown <- setdiff(names(x), c(required, optional))
  if (length(unknown)) {
    stop(what, " contains unsupported fields: ",
         paste(unknown, collapse = ", "), ".", call. = FALSE)
  }
  x
}

.mtblr_summary_public_csr_payload <- function(csr, marker_ids, diagonal,
                                              approximation) {
  if (!is.list(csr)) stop("CSR data must be a named list.", call. = FALSE)
  if (all(c("row_ptr", "column_index", "values", "diagonal", "complete") %in%
          names(csr))) {
    return(list(
      row_ptr = as.numeric(csr$row_ptr),
      column_index = as.integer(csr$column_index),
      values = as.numeric(csr$values), diagonal = as.numeric(csr$diagonal),
      complete = isTRUE(csr$complete)))
  }
  required <- c("row_ptr", "col_idx", "values")
  if (length(setdiff(required, names(csr)))) {
    stop("CSR data must be returned by sparseLD_read_CSR() or follow the Phase 2 explicit CSR contract.",
         call. = FALSE)
  }
  m <- length(marker_ids)
  diagonal <- as.numeric(diagonal)
  if (length(diagonal) != m || any(!is.finite(diagonal)) ||
      any(diagonal <= 0)) {
    stop("Each public CSR resource requires one finite positive cross-product diagonal per marker.",
         call. = FALSE)
  }
  row_ptr <- as.numeric(csr$row_ptr)
  columns <- as.integer(csr$col_idx)
  values <- as.numeric(csr$values)
  if (length(row_ptr) != m + 1L || row_ptr[[1L]] != 0 ||
      any(!is.finite(row_ptr)) || any(row_ptr != floor(row_ptr)) ||
      any(diff(row_ptr) < 0) || tail(row_ptr, 1L) != length(values) ||
      length(columns) != length(values) || anyNA(columns) ||
      any(!is.finite(values))) {
    stop("CSR row pointers, columns, or values are invalid.", call. = FALSE)
  }
  index_base <- as.integer(csr$index_base %||% 1L)
  if (!index_base %in% c(0L, 1L)) {
    stop("CSR index_base must be zero or one.", call. = FALSE)
  }
  columns <- columns + if (index_base == 0L) 1L else 0L
  if (any(columns < 1L | columns > m)) {
    stop("CSR column indices are outside the resource marker range.",
         call. = FALSE)
  }
  correlation_scale <- identical(csr$value %||% NULL, "r") ||
    identical(csr$ld_normalization %||% NULL, "sqrt_xx") ||
    identical(csr$diag %||% NULL, "implicit_1")
  upper <- isTRUE(csr$upper_triangle)
  row_columns <- vector("list", m)
  row_values <- vector("list", m)
  add <- function(row, column, value) {
    row_columns[[row]] <<- c(row_columns[[row]], column)
    row_values[[row]] <<- c(row_values[[row]], value)
  }
  for (row in seq_len(m)) {
    first <- row_ptr[[row]] + 1L
    last <- row_ptr[[row + 1L]]
    if (first <= last) for (position in first:last) {
      column <- columns[[position]]
      if (column == row) next
      value <- values[[position]]
      if (correlation_scale) {
        value <- value * sqrt(diagonal[[row]] * diagonal[[column]])
      }
      add(row, column, value)
      if (upper) add(column, row, value)
    }
  }
  for (row in seq_len(m)) add(row, row, diagonal[[row]])
  for (row in seq_len(m)) {
    order <- order(row_columns[[row]], method = "radix")
    row_columns[[row]] <- row_columns[[row]][order]
    row_values[[row]] <- row_values[[row]][order]
  }
  list(
    row_ptr = c(0, cumsum(lengths(row_columns))),
    column_index = as.integer(unlist(row_columns, use.names = FALSE)),
    values = as.numeric(unlist(row_values, use.names = FALSE)),
    diagonal = diagonal,
    complete = identical(approximation, "exact_declared_operator"))
}

.mtblr_summary_public_resource <- function(descriptor, representation) {
  common <- c(
    "resource_id", "marker_ids", "alleles", "genotype_coding",
    "centering", "standardization", "operator_scale", "approximation",
    "provenance")
  if (identical(representation, "csr")) {
    descriptor <- .mtblr_summary_descriptor(
      descriptor, "CSR resource descriptor",
      c(common, "csr", "diagonal"))
    marker_ids <- .blr_ids(descriptor$marker_ids,
                           "CSR resource marker_ids")
    payload <- .mtblr_summary_public_csr_payload(
      descriptor$csr, marker_ids, descriptor$diagonal,
      descriptor$approximation)
    return(.blr_new_operator_resource(
      descriptor$resource_id, "csr", marker_ids, descriptor$alleles,
      descriptor$genotype_coding, descriptor$centering,
      descriptor$standardization, descriptor$operator_scale,
      .blr_new_operator_storage_ref("csr_cross_product", payload),
      block_eigen = NULL, approximation = descriptor$approximation,
      provenance = descriptor$provenance))
  }
  descriptor <- .mtblr_summary_descriptor(
    descriptor, "block-eigen resource descriptor", c(common, "blocks"))
  marker_ids <- .blr_ids(descriptor$marker_ids,
                         "block-eigen resource marker_ids")
  if (!is.list(descriptor$blocks) || !length(descriptor$blocks)) {
    stop("block-eigen resources require a nonempty blocks list.",
         call. = FALSE)
  }
  blocks <- lapply(seq_along(descriptor$blocks), function(index) {
    block <- .mtblr_summary_descriptor(
      descriptor$blocks[[index]], paste0("block-eigen block ", index),
      c("marker_ids", "eigenvectors", "eigenvalues"), "block_id")
    ids <- .blr_ids(block$marker_ids, "block-eigen block marker_ids")
    marker_indices <- match(ids, marker_ids)
    if (anyNA(marker_indices)) {
      stop("Every block marker must belong to its resource.", call. = FALSE)
    }
    vectors <- as.matrix(block$eigenvectors)
    values <- as.numeric(block$eigenvalues)
    list(
      block_id = block$block_id %||% paste0("block", index),
      marker_indices = as.integer(marker_indices), eigenvectors = vectors,
      eigenvalues = values, retained_rank = as.integer(ncol(vectors)))
  })
  retained <- any(vapply(blocks, function(block) {
    block$retained_rank < length(block$marker_indices)
  }, logical(1)))
  type <- if (retained) "retained_rank_block_eigen" else
    "full_rank_block_eigen"
  if (retained && identical(descriptor$approximation,
                            "exact_declared_block_diagonal_operator")) {
    stop("A retained-rank block-eigen resource must declare an approximation.",
         call. = FALSE)
  }
  .blr_new_operator_resource(
    descriptor$resource_id, type, marker_ids, descriptor$alleles,
    descriptor$genotype_coding, descriptor$centering,
    descriptor$standardization, descriptor$operator_scale,
    .blr_new_operator_storage_ref("block_eigen", NULL),
    block_eigen = list(blocks = blocks),
    approximation = descriptor$approximation,
    provenance = descriptor$provenance)
}

.mtblr_summary_public_collection <- function(
    providers, operator_resources, global_marker_ids, global_alleles,
    trait_ids, representation) {
  trait_ids <- .blr_ids(trait_ids, "trait_ids")
  global_marker_ids <- .blr_ids(global_marker_ids, "global_marker_ids")
  global_map <- .blr_new_global_marker_map(global_marker_ids, global_alleles)
  if (!is.list(providers) || !length(providers) ||
      !is.list(operator_resources) || !length(operator_resources)) {
    stop("providers and operator_resources must be nonempty lists.",
         call. = FALSE)
  }
  provider_required <- c(
    "provider_id", "trait_id", "operator_resource_id", "score",
    "sample_size", "residual_scale", "likelihood_regime", "effect_scale")
  provider_optional <- c("population", "overlap_group", "provenance")
  providers <- lapply(seq_along(providers), function(index) {
    provider <- .mtblr_summary_descriptor(
      providers[[index]], paste0("provider descriptor ", index),
      provider_required, provider_optional)
    if (!identical(provider$likelihood_regime, "independent_summary")) {
      stop("Every provider must declare likelihood_regime = 'independent_summary'.",
           call. = FALSE)
    }
    if (!is.null(provider$overlap_group)) {
      stop("Declared sample overlap requires an overlap-aware likelihood, which is not available.",
           call. = FALSE)
    }
    provider
  })
  resources <- lapply(operator_resources, .mtblr_summary_public_resource,
                      representation = representation)
  resource_ids <- vapply(resources, `[[`, character(1), "resource_id")
  if (anyDuplicated(resource_ids)) {
    stop("operator resource IDs must be unique.", call. = FALSE)
  }
  names(resources) <- resource_ids
  provider_objects <- lapply(providers, function(provider) {
    id <- .blr_character_scalar(provider$provider_id, "provider_id")
    trait <- .blr_character_scalar(provider$trait_id, "trait_id")
    if (!trait %in% trait_ids) stop("provider trait_id is undeclared.",
                                   call. = FALSE)
    resource_id <- .blr_character_scalar(
      provider$operator_resource_id, "operator_resource_id")
    resource <- resources[[resource_id]]
    if (is.null(resource)) stop("provider references an unknown resource.",
                                call. = FALSE)
    score <- provider$score
    if (!is.numeric(score) || length(score) != length(resource$marker_ids) ||
        any(!is.finite(score)) || is.null(names(score)) ||
        !identical(names(score), resource$marker_ids)) {
      stop("provider score must be a finite named vector in resource marker order.",
           call. = FALSE)
    }
    sample_size <- .blr_scalar_whole(
      provider$sample_size, "provider sample_size", 1, 2^53)
    residual_scale <- provider$residual_scale
    if (!is.numeric(residual_scale) || length(residual_scale) != 1L ||
        !is.finite(residual_scale) || residual_scale <= 0) {
      stop("provider residual_scale must be one finite positive value.",
           call. = FALSE)
    }
    local_to_global <- stats::setNames(
      match(resource$marker_ids, global_marker_ids), resource$marker_ids)
    if (anyNA(local_to_global)) {
      stop("Every provider marker must occur in global_marker_ids.",
           call. = FALSE)
    }
    .blr_new_likelihood_provider(
      id, trait, resource_id, local_to_global,
      sufficient_statistics = list(
        score = matrix(as.numeric(score), ncol = 1L,
                       dimnames = list(marker = resource$marker_ids,
                                       trait = trait)),
        residual_scale = as.numeric(residual_scale)),
      sample_size = stats::setNames(sample_size, trait),
      likelihood_regime = "independent_summary",
      residual_contract = "fixed_provider_residual_scale",
      population = provider$population %||% NA_character_,
      effect_scale = .blr_character_scalar(
        provider$effect_scale, "provider effect_scale"),
      overlap_group = NULL, provenance = provider$provenance %||% list())
  })
  provider_ids <- vapply(provider_objects, `[[`, character(1), "provider_id")
  if (anyDuplicated(provider_ids)) stop("provider IDs must be unique.",
                                        call. = FALSE)
  names(provider_objects) <- provider_ids
  .blr_new_provider_collection(
    global_map, resources, provider_objects,
    likelihood_regime = "independent_summary",
    analysis_mode = "joint_multitrait")
}

.blr_phase6b_promote_summary_raw <- function(raw, interface, operator) {
  validate_blr_raw_v2(raw)
  raw$input$schema$compatibility_id <- paste0(
    "phase6b-public-summary-cheng-mt-", operator,
    ";seed=unified_fnv_splitmix_v1;retention=postburn_divisible_v1")
  raw$diagnostics$qualification$status <- "publicly_supported"
  raw$diagnostics$qualification$implementation <-
    "general_t_independent_summary_cheng_mt_bayesc"
  raw$diagnostics$qualification$public_interface <- interface
  raw$schema$source_schema$name <- paste0(
    "general_t_summary_cheng_mt_bayesc_", operator)
  raw$schema$migration$status <- "public_phase6b"
  raw$schema$migration$legacy_mt_conversion <- FALSE
  validate_blr_raw_v2(raw)
  raw
}

.mtblr_summary_public_fit <- function(
    providers, operator_resources, global_marker_ids, global_alleles,
    trait_ids, representation, interface, method,
    initial_marker_covariance,
    marker_covariance_prior_degrees_of_freedom,
    marker_covariance_prior_scale,
    initial_activity_pattern_probability,
    activity_pattern_dirichlet_prior,
    nit, nburn, nthin, seed, nchains, ncores, chain_seeds,
    keep_chains, convergence, convergence_control, keep_traces,
    memory_limit_bytes) {
  if (!identical(method, "bayesc")) {
    stop(interface, " supports only method = 'bayesc' under the corrected Cheng model.",
         call. = FALSE)
  }
  convergence <- match.arg(convergence, c("core", "none"))
  if (!is.null(convergence_control)) {
    stop(interface, " requires convergence_control = NULL.", call. = FALSE)
  }
  keep_chains <- .blr_logical_scalar(keep_chains, "keep_chains")
  keep_traces <- .blr_logical_scalar(keep_traces, "keep_traces")
  if (identical(convergence, "none") && keep_traces) {
    stop("convergence = 'none' requires keep_traces = FALSE.", call. = FALSE)
  }
  collection <- .mtblr_summary_public_collection(
    providers, operator_resources, global_marker_ids, global_alleles,
    trait_ids, representation)
  raw <- .blr_cheng_mt_bayesc_summary_qualification(
    collection = collection, trait_ids = trait_ids,
    initial_marker_covariance = initial_marker_covariance,
    marker_covariance_prior_df =
      marker_covariance_prior_degrees_of_freedom,
    marker_covariance_prior_scale = marker_covariance_prior_scale,
    initial_activity_pattern_probability =
      initial_activity_pattern_probability,
    activity_pattern_dirichlet_prior = activity_pattern_dirichlet_prior,
    update_marker_covariance = TRUE,
    update_activity_pattern_probability = TRUE,
    burn_in_iterations = nburn, sampling_iterations = nit,
    thin_interval = nthin, chains = nchains, cores = ncores, seed = seed,
    chain_seeds = chain_seeds,
    keep_traces = keep_traces && identical(convergence, "core"),
    memory_limit_bytes = memory_limit_bytes)
  operator <- if (identical(representation, "csr")) "csr" else
    "block_eigen"
  raw <- .blr_phase6b_promote_summary_raw(raw, interface, operator)
  .blr_format_cheng_mt_raw_v2(
    raw, keep_chains = keep_chains, operator = operator)
}
