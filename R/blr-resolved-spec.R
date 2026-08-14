.blr_exact_names <- function(x, what, allowed = NULL) {
  if (!is.list(x) || is.data.frame(x)) {
    stop(what, " must be a list.", call. = FALSE)
  }
  if (!length(x)) return(invisible(character()))
  nm <- names(x)
  if (is.null(nm) || length(nm) != length(x) || anyNA(nm) ||
      any(!nzchar(nm)) || anyDuplicated(nm)) {
    stop(what, " must have unique, nonempty, non-NA names.", call. = FALSE)
  }
  if (!is.null(allowed)) {
    unknown <- setdiff(nm, allowed)
    if (length(unknown)) {
      stop("Unsupported argument(s). Unknown argument(s); unused argument(s) in ",
           what, ": ",
           paste(unknown, collapse = ", "), ".", call. = FALSE)
    }
  }
  invisible(nm)
}

.blr_exact_fields <- function(x, what, fields) {
  names <- .blr_exact_names(x, what, fields)
  missing <- setdiff(fields, names)
  if (length(missing)) {
    stop(what, " is missing required field(s): ",
         paste(missing, collapse = ", "), ".", call. = FALSE)
  }
  invisible(names)
}

.blr_capture_forwarded_args <- function(args, accepted, aliases = character(),
                                        what = "forwarded arguments") {
  .blr_exact_names(args, what)
  if (!is.character(accepted) || anyNA(accepted) || any(!nzchar(accepted)) ||
      anyDuplicated(accepted)) {
    stop("accepted must contain unique, nonempty argument names.", call. = FALSE)
  }
  if (length(aliases)) {
    if (is.null(names(aliases)) || anyNA(names(aliases)) ||
        any(!nzchar(names(aliases))) || anyDuplicated(names(aliases)) ||
        anyNA(aliases) || any(!nzchar(aliases))) {
      stop("aliases must be a uniquely named character vector.", call. = FALSE)
    }
    hit <- intersect(names(args), names(aliases))
    if (length(hit)) {
      target <- unname(aliases[hit])
      collision <- target %in% names(args) | anyDuplicated(target)
      if (any(collision)) {
        stop("Alias resolution would create duplicate argument(s): ",
             paste(target[collision], collapse = ", "), ".", call. = FALSE)
      }
      names(args)[match(hit, names(args))] <- target
      attr(args, "migration_actions") <- stats::setNames(target, hit)
    }
  }
  .blr_exact_names(args, what, accepted)
  args
}

.blr_validate_exact_public_call <- function(call, definition,
                                            what = "public BLR call") {
  supplied <- as.list(call)[-1L]
  supplied_names <- names(supplied)
  if (is.null(supplied_names)) return(invisible(TRUE))
  named <- !is.na(supplied_names) & nzchar(supplied_names)
  duplicated_names <- unique(supplied_names[named][duplicated(
    supplied_names[named])])
  if (length(duplicated_names)) {
    stop(what, " has duplicated argument `", duplicated_names[[1L]], "`.",
         call. = FALSE)
  }
  if (anyNA(supplied_names)) {
    stop(what, " must use unique, nonempty argument names.", call. = FALSE)
  }
  formal_names <- setdiff(names(formals(definition)), "...")
  abbreviated <- supplied_names[named & !supplied_names %in% formal_names &
    vapply(supplied_names, function(name) {
      nzchar(name) && any(startsWith(formal_names, name))
    }, logical(1))]
  if (length(abbreviated)) {
    stop(what, " requires exact argument names; abbreviated argument(s): ",
         paste(unique(abbreviated), collapse = ", "), ".", call. = FALSE)
  }
  invisible(TRUE)
}

.blr_scalar_whole <- function(x, what, lower = 0, upper = Inf) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x) ||
      x != floor(x) || x < lower || x > upper) {
    stop(what, " must be one finite integer-compatible value in [",
         lower, ", ", upper, "].", call. = FALSE)
  }
  as.numeric(x)
}

.blr_ids <- function(x, what) {
  if (!is.character(x) || !length(x) || anyNA(x) || any(!nzchar(x)) ||
      anyDuplicated(x)) {
    stop(what, " must contain unique, nonempty, non-NA character IDs.",
         call. = FALSE)
  }
  x
}

.blr_character_scalar <- function(x, what, allowed = NULL) {
  if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x) ||
      (!is.null(allowed) && !x %in% allowed)) {
    suffix <- if (is.null(allowed)) "." else paste0(
      "; expected one of: ", paste(allowed, collapse = ", "), ".")
    stop(what, " must be one nonempty character value", suffix,
         call. = FALSE)
  }
  x
}

.blr_validate_allele_table <- function(x, marker_ids, what,
                                       require_coding = FALSE) {
  required <- c("marker_id", "effect", "other",
                if (require_coding) "coding")
  if (!is.data.frame(x) || !all(required %in% names(x)) ||
      nrow(x) != length(marker_ids) ||
      !identical(as.character(x$marker_id), marker_ids)) {
    stop(what, " must have exactly one row in declared marker order.",
         call. = FALSE)
  }
  for (field in setdiff(required, "marker_id")) {
    value <- x[[field]]
    if (!is.character(value) || length(value) != length(marker_ids)) {
      stop(what, "$", field, " has an invalid allele contract.",
           call. = FALSE)
    }
    present <- !is.na(value)
    if (any(!nzchar(value[present]))) {
      stop(what, "$", field, " must be nonempty when present.",
           call. = FALSE)
    }
  }
  if (require_coding && anyNA(x$coding)) {
    stop(what, "$coding must be nonmissing.", call. = FALSE)
  }
  if (any(xor(is.na(x$effect), is.na(x$other)))) {
    stop(what, " must declare both alleles or mark both unavailable.",
         call. = FALSE)
  }
  known <- !is.na(x$effect)
  if (any(x$effect[known] == x$other[known])) {
    stop(what, " effect and other alleles must differ when declared.",
         call. = FALSE)
  }
  invisible(TRUE)
}

.blr_validate_spd_matrix <- function(x, trait_ids, what) {
  if (!is.matrix(x) || !is.numeric(x) ||
      !identical(dim(x), rep.int(length(trait_ids), 2L)) ||
      any(!is.finite(x)) ||
      !isTRUE(all.equal(x, t(x), tolerance = 0)) ||
      inherits(try(chol(x), silent = TRUE), "try-error")) {
    stop(what, " must be a finite symmetric positive-definite matrix.",
         call. = FALSE)
  }
  if (!is.null(dimnames(x)) &&
      !identical(dimnames(x), list(trait_ids, trait_ids))) {
    stop(what, " dimnames must follow declared trait order.", call. = FALSE)
  }
  invisible(TRUE)
}

.blr_validate_covariance_prior <- function(x, trait_ids, what) {
  if (is.null(x)) return(invisible(TRUE))
  .blr_exact_fields(x, what, c(
    "degrees_of_freedom", "scale", "fixed_value", "sampled"))
  if (!is.logical(x$sampled) || length(x$sampled) != 1L ||
      is.na(x$sampled)) {
    stop(what, "$sampled must be one nonmissing logical value.", call. = FALSE)
  }
  if (isTRUE(x$sampled)) {
    if (!is.numeric(x$degrees_of_freedom) ||
        length(x$degrees_of_freedom) != 1L ||
        is.na(x$degrees_of_freedom) || !is.finite(x$degrees_of_freedom) ||
        x$degrees_of_freedom <= length(trait_ids) - 1L) {
      stop(what, "$degrees_of_freedom does not define a proper inverse-Wishart prior.",
           call. = FALSE)
    }
    .blr_validate_spd_matrix(x$scale, trait_ids, paste0(what, "$scale"))
  } else if (!is.null(x$degrees_of_freedom) || !is.null(x$scale)) {
    stop(what, " cannot carry inverse-Wishart parameters when sampled is FALSE.",
         call. = FALSE)
  }
  if (!is.null(x$fixed_value)) {
    .blr_validate_spd_matrix(x$fixed_value, trait_ids,
                             paste0(what, "$fixed_value"))
  }
  invisible(TRUE)
}

.blr_validate_declared_prior <- function(x, what) {
  if (is.null(x)) return(invisible(TRUE))
  .blr_exact_names(x, what)
  validate_value <- function(value, path) {
    if (is.factor(value)) value <- as.character(value)
    if (is.null(value) || is.character(value) || is.logical(value)) {
      if (is.atomic(value) && anyNA(value)) {
        stop(path, " cannot contain missing values.", call. = FALSE)
      }
      return(invisible(TRUE))
    }
    if (is.numeric(value)) {
      if (anyNA(value) || any(!is.finite(value))) {
        stop(path, " must contain finite numeric values.", call. = FALSE)
      }
      return(invisible(TRUE))
    }
    if (is.list(value) && !is.data.frame(value)) {
      if (length(value) && !is.null(names(value))) .blr_exact_names(value, path)
      labels <- names(value) %||% as.character(seq_along(value))
      for (index in seq_along(value)) validate_value(
        value[[index]], paste0(path, "$", labels[[index]]))
      return(invisible(TRUE))
    }
    stop(path, " has an unsupported prior value type.", call. = FALSE)
  }
  for (name in names(x)) validate_value(x[[name]], paste0(what, "$", name))
  invisible(TRUE)
}

.blr_legacy_retained_indices <- function(nit, nthin) {
  nit <- as.integer(.blr_scalar_whole(nit, "sampling_iterations", 1,
                                      .Machine$integer.max))
  nthin <- as.integer(.blr_scalar_whole(nthin, "thin_interval", 1,
                                        .Machine$integer.max))
  seq.int(1L, nit, by = nthin)
}

.blr_target_retained_indices <- function(nit, nthin) {
  nit <- as.integer(.blr_scalar_whole(nit, "sampling_iterations", 1,
                                      .Machine$integer.max))
  nthin <- as.integer(.blr_scalar_whole(nthin, "thin_interval", 1,
                                        .Machine$integer.max))
  index <- seq_len(nit)
  index[index %% nthin == 0L]
}

.blr_valid_execution <- function(analysis_mode, execution_mode,
                                 parallelization) {
  allowed <- list(
    single_trait = c("none", "chains"),
    independent_traits = c("none", "chains", "traits", "trait_chains"),
    joint_multitrait = c("none", "chains")
  )
  if (!analysis_mode %in% names(allowed) ||
      !execution_mode %in% c("serial", "parallel") ||
      !parallelization %in% c("none", "chains", "traits", "trait_chains")) {
    return(FALSE)
  }
  if (identical(execution_mode, "serial")) {
    return(identical(parallelization, "none"))
  }
  parallelization %in% setdiff(allowed[[analysis_mode]], "none")
}

.blr_validate_task_seeds <- function(x, analysis_mode, trait_ids, chains) {
  valid_value <- is.numeric(x) && length(x) > 0L && !anyNA(x) &&
    all(is.finite(x)) && all(x == floor(x)) && all(x >= 0) &&
    all(x <= 4294967295)
  chain_ids <- paste0("chain", seq_len(chains))
  if (identical(analysis_mode, "independent_traits")) {
    shape <- identical(dim(x), as.integer(c(length(trait_ids), chains))) &&
      identical(dimnames(x), list(trait = trait_ids, chain = chain_ids))
  } else {
    shape <- is.null(dim(x)) && length(x) == chains &&
      identical(names(x), chain_ids)
  }
  if (!valid_value || !shape) {
    stop("mcmc$task_seeds has invalid uint32 values, shape, or names.",
         call. = FALSE)
  }
  invisible(TRUE)
}

.blr_validate_resource_provider_contract <- function(data) {
  .blr_ids(data$trait_ids, "data$trait_ids")
  markers <- .blr_ids(data$global_markers, "data$global_markers")
  .blr_validate_allele_table(data$global_alleles, markers,
                             "data$global_alleles", require_coding = TRUE)
  regimes <- c("common_sample", "independent_summary", "overlap_aware")
  .blr_character_scalar(data$likelihood_regime, "data$likelihood_regime",
                        regimes)
  if (!is.null(data$statistical_regions)) {
    regions <- data$statistical_regions
    if (!is.character(regions) || length(regions) != length(markers) ||
        anyNA(regions) || any(!nzchar(regions)) ||
        !identical(names(regions), markers)) {
      stop("data$statistical_regions must be a marker-named nonmissing region assignment.",
           call. = FALSE)
    }
  }
  .blr_exact_names(data$operator_resources, "data$operator_resources")
  .blr_exact_names(data$providers, "data$providers")
  .blr_exact_names(data$provider_maps, "data$provider_maps")
  if (!length(data$operator_resources) || !length(data$providers)) {
    stop("At least one operator resource and likelihood provider are required.",
         call. = FALSE)
  }
  resource_names <- names(data$operator_resources)
  resource_ids <- unname(vapply(data$operator_resources, function(resource) {
    .blr_exact_fields(resource, "operator resource", c(
      "resource_id", "operator_type", "marker_ids", "alleles",
      "genotype_coding", "centering", "standardization", "operator_scale",
      "storage", "block_eigen", "approximation", "provenance"))
    .blr_ids(resource$marker_ids, "operator resource marker_ids")
    if (!resource$operator_type %in% c(
        "bed", "dense", "dense_cross_product", "csr", "block_eigen",
        "full_rank_block_eigen", "retained_rank_block_eigen")) {
      stop("An operator resource has an unsupported operator_type.",
           call. = FALSE)
    }
    .blr_validate_allele_table(resource$alleles, resource$marker_ids,
                               "operator resource alleles")
    .blr_character_scalar(resource$genotype_coding,
                          "operator resource genotype_coding")
    as.character(resource$resource_id)
  }, character(1)))
  if (anyDuplicated(resource_ids) || !identical(resource_names, resource_ids)) {
    stop("Operator resources must be named by unique resource_id values.",
         call. = FALSE)
  }
  provider_ids <- unname(vapply(data$providers, function(provider) {
    .blr_exact_fields(provider, "likelihood provider", c(
      "provider_id", "trait_ids", "operator_resource_id", "local_to_global",
      "sufficient_statistics", "sample_size", "likelihood_regime",
      "residual_contract", "population", "effect_scale", "overlap_group",
      "provenance"))
    .blr_ids(provider$trait_ids, "provider$trait_ids")
    .blr_character_scalar(provider$likelihood_regime,
                          "provider$likelihood_regime", regimes)
    if (!all(provider$trait_ids %in% data$trait_ids)) {
      stop("A provider references an unknown trait ID.", call. = FALSE)
    }
    if (!provider$operator_resource_id %in% resource_ids) {
      stop("A provider references an unknown operator resource.", call. = FALSE)
    }
    resource <- data$operator_resources[[provider$operator_resource_id]]
    map <- provider$local_to_global
    if (!is.numeric(map) || length(map) != length(resource$marker_ids) ||
        anyNA(map) || any(map != floor(map)) || any(map < 1L) ||
        any(map > length(markers)) ||
        !identical(names(map), resource$marker_ids)) {
      stop("A provider local_to_global map is invalid or out of local order.",
           call. = FALSE)
    }
    if (!identical(unname(resource$marker_ids), unname(markers[map]))) {
      stop("A provider local_to_global map does not preserve marker identity.",
           call. = FALSE)
    }
    global_alleles <- data$global_alleles[map, , drop = FALSE]
    for (field in c("effect", "other")) {
      if (!identical(as.character(resource$alleles[[field]]),
                     as.character(global_alleles[[field]]))) {
        stop("Operator-resource alleles disagree with the global marker map.",
             call. = FALSE)
      }
    }
    sample_size <- provider$sample_size
    if (!is.numeric(sample_size) ||
        length(sample_size) != length(provider$trait_ids) ||
        anyNA(sample_size) || any(!is.finite(sample_size)) ||
        any(sample_size <= 0) ||
        !identical(names(sample_size), provider$trait_ids)) {
      stop("A provider sample_size must be a finite positive trait-named vector.",
           call. = FALSE)
    }
    as.character(provider$provider_id)
  }, character(1)))
  if (anyDuplicated(provider_ids) ||
      !identical(names(data$providers), provider_ids) ||
      !identical(names(data$provider_maps), provider_ids)) {
    stop("Providers and provider_maps must be named by unique provider IDs.",
         call. = FALSE)
  }
  for (id in provider_ids) {
    if (!identical(data$provider_maps[[id]],
                   data$providers[[id]]$local_to_global)) {
      stop("data$provider_maps does not match provider local_to_global maps.",
           call. = FALSE)
    }
  }
  represented_traits <- unique(unlist(lapply(
    data$providers, `[[`, "trait_ids"), use.names = FALSE))
  if (!setequal(represented_traits, data$trait_ids)) {
    stop("Every resolved trait must be represented by a likelihood provider.",
         call. = FALSE)
  }
  provider_regimes <- vapply(data$providers, `[[`, character(1),
                               "likelihood_regime")
  if (!all(provider_regimes == data$likelihood_regime)) {
    stop("Provider likelihood regimes must match data$likelihood_regime.",
         call. = FALSE)
  }
  if (identical(data$likelihood_regime, "independent_summary") &&
      any(vapply(data$providers, function(x) length(x$trait_ids) != 1L,
                 logical(1)))) {
    stop("independent_summary providers must each own exactly one trait.",
         call. = FALSE)
  }
  if (identical(data$analysis_mode, "joint_multitrait") &&
      identical(data$likelihood_regime, "common_sample") &&
      !any(vapply(data$providers, function(x) length(x$trait_ids) > 1L,
                  logical(1)))) {
    stop("A common-sample joint_multitrait likelihood requires a non-factorized multi-trait provider.",
         call. = FALSE)
  }
  invisible(TRUE)
}

validate_blr_resolved_spec <- function(spec) {
  .blr_exact_names(spec, "blr_resolved_spec")
  expected <- c("schema", "data", "model", "prior", "mcmc", "compute",
                "output")
  if (!identical(names(spec), expected)) {
    stop("blr_resolved_spec has invalid namespace names or order.", call. = FALSE)
  }
  .blr_exact_fields(spec$schema, "schema", c(
    "name", "version", "compatibility_id", "seed_contract_version",
    "retention_contract_version", "dimension_contract_version"))
  schema_versions <- c(
    version = spec$schema$version,
    seed = spec$schema$seed_contract_version,
    retention = spec$schema$retention_contract_version,
    dimension = spec$schema$dimension_contract_version)
  if (!identical(spec$schema$name, "blr_resolved_spec") ||
      !is.character(spec$schema$compatibility_id) ||
      length(spec$schema$compatibility_id) != 1L ||
      is.na(spec$schema$compatibility_id) ||
      !nzchar(spec$schema$compatibility_id) ||
      !is.numeric(schema_versions) || length(schema_versions) != 4L ||
      anyNA(schema_versions) ||
      any(!is.finite(schema_versions)) ||
      any(schema_versions != floor(schema_versions)) ||
      schema_versions[["version"]] != 1L ||
      !schema_versions[["seed"]] %in% c(0L, 1L) ||
      !schema_versions[["retention"]] %in% c(0L, 1L) ||
      schema_versions[["dimension"]] != 1L) {
    stop("Unsupported blr_resolved_spec schema or dimension version.",
         call. = FALSE)
  }
  .blr_exact_fields(spec$data, "data", c(
    "analysis_mode", "trait_ids", "global_markers", "global_alleles",
    "operator_resources", "providers", "provider_maps",
    "likelihood_regime", "statistical_regions"))
  if (!spec$data$analysis_mode %in% c(
      "single_trait", "independent_traits", "joint_multitrait")) {
    stop("data$analysis_mode is invalid.", call. = FALSE)
  }
  if (identical(spec$data$analysis_mode, "single_trait") &&
      length(spec$data$trait_ids) != 1L) {
    stop("single_trait requires exactly one trait ID.", call. = FALSE)
  }
  if (!identical(spec$data$analysis_mode, "single_trait") &&
      length(spec$data$trait_ids) < 2L) {
    stop(spec$data$analysis_mode, " requires at least two trait IDs.",
         call. = FALSE)
  }
  .blr_validate_resource_provider_contract(spec$data)
  .blr_exact_fields(spec$model, "model", c(
    "family", "state_space", "null_state_index",
    "effect_storage_convention", "probability_policy", "marker_scale_policy",
    "marker_covariance_policy", "residual_policy", "update_order_version"))
  .blr_character_scalar(spec$model$family, "model$family",
                        c("bayesc", "bayesr", "bayesrc"))
  .blr_character_scalar(spec$model$effect_storage_convention,
                        "model$effect_storage_convention",
                        c("realised", "base_latent", "scaled_latent"))
  .blr_character_scalar(spec$model$probability_policy,
                        "model$probability_policy",
                        c("global", "fixed_marker", "learned_logistic",
                          "group", "annotation_probit_stick"))
  .blr_character_scalar(spec$model$marker_scale_policy,
                        "model$marker_scale_policy",
                        c("unit", "component", "maf_s",
                          "component_maf_s", "annotation_log_variance"))
  .blr_character_scalar(spec$model$marker_covariance_policy,
                        "model$marker_covariance_policy",
                        c("traitwise_scalar", "global_matrix", "regional",
                          "covariance_templates"))
  .blr_character_scalar(spec$model$residual_policy,
                        "model$residual_policy",
                        c("scalar", "global_projected_legacy", "gctb_block",
                          "fixed_block", "fixed_full", "sampled_full",
                          "diagonal", "overlap_aware"))
  .blr_scalar_whole(spec$model$update_order_version,
                    "model$update_order_version", 1L, 1L)
  if (!is.character(spec$model$state_space) ||
      !length(spec$model$state_space) || anyNA(spec$model$state_space) ||
      any(!nzchar(spec$model$state_space)) ||
      anyDuplicated(spec$model$state_space)) {
    stop("model$state_space must contain ordered, unique state IDs.",
         call. = FALSE)
  }
  null_index <- .blr_scalar_whole(
    spec$model$null_state_index, "model$null_state_index", 1,
    length(spec$model$state_space))
  .blr_exact_fields(spec$prior, "prior", c(
    "probability", "component_multipliers", "marker_multipliers",
    "scalar_variance", "marker_covariance", "residual_covariance",
    "annotation"))
  if (!is.list(spec$prior$probability)) {
    stop("prior$probability must be an explicitly named list.", call. = FALSE)
  }
  .blr_exact_names(spec$prior$probability, "prior$probability")
  .blr_validate_declared_prior(spec$prior$probability, "prior$probability")
  probability_names <- c(names(spec$prior$probability),
    if (is.list(spec$prior$probability$values))
      names(spec$prior$probability$values) else character())
  if (any(c("pi", "pis", "pim") %in% probability_names)) {
    stop("prior$probability contains an ambiguous probability name.",
         call. = FALSE)
  }
  for (field in intersect(c("alpha", "beta", "shape", "strength"),
                          names(spec$prior$probability))) {
    value <- spec$prior$probability[[field]]
    if (!is.numeric(value) || length(value) != 1L || is.na(value) ||
        !is.finite(value) || value <= 0) {
      stop("prior$probability$", field,
           " must be one finite positive value.", call. = FALSE)
    }
  }
  probability_values <- spec$prior$probability$values
  if (is.list(probability_values)) {
    for (field in intersect(
        c("state_probability_initial", "component_probability_initial"),
        names(probability_values))) {
      value <- probability_values[[field]]
      if (!is.numeric(value) || length(value) != length(spec$model$state_space) ||
          anyNA(value) || any(!is.finite(value)) || any(value < 0) ||
          abs(sum(value) - 1) > 1e-10) {
        stop("prior$probability$values$", field,
             " must be a finite simplex in declared state order.",
             call. = FALSE)
      }
    }
    if ("annotation_stick_probability_initial" %in%
        names(probability_values)) {
      value <- probability_values$annotation_stick_probability_initial
      if (!is.numeric(value) ||
          length(value) != length(spec$model$state_space) - 1L ||
          anyNA(value) || any(!is.finite(value)) ||
          any(value <= 0 | value >= 1)) {
        stop("prior$probability$values$annotation_stick_probability_initial must contain one finite value in (0, 1) per nonterminal stick.",
             call. = FALSE)
      }
    }
  }
  if (!is.null(spec$prior$component_multipliers)) {
    gamma <- spec$prior$component_multipliers
    if (!is.numeric(gamma) || length(gamma) != length(spec$model$state_space) ||
        anyNA(gamma) || any(!is.finite(gamma)) || any(gamma < 0) ||
        gamma[[null_index]] != 0 ||
        !identical(names(gamma), spec$model$state_space)) {
      stop("prior$component_multipliers must be a state-named finite nonnegative vector with zero at the null state.",
           call. = FALSE)
    }
  }
  q <- spec$prior$marker_multipliers
  if (!is.null(q) && (!is.numeric(q) || length(q) != length(spec$data$global_markers) ||
      anyNA(q) || any(!is.finite(q)) || any(q <= 0) ||
      !identical(names(q), spec$data$global_markers))) {
    stop("prior$marker_multipliers must be a marker-named finite positive vector.",
         call. = FALSE)
  }
  .blr_validate_declared_prior(spec$prior$scalar_variance,
                               "prior$scalar_variance")
  if (!is.null(spec$prior$scalar_variance)) {
    scalar_values <- unlist(spec$prior$scalar_variance, recursive = TRUE,
                            use.names = TRUE)
    numeric_values <- suppressWarnings(as.numeric(scalar_values))
    names(numeric_values) <- names(scalar_values)
    positive <- grepl("shape|degrees_of_freedom", names(numeric_values))
    nonnegative <- grepl("scale", names(numeric_values))
    if (any(positive & (!is.finite(numeric_values) | numeric_values <= 0)) ||
        any(nonnegative & (!is.finite(numeric_values) | numeric_values < 0))) {
      stop("prior$scalar_variance shape and degrees-of-freedom values must be positive, and scale values must be nonnegative.",
           call. = FALSE)
    }
  }
  .blr_validate_covariance_prior(spec$prior$marker_covariance,
                                 spec$data$trait_ids,
                                 "prior$marker_covariance")
  .blr_validate_covariance_prior(spec$prior$residual_covariance,
                                 spec$data$trait_ids,
                                 "prior$residual_covariance")
  if (identical(spec$model$marker_covariance_policy, "global_matrix") &&
      is.null(spec$prior$marker_covariance)) {
    stop("global_matrix marker covariance requires prior$marker_covariance.",
         call. = FALSE)
  }
  if (identical(spec$model$residual_policy, "fixed_full") &&
      (is.null(spec$prior$residual_covariance) ||
       is.null(spec$prior$residual_covariance$fixed_value) ||
       isTRUE(spec$prior$residual_covariance$sampled))) {
    stop("fixed_full residual policy requires one fixed residual covariance.",
         call. = FALSE)
  }
  if (identical(spec$model$residual_policy, "sampled_full") &&
      (is.null(spec$prior$residual_covariance) ||
       !isTRUE(spec$prior$residual_covariance$sampled))) {
    stop("sampled_full residual policy requires a sampled covariance prior.",
         call. = FALSE)
  }
  .blr_validate_declared_prior(spec$prior$annotation, "prior$annotation")
  .blr_exact_fields(spec$mcmc, "mcmc", c(
    "burn_in_iterations", "sampling_iterations", "thin_interval",
    "retained_draws", "retained_transition_indices", "chains", "seed",
    "task_seeds", "update_flags"))
  burn <- .blr_scalar_whole(spec$mcmc$burn_in_iterations,
                            "mcmc$burn_in_iterations", 0,
                            .Machine$integer.max)
  nit <- .blr_scalar_whole(spec$mcmc$sampling_iterations,
                           "mcmc$sampling_iterations", 1,
                           .Machine$integer.max)
  thin <- .blr_scalar_whole(spec$mcmc$thin_interval,
                            "mcmc$thin_interval", 1,
                            .Machine$integer.max)
  chains <- as.integer(.blr_scalar_whole(spec$mcmc$chains, "mcmc$chains", 1,
                                         .Machine$integer.max))
  .blr_scalar_whole(spec$mcmc$seed, "mcmc$seed", 0, 4294967295)
  retention_version <- as.integer(spec$schema$retention_contract_version)
  expected_index <- if (retention_version == 1L) {
    .blr_target_retained_indices(nit, thin)
  } else if (retention_version == 0L) {
    .blr_legacy_retained_indices(nit, thin)
  } else {
    stop("Unsupported retention-contract version.", call. = FALSE)
  }
  retained_indices <- spec$mcmc$retained_transition_indices
  retained_draws <- spec$mcmc$retained_draws
  if (!is.numeric(retained_indices) || anyNA(retained_indices) ||
      any(!is.finite(retained_indices)) ||
      any(retained_indices != floor(retained_indices)) ||
      !is.numeric(retained_draws) || length(retained_draws) != 1L ||
      is.na(retained_draws) || !is.finite(retained_draws) ||
      retained_draws != floor(retained_draws) || retained_draws < 0 ||
      !identical(as.integer(retained_indices),
                 as.integer(expected_index)) ||
      !identical(as.integer(retained_draws),
                 as.integer(length(expected_index)))) {
    stop("MCMC retained indices or draw count do not match the declared retention contract.",
         call. = FALSE)
  }
  if (!length(expected_index) && length(spec$output$retained_parameters)) {
    stop("Retained parameters were requested but no transition is retained.",
         call. = FALSE)
  }
  .blr_exact_names(spec$mcmc$update_flags, "mcmc$update_flags")
  if (!all(vapply(spec$mcmc$update_flags, function(x) {
    is.logical(x) && length(x) == 1L && !is.na(x)
  }, logical(1)))) stop("mcmc$update_flags must contain logical scalars.", call. = FALSE)
  .blr_validate_task_seeds(spec$mcmc$task_seeds, spec$data$analysis_mode,
                           spec$data$trait_ids, chains)
  .blr_exact_fields(spec$compute, "compute", c(
    "execution_mode", "parallelization", "cores", "scheduler_version",
    "memory_limit_bytes", "operator_numerical_controls"))
  if (!.blr_valid_execution(spec$data$analysis_mode,
                            spec$compute$execution_mode,
                            spec$compute$parallelization)) {
    stop("Invalid analysis/execution combination.", call. = FALSE)
  }
  .blr_scalar_whole(spec$compute$cores, "compute$cores", 1,
                    .Machine$integer.max)
  .blr_scalar_whole(spec$compute$scheduler_version,
                    "compute$scheduler_version", 0,
                    .Machine$integer.max)
  memory_limit <- spec$compute$memory_limit_bytes
  if (!is.null(memory_limit) &&
      (!is.numeric(memory_limit) || length(memory_limit) != 1L ||
       is.na(memory_limit) || is.nan(memory_limit) ||
       memory_limit < 0 || identical(memory_limit, -Inf))) {
    stop("compute$memory_limit_bytes must be NULL, one finite nonnegative value, or positive Inf.",
         call. = FALSE)
  }
  .blr_exact_names(spec$compute$operator_numerical_controls,
                   "compute$operator_numerical_controls")
  .blr_exact_fields(spec$output, "output", c(
    "posterior_summaries", "retained_parameters", "effect_draw_policy",
    "state_draw_policy", "convergence_policy", "derived_quantities",
    "preserve_chains", "memory_estimate_bytes"))
  if (!is.logical(spec$output$posterior_summaries) ||
      length(spec$output$posterior_summaries) != 1L ||
      is.na(spec$output$posterior_summaries) ||
      !is.character(spec$output$retained_parameters) ||
      !is.character(spec$output$effect_draw_policy) ||
      length(spec$output$effect_draw_policy) != 1L ||
      !is.character(spec$output$state_draw_policy) ||
      length(spec$output$state_draw_policy) != 1L ||
      !is.character(spec$output$derived_quantities) ||
      !is.logical(spec$output$preserve_chains) ||
      length(spec$output$preserve_chains) != 1L ||
      is.na(spec$output$preserve_chains) ||
      !is.numeric(spec$output$memory_estimate_bytes) ||
      length(spec$output$memory_estimate_bytes) != 1L ||
      (!is.na(spec$output$memory_estimate_bytes) &&
       (!is.finite(spec$output$memory_estimate_bytes) ||
        spec$output$memory_estimate_bytes < 0))) {
    stop("output contains an invalid policy or memory estimate.", call. = FALSE)
  }
  .blr_exact_fields(spec$output$convergence_policy,
                    "output$convergence_policy", c("mode", "quantities"))
  .blr_character_scalar(spec$output$convergence_policy$mode,
                        "output$convergence_policy$mode",
                        c("legacy_wrapper", "none", "auto", "core",
                          "extended"))
  if (!is.character(spec$output$convergence_policy$quantities) ||
      anyNA(spec$output$convergence_policy$quantities) ||
      any(!nzchar(spec$output$convergence_policy$quantities)) ||
      anyDuplicated(spec$output$convergence_policy$quantities)) {
    stop("output$convergence_policy$quantities must contain unique nonmissing names.",
         call. = FALSE)
  }
  invisible(TRUE)
}

new_blr_resolved_spec <- function(schema, data, model, prior, mcmc, compute,
                                  output) {
  spec <- list(schema = schema, data = data, model = model, prior = prior,
               mcmc = mcmc, compute = compute, output = output)
  validate_blr_resolved_spec(spec)
  class(spec) <- c("blr_resolved_spec_v1", "blr_resolved_spec", "list")
  spec
}

.blr_legacy_task_seed_table <- function(chain, trait_ids) {
  seeds <- as.numeric(.blr_st_task_seeds(chain, length(trait_ids)))
  chain_ids <- paste0("chain", seq_len(chain$nchains))
  if (length(trait_ids) == 1L) {
    names(seeds) <- chain_ids
    return(seeds)
  }
  matrix(seeds, nrow = length(trait_ids), ncol = chain$nchains, byrow = TRUE,
         dimnames = list(trait = trait_ids, chain = chain_ids))
}

.blr_legacy_resource <- function(operator, marker_ids) {
  type <- switch(operator, packed_bed = "bed", csr = "csr",
                 block_eigen = "block_eigen", operator)
  list(
    resource_id = paste0(type, "_legacy_resource"), operator_type = type,
    marker_ids = marker_ids,
    alleles = data.frame(marker_id = marker_ids, effect = NA_character_,
                         other = NA_character_, stringsAsFactors = FALSE),
    genotype_coding = "legacy_wrapper_declared",
    centering = "legacy_wrapper_declared",
    standardization = "legacy_wrapper_declared",
    operator_scale = if (operator == "packed_bed") "individual_genotypes" else
      "cross_product",
    storage = list(kind = "legacy_immutable_view", payload = NULL),
    block_eigen = if (operator == "block_eigen")
      list(contract = "legacy_provider_operator") else NULL,
    approximation = if (operator == "block_eigen") "declared_approximation" else
      "legacy_route_contract",
    provenance = list(source = "maintained Phase 1 compatibility wrapper")
  )
}

resolve_blr_spec_from_wrapper <- function(model, operator, trait_ids, marker_ids,
                                          chain, sample_sizes,
                                          probability_policy = "global",
                                          marker_scale_policy = "unit",
                                          residual_policy = "scalar",
                                          component_multipliers = NULL,
                                          probability_prior = NULL,
                                          scalar_variance_prior = NULL,
                                          annotation_prior = NULL,
                                          update_flags = list(marker_effects = TRUE,
                                                              residual_variance = TRUE,
                                                              probability = TRUE),
                                          numerical_controls = list(),
                                          migration_actions = character()) {
  trait_ids <- .blr_ids(trait_ids, "trait_ids")
  marker_ids <- .blr_ids(marker_ids, "marker_ids")
  if (!is.numeric(sample_sizes) ||
      !length(sample_sizes) %in% c(1L, length(trait_ids)) ||
      anyNA(sample_sizes) || any(!is.finite(sample_sizes)) ||
      any(sample_sizes <= 0)) {
    stop("sample_sizes must contain one finite positive value per trait or one shared value.",
         call. = FALSE)
  }
  sample_sizes <- rep_len(as.numeric(sample_sizes), length(trait_ids))
  names(sample_sizes) <- trait_ids
  analysis_mode <- if (length(trait_ids) == 1L) "single_trait" else
    "independent_traits"
  resource <- .blr_legacy_resource(operator, marker_ids)
  resource_id <- resource$resource_id
  resources <- stats::setNames(list(resource), resource_id)
  providers <- lapply(seq_along(trait_ids), function(index) {
    id <- paste0("provider_", index)
    list(
      provider_id = id, trait_ids = trait_ids[index],
      operator_resource_id = resource_id,
      local_to_global = stats::setNames(seq_along(marker_ids), marker_ids),
      sufficient_statistics = list(source = "legacy_wrapper", payload = NULL),
      sample_size = sample_sizes[trait_ids[index]],
      likelihood_regime = if (operator == "packed_bed") "common_sample" else
        "independent_summary",
      residual_contract = "traitwise_scalar",
      population = NA_character_, effect_scale = "legacy_declared",
      overlap_group = NULL,
      provenance = list(source = "maintained Phase 1 compatibility wrapper")
    )
  })
  names(providers) <- vapply(providers, `[[`, character(1), "provider_id")
  retained <- .blr_legacy_retained_indices(chain$nit, chain$nthin)
  tasks <- length(trait_ids) * chain$nchains
  parallelization <- if (tasks <= 1L || chain$ncores <= 1L) "none" else if (
    length(trait_ids) > 1L && chain$nchains > 1L) "trait_chains" else if (
      length(trait_ids) > 1L) "traits" else "chains"
  execution <- if (parallelization == "none") "serial" else "parallel"
  family <- switch(model,
    bayesc = "bayesc", sbayesc = "bayesc",
    bayesr = "bayesr", sbayesr = "bayesr",
    bayesrc = "bayesrc", sbayesrc = "bayesrc",
    stop("Unsupported ST model for Phase 1 resolution: ", model,
         call. = FALSE))
  if (family == "bayesc") {
    states <- c("null", "active")
    component_multipliers <- NULL
  } else {
    component_multipliers <- component_multipliers %||% c(0, 0.01, 0.1, 1)
    if (!is.numeric(component_multipliers) ||
        length(component_multipliers) < 2L ||
        anyNA(component_multipliers) || any(!is.finite(component_multipliers)) ||
        component_multipliers[[1L]] != 0 ||
        any(component_multipliers[-1L] <= 0)) {
      stop("component_multipliers must define a finite null-plus-active BayesR mixture.",
           call. = FALSE)
    }
    states <- if (family == "bayesrc") {
      paste0("gamma_", formatC(component_multipliers, format = "f", digits = 2L))
    } else {
      paste0("component_", seq_along(component_multipliers) - 1L)
    }
    names(component_multipliers) <- states
  }
  spec <- new_blr_resolved_spec(
    schema = list(
      name = "blr_resolved_spec", version = 1L,
      compatibility_id = paste0("phase1-legacy-st-v1;seed=legacy_st_arithmetic_v1;retention=",
                                if (operator == "packed_bed") "st_bed_v1" else "st_scalar_v1",
                                if (length(migration_actions)) paste0(";aliases=", paste(migration_actions, collapse = ",")) else ""),
      seed_contract_version = 0L, retention_contract_version = 0L,
      dimension_contract_version = 1L),
    data = list(
      analysis_mode = analysis_mode, trait_ids = trait_ids,
      global_markers = marker_ids,
      global_alleles = data.frame(marker_id = marker_ids,
                                  effect = NA_character_, other = NA_character_,
                                  coding = "legacy_wrapper_declared",
                                  stringsAsFactors = FALSE),
      operator_resources = resources, providers = providers,
      provider_maps = lapply(providers, `[[`, "local_to_global"),
      likelihood_regime = if (operator == "packed_bed") "common_sample" else
        "independent_summary", statistical_regions = NULL),
    model = list(
      family = family, state_space = states, null_state_index = 1L,
      effect_storage_convention = "realised",
      probability_policy = probability_policy,
      marker_scale_policy = marker_scale_policy,
      marker_covariance_policy = "traitwise_scalar",
      residual_policy = residual_policy, update_order_version = 1L),
    prior = list(
      probability = list(contract = "legacy_wrapper_resolved",
                         values = probability_prior),
      component_multipliers = component_multipliers,
      marker_multipliers = stats::setNames(rep(1, length(marker_ids)), marker_ids),
      scalar_variance = list(contract = "legacy_wrapper_resolved",
                             values = scalar_variance_prior),
      marker_covariance = NULL, residual_covariance = NULL,
      annotation = annotation_prior),
    mcmc = list(
      burn_in_iterations = chain$nburn,
      sampling_iterations = chain$nit, thin_interval = chain$nthin,
      retained_draws = length(retained),
      retained_transition_indices = retained, chains = chain$nchains,
      seed = chain$seed, task_seeds = .blr_legacy_task_seed_table(chain, trait_ids),
      update_flags = update_flags),
    compute = list(
      execution_mode = execution, parallelization = parallelization,
      cores = chain$ncores, scheduler_version = 0L,
      memory_limit_bytes = NULL,
      operator_numerical_controls = numerical_controls),
    output = list(
      posterior_summaries = TRUE, retained_parameters = character(),
      effect_draw_policy = "none", state_draw_policy = "none",
      convergence_policy = list(mode = "legacy_wrapper",
                                quantities = character()),
      derived_quantities = character(), preserve_chains = TRUE,
      memory_estimate_bytes = 0)
  )
  spec
}

.blr_complete_legacy_spec_from_fit <- function(spec, fit) {
  validate_blr_resolved_spec(spec)
  input <- fit$input %||% list()
  map_present <- function(mapping) {
    source <- names(mapping)
    hit <- source[source %in% names(input)]
    hit <- hit[!vapply(input[hit], is.null, logical(1))]
    values <- input[hit]
    names(values) <- unname(mapping[hit])
    values
  }
  probability_map <- c(
    pi_init = "non_null_probability_initial",
    pi_vb_init = "variance_calibration_non_null_probability_initial",
    pi_prior_mean = "non_null_probability_prior_mean",
    pi_prior_strength = "non_null_probability_prior_strength",
    pi_prior_a = "non_null_probability_prior_shape_active",
    pi_prior_b = "non_null_probability_prior_shape_null")
  probability_values <- map_present(probability_map)
  if (!is.null(input$pi)) {
    name <- switch(spec$model$family,
      bayesc = "state_probability_initial",
      bayesr = "component_probability_initial",
      bayesrc = if (length(input$pi) == length(spec$model$state_space))
        "component_probability_initial" else NULL)
    if (!is.null(name)) {
      probability_values[[name]] <- if (
        spec$model$family == "bayesc" && length(input$pi) == 1L) {
        stats::setNames(c(1 - as.numeric(input$pi), as.numeric(input$pi)),
                        spec$model$state_space)
      } else input$pi
    }
  }
  spec$prior$probability$values <- probability_values
  spec$prior$probability$source_arguments <- names(probability_map)[
    names(probability_map) %in% names(input)]
  scalar_map <- c(
    h2 = "heritability_calibration",
    nub = "marker_variance_prior_degrees_of_freedom",
    nue = "residual_variance_prior_degrees_of_freedom",
    B = "marker_variance_prior_scale",
    E = "residual_variance_prior_scale",
    ssb_prior = "marker_variance_prior_scale_by_trait",
    sse_prior = "residual_variance_prior_scale_by_trait")
  spec$prior$scalar_variance$values <- map_present(scalar_map)
  spec$prior$scalar_variance$source_arguments <- names(scalar_map)[
    names(scalar_map) %in% names(input)]
  annotation_map <- c(
    annotation_model = "model",
    annotation_policy = "probability_policy",
    effect_scale_policy = "marker_scale_policy",
    effect_scale = "marker_scale_contract",
    A = "annotation_matrix",
    pi_marker = "fixed_marker_inclusion_probabilities",
    vb_multiplier = "marker_variance_multipliers",
    use_pi_marker = "uses_fixed_marker_inclusion_probabilities",
    use_vb_multiplier = "uses_marker_variance_multipliers",
    group = "group_assignments",
    group_names = "group_ids",
    group_size = "group_sizes")
  spec$prior$annotation <- c(
    spec$prior$annotation %||% list(),
    map_present(annotation_map))

  sample_size <- input$n %||% input$N %||% NULL
  if (is.numeric(sample_size) && length(sample_size) %in%
      c(1L, length(spec$data$trait_ids)) && !anyNA(sample_size) &&
      all(is.finite(sample_size)) && all(sample_size > 0)) {
    sample_size <- rep_len(as.numeric(sample_size), length(spec$data$trait_ids))
    names(sample_size) <- spec$data$trait_ids
    for (provider_id in names(spec$data$providers)) {
      provider_traits <- spec$data$providers[[provider_id]]$trait_ids
      spec$data$providers[[provider_id]]$sample_size <- sample_size[provider_traits]
    }
  }
  estimated_bytes <- fit$memory_estimate$estimated_total_bytes %||% NA_real_
  spec$output$memory_estimate_bytes <- if (
    is.numeric(estimated_bytes) && length(estimated_bytes) == 1L &&
      !is.na(estimated_bytes) && is.finite(estimated_bytes) &&
      estimated_bytes >= 0) as.numeric(estimated_bytes) else NA_real_
  validate_blr_resolved_spec(spec)
  spec
}

.new_blr_resolved_spec <- new_blr_resolved_spec
.validate_blr_resolved_spec <- validate_blr_resolved_spec
.resolve_blr_spec_from_wrapper <- resolve_blr_spec_from_wrapper
