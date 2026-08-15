.blr_raw_axis_contract <- list(
  realised_effects = c("draw", "chain", "marker", "trait"),
  latent_effects = c("draw", "chain", "marker", "trait"),
  scaled_effects = c("draw", "chain", "marker", "trait"),
  independent_trait_states = c("draw", "chain", "marker", "trait"),
  joint_states = c("draw", "chain", "marker"),
  traitwise_activity = c("draw", "chain", "marker", "trait"),
  realised_effect_mean = c("marker", "trait"),
  latent_effect_mean = c("marker", "trait"),
  scaled_effect_mean = c("marker", "trait"),
  pips = c("marker", "trait"),
  traitwise_state_probabilities = c("marker", "trait", "state"),
  joint_state_probabilities = c("marker", "joint_state"),
  activity_pattern_probabilities = c("marker", "activity_pattern"),
  traitwise_component_assignment_probabilities = c("marker", "trait", "component"),
  joint_component_assignment_probabilities = c("marker", "component"),
  traitwise_probability_parameter_mean = c("trait", "state"),
  traitwise_component_probability_parameter_mean = c("trait", "component"),
  traitwise_probability_parameters = c("draw", "chain", "trait", "state"),
  joint_probability_parameters = c("draw", "chain", "joint_state"),
  activity_pattern_parameters = c("draw", "chain", "activity_pattern"),
  traitwise_component_probability_parameters = c("draw", "chain", "trait", "component"),
  joint_component_probability_parameters = c("draw", "chain", "component"),
  marker_variance = c("draw", "chain", "trait"),
  residual_variance = c("draw", "chain", "trait"),
  marker_covariance = c("draw", "chain", "trait_row", "trait_col"),
  residual_covariance = c("draw", "chain", "trait_row", "trait_col"),
  regional_marker_covariance = c("draw", "chain", "region", "trait_row", "trait_col"),
  predictions = c("draw", "chain", "observation", "trait")
)

.blr_st_raw_capture_state <- new.env(parent = emptyenv())
.blr_st_raw_capture_state$depth <- 0L

.blr_begin_st_raw_capture <- function() {
  previous <- .blr_st_raw_capture_state$depth
  .blr_st_raw_capture_state$depth <- previous + 1L
  previous
}

.blr_end_st_raw_capture <- function(previous) {
  .blr_st_raw_capture_state$depth <- as.integer(previous)
  invisible(NULL)
}

.blr_st_raw_capture_active <- function() {
  .blr_st_raw_capture_state$depth > 0L
}

.blr_make_array <- function(value, axes) {
  if (!is.list(axes) || is.null(names(axes)) || anyNA(names(axes)) ||
      any(!nzchar(names(axes))) || anyDuplicated(names(axes)) ||
      any(vapply(axes, function(x) !is.character(x) || anyNA(x) ||
        any(!nzchar(x)) || anyDuplicated(x), logical(1)))) {
    stop("Array axes must be a uniquely named list of stable character IDs.",
         call. = FALSE)
  }
  expected <- prod(lengths(axes))
  if (!is.atomic(value) || length(value) != expected) {
    stop("Array values do not match the contracted axis sizes.", call. = FALSE)
  }
  out <- array(value, dim = as.integer(lengths(axes)), dimnames = unname(axes))
  attr(out, "dim_axis_names") <- names(axes)
  out
}

.blr_validate_axis_array <- function(x, axes, field, finite = TRUE) {
  if (!is.array(x) || !identical(attr(x, "dim_axis_names"), axes) ||
      length(dim(x)) != length(axes) || is.null(dimnames(x)) ||
      length(dimnames(x)) != length(axes) ||
      any(vapply(dimnames(x), is.null, logical(1))) ||
      !identical(names(dimnames(x)), rep("", length(axes)))) {
    # R arrays intentionally keep axis names in dim_axis_names; dimnames are
    # the stable IDs and are un-named to avoid two competing axis sources.
    if (!is.array(x) || !identical(attr(x, "dim_axis_names"), axes) ||
        length(dim(x)) != length(axes) || is.null(dimnames(x)) ||
        any(vapply(dimnames(x), is.null, logical(1)))) {
      stop(field, " must have axes ", paste(axes, collapse = " x "),
           " with stable dimnames.", call. = FALSE)
    }
  }
  if (finite && (!is.numeric(x) || any(!is.finite(x)))) {
    stop(field, " must contain finite numeric values.", call. = FALSE)
  }
  invisible(TRUE)
}

.blr_required_raw_fields <- list(
  posterior = c(
    "realised_effect_mean", "latent_effect_mean", "scaled_effect_mean", "pips",
    "traitwise_state_probabilities", "joint_state_probabilities",
    "activity_pattern_probabilities",
    "traitwise_component_assignment_probabilities",
    "joint_component_assignment_probabilities", "marker_covariance_mean",
    "marker_variance_mean", "residual_covariance_mean",
    "residual_variance_mean", "uncertainty"),
  draws = c(
    "realised_effects", "latent_effects", "scaled_effects",
    "independent_trait_states", "joint_states", "traitwise_activity",
    "traitwise_probability_parameters", "joint_probability_parameters",
    "activity_pattern_parameters",
    "traitwise_component_probability_parameters",
    "joint_component_probability_parameters", "marker_covariance",
    "residual_covariance", "marker_variance", "residual_variance",
    "regional_marker_covariance", "convergence"),
  final = c(
    "realised_effects", "latent_effects", "scaled_effects",
    "independent_trait_states", "joint_states",
    "traitwise_probability_parameters", "joint_probability_parameters",
    "activity_pattern_parameters",
    "traitwise_component_probability_parameters",
    "joint_component_probability_parameters", "marker_covariance",
    "residual_covariance", "marker_variance", "residual_variance",
    "rng_continuation"),
  derived = c("predictions", "genetic_variance", "genomic_covariance",
              "operator_relative_quadratics", "descriptive_bilinear_forms"),
  diagnostics = c("convergence", "acceptance", "runtime", "memory", "workers",
                  "numerical_safeguards", "approximation_warnings"),
  provenance = c("package_version", "git_sha", "dirty_build", "compiler",
                 "operator_resources", "marker_alignment",
                 "seed_contract_version", "task_seeds", "timestamp")
)

.blr_fill_required_nulls <- function(x, fields) {
  if (is.null(x)) x <- list()
  .blr_exact_names(x, "blr_raw namespace")
  out <- stats::setNames(vector("list", length(fields)), fields)
  present <- intersect(fields, names(x))
  out[present] <- x[present]
  extras <- setdiff(names(x), fields)
  if (length(extras)) out <- c(out, x[extras])
  out
}

new_blr_raw_v2 <- function(model, input, posterior, draws, final,
                           derived = list(), diagnostics = list(),
                           provenance = list(),
                           compatibility_id = "phase1-r-v2",
                           source_schema = NULL,
                           migration = NULL) {
  validate_blr_resolved_spec(input)
  posterior <- .blr_fill_required_nulls(posterior,
                                        .blr_required_raw_fields$posterior)
  draws <- .blr_fill_required_nulls(draws, .blr_required_raw_fields$draws)
  final <- .blr_fill_required_nulls(final, .blr_required_raw_fields$final)
  derived <- .blr_fill_required_nulls(derived, .blr_required_raw_fields$derived)
  diagnostics <- .blr_fill_required_nulls(
    diagnostics, .blr_required_raw_fields$diagnostics)
  provenance <- .blr_fill_required_nulls(
    provenance, .blr_required_raw_fields$provenance)
  raw <- list(
    schema = list(
      name = "blr_raw", version = 2L,
      compatibility_id = compatibility_id,
      dimension_contract_version = 1L,
      seed_contract_version = input$schema$seed_contract_version,
      retention_contract_version = input$schema$retention_contract_version,
      source_schema = source_schema, migration = migration),
    model = model, input = input, posterior = posterior, draws = draws,
    final = final, derived = derived, diagnostics = diagnostics,
    provenance = provenance)
  validate_blr_raw_v2(raw)
  class(raw) <- c("blr_raw_v2", "blr_raw", "list")
  raw
}

is_blr_raw_v2 <- function(x) {
  is.list(x) && is.list(x$schema) && identical(x$schema$name, "blr_raw") &&
    identical(as.integer(x$schema$version), 2L)
}

.blr_probability_array <- function(x, field, axis) {
  if (is.null(x)) return(invisible(TRUE))
  if (any(x < -1e-12 | x > 1 + 1e-12)) {
    stop(field, " must lie in [0, 1].", call. = FALSE)
  }
  sums <- apply(x, setdiff(seq_along(dim(x)), axis), sum)
  if (any(abs(sums - 1) > 1e-8)) {
    stop(field, " must normalize over its probability axis.", call. = FALSE)
  }
  invisible(TRUE)
}

.blr_validate_covariance_array <- function(x, field) {
  if (is.null(x)) return(invisible(TRUE))
  d <- dim(x)
  nr <- d[length(d) - 1L]
  if (nr != d[length(d)]) stop(field, " must be square.", call. = FALSE)
  leading <- d[-c(length(d) - 1L, length(d))]
  index <- if (length(leading)) {
    arrayInd(seq_len(prod(leading)), leading)
  } else {
    matrix(integer(), nrow = 1L, ncol = 0L)
  }
  for (i in seq_len(nrow(index))) {
    subscripts <- c(as.list(index[i, ]), list(TRUE, TRUE), list(drop = TRUE))
    current <- do.call(`[`, c(list(x), subscripts))
    current <- matrix(current, nrow = nr, ncol = nr)
    if (any(!is.finite(current)) ||
        !isTRUE(all.equal(current, t(current), tolerance = 1e-10)) ||
        inherits(try(chol(current), silent = TRUE), "try-error")) {
      stop(field, " contains a non-finite, asymmetric, or non-positive-definite draw.",
           call. = FALSE)
    }
  }
  invisible(TRUE)
}

.blr_validate_trait_summary <- function(x, traits, field,
                                        strictly_positive = FALSE) {
  if (is.null(x)) return(invisible(TRUE))
  if (!is.numeric(x) || length(x) != length(traits) || any(!is.finite(x)) ||
      !identical(names(x), traits) || any(if (strictly_positive) x <= 0 else x < 0)) {
    stop(field, " must be a finite trait-named ",
         if (strictly_positive) "positive" else "nonnegative", " vector.",
         call. = FALSE)
  }
  invisible(TRUE)
}

.blr_validate_discrete_states <- function(x, field, maximum) {
  if (is.null(x)) return(invisible(TRUE))
  if (!is.numeric(x) || any(!is.finite(x)) || any(x != floor(x)) ||
      any(x < 0) || any(x > maximum)) {
    stop(field, " must contain zero-based declared state indices.",
         call. = FALSE)
  }
  invisible(TRUE)
}

.blr_validate_raw_axis_ids <- function(raw) {
  input <- raw$input
  regions <- input$data$statistical_regions
  region_ids <- if (is.null(regions)) character() else unique(as.character(regions))
  state_ids <- input$model$state_space
  component_ids <- if (input$model$family %in% c("bayesr", "bayesrc"))
    state_ids else character()
  joint_ids <- if (identical(input$data$analysis_mode, "joint_multitrait"))
    state_ids else character()
  ids <- list(
    draw = paste0("draw", seq_len(input$mcmc$retained_draws)),
    chain = paste0("chain", seq_len(input$mcmc$chains)),
    marker = input$data$global_markers,
    trait = input$data$trait_ids,
    trait_row = input$data$trait_ids,
    trait_col = input$data$trait_ids,
    state = state_ids,
    joint_state = joint_ids,
    component = component_ids,
    activity_pattern = joint_ids,
    region = region_ids,
    provider = names(input$data$providers))
  check_group <- function(group, group_name) {
    for (field in names(group)) {
      x <- group[[field]]
      if (!is.array(x)) next
      axes <- attr(x, "dim_axis_names")
      if (is.null(axes) || length(axes) != length(dim(x))) {
        stop(group_name, "$", field,
             " must declare one axis name per dimension.", call. = FALSE)
      }
      for (axis in axes) {
        if (identical(axis, "observation")) {
          observation_ids <- input$data$observation_ids %||% NULL
          if (is.null(observation_ids)) {
            stop(group_name, "$", field,
                 " uses an observation axis without declared observation IDs.",
                 call. = FALSE)
          }
          expected <- .blr_ids(observation_ids, "observation IDs")
        } else {
          expected <- ids[[axis]]
        }
        if (is.null(expected) || !length(expected)) {
          stop(group_name, "$", field, " uses undeclared axis ", axis, ".",
               call. = FALSE)
        }
        location <- match(axis, axes)
        if (!identical(dimnames(x)[[location]], expected)) {
          stop("Raw axis IDs do not match the resolved specification for ",
               group_name, "$", field, "$", axis, ".", call. = FALSE)
        }
      }
    }
  }
  check_group(raw$posterior, "posterior")
  check_group(raw$draws, "draws")
  check_group(raw$final, "final")
  check_group(raw$derived, "derived")
  check_group(raw$diagnostics, "diagnostics")
  invisible(TRUE)
}

.blr_validate_raw_provenance <- function(raw) {
  provenance <- raw$provenance
  .blr_exact_fields(provenance, "blr_raw$provenance",
                    .blr_required_raw_fields$provenance)
  character_or_null <- function(x, what) {
    if (!is.null(x) && (!is.character(x) || length(x) != 1L ||
                        is.na(x) || !nzchar(x))) {
      stop(what, " must be NULL or one nonempty character value.",
           call. = FALSE)
    }
  }
  .blr_character_scalar(provenance$package_version,
                        "provenance$package_version")
  if (!is.null(provenance$git_sha) &&
      (!is.character(provenance$git_sha) || length(provenance$git_sha) != 1L ||
       is.na(provenance$git_sha) ||
       !grepl("^[0-9a-fA-F]{7,40}$", provenance$git_sha))) {
    stop("provenance$git_sha must be NULL or a 7-40 digit hexadecimal Git identifier.",
         call. = FALSE)
  }
  if (!is.null(provenance$dirty_build) &&
      (!is.logical(provenance$dirty_build) ||
       length(provenance$dirty_build) != 1L ||
       is.na(provenance$dirty_build))) {
    stop("provenance$dirty_build must be NULL or one nonmissing logical value.",
         call. = FALSE)
  }
  character_or_null(provenance$compiler, "provenance$compiler")
  character_or_null(provenance$timestamp, "provenance$timestamp")
  if (!identical(provenance$operator_resources,
                 raw$input$data$operator_resources) ||
      !identical(provenance$marker_alignment, raw$input$data$provider_maps)) {
    stop("Provenance operator resources or marker alignment disagree with input.",
         call. = FALSE)
  }
  if (!identical(provenance$seed_contract_version,
                 raw$schema$seed_contract_version)) {
    stop("Provenance seed-contract version disagrees with schema.",
         call. = FALSE)
  }
  .blr_validate_task_seeds(provenance$task_seeds,
                           raw$input$data$analysis_mode,
                           raw$input$data$trait_ids,
                           raw$input$mcmc$chains)
  if (!identical(provenance$task_seeds, raw$input$mcmc$task_seeds)) {
    stop("Provenance task seeds disagree with the resolved input.",
         call. = FALSE)
  }
  invisible(TRUE)
}

validate_blr_raw_v2 <- function(raw) {
  .blr_exact_names(raw, "blr_raw")
  expected <- c("schema", "model", "input", "posterior", "draws", "final",
                "derived", "diagnostics", "provenance")
  if (!identical(names(raw), expected)) {
    stop("blr_raw has invalid namespace names or order.", call. = FALSE)
  }
  if (!is.list(raw$schema) || !identical(raw$schema$name, "blr_raw")) {
    stop("Unsupported blr_raw schema or dimension version.", call. = FALSE)
  }
  schema_fields <- c(
    "name", "version", "compatibility_id", "dimension_contract_version",
    "seed_contract_version", "retention_contract_version",
    "source_schema", "migration")
  .blr_exact_fields(raw$schema, "blr_raw$schema", schema_fields)
  .blr_scalar_whole(raw$schema$version, "blr_raw$schema$version", 2L, 2L)
  .blr_scalar_whole(raw$schema$dimension_contract_version,
                    "blr_raw$schema$dimension_contract_version", 1L, 1L)
  .blr_character_scalar(raw$schema$compatibility_id,
                        "blr_raw$schema$compatibility_id",
                        c("phase0-v2", "phase1-r-v2"))
  .blr_scalar_whole(raw$schema$seed_contract_version,
                    "blr_raw$schema$seed_contract_version", 0L, 1L)
  .blr_scalar_whole(raw$schema$retention_contract_version,
                    "blr_raw$schema$retention_contract_version", 0L, 1L)
  if (!is.null(raw$schema$source_schema) &&
      (!is.list(raw$schema$source_schema) ||
       is.data.frame(raw$schema$source_schema))) {
    stop("blr_raw$schema$source_schema must be NULL or a named descriptor.",
         call. = FALSE)
  }
  if (!is.null(raw$schema$source_schema)) {
    .blr_exact_names(raw$schema$source_schema,
                     "blr_raw$schema$source_schema")
  }
  if (!is.null(raw$schema$migration) &&
      (!is.list(raw$schema$migration) || is.data.frame(raw$schema$migration))) {
    stop("blr_raw$schema$migration must be NULL or a named descriptor.",
         call. = FALSE)
  }
  if (!is.null(raw$schema$migration)) {
    .blr_exact_names(raw$schema$migration, "blr_raw$schema$migration")
  }
  validate_blr_resolved_spec(raw$input)
  .blr_exact_fields(raw$model, "blr_raw$model", c(
    "analysis_mode", "family", "state_space", "null_state_index",
    "effect_storage_convention", "probability_policy",
    "marker_scale_policy", "marker_covariance_policy", "residual_policy",
    "update_order_version"))
  if (!identical(raw$model$analysis_mode, raw$input$data$analysis_mode) ||
      !identical(raw$model[names(raw$input$model)], raw$input$model) ||
      !identical(raw$schema$seed_contract_version,
                 raw$input$schema$seed_contract_version) ||
      !identical(raw$schema$retention_contract_version,
                 raw$input$schema$retention_contract_version)) {
    stop("blr_raw model or contract versions disagree with input.",
         call. = FALSE)
  }
  .blr_validate_raw_axis_ids(raw)
  forbidden <- c("pi", "pis", "pim", "state_probabilities",
                 "pattern_probabilities")
  for (group in c("posterior", "draws", "final")) {
    if (!is.list(raw[[group]])) stop("blr_raw$", group, " must be a list.",
                                     call. = FALSE)
    if (any(forbidden %in% names(raw[[group]]))) {
      stop("blr_raw v2 contains an ambiguous probability field.", call. = FALSE)
    }
  }
  for (group in names(.blr_required_raw_fields)) {
    if (!is.list(raw[[group]]) ||
        !all(.blr_required_raw_fields[[group]] %in% names(raw[[group]]))) {
      stop("blr_raw$", group, " is missing required-present fields.",
           call. = FALSE)
    }
  }
  for (field in intersect(names(.blr_raw_axis_contract), names(raw$posterior))) {
    if (!is.null(raw$posterior[[field]])) {
      .blr_validate_axis_array(raw$posterior[[field]],
                               .blr_raw_axis_contract[[field]],
                               paste0("posterior$", field))
    }
  }
  for (field in intersect(names(.blr_raw_axis_contract), names(raw$draws))) {
    if (!is.null(raw$draws[[field]])) {
      .blr_validate_axis_array(raw$draws[[field]],
                               .blr_raw_axis_contract[[field]],
                               paste0("draws$", field),
                               finite = !field %in% c("independent_trait_states",
                                                      "joint_states"))
    }
  }
  for (field in intersect(names(.blr_raw_axis_contract), names(raw$derived))) {
    if (!is.null(raw$derived[[field]])) {
      .blr_validate_axis_array(raw$derived[[field]],
                               .blr_raw_axis_contract[[field]],
                               paste0("derived$", field))
    }
  }
  final_axes <- list(
    realised_effects = c("chain", "marker", "trait"),
    latent_effects = c("chain", "marker", "trait"),
    scaled_effects = c("chain", "marker", "trait"),
    independent_trait_states = c("chain", "marker", "trait"),
    joint_states = c("chain", "marker"),
    traitwise_probability_parameters = c("chain", "trait", "state"),
    joint_probability_parameters = c("chain", "joint_state"),
    activity_pattern_parameters = c("chain", "activity_pattern"),
    traitwise_component_probability_parameters = c("chain", "trait", "component"),
    joint_component_probability_parameters = c("chain", "component"),
    marker_variance = c("chain", "trait"),
    genetic_variance = c("chain", "trait"),
    residual_variance = c("chain", "trait"),
    marker_covariance = c("chain", "trait_row", "trait_col"),
    residual_covariance = c("chain", "trait_row", "trait_col"))
  for (field in intersect(names(final_axes), names(raw$final))) {
    if (!is.null(raw$final[[field]])) {
      .blr_validate_axis_array(raw$final[[field]], final_axes[[field]],
                               paste0("final$", field),
                               finite = !field %in% c("independent_trait_states",
                                                      "joint_states"))
    }
  }
  if (!is.null(raw$posterior$pips) &&
      any(raw$posterior$pips < -1e-12 | raw$posterior$pips > 1 + 1e-12)) {
    stop("posterior$pips must lie in [0, 1].", call. = FALSE)
  }
  {
    traits <- raw$input$data$trait_ids
    .blr_validate_trait_summary(raw$posterior$marker_variance_mean, traits,
                                "posterior$marker_variance_mean")
    .blr_validate_trait_summary(raw$posterior$residual_variance_mean, traits,
                                "posterior$residual_variance_mean", TRUE)
    .blr_validate_covariance_array(raw$posterior$marker_covariance_mean,
                                   "posterior$marker_covariance_mean")
    .blr_validate_covariance_array(raw$posterior$residual_covariance_mean,
                                   "posterior$residual_covariance_mean")
    maximum_state <- length(raw$input$model$state_space) - 1L
    .blr_validate_discrete_states(raw$draws$independent_trait_states,
                                  "draws$independent_trait_states",
                                  maximum_state)
    .blr_validate_discrete_states(raw$draws$joint_states,
                                  "draws$joint_states", maximum_state)
    .blr_validate_discrete_states(raw$final$independent_trait_states,
                                  "final$independent_trait_states",
                                  maximum_state)
    .blr_validate_discrete_states(raw$final$joint_states,
                                  "final$joint_states", maximum_state)
    .blr_validate_discrete_states(raw$draws$traitwise_activity,
                                  "draws$traitwise_activity", 1L)
    .blr_probability_array(raw$posterior$traitwise_state_probabilities,
                           "posterior$traitwise_state_probabilities", 3L)
    .blr_probability_array(raw$posterior$joint_state_probabilities,
                           "posterior$joint_state_probabilities", 2L)
    .blr_probability_array(raw$posterior$activity_pattern_probabilities,
                           "posterior$activity_pattern_probabilities", 2L)
    .blr_probability_array(
      raw$posterior$traitwise_component_assignment_probabilities,
      "posterior$traitwise_component_assignment_probabilities", 3L)
    .blr_probability_array(
      raw$posterior$joint_component_assignment_probabilities,
      "posterior$joint_component_assignment_probabilities", 2L)
    .blr_probability_array(raw$draws$traitwise_probability_parameters,
                           "draws$traitwise_probability_parameters", 4L)
    .blr_probability_array(raw$draws$joint_probability_parameters,
                           "draws$joint_probability_parameters", 3L)
    .blr_probability_array(raw$draws$activity_pattern_parameters,
                           "draws$activity_pattern_parameters", 3L)
    .blr_probability_array(
      raw$draws$traitwise_component_probability_parameters,
      "draws$traitwise_component_probability_parameters", 4L)
    .blr_probability_array(
      raw$draws$joint_component_probability_parameters,
      "draws$joint_component_probability_parameters", 3L)
    .blr_probability_array(raw$final$traitwise_probability_parameters,
                           "final$traitwise_probability_parameters", 3L)
    .blr_probability_array(raw$final$joint_probability_parameters,
                           "final$joint_probability_parameters", 2L)
    .blr_probability_array(raw$final$activity_pattern_parameters,
                           "final$activity_pattern_parameters", 2L)
    .blr_probability_array(
      raw$final$traitwise_component_probability_parameters,
      "final$traitwise_component_probability_parameters", 3L)
    .blr_probability_array(
      raw$final$joint_component_probability_parameters,
      "final$joint_component_probability_parameters", 2L)
    for (field in c("marker_covariance", "residual_covariance",
                    "regional_marker_covariance")) {
      .blr_validate_covariance_array(raw$draws[[field]], paste0("draws$", field))
    }
    for (field in c("marker_covariance", "residual_covariance")) {
      .blr_validate_covariance_array(raw$final[[field]], paste0("final$", field))
    }
  }
  .blr_validate_raw_provenance(raw)
  invisible(TRUE)
}

.blr_provenance_cache <- new.env(parent = emptyenv())

.blr_cached_provenance <- function(refresh = FALSE) {
  if (!refresh && exists("value", envir = .blr_provenance_cache,
                         inherits = FALSE)) {
    return(get("value", envir = .blr_provenance_cache, inherits = FALSE))
  }
  version <- tryCatch(as.character(utils::packageVersion("sblr")),
                      error = function(e) {
                        desc <- tryCatch(read.dcf("DESCRIPTION"),
                                         error = function(e2) NULL)
                        if (is.null(desc)) "unknown" else
                          unname(desc[1L, "Version"])
                      })
  sha_value <- Sys.getenv("SBLR_BUILD_GIT_SHA", unset = "")
  sha <- if (nzchar(sha_value) &&
             grepl("^[0-9a-fA-F]{7,40}$", sha_value)) sha_value else NULL
  dirty_value <- Sys.getenv("SBLR_BUILD_GIT_DIRTY", unset = "")
  dirty <- if (dirty_value %in% c("1", "true", "TRUE")) TRUE else if (
    dirty_value %in% c("0", "false", "FALSE")) FALSE else NULL
  value <- list(package_version = version, git_sha = NULL,
                dirty_build = NULL, compiler = NULL, timestamp = NULL)
  value["git_sha"] <- list(sha)
  value["dirty_build"] <- list(dirty)
  assign("value", value, envir = .blr_provenance_cache)
  value
}

.blr_st_component_probabilities <- function(raw, markers, traits) {
  value <- raw$component$prob
  if (!is.list(value) || length(value) != length(traits)) return(NULL)
  components <- raw$component$names %||%
    paste0("component_", seq_len(ncol(as.matrix(value[[1L]]))) - 1L)
  out <- array(NA_real_, c(length(markers), length(traits), length(components)),
               dimnames = list(markers, traits, components))
  for (trait in seq_along(traits)) out[, trait, ] <- as.matrix(value[[trait]])
  attr(out, "dim_axis_names") <- c("marker", "trait", "component")
  out
}

.blr_validate_st_v1_conversion_registry <- function(raw, spec) {
  raw_family <- switch(as.character(raw$meta$model),
    bayesc = "bayesc", sbayesc = "bayesc",
    bayesr = "bayesr", sbayesr = "bayesr",
    bayesrc = "bayesrc", sbayesrc = "bayesrc",
    NA_character_)
  if (is.na(raw_family) || !identical(raw_family, spec$model$family)) {
    stop("Legacy raw model family does not match the registered converter specification.",
         call. = FALSE)
  }
  resource_types <- unique(vapply(
    spec$data$operator_resources, `[[`, character(1), "operator_type"))
  if (length(resource_types) != 1L) {
    stop("The Phase 1 ST converter requires one declared operator-resource type.",
         call. = FALSE)
  }
  backend <- as.character(raw$meta$backend)
  backend_ok <- switch(resource_types,
    bed = grepl("bed", backend, fixed = TRUE),
    csr = grepl("csr", backend, fixed = TRUE),
    block_eigen = grepl("csr|block", backend),
    retained_rank_block_eigen = grepl("csr|block", backend),
    full_rank_block_eigen = grepl("csr|block", backend),
    FALSE)
  if (!isTRUE(backend_ok)) {
    stop("Legacy raw backend does not match the registered operator converter.",
         call. = FALSE)
  }
  expected_prior <- if (spec$model$family == "bayesrc") {
    "annotation_component"
  } else if (identical(spec$model$marker_scale_policy,
                       "annotation_log_variance")) {
    "annotation_log_variance"
  } else if (spec$model$family == "bayesr") {
    "component"
  } else {
    switch(spec$model$probability_policy,
      global = "global", fixed_marker = "marker",
      learned_logistic = "annotation", group = "group",
      annotation_probit_stick = "annotation_component", NA_character_)
  }
  if (is.na(expected_prior) ||
      !identical(as.character(raw$meta$prior_type), expected_prior)) {
    stop("Legacy raw probability policy does not match the registered converter specification.",
         call. = FALSE)
  }
  residual_policies <- if (resource_types %in% c(
      "block_eigen", "retained_rank_block_eigen", "full_rank_block_eigen")) {
    c("global_projected_legacy", "gctb_block", "fixed_block")
  } else {
    "scalar"
  }
  if (!identical(spec$model$marker_covariance_policy, "traitwise_scalar") ||
      !spec$model$residual_policy %in% residual_policies ||
      !identical(spec$model$effect_storage_convention, "realised")) {
    stop("The Phase 1 ST converter supports only registered scalar covariance/residual and realised-effect policies.",
         call. = FALSE)
  }
  invisible(TRUE)
}

convert_stblr_raw_v1_to_blr_raw_v2 <- function(raw, spec,
                                                trait_ids = spec$data$trait_ids,
                                                marker_ids = spec$data$global_markers,
                                                convergence = NULL,
                                                formatted_pips = NULL) {
  .validate_stblr_raw(raw)
  validate_blr_resolved_spec(spec)
  .blr_validate_st_v1_conversion_registry(raw, spec)
  if (identical(spec$data$analysis_mode, "joint_multitrait")) {
    stop("Current MT covariance output is not convertible to sampled blr_raw v2 V_b.",
         call. = FALSE)
  }
  if (spec$mcmc$chains != 1L) {
    stop(paste0(
      "stblr_raw v1 does not expose every chain's final effect vector; ",
      "multi-chain results remain on the legacy schema until the native result ",
      "contract is extended."), call. = FALSE)
  }
  if (raw$meta$m != length(marker_ids) || raw$meta$nt != length(trait_ids)) {
    stop("Legacy raw dimensions do not match the resolved marker and trait IDs.",
         call. = FALSE)
  }
  marker_axes <- list(marker = marker_ids, trait = trait_ids)
  effect_mean <- .blr_make_array(as.numeric(raw$marker$bm), marker_axes)
  pips <- .blr_make_array(as.numeric(raw$marker$dm), marker_axes)
  if (!is.null(formatted_pips)) {
    candidate <- as.matrix(formatted_pips)
    if (!identical(dim(candidate), dim(pips)) ||
        !identical(dimnames(candidate), dimnames(pips)) ||
        !isTRUE(all.equal(unname(candidate), unname(pips),
                          tolerance = 1e-12, check.attributes = FALSE))) {
      stop("Formatted legacy PIPs disagree with the canonical raw PIPs.",
           call. = FALSE)
    }
    pips[] <- candidate
  }
  family <- spec$model$family
  component_prob <- if (family %in% c("bayesr", "bayesrc")) {
    .blr_st_component_probabilities(raw, marker_ids, trait_ids)
  } else NULL
  state_prob <- if (family == "bayesc") {
    value <- array(0, c(length(marker_ids), length(trait_ids), 2L),
                   dimnames = list(marker_ids, trait_ids,
                                   spec$model$state_space))
    value[, , 2L] <- pips
    value[, , 1L] <- 1 - pips
    attr(value, "dim_axis_names") <- c("marker", "trait", "state")
    value
  } else NULL
  postburn <- seq.int(raw$meta$nburn + 1L,
                      raw$meta$nburn + raw$meta$nit)
  retained <- raw$meta$nburn + spec$mcmc$retained_transition_indices
  trace_matrix <- function(x) if (is.null(x)) NULL else as.matrix(x)
  trait_mean <- function(x) {
    if (is.null(x)) return(NULL)
    stats::setNames(colMeans(trace_matrix(x)[postburn, , drop = FALSE]), trait_ids)
  }
  trait_diagonal <- function(x) {
    if (is.null(x)) return(NULL)
    stats::setNames(diag(as.matrix(x)), trait_ids)
  }
  variance_draw <- function(x) {
    if (is.null(x)) return(NULL)
    values <- trace_matrix(x)[retained, , drop = FALSE]
    .blr_make_array(as.numeric(values), list(
      draw = paste0("draw", seq_len(nrow(values))), chain = "chain1",
      trait = trait_ids))
  }
  final_effect <- .blr_make_array(as.numeric(raw$marker$b), list(
    chain = "chain1", marker = marker_ids, trait = trait_ids))
  final_state <- .blr_make_array(as.numeric(raw$marker$state), list(
    chain = "chain1", marker = marker_ids, trait = trait_ids))
  final_scalar <- function(x) {
    if (is.null(x)) return(NULL)
    .blr_make_array(diag(as.matrix(x)), list(chain = "chain1", trait = trait_ids))
  }
  pi_final <- as.matrix(raw$pi$final)
  pi_mean <- as.matrix(raw$pi$mean)
  probability_names <- if (family == "bayesc") {
    spec$model$state_space
  } else {
    raw$component$names %||% spec$model$state_space
  }
  if (length(probability_names) != ncol(pi_final)) {
    stop("Legacy probability parameters do not match the resolved state space.",
         call. = FALSE)
  }
  probability_axis <- if (family == "bayesc") "state" else "component"
  probability_axes <- list(chain = "chain1", trait = trait_ids)
  probability_axes[[probability_axis]] <- probability_names
  final_probability <- .blr_make_array(as.numeric(pi_final), probability_axes)
  probability_draws <- NULL
  if (family == "bayesc" && !is.null(raw$trace$pis)) {
    active <- trace_matrix(raw$trace$pis)[retained, , drop = FALSE]
    probability_draws <- array(
      NA_real_, c(nrow(active), 1L, length(trait_ids), 2L),
      dimnames = list(
        paste0("draw", seq_len(nrow(active))), "chain1", trait_ids,
        probability_names))
    probability_draws[, 1L, , 1L] <- 1 - active
    probability_draws[, 1L, , 2L] <- active
    attr(probability_draws, "dim_axis_names") <-
      c("draw", "chain", "trait", "state")
  }
  posterior <- list(
    realised_effect_mean = effect_mean, latent_effect_mean = NULL,
    scaled_effect_mean = NULL, pips = pips,
    traitwise_state_probabilities = state_prob,
    joint_state_probabilities = NULL, activity_pattern_probabilities = NULL,
    traitwise_component_assignment_probabilities = component_prob,
    joint_component_assignment_probabilities = NULL,
    marker_covariance_mean = NULL,
    marker_variance_mean = trait_diagonal(raw$variance$covb),
    residual_covariance_mean = NULL,
    residual_variance_mean = trait_diagonal(raw$variance$cove),
    uncertainty = NULL,
    traitwise_probability_parameter_mean = if (family == "bayesc")
      .blr_make_array(as.numeric(pi_mean), list(
        trait = trait_ids, state = probability_names)) else NULL,
    traitwise_component_probability_parameter_mean = if (family != "bayesc")
      .blr_make_array(as.numeric(pi_mean), list(
        trait = trait_ids, component = probability_names)) else NULL)
  draws <- list(
    realised_effects = NULL, latent_effects = NULL, scaled_effects = NULL,
    independent_trait_states = NULL, joint_states = NULL,
    traitwise_activity = NULL,
    traitwise_probability_parameters = probability_draws,
    joint_probability_parameters = NULL, activity_pattern_parameters = NULL,
    traitwise_component_probability_parameters = NULL,
    joint_component_probability_parameters = NULL,
    marker_covariance = NULL, residual_covariance = NULL,
    marker_variance = variance_draw(raw$trace$vbs),
    residual_variance = variance_draw(raw$trace$ves),
    regional_marker_covariance = NULL, convergence = convergence)
  final <- list(
    realised_effects = final_effect, latent_effects = NULL,
    scaled_effects = NULL, independent_trait_states = final_state,
    joint_states = NULL,
    traitwise_probability_parameters = if (family == "bayesc")
      final_probability else NULL,
    joint_probability_parameters = NULL, activity_pattern_parameters = NULL,
    traitwise_component_probability_parameters = if (family != "bayesc")
      final_probability else NULL,
    joint_component_probability_parameters = NULL, marker_covariance = NULL,
    residual_covariance = NULL,
    marker_variance = final_scalar(raw$variance$vb),
    genetic_variance = final_scalar(raw$variance$vg),
    residual_variance = final_scalar(raw$variance$ve), rng_continuation = NULL)
  legacy_iteration <- list(
    marker_variance = trace_matrix(raw$trace$vbs),
    genetic_variance = trace_matrix(raw$trace$vgs),
    residual_variance = trace_matrix(raw$trace$ves),
    le_genetic_contribution = trace_matrix(raw$trace$vle),
    ld_genetic_contribution = trace_matrix(raw$trace$vld),
    non_null_probability_parameter = trace_matrix(raw$trace$pis),
    includes_burn_in = TRUE,
    iteration_contract = "legacy_full_transition_trace_v1")
  provenance <- c(.blr_cached_provenance(), list(
    operator_resources = spec$data$operator_resources,
    marker_alignment = spec$data$provider_maps,
    seed_contract_version = spec$schema$seed_contract_version,
    task_seeds = spec$mcmc$task_seeds))
  new_blr_raw_v2(
    model = c(list(analysis_mode = spec$data$analysis_mode), spec$model),
    input = spec, posterior = posterior, draws = draws, final = final,
    derived = list(
      predictions = NULL,
      genetic_variance = trait_diagonal(raw$variance$covg),
      genomic_covariance = NULL,
      operator_relative_quadratics = list(
        le_mean = trait_mean(raw$trace$vle),
        ld_mean = trait_mean(raw$trace$vld)),
      descriptive_bilinear_forms = NULL,
      legacy_iteration_quantities = legacy_iteration),
    diagnostics = list(
      convergence = convergence, acceptance = NULL,
      runtime = raw$diagnostics[c("seconds_mean", "seconds_max")],
      memory = NULL, workers = raw$diagnostics$workers %||% NULL,
      numerical_safeguards = NULL,
      approximation_warnings = NULL, native = raw$diagnostics),
    provenance = provenance,
    source_schema = list(name = "stblr_raw", version = 1L,
                         backend = raw$meta$backend),
    migration = list(converter = "stblr_raw_v1_to_blr_raw_v2_v1",
                     retention = spec$schema$compatibility_id))
}

convert_blr_raw_v1_to_v2 <- function(raw, spec, ...) {
  if (.is_stblr_raw(raw)) {
    return(convert_stblr_raw_v1_to_blr_raw_v2(raw, spec, ...))
  }
  if (is.list(raw) && is.list(raw$schema) &&
      identical(raw$schema$class, "mtblr_raw") &&
      identical(as.integer(raw$schema$version), 1L)) {
    stop(paste0(
      "mtblr_raw v1 uses the current covariance hybrid and cannot be ",
      "converted into schema-v2 sampled marker covariance; rerun under the ",
      "future corrected joint MT model."), call. = FALSE)
  }
  stop("No registered converter for the supplied raw schema and version.",
       call. = FALSE)
}

.blr_apply_v2_st_aliases <- function(fit, raw) {
  replace_values <- function(template, value, field) {
    if (is.null(value)) return(template)
    if (is.null(template) || length(template) != length(value)) {
      stop("Formatted alias ", field,
           " is incompatible with its schema-v2 source.", call. = FALSE)
    }
    template[] <- as.numeric(value)
    template
  }
  diagonal_matrix <- function(value, trait_ids) {
    if (is.null(value)) return(NULL)
    out <- diag(as.numeric(value), nrow = length(trait_ids))
    dimnames(out) <- list(trait_ids, trait_ids)
    out
  }
  fit$bm <- raw$posterior$realised_effect_mean
  fit$dm <- raw$posterior$pips
  attr(fit$bm, "dim_axis_names") <- NULL
  attr(fit$dm, "dim_axis_names") <- NULL
  fit$b <- raw$final$realised_effects[1L, , , drop = TRUE]
  fit$d <- raw$final$independent_trait_states[1L, , , drop = TRUE]
  if (is.null(dim(fit$b))) fit$b <- matrix(fit$b, ncol = 1L)
  if (is.null(dim(fit$d))) fit$d <- matrix(fit$d, ncol = 1L)
  dimnames(fit$b) <- dimnames(fit$bm)
  dimnames(fit$d) <- dimnames(fit$dm)
  legacy <- raw$derived$legacy_iteration_quantities
  fit$vbs <- replace_values(fit$vbs, legacy$marker_variance, "vbs")
  fit$vgs <- replace_values(fit$vgs, legacy$genetic_variance, "vgs")
  fit$ves <- replace_values(fit$ves, legacy$residual_variance, "ves")
  fit$vle <- replace_values(
    fit$vle, legacy$le_genetic_contribution, "vle")
  fit$vld <- replace_values(
    fit$vld, legacy$ld_genetic_contribution, "vld")
  fit$pi_trace <- replace_values(
    fit$pi_trace, legacy$non_null_probability_parameter, "pi_trace")
  posterior_probability <- if (raw$model$family == "bayesc") {
    raw$posterior$traitwise_probability_parameter_mean
  } else {
    raw$posterior$traitwise_component_probability_parameter_mean
  }
  final_probability <- if (raw$model$family == "bayesc") {
    raw$final$traitwise_probability_parameters
  } else {
    raw$final$traitwise_component_probability_parameters
  }
  fit$pi_mean <- replace_values(fit$pi_mean, posterior_probability, "pi_mean")
  fit$pi_final <- replace_values(fit$pi_final, final_probability, "pi_final")
  if (!is.null(raw$posterior$traitwise_component_assignment_probabilities)) {
    values <- raw$posterior$traitwise_component_assignment_probabilities
    for (trait in seq_along(fit$component_probabilities)) {
      fit$component_probabilities[[trait]] <- replace_values(
        fit$component_probabilities[[trait]], values[, trait, ],
        "component_probabilities")
    }
  }
  traits <- raw$input$data$trait_ids
  fit$cov_b_mean <- replace_values(
    fit$cov_b_mean,
    diagonal_matrix(raw$posterior$marker_variance_mean, traits),
    "cov_b_mean")
  fit$cov_g_mean <- replace_values(
    fit$cov_g_mean,
    diagonal_matrix(raw$derived$genetic_variance, traits),
    "cov_g_mean")
  fit$cov_e_mean <- replace_values(
    fit$cov_e_mean,
    diagonal_matrix(raw$posterior$residual_variance_mean, traits),
    "cov_e_mean")
  fit$cov_b_final <- replace_values(
    fit$cov_b_final,
    diagonal_matrix(raw$final$marker_variance, traits),
    "cov_b_final")
  fit$cov_g_final <- replace_values(
    fit$cov_g_final,
    diagonal_matrix(raw$final$genetic_variance, traits),
    "cov_g_final")
  fit$cov_e_final <- replace_values(
    fit$cov_e_final,
    diagonal_matrix(raw$final$residual_variance, traits),
    "cov_e_final")
  attr(fit, "blr_raw") <- raw
  fit$diagnostics$schema_v2_alias_sources <- list(
    bm = "posterior$realised_effect_mean", dm = "posterior$pips",
    b = "final$realised_effects", d = "final$independent_trait_states",
    vbs = "derived$legacy_iteration_quantities$marker_variance",
    vgs = "derived$legacy_iteration_quantities$genetic_variance",
    ves = "derived$legacy_iteration_quantities$residual_variance",
    vle = "derived$legacy_iteration_quantities$le_genetic_contribution",
    vld = "derived$legacy_iteration_quantities$ld_genetic_contribution",
    pi_trace = paste0(
      "derived$legacy_iteration_quantities$non_null_probability_parameter; ",
      "legacy full transition trace including burn-in"),
    pi_mean = if (raw$model$family == "bayesc")
      "posterior$traitwise_probability_parameter_mean" else
      "posterior$traitwise_component_probability_parameter_mean",
    pi_final = if (raw$model$family == "bayesc")
      "final$traitwise_probability_parameters" else
      "final$traitwise_component_probability_parameters",
    component_probabilities =
      "posterior$traitwise_component_assignment_probabilities",
    cov_b_mean = "posterior$marker_variance_mean",
    cov_g_mean = "derived$genetic_variance",
    cov_e_mean = "posterior$residual_variance_mean",
    cov_b_final = "final$marker_variance",
    cov_g_final = "final$genetic_variance",
    cov_e_final = "final$residual_variance")
  fit
}

.blr_phase1_finalize_st <- function(fit, chain, model, operator) {
  spec <- attr(fit, "blr_resolved_spec", exact = TRUE)
  if (is.null(spec)) {
    spec <- resolve_blr_spec_from_wrapper(
      model, operator, colnames(fit$bm), rownames(fit$bm), chain,
      sample_sizes = fit$input$n %||% fit$input$n_used,
      probability_policy = fit$input$probability_policy %||% "global",
      marker_scale_policy = fit$input$effect_scale_policy %||% "unit",
      update_flags = list(
        marker_effects = isTRUE(fit$input$updateB %||% TRUE),
        residual_variance = isTRUE(fit$input$updateE %||% TRUE),
        probability = isTRUE(fit$input$updatePi %||% TRUE)))
  }
  spec <- .blr_complete_legacy_spec_from_fit(spec, fit)
  attr(fit, "blr_resolved_spec") <- spec
  legacy <- attr(fit, "stblr_raw_v1", exact = TRUE)
  attr(fit, "stblr_raw_v1") <- NULL
  if (is.null(legacy)) {
    fit$diagnostics$schema_v2_migration <- list(
      status = "legacy", reason = "The canonical native stblr_raw v1 object was not retained by this adapter.")
    return(fit)
  }
  converted <- tryCatch(
    convert_stblr_raw_v1_to_blr_raw_v2(
      legacy, spec, convergence = fit$convergence_traces %||% NULL,
      formatted_pips = fit$dm),
    error = identity)
  if (inherits(converted, "error")) {
    fit$diagnostics$schema_v2_migration <- list(
      status = "legacy", reason = conditionMessage(converted),
      source_schema = "stblr_raw_v1")
    return(fit)
  }
  fit <- .blr_apply_v2_st_aliases(fit, converted)
  fit$diagnostics$schema_v2_migration <- list(
    status = "converted", source_schema = "stblr_raw_v1",
    target_schema = "blr_raw_v2",
    converter = "stblr_raw_v1_to_blr_raw_v2_v1")
  fit
}

.new_blr_raw_v2 <- new_blr_raw_v2
.validate_blr_raw_v2 <- validate_blr_raw_v2
.is_blr_raw_v2 <- is_blr_raw_v2
.convert_stblr_raw_v1_to_v2 <- convert_stblr_raw_v1_to_blr_raw_v2
.convert_blr_raw_v1_to_v2 <- convert_blr_raw_v1_to_v2
