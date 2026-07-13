# Internal Phase 1 contracts for the future unified BLR framework.
# These helpers are deliberately not connected to a public fitting route.

.blr_scalar_integer <- function(x, field, minimum = NULL) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x) ||
      x != floor(x) || x < (-.Machine$integer.max - 1) ||
      x > .Machine$integer.max) {
    stop(field, " must be a finite integer scalar.", call. = FALSE)
  }
  value <- as.integer(x)
  if (!is.null(minimum) && value < minimum) {
    stop(field, " must be >= ", minimum, ".", call. = FALSE)
  }
  value
}

.blr_scalar_flag <- function(x, field) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    stop(field, " must be TRUE or FALSE.", call. = FALSE)
  }
  x
}

.blr_scalar_tag <- function(x, field) {
  if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x)) {
    stop(field, " must be a non-empty character scalar.", call. = FALSE)
  }
  x
}

.blr_validate_ids <- function(x, expected, field) {
  if (!is.character(x) || length(x) != expected || anyNA(x) ||
      any(!nzchar(x))) {
    stop(field, " must contain exactly ", expected,
         " non-empty character IDs.", call. = FALSE)
  }
  if (anyDuplicated(x)) {
    stop(field, " must not contain duplicate IDs.", call. = FALSE)
  }
  x
}

.blr_data_spec <- function(n_markers, n_traits, marker_ids, trait_ids,
                           sample_size,
                           scaling = "standardized_genotype",
                           resource_id) {
  list(
    representation = "csr",
    design = "independent_traits",
    n_markers = n_markers,
    n_traits = n_traits,
    marker_ids = marker_ids,
    trait_ids = trait_ids,
    sample_size = sample_size,
    scaling = scaling,
    csr = list(
      resource_id = resource_id,
      marker_count = n_markers,
      shared_read_only = TRUE,
      per_chain_data = FALSE,
      lifetime_exceeds_chains = TRUE
    )
  )
}

.blr_model_spec <- function() {
  list(
    kernel = "scalar",
    family = "bayesc",
    state = "binary",
    probability = "global_binary",
    scale = "unit",
    trait_covariance = "scalar_independent",
    residual_covariance = "scalar_independent"
  )
}

.blr_mcmc_control <- function(nit, nburn, nthin = 1L, nchains = 1L,
                              ncores = 1L, seed = 1L,
                              chain_seeds = NULL) {
  list(
    nit = nit,
    nburn = nburn,
    nthin = nthin,
    nchains = nchains,
    ncores = ncores,
    seed = seed,
    chain_seeds = chain_seeds
  )
}

.blr_output_spec <- function(marker_mean = TRUE, marker_pip = TRUE,
                             parameter_traces = TRUE, final_state = TRUE,
                             keep_chain_summaries = FALSE) {
  list(
    marker_mean = marker_mean,
    marker_pip = marker_pip,
    parameter_traces = parameter_traces,
    final_state = final_state,
    keep_chain_summaries = keep_chain_summaries
  )
}

.blr_execution_spec <- function() {
  list(
    operator = "csr",
    backend_reference = "stblr_cpg_omp_csr",
    scheduled = FALSE
  )
}

.blr_resolved_spec <- function(n_markers, n_traits, marker_ids, trait_ids,
                               sample_size, resource_id, nit, nburn,
                               nthin = 1L, nchains = 1L, ncores = 1L,
                               seed = 1L, chain_seeds = NULL,
                               scaling = "standardized_genotype",
                               keep_chain_summaries = FALSE) {
  spec <- list(
    schema = list(name = "blr_resolved_spec", version = 1L),
    data = .blr_data_spec(
      n_markers = n_markers,
      n_traits = n_traits,
      marker_ids = marker_ids,
      trait_ids = trait_ids,
      sample_size = sample_size,
      scaling = scaling,
      resource_id = resource_id
    ),
    model = .blr_model_spec(),
    mcmc = .blr_mcmc_control(
      nit = nit,
      nburn = nburn,
      nthin = nthin,
      nchains = nchains,
      ncores = ncores,
      seed = seed,
      chain_seeds = chain_seeds
    ),
    output = .blr_output_spec(
      keep_chain_summaries = keep_chain_summaries
    ),
    execution = .blr_execution_spec()
  )
  .validate_blr_resolved_spec(spec)
}

.validate_blr_resolved_spec <- function(spec) {
  required <- c("schema", "data", "model", "mcmc", "output", "execution")
  if (!is.list(spec) || !identical(names(spec), required)) {
    stop("spec must be a resolved BLR specification with fields: ",
         paste(required, collapse = ", "), ".", call. = FALSE)
  }

  if (!identical(spec$schema$name, "blr_resolved_spec")) {
    stop("schema$name must be 'blr_resolved_spec'.", call. = FALSE)
  }
  if (.blr_scalar_integer(spec$schema$version, "schema$version") != 1L) {
    stop("schema$version must be 1.", call. = FALSE)
  }

  n_markers <- .blr_scalar_integer(spec$data$n_markers,
                                    "data$n_markers", 1L)
  n_traits <- .blr_scalar_integer(spec$data$n_traits,
                                  "data$n_traits", 1L)
  .blr_validate_ids(spec$data$marker_ids, n_markers, "data$marker_ids")
  .blr_validate_ids(spec$data$trait_ids, n_traits, "data$trait_ids")
  if (!is.numeric(spec$data$sample_size) ||
      length(spec$data$sample_size) != n_traits ||
      anyNA(spec$data$sample_size) || any(!is.finite(spec$data$sample_size)) ||
      any(spec$data$sample_size != floor(spec$data$sample_size)) ||
      any(spec$data$sample_size <= 0) ||
      any(spec$data$sample_size > .Machine$integer.max)) {
    stop("data$sample_size must contain one positive integer per trait.",
         call. = FALSE)
  }
  if (!identical(.blr_scalar_tag(spec$data$representation,
                                 "data$representation"), "csr")) {
    stop("data$representation must be 'csr' in Phase 1.", call. = FALSE)
  }
  if (!identical(.blr_scalar_tag(spec$data$design, "data$design"),
                 "independent_traits")) {
    stop("data$design must be 'independent_traits' in Phase 1.",
         call. = FALSE)
  }
  if (!identical(.blr_scalar_tag(spec$data$scaling, "data$scaling"),
                 "standardized_genotype")) {
    stop("data$scaling must be 'standardized_genotype' in Phase 1.",
         call. = FALSE)
  }

  csr <- spec$data$csr
  if (!is.list(csr)) stop("data$csr must be a list.", call. = FALSE)
  .blr_scalar_tag(csr$resource_id, "data$csr$resource_id")
  if (.blr_scalar_integer(csr$marker_count, "data$csr$marker_count", 1L) !=
      n_markers) {
    stop("data$csr$marker_count must equal data$n_markers.", call. = FALSE)
  }
  if (!isTRUE(.blr_scalar_flag(csr$shared_read_only,
                               "data$csr$shared_read_only"))) {
    stop("data$csr$shared_read_only must be TRUE.", call. = FALSE)
  }
  if (isTRUE(.blr_scalar_flag(csr$per_chain_data,
                              "data$csr$per_chain_data"))) {
    stop("data$csr$per_chain_data must be FALSE.", call. = FALSE)
  }
  if (!isTRUE(.blr_scalar_flag(csr$lifetime_exceeds_chains,
                               "data$csr$lifetime_exceeds_chains"))) {
    stop("data$csr$lifetime_exceeds_chains must be TRUE.", call. = FALSE)
  }

  supported_model <- c(
    kernel = "scalar",
    family = "bayesc",
    state = "binary",
    probability = "global_binary",
    scale = "unit",
    trait_covariance = "scalar_independent",
    residual_covariance = "scalar_independent"
  )
  for (field in names(supported_model)) {
    value <- .blr_scalar_tag(spec$model[[field]], paste0("model$", field))
    if (!identical(value, unname(supported_model[[field]]))) {
      stop("model$", field, " must be '", supported_model[[field]],
           "' in Phase 1.", call. = FALSE)
    }
  }

  for (field in c("nit", "nthin", "nchains", "ncores")) {
    .blr_scalar_integer(spec$mcmc[[field]], paste0("mcmc$", field), 1L)
  }
  .blr_scalar_integer(spec$mcmc$nburn, "mcmc$nburn", 0L)
  .blr_scalar_integer(spec$mcmc$seed, "mcmc$seed")
  if (!is.null(spec$mcmc$chain_seeds)) {
    seeds <- spec$mcmc$chain_seeds
    if (!is.numeric(seeds) || length(seeds) != spec$mcmc$nchains ||
        anyNA(seeds) || any(!is.finite(seeds)) ||
        any(seeds != floor(seeds)) ||
        any(seeds < (-.Machine$integer.max - 1)) ||
        any(seeds > .Machine$integer.max)) {
      stop("mcmc$chain_seeds must be NULL or contain one integer per chain.",
           call. = FALSE)
    }
  }

  for (field in c("marker_mean", "marker_pip", "parameter_traces",
                  "final_state", "keep_chain_summaries")) {
    .blr_scalar_flag(spec$output[[field]], paste0("output$", field))
  }

  if (!identical(.blr_scalar_tag(spec$execution$operator,
                                 "execution$operator"), "csr")) {
    stop("execution$operator must be 'csr' in Phase 1.", call. = FALSE)
  }
  if (!identical(.blr_scalar_tag(spec$execution$backend_reference,
                                 "execution$backend_reference"),
                 "stblr_cpg_omp_csr")) {
    stop("execution$backend_reference must be 'stblr_cpg_omp_csr' in Phase 1.",
         call. = FALSE)
  }
  if (isTRUE(.blr_scalar_flag(spec$execution$scheduled,
                              "execution$scheduled"))) {
    stop("execution$scheduled must be FALSE; scheduled BayesC is unsupported in Phase 1.",
         call. = FALSE)
  }

  spec
}

.blr_validate_spec_cpp <- function(spec) {
  blr_phase1_validate_spec_cpp(.validate_blr_resolved_spec(spec))
}

.blr_validate_result_dimensions_cpp <- function(dimensions) {
  blr_phase1_validate_result_dimensions_cpp(dimensions)
}
