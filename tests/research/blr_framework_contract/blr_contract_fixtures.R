# Standalone base-R fixtures for the proposed unified BLR Phase 0 contracts.
# This file intentionally has no package or testthat dependency.

.blr_stop <- function(...) stop(..., call. = FALSE)

.blr_named_list <- function(x, what) {
  if (!is.list(x)) .blr_stop(what, " must be a list.")
  nm <- names(x)
  if (is.null(nm) || length(nm) != length(x) || anyNA(nm) || any(!nzchar(nm)) ||
      anyDuplicated(nm)) {
    .blr_stop(what, " must have unique, nonempty, non-NA names.")
  }
  x
}

.blr_scalar_integer <- function(x, what, lower = 0) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x) ||
      x != floor(x) || x < lower) {
    .blr_stop(what, " must be a finite integer scalar >= ", lower, ".")
  }
  as.numeric(x)
}

blr_retained_indices_v1 <- function(burn_in_iterations, sampling_iterations,
                                    thin_interval, retained_requested = TRUE) {
  burn <- .blr_scalar_integer(burn_in_iterations, "burn_in_iterations", 0)
  sampling <- .blr_scalar_integer(sampling_iterations, "sampling_iterations", 1)
  thin <- .blr_scalar_integer(thin_interval, "thin_interval", 1)
  post_burn <- seq_len(sampling)
  retained_post_burn <- post_burn[post_burn %% thin == 0]
  if (isTRUE(retained_requested) && !length(retained_post_burn)) {
    .blr_stop("retained output was requested but the retention rule keeps no draws.")
  }
  list(
    post_burn = retained_post_burn,
    absolute_transition = burn + retained_post_burn,
    retained_draws = length(retained_post_burn)
  )
}

# Unsigned 64-bit arithmetic represented by four little-endian base-2^16 limbs.
.u64_base <- 65536

.u64 <- function(limbs = numeric(4)) {
  if (length(limbs) != 4L || any(!is.finite(limbs)) || any(limbs < 0) ||
      any(limbs >= .u64_base) || any(limbs != floor(limbs))) {
    .blr_stop("invalid uint64 limbs.")
  }
  as.numeric(limbs)
}

.u64_from_uint32 <- function(x) {
  x <- .blr_scalar_integer(x, "uint32", 0)
  if (x > 4294967295) .blr_stop("uint32 must be <= 2^32 - 1.")
  .u64(c(x %% .u64_base, floor(x / .u64_base), 0, 0))
}

.u64_from_hex <- function(x) {
  if (!is.character(x) || length(x) != 1L ||
      !grepl("^[0-9A-Fa-f]{16}$", x)) .blr_stop("hex uint64 must have 16 digits.")
  chunks <- substring(x, c(13, 9, 5, 1), c(16, 12, 8, 4))
  .u64(strtoi(chunks, base = 16L))
}

.u64_hex <- function(x) {
  x <- .u64(x)
  paste0(sprintf("%04x", rev(as.integer(x))), collapse = "")
}

.xor16 <- function(a, b) {
  a <- as.integer(a); b <- as.integer(b)
  lo <- bitwXor(bitwAnd(a, 255L), bitwAnd(b, 255L))
  hi <- bitwXor(bitwShiftR(a, 8L), bitwShiftR(b, 8L))
  as.numeric(lo + 256L * hi)
}

.u64_xor <- function(a, b) .u64(mapply(.xor16, .u64(a), .u64(b)))

.u64_add <- function(a, b) {
  a <- .u64(a); b <- .u64(b); out <- numeric(4); carry <- 0
  for (i in seq_len(4)) {
    total <- a[[i]] + b[[i]] + carry
    out[[i]] <- total %% .u64_base
    carry <- floor(total / .u64_base)
  }
  .u64(out)
}

.u64_mul <- function(a, b) {
  a <- .u64(a); b <- .u64(b); accum <- numeric(8)
  for (i in seq_len(4)) for (j in seq_len(4)) accum[[i + j - 1L]] <-
    accum[[i + j - 1L]] + a[[i]] * b[[j]]
  for (i in seq_len(7)) {
    carry <- floor(accum[[i]] / .u64_base)
    accum[[i]] <- accum[[i]] %% .u64_base
    accum[[i + 1L]] <- accum[[i + 1L]] + carry
  }
  .u64(accum[seq_len(4)] %% .u64_base)
}

.u64_rshift <- function(x, n) {
  x <- .u64(x); n <- .blr_scalar_integer(n, "shift", 0)
  if (n >= 64) return(.u64())
  words <- floor(n / 16); bits <- n %% 16; out <- numeric(4)
  for (i in seq_len(4)) {
    source <- i + words
    if (source <= 4) {
      out[[i]] <- floor(x[[source]] / 2^bits)
      if (bits > 0 && source < 4) {
        out[[i]] <- out[[i]] + (x[[source + 1L]] %% 2^bits) * 2^(16 - bits)
      }
    }
  }
  .u64(out)
}

.fnv1a64 <- function(text) {
  if (!is.character(text) || length(text) != 1L || is.na(text) || !nzchar(text)) {
    .blr_stop("trait identity must be one nonempty UTF-8 string.")
  }
  hash <- .u64_from_hex("cbf29ce484222325")
  prime <- .u64_from_hex("00000100000001b3")
  bytes <- as.integer(charToRaw(enc2utf8(text)))
  for (byte in bytes) {
    operand <- .u64(c(byte, 0, 0, 0))
    hash <- .u64_mul(.u64_xor(hash, operand), prime)
  }
  hash
}

.splitmix64 <- function(x) {
  z <- .u64_add(.u64(x), .u64_from_hex("9e3779b97f4a7c15"))
  z <- .u64_mul(.u64_xor(z, .u64_rshift(z, 30)),
                .u64_from_hex("bf58476d1ce4e5b9"))
  z <- .u64_mul(.u64_xor(z, .u64_rshift(z, 27)),
                .u64_from_hex("94d049bb133111eb"))
  .u64_xor(z, .u64_rshift(z, 31))
}

.blr_seed_v1 <- function(user_seed, analysis_mode, trait_id, chain_index) {
  mode_code <- c(single_trait = 1, independent_traits = 2, joint_multitrait = 3)
  if (!analysis_mode %in% names(mode_code)) .blr_stop("invalid analysis_mode.")
  user_seed <- .blr_scalar_integer(user_seed, "seed", 0)
  if (user_seed > 4294967295) .blr_stop("seed must be <= 2^32 - 1.")
  chain_index <- .blr_scalar_integer(chain_index, "chain_index", 0)
  if (chain_index > 4294967295) .blr_stop("chain_index must be <= 2^32 - 1.")
  identity <- switch(
    analysis_mode,
    single_trait = "sblr:single_trait",
    independent_traits = trait_id,
    joint_multitrait = "sblr:joint_multitrait"
  )
  x <- .u64_from_uint32(user_seed)
  x <- .splitmix64(.u64_xor(x, .u64_from_uint32(1)))
  x <- .splitmix64(.u64_xor(x, .u64_from_uint32(unname(mode_code[[analysis_mode]]))))
  x <- .splitmix64(.u64_xor(x, .fnv1a64(identity)))
  x <- .splitmix64(.u64_xor(x, .u64_from_uint32(chain_index)))
  folded <- .u64_xor(x, .u64_rshift(x, 32))
  folded[[1]] + .u64_base * folded[[2]]
}

blr_task_seeds_v1 <- function(seed, analysis_mode, trait_ids, chains) {
  chains <- as.integer(.blr_scalar_integer(chains, "chains", 1))
  if (!is.character(trait_ids) || !length(trait_ids) || anyNA(trait_ids) ||
      any(!nzchar(trait_ids)) || anyDuplicated(trait_ids)) {
    .blr_stop("trait_ids must be unique, nonempty strings.")
  }
  chain_ids <- paste0("chain", seq_len(chains))
  if (analysis_mode == "independent_traits") {
    out <- matrix(0, length(trait_ids), chains,
                  dimnames = list(trait = trait_ids, chain = chain_ids))
    for (t in seq_along(trait_ids)) for (ch in seq_len(chains)) {
      out[t, ch] <- .blr_seed_v1(seed, analysis_mode, trait_ids[[t]], ch - 1L)
    }
    return(out)
  }
  if (!analysis_mode %in% c("single_trait", "joint_multitrait")) {
    .blr_stop("invalid analysis_mode.")
  }
  if (analysis_mode == "single_trait" && length(trait_ids) != 1L) {
    .blr_stop("single_trait requires exactly one trait ID.")
  }
  out <- vapply(seq_len(chains), function(ch) {
    .blr_seed_v1(seed, analysis_mode, trait_ids[[1]], ch - 1L)
  }, numeric(1))
  names(out) <- chain_ids
  out
}

blr_seed_reference_vectors_v1 <- function() {
  data.frame(
    user_seed = c(0, 0, 0, 0, 0, 17, 17, 17),
    analysis_mode = c("single_trait", "single_trait", "independent_traits",
                      "independent_traits", "joint_multitrait", "single_trait",
                      "independent_traits", "joint_multitrait"),
    trait_id = c("traitA", "traitA", "traitA", "traitB", "traitA", "traitA",
                 "traitB", "traitA"),
    chain_index = c(0, 1, 0, 1, 0, 0, 1, 1),
    native_seed = c(830191578, 160141543, 226943096, 286956759,
                    3100589946, 3397578794, 1132619387, 3700933392),
    stringsAsFactors = FALSE
  )
}

blr_valid_execution <- function(analysis_mode, execution_mode, parallelization) {
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
  if (execution_mode == "serial") return(parallelization == "none")
  parallelization %in% setdiff(allowed[[analysis_mode]], "none")
}

blr_operator_resource <- function(resource_id, operator_type, marker_ids,
                                  approximation = "exact_declared_operator") {
  list(
    resource_id = resource_id,
    operator_type = operator_type,
    marker_ids = marker_ids,
    alleles = data.frame(marker_id = marker_ids, effect = "A", other = "C"),
    genotype_coding = "effect_allele_count",
    centering = "declared",
    standardization = "declared",
    operator_scale = "cross_product",
    storage = list(kind = "fixture", payload = NULL),
    block_eigen = NULL,
    approximation = approximation,
    provenance = list(source = "Phase 0 fixture")
  )
}

blr_provider <- function(provider_id, trait_ids, resource_id, marker_ids,
                         global_markers, likelihood_regime, sample_size) {
  map <- match(marker_ids, global_markers)
  list(
    provider_id = provider_id,
    trait_ids = trait_ids,
    operator_resource_id = resource_id,
    local_to_global = setNames(map, marker_ids),
    sufficient_statistics = list(score = matrix(0, length(marker_ids), length(trait_ids),
      dimnames = list(marker = marker_ids, trait = trait_ids))),
    sample_size = setNames(rep(sample_size, length(trait_ids)), trait_ids),
    likelihood_regime = likelihood_regime,
    residual_contract = if (length(trait_ids) > 1L) "full_common_sample" else "marginal",
    population = "fixture_population",
    effect_scale = "declared_common_scale",
    overlap_group = NULL,
    provenance = list(source = "Phase 0 fixture")
  )
}

blr_validate_resources_providers <- function(resources, providers, global_markers,
                                             analysis_mode) {
  .blr_named_list(resources, "operator_resources")
  .blr_named_list(providers, "providers")
  ids <- vapply(resources, `[[`, character(1), "resource_id")
  if (anyDuplicated(ids) || !identical(names(resources), unname(ids))) .blr_stop("resource IDs must be unique and named identically.")
  pids <- vapply(providers, `[[`, character(1), "provider_id")
  if (anyDuplicated(pids) || !identical(names(providers), unname(pids))) .blr_stop("provider IDs must be unique and named identically.")
  for (p in providers) {
    if (!length(p$trait_ids) || anyNA(p$trait_ids) || any(!nzchar(p$trait_ids))) .blr_stop("providers require nonempty trait_ids.")
    if (!p$operator_resource_id %in% ids) .blr_stop("provider references an unknown operator resource.")
    resource <- resources[[p$operator_resource_id]]
    if (!identical(names(p$local_to_global), resource$marker_ids)) .blr_stop("provider map must follow resource-local marker order.")
    if (anyNA(p$local_to_global) || any(p$local_to_global < 1L) || any(p$local_to_global > length(global_markers))) .blr_stop("invalid local-to-global marker map.")
  }
  if (analysis_mode == "joint_multitrait") {
    joint <- vapply(providers, function(p) length(p$trait_ids) > 1L, logical(1))
    if (!any(joint)) .blr_stop("common-sample joint_multitrait data require a non-factorized multi-trait provider.")
  }
  TRUE
}

blr_spec_fixture <- function(analysis_mode = c("single_trait", "independent_traits", "joint_multitrait"),
                             residual_policy = c("fixed_full", "sampled_full")) {
  analysis_mode <- match.arg(analysis_mode); residual_policy <- match.arg(residual_policy)
  trait_ids <- switch(analysis_mode, single_trait = "trait1", c("trait1", "trait2"))
  global_markers <- c("rs1", "rs2", "rs3")
  bed <- blr_operator_resource("bed_shared", "bed", global_markers)
  resources <- list(bed_shared = bed)
  if (analysis_mode == "single_trait") {
    providers <- list(p1 = blr_provider("p1", trait_ids, "bed_shared", global_markers,
                                        global_markers, "common_sample", 100))
  } else if (analysis_mode == "independent_traits") {
    providers <- list(
      p1 = blr_provider("p1", "trait1", "bed_shared", global_markers, global_markers, "independent", 100),
      p2 = blr_provider("p2", "trait2", "bed_shared", global_markers, global_markers, "independent", 100)
    )
  } else {
    providers <- list(p_joint = blr_provider("p_joint", trait_ids, "bed_shared", global_markers,
                                              global_markers, "common_sample", 100))
  }
  burn_in <- if (analysis_mode == "single_trait") 1L else 2L
  sampling <- if (analysis_mode == "single_trait") 2L else 4L
  retained <- blr_retained_indices_v1(burn_in, sampling, 2)
  chains <- if (analysis_mode == "single_trait") 1 else 2
  task_seeds <- blr_task_seeds_v1(0, analysis_mode, trait_ids, chains)
  list(
    schema = list(name = "blr_resolved_spec", version = 1L,
                  compatibility_id = "phase0-v1", seed_contract_version = 1L,
                  retention_contract_version = 1L, dimension_contract_version = 1L),
    data = list(analysis_mode = analysis_mode, trait_ids = trait_ids,
                global_markers = global_markers,
                global_alleles = data.frame(marker_id = global_markers, effect = "A", other = "C"),
                operator_resources = resources, providers = providers,
                provider_maps = lapply(providers, `[[`, "local_to_global"),
                likelihood_regime = if (analysis_mode == "independent_traits") "independent_summary" else "common_sample",
                statistical_regions = NULL),
    model = list(family = "bayesc", state_space = c("null", "active"), null_state_index = 1L,
                 effect_storage_convention = "realised",
                 probability_policy = "global_inclusion",
                 marker_scale_policy = "unit", marker_covariance_policy = if (length(trait_ids) == 1L) "scalar" else "global_matrix",
                 residual_policy = residual_policy, update_order_version = 1L),
    prior = list(probability = list(alpha = 1, beta = 1), component_multipliers = NULL,
                 marker_multipliers = setNames(rep(1, length(global_markers)), global_markers),
                 scalar_variance = if (length(trait_ids) == 1L) list(shape = 2, scale = 1) else NULL,
                 marker_covariance = if (length(trait_ids) > 1L) list(degrees_of_freedom = 4, scale = diag(2), sampled = TRUE) else NULL,
                 residual_covariance = list(degrees_of_freedom = if (residual_policy == "sampled_full") 4 else NULL,
                                            scale = if (residual_policy == "sampled_full") diag(2) else NULL,
                                            fixed_value = if (residual_policy == "fixed_full") diag(length(trait_ids)) else NULL,
                                            sampled = residual_policy == "sampled_full"),
                 annotation = NULL),
    mcmc = list(burn_in_iterations = burn_in, sampling_iterations = sampling, thin_interval = 2L,
                retained_draws = retained$retained_draws,
                retained_transition_indices = retained$post_burn,
                chains = chains, seed = 0, task_seeds = task_seeds,
                update_flags = list(marker_effects = TRUE, marker_covariance = TRUE,
                                    residual_covariance = residual_policy == "sampled_full")),
    compute = list(execution_mode = "serial", parallelization = "none", cores = 1L,
                   scheduler_version = 1L, memory_limit_bytes = NULL,
                   operator_numerical_controls = list()),
    output = list(posterior_summaries = TRUE, retained_parameters = c("marker_covariance"),
                  effect_draw_policy = "full", state_draw_policy = "full",
                  convergence_policy = "core", derived_quantities = character(),
                  preserve_chains = TRUE, memory_estimate_bytes = 0)
  )
}

blr_summary_spec_fixture <- function() {
  spec <- blr_spec_fixture("independent_traits", "fixed_full")
  markers <- spec$data$global_markers
  dense <- blr_operator_resource("dense_t1", "dense_cross_product", markers[1:2])
  csr <- blr_operator_resource("csr_t2", "csr", markers[2:3])
  spec$data$operator_resources <- list(dense_t1 = dense, csr_t2 = csr)
  spec$data$providers <- list(
    p1 = blr_provider("p1", "trait1", "dense_t1", markers[1:2], markers, "independent_summary", 120),
    p2 = blr_provider("p2", "trait2", "csr_t2", markers[2:3], markers, "independent_summary", 85)
  )
  spec$data$provider_maps <- lapply(spec$data$providers, `[[`, "local_to_global")
  spec$data$likelihood_regime <- "independent_summary"
  spec
}

.blr_array <- function(dim, dimnames, value = 0) {
  if (length(dim) != length(dimnames) || is.null(names(dimnames))) .blr_stop("array dimensions require named dimnames.")
  array(value, dim = unname(dim), dimnames = unname(dimnames)) |>
    structure(dim_axis_names = names(dimnames))
}

.blr_identity <- function(ids) {
  out <- diag(length(ids))
  dimnames(out) <- list(ids, ids)
  attr(out, "dim_axis_names") <- c("trait_row", "trait_col")
  out
}

blr_raw_fixture <- function(analysis_mode = c("single_trait", "independent_traits", "joint_multitrait")) {
  analysis_mode <- match.arg(analysis_mode)
  spec <- blr_spec_fixture(analysis_mode, if (analysis_mode == "joint_multitrait") "sampled_full" else "fixed_full")
  traits <- spec$data$trait_ids; markers <- spec$data$global_markers
  chains <- paste0("chain", seq_len(spec$mcmc$chains)); draws <- paste0("draw", seq_len(spec$mcmc$retained_draws))
  states <- c("null", "active"); patterns <- if (analysis_mode == "joint_multitrait") c("00", "10", "01", "11") else character()
  d_eff <- list(draw = draws, chain = chains, marker = markers, trait = traits)
  effect_draws <- .blr_array(vapply(d_eff, length, integer(1)), d_eff)
  cov_dims <- list(draw = draws, chain = chains, trait_row = traits, trait_col = traits)
  final_eff_dims <- list(chain = chains, marker = markers, trait = traits)
  posterior <- list(
    realised_effect_mean = .blr_array(c(length(markers), length(traits)), list(marker = markers, trait = traits)),
    latent_effect_mean = NULL,
    scaled_effect_mean = NULL,
    pips = .blr_array(c(length(markers), length(traits)), list(marker = markers, trait = traits)),
    traitwise_state_probabilities = if (analysis_mode != "joint_multitrait") .blr_array(c(length(markers), length(traits), 2L), list(marker = markers, trait = traits, state = states)) else NULL,
    joint_state_probabilities = if (analysis_mode == "joint_multitrait") .blr_array(c(length(markers), length(patterns)), list(marker = markers, joint_state = patterns)) else NULL,
    activity_pattern_probabilities = if (analysis_mode == "joint_multitrait") .blr_array(c(length(markers), length(patterns)), list(marker = markers, activity_pattern = patterns)) else NULL,
    traitwise_component_assignment_probabilities = NULL,
    joint_component_assignment_probabilities = NULL,
    marker_covariance_mean = if (length(traits) > 1L) .blr_identity(traits) else NULL,
    marker_variance_mean = if (length(traits) == 1L) setNames(1, traits) else NULL,
    residual_covariance_mean = if (length(traits) > 1L) .blr_identity(traits) else NULL,
    residual_variance_mean = if (length(traits) == 1L) setNames(1, traits) else NULL,
    uncertainty = NULL
  )
  draws_list <- list(
    realised_effects = effect_draws,
    latent_effects = NULL,
    scaled_effects = NULL,
    independent_trait_states = if (analysis_mode != "joint_multitrait") .blr_array(vapply(d_eff, length, integer(1)), d_eff) else NULL,
    joint_states = if (analysis_mode == "joint_multitrait") .blr_array(c(length(draws), length(chains), length(markers)), list(draw = draws, chain = chains, marker = markers)) else NULL,
    traitwise_activity = effect_draws,
    traitwise_probability_parameters = if (analysis_mode != "joint_multitrait") .blr_array(c(length(draws), length(chains), length(traits), 2L), list(draw = draws, chain = chains, trait = traits, state = states)) else NULL,
    joint_probability_parameters = if (analysis_mode == "joint_multitrait") .blr_array(c(length(draws), length(chains), length(patterns)), list(draw = draws, chain = chains, joint_state = patterns)) else NULL,
    activity_pattern_parameters = if (analysis_mode == "joint_multitrait") .blr_array(c(length(draws), length(chains), length(patterns)), list(draw = draws, chain = chains, activity_pattern = patterns)) else NULL,
    traitwise_component_probability_parameters = NULL,
    joint_component_probability_parameters = NULL,
    marker_covariance = if (length(traits) > 1L) .blr_array(vapply(cov_dims, length, integer(1)), cov_dims) else NULL,
    residual_covariance = if (length(traits) > 1L) .blr_array(vapply(cov_dims, length, integer(1)), cov_dims) else NULL,
    marker_variance = if (length(traits) == 1L) .blr_array(c(length(draws), length(chains), 1L), list(draw = draws, chain = chains, trait = traits)) else NULL,
    residual_variance = if (length(traits) == 1L) .blr_array(c(length(draws), length(chains), 1L), list(draw = draws, chain = chains, trait = traits)) else NULL,
    regional_marker_covariance = NULL,
    convergence = NULL
  )
  final <- list(
    realised_effects = .blr_array(vapply(final_eff_dims, length, integer(1)), final_eff_dims),
    latent_effects = NULL,
    scaled_effects = NULL,
    independent_trait_states = if (analysis_mode != "joint_multitrait") .blr_array(vapply(final_eff_dims, length, integer(1)), final_eff_dims) else NULL,
    joint_states = if (analysis_mode == "joint_multitrait") .blr_array(c(length(chains), length(markers)), list(chain = chains, marker = markers)) else NULL,
    traitwise_probability_parameters = NULL,
    joint_probability_parameters = NULL,
    activity_pattern_parameters = NULL,
    traitwise_component_probability_parameters = NULL,
    joint_component_probability_parameters = NULL,
    marker_covariance = if (length(traits) > 1L) .blr_array(c(length(chains), length(traits), length(traits)), list(chain = chains, trait_row = traits, trait_col = traits)) else NULL,
    residual_covariance = if (length(traits) > 1L) .blr_array(c(length(chains), length(traits), length(traits)), list(chain = chains, trait_row = traits, trait_col = traits)) else NULL,
    marker_variance = if (length(traits) == 1L) .blr_array(c(length(chains), 1L), list(chain = chains, trait = traits)) else NULL,
    residual_variance = if (length(traits) == 1L) .blr_array(c(length(chains), 1L), list(chain = chains, trait = traits)) else NULL,
    rng_continuation = NULL
  )
  list(
    schema = list(name = "blr_raw", version = 2L, compatibility_id = "phase0-v2",
                  dimension_contract_version = 1L),
    model = spec$model,
    input = spec,
    posterior = posterior,
    draws = draws_list,
    final = final,
    derived = list(predictions = NULL, genetic_variance = NULL, genomic_covariance = NULL,
                   operator_relative_quadratics = NULL, descriptive_bilinear_forms = NULL),
    diagnostics = list(convergence = NULL, acceptance = NULL, runtime = NULL, memory = NULL,
                       workers = NULL, numerical_safeguards = NULL, approximation_warnings = NULL),
    provenance = list(package_version = "fixture", git_sha = NULL, dirty_build = NULL,
                      compiler = NULL, operator_resources = spec$data$operator_resources,
                      marker_alignment = spec$data$provider_maps,
                      seed_contract_version = 1L, task_seeds = spec$mcmc$task_seeds,
                      timestamp = NULL)
  )
}

blr_validate_spec_fixture <- function(spec) {
  .blr_named_list(spec, "blr_resolved_spec")
  expected <- c("schema", "data", "model", "prior", "mcmc", "compute", "output")
  if (!identical(names(spec), expected)) .blr_stop("blr_resolved_spec envelope names are invalid.")
  if (!blr_valid_execution(spec$data$analysis_mode, spec$compute$execution_mode,
                           spec$compute$parallelization)) .blr_stop("invalid analysis/execution combination.")
  retained <- blr_retained_indices_v1(spec$mcmc$burn_in_iterations,
                                      spec$mcmc$sampling_iterations,
                                      spec$mcmc$thin_interval)
  if (!identical(spec$mcmc$retained_transition_indices, retained$post_burn) ||
      spec$mcmc$retained_draws != retained$retained_draws) .blr_stop("retention fields are inconsistent.")
  expected_seed_dim <- if (spec$data$analysis_mode == "independent_traits") as.integer(c(length(spec$data$trait_ids), spec$mcmc$chains)) else as.integer(spec$mcmc$chains)
  if (!identical(dim(spec$mcmc$task_seeds), if (length(expected_seed_dim) == 1L) NULL else expected_seed_dim) ||
      length(spec$mcmc$task_seeds) != prod(expected_seed_dim) || anyNA(spec$mcmc$task_seeds) ||
      any(!is.finite(spec$mcmc$task_seeds)) || any(spec$mcmc$task_seeds < 0) ||
      any(spec$mcmc$task_seeds > 4294967295) || any(spec$mcmc$task_seeds != floor(spec$mcmc$task_seeds))) {
    .blr_stop("task_seeds have invalid shape or uint32 values.")
  }
  chain_names <- paste0("chain", seq_len(spec$mcmc$chains))
  if (spec$data$analysis_mode == "independent_traits") {
    if (!identical(dimnames(spec$mcmc$task_seeds),
                   list(trait = spec$data$trait_ids, chain = chain_names))) {
      .blr_stop("task_seeds have invalid trait or chain dimnames.")
    }
  } else if (!identical(names(spec$mcmc$task_seeds), chain_names)) {
    .blr_stop("task_seeds have invalid chain names.")
  }
  blr_validate_resources_providers(spec$data$operator_resources, spec$data$providers,
                                   spec$data$global_markers, spec$data$analysis_mode)
  trait_count <- length(spec$data$trait_ids)
  check_covariance_prior <- function(x, what) {
    if (is.null(x) || !isTRUE(x$sampled)) return(TRUE)
    if (!is.numeric(x$degrees_of_freedom) || length(x$degrees_of_freedom) != 1L ||
        !is.finite(x$degrees_of_freedom) || x$degrees_of_freedom <= trait_count - 1L) {
      .blr_stop(what, " degrees_of_freedom must define a proper inverse-Wishart distribution.")
    }
    if (!is.matrix(x$scale) || !identical(dim(x$scale), c(trait_count, trait_count)) ||
        any(!is.finite(x$scale)) || !isTRUE(all.equal(x$scale, t(x$scale), tolerance = 0)) ||
        inherits(try(chol(x$scale), silent = TRUE), "try-error")) {
      .blr_stop(what, " scale must be finite, symmetric, and positive definite.")
    }
    TRUE
  }
  check_covariance_prior(spec$prior$marker_covariance, "marker covariance")
  check_covariance_prior(spec$prior$residual_covariance, "residual covariance")
  TRUE
}

blr_validate_raw_fixture <- function(raw) {
  .blr_named_list(raw, "blr_raw")
  expected <- c("schema", "model", "input", "posterior", "draws", "final",
                "derived", "diagnostics", "provenance")
  if (!identical(names(raw), expected)) .blr_stop("blr_raw envelope names are invalid.")
  forbidden <- c("pi", "pis", "pim", "state_probabilities", "pattern_probabilities")
  if (any(forbidden %in% names(raw$posterior)) || any(forbidden %in% names(raw$draws))) .blr_stop("ambiguous probability field found.")
  arrays <- c(Filter(is.array, raw$posterior), Filter(is.array, raw$draws), Filter(is.array, raw$final))
  for (x in arrays) {
    if (is.null(attr(x, "dim_axis_names")) || length(dim(x)) != length(dimnames(x)) ||
        any(vapply(dimnames(x), is.null, logical(1)))) .blr_stop("raw array lacks fixed named axes or dimnames.")
  }
  required_present <- list(
    posterior = c("latent_effect_mean", "scaled_effect_mean",
                  "traitwise_component_assignment_probabilities",
                  "joint_component_assignment_probabilities", "uncertainty"),
    draws = c("latent_effects", "scaled_effects",
              "traitwise_component_probability_parameters",
              "joint_component_probability_parameters",
              "regional_marker_covariance", "convergence"),
    final = c("traitwise_probability_parameters", "joint_probability_parameters",
              "activity_pattern_parameters",
              "traitwise_component_probability_parameters",
              "joint_component_probability_parameters", "rng_continuation")
  )
  for (group in names(required_present)) if (!all(required_present[[group]] %in% names(raw[[group]]))) .blr_stop("required-present-NULL field is absent.")
  TRUE
}
