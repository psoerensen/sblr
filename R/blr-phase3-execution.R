.blr_u64_base <- 65536

.blr_u64 <- function(x = numeric(4L)) {
  if (!is.numeric(x) || length(x) != 4L || anyNA(x) ||
      any(!is.finite(x)) || any(x < 0) || any(x >= .blr_u64_base) ||
      any(x != floor(x))) {
    stop("Invalid uint64 limb representation.", call. = FALSE)
  }
  as.numeric(x)
}

.blr_u64_from_uint32 <- function(x) {
  x <- .blr_scalar_whole(x, "uint32", 0, 4294967295)
  .blr_u64(c(x %% .blr_u64_base, floor(x / .blr_u64_base), 0, 0))
}

.blr_u64_from_hex <- function(x) {
  if (!is.character(x) || length(x) != 1L || is.na(x) ||
      !grepl("^[0-9A-Fa-f]{16}$", x)) {
    stop("A uint64 hexadecimal constant must contain 16 digits.",
         call. = FALSE)
  }
  chunks <- substring(x, c(13L, 9L, 5L, 1L), c(16L, 12L, 8L, 4L))
  .blr_u64(strtoi(chunks, base = 16L))
}

.blr_xor16 <- function(a, b) {
  a <- as.integer(a)
  b <- as.integer(b)
  lo <- bitwXor(bitwAnd(a, 255L), bitwAnd(b, 255L))
  hi <- bitwXor(bitwShiftR(a, 8L), bitwShiftR(b, 8L))
  as.numeric(lo + 256L * hi)
}

.blr_u64_xor <- function(a, b) {
  .blr_u64(mapply(.blr_xor16, .blr_u64(a), .blr_u64(b)))
}

.blr_u64_add <- function(a, b) {
  a <- .blr_u64(a)
  b <- .blr_u64(b)
  out <- numeric(4L)
  carry <- 0
  for (index in seq_len(4L)) {
    total <- a[[index]] + b[[index]] + carry
    out[[index]] <- total %% .blr_u64_base
    carry <- floor(total / .blr_u64_base)
  }
  .blr_u64(out)
}

.blr_u64_mul <- function(a, b) {
  a <- .blr_u64(a)
  b <- .blr_u64(b)
  accum <- numeric(8L)
  for (left in seq_len(4L)) for (right in seq_len(4L)) {
    position <- left + right - 1L
    accum[[position]] <- accum[[position]] + a[[left]] * b[[right]]
  }
  for (index in seq_len(7L)) {
    carry <- floor(accum[[index]] / .blr_u64_base)
    accum[[index]] <- accum[[index]] %% .blr_u64_base
    accum[[index + 1L]] <- accum[[index + 1L]] + carry
  }
  .blr_u64(accum[seq_len(4L)] %% .blr_u64_base)
}

.blr_u64_rshift <- function(x, bits) {
  x <- .blr_u64(x)
  bits <- .blr_scalar_whole(bits, "uint64 shift", 0, 64)
  if (bits >= 64) return(.blr_u64())
  words <- floor(bits / 16)
  remainder <- bits %% 16
  out <- numeric(4L)
  for (index in seq_len(4L)) {
    source <- index + words
    if (source <= 4L) {
      out[[index]] <- floor(x[[source]] / 2^remainder)
      if (remainder > 0 && source < 4L) {
        out[[index]] <- out[[index]] +
          (x[[source + 1L]] %% 2^remainder) * 2^(16 - remainder)
      }
    }
  }
  .blr_u64(out)
}

.blr_fnv1a64_utf8 <- function(identity) {
  identity <- .blr_character_scalar(identity, "logical task identity")
  hash <- .blr_u64_from_hex("cbf29ce484222325")
  prime <- .blr_u64_from_hex("00000100000001b3")
  for (byte in as.integer(charToRaw(enc2utf8(identity)))) {
    hash <- .blr_u64_mul(
      .blr_u64_xor(hash, .blr_u64(c(byte, 0, 0, 0))), prime)
  }
  hash
}

.blr_splitmix64 <- function(x) {
  z <- .blr_u64_add(.blr_u64(x), .blr_u64_from_hex("9e3779b97f4a7c15"))
  z <- .blr_u64_mul(
    .blr_u64_xor(z, .blr_u64_rshift(z, 30)),
    .blr_u64_from_hex("bf58476d1ce4e5b9"))
  z <- .blr_u64_mul(
    .blr_u64_xor(z, .blr_u64_rshift(z, 27)),
    .blr_u64_from_hex("94d049bb133111eb"))
  .blr_u64_xor(z, .blr_u64_rshift(z, 31))
}

.blr_seed_v1 <- function(user_seed, analysis_mode, trait_id, chain_index) {
  mode_codes <- c(single_trait = 1, independent_traits = 2,
                  joint_multitrait = 3)
  analysis_mode <- .blr_character_scalar(
    analysis_mode, "analysis_mode", names(mode_codes))
  user_seed <- .blr_scalar_whole(user_seed, "seed", 0, 4294967295)
  chain_index <- .blr_scalar_whole(
    chain_index, "zero-based chain index", 0, 4294967295)
  identity <- switch(
    analysis_mode,
    single_trait = "sblr:single_trait",
    independent_traits = .blr_character_scalar(trait_id, "trait ID"),
    joint_multitrait = "sblr:joint_multitrait")
  x <- .blr_u64_from_uint32(user_seed)
  x <- .blr_splitmix64(.blr_u64_xor(x, .blr_u64_from_uint32(1)))
  x <- .blr_splitmix64(.blr_u64_xor(
    x, .blr_u64_from_uint32(unname(mode_codes[[analysis_mode]]))))
  x <- .blr_splitmix64(.blr_u64_xor(x, .blr_fnv1a64_utf8(identity)))
  x <- .blr_splitmix64(.blr_u64_xor(
    x, .blr_u64_from_uint32(chain_index)))
  folded <- .blr_u64_xor(x, .blr_u64_rshift(x, 32))
  folded[[1L]] + .blr_u64_base * folded[[2L]]
}

.blr_logical_task_plan <- function(analysis_mode, trait_ids, chains) {
  analysis_mode <- .blr_character_scalar(
    analysis_mode, "analysis_mode",
    c("single_trait", "independent_traits", "joint_multitrait"))
  trait_ids <- .blr_ids(trait_ids, "trait_ids")
  chains <- as.integer(.blr_scalar_whole(
    chains, "chains", 1, .Machine$integer.max))
  if (analysis_mode == "single_trait" && length(trait_ids) != 1L) {
    stop("single_trait requires exactly one trait ID.", call. = FALSE)
  }
  if (analysis_mode == "independent_traits") {
    trait <- rep(trait_ids, each = chains)
    chain <- rep(seq_len(chains) - 1L, times = length(trait_ids))
  } else {
    trait <- rep(NA_character_, chains)
    chain <- seq_len(chains) - 1L
  }
  task_id <- if (analysis_mode == "independent_traits") {
    paste0("trait:", trait, "|chain:", chain)
  } else paste0("chain:", chain)
  data.frame(
    task_index = seq_along(task_id) - 1L,
    task_id = task_id,
    trait_id = trait,
    chain_index = chain,
    stringsAsFactors = FALSE)
}

.blr_chain_seed_base_uint32 <- function(chain_seeds, chains) {
  if (is.null(chain_seeds)) return(NULL)
  if (!is.numeric(chain_seeds) || length(chain_seeds) != chains ||
      anyNA(chain_seeds) || any(!is.finite(chain_seeds)) ||
      any(chain_seeds != floor(chain_seeds)) ||
      any(chain_seeds < -2147483648 | chain_seeds > 4294967295)) {
    stop(paste0(
      "chain_seeds must supply one signed-int32 or uint32-compatible base ",
      "seed per chain."), call. = FALSE)
  }
  ifelse(chain_seeds < 0, chain_seeds + 4294967296, chain_seeds)
}

.blr_task_seeds_v1 <- function(seed, analysis_mode, trait_ids, chains,
                               chain_seeds = NULL) {
  plan <- .blr_logical_task_plan(analysis_mode, trait_ids, chains)
  bases <- .blr_chain_seed_base_uint32(chain_seeds, chains)
  final <- vapply(seq_len(nrow(plan)), function(index) {
    chain <- plan$chain_index[[index]]
    base <- if (is.null(bases)) seed else bases[[chain + 1L]]
    .blr_seed_v1(base, analysis_mode, plan$trait_id[[index]], chain)
  }, numeric(1))
  chain_ids <- paste0("chain", seq_len(chains))
  if (analysis_mode == "independent_traits") {
    return(matrix(
      final, nrow = length(trait_ids), ncol = chains, byrow = TRUE,
      dimnames = list(trait = trait_ids, chain = chain_ids)))
  }
  stats::setNames(final, chain_ids)
}

.blr_task_seed_vector <- function(spec) {
  validate_blr_resolved_spec(spec)
  seeds <- spec$mcmc$task_seeds
  if (spec$data$analysis_mode == "independent_traits") {
    return(as.numeric(t(seeds)))
  }
  as.numeric(seeds)
}

.blr_retention_plan <- function(burn_in_iterations, sampling_iterations,
                                thin_interval, contract_version = 1L,
                                retained_requested = TRUE) {
  burn <- as.integer(.blr_scalar_whole(
    burn_in_iterations, "burn_in_iterations", 0, .Machine$integer.max))
  sampling <- as.integer(.blr_scalar_whole(
    sampling_iterations, "sampling_iterations", 1, .Machine$integer.max))
  thin <- as.integer(.blr_scalar_whole(
    thin_interval, "thin_interval", 1, .Machine$integer.max))
  version <- as.integer(.blr_scalar_whole(
    contract_version, "retention-contract version", 0, 1))
  post_burn <- if (version == 1L) {
    index <- seq_len(sampling)
    index[index %% thin == 0L]
  } else seq.int(1L, sampling, by = thin)
  if (isTRUE(retained_requested) && !length(post_burn)) {
    stop("Retained output was requested but the retention rule keeps no draws.",
         call. = FALSE)
  }
  list(
    contract_version = version,
    post_burn = as.integer(post_burn),
    absolute_transition = as.integer(burn + post_burn),
    retained_draws = as.integer(length(post_burn)))
}

.blr_convergence_iteration_plan <- function(burn_in_iterations,
                                             sampling_iterations) {
  burn <- as.integer(.blr_scalar_whole(
    burn_in_iterations, "burn_in_iterations", 0, .Machine$integer.max))
  sampling <- as.integer(.blr_scalar_whole(
    sampling_iterations, "sampling_iterations", 1, .Machine$integer.max))
  list(
    post_burn = seq_len(sampling),
    absolute_transition = burn + seq_len(sampling),
    checkpoint = "completed_post_burn_iteration",
    thinning = "unthinned",
    rng_draws = 0L)
}

.blr_native_execution_contract <- function(spec) {
  validate_blr_resolved_spec(spec)
  if (spec$schema$seed_contract_version == 0L &&
      spec$schema$retention_contract_version == 0L &&
      spec$compute$scheduler_version == 0L) return(NULL)
  if (!identical(as.integer(spec$schema$seed_contract_version), 1L) ||
      !identical(as.integer(spec$schema$retention_contract_version), 1L) ||
      !identical(as.integer(spec$compute$scheduler_version), 1L)) {
    stop("Phase 3 activation requires seed, retention, and scheduler version 1.",
         call. = FALSE)
  }
  plan <- .blr_logical_task_plan(
    spec$data$analysis_mode, spec$data$trait_ids, spec$mcmc$chains)
  list(
    seed_contract_version = 1L,
    retention_contract_version = 1L,
    scheduler_version = 1L,
    task_seeds = .blr_task_seed_vector(spec),
    task_ids = plan$task_id,
    retained_transition_indices =
      as.integer(spec$mcmc$retained_transition_indices))
}
