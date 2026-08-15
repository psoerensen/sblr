st_bayesc_csr_reference_metadata <- list(
  starting_commit = "7df2d1ddec1c5f45be03e399dc7afb419afd56d8",
  r_version = "R 4.4.1",
  toolchain = "Rtools44 GNU C++17",
  fixture = paste(
    "three-marker CSR with one 0.4 correlation between m1 and m2;",
    "wy trait1 = c(25, 15, 4), trait2 = c(12, 20, 6);",
    "ww = 100; yy = 100; n = 100"
  ),
  seed = 31L,
  iterations = 8L,
  burn_in = 2L,
  thinning = 1L,
  expected_schema_version = 1L
)

st_bayesc_csr_reference_configurations <- list(
  one_trait_one_chain_one_core = list(
    traits = 1L, chains = 1L, cores = 1L, keep_chains = FALSE,
    chain_seeds = NULL, maf_effect_s = NULL
  ),
  one_trait_two_chains_one_core = list(
    traits = 1L, chains = 2L, cores = 1L, keep_chains = FALSE,
    chain_seeds = NULL, maf_effect_s = NULL
  ),
  one_trait_two_chains_two_cores = list(
    traits = 1L, chains = 2L, cores = 2L, keep_chains = FALSE,
    chain_seeds = NULL, maf_effect_s = NULL
  ),
  multiple_traits = list(
    traits = 2L, chains = 2L, cores = 2L, keep_chains = FALSE,
    chain_seeds = NULL, maf_effect_s = NULL
  ),
  explicit_chain_seeds = list(
    traits = 1L, chains = 2L, cores = 1L, keep_chains = FALSE,
    chain_seeds = c(401L, 402L), maf_effect_s = NULL
  ),
  keep_chains = list(
    traits = 1L, chains = 2L, cores = 1L, keep_chains = TRUE,
    chain_seeds = c(401L, 402L), maf_effect_s = NULL
  ),
  fixed_maf_effect_s = list(
    traits = 1L, chains = 1L, cores = 1L, keep_chains = FALSE,
    chain_seeds = NULL, maf_effect_s = -0.5
  )
)

# Captured from the pre-extraction DLL on 2026-08-09. Timing and temporary LD
# prefixes are normalized by st_bayesc_csr_reference_md5(). These hashes freeze
# the complete accessible raw/formatted trajectory and schema across ordinary,
# multichain, explicit-seed, retained-chain, and fixed-scale routes.
st_bayesc_csr_pre_engine_extraction_hashes <- list(
  one_trait_one_chain_one_core = c(
    raw = "1bdbead1dbdc2f3b5f9c3576d2f78669",
    fit = "fb01c4ef85b01cfecef4c7b0297a2ac2"
  ),
  one_trait_two_chains_one_core = c(
    raw = "6c5667d2ab9a1aa9fdd1e4d3bd289122",
    fit = "7e546722842b9d59f283a4a261cd8522"
  ),
  one_trait_two_chains_two_cores = c(
    raw = "6c5667d2ab9a1aa9fdd1e4d3bd289122",
    fit = "34fb0eea180f5378f0cc6ae4773a29a9"
  ),
  multiple_traits = c(
    raw = "e63936437a33fd3f2b479401a2520fac",
    fit = "a9ba3f3a5e1a38bd6a5de6c2e45029fc"
  ),
  explicit_chain_seeds = c(
    raw = "fc9847ed816c0cc2e41b1ded38dccc0a",
    fit = "68bd9a521469fab766640478eb76d9f2"
  ),
  keep_chains = c(
    raw = "25617ccc7903fe7e28116871593cb495",
    fit = "b258ac45a7d4d30914305b210047d130"
  ),
  fixed_maf_effect_s = c(
    raw = "1c50036f9bca4b12456d97d333bbedc0",
    fit = "8de361541b489facc3ae62d5dfeaa631"
  )
)

st_bayesc_csr_reference_write_csr <- function() {
  prefix <- tempfile("blr_st_bayesc_csr_reference_csr_")
  sblr:::.stblr_write_uint64_file(
    paste0(prefix, ".row_ptr.u64.bin"), c(0, 1, 1, 1)
  )
  sblr:::.stblr_write_uint32_file(
    paste0(prefix, ".col_idx.u32.0based.bin"), 1L
  )
  writeBin(
    0.4, paste0(prefix, ".values.f32.bin"), size = 4, endian = "little"
  )
  writeLines(c(
    "format=sparse_ld_csr", "storage=streamed_upper_triangle",
    "n_bed=NA", "n_used=NA", "n_samples_used=NA", "n_variants=3",
    "nnz=1", "triangle=upper", "diagonal=implicit_1",
    paste0("row_ptr_file=", prefix, ".row_ptr.u64.bin"),
    paste0("col_idx_file=", prefix, ".col_idx.u32.0based.bin"),
    paste0("values_file=", prefix, ".values.f32.bin"),
    "row_ptr_type=uint64", "col_idx_type=uint32", "values_type=float32",
    "index_base=0", "value=r"
  ), paste0(prefix, ".meta.txt"))
  prefix
}

st_bayesc_csr_reference_inputs <- function(traits) {
  marker_ids <- c("m1", "m2", "m3")
  trait_ids <- paste0("trait", seq_len(traits))
  wy_values <- list(c(25, 15, 4), c(12, 20, 6))[seq_len(traits)]
  wy <- lapply(wy_values, stats::setNames, nm = marker_ids)
  ww <- replicate(
    traits, stats::setNames(rep(100, 3), marker_ids), simplify = FALSE
  )
  names(wy) <- names(ww) <- trait_ids
  stats <- list(
    wy = wy, ww = ww, yy = stats::setNames(rep(100, traits), trait_ids),
    n = 100L, m = 3L, marker_names = marker_ids, trait_names = trait_ids
  )
  glist <- list(
    rsidsLD = list(marker_ids),
    rsids = list(c("m3", "m1", "m2")),
    maf = list(c(0.40, 0.05, 0.20))
  )
  list(stats = stats, glist = glist)
}

st_bayesc_csr_reference_native <- function(config, prefix, inputs,
                                            execution_contract = NULL) {
  stats <- inputs$stats
  nt <- length(stats$yy)
  m <- stats$m
  vy <- as.numeric(stats$yy) / (stats$n - 1)
  B <- diag((vy * 0.3) / (m * 0.5), nt, nt)
  E <- diag(vy * 0.7, nt, nt)
  ssb <- diag(((4 - 2) / 4) * (vy * 0.3) / (m * 0.5), nt, nt)
  sse <- diag(((4 - 2) / 4) * (vy * 0.7), nt, nt)
  selection <- sblr:::.stblr_prepare_csr_bayesc_maf_effect_s(
    maf_effect_s = config$maf_effect_s,
    Glist = inputs$glist,
    m = m,
    scheduled = FALSE,
    return_log_h = FALSE
  )
  sblr:::stblr_cpg_omp_csr(
    wy = stats$wy, ww = stats$ww, yy = stats$yy,
    b_init = replicate(nt, rep(0, m), simplify = FALSE),
    d_init = replicate(nt, rep(0, m), simplify = FALSE),
    use_d_init = FALSE, r_init = stats$wy, use_r_init = FALSE,
    rebuild_r_before_updateE = FALSE, ld_prefix = prefix,
    B = B, E = E,
    ssb_prior = split(ssb, rep(seq_len(nt), each = nt)),
    sse_prior = split(sse, rep(seq_len(nt), each = nt)),
    pi = c(0.5, 0.5), nub = 4, nue = 4,
    updateB = FALSE, updateE = FALSE, updatePi = TRUE, adjE = 0.9,
    n = rep(stats$n, nt), nit = 8L, nburn = 2L, nthin = 1L,
    pi_prior_a = 1, pi_prior_b = 1,
    ncores = config$cores, seed = 31L, nchains = config$chains,
    keep_chains = config$keep_chains,
    chain_seeds = if (is.null(config$chain_seeds)) integer() else config$chain_seeds,
    updateLDswap = FALSE,
    maf_effect_s_prior_scale = selection$prior_scale,
    execution_contract = execution_contract
  )
}

st_bayesc_csr_reference_formatted <- function(config, prefix, inputs) {
  args <- list(
    stats = inputs$stats,
    Glist = inputs$glist,
    ld_prefix = prefix,
    method = "sbayesc",
    pi_init = 0.5,
    pi_prior_mean = 0.5,
    pi_prior_strength = 2,
    updateB = FALSE,
    updateE = FALSE,
    updatePi = TRUE,
    nit = 8L,
    nburn = 2L,
    nthin = 1L,
    seed = 31L,
    nchains = config$chains,
    keep_chains = config$keep_chains,
    chain_seeds = config$chain_seeds,
    ncores = config$cores,
    updateLDswap = FALSE,
    scheduled = FALSE
  )
  if (!is.null(config$maf_effect_s)) {
    args$maf_effect_s <- config$maf_effect_s
    args$allow_reference_maf_for_maf_effect_s <- TRUE
  }
  do.call(stblr_csr, args)
}

st_bayesc_csr_reference_normalize <- function(x) {
  normalize <- function(value) {
    if (!is.list(value)) return(value)
    for (name in names(value)) {
      if (name %in% c("seconds_mean", "seconds_max")) {
        value[[name]] <- rep(0, length(value[[name]]))
      } else if (identical(name, "ld_prefix")) {
        value[[name]] <- "<fixture-csr-prefix>"
      } else {
        value[[name]] <- normalize(value[[name]])
      }
    }
    value
  }
  normalize(x)
}

st_bayesc_csr_reference_md5 <- function(x) {
  path <- tempfile("blr_st_bayesc_csr_reference_", fileext = ".rds")
  on.exit(unlink(path), add = TRUE)
  saveRDS(st_bayesc_csr_reference_normalize(x), path, version = 3, compress = FALSE)
  unname(tools::md5sum(path))
}

st_bayesc_csr_reference_objects <- function(config) {
  prefix <- st_bayesc_csr_reference_write_csr()
  inputs <- st_bayesc_csr_reference_inputs(config$traits)
  list(
    raw = st_bayesc_csr_reference_native(config, prefix, inputs),
    fit = st_bayesc_csr_reference_formatted(config, prefix, inputs)
  )
}

st_bayesc_csr_reference_hashes <- function() {
  lapply(st_bayesc_csr_reference_configurations, function(config) {
    objects <- st_bayesc_csr_reference_objects(config)
    c(raw = st_bayesc_csr_reference_md5(objects$raw),
      fit = st_bayesc_csr_reference_md5(objects$fit))
  })
}
