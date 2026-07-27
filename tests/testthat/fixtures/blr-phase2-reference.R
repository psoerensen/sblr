phase2_reference_metadata <- list(
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

phase2_reference_configurations <- list(
  one_trait_one_chain_one_core = list(
    traits = 1L, chains = 1L, cores = 1L, keep_chains = FALSE,
    chain_seeds = NULL, selection_s = NULL
  ),
  one_trait_two_chains_one_core = list(
    traits = 1L, chains = 2L, cores = 1L, keep_chains = FALSE,
    chain_seeds = NULL, selection_s = NULL
  ),
  one_trait_two_chains_two_cores = list(
    traits = 1L, chains = 2L, cores = 2L, keep_chains = FALSE,
    chain_seeds = NULL, selection_s = NULL
  ),
  multiple_traits = list(
    traits = 2L, chains = 2L, cores = 2L, keep_chains = FALSE,
    chain_seeds = NULL, selection_s = NULL
  ),
  explicit_chain_seeds = list(
    traits = 1L, chains = 2L, cores = 1L, keep_chains = FALSE,
    chain_seeds = c(401L, 402L), selection_s = NULL
  ),
  keep_chains = list(
    traits = 1L, chains = 2L, cores = 1L, keep_chains = TRUE,
    chain_seeds = c(401L, 402L), selection_s = NULL
  ),
  fixed_selection_s = list(
    traits = 1L, chains = 1L, cores = 1L, keep_chains = FALSE,
    chain_seeds = NULL, selection_s = -0.5
  )
)

phase2_reference_write_csr <- function() {
  prefix <- tempfile("blr_phase2_reference_csr_")
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

phase2_reference_inputs <- function(traits) {
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

phase2_reference_native <- function(config, prefix, inputs) {
  stats <- inputs$stats
  nt <- length(stats$yy)
  m <- stats$m
  vy <- as.numeric(stats$yy) / (stats$n - 1)
  B <- diag((vy * 0.3) / (m * 0.5), nt, nt)
  E <- diag(vy * 0.7, nt, nt)
  ssb <- diag(((4 - 2) / 4) * (vy * 0.3) / (m * 0.5), nt, nt)
  sse <- diag(((4 - 2) / 4) * (vy * 0.7), nt, nt)
  selection <- sblr:::.stblr_prepare_csr_bayesc_selection_s(
    selection_s = config$selection_s,
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
    selection_s_prior_scale = selection$prior_scale
  )
}

phase2_reference_formatted <- function(config, prefix, inputs) {
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
  if (!is.null(config$selection_s)) {
    args$selection_s <- config$selection_s
    args$allow_reference_maf_for_selection_s <- TRUE
  }
  do.call(stblr_csr, args)
}

phase2_reference_normalize <- function(x) {
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

phase2_reference_md5 <- function(x) {
  path <- tempfile("blr_phase2_reference_", fileext = ".rds")
  on.exit(unlink(path), add = TRUE)
  saveRDS(phase2_reference_normalize(x), path, version = 3, compress = FALSE)
  unname(tools::md5sum(path))
}

phase2_reference_objects <- function(config) {
  prefix <- phase2_reference_write_csr()
  inputs <- phase2_reference_inputs(config$traits)
  list(
    raw = phase2_reference_native(config, prefix, inputs),
    fit = phase2_reference_formatted(config, prefix, inputs)
  )
}

phase2_reference_hashes <- function() {
  lapply(phase2_reference_configurations, function(config) {
    objects <- phase2_reference_objects(config)
    c(raw = phase2_reference_md5(objects$raw),
      fit = phase2_reference_md5(objects$fit))
  })
}
