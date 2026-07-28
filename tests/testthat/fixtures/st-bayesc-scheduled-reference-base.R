st_bayesc_scheduled_base_starting_commit <- "ab18b9c"

st_bayesc_scheduled_base_prefix <- function() {
  prefix <- tempfile("st_bayesc_scheduled_base_scheduled_")
  row_ptr <- c(0, 2, 3, 4, 4, 5, 5)
  col_idx <- c(1L, 2L, 2L, 3L, 5L)
  values <- c(.65, -.25, .45, .55, -.35)
  sblr:::.stblr_write_uint64_file(paste0(prefix, ".row_ptr.u64.bin"), row_ptr)
  sblr:::.stblr_write_uint32_file(paste0(prefix, ".col_idx.u32.0based.bin"), col_idx)
  writeBin(values, paste0(prefix, ".values.f32.bin"), size = 4, endian = "little")
  writeLines(c(
    "format=sparse_ld_csr", "storage=streamed_upper_triangle", "n_bed=NA",
    "n_used=NA", "n_samples_used=NA", "n_variants=6", "nnz=5",
    "triangle=upper", "diagonal=implicit_1",
    paste0("row_ptr_file=", prefix, ".row_ptr.u64.bin"),
    paste0("col_idx_file=", prefix, ".col_idx.u32.0based.bin"),
    paste0("values_file=", prefix, ".values.f32.bin"),
    "row_ptr_type=uint64", "col_idx_type=uint32", "values_type=float32",
    "index_base=0", "value=r"
  ), paste0(prefix, ".meta.txt"))
  prefix
}

st_bayesc_scheduled_base_stats <- function() {
  markers <- paste0("m", 1:6)
  list(
    wy = list(trait1 = stats::setNames(c(4, -2, .25, 1.5, -.1, .8), markers)),
    ww = list(trait1 = stats::setNames(rep(80, 6), markers)),
    yy = stats::setNames(80, "trait1"), n = 80L, m = 6L,
    marker_names = markers, trait_names = "trait1"
  )
}

st_bayesc_scheduled_base_configs <- list(
  dense_one = list(nchains=1L,ncores=1L,seeds=NULL,full=1L,base=1L,max=1L,
                   threshold=0,lifetime=0L,burnin_only=FALSE,wakeup=FALSE),
  skip_two_one = list(nchains=2L,ncores=1L,seeds=c(1101L,1102L),full=3L,
                      base=2L,max=7L,threshold=.12,lifetime=2L,
                      burnin_only=FALSE,wakeup=TRUE),
  skip_two_two = list(nchains=2L,ncores=2L,seeds=c(1201L,1202L),full=4L,
                      base=3L,max=9L,threshold=.08,lifetime=3L,
                      burnin_only=TRUE,wakeup=TRUE)
)

st_bayesc_scheduled_base_normalize <- function(x) {
  if (!is.list(x)) return(x)
  nms <- names(x)
  for (i in seq_along(x)) {
    nm <- if (is.null(nms)) NA_character_ else nms[[i]]
    if (is.na(nm) || !nzchar(nm)) x[i] <- list(st_bayesc_scheduled_base_normalize(x[[i]]))
    else if (nm %in% c("seconds", "seconds_mean", "seconds_max")) x[[nm]][] <- 0
    else if (identical(nm, "ld_prefix")) x[[nm]] <- "<fixture>"
    else x[nm] <- list(st_bayesc_scheduled_base_normalize(x[[nm]]))
  }
  x
}

st_bayesc_scheduled_base_run <- function(cfg) {
  captured <- NULL
  ns <- asNamespace("sblr")
  suppressMessages(trace(".as_stblr_fit", where=ns,
    tracer=quote(assign(".st_bayesc_scheduled_base_raw_capture", raw, envir=.GlobalEnv)),
    print=FALSE))
  on.exit(suppressMessages(untrace(".as_stblr_fit", where=ns)), add=TRUE)
  fit <- sblr::stblr_csr(
    stats=st_bayesc_scheduled_base_stats(), ld_prefix=st_bayesc_scheduled_base_prefix(), scheduled=TRUE,
    pi_init=.35, pi_prior_mean=.35, pi_prior_strength=3,
    updateB=FALSE, updateE=FALSE, updatePi=FALSE,
    nit=8L, nburn=2L, nthin=1L, seed=1001L,
    nchains=cfg$nchains, ncores=cfg$ncores, keep_chains=FALSE,
    convergence="none",
    chain_seeds=cfg$seeds, full_sweep_every=cfg$full,
    null_skip_base=cfg$base, null_skip_max=cfg$max,
    candidate_threshold=cfg$threshold, candidate_lifetime=cfg$lifetime,
    skip_nulls_burnin_only=cfg$burnin_only,
    wakeup_ld_neighbors=cfg$wakeup, wakeup_diff_threshold=.02,
    wakeup_max_neighbors=2L, updateLDswap=FALSE
  )
  captured <- get(".st_bayesc_scheduled_base_raw_capture", envir=.GlobalEnv)
  rm(".st_bayesc_scheduled_base_raw_capture", envir=.GlobalEnv)
  list(raw=st_bayesc_scheduled_base_normalize(captured), fit=st_bayesc_scheduled_base_normalize(fit))
}

st_bayesc_scheduled_base_metadata <- function(name,cfg) list(
  starting_commit=st_bayesc_scheduled_base_starting_commit, model="scheduled CSR BayesC",
  configuration=name, marker_count=6L, trait_count=1L, chains=cfg$nchains,
  cores=cfg$ncores, chain_seeds=cfg$seeds, iterations=8L, burnin=2L,
  thinning=1L, controls=cfg, schema="stblr_raw_v1",
  reference_mode="fresh R process", compiler="Rtools44 GCC 13.2 C++17"
)
