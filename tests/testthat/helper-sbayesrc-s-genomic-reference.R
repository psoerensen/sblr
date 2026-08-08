# Compact internal-only Phase-4B CSR qualification helpers.

.sbs4b_csr_prefix <- function(marker_count) {
  prefix <- tempfile("sbs4b_csr_")
  sblr:::.stblr_write_uint64_file(
    paste0(prefix, ".row_ptr.u64.bin"), rep(0, marker_count + 1L)
  )
  file.create(paste0(prefix, ".col_idx.u32.0based.bin"))
  file.create(paste0(prefix, ".values.f32.bin"))
  writeLines(c(
    "format=sparse_ld_csr", "storage=streamed_upper_triangle",
    "n_bed=NA", "n_used=NA", "n_samples_used=NA",
    paste0("n_variants=", marker_count), "nnz=0", "triangle=upper",
    "diagonal=implicit_1", paste0("row_ptr_file=", prefix, ".row_ptr.u64.bin"),
    paste0("col_idx_file=", prefix, ".col_idx.u32.0based.bin"),
    paste0("values_file=", prefix, ".values.f32.bin"),
    "row_ptr_type=uint64", "col_idx_type=uint32", "values_type=float32",
    "index_base=0", "value=r"
  ), paste0(prefix, ".meta.txt"))
  prefix
}

.sbs4b_fixture <- function(marker_count = 80L, seed = 20270899L) {
  set.seed(seed)
  annotation <- cbind(
    enriched = stats::rbinom(marker_count, 1L, 0.25),
    continuous = as.numeric(scale(stats::rnorm(marker_count))),
    null = as.numeric(scale(stats::rnorm(marker_count)))
  )
  A <- cbind(Intercept = 1, annotation)
  gamma <- c(0, 0.03, 0.2, 1)
  component <- rep(0:3, length.out = marker_count)
  beta <- numeric(marker_count)
  beta[component > 0L] <- stats::rnorm(sum(component > 0L), 0, 0.16)
  beta[component == 3L] <- beta[component == 3L] +
    0.18 * (2 * annotation[component == 3L, 1L] - 1)
  sample_size <- 120L
  score <- sample_size * beta + stats::rnorm(marker_count, 0, 1.5)
  intercept_prior_specification <- sblr:::.sbayesrc_resolve_intercept_prior(
    rep(0.25, 4L), annotation_intercept_prior = list(
      mean = "initial_mixture", sd = 1
    )
  )
  intercept_prior <- intercept_prior_specification$native
  intercept_prior <- rbind(
    intercept_prior,
    update_sigmaSqAlpha = rep(1, 3L),
    allocation_updates_per_cycle = rep(1, 3L),
    annotation_updates_per_cycle = rep(1, 3L)
  )
  list(
    wy = list(score), ww = list(rep(sample_size, marker_count)),
    yy = sample_size, b_init = list(beta), comp_init = list(component),
    r_init = list(score - sample_size * beta), prefix = .sbs4b_csr_prefix(marker_count),
    A = A, annotation = annotation, gamma = gamma,
    alpha_init = matrix(0, 4L, 3L), delta_init = c(1L, 1L, 0L),
    tau2_init = rep(0.8, 3L), intercept_prior = intercept_prior,
    intercept_prior_specification = intercept_prior_specification,
    B = matrix(0.08, 1L, 1L), E = matrix(1, 1L, 1L),
    ssb_prior = list(0.08), sse_prior = list(1), n = sample_size
  )
}

.sbs4b_run <- function(fixture, seed = 20270900L, nit = 800L,
                       nburn = 200L, fixed_delta = integer(),
                       updateB = FALSE, updateE = FALSE,
                       initial_delta = fixture$delta_init,
                       update_hierarchy = TRUE,
                       update_pi_A = TRUE, update_tau2 = TRUE,
                       hierarchy_sweeps = 1L, genomic_sweeps = 1L) {
  intercept_prior <- fixture$intercept_prior
  intercept_prior["allocation_updates_per_cycle", ] <- genomic_sweeps
  intercept_prior["annotation_updates_per_cycle", ] <- hierarchy_sweeps
  .st_sbayesrc_selection_csr(
    fixture$wy, fixture$ww, fixture$yy, fixture$b_init,
    fixture$comp_init, TRUE, fixture$r_init, TRUE, fixture$prefix,
    fixture$B, fixture$E, fixture$ssb_prior, fixture$sse_prior,
    fixture$A, fixture$gamma, fixture$alpha_init,
    as.integer(initial_delta), 0.35, fixture$tau2_init,
    1, 1, 3, 1.6, as.integer(fixed_delta), update_hierarchy,
    update_pi_A, update_tau2, intercept_prior,
    1e-12, 4, 4, updateB, updateE, 0.9, fixture$n,
    nit, nburn, 1L, 1L, seed, 1L, seed
  )
}

.sbs4b_run_standard <- function(fixture, seed = 20270900L, nit = 800L,
                                nburn = 200L, update_hierarchy = TRUE) {
  stblr_cpg_omp_csr_sbayesrc(
    fixture$wy, fixture$ww, fixture$yy, fixture$b_init,
    fixture$comp_init, TRUE, fixture$r_init, TRUE, FALSE, fixture$prefix,
    fixture$B, fixture$E, fixture$ssb_prior, fixture$sse_prior,
    fixture$A, fixture$gamma, fixture$alpha_init, fixture$tau2_init,
    fixture$intercept_prior, 6, 3.2, 1e-12, 4, 4,
    update_hierarchy, FALSE, FALSE, 1L, 0.9, fixture$n, nit, nburn, 1L, 1L,
    seed, 1L, TRUE, seed, FALSE, 0, 0.8, 50L, 0L,
    NULL, FALSE, 0, c(-3, 2), 0.35, NULL, integer(), TRUE,
    FALSE, FALSE, TRUE
  )
}
