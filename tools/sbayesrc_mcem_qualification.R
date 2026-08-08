# Development qualification for internal Phase-5A/5B SBayesRC-EM paths.
# Not part of the supported public API.

if (!file.exists("local_reference/sbayesrc_mcem_exact_validation.R")) {
  stop("Run from the sblr repository root with local_reference available.")
}

devtools::load_all(".", quiet = TRUE)
source("local_reference/sbayesrc_mcem_exact_validation.R", local = FALSE)
source("tests/testthat/helper-sbayesrc-mcem-reference.R", local = FALSE)

write_csr <- function(correlation) {
  prefix <- tempfile("mcem_phase5a_")
  marker_count <- nrow(correlation)
  edge <- which(upper.tri(correlation) & correlation != 0, arr.ind = TRUE)
  edge <- edge[order(edge[, 1L], edge[, 2L]), , drop = FALSE]
  row_ptr <- c(0, cumsum(tabulate(edge[, 1L], nbins = marker_count)))
  sblr:::.stblr_write_uint64_file(
    paste0(prefix, ".row_ptr.u64.bin"), row_ptr
  )
  sblr:::.stblr_write_uint32_file(
    paste0(prefix, ".col_idx.u32.0based.bin"), edge[, 2L] - 1L
  )
  writeBin(as.numeric(correlation[edge]), paste0(prefix, ".values.f32.bin"),
           size = 4L, endian = .Platform$endian)
  writeLines(c(
    "format=sparse_ld_csr", "storage=streamed_upper_triangle",
    "n_bed=NA", "n_used=NA", "n_samples_used=NA",
    paste0("n_variants=", marker_count), paste0("nnz=", nrow(edge)),
    "triangle=upper", "diagonal=implicit_1",
    paste0("row_ptr_file=", prefix, ".row_ptr.u64.bin"),
    paste0("col_idx_file=", prefix, ".col_idx.u32.0based.bin"),
    paste0("values_file=", prefix, ".values.f32.bin"),
    "row_ptr_type=uint64", "col_idx_type=uint32", "values_type=float32",
    "index_base=0", "value=r"
  ), paste0(prefix, ".meta.txt"))
  prefix
}

cleanup_csr <- function(prefix) {
  unlink(paste0(prefix, c(
    ".row_ptr.u64.bin", ".col_idx.u32.0based.bin", ".values.f32.bin",
    ".meta.txt"
  )))
}

run_native_mcem <- function(data, alpha_start, retained, burn, max_outer, seed) {
  diagonal <- diag(data$XtX)
  correlation <- data$XtX / sqrt(base::outer(diagonal, diagonal))
  diag(correlation) <- 1
  prefix <- write_csr(correlation)
  on.exit(cleanup_csr(prefix), add = TRUE)
  intercept_prior <- rbind(
    type = rep(0, ncol(alpha_start)),
    mean = data$intercept_prior_mean,
    precision = rep(1, ncol(alpha_start))
  )
  marker_count <- ncol(data$X)
  sblr:::.stblr_mcem_sbayesrc_csr(
    wy = list(data$Xty), ww = list(diagonal), yy = data$yty,
    b_init = list(rep(0, marker_count)),
    comp_init = list(rep(0, marker_count)), r_init = list(data$Xty),
    ld_prefix = prefix, B = matrix(data$sigma2_beta, 1L, 1L),
    E = matrix(data$sigma2_e, 1L, 1L),
    ssb_prior = list(data$sigma2_beta), sse_prior = list(data$sigma2_e),
    A = data$A, gamma = data$gamma, alpha_init = alpha_start,
    sigmaSqAlpha_init = rep(1, ncol(alpha_start)),
    intercept_prior_resolved = intercept_prior, n = nrow(data$X),
    inner_sweeps = retained, inner_burn = burn,
    final_sweeps = max(retained, 1500L), final_burn = max(burn, 500L),
    damping = 0.5, min_outer = 3L, max_outer = max_outer,
    seed = seed, ncores = 1L
  )
}

quick <- identical(Sys.getenv("SBLR_MCEM_QUICK"), "1")

step1_data <- simulate_step1_case(2026080801L)
step1_starts <- make_starting_alphas(step1_data$A, step1_data$baseline_pi)
step1_log_bf <- independent_log_bf(step1_data)
step1_direct <- direct_optimize_orthogonal(
  step1_data, step1_starts, step1_log_bf
)$best
step1_fit <- run_native_mcem(
  step1_data, step1_starts$baseline,
  retained = if (quick) 300L else 1000L,
  burn = if (quick) 100L else 500L,
  max_outer = if (quick) 12L else 50L,
  seed = 20271301L
)
step1_exact_responsibility <- orthogonal_exact_responsibilities(
  step1_fit$mcem$alpha_map, step1_data, step1_log_bf
)

step2_data <- simulate_step2_case(2026080802L)
step2_starts <- make_starting_alphas(step2_data$A, step2_data$baseline_pi)
step2_exact_model <- precompute_exact_ld_model(step2_data)
step2_direct <- direct_optimize_exact_ld(
  step2_data, step2_starts, step2_exact_model
)$best
step2_fit <- run_native_mcem(
  step2_data, step2_starts$baseline,
  retained = if (quick) 600L else 2800L,
  burn = if (quick) 200L else 1200L,
  max_outer = if (quick) 15L else 50L,
  seed = 20271302L
)
step2_exact_responsibility <- exact_ld_responsibilities(
  step2_fit$mcem$alpha_map, step2_data, step2_exact_model
)

summarize <- function(label, fit, direct, exact_responsibility) {
  data.frame(
    case = label,
    n_outer = fit$mcem$n_outer,
    converged = fit$mcem$converged,
    alpha_rmse = sqrt(mean((fit$mcem$alpha_map - direct$alpha)^2)),
    alpha_max_abs = max(abs(fit$mcem$alpha_map - direct$alpha)),
    prior_rmse = sqrt(mean((fit$mcem$component_prior - direct$prior)^2)),
    responsibility_rmse = sqrt(mean((
      fit$mcem$final_responsibilities - exact_responsibility
    )^2)),
    max_delta_alpha = tail(fit$mcem$history$summary$max_delta_alpha, 1L),
    max_delta_prior = tail(fit$mcem$history$summary$max_delta_prior, 1L),
    stringsAsFactors = FALSE
  )
}

summary <- rbind(
  summarize("orthogonal", step1_fit, step1_direct, step1_exact_responsibility),
  summarize("correlated_ld", step2_fit, step2_direct, step2_exact_responsibility)
)
print(summary, digits = 6, row.names = FALSE)

output_dir <- file.path("results", "local", "sbayesrc_mcem", "phase5A")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(summary, file.path(output_dir, "reference_parity.csv"), row.names = FALSE)
saveRDS(list(summary = summary, step1 = step1_fit$mcem,
             step2 = step2_fit$mcem),
        file.path(output_dir, "reference_parity.rds"))

# Phase 5B block-eigen gates. The shared deterministic fixture is generated
# once; only backend and qualified hyperparameter update modes vary.
block_fixture <- .mcem_block_fixture(
  seed = 2026081501L,
  marker_count = if (quick) 10L else 16L,
  sample_count = if (quick) 70L else 120L
)
block_baseline <- matrix(c(-0.6, -0.5, 0, 0), 2L, 2L)
block_dispersed <- matrix(c(0.2, 0.1, -0.6, 0.5), 2L, 2L)
block_retained <- if (quick) 300L else 900L
block_burn <- if (quick) 100L else 350L
block_outer <- if (quick) 12L else 35L

run_fixture_csr <- function(alpha_start, seed) {
  sblr:::.stblr_mcem_sbayesrc_csr(
    wy = block_fixture$stats$wy, ww = block_fixture$stats$ww,
    yy = block_fixture$stats$yy,
    b_init = list(rep(0, nrow(block_fixture$A))),
    comp_init = list(rep(0L, nrow(block_fixture$A))),
    r_init = block_fixture$stats$wy, ld_prefix = block_fixture$prefix,
    B = block_fixture$B, E = block_fixture$E,
    ssb_prior = block_fixture$ssb_prior,
    sse_prior = block_fixture$sse_prior,
    A = block_fixture$A, gamma = block_fixture$gamma,
    alpha_init = alpha_start, sigmaSqAlpha_init = c(1, 1),
    intercept_prior_resolved = block_fixture$intercept_prior,
    n = block_fixture$stats$n, inner_sweeps = block_retained,
    inner_burn = block_burn, final_sweeps = block_retained + 200L,
    final_burn = block_burn, max_outer = block_outer,
    seed = seed, ncores = 1L
  )
}

run_fixture_block <- function(alpha_start, seed, updateB = FALSE,
                              updateE = FALSE, retained = block_retained) {
  .mcem_run_block(
    block_fixture, alpha_start, seed, updateB = updateB, updateE = updateE,
    inner_sweeps = retained, inner_burn = block_burn,
    max_outer = block_outer
  )
}

csr_fixed <- run_fixture_csr(block_baseline, 20271501L)
block_fixed <- run_fixture_block(block_baseline, 20271502L)
block_fixed_start <- run_fixture_block(block_dispersed, 20271503L)
block_e <- run_fixture_block(block_baseline, 20271504L, updateE = TRUE)
block_e_start <- run_fixture_block(block_dispersed, 20271505L, updateE = TRUE)
block_b <- run_fixture_block(block_baseline, 20271506L, updateB = TRUE)
block_b_start <- run_fixture_block(block_dispersed, 20271507L, updateB = TRUE)
block_be <- run_fixture_block(
  block_baseline, 20271508L, updateB = TRUE, updateE = TRUE
)
block_be_start <- run_fixture_block(
  block_dispersed, 20271509L, updateB = TRUE, updateE = TRUE
)

gate_row <- function(gate, fit, reference, start_fit, updateB, updateE,
                     alpha_limit, prior_limit, rb_limit) {
  tail_mean <- function(x) mean(tail(x, min(10L, length(x))))
  alpha_difference <- max(abs(fit$mcem$alpha_map - reference$mcem$alpha_map))
  prior_difference <- max(abs(
    fit$mcem$component_prior - reference$mcem$component_prior
  ))
  rb_difference <- max(abs(
    fit$mcem$last_estep_responsibilities -
      reference$mcem$last_estep_responsibilities
  ))
  start_difference <- max(abs(
    fit$mcem$alpha_map - start_fit$mcem$alpha_map
  ))
  B_tail <- tail_mean(fit$mcem$history$summary$B)
  B_start_tail <- tail_mean(start_fit$mcem$history$summary$B)
  E_tail <- tail_mean(fit$mcem$history$summary$E)
  E_start_tail <- tail_mean(start_fit$mcem$history$summary$E)
  data.frame(
    gate = gate, backend = fit$mcem$backend,
    B_mode = if (updateB) "learned" else "fixed",
    E_mode = if (updateE) "learned_gctb_allMixVe" else "fixed_gctb_fixVe",
    sigmaSqAlpha_mode = fit$mcem$sigmaSqAlpha_mode,
    converged = fit$mcem$converged, n_outer = fit$mcem$n_outer,
    alpha_reference_max = alpha_difference,
    prior_reference_max = prior_difference,
    rb_reference_max = rb_difference,
    start_alpha_max = start_difference,
    B_final = fit$mcem$genomic_hyperparameters$B_final[1L, 1L],
    E_final = fit$mcem$genomic_hyperparameters$E_final[1L, 1L],
    B_tail_mean = B_tail, B_start_tail_mean = B_start_tail,
    E_tail_mean = E_tail, E_start_tail_mean = E_start_tail,
    B_final_block_mean = mean(as.numeric(fit$genomic$trace$vbs)),
    B_start_final_block_mean = mean(as.numeric(start_fit$genomic$trace$vbs)),
    E_final_block_mean = mean(as.numeric(fit$genomic$trace$ves)),
    E_start_final_block_mean = mean(as.numeric(start_fit$genomic$trace$ves)),
    pass = alpha_difference < alpha_limit &&
      prior_difference < prior_limit && rb_difference < rb_limit &&
      start_difference < 0.18,
    stringsAsFactors = FALSE
  )
}

gate_summary <- rbind(
  gate_row("G1", block_fixed, csr_fixed, block_fixed_start,
           FALSE, FALSE, 0.12, 0.04, 0.08),
  gate_row("G2", block_e, block_fixed, block_e_start,
           FALSE, TRUE, 0.35, 0.12, 0.20),
  gate_row("G3_B", block_b, block_fixed, block_b_start,
           TRUE, FALSE, 0.35, 0.12, 0.20),
  gate_row("G3_BE", block_be, block_fixed, block_be_start,
           TRUE, TRUE, 0.40, 0.15, 0.22)
)
print(gate_summary, digits = 5, row.names = FALSE)

inner_lengths <- if (quick) c(150L, 300L, 600L) else c(250L, 500L, 1000L)
inner_fits <- lapply(seq_along(inner_lengths), function(index) {
  run_fixture_block(
    block_baseline, 20271600L + index, updateB = TRUE, updateE = TRUE,
    retained = inner_lengths[index]
  )
})
inner_reference <- inner_fits[[length(inner_fits)]]
inner_summary <- do.call(rbind, lapply(seq_along(inner_fits), function(index) {
  fit <- inner_fits[[index]]
  data.frame(
    retained = inner_lengths[index], converged = fit$mcem$converged,
    n_outer = fit$mcem$n_outer,
    alpha_max_vs_long = max(abs(
      fit$mcem$alpha_map - inner_reference$mcem$alpha_map
    )),
    prior_max_vs_long = max(abs(
      fit$mcem$component_prior - inner_reference$mcem$component_prior
    )),
    B_final = fit$mcem$genomic_hyperparameters$B_final[1L, 1L],
    E_final = fit$mcem$genomic_hyperparameters$E_final[1L, 1L],
    stringsAsFactors = FALSE
  )
}))
print(inner_summary, digits = 5, row.names = FALSE)

output_dir_5b <- file.path("results", "local", "sbayesrc_mcem", "phase5B")
dir.create(output_dir_5b, recursive = TRUE, showWarnings = FALSE)
write.csv(gate_summary, file.path(output_dir_5b, "gate_summary.csv"),
          row.names = FALSE)
write.csv(inner_summary, file.path(output_dir_5b, "inner_length_summary.csv"),
          row.names = FALSE)
saveRDS(
  list(gates = gate_summary, inner_length = inner_summary),
  file.path(output_dir_5b, "qualification_summary.rds")
)
.mcem_cleanup_block_fixture(block_fixture)
