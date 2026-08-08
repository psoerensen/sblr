# Development qualification for the internal Phase-5A MCEM-SBayesRC path.
# Not part of the supported public API.

if (!file.exists("local_reference/sbayesrc_mcem_exact_validation.R")) {
  stop("Run from the sblr repository root with local_reference available.")
}

devtools::load_all(".", quiet = TRUE)
source("local_reference/sbayesrc_mcem_exact_validation.R", local = FALSE)

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
