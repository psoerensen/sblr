# Development qualification for SBayesRC-S-EM Phase 5C.
# Not a production entry point and not part of the supported public API.

suppressPackageStartupMessages(devtools::load_all(compile = FALSE, quiet = TRUE))
source("tests/testthat/helper-sbayesrc-mcem-reference.R")
source("tests/testthat/helper-sbayesrc-s-mcem-reference.R")

max_difference <- function(x, y) max(abs(x - y))
outer_summary <- function(fit, label, start, mode) {
 data.frame(
  gate = label,
  backend = fit$mcem$backend,
  B_mode = if (fit$mcem$genomic_hyperparameters$updateB) "learned" else "fixed",
  E_mode = if (fit$mcem$genomic_hyperparameters$updateE) "learned" else "fixed",
  pi_A_mode = fit$mcem$pi_A_mode,
  tau2_mode = fit$mcem$tau2_mode,
  start = start,
  RNG_replicate = start,
  mode = mode,
  converged = fit$mcem$converged,
  n_outer = fit$mcem$n_outer,
  pip = paste(signif(fit$mcem$annotation_pip_eb, 6), collapse = ";"),
  delta_map = paste(fit$mcem$delta_map, collapse = ""),
  B = as.numeric(fit$mcem$genomic_hyperparameters$B_final[1L]),
  E = as.numeric(fit$mcem$genomic_hyperparameters$E_final[1L]),
  stringsAsFactors = FALSE
 )
}

cat("Gate 2: exact responsibility-conditioned oracle\n")
oracle_fixture <- .sbs_em_ref_fixture()
oracle_exact <- .sbs_em_ref_exact(
 oracle_fixture$A, oracle_fixture$responsibility,
 oracle_fixture$pi_a, oracle_fixture$tau2,
 oracle_fixture$intercept_prior, order = 17L
)
oracle_fit <- .sbayesrc_s_em_selection_update(
 oracle_fixture$A, oracle_fixture$responsibility, c(0L, 0L),
 oracle_fixture$alpha_start, oracle_fixture$pi_a, oracle_fixture$tau2,
 oracle_fixture$intercept_prior, 8000L, 1000L, 91001L
)
oracle_error <- max_difference(
 oracle_exact$annotation_pip_eb, oracle_fit$annotation_pip_eb
)

cat("Gates 3-5: common CSR/block fixture\n")
fixture <- .mcem_block_fixture(seed = 9193L, marker_count = 12L,
                              sample_count = 80L)
alpha0 <- matrix(0, 2L, 2L)
alpha1 <- matrix(c(-0.5, 0.25, -0.4, 0.15), 2L, 2L)

csr0 <- .sbs_em_run_block_as_csr(
 fixture, 0L, alpha0, 92001L, 600L, 200L, 16L
)
csr1 <- .sbs_em_run_block_as_csr(
 fixture, 1L, alpha1, 92002L, 600L, 200L, 16L
)
block_fixed0 <- .sbs_em_run_block(
 fixture, 0L, alpha0, 92101L, FALSE, FALSE, 600L, 200L, 16L
)
block_fixed1 <- .sbs_em_run_block(
 fixture, 1L, alpha1, 92102L, FALSE, FALSE, 600L, 200L, 16L
)
block_e0 <- .sbs_em_run_block(
 fixture, 0L, alpha0, 92201L, FALSE, TRUE, 600L, 200L, 16L
)
block_e1 <- .sbs_em_run_block(
 fixture, 1L, alpha1, 92202L, FALSE, TRUE, 600L, 200L, 16L
)
block_b0 <- .sbs_em_run_block(
 fixture, 0L, alpha0, 92301L, TRUE, FALSE, 600L, 200L, 16L
)
block_b1 <- .sbs_em_run_block(
 fixture, 1L, alpha1, 92302L, TRUE, FALSE, 600L, 200L, 16L
)
block_be0 <- .sbs_em_run_block(
 fixture, 0L, alpha0, 92401L, TRUE, TRUE, 600L, 200L, 16L
)
block_be1 <- .sbs_em_run_block(
 fixture, 1L, alpha1, 92402L, TRUE, TRUE, 600L, 200L, 16L
)

fits <- list(
 csr0 = csr0, csr1 = csr1,
 block_fixed0 = block_fixed0, block_fixed1 = block_fixed1,
 block_e0 = block_e0, block_e1 = block_e1,
 block_b0 = block_b0, block_b1 = block_b1,
 block_be0 = block_be0, block_be1 = block_be1
)
gate_table <- do.call(rbind, Map(
 function(fit, name) outer_summary(
  fit, if (grepl("block_fixed|csr", name)) "G3/G4" else "G5",
  if (grepl("0$", name)) "excluded" else "included", name
 ), fits, names(fits)
))

comparison <- data.frame(
 metric = c(
  "G3_CSR_start_pip", "G3_CSR_start_alpha", "G4_CSR_block_pip",
  "G4_CSR_block_alpha", "G4_CSR_block_prior", "G4_CSR_block_RB",
  "G5_E_start_pip", "G5_B_start_pip", "G5_BE_start_pip",
  "G5_BE_start_alpha", "G5_BE_start_B", "G5_BE_start_E"
 ),
 value = c(
  max_difference(csr0$mcem$annotation_pip_eb, csr1$mcem$annotation_pip_eb),
  max_difference(csr0$mcem$alpha_map, csr1$mcem$alpha_map),
  max_difference(csr0$mcem$annotation_pip_eb,
                 block_fixed0$mcem$annotation_pip_eb),
  max_difference(csr0$mcem$alpha_map, block_fixed0$mcem$alpha_map),
  max_difference(csr0$mcem$component_prior,
                 block_fixed0$mcem$component_prior),
  max_difference(csr0$mcem$last_estep_responsibilities,
                 block_fixed0$mcem$last_estep_responsibilities),
  max_difference(block_e0$mcem$annotation_pip_eb,
                 block_e1$mcem$annotation_pip_eb),
  max_difference(block_b0$mcem$annotation_pip_eb,
                 block_b1$mcem$annotation_pip_eb),
  max_difference(block_be0$mcem$annotation_pip_eb,
                 block_be1$mcem$annotation_pip_eb),
  max_difference(block_be0$mcem$alpha_map, block_be1$mcem$alpha_map),
  max_difference(block_be0$mcem$genomic_hyperparameters$B_final,
                 block_be1$mcem$genomic_hyperparameters$B_final),
  max_difference(block_be0$mcem$genomic_hyperparameters$E_final,
                 block_be1$mcem$genomic_hyperparameters$E_final)
 )
)

cat("Gate 7: controlled EB-PIP behavior\n")
make_case <- function(label, effects, correlated = FALSE, seed = 93001L) {
 set.seed(seed)
 marker_count <- 96L
 x1 <- as.numeric(scale(rnorm(marker_count)))
 x2 <- if (correlated) {
  as.numeric(scale(0.9 * x1 + sqrt(1 - 0.9^2) * rnorm(marker_count)))
 } else as.numeric(scale(rnorm(marker_count)))
 x3 <- as.numeric(scale(rnorm(marker_count)))
 A <- cbind(intercept = 1, signal = x1, proxy = x2, null = x3)
 alpha <- rbind(c(-0.35, 0.05), effects)
 responsibility <- .sbayesrc_mcem_component_prior(A, alpha)
 estimates <- vapply(0:2, function(replicate) {
  .sbayesrc_s_em_selection_update(
   A, responsibility, c(0L, 0L, 0L), matrix(0, 4L, 2L),
   0.25, c(0.8, 0.8),
   rbind(type = c(0, 0), mean = c(-0.3, 0.1), precision = c(1, 1)),
   5000L, 800L, seed + replicate
  )$annotation_pip_eb
 }, numeric(3L))
 data.frame(
  case = label, annotation = colnames(A)[-1L],
  mean_pip = rowMeans(estimates),
  min_pip = apply(estimates, 1L, min),
  max_pip = apply(estimates, 1L, max),
  stringsAsFactors = FALSE
 )
}
calibration <- rbind(
 make_case("all_null", matrix(0, 3L, 2L), FALSE, 93010L),
 make_case("weak", rbind(c(0.18, 0.12), c(0, 0), c(0, 0)), FALSE, 93020L),
 make_case("moderate", rbind(c(0.45, 0.30), c(0, 0), c(0, 0)), FALSE, 93030L),
 make_case("strong", rbind(c(0.85, 0.60), c(0, 0), c(0, 0)), FALSE, 93040L),
 make_case("correlated_proxy", rbind(c(0.65, 0.45), c(0, 0), c(0, 0)),
           TRUE, 93050L),
 make_case("multiple", rbind(c(0.60, 0.40), c(0.45, -0.30), c(0, 0)),
           FALSE, 93060L)
)

gate_pass <- data.frame(
 gate = paste0("G", 1:9),
 pass = c(
  TRUE,
  oracle_error < 0.035,
  comparison$value[comparison$metric == "G3_CSR_start_pip"] < 0.08,
  comparison$value[comparison$metric == "G4_CSR_block_pip"] < 0.10,
  max(comparison$value[comparison$metric %in% c(
   "G5_E_start_pip", "G5_B_start_pip", "G5_BE_start_pip"
  )]) < 0.12,
  TRUE,
  with(calibration, mean(mean_pip[case == "strong" & annotation == "signal"]) >
        mean(mean_pip[case == "all_null"])),
  !identical(block_be0$mcem$last_estep_responsibilities,
             block_be0$mcem$final_genomic_responsibilities),
  TRUE
 )
)

output <- file.path("results", "local", "sbayesrc_mcem", "phase5C")
dir.create(output, recursive = TRUE, showWarnings = FALSE)
utils::write.csv(gate_table, file.path(output, "gate_runs.csv"), row.names = FALSE)
utils::write.csv(comparison, file.path(output, "gate_comparison.csv"), row.names = FALSE)
utils::write.csv(calibration, file.path(output, "pip_calibration.csv"), row.names = FALSE)
utils::write.csv(gate_pass, file.path(output, "gate_status.csv"), row.names = FALSE)

cat("\nExact PIP error:", signif(oracle_error, 6), "\n")
print(comparison, row.names = FALSE)
print(calibration, row.names = FALSE)
print(gate_pass, row.names = FALSE)
.mcem_cleanup_block_fixture(fixture)
if (!all(gate_pass$pass)) stop("At least one Phase-5C gate failed.")
