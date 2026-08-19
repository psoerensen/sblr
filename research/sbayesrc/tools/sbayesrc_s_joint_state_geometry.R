# Development/reference implementation.
# Compact Phase-4C geometry audit for the validated internal CSR backend.
# Not a supported public sampler or API.

devtools::load_all(quiet = TRUE)
source(file.path("tests", "testthat", "helper-sbayesrc-s-genomic-reference.R"))

output_dir <- file.path(
  "results", "local", "sbayesrc_s_reference", "phase4C"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
fixture <- .sbs4b_fixture(160L, 20270930L)
initials <- list(c(0L, 0L, 0L), c(1L, 1L, 1L),
                 c(1L, 0L, 1L), c(0L, 1L, 0L))

run_configuration <- function(name, fixed_delta = integer()) {
  fits <- lapply(seq_along(initials), function(i) .sbs4b_run(
    fixture, 20271200L + i, 1800L, 400L,
    fixed_delta = fixed_delta, updateB = TRUE, updateE = FALSE,
    initial_delta = if (length(fixed_delta)) fixed_delta else initials[[i]]
  ))
  do.call(rbind, lapply(seq_along(fits), function(i) {
    trace <- fits[[i]]$chains[[1L]][[1L]]$convergence_trace
    alpha <- as.matrix(trace$alpha)
    delta <- as.matrix(trace$annotation_delta)
    active <- as.numeric(trace$realized_active_count)
    expected_active <- vapply(seq_len(nrow(alpha)), function(draw) {
      sum(stats::pnorm(drop(fixture$A %*% alpha[draw, 1:4])))
    }, numeric(1L))
    delta_switch <- c(FALSE, rowSums(abs(delta[-1L, , drop = FALSE] -
                                      delta[-nrow(delta), , drop = FALSE])) > 0)
    expected_jump <- c(NA_real_, abs(diff(expected_active)))
    active_jump <- c(NA_real_, abs(diff(active)))
    data.frame(
      configuration = name, chain = i,
      mean_expected_active = mean(expected_active),
      mean_realized_active = mean(active),
      expected_realized_correlation = stats::cor(expected_active, active),
      mean_expected_jump_at_delta_switch = if (any(delta_switch))
        mean(expected_jump[delta_switch], na.rm = TRUE) else NA_real_,
      mean_active_jump_at_delta_switch = if (any(delta_switch))
        mean(active_jump[delta_switch], na.rm = TRUE) else NA_real_,
      mean_expected_jump_otherwise = mean(expected_jump[!delta_switch], na.rm = TRUE),
      mean_active_jump_otherwise = mean(active_jump[!delta_switch], na.rm = TRUE),
      intercept_1 = mean(alpha[, 1L]),
      alpha_1_1 = mean(alpha[, 2L]),
      alpha_2_1 = mean(alpha[, 3L]),
      alpha_3_1 = mean(alpha[, 4L])
    )
  }))
}

geometry <- rbind(
  run_configuration("C1_fixed_delta_truth", c(1L, 1L, 0L)),
  run_configuration("C5_full")
)
write.csv(geometry, file.path(output_dir, "state_dependence.csv"), row.names = FALSE)
print(geometry)
