# Development/reference implementation.
# SBayesRC-S Phase 4B internal CSR genomic qualification.
# Not a supported public sampler or API.

devtools::load_all(quiet = TRUE)
source(file.path("tests", "testthat", "helper-sbayesrc-s-genomic-reference.R"))

annotation_summary <- function(fit) fit$chains[[1L]][[1L]]$annotation
component_summary <- function(fit) fit$component$prob[[1L]]

started <- proc.time()[["elapsed"]]
output_dir <- file.path("results", "local", "sbayesrc_s_reference", "phase4B")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
fixture <- .sbs4b_fixture(160L, 20270930L)
initials <- list(c(0L, 0L, 0L), c(1L, 1L, 1L),
                 c(1L, 0L, 1L), c(0L, 1L, 0L))
fits <- tryCatch(
  lapply(seq_along(initials), function(i) {
    .sbs4b_run(
      fixture, 20270930L + i, 1800L, 400L,
      initial_delta = initials[[i]], updateB = TRUE, updateE = FALSE
    )
  }),
  error = identity
)
if (inherits(fits, "error")) {
  blocked <- list(
    pass = FALSE,
    decision = "SBS4B-R3",
    route = "CSR summary-statistic BayesRC engine",
    marker_count = nrow(fixture$A),
    failure = conditionMessage(fits),
    interpretation = paste(
      "The genomic allocation state can contain an empty later stick.",
      "The validated Phase-3 flat intercept then has an improper conditional;",
      "no fallback or prior substitution was introduced."
    ),
    runtime_seconds = proc.time()[["elapsed"]] - started
  )
  saveRDS(blocked, file.path(output_dir, "phase4b_csr_qualification.rds"))
  writeLines(capture.output(str(blocked)),
             file.path(output_dir, "phase4b_blocked.txt"))
  print(blocked)
  stop("Phase 4B stopped at the empty-eligible-stick support boundary")
}
annotation <- lapply(fits, annotation_summary)
chain_pip <- do.call(rbind, lapply(annotation, function(x) {
  as.numeric(x$annotation_pip)
}))
colnames(chain_pip) <- colnames(fixture$annotation)
pooled_pip <- colMeans(chain_pip)
pip_range <- apply(chain_pip, 2L, function(x) diff(range(x)))
switches <- do.call(rbind, lapply(annotation, function(x) {
  as.numeric(x$annotation_switches)
}))
pi_A <- vapply(annotation, function(x) x$annotation_pi_A[[1L]], numeric(1L))
tau2 <- do.call(rbind, lapply(annotation, function(x) as.numeric(x$annotation_tau2)))
included <- vapply(annotation, function(x) {
  x$annotation_included_mean[[1L]]
}, numeric(1L))
active_count <- vapply(fits, function(x) sum(x$marker$dm[, 1L]), numeric(1L))
component_occupancy <- do.call(rbind, lapply(fits, function(x) {
  colSums(component_summary(x))
}))

# Fixed-hierarchy standard bridge already has exact automated coverage.  This
# longer comparison records the SNP-level spread between the learned selection
# hierarchy and the frozen continuous-alpha baseline without claiming equal
# posteriors.
standard <- .sbs4b_run_standard(
  fixture, 20270950L, 1800L, 400L, update_hierarchy = FALSE
)
mean_selection_bm <- Reduce(`+`, lapply(fits, function(x) x$marker$bm)) /
  length(fits)
mean_selection_dm <- Reduce(`+`, lapply(fits, function(x) x$marker$dm)) /
  length(fits)
snp_safeguards <- c(
  beta_correlation = stats::cor(
    as.numeric(mean_selection_bm), as.numeric(standard$marker$bm)
  ),
  pip_correlation = stats::cor(
    as.numeric(mean_selection_dm), as.numeric(standard$marker$dm)
  ),
  max_chain_beta_spread = max(apply(
    do.call(cbind, lapply(fits, function(x) x$marker$bm[, 1L])), 1L,
    function(x) diff(range(x))
  )),
  max_chain_pip_spread = max(apply(
    do.call(cbind, lapply(fits, function(x) x$marker$dm[, 1L])), 1L,
    function(x) diff(range(x))
  ))
)

# Correlated signal/proxy stress: exchangeability is not expected because the
# proxy is noisy, but both columns must remain index-coherent and explored.
proxy_fixture <- .sbs4b_fixture(160L, 20270960L)
proxy_fixture$annotation[, 3L] <- as.numeric(scale(
  0.92 * proxy_fixture$annotation[, 2L] +
    sqrt(1 - 0.92^2) * proxy_fixture$annotation[, 3L]
))
proxy_fixture$A[, 4L] <- proxy_fixture$annotation[, 3L]
proxy_fits <- list(
  .sbs4b_run(proxy_fixture, 20270961L, 1800L, 400L,
             initial_delta = c(0L, 0L, 1L), updateB = TRUE),
  .sbs4b_run(proxy_fixture, 20270962L, 1800L, 400L,
             initial_delta = c(0L, 1L, 0L), updateB = TRUE)
)
proxy_pip <- do.call(rbind, lapply(proxy_fits, function(x) {
  as.numeric(annotation_summary(x)$annotation_pip)
}))
proxy_switches <- do.call(rbind, lapply(proxy_fits, function(x) {
  as.numeric(annotation_summary(x)$annotation_switches)
}))

runtime_seconds <- proc.time()[["elapsed"]] - started
result <- list(
  route = "CSR summary-statistic BayesRC engine",
  marker_count = nrow(fixture$A),
  iterations = 1800L,
  burnin = 400L,
  chain_pip = chain_pip,
  pooled_pip = pooled_pip,
  pip_range = pip_range,
  switches = switches,
  pi_A_chain_mean = pi_A,
  tau2_chain_mean = tau2,
  included_chain_mean = included,
  active_count_chain_mean = active_count,
  component_occupancy = component_occupancy,
  snp_safeguards = snp_safeguards,
  proxy_pip = proxy_pip,
  proxy_switches = proxy_switches,
  runtime_seconds = runtime_seconds
)
result$pass <- isTRUE(
  max(pip_range) <= 0.10 &&
    all(rowSums(switches) > 0) &&
    all(is.finite(c(pi_A, tau2, included, active_count,
                    component_occupancy, snp_safeguards, proxy_pip))) &&
    snp_safeguards[["beta_correlation"]] >= 0.90 &&
    snp_safeguards[["pip_correlation"]] >= 0.80 &&
    all(rowSums(proxy_switches) > 0)
)

saveRDS(result, file.path(output_dir, "phase4b_csr_qualification.rds"))
write.csv(data.frame(annotation = names(pooled_pip), pooled_pip,
                     chain_range = pip_range),
          file.path(output_dir, "phase4b_annotation_pips.csv"), row.names = FALSE)
write.csv(data.frame(metric = names(snp_safeguards),
                     value = unname(snp_safeguards)),
          file.path(output_dir, "phase4b_snp_safeguards.csv"), row.names = FALSE)
print(result)
if (!result$pass) stop("Phase 4B qualification gates did not pass")
